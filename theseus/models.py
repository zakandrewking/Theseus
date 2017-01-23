# -*- coding: utf-8 -*-

import cobra
import cobra.io
from cobra.core.Formula import Formula
import os
from os.path import join, abspath, dirname
import re
import pickle
from six import iterkeys

try:
    from cobrame import MetabolicReaction, StoichiometricData
    from cobrame.solve.symbolic import compile_expressions
except ImportError:
    pass

data_path = join(abspath(dirname(__file__)), 'data')

def get_model_list():
    """Get the models that are available, as SBML, in data/models"""
    return [x.replace('.xml', '').replace('.mat', '').replace('.json', '') for x in
            os.listdir(join(data_path, 'models'))
            if '.xml' in x or '.mat' in x or '.json' in x]

def check_for_model(name):
    """Check for model, case insensitive, and ignore periods and underscores"""
    def min_name(n):
        return n.lower().replace('.','').replace(' ','').replace('_','')
    for x in get_model_list():
        if min_name(name)==min_name(x):
            return x
    return None

# the regex to separate the base id, the chirality ('_L') and the compartment ('_c')
reg = re.compile(r'(.*?)(?:(.*[^_])_([LDSR]))?[_\(\[]([a-z])[_\)\]]?$')
def id_for_new_id_style(old_id, is_metabolite=False, new_id_style='cobrapy'):
    """ Get the new style id"""

    def join_parts(the_id, the_compartment):
        if (new_id_style.lower()=='cobrapy'):
            if the_compartment:
                the_id = the_id+'_'+the_compartment
            the_id = the_id.replace('-', '__')
        elif (new_id_style.lower()=='simpheny'):
            if the_compartment and is_metabolite:
                the_id = the_id+'['+the_compartment+']'
            elif the_compartment:
                the_id = the_id+'('+the_compartment+')'
            the_id = the_id.replace('__', '-')
        else:
            raise Exception('Invalid id style')
        return the_id

    # separate the base id, the chirality ('_L') and the compartment ('_c')
    m = reg.match(old_id)
    if m is None:
        # still change the underscore/dash
        new_id = join_parts(old_id, None)
    elif m.group(2) is None:
        new_id = join_parts(m.group(1), m.group(4))
    else:
        # if the chirality is not joined by two underscores, then fix that
        a = "__".join(m.groups()[1:3])
        new_id = join_parts(a, m.group(4))

    # deal with inconsistent notation of (sec) vs. [sec] in iJO1366 versions
    new_id = new_id.replace('[sec]', '_sec_').replace('(sec)', '_sec_')

    return new_id

def convert_ids(model, new_id_style):
    """Converts metabolite and reaction ids to the new style. Style options:

    cobrapy: EX_lac__L_e
    simpheny: EX_lac-L(e)

    """
    # loop through the ids:

    # this code comes from cobra.io.sbml
    # legacy_ids add special characters to the names again
    for metabolite in model.metabolites:
        metabolite.id = fix_legacy_id(metabolite.id, use_hyphens=False)
    model.metabolites._generate_index()
    for reaction in model.reactions:
        reaction.id = fix_legacy_id(reaction.id, use_hyphens=False)
    model.reactions._generate_index()
    # remove boundary metabolites (end in _b and only present in exchanges) . Be
    # sure to loop through a static list of ids so the list does not get
    # shorter as the metabolites are deleted
    for metabolite_id in [str(x) for x in model.metabolites]:
        metabolite = model.metabolites.get_by_id(metabolite_id)
        if not metabolite.id.endswith("_b"):
            continue
        for reaction in list(metabolite._reaction):
            if reaction.id.startswith("EX_"):
                metabolite.remove_from_model()
                break
    model.metabolites._generate_index()

    # separate ids and compartments, and convert to the new_id_style
    for reaction in model.reactions:
        reaction.id = id_for_new_id_style(reaction.id, new_id_style=new_id_style)
    model.reactions._generate_index()
    for metabolite in model.metabolites:
        metabolite.id = id_for_new_id_style(metabolite.id, is_metabolite=True, new_id_style=new_id_style)
    model.metabolites._generate_index()

    return model

def load_model_me(unmodified_me=False):
    import cloudpickle

    def me_no_glucose_ex(me):
        me.reactions.get_by_id('EX_glc__D_e').lower_bound = 0
        return me

    with open(join(data_path, 'models', 'prototype_67.pickle'), 'rb') as f:
        me = pickle.load(f)

    return me if unmodified_me else me_no_glucose_ex(me)

def load_model(name, id_style='cobrapy', unmodified_me=False):
    """Load a model, and give it a particular id style"""

    if name == 'ME':
        me = load_model_me(unmodified_me)
        me.id = name
        return me

    # check for model
    name = check_for_model(name)
    if not name:
        raise Exception('Could not find model')

    # load the model pickle, or, if not, the sbml
    try:
        with open(join(data_path, 'model_pickles', name+'.pickle'), 'rb') as f:
            model = pickle.load(f)
    except:
        try:
            model = cobra.io.load_matlab_model(join(data_path, 'models', name+'.mat'))
        except:
            try:
                model = cobra.io.read_sbml_model(join(data_path, 'models', name+'.xml'))
            except:
                model = cobra.io.load_json_model(join(data_path, 'models', name+'.json'))
        with open(join(data_path, 'model_pickles', name+'.pickle'), 'w') as f:
            pickle.dump(model, f)

    # convert the ids
    model = convert_ids(model, id_style)

    # extract metabolite formulas from names (e.g. for iAF1260)
    model = get_formulas_from_names(model)

    # turn off carbon sources
    model = turn_off_carbon_sources(model)

    return model

def get_formulas_from_names(model):
    reg = re.compile(r'.*_([A-Za-z0-9]+)$')
    for metabolite in model.metabolites:
        if (metabolite.formula is not None
            and str(metabolite.formula).strip() != ''): continue
        m = reg.match(metabolite.name)
        if m:
            metabolite.formula = Formula(m.group(1))
    return model

def turn_off_carbon_sources(model):
    for reaction in model.reactions:
        if not reaction.id.startswith('EX_'): continue
        if carbons_for_exchange_reaction(reaction) > 0:
            reaction.lower_bound = 0
    return model

def setup_model(model, substrate_reactions, aerobic=True, sur=10, max_our=10):
    """Set up the model with environmntal parameters.

    model: a cobra model
    substrate_reactions: A single reaction id, list of reaction ids, or dictionary with reaction
    ids as keys and max substrate uptakes as keys. If a list or single id is
    given, then each substrate will be limited to /sur/
    aerobic: True or False
    sur: substrate uptake rate. Ignored if substrate_reactions is a dictionary.
    max_our: Max oxygen uptake rate.
    id_style: 'cobrapy' or 'simpheny'.

    """

    if isinstance(substrate_reactions, dict):
        for r, v in substrate_reactions.iteritems():
            model.reactions.get_by_id(r).lower_bound = -abs(v)
    elif isinstance(substrate_reactions, list):
        for r in substrate_reactions:
            model.reactions.get_by_id(r).lower_bound = -abs(sur)
    elif isinstance(substrate_reactions, str):
        model.reactions.get_by_id(substrate_reactions).lower_bound = -abs(sur)
    else: raise Exception('bad substrate_reactions argument')

    o2 = 'EX_o2_e'
    if aerobic:
        model.reactions.get_by_id(o2).lower_bound = -abs(max_our)
    else:
        model.reactions.get_by_id(o2).lower_bound = 0

    # model specific setup
    if str(model) == 'iJO1366' and aerobic is False:
        for r in ['CAT', 'SPODM', 'SPODMpp']:
            model.reactions.get_by_id(r).lower_bound = 0
            model.reactions.get_by_id(r).upper_bound = 0

    elif str(model) == 'iJR904':
        model.objective = {model.reactions.get_by_id('BIOMASS_Ecoli'): 1}

    elif str(model) == 'iMM904' and aerobic is False:
        necessary_ex = ['EX_ergst_e', 'EX_zymst_e', 'EX_hdcea_e', 'EX_ocdca_e',
                        'EX_ocdcea_e', 'EX_ocdcya_e']
        for r in necessary_ex:
            rxn = model.reactions.get_by_id(r)
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000

    return model

def turn_on_subsystem(model, subsytem):
    raise NotImplementedError()
    for reaction in model.reactions:
        if reaction.subsystem.strip('_') == subsytem.strip('_'):
            reaction.lower_bound = -1000 if reaction.reversibility else 0
            reaction.upper_bound = 1000
    return model

def carbons_for_exchange_reaction(reaction):
    if len(reaction._metabolites) > 1:
        raise Exception('%s not an exchange reaction' % str(reaction))

    metabolite = next(iterkeys(reaction._metabolites))
    try:
        return metabolite.elements['C']
    except KeyError:
        return 0

def add_me_reaction(model, reaction_id, stoichiometry, bounds=(-1000.0, 1000.0),
                    keff=65.0):
    """Add a new reaction to the ME model.

    stoichiometry: {metabolite_id: coefficient}

    """
    data = StoichiometricData(reaction_id, model)
    data.lower_bound, data.upper_bound = bounds
    data._stoichiometry = stoichiometry

    for rev_str, reverse in ('FWD', False), ('REV', True):
        reaction = MetabolicReaction('%s_%s_CPLX_dummy' % (reaction_id, rev_str))
        reaction.keff = keff
        reaction.stoichiometric_data = data
        reaction.reverse = reverse
        reaction.complex_data = model.complex_data.CPLX_dummy
        model.add_reaction(reaction)
        reaction.update()

def add_pathway(model, new_metabolites, new_reactions, subsystems, bounds,
                check_mass_balance=False, check_charge_balance=False,
                ignore_repeats=False, recompile_expressions=True):
    """Add a pathway to the model. Reversibility defaults to reversible (1).

    check_charge_balance: Only works if check_mass_balance is True.

    recompile_expressions: If True, then recompile_expressions when new ME
    reactions are added.


    new_metabolites: e.g. { 'ggpp_c': {'formula': 'C20H33O7P2', 'name': 'name'},
                            'phyto_c': {'formula': 'C40H64'}},
                            'lyco_c': {'formula': 'C40H56'},
                            'lyco_e': {'formula': 'C40H56'} }
    new_reactions: e.g. { 'FPS': { 'ipdp_c': -2,
                                   'ppi_c': 1,
                                   'grdp_c': 1 },
                          'CRTE': { 'ipdp_c': -1,
                                    'frdp_c': -1,
                                    'ggpp_c': 1,
                                    'ppi_c': 1 } }
    subsystems: e.g. { 'FPS': 'Lycopene production',
                       'CRTE': 'Lycopene production' }
    bound: e.g. { 'FPS': (0, 0),
                  'CRTE': (0, 1000) }

    """

    for k, v in new_metabolites.iteritems():
        m = cobra.Metabolite(id=k,
                             formula=v.get('formula', None),
                             name=v.get('name', None))
        m.charge = v.get('charge', None)
        try:
            model.add_metabolites([m])
        except Exception as err:
            if (not ignore_repeats or
                "already in the model" not in str(err)):
                raise(err)

    has_new_me_reactions = False
    for name, mets in new_reactions.iteritems():
        if bounds and (name in bounds):
            r_bounds = bounds[name]
        else:
            r_bounds = (-1000, 1000)

        if model.id == 'ME' and not name.startswith('EX_'):
            # me reaction
            if ignore_repeats and ('%_FWD_CPLX_dummy' % name in model.reactions or
                                   '%_REV_CPLX_dummy' % name in model.reactions):
                continue
            print('Adding ME reactions for %s' % name)
            add_me_reaction(model, name, mets, r_bounds)
            has_new_me_reactions = True
        else:
            # m reaction
            if ignore_repeats and name in model.reactions:
                continue
            r = cobra.Reaction(name)
            m_obj = {}
            for k, v in mets.iteritems():
                m_obj[model.metabolites.get_by_id(k)] = v
            r.add_metabolites(m_obj)
            r.lower_bound, r.upper_bound = r_bounds
            if subsystems and (name in subsystems):
                r.subsystem = subsystems[name]
            model.add_reaction(r)

            # mass balance
            if check_mass_balance and 'EX_' not in name:
                balance = model.reactions.get_by_id(name).check_mass_balance()
                if not check_charge_balance and 'charge' in balance:
                    del balance['charge']
                if len(balance) > 0:
                    raise Exception('Bad balance: %s' % str(balance))

    # recompile the expressions
    if has_new_me_reactions and recompile_expressions:
        print('Recompiling expressions')
        model.expressions = compile_expressions(model)

    return model


def fix_legacy_id(id, use_hyphens=False):
    id = id.replace('_DASH_', '__')
    id = id.replace('_FSLASH_', '/')
    id = id.replace('_BSLASH_', "\\")
    id = id.replace('_LPAREN_', '(')
    id = id.replace('_LSQBKT_', '[')
    id = id.replace('_RSQBKT_', ']')
    id = id.replace('_RPAREN_', ')')
    id = id.replace('_COMMA_', ',')
    id = id.replace('_PERIOD_', '.')
    id = id.replace('_APOS_', "'")
    id = id.replace('&amp;', '&')
    id = id.replace('&lt;', '<')
    id = id.replace('&gt;', '>')
    id = id.replace('&quot;', '"')
    if use_hyphens:
        id = id.replace('__', '-')
    else:
        id = id.replace("-", "__")
    return id
