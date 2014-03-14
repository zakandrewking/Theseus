# -*- coding: utf-8 -*-
import cobra
import cobra.io
from cobra.io.sbml import fix_legacy_id
import os
from os.path import join, abspath, dirname
import re
import cPickle as pickle

data_path = join(abspath(dirname(__file__)), 'data')

def get_model_list():
    """Get the models that are available, as SBML, in data/models"""
    return [x.replace('.xml', '') for x in
            os.listdir(join(data_path, 'models'))
            if '.xml' in x]

def check_for_model(name):
    """Check for model, case insensitive, and ignore periods and underscores"""
    def min_name(n):
        return n.lower().replace('.','').replace(' ','').replace('_','')
    for x in get_model_list():
        if min_name(name)==min_name(x):
            return x
    return None

def convert_ids(model, new_id_style):
    """Converts metabolite and reaction ids to the new style. Style options:

    cobrapy: EX_lac__L_e
    simpheny: EX_lac-L(e)

    """

    # the regex to separate the base id, the chirality ('_L') and the compartment ('_c')
    reg = re.compile(r'(.*?)(?:(.*[^_])_([LD]))?[_\(\[]([a-z])[_\)\]]?$')
    
    def id_for_new_id_style(old_id, is_metabolite=False):
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
        return new_id

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
        reaction.id = id_for_new_id_style(reaction.id)
    model.reactions._generate_index()
    for metabolite in model.metabolites:
        metabolite.id = id_for_new_id_style(metabolite.id, is_metabolite=True)
    model.metabolites._generate_index()
    
    return model

def load_model(name, id_style='cobrapy'):
    """Load a model, and give it a particular id style"""

    # check for model
    name = check_for_model(name)
    if not name:
        raise Exception('Could not find model')

    # load the model pickle, or, if not, the sbml
    try:
        with open(join(data_path, 'model_pickles', name+'.pickle'), 'r') as f:
            model = pickle.load(f)
    except:
        model = cobra.io.read_sbml_model(join(data_path, 'models', name+'.xml'))
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
        if metabolite.formula!='': continue
        m = reg.match(metabolite.name)
        if m:
            metabolite.formula = m.group(1)
    return model

def turn_off_carbon_sources(model):
    for reaction in model.reactions:
        if 'EX_' not in str(reaction): continue
        if carbons_for_exchange_reaction(reaction) > 0:
            reaction.lower_bound = 0
    return model

def setup_model(model, substrate_reactions, aerobic=True, sur=10, max_our=10, id_style='cobrapy'):
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
    if id_style=='cobrapy': o2 = 'EX_o2_e'
    elif id_style=='simpheny': o2 = 'EX_o2(e)'
    else: raise Exception('Invalid id_style')

    if isinstance(substrate_reactions, dict):
        for r, v in substrate_reactions.iteritems():
            model.reactions.get_by_id(r).lower_bound = -abs(v)
    elif isinstance(substrate_reactions, list):
        for r in substrate_reactions:
            model.reactions.get_by_id(r).lower_bound = -abs(sur)
    elif isinstance(substrate_reactions, str):                
        model.reactions.get_by_id(substrate_reactions).lower_bound = -abs(sur)
    else: raise Exception('bad substrate_reactions argument')
        
    if aerobic:
        model.reactions.get_by_id(o2).lower_bound = -abs(max_our)
    else:
        model.reactions.get_by_id(o2).lower_bound = 0
        
    # model specific setup
    if str(model)=='iJO1366' and aerobic==False:
        for r in ['CAT', 'SPODM', 'SPODMpp']:
            model.reactions.get_by_id(r).lower_bound = 0
            model.reactions.get_by_id(r).upper_bound = 0

    # TODO hydrogen reaction for ijo
            
    if str(model)=='iMM904' and aerobic==False:
        necessary_ex = ['EX_ergst(e)', 'EX_zymst(e)', 'EX_hdcea(e)',
                        'EX_ocdca(e)', 'EX_ocdcea(e)', 'EX_ocdcya(e)']
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
    
    metabolite = reaction._metabolites.iterkeys().next()
    match = re.match(r'C([0-9]+)', str(metabolite.formula))
    try:
        return int(match.group(1))
    except AttributeError:
        return 0
