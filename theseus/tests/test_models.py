from theseus.models import *

import cobra
import os
import pytest

def test_get_model_list():
    model_list = get_model_list()
    assert 'iJO1366' in model_list
    assert 'iAF1260' in model_list
    assert 'E coli core' in model_list

def test_check_for_model():
    assert check_for_model('e_coli_core')=='E coli core'
    assert check_for_model('E. coli core')=='E coli core'

def test_load_model():
    model = load_model('iJO1366')
    assert isinstance(model, cobra.core.Model)

    model = load_model('RECON1')
    assert isinstance(model, cobra.core.Model)


def test_load_model_json():
    model = load_model('iML1503')
    assert isinstance(model, cobra.core.Model)


def test_id_for_new_id_style():
    assert(id_for_new_id_style('EX_glc(r)', is_metabolite=False)=='EX_glc_r')
    assert(id_for_new_id_style('glucose[r]', is_metabolite=True)=='glucose_r')

def test_convert_ids():
    for model_name in 'iJO1366', 'iAF1260', 'E coli core':
        print "\n"
        print model_name
        model = load_model(model_name)

        # cobrapy style
        model = convert_ids(model, new_id_style='cobrapy')
        print 'cobrapy ids'
        print [str(x) for x in model.reactions if 'lac' in str(x)]
        assert 'EX_lac__D_e' in [str(x) for x in model.reactions]
        assert ['-' not in str(x) for x in model.reactions]
        print [str(x) for x in model.metabolites if 'lac' in str(x)]
        assert 'lac__D_e' in [str(x) for x in model.metabolites]
        assert ['-' not in str(x) for x in model.metabolites]

        # simpheny style
        model = convert_ids(model, new_id_style='simpheny')
        print 'simpheny ids'
        print [str(x) for x in model.reactions if 'lac' in str(x)]
        assert 'EX_lac-D(e)' in [str(x) for x in model.reactions]
        assert ['__' not in str(x) for x in model.reactions]
        print [str(x) for x in model.metabolites if 'lac' in str(x)]
        assert 'lac-D[e]' in [str(x) for x in model.metabolites]
        assert ['__' not in str(x) for x in model.metabolites]

def test_get_formulas_from_names():
    model = load_model('iAF1260')
    assert model.metabolites.get_by_id('acald_c').formula.elements['C'] == 2
    assert model.metabolites.get_by_id('acald_c').formula.elements['H'] == 4
    assert model.metabolites.get_by_id('acald_c').formula.elements['O'] == 1

def test_turn_off_carbon_sources():
    for model_name in 'iJO1366', 'iAF1260', 'E coli core':
        print model_name
        model = load_model(model_name)
        model.reactions.get_by_id('EX_glc_e').lower_bound = -100
        model.reactions.get_by_id('EX_ac_e').lower_bound = -100
        model = turn_off_carbon_sources(model)
        assert model.reactions.get_by_id('EX_glc_e').lower_bound==0
        assert model.reactions.get_by_id('EX_ac_e').lower_bound==0

def test_setup_model():
    model = load_model('iJO1366')
    model = setup_model(model, 'EX_glc_e', aerobic=False, sur=18)
    assert model.reactions.get_by_id('EX_glc_e').lower_bound==-18
    assert model.reactions.get_by_id('EX_o2_e').lower_bound==0
    assert model.reactions.get_by_id('CAT').upper_bound==0
    assert model.reactions.get_by_id('SPODM').upper_bound==0

    model = setup_model(model, ['EX_glc_e', 'EX_xyl__D_e'], sur=-18)
    assert model.reactions.get_by_id('EX_glc_e').lower_bound==-18
    assert model.reactions.get_by_id('EX_xyl__D_e').lower_bound==-18

    model = setup_model(model, {'EX_glc_e':-5, 'EX_xyl__D_e':5}, sur=999)
    assert model.reactions.get_by_id('EX_glc_e').lower_bound==-5
    assert model.reactions.get_by_id('EX_xyl__D_e').lower_bound==-5

def test_turn_on_subsystem():
    with pytest.raises(NotImplementedError):
        turn_on_subsystem(None, None)

def test_carbons_for_exchange_reaction():
    model = load_model('iJO1366')
    assert carbons_for_exchange_reaction(model.reactions.get_by_id('EX_glc_e'))==6
    assert carbons_for_exchange_reaction(model.reactions.get_by_id('EX_ac_e'))==2
    assert carbons_for_exchange_reaction(model.reactions.get_by_id('EX_for_e'))==1

def test_add_pathway():
    new = [ { '1poh_c': {'formula': 'C3H8O', 'name': '1-propanol'} },
            { '2OBUTDC': {'2obut_c': -1, 'h_c': -1, 'ppal_c': 1, 'co2_c': 1},
              '1PDH': {'ppal_c': -1, 'nadh_c': -1, 'h_c': -1, '1poh_c': 1, 'nad_c': 1},
              'EX_1poh_e': {'1poh_c': -1} },
            { 'EX_1poh_e': '1-propanol production' },
            { 'EX_1poh_e': (0, 1000) } ]

    m = load_model('iJO1366')
    model = add_pathway(m.copy(), *new, check_mass_balance=True)
    assert isinstance(model.metabolites.get_by_id(new[0].keys()[0]), cobra.Metabolite)
    reaction = model.reactions.get_by_id('EX_1poh_e')
    assert isinstance(model.reactions.get_by_id(new[1].keys()[0]), cobra.Reaction)
    assert reaction.reversibility == False
    assert reaction.upper_bound == 1000
    assert reaction.lower_bound == 0

    # test metabolite name
    assert model.metabolites.get_by_id('1poh_c').name == '1-propanol'
    # test mass balance
    assert len(model.reactions.get_by_id('1PDH').check_mass_balance()) == 0

    # test ignore repeats
    model = add_pathway(m.copy(), *new)
    with pytest.raises(Exception):
        model = add_pathway(model, *new)
    model = add_pathway(model, *new, ignore_repeats=True)
