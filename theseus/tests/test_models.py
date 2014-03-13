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

def test_turn_on_subsystem():
    with pytest.raises(NotImplementedError):
        turn_on_subsystem(None, None)

def test_carbons_for_exchange_reaction():
    model = load_model('iJO1366')
    assert carbons_for_exchange_reaction(model.reactions.get_by_id('EX_glc_e'))==6
