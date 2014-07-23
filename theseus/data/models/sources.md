# iAF1260
file: Ec_iAF1260.xml
publication: 10.1038/msb4100155

# iAF1260b
file: Ec_iAF1260b.xml
publication: 10.1016/j.ymben.2009.10.003
generated from iAF1260, using the following code:

```python
new_reactions = { 'DHORDfum': {'dhor__S_c': -1, 'fum_c': -1, 'orot_c': 1, 'succ_c': 1},
                  'MALt3pp': {'mal__L_c': -1, 'h_p': -1, 'h_c': 1, 'mal__L_p': 1},
                  'ALAt2rpp': {'ala__L_p': -1, 'h_p': -1, 'ala__L_c': 1, 'h_c': 1},
                  'GLYt2rpp': {'gly_p': -1, 'h_p': -1, 'gly_c': 1, 'h_c': 1},
				               'CITt3pp': {'cit_c': -1, 'h_p': -1, 'h_c': 1, 'cit_p': 1},
                  'ASPt2rpp': {'asp__L_p': -1, 'h_p': -1, 'asp__L_c': 1, 'h_c': 1}
                 }
bounds = {'DHORDfum': (0, 1000),
          'MALt3pp': (0, 1000),
          'CITt3pp': (0, 1000)}
iAF1260b = theseus.add_pathway(iAF1260.copy(), {}, new_reactions, {}, bounds, check_mass_balance=True)
```						 

# iJO1366
file: msb201165-sup-0003.xml
publication: 10.1038/msb.2011.65

# E coli core
file: e_coli_core.xml
publication: http://systemsbiology.ucsd.edu/Downloads

# RECON1
file: RECON1.xml
source: bigg.ucsd.edu

# iMM904
file: iMM904.xml
publication: 10.1186/1752-0509-3-37

# iJR904
file: iJO904.xml
source: bigg.ucsd.edu
TODO: fix model.id

# iJE660b
file: iJE660b.mat
source: Simpheny database
NOTE: currently has not biomass function