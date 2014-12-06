from sys import argv, exit
from cobra.io import save_json_model
from theseus import load_model

try:
    name = argv[1]
except KeyError:
    print('Usage python -m theseus.save_json iJO1366')
    exit()

model = load_model(name)
    
save_json_model(model, name + '.json')
