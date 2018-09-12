#!/opt/anaconda3/bin/python

from pymatgen import MPRester

api = MPRester("fB610TDF3LSwxiN9")
results = api.query({"elements": {'$all': ['S', 'O']}, 'nelements': {'$lte': 6}, 'elasticity': {'$ne': None}}, properties=['pretty_formula', 'elasticity'])
N = len(results)

print("Le nombre d'éléments contenant S et O est : " + str(N))
