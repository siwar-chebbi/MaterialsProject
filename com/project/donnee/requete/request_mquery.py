#!/opt/anaconda3/bin/python

from pymatgen import MPRester
from pymatgen.core.periodic_table import Element
api = MPRester("fB610TDF3LSwxiN9")


actinoid = list()
alkali = list()
alkaline = list()
chalcogen = list()
halogen = list()
lanthanoid = list()
metalloid = list()
noble_gas = list()
quadrupolar = list()
rare_earth_metal = list()
transition_metal = list()
tableau = {'actinoid': actinoid, 'alkali': alkali, 'alkaline': alkaline, 'chalcogen': chalcogen, 'halogen': halogen, 'lanthanoid': lanthanoid, 'metalloid': metalloid, 'noble_gas': noble_gas, 'quadrupolar': quadrupolar, 'rare_earth_metal': rare_earth_metal, 'transition_metal': transition_metal}

for element in Element:
    if element.is_actinoid : actinoid.append(element.value)
    if element.is_alkali : alkali.append(element.value)
    if element.is_alkaline : alkaline.append(element.value)
    if element.is_chalcogen : chalcogen.append(element.value)
    if element.is_halogen : halogen.append(element.value)
    if element.is_lanthanoid : lanthanoid.append(element.value)
    if element.is_metalloid : metalloid.append(element.value)
    if element.is_noble_gas : noble_gas.append(element.value)
    if element.is_quadrupolar : quadrupolar.append(element.value)
    if element.is_rare_earth_metal : rare_earth_metal.append(element.value)
    if element.is_transition_metal : transition_metal.append(element.value)
def print_group (group):
    print("les éléments de " + str.upper(group) + " sont : "+ str(tableau.get(group)) + "\n" )


def request (group):
    results = api.query({'elements': {'$in': tableau.get(group)}, 'nelements': {'$lte': 6}, 'elasticity': {'$ne': None}}, properties=['pretty_formula', 'elasticity'])
    N = len(results)
    print("\nLe nombre d'élément de la famille " + str(group) + " est de : " + str(N))
    return N

for key in tableau.keys():
    print_group (key)
somme=0
for key in tableau.keys():
    request(key)
    somme=somme+request(key)
print(str(somme))