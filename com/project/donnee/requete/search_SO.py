#!/opt/anaconda3/bin/python

from pymatgen import MPRester

api = MPRester("fB610TDF3LSwxiN9")

compteur = 0

for i in range(1, 7):
    entries = api.get_entries({"nelements": i}, property_data=['elasticity'])
    for entry in entries:
        if entry.data["elasticity"]:
            oxygen = False
            soufre = False
            for element in entry.composition.elements:
                if element.name == "O":
                    oxygen = True
                if element.name == "S":
                    soufre = True
            if oxygen and soufre:
                compteur = compteur + 1

print("Le nombre d'éléments contenant S et O est : " + str(compteur))
