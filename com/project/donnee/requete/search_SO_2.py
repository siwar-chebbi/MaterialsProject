#!/opt/anaconda3/bin/python

from pymatgen import MPRester

api = MPRester("eDCEK5m9WVjmajp7e8af")

compteur = 0

for i in range(1, 7):
    # Nous pouvons mettre (2,7) puisque nous recherchons 2 composés
    entries = api.get_entries({"nelements": i}, property_data=['elasticity'])
    for entry in entries:
        if entry.data["elasticity"]:
            oxygenFound = False
            soufreFound = False
            position = 0
            while position < len(entry.composition.elements) and not (oxygenFound and soufreFound):
                # element = entry.composition.elements[position]
                if entry.composition.elements[position].name == "O":  # element.name
                    oxygenFound = True
                if entry.composition.elements[position].name == "S":  # element.name
                    soufreFound = True
                position += 1

            if oxygenFound and soufreFound:
                compteur = compteur + 1

print("Le nombre d'éléments contenant S et O est : " + str(compteur))
