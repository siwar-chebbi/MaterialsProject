#!/opt/anaconda3/bin/python

from pymatgen import MPRester

api = MPRester("eDCEK5m9WVjmajp7e8af")

compteur = 0
compteuri = 0
texte = ""

for i in range(1, 6):
    compteuri = 0
    entries = api.get_entries({"nelements": i}, property_data=["elasticity"])
    for entry in entries:
        if entry.data["elasticity"]:
            compteuri = compteuri + 1
            compteur = compteur + 1

    texte += "Number of elements containing " + str(i) + " atome(s) : " + str(compteuri) + "\n"

texte += "\nThe total number of all elements : " + str(compteur)
print(texte)
