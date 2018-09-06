#!/opt/anaconda3/bin/python

from pymatgen import MPRester

api = MPRester("fB610TDF3LSwxiN9")

compteur = 0
compteuri = 0
texte = ""


for i in range(1, 7):
    compteuri = 0
    entries = api.get_entries({"nelements": i}, property_data=["elasticity"])
    for entry in entries:
        if entry.data["elasticity"]:
            compteuri = compteuri + 1
            compteur = compteur + 1

    texte += "Le nombre de molécules contenant " + str(i) + " atome(s) : " + str(compteuri) + "\n"


texte += "\nLe nombre total d'éléments : " + str(compteur)
print(texte)






