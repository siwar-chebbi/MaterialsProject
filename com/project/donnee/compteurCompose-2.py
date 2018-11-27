#!/opt/anaconda3/bin/python

from pymatgen import MPRester

api = MPRester("eDCEK5m9WVjmajp7e8af")

texte = ""


for i in range(1, 7):
    entries = api.query({"nelements": i, "elasticity": {"$exists": True}}, ["material_id"])

    texte += "Number of elements containing " + str(i) + " atome(s) : " + str(len(api.query({"nelements": i, "elasticity": {"$exists": True}}, ["material_id"]))) + "\n"

print(texte)






