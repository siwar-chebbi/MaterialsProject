#!/opt/anaconda3/bin/python

from pymatgen import MPRester

api = MPRester("eDCEK5m9WVjmajp7e8af")


texte = ""

for i in range(1, 8):
    entries = api.query({"nelements": i, "elasticity": {"$exists": True}}, ["task_id"])
    a = len(entries)
    texte += "Number of elements containing " + str(i) + " atome(s) : " + str(a) + "\n"

B = str(len(api.query({"nelements": {"$in": [1,2,3,4,5,6,7]}, "elasticity": {"$exists": True}}, ["task_id"])))


print(texte)
print(B)

#     a = len(entries)
#     texte += "Number of elements containing " + str(i) + " atome(s) : " + str(a) + "\n"
#
# B = str(len(api.query({"nelements": {"$in": [1,2,3,4,5,6,7]}, "elasticity": {"$exists": True}}, ["task_id"])))
#
#
print(str(len(entries)))
# print(B)