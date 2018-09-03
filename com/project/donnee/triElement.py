#!/opt/anaconda3/bin/python
from collections import OrderedDict

from pymatgen import MPRester

api = MPRester("fB610TDF3LSwxiN9")
tableau = {'Ac': 0, 'Al': 0, 'Am': 0, 'Sb': 0, 'Ag': 0, 'Ar': 0, 'As': 0, 'At': 0, 'N': 0, 'Ba': 0, 'Bk': 0, 'Be': 0, 'Bi': 0, 'Bh': 0, 'B': 0, 'Br': 0, 'Cd': 0, 'Ca': 0, 'Cf': 0, 'C': 0, 'Ce': 0, 'Cs': 0, 'Cl': 0, 'Cr': 0, 'Co': 0, 'Cu': 0, 'Cm': 0, 'Ds': 0, 'Db': 0, 'Dy': 0, 'Es': 0, 'Er': 0, 'Sn': 0, 'Eu': 0, 'Fe': 0, 'Fm': 0, 'F': 0, 'Fr': 0, 'Gd': 0, 'Ga': 0, 'Ge': 0, 'Hf': 0, 'Hs': 0, 'He': 0, 'Ho': 0, 'H': 0, 'In': 0, 'I': 0, 'Ir': 0, 'Kr': 0, 'La': 0, 'Lr': 0, 'Li': 0, 'Lu': 0, 'Mg': 0, 'Mn': 0, 'Mt': 0, 'Md': 0, 'Hg': 0, 'Mo': 0, 'Nd': 0, 'Ne': 0, 'Np': 0, 'Ni': 0, 'Nb': 0, 'No': 0, 'Os': 0, 'Au': 0, 'O': 0, 'Pd': 0, 'P': 0, 'Pt': 0, 'Pb': 0, 'Pu': 0, 'Po': 0, 'K': 0, 'Pr': 0, 'Pm': 0, 'Pa': 0, 'Ra': 0, 'Rn': 0, 'Re': 0, 'Rh': 0, 'Rb': 0, 'Ru': 0, 'Rf': 0, 'Sm': 0, 'Sc': 0, 'Sg': 0, 'Se': 0, 'Si': 0, 'Na': 0, 'Sr': 0, 'S': 0, 'Ta': 0, 'Tc': 0, 'Te': 0, 'Tb': 0, 'Tl': 0, 'Th': 0, 'Tm': 0, 'Ti': 0, 'W': 0, 'Uub': 0, 'Uuh': 0, 'Uuo': 0, 'Uup': 0, 'Uuq': 0, 'Uus': 0, 'Uut': 0, 'Uuu': 0, 'U': 0, 'V': 0, 'Xe': 0, 'Yb': 0, 'Y': 0, 'Zn': 0, 'Zr': 0}

# entries = api.get_entries_in_chemsys(['Ca', 'C', 'O'], property_data=['elasticity'])
texte = ""
#testgit

for i in range(1, 7):
    entries = api.get_entries({"nelements": i}, property_data=['elasticity'])
    for entry in entries:
        if entry.data["elasticity"]:

            for element in entry.composition.elements:

                tableau[element.name] = tableau.get(element.name)+1


tableauTrie = OrderedDict(sorted(tableau.items(), reverse=True,  key=lambda t: t[1]))
for elementkey in tableauTrie:
    texte += elementkey + " : " + str(tableauTrie.get(elementkey)) + "\n"
print("La liste des éléments est : \n" + texte)

# print("Le nombre d'élément est : " + compteur)





