from pymatgen import MPRester
import pandas as pd
import numpy as np
from pymatgen.core.periodic_table import Element

api = MPRester("eDCEK5m9WVjmajp7e8af")

propsTableauCritere = ['material_id', 'pretty_formula', 'energy', 'energy_per_atom', 'density',
                       'formation_energy_per_atom', 'volume']

propsTableau = ['elasticity.poisson_ratio', ]

all_elements = dict()


# Recuperer tous les elements du tableau périodique
def get_all_elements():
    for element in Element:
        try:
            all_elements[element.value] = \
                api._make_request("/element/%s/tasks/isolated_atom" % element.value, mp_decode=False)[0]
        except:
            pass
            print("Problème de recuperation de l'élément atomique  : " + str(element))


# Import du fichier csv initlial
def importer(fichier):
    return pd.read_csv(fichier, index_col=0)


# material Ids du tableau a remlir dans get_calculated_properties
tableau_material_ids = []


def get_calculated_properties(entries_var):
    i = 0
    tableau = np.zeros(shape=(entries_var.__len__(), 2))
    for entry in entries_var:
        tableau_material_ids.append(entry.data['material_id'])
        try:
            tableau[i, 0] = get_cohesive_energy(entry, False)
        except:
            pass
            print("Erreur calcul get_cohesive_energy de  : " + str(entry.data['material_id']))
            tableau[i, 0] = -1

        try:
            tableau[i, 1] = entry.composition.average_electroneg
        except:
            pass
            print("Erreur calcul average_electroneg de  : " + str(entry.data['material_id']))
            tableau[i, 0] = -2
        i = i + 1

    return pd.DataFrame(tableau, index=tableau_material_ids, columns=['cohesive_energy', 'average_electroneg'])


def export_additional_properties(data1, data2, fichier):
    mergedDf = data1.merge(data2, left_index=True, right_index=True)
    mergedDf.to_csv(fichier, index=True)


def get_cohesive_energy(entry, per_atom):
    ebulk = entry.energy / entry.composition.get_integer_formula_and_factor()[1]
    comp_dict = entry.composition.reduced_composition.as_dict()
    isolated_atom_e_sum, n = 0, 0
    for atom in comp_dict.keys():
        e = all_elements.get(atom)
        isolated_atom_e_sum += e['output']["final_energy_per_atom"] * comp_dict[atom]
        n += comp_dict[atom]
    ecoh_per_formula = isolated_atom_e_sum - ebulk
    return ecoh_per_formula / n if per_atom else ecoh_per_formula


# 1- On recupere toutes les lignes du fichiers csv
data_from_cv = importer("elasticElate_ALL_revisionArt_without_Zero.csv")

print("\nNombre de tous les éléments dans le fichier csv = {}\n".format(data_from_cv.shape[0]))

# 2- Lister tous les material Ids
Materials_Ids = tuple(list(data_from_cv.index.values))

# 3- Recuperation de tous les entries avec les propriétés correspondantes
entries = api.get_entries({'material_id': {'$in': Materials_Ids}}, property_data=propsTableauCritere)
# materials = api.query(criteria={'material_id': {'$in': Materials_Ids}}, properties=propsTableauCritere)

# 4- Recupération des propriétés des éléments du tableau periodique
get_all_elements()

# 5-Calcul des propriétés à completer au fichier csv
data_additional_prop = get_calculated_properties(entries)

# 6- Generation du nouveau fichier csv contenant toutes les proprietes
export_additional_properties(data_from_cv, data_additional_prop, 'test.csv')
