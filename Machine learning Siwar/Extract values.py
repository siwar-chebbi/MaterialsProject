from pymatgen import MPRester
from pymatgen import Composition

api = MPRester("eDCEK5m9WVjmajp7e8af")
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.backends.backend_pdf
import matplotlib.cm as cm
import matplotlib.colors as color
import math

propsTableauCritere = ['pretty_formula', 'energy', 'energy_per_atom', 'density', 'formation_energy_per_atom']
propsTableau = ['elasticity.poisson_ratio', ]


def importer(fichier):
    return pd.read_csv(fichier)


data = importer("elasticElate_ALL_revisionArt_without_Zero.csv")

print("\nNombre de tous les éléments dans le fichier csv = {}\n".format(data.shape[0]))

Materials_ID = []
for material in data.iloc[:, 0]:
    Materials_ID.append(material)

print(len(Materials_ID))





materials = api.query(criteria={'material_id': {'$in': Materials_ID}}, properties= propsTableauCritere)

def recup(materials):
    j = 0
    tableau = np.zeros(shape=(lin, col))
    # elements = []

    for material in materials:
        i = 0
        for prop in propsTableau:
            tableau[i, j] = material.get(prop)
            i = i + 1
        j = j + 1
    return tableau


def export(donnees, ligne, nomColonnes, fichier):
    my_df = pd.DataFrame(donnees)
    my_df.index = ligne
    my_df.to_csv(fichier, index=ligne, header=nomColonnes)

resultat = recup(materials)
export(resultat.transpose(), Materials_ID, propsTableau, "test.csv")')






comp = Composition('C')
print(comp.average_electroneg)
