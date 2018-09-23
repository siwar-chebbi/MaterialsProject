from pymatgen import MPRester
from project.elate import elastic
import numpy as np
import pandas as pd

api = MPRester("fB610TDF3LSwxiN9")

propsTableau = ['material_id','pretty_formula',"elasticity.elastic_tensor"]
composes = ['S', 'O']
critere = {"nelements": {'$lte': 6}, 'elements': {'$all': composes}, "elasticity": {'$ne': None},
           "elasticity.G_Reuss": {'$gte': 0}, "elasticity.G_Voigt": {'$gte': 0},
           "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0}, "elasticity.K_Reuss": {'$gte': 0},
           "elasticity.K_Voigt": {'$gte': 0}, "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0}}

materials = api.query(criteria=critere, properties=propsTableau)


# test= elastic.ELATE_MaterialsProject("mp-2133")
# matrix = elastic_tensor

def generateElas(matrix):
    try:
        elas = elastic.Elastic(matrix)
    except ValueError as e:
        print('Error: Invalid stiffness matrix: ')
        print(e.args[0])

    if elas.isOrthorhombic():
       elas = elastic.ElasticOrtho(elas)
    return elas


def calculMinLC(elas):
 return elastic.minimize(elas.LC, 2)

def calculMaxLC(elas):
 return elastic.maximize(elas.LC, 2)


propsDisplay =["minLC", "maxLC"]

lin = len(propsDisplay)
col = len(materials)
elements = []
materialIds = []

def recup(materials):
    j = 0
    tableau = np.zeros(shape=(lin, col))

    for material in materials:
        elements.append(material.get('pretty_formula'))
        materialIds.append(material.get('material_id'))
        matrix = material.get('elasticity.elastic_tensor')
        elastElement= generateElas(matrix)
        minLC = calculMinLC(elastElement)[1]
        maxLC = calculMaxLC(elastElement)[1]

        #i = 0
        # for prop in propsDisplay:
        #     tableau[i, j] = globals()[prop]
        #     i = i + 1
        tableau[0, j] = minLC
        tableau[1, j] = maxLC
        j = j + 1
    return tableau

def export (donnees,ligne,nomColonnes,fichier):
  my_df = pd.DataFrame(donnees)
  my_df.to_csv(fichier, index=ligne, header=nomColonnes)


resultat = recup(materials)

export(resultat, propsDisplay, materialIds, "elastic.csv")

