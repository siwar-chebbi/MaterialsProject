from pymatgen import MPRester
from project.elate import elastic
import numpy as np
import pandas as pd

api = MPRester("fB610TDF3LSwxiN9")

propsTableau = ['material_id','pretty_formula',"elasticity.elastic_tensor"]
composes = ['S', 'O']
#critere1: tous les elements elastiques contenant les composes S,O
critere1 = {"nelements": {'$lte': 6}, 'elements': {'$all': composes}, "elasticity": {'$ne': None},
           "elasticity.G_Reuss": {'$gte': 0}, "elasticity.G_Voigt": {'$gte': 0},
           "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0}, "elasticity.K_Reuss": {'$gte': 0},
           "elasticity.K_Voigt": {'$gte': 0}, "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0}}

#critere2: tous les elements elastiques
critere2 = {"nelements": {'$lte': 6}, "elasticity": {'$ne': None}, "elasticity.G_Reuss": {'$gte': 0},
            "elasticity.G_Voigt": {'$gte': 0}, "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0},
            "elasticity.K_Reuss": {'$gte': 0}, "elasticity.K_Voigt": {'$gte': 0},
            "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0}}


materials = api.query(criteria=critere2, properties=propsTableau)


# test= elastic.ELATE_MaterialsProject("mp-2133")
# matrix = elastic_tensor

def generateElas(matrix):
    try:
        elas = elastic.Elastic(matrix)
    except ValueError as e:
        print('Error: Invalid stiffness matrix: ')
        print(e.args[0])
        return None

    if elas.isOrthorhombic():
        elas = elastic.ElasticOrtho(elas)
    return elas


def calculMinLC(elas):
 return elastic.minimize(elas.LC, 2)

def calculMaxLC(elas):
 return elastic.maximize(elas.LC, 2)


propsDisplay =["minLC", "maxLC"]

col = len(propsDisplay)
lin = len(materials)
elements = []
materialIds =[]
materialNonConforme=[]

def recup(materials):
    i = 0
    tableau = np.zeros(shape=(lin, col))
    for material in materials:
        elements.append(material.get('pretty_formula'))
        materialIds.append(material.get('material_id'))
        matrix = material.get('elasticity.elastic_tensor')
        elastElement= generateElas(matrix)
        if elastElement:
             minLC = calculMinLC(elastElement)[1]
             maxLC = calculMaxLC(elastElement)[1]
        else:
            materialNonConforme.append(material.get('material_id'))
            minLC = -1000
            maxLC = -1000

        #i = 0
        # for prop in propsDisplay:
        #     tableau[i, j] = globals()[prop]
        #     i = i + 1
        tableau[i, 0] = minLC
        tableau[i, 1] = maxLC
        i = i + 1
    return tableau

def export (donnees,ligne,nomColonnes,fichier):
  my_df = pd.DataFrame(donnees)
  my_df.index = ligne
  my_df.to_csv(fichier, index=ligne, header=nomColonnes)


resultat = recup(materials)

export(resultat, materialIds, propsDisplay, "elastic.csv")
print("materials non conformes:\n" + str(materialNonConforme))

