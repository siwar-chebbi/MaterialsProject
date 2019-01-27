from pymatgen import MPRester
from com.project.elate import elastic
import numpy as np
import pandas as pd

#######################################NE PAS PRENDRE LES VALEURS -1000 ET -1500 DES MATERIAUX NON CONFORMES
api = MPRester("eDCEK5m9WVjmajp7e8af")

propsTableau = ['material_id', 'pretty_formula', "elasticity.elastic_tensor", "elasticity.G_Voigt_Reuss_Hill",
                "elasticity.K_Voigt_Reuss_Hill", "elasticity.poisson_ratio"]
composes = ['S', 'O']
# critere1: tous les elements elastiques contenant les composes S,O
critere1 = {"nelements": {'$lte': 6}, 'elements': {'$all': composes}, "elasticity": {'$ne': None},
            "elasticity.G_Reuss": {'$gte': 0}, "elasticity.G_Voigt": {'$gte': 0},
            "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0}, "elasticity.K_Reuss": {'$gte': 0},
            "elasticity.K_Voigt": {'$gte': 0}, "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0}}

# critere2: tous les elements elastiques
critere2 = {"nelements": {'$lte': 6}, "elasticity": {'$ne': None}, "elasticity.poisson_ratio": {'$gt': 0},
            "elasticity.G_Reuss": {'$gte': 0, '$lte': 1000},
            "elasticity.G_Voigt": {'$gte': 0, '$lte': 1000}, "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0, '$lte': 1000},
            "elasticity.K_Reuss": {'$gte': 0, '$lte': 1000}, "elasticity.K_Voigt": {'$gte': 0, '$lte': 1000},
            "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0, '$lte': 1000}}

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


def calculMinNu(elas):
    return elastic.minimize(elas.Poisson, 3)


def calculMaxNu(elas):
    return elastic.maximize(elas.Poisson, 3)


def calculMinYoung(elas):
    return elastic.minimize(elas.Young, 2)


def calculMaxYoung(elas):
    return elastic.maximize(elas.Young, 2)


def calculMinG(elas):
    return elastic.minimize(elas.shear, 3)


def calculMaxG(elas):
    return elastic.maximize(elas.shear, 3)


propsDisplay = ["minLC", "maxLC", "minNu", "maxNu", "K_Voigt_Reuss_Hill", "Emin", "Emax", "Gmin", "Gmax","elasticity.poisson_ratio"]

col = len(propsDisplay)
lin = len(materials)
newlin = lin
elements = []
materialIds = []
materialNonConformeEigenvalNegative = []
materialNonConformeMatSinguliere = []


def recup(materials):
    i = 0
    tableau = np.zeros(shape=(lin, col))
    for material in materials:
        print(str(i + 1), "-", str(material.get('material_id')), ":", material.get('pretty_formula'))

        matrix = material.get('elasticity.elastic_tensor')
        elastElement = generateElas(matrix)
        if elastElement:
            eigenval = sorted(np.linalg.eig(elastElement.CVoigt)[0])
            if eigenval[0] > 0:
                elements.append(material.get('pretty_formula'))
                materialIds.append(material.get('material_id'))
                tableau[i, 0] = calculMinLC(elastElement)[1]
                tableau[i, 1] = calculMaxLC(elastElement)[1]
                tableau[i, 2] = calculMinNu(elastElement)[1]
                tableau[i, 3] = calculMaxNu(elastElement)[1]
                tableau[i, 4] = material.get("elasticity.K_Voigt_Reuss_Hill")
                tableau[i, 5] = calculMinYoung(elastElement)[1]
                tableau[i, 6] = calculMaxYoung(elastElement)[1]
                tableau[i, 7] = calculMinG(elastElement)[1]
                tableau[i, 8] = calculMaxG(elastElement)[1]
                tableau[i, 9] = material.get("elasticity.poisson_ratio")

                i = i + 1
            else:
                materialNonConformeEigenvalNegative.append(material.get('material_id'))
        else:
            materialNonConformeMatSinguliere.append(material.get('material_id'))

    nb_ligne_supp = lin - len(materialIds)
    if nb_ligne_supp > 0:
        tableau = tableau[:-nb_ligne_supp, :]
    return tableau


def export(donnees, ligne, nomColonnes, fichier):
    my_df = pd.DataFrame(donnees)
    my_df.index = ligne
    my_df.to_csv(fichier, index=ligne, header=nomColonnes)


resultat = recup(materials)

export(resultat, materialIds, propsDisplay, "elasticRatioNega.csv")

print("materials non conformes, eigenVal negative:\n" + str(materialNonConformeEigenvalNegative))
print("materials non conformes, matrice singuliere:\n" + str(materialNonConformeMatSinguliere))
print("nbre de materials non conformes, eigenVal negative:\n" + str(materialNonConformeEigenvalNegative.__len__()))
print("nbre de materials non conformes, matrice singuliere:\n" + str(materialNonConformeMatSinguliere.__len__()))
