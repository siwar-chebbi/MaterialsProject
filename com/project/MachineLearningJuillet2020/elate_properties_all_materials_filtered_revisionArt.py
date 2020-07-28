from pymatgen import MPRester
from MaterialsProject.com.project.elate import elastic
import numpy as np
import pandas as pd

#######################################NE PAS PRENDRE LES VALEURS -1000 ET -1500 DES MATERIAUX NON CONFORMES
api = MPRester("78OAi0lR9kdkyiAi")

propsTableau = ['material_id', 'pretty_formula', "elasticity.elastic_tensor", 'elasticity.poisson_ratio',
                'elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill',
                'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']
composes = ['S', 'O']
# critere1: tous les elements elastiques contenant les composes S,O
critere1 = {"nelements": {'$lte': 6}, 'elements': {'$all': composes}, "elasticity": {'$ne': None},
            "elasticity.G_Reuss": {'$gte': 0}, "elasticity.G_Voigt": {'$gte': 0},
            "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0}, "elasticity.K_Reuss": {'$gte': 0},
            "elasticity.K_Voigt": {'$gte': 0}, "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0}}

# critere3: tous les elements elastiques
critere3 = {"nelements": {'$lte': 6}, "elasticity": {'$ne': None},
            "elasticity.G_Reuss": {'$gte': 0, '$lte': 1000},
            "elasticity.G_Voigt": {'$gte': 0, '$lte': 1000}, "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0, '$lte': 1000},
            "elasticity.K_Reuss": {'$gte': 0, '$lte': 1000}, "elasticity.K_Voigt": {'$gte': 0, '$lte': 1000},
            "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0, '$lte': 1000}}

critere3EXP = {"nelements": {'$lte': 6}, "elasticity": {'$ne': None}, 'icsd_ids.0': {'$exists': True},
               "elasticity.G_Reuss": {'$gte': 0, '$lte': 1000},
               "elasticity.G_Voigt": {'$gte': 0, '$lte': 1000},
               "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0, '$lte': 1000},
               "elasticity.K_Reuss": {'$gte': 0, '$lte': 1000}, "elasticity.K_Voigt": {'$gte': 0, '$lte': 1000},
               "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0, '$lte': 1000}}

critere3HYP = {"nelements": {'$lte': 6}, "elasticity": {'$ne': None}, 'icsd_ids.0': {'$exists': False},
               "elasticity.G_Reuss": {'$gte': 0, '$lte': 1000},
               "elasticity.G_Voigt": {'$gte': 0, '$lte': 1000},
               "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0, '$lte': 1000},
               "elasticity.K_Reuss": {'$gte': 0, '$lte': 1000}, "elasticity.K_Voigt": {'$gte': 0, '$lte': 1000},
               "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0, '$lte': 1000}}
# critere4: tous les elements elastiques ratio de poisson negatif

critere4PoissonNeg = {"nelements": {'$lte': 6}, "elasticity": {'$ne': None}, "elasticity.poisson_ratio": {'$lt': 0},
            "elasticity.G_Reuss": {'$gte': 0, '$lte': 1000},
            "elasticity.G_Voigt": {'$gte': 0, '$lte': 1000}, "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0, '$lte': 1000},
            "elasticity.K_Reuss": {'$gte': 0, '$lte': 1000}, "elasticity.K_Voigt": {'$gte': 0, '$lte': 1000},
            "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0, '$lte': 1000}}

critere4PoissonNegEXP = {"nelements": {'$lte': 6}, "elasticity": {'$ne': None}, "elasticity.poisson_ratio": {'$lt': 0},
               'icsd_ids.0': {'$exists': True},
               "elasticity.G_Reuss": {'$gte': 0, '$lte': 1000},
               "elasticity.G_Voigt": {'$gte': 0, '$lte': 1000},
               "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0, '$lte': 1000},
               "elasticity.K_Reuss": {'$gte': 0, '$lte': 1000}, "elasticity.K_Voigt": {'$gte': 0, '$lte': 1000},
               "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0, '$lte': 1000}}


materials = api.query(criteria=critere4PoissonNeg, properties=propsTableau)


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


propsDisplay_old = ["minLC", "maxLC", "minNu", "maxNu", "K_Voigt_Reuss_Hill", "Emin", "Emax", "Gmin", "Gmax",
                    "elasticity.poisson_ratio"]
propsDisplay = ["minLC", "maxLC", "minNu", "maxNu", "Emin", "Emax", "Gmin", "Gmax",
                'elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill',
                'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill', 'elasticity.poisson_ratio']

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
                tableau[i, 4] = calculMinYoung(elastElement)[1]
                tableau[i, 5] = calculMaxYoung(elastElement)[1]
                tableau[i, 6] = calculMinG(elastElement)[1]
                tableau[i, 7] = calculMaxG(elastElement)[1]
                tableau[i, 8] = material.get('elasticity.G_Reuss')
                tableau[i, 9] = material.get('elasticity.G_Voigt')
                tableau[i, 10] = material.get('elasticity.G_Voigt_Reuss_Hill')
                tableau[i, 11] = material.get('elasticity.K_Reuss')
                tableau[i, 12] = material.get('elasticity.K_Voigt')
                tableau[i, 13] = material.get("elasticity.K_Voigt_Reuss_Hill")
                tableau[i, 14] = material.get("elasticity.poisson_ratio")
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

def importer(fichier):
    return pd.read_csv(fichier, index_col=0)

def export_equal_0(file_name):
    data = importer(file_name)
    extract_data2 = data[(data['elasticity.G_Reuss'] != 0) & (data['elasticity.G_Voigt'] != 0) &
                         (data['elasticity.G_Voigt_Reuss_Hill'] != 0) &
                         (data['elasticity.K_Reuss'] != 0) &
                         (data['elasticity.K_Voigt'] != 0) &
                         (data['elasticity.K_Voigt_Reuss_Hill'] != 0)]
    extract_data3 = data[(data['elasticity.G_Reuss'] == 0) | (data['elasticity.G_Voigt'] == 0) |
                         (data['elasticity.G_Voigt_Reuss_Hill'] == 0) |
                         (data['elasticity.K_Reuss'] == 0) |
                         (data['elasticity.K_Voigt'] == 0) |
                         (data['elasticity.K_Voigt_Reuss_Hill'] == 0)]
    extract_data2.to_csv("elasticElate_All_without_Zero_JUL2020.csv")
    extract_data3.to_csv("elasticElate_All_with_Zero_JUL2020.csv")

resultat = recup(materials)

export(resultat, materialIds, propsDisplay, "elasticElate_PoissonNeg_JUL2020.csv")

#export_equal_0("elasticElate_ALL_revisionArt_PoissonNega_EXP.csv")

print("materials non conformes, eigenVal negative:\n" + str(materialNonConformeEigenvalNegative))
print("materials non conformes, matrice singuliere:\n" + str(materialNonConformeMatSinguliere))
print("nbre de materials non conformes, eigenVal negative:\n" + str(materialNonConformeEigenvalNegative.__len__()))
print("nbre de materials non conformes, matrice singuliere:\n" + str(materialNonConformeMatSinguliere.__len__()))
