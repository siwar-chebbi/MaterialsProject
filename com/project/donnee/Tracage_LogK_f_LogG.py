#!/opt/anaconda3/bin/python
from typing import List
from pandas._libs import properties
from pymatgen import MPRester
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as color
import matplotlib.backends.backend_pdf
import math


import numpy as np

api = MPRester("fB610TDF3LSwxiN9")

composes = ['S', 'O']

propsTableau = ['elasticity.poisson_ratio', 'elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill',
                'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']

propsTableauG = ['elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill']
propsTableauK = ['elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']

critere1 = {"nelements": {'$lte': 6}, 'elements': {'$all': composes}, "elasticity": {'$ne': None},
            "elasticity.G_Reuss": {'$gte': 0}, "elasticity.G_Voigt": {'$gte': 0},
            "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0}, "elasticity.K_Reuss": {'$gte': 0},
            "elasticity.K_Voigt": {'$gte': 0}, "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0}}
materials = api.query(criteria=critere1, properties=propsTableau)

lin = len(propsTableau)
col = len(materials)

def recup(materials):
    j = 0
    tableau = np.zeros(shape=(lin, col))
    elements = []

    for material in materials:

        elements.append(material.get('pretty_formula'))
        i = 0
        for prop in propsTableau:
            tableau[i, j] = material.get(prop)
            i = i + 1
        j = j + 1
    return tableau


def logTableau(tableauSource, propsTabExtraite):
    linExt = len(propsTabExtraite)
    i = 0
    tableauLog = np.zeros(shape=(linExt, col))
    for prop in propsTableau:
        index = propsTableau.index(prop)
        if prop in propsTabExtraite:
            for j in range(0, col):
                tableauLog[i, j] = math.log10(tableauSource[index, j])
            i = i + 1
    return tableauLog


def drawTable(tableau1, propsTableau1, minX, maxX, tableau2, propsTableau2, minY, maxY, pdffile):
    pdf = matplotlib.backends.backend_pdf.PdfPages(pdffile)
    for prop1 in propsTableau1:
        for prop2 in propsTableau2:
            if prop1 != prop2:
                x = tableau1[propsTableau1.index(prop1), :]
                y = tableau2[propsTableau2.index(prop2), :]
                area = 5  # 0 to 15 point radii
                plt.scatter(x, y, s=area, c=poisson, cmap=cm.get_cmap('seismic'), norm=normalize, alpha=1)
                plt.xlim(minX * 1.1, maxX * 1.1)
                plt.ylim(minY * 1.1, maxY * 1.1)
                plt.xlabel(prop1)
                plt.ylabel(prop2)
                plt.title(str(prop2) + ' versus ' + str(prop1))
                plt.colorbar()
                pdf.savefig()
                plt.close()
    pdf.close()

#Execution des fonction
#recuperation du tableau contenant les valeurs correspondantes au différents matériaux
resultat = recup(materials)
#calcul et recuperation des logs du tableau selon les proprietes K ou G avec determination des Max et Min (pour determiner min et max des echelles)

logTableauK = logTableau(resultat, propsTableauK)
maxK = logTableauK.max()
minK = logTableauK.min()

logTableauG = logTableau(resultat, propsTableauG)
maxG = logTableauG.max()
minG = logTableauG.min()

poisson = resultat[propsTableau.index('elasticity.poisson_ratio'), :]

normalize = color.Normalize(vmin=min(poisson), vmax=max(poisson))

drawTable(logTableauK, propsTableauK,  minK, maxK, logTableauG, propsTableauG, minG, maxG, "G_f_K.pdf")
drawTable(logTableauG, propsTableauG, minG, maxG, logTableauK, propsTableauK, minK, maxK, "K_f_G.pdf")


