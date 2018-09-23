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

propsTableauCritere = ['pretty_formula', 'elasticity.poisson_ratio', 'elasticity.G_Reuss', 'elasticity.G_Voigt',
                       'elasticity.G_Voigt_Reuss_Hill',
                       'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']


propsPlot = ['elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill', 'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']

critere = {"nelements": {'$lte': 6}, 'elements': {'$all': composes}, "elasticity": {'$ne': None},
           "elasticity.G_Reuss": {'$gte': 0}, "elasticity.G_Voigt": {'$gte': 0},
           "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0}, "elasticity.K_Reuss": {'$gte': 0},
           "elasticity.K_Voigt": {'$gte': 0}, "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0}}

materials = api.query(criteria=critere, properties=propsTableauCritere)

lin = len(propsPlot)
col = len(materials)

def recup(materials):
    j = 0
    tableau = np.zeros(shape=(lin, col))
    #elements = []

    for material in materials:

        #elements.append(material.get('pretty_formula'))
        i = 0
        for prop in propsPlot:
            tableau[i, j] = material.get(prop)
            i = i + 1
        j = j + 1
    return tableau



def logTableau(tableauSource):
    logTab = np.zeros(shape=(lin, col))
    for j in range(0, col):
        for i in range(0,lin):
            logTab[i,j] = math.log10(tableauSource[i,j])
    return logTab


def drawTable(tableau1, propsTableau, pdffile):
    pdf = matplotlib.backends.backend_pdf.PdfPages(pdffile)
    for prop1 in propsTableau:
        for prop2 in propsTableau:
            if prop1 != prop2:
                x = tableau1[propsTableau.index(prop1), :]
                y = tableau1[propsTableau.index(prop2), :]
                #area = 5  # 0 to 15 point radii
                plt.scatter(x, y)
                #plt.xlim(minX, maxX * 1.1)
                #plt.ylim(minY, maxY * 1.1)
                #Ne pas afficher le mot "elasticity. (suppression 11 caractères)
                plt.xlabel(prop1[11:])
                plt.ylabel(prop2[11:])
                plt.title(str(prop2[11:]) + ' versus ' + str(prop1[11:]))
                #plt.colorbar()
                pdf.savefig()
                plt.close()
    pdf.close()

resultat = recup(materials)


resultat2 = logTableau(resultat)

drawTable(resultat2, propsPlot, "LOGsiwar.pdf")






