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

propsTableau = ['elasticity.poisson_ratio', 'elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill',
                'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']


critere = {"nelements": {'$lte': 6}, 'elements': {'$all': composes}, "elasticity": {'$ne': None},
           "elasticity.G_Reuss": {'$gte': 0}, "elasticity.G_Voigt": {'$gte': 0},
           "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0}, "elasticity.K_Reuss": {'$gte': 0},
           "elasticity.K_Voigt": {'$gte': 0}, "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0}}

materials = api.query(criteria=critere, properties=propsTableauCritere)

lin = len(propsTableau)
col = len(materials)
elements = []
def recup(materials):
    j = 0
    tableau = np.zeros(shape=(lin, col))

    for material in materials:

        elements.append(material.get('pretty_formula'))
        i = 0
        for prop in propsTableau:
            tableau[i, j] = material.get(prop)
            i = i + 1
        j = j + 1
    return tableau


def drawTable(tableauSource, propsTableauToPlot, pdffile):
    pdf = matplotlib.backends.backend_pdf.PdfPages(pdffile)
    minY = 1000
    maxY = -1000

    for prop in propsTableauToPlot:
        couleur = cm((1 + propsTableauToPlot.index(prop)) / (len(propsTableauToPlot) + 1))
        minY = min(minY, (tableauSource[propsTableau.index(prop), :]).min())
        epaisseur = 0.2
        abscisses = []
        for i in range(0, col):
            abscisses.append(i+1 + propsTableauToPlot.index(prop)*epaisseur-int(len(propsTableauToPlot)*epaisseur/2))

        maxY = max(maxY, (tableauSource[propsTableau.index(prop), :]).max())
        plt.bar(abscisses, tableauSource[propsTableau.index(prop), :], color=couleur, width=epaisseur, label=prop[11:], align='center')  # bar est le defaut
    #maxY = 3
    plt.ylim(minY, maxY*1.01)
    plt.ylabel('valeurs')
    plt.xlabel('elements')
    plt.xticks(range(col), elements, rotation='vertical', fontsize=7)
    plt.title('Graphiques Elements')
    plt.legend()
    pdf.savefig()
    plt.close()
    pdf.close()


#Execution des fonction
#recuperation du tableau contenant les valeurs correspondantes au différents matériaux
resultat = recup(materials)
#calcul et recuperation des logs du tableau selon les proprietes K ou G avec determination des Max et Min (pour determiner min et max des echelles)



cm = cm.get_cmap('gist_rainbow')
#propsToPlot = ['elasticity.G_Reuss','elasticity.G_Voigt','elasticity.G_Voigt_Reuss_Hill','elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']
propsToPlot = ['elasticity.G_Reuss']

drawTable(resultat, propsToPlot, "bar.pdf")



