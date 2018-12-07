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

api = MPRester("eDCEK5m9WVjmajp7e8af")

composes = ['S', 'O']

propsTableauCritere = ['pretty_formula', 'elasticity.poisson_ratio', 'elasticity.G_Reuss', 'elasticity.G_Voigt',
                       'elasticity.G_Voigt_Reuss_Hill',
                       'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']

propsTableau = ['elasticity.poisson_ratio', 'elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill',
                'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']
propsPlotLabel = [u'$Poisson\u2000ratio$', u'$G_{Reuss} (GPa)$', u'$G_{Voigt}(GPa)$', u'$G_{Voigt\u2000Reuss\u2000Hill}(GPa)$', u'$K_{Reuss}(GPa)$', '$K_{Voigt}(GPa)$', u'$K_{Voigt\u2000Reuss\u2000Hill}(GPa)$']

critere4 = {"nelements": {'$lte': 6}, "elasticity": {'$ne': None}, "elasticity.G_Reuss": {'$gte': 0, '$lte': 1000},
                "elasticity.G_Voigt": {'$gte': 0, '$lte': 1000}, "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0, '$lte': 1000},
                "elasticity.K_Reuss": {'$gte': 0, '$lte': 1000}, "elasticity.K_Voigt": {'$gte': 0, '$lte': 1000},
                "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0, '$lte': 1000}}


critere = {"nelements": {'$lte': 6}, 'elements': {'$all': composes}, "elasticity": {'$ne': None},
           "elasticity.G_Reuss": {'$gte': 0}, "elasticity.G_Voigt": {'$gte': 0},
           "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0}, "elasticity.K_Reuss": {'$gte': 0},
           "elasticity.K_Voigt": {'$gte': 0}, "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0}}

materials = api.query(criteria=critere4, properties=propsTableauCritere)

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
    dataToPlot = []
    tableauLabel=[]
    couleur=[]
    minY = 1000
    maxY = -1000
    for prop in propsTableauToPlot:
        dataToPlot.append(tableauSource[propsTableau.index(prop), :])
        minY = min(minY, (tableauSource[propsTableau.index(prop), :]).min())
        #maxY=10
        maxY = max(maxY, (tableauSource[propsTableau.index(prop), :]).max())
        tableauLabel.append(propsPlotLabel[propsTableau.index(prop)])
        couleur.append(cm((1+propsTableauToPlot.index(prop))/(len(propsTableauToPlot)+1)))
    #http://www.python-simple.com/python-matplotlib/histogram.php
    nbIntervalle=50
    pas = (maxY-minY)/nbIntervalle
    bins = []
    for i in range(0, nbIntervalle):
        bins.append(minY + i*pas)
    bins.append(maxY)

    plt.hist(dataToPlot, bins=bins, color=couleur, edgecolor="black", lw=1, label=tableauLabel, histtype='bar')  # bar est le defaut
#plt.ylim(minY, maxY)
    plt.ylabel('Nombre of elements')
    #plt.xlabel('propriete')
    #plt.title('Histogramme')
    plt.legend()
    pdf.savefig()
    plt.close()
    pdf.close()


#Execution des fonction
#recuperation du tableau contenant les valeurs correspondantes au différents matériaux
resultat = recup(materials)
#calcul et recuperation des logs du tableau selon les proprietes K ou G avec determination des Max et Min (pour determiner min et max des echelles)


cm = cm.get_cmap('gist_rainbow')
propsToPlot = ['elasticity.G_Voigt_Reuss_Hill']
drawTable(resultat, propsToPlot, "histogrammeGVRH.pdf")

propsToPlot2 = ['elasticity.poisson_ratio']
drawTable(resultat, propsToPlot2, "histogrammeRatioGVRH.pdf")

propsToPlot3 = ['elasticity.G_Reuss']
drawTable(resultat, propsToPlot3, "histogrammeGReuss.pdf")


propsToPlot4 = ['elasticity.G_Voigt']
drawTable(resultat, propsToPlot4, "histogrammeGVoigt.pdf")


propsToPlot5 = ['elasticity.K_Reuss']
drawTable(resultat, propsToPlot5, "histogrammeKReuss.pdf")


propsToPlot6 = ['elasticity.K_Voigt']
drawTable(resultat, propsToPlot6, "histogrammeKVoigt.pdf")


propsToPlot7 = ['elasticity.K_Voigt_Reuss_Hill']
drawTable(resultat, propsToPlot7, "histogrammeKVRH.pdf")