import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.backends.backend_pdf
import matplotlib.cm as cm
import matplotlib.colors as color
import math

propsTableau = ["minLC", "maxLC", "minNu", "maxNu", "Emin", "Emax", "Gmin", "Gmax",
                'elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill',
                'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill', 'elasticity.poisson_ratio']

propsPlotLabel = [u'$LC_{min} (GPa)$', u'$LC_{max}(GPa)$', u'$\mu_{min}(GPa)$', u'$\mu_{max}(GPa)$', u'$E_{min}(GPa)$',
                  '$E_{max}(GPa)$', u'$G_{min}(GPa)$', '$G_{max}(GPa)$',
                  u'$G_{Reuss} (GPa)$', u'$G_{Voigt}(GPa)$', u'$G_{Voigt\u2000Reuss\u2000Hill}(GPa)$',
                  u'$K_{Reuss}(GPa)$', '$K_{Voigt}(GPa)$', u'$K_{Voigt\u2000Reuss\u2000Hill}(GPa)$',
                  "Poisson's ratio"]


# propsTableau = ['elasticity.poisson_ratio', 'elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill', 'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']
# propsPlotLabel = [u'$Poisson\u2000ratio$', u'$G_{Reuss} (GPa)$', u'$G_{Voigt}(GPa)$',u'$G_{Voigt\u2000Reuss\u2000Hill}(GPa)$', u'$K_{Reuss}(GPa)$', '$K_{Voigt}(GPa)$',u'$K_{Voigt\u2000Reuss\u2000Hill}(GPa)$']


# pdf = matplotlib.backends.backend_pdf.PdfPages("elastic_property_from_MP_DB_HIST_HYP4187.pdf")


def importer(fichier):
    return pd.read_csv(fichier)


data = importer("elasticElate_ALL_revisionArt_without_Zero.csv")
data.head()


def drawTable(propsTableauToPlot, pdffile, prop_graph, prop_scale):
    pdf = matplotlib.backends.backend_pdf.PdfPages(pdffile)
    # dataToPlot = []
    # tableauLabel=[]
    # couleur=[]
    minY = 1000
    maxY = -1000
    for prop in propsTableauToPlot:
        dataToPlot = data[propsTableauToPlot].get_values()
        minY = min(minY, (data[propsTableauToPlot].get_values()).min())
        # maxY=10
        maxY = max(maxY, (data[propsTableauToPlot].get_values()).max())
        tableauLabel = propsPlotLabel[propsTableau.index(prop)]
        # couleur = cm((1 + propsTableauToPlot.index(prop)) / (len(propsTableauToPlot) + 1))
    # http://www.python-simple.com/python-matplotlib/histogram.php
    nbIntervalle = 50  # 45 avec ratio poisson ; 50 pour K et G
    pas = (maxY - minY) / nbIntervalle
    bins = []
    for i in range(0, nbIntervalle):
        bins.append(minY + i * pas)
    bins.append(maxY)
    if prop_graph == "G":
        if prop_scale == "log":
            plt.hist(dataToPlot, bins=np.logspace(np.log10(1), np.log10(1000), 50), color="green", edgecolor="black",
                     lw=1,
                     label=tableauLabel, histtype='bar')  # bar est le defaut
            plt.gca().set_xscale("log")
            plt.axvspan(10, 100, facecolor='r', alpha=0.5)
        else:
            plt.hist(dataToPlot, bins=bins, color="green", edgecolor="black", lw=1, label=tableauLabel,
                     histtype='bar')
    elif prop_graph == "P":
        plt.hist(dataToPlot, bins=bins, color="green", edgecolor="black", lw=1, label=tableauLabel,
                 histtype='bar')  # bar est le
    elif prop_graph == "K":
        if prop_scale == "log":
            plt.hist(dataToPlot, bins=np.logspace(np.log10(1), np.log10(1000), 50), color="green", edgecolor="black",
                     lw=1,
                     label=tableauLabel, histtype='bar')  # bar est le defaut
            plt.axvspan(1, 10, facecolor='b', alpha=0.1)
            plt.axvspan(10, 100, facecolor='r', alpha=0.5)
            plt.axvspan(100, 300, facecolor='g', alpha=0.5)
            plt.axvspan(300, 1000, facecolor='b', alpha=0.1)
            plt.gca().set_xscale("log")
        else:
            plt.hist(dataToPlot, bins=bins, color="green", edgecolor="black", lw=1, label=tableauLabel,
                     histtype='bar')

    else:
        pass

    # plt.ylim(minY, maxY)
    plt.ylabel('Number of structures')
    # plt.xlabel('propriete')
    # plt.title('Histogramme')
    plt.legend()
    pdf.savefig()
    plt.close()
    pdf.close()


# Execution des fonction
# recuperation du tableau contenant les valeurs correspondantes au différents matériaux
# calcul et recuperation des logs du tableau selon les proprietes K ou G avec determination des Max et Min (pour determiner min et max des echelles)


cm = cm.get_cmap('gist_rainbow')
propsToPlot = ['elasticity.G_Voigt_Reuss_Hill']
drawTable(propsToPlot, "histogrammeGVRH_LOG_ALL.pdf", "G", "log")
drawTable(propsToPlot, "histogrammeGVRH_ALL.pdf", "G", "")
propsToPlot3 = ['elasticity.G_Reuss']
drawTable(propsToPlot3, "histogrammeGReuss_LOG_ALL.pdf", "G", "log")
drawTable(propsToPlot3, "histogrammeGReuss_ALL.pdf", "G", "")
propsToPlot4 = ['elasticity.G_Voigt']
drawTable(propsToPlot4, "histogrammeGVoigt_LOG_ALL.pdf", "G", "log")
drawTable(propsToPlot4, "histogrammeGVoigt_ALL.pdf", "G", "")
propsToPlot5 = ['elasticity.K_Reuss']
drawTable(propsToPlot5, "histogrammeKReuss_LOG_ALL.pdf", "K", "log")
drawTable(propsToPlot5, "histogrammeKReuss_ALL.pdf", "K", "")
propsToPlot6 = ['elasticity.K_Voigt']
drawTable(propsToPlot6, "histogrammeKVoigt_LOG_ALL.pdf", "K", "log")
drawTable(propsToPlot6, "histogrammeKVoigt_ALL.pdf", "K", "")
propsToPlot7 = ['elasticity.K_Voigt_Reuss_Hill']
drawTable(propsToPlot7, "histogrammeKVRH_LOG_ALL.pdf", "K", "log")
drawTable(propsToPlot7, "histogrammeKVRH_ALL.pdf", "K", "")
# propsToPlot2 = ['elasticity.poisson_ratio']
# drawTable(propsToPlot2, "histogrammeRatio_ALL.pdf", "P", "")
