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
                  u'$K_{Reuss}(GPa)$', '$K_{Voigt}(GPa)$', u'$K_{Voigt\u2000Reuss\u2000Hill}(GPa)$', u'$Poisson\u2000ratio$']



#propsTableau = ["minLC", "maxLC", "minNu", "maxNu", "K_Voigt_Reuss_Hill", "Emin", "Emax", "Gmin", "Gmax"]

#propsPlotLabel = [u'$LC_{min} (GPa)$', u'$LC_{max}(GPa)$', u'$\mu_{min}(GPa)$', u'$\mu_{max}(GPa)$', u'$K_{Voigt\u2000Reuss\u2000Hill}(GPa)$', u'$E_{min}(GPa)$', '$E_{max}(GPa)$', u'$G_{min}(GPa)$', '$G_{max}(GPa)$']


def importer(fichier):
    return pd.read_csv(fichier)


data = importer("elasticElate_ALL_revisionArt_without_Zero.csv")
data.head()

Emax_list = data['Emax'].get_values()
Emin_list = data['Emin'].get_values()
Emax_sur_Emin = []

for x, y in zip(Emax_list, Emin_list):
    if y == 0:
        continue
    else:
        Emax_sur_Emin.append(x / y)

print(len(Emax_sur_Emin))
print(min(Emax_sur_Emin))
print(max(Emax_sur_Emin))

Emax_sur_EminLOG = [math.log10(i) for i in Emax_sur_Emin]
print(len(Emax_sur_EminLOG))
print(min(Emax_sur_EminLOG))
print(max(Emax_sur_EminLOG))

def drawTable(propsTableauToPlot, pdffile):
    pdf = matplotlib.backends.backend_pdf.PdfPages(pdffile)

    dataToPlot = Emax_sur_Emin
    #minY = min(Emax_sur_Emin)

    #maxY = max(Emax_sur_Emin)
    tableauLabel = propsTableauToPlot
    #couleur = "green"

    # http://www.python-simple.com/python-matplotlib/histogram.php
    maxY= 100
    minY= 1

    nbIntervalle = 50
    pas = (maxY - minY) / nbIntervalle
    bins = []
    for i in range(0, nbIntervalle):
        bins.append(minY + i * pas)
    bins.append(maxY)

    #plt.hist(dataToPlot, bins=np.linspace(1, 100, 100), color="green", edgecolor="black", lw=1, label=tableauLabel, histtype='bar') # bar est le defaut

    plt.hist(dataToPlot, bins=np.logspace(np.log10(1),np.log10(300), 50), color="green", edgecolor="black", lw=1, label=tableauLabel,
             histtype='bar')  # bar est le defaut

    # plt.ylim(minY, maxY)
    plt.ylabel('Number of structures')
    plt.gca().set_xscale("log")
    # plt.xlabel('propriete')
    # plt.title('Histogramme')
    plt.legend()
    pdf.savefig()
    plt.close()
    pdf.close()


cm = cm.get_cmap('gist_rainbow')
propsToPlot = ['Elastic anisotropy']
drawTable(propsToPlot, "Emax_sur_EminLOG.pdf")
