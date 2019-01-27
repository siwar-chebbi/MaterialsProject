import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.backends.backend_pdf
import matplotlib.cm as cm
import matplotlib.colors as color
import math

propsTableau = ["minLC", "maxLC", "minNu", "maxNu", "K_Voigt_Reuss_Hill", "Emin", "Emax", "Gmin", "Gmax"]

propsPlotLabel = [u'$LC_{min} (GPa)$', u'$LC_{max}(GPa)$', u'$\mu_{min}(GPa)$', u'$\mu_{max}(GPa)$',
                  u'$K_{Voigt\u2000Reuss\u2000Hill}(GPa)$', u'$E_{min}(GPa)$', '$E_{max}(GPa)$', u'$G_{min}(GPa)$',
                  '$G_{max}(GPa)$']


def importer(fichier):
    return pd.read_csv(fichier)


data = importer("elasticRatioPoissonPositive.csv")
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


def drawTable(propsTableauToPlot, pdffile):
    pdf = matplotlib.backends.backend_pdf.PdfPages(pdffile)

    dataToPlot = Emax_sur_Emin
    minY = min(Emax_sur_Emin)

    maxY = max(Emax_sur_Emin)
    tableauLabel = propsTableauToPlot
    couleur = "green"

    # http://www.python-simple.com/python-matplotlib/histogram.php
    nbIntervalle = 50

    # maxY=250
    # minY=0

    pas = (maxY - minY) / nbIntervalle
    bins = []
    for i in range(0, nbIntervalle):
        bins.append(minY + i * pas)
    bins.append(maxY)

    plt.hist(dataToPlot, bins=bins, color="green", edgecolor="black", lw=1, label=tableauLabel,
             histtype='bar')  # bar est le defaut
    # plt.hist(dataToPlot, bins=100, color='blue', edgecolor='black',lw=1,histtype='bar')
    # plt.ylim(minY, maxY)
    plt.ylabel('Nombre of elements')
    # plt.xlabel('propriete')
    # plt.title('Histogramme')
    plt.legend()
    pdf.savefig()
    plt.close()
    pdf.close()


cm = cm.get_cmap('gist_rainbow')
propsToPlot = ['Emax_sur_Emin']
drawTable(propsToPlot, "Emax_sur_Emin.pdf")
