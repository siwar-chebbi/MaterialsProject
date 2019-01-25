import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_squared_error,  r2_score
import matplotlib.backends.backend_pdf
import matplotlib.cm as cm
import matplotlib.colors as color
import math


propsTableau = ['elasticity.poisson_ratio', 'elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill',
                'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']
propsPlotLabel = [u'$Poisson\u2000ratio$', u'$G_{Reuss} (GPa)$', u'$G_{Voigt}(GPa)$', u'$G_{Voigt\u2000Reuss\u2000Hill}(GPa)$', u'$K_{Reuss}(GPa)$', '$K_{Voigt}(GPa)$', u'$K_{Voigt\u2000Reuss\u2000Hill}(GPa)$']



#pdf = matplotlib.backends.backend_pdf.PdfPages("elastic_property_from_MP_DB_HIST_HYP4187.pdf")



def importer (fichier):
    return pd.read_csv(fichier)

data=importer("elastic_property_from_MP_DB_12522.csv")
data.head()


def drawTable(propsTableauToPlot, pdffile):
    pdf = matplotlib.backends.backend_pdf.PdfPages(pdffile)
   # dataToPlot = []
   # tableauLabel=[]
  #couleur=[]
    minY = 1000
    maxY = -1000
    for prop in propsTableauToPlot:
        dataToPlot = data[propsTableauToPlot].get_values()
        minY = min(minY, (data[propsTableauToPlot].get_values()).min())
        #maxY=10
        maxY = max(maxY, (data[propsTableauToPlot].get_values()).max())
        tableauLabel = propsPlotLabel[propsTableau.index(prop)]
        couleur = cm((1+propsTableauToPlot.index(prop))/(len(propsTableauToPlot)+1))
    #http://www.python-simple.com/python-matplotlib/histogram.php
    nbIntervalle=50
    pas = (maxY-minY)/nbIntervalle
    bins = []
    for i in range(0, nbIntervalle):
        bins.append(minY + i*pas)
    bins.append(maxY)

    plt.hist(dataToPlot, bins=bins, color="green", edgecolor="black", lw=1, label=tableauLabel, histtype='bar')  # bar est le defaut
#plt.ylim(minY, maxY)
    plt.ylabel('Nombre of structures')
    #plt.xlabel('propriete')
    #plt.title('Histogramme')
    plt.legend()
    pdf.savefig()
    plt.close()
    pdf.close()


#Execution des fonction
#recuperation du tableau contenant les valeurs correspondantes au différents matériaux
#calcul et recuperation des logs du tableau selon les proprietes K ou G avec determination des Max et Min (pour determiner min et max des echelles)


cm = cm.get_cmap('gist_rainbow')
propsToPlot = ['elasticity.G_Voigt_Reuss_Hill']
drawTable(propsToPlot, "histogrammeGVRH.pdf")

propsToPlot2 = ['elasticity.poisson_ratio']
drawTable(propsToPlot2, "histogrammeRatioGVRH.pdf")

propsToPlot3 = ['elasticity.G_Reuss']
drawTable(propsToPlot3, "histogrammeGReuss.pdf")


propsToPlot4 = ['elasticity.G_Voigt']
drawTable(propsToPlot4, "histogrammeGVoigt.pdf")


propsToPlot5 = ['elasticity.K_Reuss']
drawTable(propsToPlot5, "histogrammeKReuss.pdf")


propsToPlot6 = ['elasticity.K_Voigt']
drawTable(propsToPlot6, "histogrammeKVoigt.pdf")


propsToPlot7 = ['elasticity.K_Voigt_Reuss_Hill']
drawTable(propsToPlot7, "histogrammeKVRH.pdf")
