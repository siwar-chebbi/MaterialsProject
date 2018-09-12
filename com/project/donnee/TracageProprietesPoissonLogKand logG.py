#!/opt/anaconda3/bin/python
from typing import List

from pandas._libs import properties
from pymatgen import MPRester
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as color
import os
import matplotlib.backends.backend_pdf
import math
#pdf = matplotlib.backends.backend_pdf.PdfPages("C:\\Users\\siwar\\Desktop\\image\\output.pdf")
pdf = matplotlib.backends.backend_pdf.PdfPages("t.pdf")

from pymatgen.analysis.elasticity import elastic
from pymatgen.util import plotting
import numpy as np

api = MPRester("fB610TDF3LSwxiN9")

composes = ['S', 'O']


propsTableau = ['elasticity.poisson_ratio','elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill', 'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']
critere1 = {"nelements": {'$lte': 6}, 'elements': {'$all': composes}, "elasticity": {'$ne': None}, "elasticity.G_Reuss": {'$gte': 0 }, "elasticity.G_Voigt": {'$gte': 0 }, "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0 }, "elasticity.K_Reuss": {'$gte': 0 }, "elasticity.K_Voigt": {'$gte': 0 }, "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0 }}


materials = api.query(criteria=critere1, properties=propsTableau)

propsPlot = ['elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill', 'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']



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

resultat = recup(materials)

def log():

    propsPlot1 = []

    for i in :
        propsPlot1.append(math.log(i))




Greuss = resultat[propsTableau.index('elasticity.G_Reuss'), :]

mini=min(Greuss)

poisson = resultat[propsTableau.index('elasticity.poisson_ratio'), :]
normalize = color.Normalize(vmin=min(poisson), vmax=max(poisson))
for prop1 in propsPlot1:
    for prop2 in propsPlot1:
        if prop1 != prop2:
            x = resultat[propsTableau.index(prop1), :]
            y = resultat[propsTableau.index(prop2), :]
            area = 5  # 0 to 15 point radii
            plt.scatter(x, y, s=area, c=poisson,cmap=cm.get_cmap('seismic'),  norm=normalize, alpha=1)
            #plt.xlim(x.min() * 1.1, x.max() * 1.1)
            #plt.ylim(y.min() * 1.1, y.max() * 1.1)
            plt.xlim(x.min(), 1000)
            plt.ylim(y.min(), 1000)
            plt.xlabel(prop1)
            plt.ylabel(prop2)
            plt.title(str(prop2) + ' versus ' +str(prop1))
            plt.colorbar()
            #filename= 'C:\\Users\\siwar\\Desktop\\image\\'+str(prop2) +' versus '+str(prop1)+'.pdf'
            #if os.path.isfile(filename):
            #    os.remove(filename)  # Opt.: os.system("rm "+strFile)
            pdf.savefig()
            plt.close()
pdf.close()


