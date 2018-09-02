from pandas._libs import properties
from pymatgen import MPRester
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as color
import os
import matplotlib.backends.backend_pdf

pdf = matplotlib.backends.backend_pdf.PdfPages("C:\\Users\\siwar\\Desktop\\image\\output.pdf")
# from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter

from pymatgen.analysis.elasticity import elastic
from pymatgen.util import plotting
import numpy as np

api = MPRester("fB610TDF3LSwxiN9")

#entries = api.get_entries({"nelements": {'$lte': 6, '$gte': 1},"elasticity": {'$ne': None}}, property_data=['pretty_formula','elasticity', 'elements'])
covalent = ['B', 'C', 'Si']
ionique = ['N', 'O', 'F', 'P', 'S', 'Cl', 'Sr', 'Br', 'I']

propsTableau = ['elasticity.poisson_ratio','elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill', 'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']
critere1= {"nelements": {'$gte': 1, '$lte': 6}, "elasticity": {'$ne': None}}
critere2 = {"nelements": {'$lte': 6}, 'elements': {'$in': covalent}, "elasticity": {'$ne': None}}
critere3 = {"nelements": {'$lte': 6}, 'elements': {'$in': ionique}, "elasticity": {'$ne': None}}
#critere4 = {"nelements": {'$lte': 6}, "elasticity": {'$ne': None}, "elasticity.G_Reuss": {'$gte': 0 }}
materials = api.query(criteria=critere1, properties=propsTableau)
#materials_covalent= api.query(criteria=critere2, properties=propsTableau)
#materials_ionique = api.query(criteria=critere3, properties=propsTableau)

propsPlot = ['elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill', 'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']
#propsPlot = ['elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill']




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

Greuss = resultat[propsTableau.index('elasticity.G_Reuss'), :]

mini=min(Greuss)

poisson = resultat[propsTableau.index('elasticity.poisson_ratio'), :]
normalize = color.Normalize(vmin=min(poisson), vmax=max(poisson))
for prop1 in propsPlot:
    for prop2 in propsPlot:
        if prop1 != prop2:
            x = resultat[propsTableau.index(prop1), :]
            y = resultat[propsTableau.index(prop2), :]
            area = 5  # 0 to 15 point radii
            plt.scatter(x, y, s=area, c=poisson,cmap=cm.get_cmap('plasma'),  norm=normalize, alpha=1)
            plt.xlim(x.min() * 1.1, x.max() * 1.1)
            plt.ylim(y.min() * 1.1, y.max() * 1.1)
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



#entries = api.get_entries({"nelements": {'$lte': 6, '$gte': 6},"elasticity": {'$ne': None}}, property_data=['pretty_formula','elasticity', 'elements'])


#datax = api.query(criteria={"nelements": {'$lte': 6, '$gte': 6}, "elasticity": {'$ne': None}}, properties=["pretty_formula","elasticity", "elements"])
#elasticity.G_Reuss


#datax = api.query(criteria={"nelements": {'$lte': 6, '$gte': 6}, "elasticity": {'$ne': None}}, properties=["pretty_formula","elasticity"])

#datay = api.query(criteria={"nelements": {'$lte': 1}, "elasticity": {'$ne': None}}, properties=["elasticity.K_VRH"])



#area = 5  # 0 to 15 point radii
#plt.scatter(x, y, s=area, alpha=0.5)
#plt.show()

# N = len(data)
# print(N)
#
# entries = api.get_entries({"nelements": {'$lte': 6}, "elements": {'$all': ['S', 'O']}, "elasticity": {'$ne': None}},
#                           property_data=['elasticity'])
# N = len(entries)
# print(N)
#
# x = list()
# y = list()
#
# for entry in entries:
#     if entry.data["elasticity"]:
#         x.append(entry.data["elasticity"]['G_Reuss'])
#         y.append(entry.data["elasticity"]['K_VRH'])
#
# area = 5  # 0 to 15 point radii
# plt.scatter(x, y, s=area, alpha=0.5)
# plt.show()
