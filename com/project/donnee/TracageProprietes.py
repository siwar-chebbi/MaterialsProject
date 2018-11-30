from pandas._libs import properties
from pymatgen import MPRester
import matplotlib.pyplot as plt
# from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter

from pymatgen.analysis.elasticity import elastic
from pymatgen.util import plotting
import numpy as np

api = MPRester("eDCEK5m9WVjmajp7e8af")


#datax = api.query(criteria={"nelements": {'$lte': 6 ,'$gte': 1 }, "elements": {'$all': ['S','O']}, "elasticity": {'$ne': None}}, properties=['pretty_formula','elasticity.G_Reuss', 'elasticity.G_VRH', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill', 'elasticity.K_Reuss', 'elasticity.K_VRH', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill'])

materials = api.query(criteria={"nelements": {'$lte': 6 ,'$gte': 1}, "elasticity": {'$ne': None}}, properties=['pretty_formula','elasticity.G_Reuss', 'elasticity.G_VRH', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill', 'elasticity.K_Reuss', 'elasticity.K_VRH', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill'])

props = ['elasticity.G_Reuss', 'elasticity.G_VRH', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill', 'elasticity.K_Reuss', 'elasticity.K_VRH', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']
lin = len(props)
col = len(materials)



def recup(materials):

    j = 0
    tableau = np.zeros(shape=(lin, col))
    elements = []

    for material in materials:

        elements.append(material.get('pretty_formula'))
        i = 0
        for prop in props:
            tableau[i, j] = material.get(prop)
            i = i + 1
        j = j + 1
    return tableau

resultat = recup(materials)


for prop1 in props:
    for prop2 in props:
        if prop1 != prop2:
            x = resultat[props.index(prop1), :]
            y = resultat[props.index(prop2), :]
            area = 5  # 0 to 15 point radii
            plt.scatter(x, y, s=area, alpha=0.5)
            plt.xlim(x.min() * 1.1, x.max() * 1.1)
            plt.ylim(y.min() * 1.1, y.max() * 1.1)
            plt.xlabel(prop1)
            plt.ylabel(prop2)
            plt.title(str(prop2) + ' versus ' +str(prop1))
            plt.show()
