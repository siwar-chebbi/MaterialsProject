from pymatgen import MPRester
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as color
import matplotlib.backends.backend_pdf
import numpy as np

pdf = matplotlib.backends.backend_pdf.PdfPages("output.pdf")
api = MPRester("fB610TDF3LSwxiN9")


#Proprietes utilisees dans la requete
propsTableauCritere = ['pretty_formula', 'elasticity.poisson_ratio', 'elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill', 'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']

#Proprietes utilisees dans la generation du tableau
propsTableau = ['elasticity.poisson_ratio', 'elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill', 'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']

#Proprietes utilisees dans le tracage des graphes
propsPlot = ['elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill', 'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']


critere1 = {"nelements": {'$lte': 6}, "elasticity": {'$ne': None}, "elasticity.G_Reuss": {'$gte': 0 }, "elasticity.G_Voigt": {'$gte': 0 }, "elasticity.G_Voigt_Reuss_Hill": {'$gte': 0 }, "elasticity.K_Reuss": {'$gte': 0 }, "elasticity.K_Voigt": {'$gte': 0 }, "elasticity.K_Voigt_Reuss_Hill": {'$gte': 0 }}



#requete
materials = api.query(criteria=critere1, properties=propsTableauCritere)

#lin= len(propsTableauCritere)
#dimensions du tableau
lin = len(propsTableau)
col = len(materials)

#generation des valeurs correspondantes aux proprietes des elements
def recup(materials):
    j = 0
    tableau = np.zeros(shape=(lin, col))
    #elements = []

    for material in materials:

     #  elements.append(material.get('pretty_formula'))
        i = 0
        for prop in propsTableau:
            tableau[i, j] = material.get(prop)
            i = i + 1
        j = j + 1
    return tableau

resultat = recup(materials)

#poisson = resultat[propsTableau.index('elasticity.poisson_ratio'), :]
#normalize = color.Normalize(vmin=min   (poisson), vmax=max(poisson))
for prop in propsTableau:
    #for prop2 in propsPlot:
        #if prop1 != prop2:
    data = resultat[:, propsTableau.index(prop)]
    #data = resultat[propsTableau.index(prop), :]
    #y = resultat[propsTableau.index(prop2), :]
    #area = 5  # 0 to 15 point radii
    plt.hist(data)
    #plt.scatter(x, y, s=area, c=poisson,cmap=cm.get_cmap('seismic'),  norm=normalize, alpha=1)
    #plt.xlim(x.min(), x.max() * 1.1)
    #plt.ylim(y.min(), y.max() * 1.1)
    #plt.xlabel(prop1[11:])
    #plt.ylabel(prop2[11:])
    plt.title(str(prop[11:]))
    #plt.colorbar()
    #filename= 'C:\\Users\\siwar\\Desktop\\image\\'+str(prop2) +' versus '+str(prop1)+'.pdf'
    #if os.path.isfile(filename):
    #    os.remove(filename)  # Opt.: os.system("rm "+strFile)
    pdf.savefig()
    plt.close()
pdf.close()



