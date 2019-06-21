import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.backends.backend_pdf

# correspond a elate_properties_all_materials_with_crashes_cases.py
# propsDisplay = ["minLC", "maxLC", "minNu", "maxNu", "G_Voigt_Reuss_Hill", "K_Voigt_Reuss_Hill"]

# correspond a elate_properties_all_materials_filtered.py
propsDisplay = ['elasticity.poisson_ratio']

propsPlotLabel = [u'$G_{Reuss} (GPa)$', u'$G_{Voigt}(GPa)$', u'$G_{Voigt\u2000Reuss\u2000Hill}(GPa)$',
                  u'$K_{Reuss}(GPa)$', '$K_{Voigt}(GPa)$', u'$K_{Voigt\u2000Reuss\u2000Hill}(GPa)$']

pdf = matplotlib.backends.backend_pdf.PdfPages("RatioPoissonAnisotropy.pdf")


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

data_X = Emax_sur_Emin

print(len(data_X))
# data_X.shape
data_X = np.vstack(data_X)

data_Y = data['elasticity.poisson_ratio'].get_values()

area = 5
fig = plt.figure()
ax1 = fig.add_subplot(111, label="log10")
ax1.set_xlim(0.80, 300)
ax1.set_ylim(data_Y.min(), data_Y.max())
ax1.set_xlabel("Elastic anisotropy")
ax1.set_ylabel(u'$\mu$')
ax1.set_xscale("log")
ax1.scatter(data_X, data_Y, s=area, alpha=1, color="green")

ax1.axvspan(1, 2, facecolor='r', alpha=0.5)

pdf.savefig()
plt.close()
pdf.close()
