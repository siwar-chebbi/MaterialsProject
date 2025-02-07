import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.backends.backend_pdf
import matplotlib.cm as cm
import matplotlib.colors as color
import math

# correspond a elate_properties_all_materials_with_crashes_cases.py
# propsDisplay = ["minLC", "maxLC", "minNu", "maxNu", "G_Voigt_Reuss_Hill", "K_Voigt_Reuss_Hill"]

# correspond a elate_properties_all_materials_filtered.py
propsDisplay = ['minLC', 'maxLC', 'minNu', 'maxNu', 'elasticity.K_Voigt_Reuss_Hill',
       'Emin', 'Emax', 'Gmin', 'Gmax', 'elasticity.G_Reuss',
       'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill',
       'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.poisson_ratio',
       'lvpa', 'cepa', 'group1', 'atomic_mass1', 'atomicRadius1', 'rowH1A',
       'rowHn3A', 'xH4A', 'xHn4A', 'cohesive_energy', 'average_electroneg',
       'bandgap', 'density', 'formation_energy-peratom', 'e_above_hull']


propsPlotLabel = ['minLC', 'maxLC', 'minNu', 'maxNu', 'elasticity.K_Voigt_Reuss_Hill',
       'Emin', 'Emax', 'Gmin', 'Gmax', 'elasticity.G_Reuss',
       'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill',
       'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.poisson_ratio',
       'lvpa', 'cepa', 'group1', 'atomic_mass1', 'atomicRadius1', 'rowH1A',
       'rowHn3A', 'xH4A', 'xHn4A', 'cohesive_energy', 'average_electroneg',
       'bandgap', 'density', 'formation_energy-peratom', 'e_above_hull']

#propsDisplay = ['elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill', 'elasticity.K_Reuss',
               # 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']

#propsPlotLabel = [u'$G_{Reuss} (GPa)$', u'$G_{Voigt}(GPa)$', u'$G_{Voigt\u2000Reuss\u2000Hill}(GPa)$',
                 # u'$K_{Reuss}(GPa)$', '$K_{Voigt}(GPa)$', u'$K_{Voigt\u2000Reuss\u2000Hill}(GPa)$']

pdf = matplotlib.backends.backend_pdf.PdfPages("elastic_property_from_MP_DB12522.pdf")


def importer(fichier):
    return pd.read_csv(fichier)


data = importer("Extract_nagative_minNu_maxNu_poisson_descriptors.csv")
data.head()

#poisson = data['elasticity.poisson_ratio'].get_values()
#normalize = color.Normalize(vmin=min(poisson), vmax=max(poisson))
for prop1 in propsDisplay:
    for prop2 in propsDisplay:
        if prop1 != prop2:
            data_X = data[prop1].get_values()
            # data_X.shape
            data_X = np.vstack(data_X)

            data_Y = data[prop2].get_values()
            # data_Y.shape
            regr = linear_model.LinearRegression()

            # Train the model using the training sets
            regr.fit(data_X, data_Y)

            # Make predictions using the testing set
            data_y_pred = regr.predict(data_X)

            # The coefficients
            print('#############' + str(prop2) + ' versus ' + str(prop1) + '######################')
            print('Coefficients: \n', regr.coef_)
            # The mean squared error
            print("Mean squared error: %.2f"
                  % mean_squared_error(data_y_pred, data_Y))
            # Explained variance score: 1 is perfect prediction
            print('Variance score: %.2f \n '
                  '##############################################################' % r2_score(data_Y, data_y_pred))

            texte = 'Coefficients: ' + "{:.2f}".format(regr.coef_[0]) + ' \n Mean squared error: ' + "{:.2f}".format(
                mean_squared_error(data_y_pred, data_Y)) + ' \n Variance score R2: ' + "{:.2f}".format(
                r2_score(data_Y, data_y_pred))
            area = 5
            plt.scatter(data_X, data_Y, s=area, cmap=cm.get_cmap('seismic'), alpha=1)

            #plt.scatter(data_X, data_Y, s=area, c=poisson, cmap=cm.get_cmap('seismic'), norm=normalize, alpha=1)

            plt.plot(data_X, data_y_pred, color='black', linewidth=2)
            # plt.xlim(data_X.min(), data_X.max() * 1.1)
            # plt.ylim(data_Y.min(), data_Y.max() * 1.1)
            #plt.xlim(data_X.min(), 1000)
            #plt.ylim(data_Y.min(), 1000)
            plt.xlabel(propsPlotLabel[propsDisplay.index(prop1)])
            plt.ylabel(propsPlotLabel[propsDisplay.index(prop2)])
            #plt.colorbar()
            # plt.title(str(prop2) + ' versus ' + str(prop1))
            plt.figtext(0.5, 0.75, texte, ha="center", fontsize=7, bbox={"facecolor": "orange", "alpha": 0.5, "pad": 5})
            pdf.savefig()
            plt.close()
pdf.close()
