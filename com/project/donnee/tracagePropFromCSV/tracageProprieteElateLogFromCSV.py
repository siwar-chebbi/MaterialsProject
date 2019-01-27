import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.backends.backend_pdf
import matplotlib.cm as cm
import matplotlib.colors as color
import math


def importer(fichier):
    return pd.read_csv(fichier)


# correspond a elate_properties_all_materials_with_crashes_cases.py
# propsDisplay = ["minLC", "maxLC", "minNu", "maxNu", "G_Voigt_Reuss_Hill", "K_Voigt_Reuss_Hill"]

# correspond a elate_properties_all_materials_filtered.py
propsDisplay = ["minLC", "maxLC", "minNu", "maxNu", "K_Voigt_Reuss_Hill", "Emin", "Emax", "Gmin", "Gmax"]

propsPlotLabel = [u'$LC_{min} (GPa)$', u'$LC_{max}(GPa)$', u'$\mu_{min}(GPa)$', u'$\mu_{max}(GPa)$',
                  u'$K_{Voigt\u2000Reuss\u2000Hill}(GPa)$', u'$E_{min}(GPa)$', '$E_{max}(GPa)$', u'$G_{min}(GPa)$',
                  '$G_{max}(GPa)$']

# fichiers input (csv) et output (pdf)
data = importer("elasticRatioPoissonPositive.csv")
data.head()
pdf = matplotlib.backends.backend_pdf.PdfPages("elasticRatioPoissonPositive.pdf")

# valeurs poisson
poisson = data['elasticity.poisson_ratio'].get_values()

for prop1 in propsDisplay:
    for prop2 in propsDisplay:
        if prop1 != prop2:
            # X et Y
            data_X = data[prop1].get_values()
            data_Y = data[prop2].get_values()
            # X et Y sans les 0 (pb log10)
            cleaned_x = []
            cleaned_y = []
            cleaned_poisson = []
            for x, y, z in zip(data_X, data_Y, poisson):
                if x == '-inf' or y == '-inf' in y or x <= 0 or y <= 0:
                    continue
                else:
                    cleaned_x.append(x)
                    cleaned_y.append(y)
                    cleaned_poisson.append(z)
            # log10 de X et Y
            data_X_log = np.vstack(np.log10(cleaned_x))
            data_Y_log = np.log10(cleaned_y)
            # regession lineaire de log10(y) =f(log10(x))
            regr = linear_model.LinearRegression()
            regr.fit(data_X_log, data_Y_log)
            data_y_pred = regr.predict(sorted(data_X_log))

            # log10(y) = a*log10(x) + b
            # Explained variance score: 1 is perfect prediction
            # The mean squared error
            print('############# ' + str(prop2) + ' versus ' + str(prop1) + ' ######################')
            print(
                'log10 =  {:.2f} * log10(x) + {:.2f} \nVariance score: {:.2f} \nMean squared error:{:.2f} \nNombre de points: {:d}\n'.format(
                    regr.coef_[0], regr.intercept_, r2_score(data_Y_log, data_y_pred),
                    mean_squared_error(data_y_pred, data_Y_log), len(cleaned_x)))

            # texte dans le graphe
            texte = u'$log_{10}($' + propsPlotLabel[propsDisplay.index(prop2)] + ') = ' + "{:.2f}".format(
                regr.coef_[0]) + ' * ' + u'$log_{10}($' + propsPlotLabel[
                        propsDisplay.index(prop1)] + ') + ' + "{:.2f}".format(
                regr.intercept_) + ' \n Variance score R2: ' + "{:.2f}".format(
                r2_score(data_Y_log, data_y_pred)) \
                # + '\n Mean squared error: ' + "{:.2f}".format( mean_squared_error(data_y_pred, data_Y_log))

            area = 5
            max_x_log = 1e3
            max_y_log = 1e3
            min_x_log = 1
            min_y_log = 1
            # 2 subplots superposees
            fig = plt.figure()
            ax1 = fig.add_subplot(111, label="log10")
            ax2 = fig.add_subplot(111, label="regression", frame_on=False)
            # subplot de tous les points
            normalize = color.Normalize(vmin=min(cleaned_poisson), vmax=max(cleaned_poisson))
            im = ax1.scatter(cleaned_x, cleaned_y, s=area, c=cleaned_poisson, cmap=cm.get_cmap('seismic'),
                             norm=normalize, alpha=1)
            ax1.set_xscale('log')
            ax1.set_yscale('log')
            ax1.set_xlim(min_x_log, max_x_log)
            ax1.set_ylim(min_y_log, max_y_log)

            # subplot regression lineaire (droite)
            ax2.plot(sorted(data_X_log), data_y_pred, color='black', linewidth=2)
            ax2.set_xlim(math.log10(min_x_log), math.log10(max_x_log))
            ax2.set_ylim(math.log10(min_y_log), math.log10(max_y_log))
            ax2.set_yticklabels([])
            ax2.set_xticklabels([])

            # colorbar des valeurs de poisson
            plt.gcf().subplots_adjust(right=0.8)
            cbar_ax = plt.gcf().add_axes([0.82, 0.125, 0.03, 0.75])
            fig.colorbar(im, cax=cbar_ax, orientation='vertical')

            ax1.set_xlabel(propsPlotLabel[propsDisplay.index(prop1)])
            ax1.set_ylabel(propsPlotLabel[propsDisplay.index(prop2)])
            plt.figtext(0.5, 0.80, texte, ha="center", fontsize=7, bbox={"facecolor": "orange", "alpha": 0.5, "pad": 5})
            pdf.savefig()
            plt.close()
pdf.close()
