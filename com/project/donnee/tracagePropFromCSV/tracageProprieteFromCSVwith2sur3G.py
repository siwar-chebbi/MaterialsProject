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

propsDisplay = ["minLC", "maxLC", "minNu", "maxNu", "Emin", "Emax", "Gmin", "Gmax",
                'elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill',
                'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']

propsPlotLabel = [u'$LC_{min} (GPa)$', u'$LC_{max}(GPa)$', u'$\mu_{min}(GPa)$', u'$\mu_{max}(GPa)$', u'$E_{min}(GPa)$',
                  '$E_{max}(GPa)$', u'$G_{min}(GPa)$', '$G_{max}(GPa)$',
                  u'$G_{Reuss} (GPa)$', u'$G_{Voigt}(GPa)$', u'$G_{Voigt\u2000Reuss\u2000Hill}(GPa)$',
                  u'$K_{Reuss}(GPa)$', '$K_{Voigt}(GPa)$', u'$K_{Voigt\u2000Reuss\u2000Hill}(GPa)$']

propsPlotLabelSansGPA = [u'$G_{Reuss} $', u'$G_{Voigt}$', u'$G_{Voigt\u2000Reuss\u2000Hill}$',
                         u'$K_{Reuss}$', '$K_{Voigt}$', u'$K_{Voigt\u2000Reuss\u2000Hill}$']


def importer(fichier):
    return pd.read_csv(fichier)


# fichiers input (csv) et output (pdf)
def generation_pdf(csv_source, pdf_file):
    data = importer(csv_source)
    data.head()
    pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_file)
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
                    if x <= 0 or y <= 0:
                        continue
                    else:
                        cleaned_x.append(x)
                        cleaned_y.append(y)
                        cleaned_poisson.append(z)

                data_X2 = np.vstack(cleaned_x)
                data_Y2 = cleaned_y
                cleaned_x2sur3 = [i * 1 / 3 for i in data_X2]
                # cleaned_x8sur3 = [i * 8 / 3 for i in data_X_log]

                # regession lineaire de log10(y) =f(log10(x))
                regr = linear_model.LinearRegression()
                regr.fit(data_X2, data_Y2)
                data_y_pred = regr.predict(sorted(data_X2))

                # y = a*x + b
                # Explained variance score: 1 is perfect prediction
                # The mean squared error
                print('############# ' + str(prop2) + ' versus ' + str(prop1) + ' ######################')
                print(
                    'y =  {:.2f} * x + {:.2f} \nVariance score: {:.2f} \nMean squared error:{:.2f} \nNombre de points: {:d}\n'.format(
                        regr.coef_[0], regr.intercept_, r2_score(data_Y2, data_y_pred),
                        mean_squared_error(data_y_pred, data_Y2), len(cleaned_x)))

                # texte dans le graphe
                # texte = u'$log_{10}($' + propsPlotLabel[propsDisplay.index(prop2)] + ') = ' + "{:.2f}".format(
                #    regr.coef_[0]) + ' * ' + u'$log_{10}($' + propsPlotLabel[
                #            propsDisplay.index(prop1)] + ') + ' + "{:.2f}".format(
                #    regr.intercept_) + ' \n Variance score R2: ' + "{:.2f}".format(
                #    r2_score(data_Y_log, data_y_pred)) \
                # + '\n Mean squared error: ' + "{:.2f}".format( mean_squared_error(data_y_pred, data_Y_log))

                area = 0.1
                # 2 subplots superposees
                fig = plt.figure()
                ax1 = fig.add_subplot(111, label="log10")
                ax2 = fig.add_subplot(111, label="regression", frame_on=False)
                # rajouter 2sur3 de G
                ax3 = fig.add_subplot(111, label="1sur3", frame_on=False)
                # ax4 = fig.add_subplot(111, label="8sur3", frame_on=False)

                # subplot de tous les points
                normalize = color.Normalize(vmin=min(cleaned_poisson), vmax=max(cleaned_poisson))
                toto = color.LinearSegmentedColormap.from_list("toto", (
                    (0.0, 0.0, 0.3), (0.0, 0.0, 1.0),
                    (0.6, 0.6, 0.6), (1.0, 0.0, 0.0),
                    (0.5, 0.0, 0.0)))
                im = ax1.scatter(cleaned_x, cleaned_y, s=area, c=cleaned_poisson, norm=normalize, alpha=3, cmap=toto)
                # im = ax1.scatter(cleaned_x, cleaned_y, s=area, c=cleaned_poisson, cmap=cm.get_cmap('coolwarm'),
                #                 norm=normalize, alpha=3)
                # ax1.set_xscale('log')
                # ax1.set_yscale('log')
                ax1.set_xlim(0, 1e3)
                ax1.set_ylim(0, 1e3)

                # subplot regression lineaire (droite)
                ax2.plot(sorted(data_X2), data_y_pred, color='black', linewidth=1)
                ax2.set_xlim(0, 1e3)
                ax2.set_ylim(0, 1e3)
                ax2.set_yticklabels([])
                ax2.set_xticklabels([])

                # subplot 1/3  (droite)
                ax3.plot(sorted(data_X2), sorted(cleaned_x2sur3), color='green', linewidth=1)
                ax3.text(1.9, 0.7, '1/3G', fontsize=6, color='green', rotation=19)
                ax3.set_xlim(0, 1e3)
                ax3.set_ylim(0, 1e3)
                ax3.set_yticklabels([])
                ax3.set_xticklabels([])

                # subplot 8/3  (droite)
                # ax4.plot(sorted(data_X_log), sorted(cleaned_x8sur3), "--", color='orange', linewidth=2)
                # ax4.text(0.78, 2.6, '8/3G', fontsize=10, color='orange', rotation= 60)
                # ax4.set_xlim(0, 3)
                # ax4.set_ylim(0, 3)
                # ax4.set_yticklabels([])
                # ax4.set_xticklabels([])

                # colorbar des valeurs de poisson
                plt.gcf().subplots_adjust(right=0.8)
                cbar_ax = plt.gcf().add_axes([0.82, 0.125, 0.03, 0.75])
                cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
                cbar.set_label("Poisson's ratio")

                ax1.set_xlabel(propsPlotLabel[propsDisplay.index(prop1)])
                ax1.set_ylabel(propsPlotLabel[propsDisplay.index(prop2)])
                # plt.figtext(0.5, 0.80, ha="center", fontsize=7, bbox={"facecolor": "orange", "alpha": 0.5, "pad": 5})
                # plt.figtext(0.5, 0.80, texte, ha="center", fontsize=7, bbox={"facecolor": "orange", "alpha": 0.5, "pad": 5})
                pdf.savefig()
                plt.close()
    pdf.close()


#generation_pdf("elasticElate_ALL_revisionArt_without_Zero.csv", "elasticElate_ALL_revisionArt_without_LOG_AND_Zero.pdf")
print("#################################################################################################\n""#################################################################################################\n")
generation_pdf("elasticElate_ALL_revisionArt_without_Zero_EXP.csv", "elasticElate_ALL_revisionArt_without_LOG_AND_Zero_EXP.pdf")
print("#################################################################################################\n""#################################################################################################\n")
generation_pdf("elasticElate_ALL_revisionArt_without_Zero_HYP.csv", "elasticElate_ALL_revisionArt_without_LOG_AND_Zero_HYP.pdf")
