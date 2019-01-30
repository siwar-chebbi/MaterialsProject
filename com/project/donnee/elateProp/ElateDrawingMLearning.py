import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.backends.backend_pdf

# correspond a elate_properties_all_materials_with_crashes_cases.py
# propsDisplay = ["minLC", "maxLC", "minNu", "maxNu", "G_Voigt_Reuss_Hill", "K_Voigt_Reuss_Hill"]

# correspond a elate_properties_all_materials_filtered.py
propsDisplay = ["minLC", "maxLC", "minNu", "maxNu", "K_Voigt_Reuss_Hill", "Emin", "Emax", "Gmin", "Gmax"]
propsPlotLabel = [u'$LC_{min} (GPa)$', u'$LC_{max}(GPa)$', u'$\mu_{min}(GPa)$', u'$\mu_{max}(GPa)$',
                  u'$K_{Voigt\u2000Reuss\u2000Hill}(GPa)$', u'$E_{min}(GPa)$', '$E_{max}(GPa)$', u'$G_{min}(GPa)$',
                  '$G_{max}(GPa)$']
pdf = matplotlib.backends.backend_pdf.PdfPages("elastic_All_EXP_AND_HYP_Filtered.pdf")


def importer(fichier):
    return pd.read_csv(fichier)


data = importer("elastic_All_EXP_AND_HYP_Filtered.csv")
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
            plt.scatter(data_X, data_Y, color='black')
            plt.plot(data_X, data_y_pred, color='blue', linewidth=3)
            plt.xlim(data_X.min(), data_X.max() * 1.1)
            plt.ylim(data_Y.min(), data_Y.max() * 1.1)
            # plt.xlim(data_X.min(), 1800)
            # plt.ylim(data_Y.min(), 1800)
            plt.xlabel(propsPlotLabel[propsDisplay.index(prop1)])
            plt.ylabel(propsPlotLabel[propsDisplay.index(prop2)])
            # plt.title(str(prop2) + ' versus ' + str(prop1))
            plt.figtext(0.5, 0.75, texte, ha="center", fontsize=7, bbox={"facecolor": "orange", "alpha": 0.5, "pad": 5})
            pdf.savefig()
            plt.close()
pdf.close()
