import sys
from sklearn.ensemble import RandomForestRegressor
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_predict, cross_val_score

# permet d'imrpimer la liste complète des valeurs
np.set_printoptions(threshold=sys.maxsize)
feature_names = np.array(
    ['eGrpHn4A', 'eGrpHn3A', 'eGrpHn2A', 'eGrpHn1A', 'eGrpH0A', 'eGrpH1A', 'eGrpH2A', 'eGrpH3A', 'eGrpH4A',
     'eGrpH0S', 'eGrpH1S', 'eAtMHn4A', 'eAtMHn3A', 'eAtMHn2A', 'eAtMHn1A', 'eAtMH0A', 'eAtMH1A', 'eAtMH2A',
     'eAtMH3A', 'eAtMH4A', 'eAtMH0S', 'eAtMH1S', 'eZHn4A', 'eZHn3A', 'eZHn2A', 'eZHn1A', 'eZH0A', 'eZH1A', 'eZH2A',
     'eZH3A', 'eZH4A', 'eZH0S', 'eZH1S', 'eRadHn4A', 'eRadHn3A', 'eRadHn2A', 'eRadHn1A', 'eRadH0A', 'eRadH1A',
     'eRadH2A', 'eRadH3A', 'eRadH4A', 'eRadH0S', 'eRadH1S', 'eRowHn4A', 'eRowHn3A', 'eRowHn2A', 'eRowHn1A',
     'eRowH0A', 'eRowH1A', 'eRowH2A', 'eRowH3A', 'eRowH4A', 'eRowH0S', 'eRowH1S', 'eXHn4A', 'eXHn3A', 'eXHn2A',
     'eXHn1A', 'eXH0A', 'eXH1A', 'eXH2A', 'eXH3A', 'eXH4A', 'eXH0S', 'eXH1S', 'eBlPtHn4A', 'eBlPtHn3A', 'eBlPtHn2A',
     'eBlPtHn1A', 'eBlPtH0A', 'eBlPtH1A', 'eBlPtH2A', 'eBlPtH3A', 'eBlPtH4A', 'eBlPtH0S', 'eBlPtH1S', 'eMlPtHn4A',
     'eMlPtHn3A', 'eMlPtHn2A', 'eMlPtHn1A', 'eMlPtH0A', 'eMlPtH1A', 'eMlPtH2A', 'eMlPtH3A', 'eMlPtH4A', 'eMlPtH0S',
     'eMlPtH1S', 'lvpa', 'cepa', 'cohesive_energy', 'average_electroneg', 'bandGap', 'rho', 'fepa', 'eah',
     'sCoorHn4A', 'sCoorHn3A', 'sCoorHn2A', 'sCoorHn1A', 'sCoorH0A', 'sCoorH1A', 'sCoorH2A', 'sCoorH3A', 'sCoorH4A',
     'sCoorH0S', 'sCoorH1S', 'sBnLnHn4AH1A', 'sBnLnHn3AH1A', 'sBnLnHn2AH1A', 'sBnLnHn1AH1A', 'sBnLnH0AH1A',
     'sBnLnH1AH1A', 'sBnLnH2AH1A', 'sBnLnH3AH1A', 'sBnLnH4AH1A', 'sBnLnH0SH1A', 'sBnLnH1SH1A', 'sRowADH0AH1A',
     'sRowADH1AH1A', 'sRowADH2AH1A', 'sRowADH3AH1A', 'sRowADH4AH1A', 'sRowADH1SH1A', 'sRowSDH1AH1A', 'sRowSDH2AH1A',
     'sRowSDH4AH1A', 'sRowSDH1SH1A', 'sGrpADH0AH1A', 'sGrpADH1AH1A', 'sGrpADH2AH1A', 'sGrpADH3AH1A', 'sGrpADH4AH1A',
     'sGrpADH1SH1A', 'sGrpSDH1AH1A', 'sGrpSDH2AH1A', 'sGrpSDH4AH1A', 'sGrpSDH1SH1A', 'sAtMADH0AH1A', 'sAtMADH1AH1A',
     'sAtMADH2AH1A', 'sAtMADH3AH1A', 'sAtMADH4AH1A', 'sAtMADH1SH1A', 'sAtMSDH1AH1A', 'sAtMSDH2AH1A', 'sAtMSDH4AH1A',
     'sAtMSDH1SH1A', 'sRadADH0AH1A', 'sRadADH1AH1A', 'sRadADH2AH1A', 'sRadADH3AH1A', 'sRadADH4AH1A', 'sRadADH1SH1A',
     'sRadSDH1AH1A', 'sRadSDH2AH1A', 'sRadSDH4AH1A', 'sRadSDH1SH1A', 'sXADH0AH1A', 'sXADH1AH1A', 'sXADH2AH1A',
     'sXADH3AH1A', 'sXADH4AH1A', 'sXADH1SH1A', 'sXSDH1AH1A', 'sXSDH2AH1A', 'sXSDH4AH1A', 'sXSDH1SH1A', 'sZADH0AH1A',
     'sZADH1AH1A', 'sZADH2AH1A', 'sZADH3AH1A', 'sZADH4AH1A', 'sZADH1SH1A', 'sZSDH1AH1A', 'sZSDH2AH1A', 'sZSDH4AH1A',
     'sZSDH1SH1A'])

data = pd.read_csv('Extract_Allvalues_descriptors_sans_angles.csv', index_col=0)

print("\nliste des propriétés dans le tableau à traiter \n")
print(str(data.columns))

y = data['elasticity.poisson_ratio']
# y = data['minNu']
# y = data['maxNu']
# without elastic properties
X = data.drop(
    ['elasticity.poisson_ratio', 'minLC', 'maxLC', 'minNu', 'maxNu', 'elasticity.K_Voigt_Reuss_Hill', 'Emin',
     'Emax', 'Gmin', 'Gmax', 'elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill',
     'elasticity.K_Reuss', 'elasticity.K_Voigt'], axis=1)

# without Numax and Numin
# X = data.drop(['elasticity.poisson_ratio', 'cepa', 'minLC',	'maxLC','elasticity.K_Voigt_Reuss_Hill',	'Emin',	'Emax',	'Gmin',	'Gmax',	'elasticity.G_Reuss',	'elasticity.G_Voigt',	'elasticity.G_Voigt_Reuss_Hill',	'elasticity.K_Reuss',	'elasticity.K_Voigt'], axis=1)


# only we kept KVRH and GVRH
# X = data.drop(['elasticity.poisson_ratio', 'cepa', 'minLC',	'maxLC', 'minNu','maxNu',	'Emin',	'Emax',	'Gmin',	'Gmax',	'elasticity.G_Reuss',	'elasticity.G_Voigt',	'elasticity.K_Reuss',	'elasticity.K_Voigt'], axis=1)

# With all elastic properties
# X = data.drop(['maxNu', 'cepa'], axis=1)
print("\ndimension du tableau des données à traiter ")
print(str(np.shape(X)))

# y = np.log10(y)


lr = linear_model.LinearRegression()

#gbr = GradientBoostingRegressor(n_estimators=1000, learning_rate=0.01, min_samples_split=2, min_samples_leaf=3,max_depth=3, max_features='sqrt', loss='ls', subsample=0.4)
# gbr = GradientBoostingRegressor(n_estimators=1000, learning_rate=0.01)
#
# gbr = AdaBoostRegressor(DecisionTreeRegressor(max_depth=4),
#                         n_estimators=300, random_state=0)
gbr = RandomForestRegressor(n_estimators=100, criterion='mse', min_samples_split=2, min_samples_leaf=3, max_depth=3, max_features='sqrt')

scores = []
iterations = 100
plt.figure(figsize=(12, 6))
predict_result = []
for i in range(0, iterations):
    print("\niteration numero : " + str(i) + " / " + str(iterations))
    cv = KFold(n_splits=3, random_state=i)
    predict = cross_val_predict(gbr, X, y, cv=cv, n_jobs=-1)
    score = cross_val_score(gbr, X, y, cv=cv, scoring='neg_mean_squared_error', n_jobs=-1)
    score = (score * -1) ** 0.5
    scores.append(score)
    if i == (iterations - 1):
        predict_result = predict

scores = np.array(scores)
print("Le resultat de prédiction de ration de poisson est :\n" + str(predict_result) + "\n")
print("\n Moyenne: \n" + str(np.mean(scores.flatten())))
print("\n Ecart-type: \n" + str(np.std(scores.flatten())))
print("\n Erreur moyenne absolue: \n" + str(mean_absolute_error(y, predict)))


# plot ratio de poisson

def tracage_ratio_poisson():

    plt.subplot(1, 2, 1)
    plt.plot([-1, 1], [-1, 1], '--', color='black')
    plt.plot(y, predict, 'o', color='b', markersize=1)
    #plt.ylabel('$\mu$ $\mathregular{_{GBR}}$ ')
    plt.ylabel('$\mu$ $\mathregular{_{RFR}}$ ')
    plt.xlabel('$\mu$ $\mathregular{_{DFT}}$ ')
    #plt.title('Nu')


def tracage_relative_importance():
    plt.subplot(1, 2, 2)
    plt.barh(pos, feature_importance[sorted_idx], align='center')

    plt.yticks(pos, "", fontsize=8)
    plt.xlabel('Relative Importance')
   # plt.title('Variable Importance')
    plt.savefig("Class_RFR_importance_withElasticProperties_avec_voronoi.pdf", format='pdf', bbox_inches="tight",
                dpi=600)


# Plot feature importance

gbr.fit(X, y)
# make importances relative to max importance
feature_importance = gbr.feature_importances_
feature_importance = 100.0 * (feature_importance / feature_importance.max())
sorted_idx = np.argsort(feature_importance)
print("\nindices des proprietes par ordre d'importance décroissant\n")
print(sorted_idx)
print("\nnom des proprietes par ordre d'importance décroissant\n")
print(str([feature_names[i] for i in sorted_idx]))
pos = np.arange(sorted_idx.shape[0]) + .5
print("\nnombre de propietes dans la liste de proprietes : " + str(len(feature_names)) + " doit être egale à " +
      str(np.shape(X)[1]))
print("\nnombre d'éléments traités dans le tableau : " + str(np.shape(X)[0]))

# =====tracage ratio de poisson
tracage_ratio_poisson()
# =====tracage relative importance
tracage_relative_importance()
