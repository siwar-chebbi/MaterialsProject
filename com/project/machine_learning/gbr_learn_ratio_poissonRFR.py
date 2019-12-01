import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_predict, cross_val_score

data = pd.read_csv('Extract_nagative_minNu_maxNu_poisson_descriptors.csv', index_col=0)

print(str(data.columns))

y = data['elasticity.poisson_ratio']

X = data.drop(['elasticity.poisson_ratio', 'cepa', 'minNu', 'maxNu'], axis=1)

print(str(np.shape(X)))
# y = np.log10(y)
lr = linear_model.LinearRegression()

gbr = RandomForestRegressor(n_estimators=100, criterion='mse', min_samples_split=2, min_samples_leaf=3,
                            max_depth=3, max_features='sqrt')

scores = []
iterations = 1
plt.figure(figsize=(12, 6))
predict_result = []
for i in range(0, iterations):
    cv = KFold(n_splits=3, random_state=i)
    predict = cross_val_predict(gbr, X, y, cv=cv, n_jobs=-1)

    score = cross_val_score(gbr, X, y, cv=cv, scoring='neg_mean_squared_error', n_jobs=-1)
    score = (score * -1) ** 0.5
    scores.append(score)
    if i == (iterations - 1):
        predict_result = predict
scores = np.array(scores)
print("Le resultat de prédiction de ration de poisson est :" + str(predict_result))

print(str(np.mean(scores.flatten())))
print(str(np.std(scores.flatten())))
print(str(mean_absolute_error(y, predict)))


# plot ratio de poisson

def tracage_ratio_poisson():
    plt.subplot(1, 2, 1)
    plt.plot([-1, 1], [-1, 1], '--', color='black')
    plt.plot(y, predict, 'o', color='b', markersize=5)
    plt.ylabel('Ratio_Poisson$\mathregular{_{RFR}}$ ')
    plt.xlabel('Ratio_Poisson$\mathregular{_{DFT}}$ ')
    plt.title('Ratio de poisson')


def tracage_relative_importance():
    plt.subplot(1, 2, 2)
    plt.barh(pos, feature_importance[sorted_idx], align='center')
    feature_names = np.array(['minLC', 'maxLC', 'elasticity.K_Voigt_Reuss_Hill',
                              'Emin', 'Emax', 'Gmin', 'Gmax', 'elasticity.G_Reuss',
                              'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill',
                              'elasticity.K_Reuss', 'elasticity.K_Voigt',
                              'lvpa', 'group1', 'atomic_mass1', 'atomicRadius1', 'rowH1A',
                              'rowHn3A', 'xH4A', 'xHn4A', 'cohesive_energy', 'average_electroneg',
                              'bandgap', 'density', 'formation_energy-peratom', 'e_above_hull'])
    print(len(feature_names))
    plt.yticks(pos, feature_names[sorted_idx], fontsize=8)
    plt.xlabel('Relative Importance')
    plt.title('Variable Importance')
    plt.savefig("Class_RFR_importanceRFRwhioutNumaxandmin.pdf", format='pdf', bbox_inches="tight", dpi=600)


# Plot feature importance

gbr.fit(X, y)
# make importances relative to max importance
feature_importance = gbr.feature_importances_
feature_importance = 100.0 * (feature_importance / feature_importance.max())
sorted_idx = np.argsort(feature_importance)
print(sorted_idx)
pos = np.arange(sorted_idx.shape[0]) + .5

print(len(X))

# =====tracage ratio de poisson
tracage_ratio_poisson()
# =====tracage relative importance
tracage_relative_importance()
