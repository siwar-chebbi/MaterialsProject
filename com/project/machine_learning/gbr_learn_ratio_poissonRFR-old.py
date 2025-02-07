import sys
import json
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from pprint import pprint
from sklearn.preprocessing import RobustScaler
import pymatgen as pmg
from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.core.structure import Structure
import scipy.stats as stats
from pymatgen.io.zeopp import ZeoCssr, ZeoVoronoiXYZ, get_voronoi_nodes, \
    get_high_accuracy_voronoi_nodes, get_void_volume_surfarea, \
    get_free_sphere_params
from sklearn.ensemble import AdaBoostRegressor, RandomForestRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import mean_squared_error
import matplotlib as mpl

data = pd.read_csv('Extract_nagative_minNu_maxNu_poisson_descriptors.csv', index_col=0)

print(str(data.columns))

y = data['elasticity.poisson_ratio']

X = data.drop(['elasticity.poisson_ratio', 'cepa'], axis=1)

print(str(np.shape(X)))

# y = np.log10(y)

from sklearn.model_selection import cross_val_predict, cross_val_score
from sklearn import linear_model

lr = linear_model.LinearRegression()
import seaborn as sns
from sklearn.model_selection import KFold

gbr = RandomForestRegressor(n_estimators=100, criterion='mse', min_samples_split=2, min_samples_leaf=3,
                                max_depth=3, max_features='sqrt')

gbr.fit(X, y)
# gbr = GradientBoostingRegressor(n_estimators=1000, learning_rate=0.01)
#
# gbr = AdaBoostRegressor(DecisionTreeRegressor(max_depth=4),
#                         n_estimators=300, random_state=0)

mpl.rcParams['pdf.fonttype'] = 42
sns.set(style="white", palette="muted")
sns.set_style("ticks")
sns.set_context("notebook", font_scale=1.8, rc={"figure.figsize": (8, 6)})
sns.set_palette("muted")
b, gr, r, p, v = sns.color_palette("muted", 5)

scores = []
iterations = 100
for i in range(0, iterations):
    cv = KFold(n_splits=3, random_state=i)
    predict = cross_val_predict(gbr, X, y, cv=cv, n_jobs=-1)
    score = cross_val_score(gbr, X, y, cv=cv, scoring='neg_mean_squared_error', n_jobs=-1)
    score = (score * -1) ** 0.5
    scores.append(score)
    if i == (iterations - 1):
        plt.plot([0, 1], [0, 1], '--', color='black')
        plt.plot(y, predict, 'o', color=b, markersize=10)

scores = np.array(scores)
# print(str(scores
print(str(np.mean(scores.flatten())))
print(str(np.std(scores.flatten())))

plt.subplot(1, 2, 1)
plt.ylabel('Ratio_Poisson$\mathregular{_{RFR}}$ ')
plt.xlabel('Ratio_Poisson$\mathregular{_{DFT}}$ ')

feature_importance = gbr.feature_importances_
feature_importance = 100.0 * (feature_importance / feature_importance.max())
sorted_idx = np.argsort(feature_importance)
pos = np.arange(sorted_idx.shape[0]) + .5
plt.subplot(1, 2, 2)
plt.barh(pos, feature_importance[sorted_idx], align='center')
plt.yticks(pos, str(np.shape(X)))
plt.xlabel('Relative Importance')
plt.title('Variable Importance')
#plt.show()
plt.savefig("Ratio_poisson_100_iterationSRFR.pdf", format='pdf', bbox_inches="tight", dpi=600)
