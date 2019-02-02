import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.backends.backend_pdf
import matplotlib.cm as cm
import matplotlib.colors as color
import math

propsTableau = ["minLC", "maxLC", "minNu", "maxNu", "K_Voigt_Reuss_Hill", "Emin", "Emax", "Gmin", "Gmax"]


def importer(fichier):
    return pd.read_csv(fichier, index_col=0)


data = importer("elasticElateALL.csv")

data.head()
print("\nNombre de tous les éléments dans le fichier csv = {}\n".format(data.shape[0]))

extract_data1 = data[(data['minNu'] < 0) | (data['elasticity.poisson_ratio'] < 0)]
extract_data1.to_csv("Extract_nagative_minNu_poisson.csv")
print("Nombre d'éléments extraits avec minNu < 0 ou  poisson_ratio < 0 = {}\n".format(extract_data1.shape[0]))

extract_data2 = data[(data['minLC'] < 0)]
extract_data2.to_csv("Extract_nagative_minLC.csv")
print("Nombre d'éléments extraits avec  minLC < 0 = {}\n".format(extract_data2.shape[0]))
