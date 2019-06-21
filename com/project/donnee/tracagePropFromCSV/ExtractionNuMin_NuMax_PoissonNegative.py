import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.backends.backend_pdf
import matplotlib.cm as cm
import matplotlib.colors as color
import math

#propsTableau_OLD = ["minLC", "maxLC", "minNu", "maxNu", "K_Voigt_Reuss_Hill", "Emin", "Emax", "Gmin", "Gmax"]

propsTableau = ["minLC", "maxLC", "minNu", "maxNu", "Emin", "Emax", "Gmin", "Gmax",
                'elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill',
                'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']


def importer(fichier):
    return pd.read_csv(fichier, index_col=0)


data = importer("elasticElate_ALL_revisionArt_without_Zero.csv")

data.head()
print("\nNombre de tous les éléments dans le fichier csv = {}\n".format(data.shape[0]))

##extract_data1 = data[(data['minNu'] < 0)]
##extract_data1.to_csv("Extract_nagative_minNu_poisson.csv")
##print("Nombre d'éléments extraits avec minNu < 0 = {}\n".format(extract_data1.shape[0]))

#extract_data2 = data[(data['minLC'] < 0)]
#extract_data2.to_csv("Extract_nagative_minLCEXP.csv")
#print("Nombre d'éléments extraits avec  minLC < 0 = {}\n".format(extract_data2.shape[0]))

#extract_data3 = data[(data['elasticity.poisson_ratio'] < 0)]
#extract_data3.to_csv("Extract_nagative_poisson_ratio.csv")
#print("Nombre d'éléments extraits avec negative poisson ratio < 0 = {}\n".format(extract_data3.shape[0]))

#extract_data4 = data[(data['minNu'] < 0) & (data['maxNu'] < 0) & (data['elasticity.poisson_ratio'] < 0)]
#extract_data4.to_csv("Extract_nagative_minNu_maxNu_poisson.csv")
#print("Nombre d'éléments extraits avec minNu < 0, maxNu < 0 poisson_ratio < 0 = {}\n".format(extract_data4.shape[0]))

#extract_data5 = data[(data['minNu'] < 0) & (data['maxNu'] > 0) & (data['elasticity.poisson_ratio'] < 0)]
#extract_data5.to_csv("Extract_nagative_minNu_poisson.csv")
#print("Nombre d'éléments extraits avec minNu < 0 and poisson_ratio < 0 = {}\n".format(extract_data5.shape[0]))

#extract_data6 = data[(data['minNu'] < 0) & (data['maxNu'] < 0)]
#extract_data6.to_csv("Extract_nagative_minAndmaxNu.csv")
#print("Nombre d'éléments extraits avec minAndmaxNu < 0 = {}\n".format(extract_data6.shape[0]))

extract_data7 = data[(data['minLC'] < 0)]
extract_data7.to_csv("Extract_nagative_minAndmaxLC.csv")
print("Nombre d'éléments extraits avec minAndmaxLC < 0 = {}\n".format(extract_data7.shape[0]))
