# Import libraries:
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import VarianceThreshold

# Read data:
df = pd.read_csv('VCP_00_bioactivity_data_3class_pIC50_pubchem_fp.csv')

# Input features:
X = df.drop('pIC50', axis = 1)
print(X)
Y = df.pIC50
print(Y)

## Examine data demension:
print(X.shape)
print(Y.shape)

## Remove low varience features:
selection = VarianceThreshold(threshold=(.8 * (1 - .8)))    
X = selection.fit_transform(X)

# Data split (80/20 ratio):
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2)
print(X_train.shape, Y_train.shape)
print(X_test.shape, Y_test.shape)

# Building a Regression Model using Random Forest
model = RandomForestRegressor(n_estimators=1000)
model.fit(X_train, Y_train)
r2 = model.score(X_test,Y_test)
print(r2)

Y_pred = model.predict(X_test)

# Scatter plot of experimental vs predicted pIC50 values:
sns.set(color_codes=True)
sns.set_style("white")

ax = sns.regplot(Y_test, Y_pred, scatter_kws={'alpha':0.4})
ax.set_xlabel('Experimental pIC50', fontsize='large', fontweight='bold')
ax.set_ylabel('Predicted pIC50', fontsize='large', fontweight='bold')
ax.set_xlim(0, 12)
ax.set_ylim(0, 12)
ax.figure.set_size_inches(5, 5)
plt.savefig('Experimental_vs_predicted_pIC50.pdf')