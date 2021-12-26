# Install lazypredict libraries: "pip install lazypredict"

# Import libraries:
import pandas as pd
import seaborn as sns
from sklearn.model_selection import train_test_split
from lazypredict.Supervised import LazyClassifier
from lazypredict.Supervised import LazyRegressor
from sklearn.feature_selection import VarianceThreshold
import matplotlib.pyplot as plt
import seaborn as sns


# Load the dataset:
df = pd.read_csv('VCP_00_bioactivity_data_3class_pIC50_pubchem_fp.csv')

# Get X and Y variables:
X = df.drop('pIC50', axis=1)
Y = df.pIC50

# Remove low variance features
selection = VarianceThreshold(threshold=(.8 * (1 - .8)))    
X = selection.fit_transform(X)

# Perform data splitting using 80/20 ratio
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

# Compare ML algorithms:

## Defines and builds the lazyclassifier:
clf = LazyRegressor(verbose=0, ignore_warnings=True, custom_metric=None)
models_train,prediction_train = clf.fit(X_train, X_train, Y_train, Y_train)
models_test,predictions_test = clf.fit(X_train, X_test, Y_train, Y_test)

## Performance table of the training set (80# subset)
print(prediction_train)

print(predictions_test)

#train["R-Squared"] = [0 if i < 0 else i for i in train.iloc[:,0] ]
plt.figure(figsize=(25, 10))
sns.set_theme(style="whitegrid")
ax = sns.barplot(y=prediction_train.index, x="R-Squared", data=prediction_train)
ax.set(xlim=(0, 1))
plt.savefig('R_squared_of_regression_models.pdf')

# Bar plot of RMSE values
plt.figure(figsize=(25, 10))
sns.set_theme(style="whitegrid")
ax = sns.barplot(y=prediction_train.index, x="RMSE", data=prediction_train)
ax.set(xlim=(0, 10))
plt.savefig('RMSE_of_regression_models.pdf')