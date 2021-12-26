# Install lazypredict libraries: "pip install lazypredict"

# Import libraries:
import pandas as pd
import seaborn as sns
from sklearn.model_selection import train_test_split
import lazypredict
from lazypredict.Supervised import LazyClassifier
from lazypredict.Supervised import LazyRegressor
from sklearn.feature_selection import VarianceThreshold



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