# Create a conda environment with rdkit using: "$ conda create -c conda-forge -n data_explore_env rdkit"
# Activate the created environment: "$ conda activate data_explore_env"

# Import libraries:
import pandas as pd
import numpy as np
from pandas.core.frame import DataFrame
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import seaborn as sns
sns.set(style='ticks')
import matplotlib.pyplot as plt

# Load bioactivity data:
df = pd.read_csv('bioactivity_preprocessed_data.csv')

# Calculate descriptors:

## The lipinski(smiles) calculate the Lipinski indicators based on the smiles notation
## of the compounds:
def lipinski(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(str(elem))
        moldata.append(mol)

    baseData = np.arange(1,1)
    i = 0
    for mol in moldata:

        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)

        row = np.array([desc_MolWt, desc_MolLogP, desc_NumHDonors, desc_NumHAcceptors])

        if i == 0:
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i + 1
    
    columnNames = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']
    descriptors = pd.DataFrame(data = baseData, columns = columnNames)

    return descriptors

## Applied function above to our dataframe
df_lipinski = lipinski(df.canonical_smiles)
print(df_lipinski) #Print out Lipinski of our chemicals
print(df) #Print out preprocessed data from bioactivity_preprocessed_data.csv

# Combine 2 dataframe above:
df_combine = pd.concat([df, df_lipinski], axis = 1)
print(df_combine) #Print out combined data

# Convert IC50 to pIC50 (or -log10(IC50)):

## Function pIC50(input) convert IC50 value in the dataframe to pIC50:
def pIC50(input):
    pIC50 = []

    for i in input['standard_value']:
        molar = i*(10**-9) #convert nM to M
        pIC50.append(-np.log10(molar))
    
    input['pIC50'] = pIC50
    x = input.drop('standard_value', 1)

    return x

## Applied the function above to our combined dataframe
df_final = pIC50(df_combine)
print(df_final.pIC50.describe()) #Check the value of pIC50 in final dataframe
print(df_final.MW.describe())

# Save the final dataframe to a .csv file:
df_final.to_csv('VCP_00_bioactivity_data_3class_pIC50.csv')

# Removing "intermediate" bioactivity class:
df_2class = df_final[df_final['bioactivity_class'] != 'intermediate']
print(df_2class) #print the new dataframe

# Exploritory data analysis:

## Frequency plot of the 2 bioactivity classes:
plt.figure(figsize=(5.5, 5.5))

sns.countplot(x = 'bioactivity_class', data = df_2class, edgecolor = 'black')

plt.xlabel('Bioactivity class', fontsize = 14, fontweight = 'bold')
plt.ylabel('Frequency', fontsize = 14, fontweight = 'bold')

plt.savefig('plot_bioactivity_class.pdf')

## Scatter plot of MW vs LogP:
plt.figure(figsize=(5.5,5.5))

sns.scatterplot(x='MW', y='LogP', hue = 'bioactivity_class', size = 'pIC50', edgecolor = 'black', alpha = 0.7, data = df_2class)

plt.xlabel('MW', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')
plt.legend(loc = 0)
plt.savefig('plot_MW_vs_LogP.pdf')

# Box plots:

## pIC50 values:
plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'pIC50', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('pIC50 value', fontsize=14, fontweight='bold')

plt.savefig('plot_pIC50.pdf')

## Stat analysis of pIC50 values
def mannwhitney(descriptor, verbose=False):
  #https://github.com/sonluongvu/code/blob/master/python/CDD_ML_Part_2_Acetylcholinesterase_Exploratory_Data_Analysis.ipynb
  from numpy.random import seed
  from numpy.random import randn
  from scipy.stats import mannwhitneyu

# seed the random number generator
  seed(1)

# actives and inactives
  selection = [descriptor, 'bioactivity_class']
  df = df_2class[selection]
  active = df[df['bioactivity_class'] == 'active']
  active = active[descriptor]

  selection = [descriptor, 'bioactivity_class']
  df = df_2class[selection]
  inactive = df[df['bioactivity_class'] == 'inactive']
  inactive = inactive[descriptor]

# compare samples
  stat, p = mannwhitneyu(active, inactive)
  #print('Statistics=%.3f, p=%.3f' % (stat, p))

# interpret
  alpha = 0.05
  if p > alpha:
    interpretation = 'Same distribution (fail to reject H0)'
  else:
    interpretation = 'Different distribution (reject H0)'
  
  results = pd.DataFrame({'Descriptor':descriptor,
                          'Statistics':stat,
                          'p':p,
                          'alpha':alpha,
                          'Interpretation':interpretation}, index=[0])
  filename = 'mannwhitneyu_' + descriptor + '.csv'
  results.to_csv(filename)

  return results

print(mannwhitney('pIC50')) #Print out stat analysis for distribution of MW in 2 bioactivity class

## MW plots:
plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'MW', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('MW', fontsize=14, fontweight='bold')

plt.savefig('plot_MW.pdf')

## Stat analysis of MW:
print(mannwhitney('MW'))

## LogP plots:
plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'LogP', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')

plt.savefig('plot_LogP.pdf')

## Stat analysis of LogP:
print(mannwhitney('LogP'))

## NumHDonnors plots:
plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'NumHDonors', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('NumHDonors', fontsize=14, fontweight='bold')

plt.savefig('plot_NumHDonors.pdf')

## Stat analysis of NumHDonnors:
print(mannwhitney('NumHDonors'))

## NumHAcceptors plots:
plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'NumHAcceptors', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('NumHAcceptors', fontsize=14, fontweight='bold')

plt.savefig('plot_NumHAcceptors.pdf')

## Stat analysis of NumHAcceptors:
print(mannwhitney('NumHAcceptors'))


# All the .csv and .pdf is saved into results.zip by using "zip -r results.zip . -i *.csv *.pdf"