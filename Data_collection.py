# ChEMBL webservice has been installed via "python3 -m pip install chembl_webresource_client" in Terminal
# importing library:
import pandas as pd
from chembl_webresource_client.new_client import new_client

# Search for target protein:
## Target search for SVIP protein:

target = new_client.target
target_query = target.search('VCP')
targets = pd.DataFrame.from_dict(target_query)
print(targets) #Print the target search

# Select and retrieve bioactivity data for "Transitional endoplasmic reticulum ATPase"

selected_target = targets.target_chembl_id[0]
print(selected_target) #Print the selected target

activity = new_client.activity
res = activity.filter(target_chembl_id = selected_target).filter(standard_type = "IC50")
df = pd.DataFrame.from_dict(res)
print(df) #Print the selected targets that have IC50 reported.


# Save the retrieved bioactivity data to a .csv file
df.standard_type.unique()
df.to_csv('bioactivity_data.csv', index = False) 

# Handling missing data:

df2 = df[df.standard_value.notna()].reset_index()
print(df2) #Print out data without including missing data

# Data preprocessing:

## Labelling compound to be active (IC50 <= 1000 nM), intermediate (1000 nM < IC50 < 10000 nM) and inactive (IC50 >= 10000):

bioactivity_class = []
for i in df2.standard_value:
    if float(i) <= 1000.0:
        bioactivity_class.append('active')
    elif float(i) <= 10000.0:
        bioactivity_class.append('intermediate')
    else:
        bioactivity_class.append('inactive')
        
## Create a dataframe contains molecule_chembl_id, canonical_smiles, standard_value and bioactivity_class:

selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
df3 = df2[selection].reset_index()
print(df3)
df3 = pd.concat([df3, pd.Series(bioactivity_class, name='bioactivity_class')], axis=1) # put the bioactivity class to the dataframe
print(df3) #Print the created dataframe



df3.to_csv('bioactivity_preprocessed_data.csv', index=False) #Save the created dataframe to a .csv file
