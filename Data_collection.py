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
print(df.head(3)) #Print the first 3 row of the selected targets that have IC50 reported.


# Save the retrieved bioactivity data to a .csv file
df.standard_type.unique()
df.to_csv('bioactivity_data.csv', index = False) 

