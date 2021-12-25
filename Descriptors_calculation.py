# Install PaDELPy-descriptor package: 'pip install padelpy'
# PaDEL descriptor is download from: 'http://www.yapcwsoft.com/dd/padeldescriptor/'


# Import libraries:
import pandas as pd

# Load bioactivity data:
df3 = pd.read_csv('VCP_00_bioactivity_data_3class_pIC50.csv')
print(df3)

# Retrieve canonical smiles and chEMBL id:
selection = ['canonical_smiles', 'molecule_chembl_id']
df3_selection = df3[selection]
df3_selection.to_csv('molecule.smi', sep='\t', index = False, header = False)

