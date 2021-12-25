# Install PaDELPy-descriptor package: 'pip install padelpy'
# PaDEL descriptor is download from: 'http://www.yapcwsoft.com/dd/padeldescriptor/'
# The .zip is extract to PaDEL-descriptor folder

# Import libraries:
import pandas as pd

# Load bioactivity data:
df3 = pd.read_csv('VCP_00_bioactivity_data_3class_pIC50.csv')
print(df3)

# Retrieve canonical smiles and chEMBL id:
selection = ['canonical_smiles', 'molecule_chembl_id']
df3_selection = df3[selection]
df3_selection.to_csv('molecule.smi', sep='\t', index = False, header = False)

# Run in terminal:"java -Xms1G -Xmx1G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro 
# -fingerprints -descriptortypes ./PaDEL-Descriptor/Descriptors.xml -dir ./ -file descriptors_output.csv"
# Put the output file to df3_X dataframe
df3_X = pd.read_csv('descriptors_output.csv')
print(df3_X)

# X data matrix
df3_X = df3_X.drop(columns='Name')
print(df3_X)

# Y variable:
df3_Y = df3['pIC50']
print(df3_Y)

# Combine X and Y:
dataset3 = pd.concat([df3_X,df3_Y], axis = 1)
print(dataset3)

# Save to .csv file:
dataset3.to_csv('VCP_00_bioactivity_data_3class_pIC50_pubchem_fp.csv', index=False)
