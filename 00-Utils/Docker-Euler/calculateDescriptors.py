import os
import pandas as pd
from tqdm import tqdm
import sqlite3
from rdkit import Chem
from rdkit.Chem import Descriptors

# Adjust working directory
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# Open SQLite table
con = sqlite3.connect("unified.db")
cur = con.cursor()

# Import data into pandas DataFrame
df = pd.read_sql_query("SELECT * FROM prod", con)

# Calculate mol objects from sequence
df['mol'] = df['seq'].apply(lambda x: Chem.MolFromSequence(x))

# Function to calculate descriptors
def calculate_descriptors(mol):
    """Calculate all descriptors for a molecule."""
    if mol is not None:
        descriptors = {desc_name: desc_func(mol) for desc_name, desc_func in Descriptors.descList}
        return descriptors
    else:
        return {desc_name: None for desc_name, desc_func in Descriptors.descList}

# Function to process DataFrame with an option to limit the number of rows for testing
def process_dataframe(df, limit=None):
    # Initialize tqdm within the apply function to show progress
    tqdm.pandas(desc="Calculating Descriptors")

    # If limit is specified, reduce the dataframe to the first 'limit' rows
    if limit is not None:
        df = df.head(limit)

    # Applying the function and converting results to DataFrame
    df_descriptors = df['mol'].progress_apply(calculate_descriptors).apply(pd.Series)

    # Concatenate the original DataFrame with the new descriptors DataFrame
    df_final = pd.concat([df, df_descriptors], axis=1)

    return df_final

# Process a limited number of rows for testing
df_final = process_dataframe(df, limit=10)  # Modify '10' to whatever number of rows you want to test

# Save the DataFrame to a new SQLite table
df_final.drop('mol', axis=1, inplace=True)  # Drop the mol column as it's not needed in the database
df_final.to_sql('prod_desc', con, if_exists='replace', index=False)

# Output a sample of the final DataFrame
print(df_final.head())

con.commit()
con.close()