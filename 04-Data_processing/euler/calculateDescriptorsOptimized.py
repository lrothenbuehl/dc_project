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

con = sqlite3.connect("unified.db")

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

# Function to process a chunk of the DataFrame
def process_chunk(chunk):
    tqdm.pandas(desc="Calculating Descriptors")
    # Applying the function and converting results to DataFrame
    descriptors_df = chunk['mol'].progress_apply(calculate_descriptors).apply(pd.Series)
    return pd.concat([chunk.drop(columns='mol'), descriptors_df], axis=1)

# Process the data in batches and save to SQLite with a limit for testing
def process_data(df, batch_size=500, limit=None):
    if limit is not None:
        df = df.head(limit)  # Limit the dataframe to the specified number of rows for testing

    num_batches = (len(df) + batch_size - 1) // batch_size  # Calculates the required number of batches
    tqdm_iterator = tqdm(range(num_batches), desc="Processing Batches")

    for i in tqdm_iterator:
        start_idx = i * batch_size
        end_idx = min((i + 1) * batch_size, len(df))
        chunk = df.iloc[start_idx:end_idx]
        processed_chunk = process_chunk(chunk)
        processed_chunk.to_sql('prod_desc', con, if_exists='append', index=False)  # Append each chunk to the table

    con.commit()
    con.close()

# Specify the limit for testing here
process_data(df, batch_size=10, limit=50)

print("Processing complete. Data has been saved to the database.")
