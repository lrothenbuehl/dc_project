import os
import sqlite3
import pandas as pd
import numpy as np
import multiprocessing as mp
from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_all_descriptors(mol):
    descriptor_names = [desc_name for desc_name, _ in Descriptors._descList]
    descriptors = {}
    for desc_name in descriptor_names:
        try:
            descriptor_value = getattr(Descriptors, desc_name)(mol)
            descriptors[desc_name] = descriptor_value
        except Exception as e:
            descriptors[desc_name] = None
            print(f"Error calculating descriptor {desc_name}: {e}")
    return descriptors

def calculate_descriptors_chunk(chunk, chunk_id):
    print(f"Starting chunk {chunk_id}")
    result = []
    for _, row in chunk.iterrows():
        mol = Chem.MolFromSequence(row['seq'])
        if mol is not None:
            descriptors = calculate_all_descriptors(mol)
            descriptors['seq'] = row['seq']
            descriptors['seq_length'] = len(row['seq'])
            result.append(descriptors)
        else:
            print(f"Invalid sequence: {row['seq']}")
            result.append({**{'seq': row['seq'], 'seq_length': len(row['seq'])},
                           **{desc_name: None for desc_name, _ in Descriptors._descList}})
    print(f"Finished chunk {chunk_id}")
    return pd.DataFrame(result)

def calculate_descriptors(df, num_threads):
    # Split the dataframe into chunks
    chunks = np.array_split(df, num_threads)

    # Create a pool of workers
    pool = mp.Pool(processes=num_threads)
    tasks = [(chunk, i) for i, chunk in enumerate(chunks)]
    result_chunks = pool.starmap(calculate_descriptors_chunk, tasks)

    # Wait for all workers to complete
    pool.close()
    pool.join()

    # Concatenate the results
    result_df = pd.concat(result_chunks, ignore_index=True)
    return result_df

def save_to_sqlite(df, db_name='unified.db', table_name='DESC'):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # Create table if it does not exist
    create_table_query = f"""
    CREATE TABLE IF NOT EXISTS {table_name} (
        id INTEGER PRIMARY KEY,
        seq TEXT,
        seq_length INTEGER,
        {', '.join([f"{col} REAL" for col in df.columns if col not in ['seq', 'seq_length']])}
    )
    """
    cursor.execute(create_table_query)

    # Insert DataFrame into SQLite table
    df.to_sql(table_name, conn, if_exists='append', index=False)

    # Commit changes and close the connection
    conn.commit()
    conn.close()

if __name__ == "__main__":
    # Get number of cores
    num_cores = os.cpu_count()
    print(f"CPU cores: {num_cores}")

    # DB import
    con = sqlite3.connect('unified.db')
    df = pd.read_sql_query("SELECT * FROM prod", con)
    con.close()
    print(f"Found {len(df['seq'])} sequences.")  

    # Calculate descriptors
    result_df = calculate_descriptors(df, num_threads=num_cores - 1)  # Leave room for main thread

    # Save the DataFrame to SQLite
    save_to_sqlite(result_df)
