import os
import sqlite3
import pandas as pd
import numpy as np
import multiprocessing as mp
from rdkit import Chem
from rdkit.Chem import Descriptors
import time
from threading import Thread

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

def save_to_sqlite(df, db_name='unified.db', table_name='DESC'):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # Create table if it does not exist
    columns = ', '.join([f"{col} REAL" for col in df.columns if col not in ['seq', 'seq_length']])
    create_table_query = f"""
    CREATE TABLE IF NOT EXISTS {table_name} (
        id INTEGER PRIMARY KEY,
        seq TEXT,
        seq_length INTEGER,
        {columns}
    )
    """
    print(create_table_query)  # Debug print to verify SQL query
    cursor.execute(create_table_query)

    # Insert DataFrame into SQLite table
    df.to_sql(table_name, conn, if_exists='append', index=False)

    # Commit changes and close the connection
    conn.commit()
    conn.close()

def calculate_descriptors_chunk(chunk, thread_id, running_threads, progress_dict):
    with running_threads.get_lock():
        running_threads.value += 1
        print(f"Thread {thread_id} started. Running threads: {running_threads.value}")

    result = []
    for idx, (_, row) in enumerate(chunk.iterrows()):
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

        # Update progress
        progress_dict[thread_id] = idx + 1

    result_df = pd.DataFrame(result)
    save_to_sqlite(result_df)  # Save the result of this chunk

    with running_threads.get_lock():
        running_threads.value -= 1
        print(f"Thread {thread_id} finished. Running threads: {running_threads.value}")

def calculate_descriptors(df, num_threads):
    # Split the dataframe into chunks
    chunks = np.array_split(df, num_threads)

    # Create a multiprocessing manager
    manager = mp.Manager()
    running_threads = manager.Value('i', 0)
    progress_dict = manager.dict({i: 0 for i in range(num_threads)})

    # Function to print progress every 10 seconds
    def print_progress():
        while running_threads.value > 0:
            time.sleep(10)
            for thread_id in range(num_threads):
                print(f"Thread {thread_id}: {progress_dict[thread_id]}/{len(chunks[thread_id])} rows processed")

    # Start the progress printing thread
    progress_thread = Thread(target=print_progress)
    progress_thread.start()

    # Create a pool of workers
    pool = mp.Pool(processes=num_threads)
    tasks = [(chunk, thread_id, running_threads, progress_dict) for thread_id, chunk in enumerate(chunks)]

    # Start the pool
    pool.starmap(calculate_descriptors_chunk, tasks)

    # Wait for all workers to complete
    pool.close()
    pool.join()

    # Stop the progress printing thread
    progress_thread.join()

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
    calculate_descriptors(df, num_threads=num_cores - 1)  # Leave room for main thread
