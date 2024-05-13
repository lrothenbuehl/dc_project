import pandas as pd
import sqlite3

# Connect to SQLite database
conn = sqlite3.connect("Z:/projects/dc_project/04-Data_processing/euler/unified_CD.db")

# Load data into a DataFrame
df = pd.read_sql_query("SELECT * FROM prod_desc", conn)

# Drop columns that contain only zeros
df = df.loc[:, (df != 0).any(axis=0)]

# Show the modified DataFrame
print(df)

# If needed, update your table in SQLite
# For example, you can replace the old table with the new DataFrame
df.to_sql('prod_desc', conn, if_exists='replace', index=False)

# Close the connection
conn.close()