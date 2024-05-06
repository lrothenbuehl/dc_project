rule all:
    input:
        "02 - Data Curation/unified.db"

rule create_database:
    input:
        "01 - Data Extraction/data_pos_neg.csv"
    output:
        "02 - Data Curation/unified.db"
    container:
        "docker://jupyter/datascience-notebook"
    notebook:
        "02 - Data Curation/DB creation.ipynb"