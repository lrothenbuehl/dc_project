rule all:
    input:
        "02 - Data Curation/unified.db"
        

rule create_database:
    input:
        "01 - Data Extraction/data_pos_neg.csv",
        "02 - Data Curation/DB creation normalized.ipynb"
    output:
        "02 - Data Curation/unified.db"
    container:
        "python:3.11"
    log:
        "logs/create_database.log"
    notebook:
        "02 - Data Curation/DB creation normalized.ipynb"


rule validate_notebook:
    input:
        "02 - Data Curation/DB creation.ipynb"
    output:
        "02 - Data Curation/DB creation normalized.ipynb"
    container:
        "python:3.11"
    log:
        "logs/validate_notebooks.log"
    script:
        "02 - Data Curation/validate notebooks.py"

    