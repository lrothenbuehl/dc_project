from datetime import datetime
from airflow.decorators import dag, task
from airflow.datasets import Dataset
from airflow.providers.papermill.operators.papermill import PapermillOperator

# Define the dataset
data_pos_neg_csv = Dataset('/mnt/01-Data_Extraction/new_data_pos_neg.csv')
unified_db = Dataset('/mnt/02-Data_Curation/unified.db')

@dag(
    dag_id='data_parsing',
    default_args={'retries': 0},
    schedule=[data_pos_neg_csv],  # This ensures the DAG runs when the dataset is modified
    start_date=datetime(2022, 10, 1),
    template_searchpath='/usr/local/airflow/include',
    catchup=False
)
def data_producer_dag():
    
    @task
    def create_SQLite_database():
        return PapermillOperator(
            task_id="create_SQLite_database",
            input_nb="/mnt/02-Data_Curation/DB_creation.ipynb",
            output_nb="/mnt/02-Data_Curation/DB_creation-executed.ipynb",
            parameters={"execution_date": "{{ execution_date }}", "virtual_env": "True"},
        ).execute({})
    
    @task(outlets=[unified_db])
    def DB_curation():
        return PapermillOperator(
            task_id="DB_curation",
            input_nb="/mnt/02-Data_Curation/DB_curation.ipynb",
            output_nb="/mnt/02-Data_Curation/DB_curation-executed.ipynb",
            parameters={"execution_date": "{{ execution_date }}", "virtual_env": "True"},
        ).execute({})
    
    db_creation = create_SQLite_database()
    db_curation = DB_curation()
    
    db_creation >> db_curation

data_producer_dag = data_producer_dag()
