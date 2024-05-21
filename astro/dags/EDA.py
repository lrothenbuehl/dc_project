from datetime import datetime
from airflow.decorators import dag, task
from airflow.datasets import Dataset
from airflow.providers.papermill.operators.papermill import PapermillOperator

# Define the dataset
unified_db = Dataset('/mnt/02-Data_Curation/unified.db')
EDA_report = Dataset('/mnt/03-Data_exploration/DataExploration-executed.ipynb')

@dag(
    dag_id='EDA',
    default_args={'retries': 0},
    schedule=[unified_db],  # This ensures the DAG runs when the dataset is modified
    start_date=datetime(2022, 10, 1),
    template_searchpath='/usr/local/airflow/include',
    catchup=False
)
def data_producer_dag():
    
    @task(outlets=[EDA_report])
    def run_EDA():
        return PapermillOperator(
            task_id="run_EDA",
            input_nb="/mnt/03-Data_exploration/DataExploration.ipynb",
            output_nb="/mnt/03-Data_exploration/DataExploration-executed.ipynb",
            parameters={"execution_date": "{{ execution_date }}", "virtual_env": "True"},
        ).execute({})
    
    run_EDA()

data_producer_dag = data_producer_dag()
