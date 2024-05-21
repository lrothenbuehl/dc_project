from datetime import datetime
from airflow.decorators import dag, task
from airflow.datasets import Dataset
from airflow.providers.papermill.operators.papermill import PapermillOperator

# Define the datasets
unified_db = Dataset('/mnt/02-Data_Curation/unified.db')
report = Dataset('/mnt/03-Data_exploration/DataExploration-completed.ipynb')

@dag(
    dag_id='Exploratory_data_analysis',
    default_args={'retries': 0},
    schedule=[unified_db],  # This ensures the DAG runs when the dataset is modified
    start_date=datetime(2022, 10, 1),
    template_searchpath='/usr/local/airflow/include',
    catchup=False
)
def data_producer_dag():
    
    @task(outlets=[report])
    def calculate_descriptive_statistics():
        return PapermillOperator(
            task_id="calculate_descriptive_statistics",
            input_nb="/mnt/03-Data_exploration/DataExploration.ipynb",
            output_nb="/mnt/03-Data_exploration/DataExploration-completed.ipynb",
            parameters={"execution_date": "{{ execution_date }}", "virtual_env": "True"},
        ).execute({})
    
    calculate_descriptive_statistics()
    

data_producer_dag = data_producer_dag()
