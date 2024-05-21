from datetime import datetime
from airflow.decorators import dag
from airflow.models.baseoperator import chain
from airflow.datasets import Dataset
from airflow.providers.papermill.operators.papermill import PapermillOperator

# Define the datasets
unified_db = Dataset('/mnt/02-Data_Curation/unified.db')
EDA_report = Dataset('/mnt/03-EDA/DataExploration-executed.ipynb')

@dag(
    dag_id='EDA',
    default_args={'retries': 0},
    schedule=[unified_db],  # This ensures the DAG runs when the dataset is modified
    start_date=datetime(2022, 10, 1),
    template_searchpath='/usr/local/airflow/include',
    catchup=False
)
def EDA_dag():

    run_EDA = PapermillOperator(
        task_id="run_EDA",
        input_nb="/mnt/03-EDA/DataExploration.ipynb",
        output_nb="/mnt/03-EDA/DataExploration-executed.ipynb",
        parameters={"execution_date": "{{ execution_date }}", "virtual_env":"True"},
        outlets=[EDA_report]
    )
    
    chain(run_EDA)

data_producer_dag = EDA_dag()
