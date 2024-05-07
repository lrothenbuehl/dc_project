# Chatgpt code ....
import requests
import json
import sys

def load_notebook(path):
    """ Load a notebook from a local file. """
    with open(path, 'r') as file:
        return json.load(file)

def save_notebook(notebook, path):
    """ Save the executed notebook to a file. """
    with open(path, 'w') as file:
        json.dump(notebook, file, indent=2)

def execute_notebook_on_server(notebook, server_url):
    """ Send the notebook to the Jupyter server for execution and return the executed notebook. """
    headers = {'Content-Type': 'application/json'}
    response = requests.post(server_url, headers=headers, data=json.dumps(notebook))
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f"Failed to execute notebook: {response.text}")

def main(notebook_path, output_path, server_url):
    """ Load, execute, and save a notebook using a remote Jupyter server. """
    # Load the notebook from the provided file path
    notebook = load_notebook(notebook_path)
    
    # Execute the notebook on the specified remote server
    executed_notebook = execute_notebook_on_server(notebook, server_url)
    
    # Save the executed notebook to the specified output path
    save_notebook(executed_notebook, output_path)
    print(f"Notebook executed and saved to {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python execute_notebook.py notebook.ipynb output.ipynb http://remote_jupyter_server/api/kernels")
    else:
        notebook_path = sys.argv[1]
        output_path = sys.argv[2]
        server_url = sys.argv[3]
        main(notebook_path, output_path, server_url)
