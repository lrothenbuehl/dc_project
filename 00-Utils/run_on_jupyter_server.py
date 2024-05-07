import requests
import json
import sys
import time

def load_notebook(path):
    """ Load a notebook from a local file. """
    with open(path, 'r') as file:
        return json.load(file)

def save_notebook(notebook, path):
    """ Save the notebook to a local file after execution. """
    with open(path, 'w') as file:
        json.dump(notebook, file, indent=2)

def create_kernel(server_url):
    """ Create a new kernel and return its ID. """
    response = requests.post(f"{server_url}/api/kernels", json={"name": "python3"})
    response.raise_for_status()  # Will raise an HTTPError if the HTTP request returned an unsuccessful status code
    return response.json()['id']

def execute_code(kernel_id, server_url, code):
    """ Execute code in the specified kernel and return the output. """
    data = {
        "code": code,
        "silent": False,
        "store_history": False,
        "user_expressions": {},
        "allow_stdin": False
    }
    headers = {'Content-Type': 'application/json'}
    response = requests.post(f"{server_url}/api/kernels/{kernel_id}/execute", headers=headers, json=data)
    response.raise_for_status()
    return response.json()

def main(notebook_path, output_path, server_url):
    """ Load, execute, and save a notebook. """
    notebook = load_notebook(notebook_path)
    kernel_id = create_kernel(server_url)

    try:
        for cell in notebook['cells']:
            if cell['cell_type'] == 'code':
                print(f"Executing cell {cell['execution_count']}...")
                output = execute_code(kernel_id, server_url, ''.join(cell['source']))
                cell['outputs'] = [output]
                cell['execution_count'] = output.get('execution_count', None)
                time.sleep(1)  # sleep to avoid rate limiting if necessary
    finally:
        # Always attempt to delete the kernel
        requests.delete(f"{server_url}/api/kernels/{kernel_id}")

    save_notebook(notebook, output_path)
    print(f"Notebook executed and saved to {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py notebook.ipynb output.ipynb http://127.0.0.1:8888")
    else:
        main(*sys.argv[1:])
