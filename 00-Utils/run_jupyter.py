import requests

def get_kernel():
    # URL for the Jupyter API to create a new kernel
    url = 'http://localhost:8888/api/kernels'

    # Data for the kernel creation
    data = {
        'name': 'python3',  # This should match the kernel name you want, usually a type like 'python3'
        'path': '/'  # The path where the kernel should be started (optional, depends on server setup)
    }

    # Send the POST request to create the kernel
    response = requests.post(url, json=data)

    # Check the response
    if response.status_code == 201:
        print("Kernel created successfully.")
        print(response.json())  # Prints the kernel ID and other details
        return response.json()['id']
    else:
        print(f"Failed to create kernel. Status code: {response.status_code}")
        print(response.text)  # Prints the error message if any

print(get_kernel())