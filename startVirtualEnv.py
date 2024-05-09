import subprocess
import os
import webbrowser
import time

cwd = os.getcwd().replace("\\","/")


def run_docker_jupyter():
    current_directory = cwd  # Get the current working directory
    docker_command = f"docker run -it -v {current_directory}:/app -d -p 8888:8888 jupyter-server"
    try:
        subprocess.run(docker_command, shell=True, check=True)
        print("Jupyter server is running.")
    except subprocess.CalledProcessError:
        print("Failed to start the Jupyter server.")

def check_docker():
    try:
        # Execute the 'docker info' command
        subprocess.run(['docker', 'info'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError:
        return False

def is_docker_image_built(image_name):
    try:
        # Run 'docker images' command
        result = subprocess.run(['docker', 'images', '--format', '{{.Repository}}'],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # Get the output and check if the image name is in there
        images = result.stdout.strip().split('\n')
        return image_name in images
    except subprocess.CalledProcessError:
        # Return False if the command fails to execute properly
        return False

    
def set_working_directory_to_script_location():
    # Get the absolute path of the script file
    script_path = os.path.abspath(__file__)
    
    # Extract the directory part from the script path
    directory = os.path.dirname(script_path)
    
    # Change the current working directory to the directory of the script
    os.chdir(directory)
    print(f"Current working directory has been set to: {os.getcwd()}")

def build_docker_image():
    try:
        # Define the Docker build command with the specific tag and Dockerfile location
        command = ['docker', 'build', '-t', 'jupyter-server', '-f', '00-Utils/Docker-jupyter-server/Dockerfile', '.']
        
        # Execute the Docker build command and wait for it to complete
        result = subprocess.run(command, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Output the results
        print("Docker image built successfully:")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        # If the command fails, output the error
        print("Failed to build Docker image:")
        print(e.stderr)

def stop_all_docker_containers():
    try:
        # Get the list of all running container IDs
        ps_command = ['docker', 'ps', '-a', '-q']
        ps_result = subprocess.run(ps_command, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        container_ids = ps_result.stdout.strip().split('\n')
        
        # If there are no containers to stop, just return
        if not container_ids or container_ids == ['']:
            print("No Docker containers are running.")
            return

        # Stop all listed containers
        stop_command = ['docker', 'stop'] + container_ids
        stop_result = subprocess.run(stop_command, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Output the results
        print("Stopped all Docker containers:")
    except subprocess.CalledProcessError as e:
        # If any command fails, output the error
        print("Failed to stop Docker containers:")
        print(e.stderr)

def open_url_in_browser(url):
    """Opens the given URL in the default web browser."""
    time.sleep(1) # give docker and server time to start
    try:
        # Attempt to open the URL in a new window, if possible
        webbrowser.open_new(url)
        print("URL successfully opened in browser.")
    except Exception as e:
        print(f"Failed to open URL: {e}")

set_working_directory_to_script_location()
stop_all_docker_containers()

if(check_docker()):
    print("Docker Engine is running. Checking if jupyter-server is already built.")
    if(is_docker_image_built("jupyter-server")):
        print("jupyter-server is already built. Executing jupyter-server")
        run_docker_jupyter()
        open_url_in_browser('http://127.0.0.1:8888/tree?')
    else:
        print("Docker image has not been built jet. Building...")
        build_docker_image()
        run_docker_jupyter()
        open_url_in_browser('http://127.0.0.1:8888/tree?')

else:
    print("Docker Engine is not running please start docker first")
    exit(0)
