import subprocess
from snakemake import checkpoint

# Check if Docker is installed
def check_docker_installed():
    try:
        result = subprocess.run(['docker', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:
            print(f"Docker is installed: {result.stdout.strip()}")
            return True
        else:
            print(f"Docker is not installed. Error: {result.stderr.strip()}")
            return False
    except FileNotFoundError:
        print("Docker is not installed or not found in the system PATH.")
        return False

# Check if Docker engine is running
def check_docker_running():
    try:
        result = subprocess.run(['docker', 'info'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:
            print("Docker engine is running.")
            return True
        else:
            print(f"Docker engine is not running. Error: {result.stderr.strip()}")
            return False
    except FileNotFoundError:
        print("Docker is not installed or not found in the system PATH.")
        return False

# Check if the Docker container exists
def check_container_exists(container_name):
    try:
        result = subprocess.run(['docker', 'images', '-q', container_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.stdout.strip():
            print(f"Docker container '{container_name}' exists.")
            return True
        else:
            print(f"Docker container '{container_name}' does not exist.")
            return False
    except FileNotFoundError:
        print("Docker is not installed or not found in the system PATH.")
        return False

# Define the checkpoint rule
checkpoint check_docker:
    output:
        touch("docker_checked.txt")
    run:
        if not check_docker_installed() or not check_docker_running():
            raise Exception("Docker is not installed or Docker engine is not running. Aborting workflow.")
        else:
            with open(output[0], 'w') as f:
                f.write("Docker is installed and running.")

# Rule that checks if the Docker container exists and builds it if it does not
rule check_and_build_container:
    input:
        checkpoint("check_docker").output[0]
    output:
        touch("container_checked.txt")
    run:
        container_name = "notebook-executor"
        if not check_container_exists(container_name):
            print(f"Building Docker container '{container_name}'...")
            result = subprocess.run(['docker', 'build', '-t', container_name, '-f', './00-Utils/Docker-ipnyb-execution/Dockerfile', '.'],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                raise Exception(f"Failed to build Docker container '{container_name}'. Error: {result.stderr.strip()}")
            print(f"Docker container '{container_name}' built successfully.")
        else:
            print(f"Docker container '{container_name}' already exists.")
        with open(output[0], 'w') as f:
            f.write(f"Docker container '{container_name}' exists or has been built.")

# Rule that generates data
rule generate_data:
    input:
        "docker_checked.txt"
    output:
        "data.txt"
    shell:
        "echo 'Generating data...' > {output}"

# Rule that processes the data
rule process_data:
    input:
        "data.txt"
    output:
        "final_output.txt"
    shell:
        "echo 'Processing data...' > {output}"

# Rule that executes a command if Docker is installed and the container exists
rule execute_command:
    input:
        "container_checked.txt"
    output:
        touch("command_executed.txt")
    shell:
        """
        echo 'Executing command...' > {output}
        # Place your actual command here
        """

# Define the final target
rule all:
    input:
        "final_output.txt",
        "command_executed.txt"
