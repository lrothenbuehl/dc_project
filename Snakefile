import os

cwd = os.getcwd().replace("\\","/")

rule all:
    input: 
        "00-Utils/Docker-images/notebook-executor.tar",
        "00-Utils/Docker-images/jupyter-server.tar",
        "serverRunning"


rule create_docker_image_notebook_executor:
    input:
        "00-Utils/Docker-ipnyb-execution/Dockerfile",
        "00-Utils/Docker-ipnyb-execution/requirements.txt",
        "00-Utils/Docker-ipnyb-execution/execute_notebook.sh"
    output:
        "00-Utils/Docker-images/notebook-executor.tar"
    log:
        "logs/setup_notebook-executor.log"
    run:
        shell("echo 'Building Docker image...'")
        shell("docker build -t notebook-executor -f 00-Utils/Docker-ipnyb-execution/Dockerfile .")
        shell("echo 'Saving Docker image to tar...'")
        shell("docker save -o 00-Utils/Docker-images/notebook-executor.tar notebook-executor")
        shell("echo 'Docker image saved.'")


rule create_docker_image_jupyter_server:
    input:
        "00-Utils/Docker-jupyter-server/Dockerfile",
        "00-Utils/Docker-jupyter-server/requirements.txt"
    output:
        "00-Utils/Docker-images/jupyter-server.tar"
    log:
        "logs/setup_jupyter-server.log"
    run:
        shell("echo 'Building Docker image...'")
        shell("docker build -t jupyter-server -f 00-Utils/Docker-jupyter-server/Dockerfile .")
        shell("echo 'Saving Docker image to tar...'")
        shell("docker save -o 00-Utils/Docker-images/jupyter-server.tar jupyter-server")
        shell("echo 'Docker image saved.'")

rule run_jupyter_server:
    input: 
        "00-Utils/Docker-images/jupyter-server.tar"
    output:
        "serverRunning"
    run:
        #shell("echo 'Unpacking jupyter-server'"),
        #shell("docker load -i 00-Utils/Docker-images/jupyter-server.tar"),
        shell("echo 'Starting jupyter-server docker instance'"),
        shell(f"docker run -it -v {cwd}:/app -d -p 8888:8888 jupyter-server"),
        shell("echo 'Starting jupyter-server: http://localhost:8888/tree?'")
        shell("echo 'running' > serverRunning")
        shell("start http://127.0.0.1:8888/tree?")