rule create_docker_images_notebook_executor:
    input:
        "00-Utils/Docker-ipnyb-execution/Dockerfile",
        "00-Utils/Docker-ipnyb-execution/requirements.txt",
        "00-Utils/Docker-ipnyb-execution/execute_notebook.sh"
    output:
        "00-Utils/Docker-ipnyb-execution/image/notebook-executor.tar"
    log:
        "logs/setup_notebook-executor.log"
    run:
        shell("echo 'Building Docker image...'")
        shell("docker build -t notebook-executor -f 00-Utils/Docker-ipnyb-execution/Dockerfile .")
        shell("echo 'Saving Docker image to tar...'")
        shell("docker save -o 00-Utils/Docker-ipnyb-execution/image/notebook-executor.tar notebook-executor")
        shell("echo 'Docker image saved.'")
