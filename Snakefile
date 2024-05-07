

rule setup_docker:
    input:
        "00 - Utils/Docker/Dockerfile",
        "00 - Utils/Docker/requirements.txt"
    log:
        "logs/setup_docker.log"
    shell:
        """
        start "Jupyter Server" /D "00 - Utils/Docker" "startDocker.bat" > {log} 2>&1
        pause
        """

    