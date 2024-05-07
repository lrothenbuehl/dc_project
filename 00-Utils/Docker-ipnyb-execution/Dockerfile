# Use an official Python runtime as a parent image
FROM python:3.11

# Set the working directory in the container
WORKDIR /app

COPY . /

# Install Jupyter
RUN pip install --no-cache-dir -r "/00-Utils/Docker-ipnyb-execution/requirements.txt"

# Copy the execution script to a specific location and make it executable
COPY 00-Utils/Docker-ipnyb-execution/execute_notebook.sh /usr/local/bin/execute_notebook.sh
RUN chmod +x /usr/local/bin/execute_notebook.sh

# Set the entrypoint to the script
ENTRYPOINT ["/usr/local/bin/execute_notebook.sh"]