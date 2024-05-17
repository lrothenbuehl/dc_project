# Use an official Python runtime as a parent image
FROM python:3.11

# Set the working directory in the container
WORKDIR /app

# Copy requirements file and script into the container
COPY /00-Utils/Docker-ipnyb-execution/requirements.txt /app/requirements.txt
COPY /00-Utils/Docker-ipnyb-execution/execute_notebook.sh /app/en.sh

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Make the script executable
RUN chmod +x /app/en.sh

# Debugging: List the contents of /app to ensure files are copied correctly
RUN sed -i 's/\r//' en.sh

# Set the entrypoint to the script
ENTRYPOINT ["bash", "./en.sh"]
