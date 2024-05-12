# Use an official Python runtime as a parent image
FROM python:3.11

# Set the working directory to /app
WORKDIR /app

RUN mkdir output

# Copy the current directory contents into the container at /app
# Ensure that the Docker context is appropriately set to include the needed files

# Install any needed packages specified in requirements.txt
# Make sure the path to requirements.txt is correct relative to the context directory
COPY calculateDescriptors.py /calculateDescriptors.py
COPY requirements.txt /app/requirements.txt

RUN pip install --no-cache-dir -r '/app/requirements.txt'

# Make port 8888 available to the world outside this container
EXPOSE 8888

# Define environment variable
ENV NAME World

# Run Jupyter Notebook when the container launches using the correct JSON array format
CMD [ "python", "/calculateDescriptors.py" ]



