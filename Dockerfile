
# Pulling the Python image
FROM	python:3.12.4

# Installing the qsprpred module
RUN 	pip install qsprpred

# Copying the source files with the trained model
COPY 	/docker_test/ /home/

# Changing the working directory
WORKDIR /home

# Defining the entrypoint
CMD ["python3"]
