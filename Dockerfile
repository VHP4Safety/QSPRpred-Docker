FROM python:3.10-slim

# Set environment variables to prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libexpat1 \
    libxext6 \
    libx11-6 \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory in the container
WORKDIR /usr/src/app

# Install any needed packages specified in requirements.txt
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt

# Copy the current directory contents into the container at /usr/src/app
COPY models models
COPY static static
COPY templates templates
COPY app.py entrypoint.sh ./

RUN chmod +x /usr/src/app/entrypoint.sh

# Define the entrypoint script
ENTRYPOINT ["/usr/src/app/entrypoint.sh"]

