# QSPRpred-Docker

This GitHub repo is set up for the development of a Docker-based service on the QSPRpred tool [github.com/CDDLeiden/QSPRpred](https://github.com/CDDLeiden/QSPRpred). In the initial stage, the service will have prediction functionalities based on imported, pre-trained model(s). 


## Example Use 

One of the authors has provided us an example 'for agonists on the thyroid hormone receptor'. Following are the steps to conduct the example. 

1. Create a Conda environment with Python >= 3.10
2. Install the qsprpred module with `pip install qsprpred`
3. Navigate into the `docker_test` directory that has the material to conduct the example
4. Run `python predict.py`

Following the steps above should give a result like `[[6.14965] [9.2567]]`. 


## Creating and using Docker Containers 

### Container for Interactive Use --> to be removed, as this functionality is included in the Flask app

Find the `Dockerfile` in the main directory. One can create the Docker image to make predictions on the terminal interactively with this file by using the following command: 

```
docker build --tag qspr_test_image .
```

Then, the container can be run with the following command: 
```
docker run -d --name qspr_test_container qspr_test_image
```

This will run the container on the system. One can then make predictions based on a SMILE interactively on the terminal. Below is an example to present this utility with two SMILEs: 
```
docker exec qspr_test_container python /usr/src/app/docker_test/predict_interactive.py "Cc1c(Cc2cc(I)c(OCC(=O)O)c(I)c2)c2c(cccc2)o1" "O=c1cnn(-c2cc(Cl)c(Oc3ccc(O)c(S(=O)(=O)N4CCc5ccccc54)c3)c(Cl)c2)c(=O)[nH]1"
```

This should give `6.1497` and `9.2567` for two SMILES used in the example above. 

## Flask Application Deployment with Docker
This guide covers the setup and usage of a Flask application for predictions via a Docker container. The application provides both a user interface for interacting with the model and an API for programmatic access.

### Building the Docker Image
To build the Docker image for the Flask application, use the following command:
```sh
docker build -t qspr_flask_image -f Dockerfile-flask .
```
### Running the Docker Container
Run the Docker container with the following command. This command also mounts a local directory to the container to make model files accessible:
```sh
docker run -d -p 5000:5000 --name qspr_flask_container -v $(pwd)/models:/usr/src/app/models  qspr_flask_image
```

### Accessing the Application
Once the container is running, navigate to http://localhost:5000 in your web browser. You will see the Flask application interface which allows you to:

- Select one or more prediction models.
- Input multiple SMILES strings either through a text box or by uploading a CSV file.
- Download the prediction results in various formats.
![QSPRpred UI](/templates/img/FlaskUIQSPRpred.png?raw=true "UI")

### Using the API
You can also interact with the Flask application via its API from your coding environment. Below are examples of how to use the API:

#### Example: JSON Output
To initiate a prediction and receive the results in JSON format, use:
```sh
curl -X POST localhost:5000/api     -H "Content-Type: application/json"     -d '{
        "smiles": ["C1=CC=CC=C1C(=O)NC2=CC=CC=C2", "CC(=O)OC1=CC=CC=C1C(=O)O"],
        "models": ["P10827_RF_Model", "P10828_RF_Model"],
        "format": "text"
    }'
```
#### Example: CSV Output
To receive the prediction results as a CSV file, use:
```sh
curl -X POST localhost:5000/api     -H "Content-Type: application/json"     -d '{
        "smiles": ["C1=CC=CC=C1C(=O)NC2=CC=CC=C2", "CC(=O)OC1=CC=CC=C1C(=O)O"],
        "models": ["P10827_RF_Model", "P10828_RF_Model"],
        "format": "csv"
    }' -o predictions.csv
```

