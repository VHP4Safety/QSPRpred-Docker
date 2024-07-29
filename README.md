# QSPRpred-Docker

This GitHub repo is set up for the development of a Docker-based service on the QSPRpred tool [github.com/CDDLeiden/QSPRpred](https://github.com/CDDLeiden/QSPRpred). In the initial stage, the service will have prediction functionalities based on imported, pre-trained model(s). 


## Example Use 

One of the authors has provided us an example 'for agonists on the thyroid hormone receptor'. Following are the steps to conduct the example. 

1. Create a Conda environment with Python >= 3.10
2. Install the qsprpred module with `pip install qsprpred`
3. Navigate into the `docker_test` directory that has the material to conduct the example
4. Run `python predict.py`

Following the steps above should give a result like `[[6.14965] [9.2567]]`. 


## Creating Docker Containers

### Container for Interactive Use

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

### Container for the Flask Application

There is another Dockerfile, `Dockerfile-flask`, along with files required for deploying a Flask app that runs the predictions. To build and run this, use the following commands:
```
docker build -t qspr_flask_image -f Dockerfile-flask .
```
The run command includes a volume so that new models can be added, and made accessible via the service. 
```
docker run -d -p 5000:5000 --name qspr_flask_container -v $(pwd)/models:/usr/src/app/models  qspr_flask_image
```


## Things to Do

- Adjust the script/container in a way that it will take input from the user. 
	- Please see `/docker_test/predict_interactive.py` for this. If it is correct, it should be taking (only one) input of SMILES from the terminal. So, if one runs `python predict_interactive.py "Cc1c(Cc2cc(I)c(OCC(=O)O)c(I)c2)c2c(cccc2)o1"`, the terminal should prompt `[[6.14965]]` which is the first output in the expected outcome. Note that there are two SMILES in the original example, `predict.py` which is why, I think, it prompts two values in the output. 
	- The `predict_interactive.py` script is now updated to receive more than one SMILEs (see [this line](https://github.com/VHP4Safety/QSPRpred-Docker/blob/a9a5abf2663c6194fb8e2815e8d9b2edf0dddb00/docker_test/predict_interactive.py#L9) for how it is set up).
- Create an UI for the "service". 
	- A new Dockerfile, `Dockerfile-flask`, has been added to the repo along with the source code for the Flask app. 
