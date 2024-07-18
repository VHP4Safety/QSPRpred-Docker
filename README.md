# QSPRpred-Docker

This GitHub repo is set up for the development of a Docker-based service on the QSPRpred tool [github.com/CDDLeiden/QSPRpred](https://github.com/CDDLeiden/QSPRpred). In the initial stage, the service will have prediction functionalities based on imported, pre-trained model(s). 


## Example Use 

One of the authors have provided us an example for agonists on the thyroid hormone receptor. Following are the steps to conduct the example. 

1. Create a Conda environment with Python >= 3.10
2. Install the qsprpred module with `pip install qsprpred`
3. Navigate into the `docker_test` directory that has the material to conduct the example
4. Run `python predict.py`

Following the steps above should give a result like `[[6.14965] [9.2567]]`. 