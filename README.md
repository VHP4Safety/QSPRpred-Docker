# QSPRpred-Docker

A Docker-based prediction service built on the QSPRpred tool [github.com/CDDLeiden/QSPRpred](https://github.com/CDDLeiden/QSPRpred). It provides a web UI and a JSON/CSV API for predicting chemical bioactivity at molecular initiating events (MIEs) from the VHP4Safety case studies (thyroid, Parkinson's disease, nephrotoxicity), using imported pre-trained models.

> ## 📦 Project-end version (v1.0.0)
>
> This is the **project-end milestone** of the QSPRpred-Docker service, developed during
> the [VHP4Safety](https://www.vhp4safety.nl/) project. Active development under the project
> has concluded, but the repository **remains open** and is not archived — the backlog of
> ideas for continuing the tool is kept in the [open issues](https://github.com/VHP4Safety/QSPRpred-Docker/issues).
>
> - **Live service:** https://qsprpred.cloud.vhp4safety.nl (available while the VHP4Safety cluster is maintained)
> - **Container image:** `ghcr.io/vhp4safety/qsprpred-docker:latest` (see [Self-hosting](#self-hosting) to run it yourself)
> - **Built on:** QSPRpred 3.1.1
> - Outstanding ideas and caveats are summarized under [Future Work](#future-work) and [Known Limitations](#known-limitations) (and tracked in the issues).
>
> <!-- DOI-BADGE -->

## Flask Application Deployment with Docker
[![Publish Docker image](https://github.com/VHP4Safety/QSPRpred-Docker/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/VHP4Safety/QSPRpred-Docker/actions/workflows/docker-publish.yml)

This guide covers the setup and usage of a Flask application for predictions via a Docker container. The application provides both a user interface for interacting with the model and an API for programmatic access.

### Building the Docker Image
To build the Docker image for the Flask application, use the following command:
```sh
docker build -t qspr_flask_image -f Dockerfile .
```
### Running the Docker Container
Run the Docker container with the following command. This command also mounts a local directory to the container to make model files accessible:
```sh
docker run -d -p 5000:5000 --name qspr_flask_container -v $(pwd)/models:/usr/src/app/models  qspr_flask_image
```

### Self-hosting
You do not need to build the image — the project-end image is published to the GitHub
Container Registry with the pretrained models, data, and QMRF documents baked in. Pull
and run it directly:
```sh
docker run -d -p 5000:5000 --name qsprpred ghcr.io/vhp4safety/qsprpred-docker:latest
```
Then open http://localhost:5000. No volume mount is required (the models live inside the
image). This is the same image that runs the live service at
https://qsprpred.cloud.vhp4safety.nl.

### Accessing the Application
Once the container is running, navigate to http://localhost:5000 in your web browser. You will see the Flask application interface which allows you to:

- Select one or more prediction models.
- Input multiple SMILES strings either through a text box or by uploading a CSV file (e.g. the [smiles_sample.csv](smiles_sample.csv)).
- Download the prediction results in various formats, or generate a report.

![QSPRpred UI](qsprpred-ui.png?raw=true "UI")

### Using the API
You can also interact with the Flask application via its API from your coding environment.
Use `GET /get_models` to list the available models and their metadata:
```sh
curl localhost:5000/get_models
```
For models that define an applicability domain, prediction responses include a
`within_applicability_domain (<model>)` field alongside the predicted value. Below are
examples of how to use the prediction endpoint:

#### Example: JSON Output
To initiate a prediction and receive the results in JSON format, use:
```sh
curl -X POST localhost:5000/api     -H "Content-Type: application/json"     -d '{
        "smiles": ["C1=CC=CC=C1C(=O)NC2=CC=CC=C2", "CC(=O)OC1=CC=CC=C1C(=O)O"],
        "models": ["TRalpha", "TRbeta"],
        "format": "json"
    }'
```
#### Example: CSV Output
To receive the prediction results as a CSV file, use:
```sh
curl -X POST localhost:5000/api     -H "Content-Type: application/json"     -d '{
        "smiles": ["C1=CC=CC=C1C(=O)NC2=CC=CC=C2", "CC(=O)OC1=CC=CC=C1C(=O)O"],
        "models": ["TRalpha", "TRbeta"],
        "format": "csv"
    }' -o predictions.csv
```

### Local Development Setup
To run the application locally without Docker:
```sh
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python app.py
```
Then access the application at http://localhost:5000.

### IFP/Docking Models
The application includes interaction fingerprint (IFP) models that use molecular docking:
- **4or2_final** - mGluR1 receptor model
- **8t6j_b_final** - mGluR5 receptor model

These models use AutoDock Vina for molecular docking and are slower than MorganFP-based models (~30-60 seconds per molecule). They are available in the "Parkinson's Disease" tab.

## Future Work
The project concluded before these enhancements were implemented. They are recorded here for anyone who continues the tool (original issue numbers in brackets):

- **Richer model metadata / provenance** — surface what each model was trained on, its container/version details, and pretrained-model information on the tiles ([#11](https://github.com/VHP4Safety/QSPRpred-Docker/issues/11), [#30](https://github.com/VHP4Safety/QSPRpred-Docker/issues/30), [#34](https://github.com/VHP4Safety/QSPRpred-Docker/issues/34)).
- **Analogue searcher** — a similarity-search UI against a model's training set ([#26](https://github.com/VHP4Safety/QSPRpred-Docker/issues/26)).
- **Model performance view** — show held-out performance metrics/plots per model ([#28](https://github.com/VHP4Safety/QSPRpred-Docker/issues/28)).
- **Structural-difference highlighting** — highlight structural differences in the nearest-neighbour table ([#40](https://github.com/VHP4Safety/QSPRpred-Docker/issues/40)).
- **Persist generated data** — save prediction outputs on the platform across sessions ([#33](https://github.com/VHP4Safety/QSPRpred-Docker/issues/33)).
- **ChEMBL lookup** — look up experimental pChEMBL values in ChEMBL for submitted compounds ([#10](https://github.com/VHP4Safety/QSPRpred-Docker/issues/10)).

## Known Limitations
- **pChEMBL is an aggregated proxy.** Predicted pChEMBL values are derived from heterogeneous assays (different assay types, conditions, and biological relevance), so they should be interpreted as a relative potency indicator, not an absolute measurement ([#35](https://github.com/VHP4Safety/QSPRpred-Docker/issues/35)).
- **Predictions are point estimates.** Models do not output confidence/prediction intervals ([#29](https://github.com/VHP4Safety/QSPRpred-Docker/issues/29)).
- **Docking models are slow and strict.** The Parkinson's Disease IFP/docking models (`4or2_final`, `8t6j_b_final`) run AutoDock Vina (~30–60 s per molecule), require each batch to contain unique structures (identical molecules submitted together are rejected), and fail a request if any molecule cannot be docked.
- **No QMRF for archived models.** The models in the "Archived" tab do not ship a QMRF document; their tiles show "QMRF not available".
