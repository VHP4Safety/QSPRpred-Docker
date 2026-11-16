# QSPRpred-Docker

A Docker-based prediction service built on the QSPRpred tool [github.com/CDDLeiden/QSPRpred](https://github.com/CDDLeiden/QSPRpred). It provides a web UI and a documented JSON API ([OpenAPI 3.1](openapi.yaml)) for predicting chemical bioactivity at molecular initiating events (MIEs) from the VHP4Safety case studies (thyroid, Parkinson's disease, nephrotoxicity), using imported pre-trained models.

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

The API is described by [`openapi.yaml`](openapi.yaml) (OpenAPI 3.1), served by the
running service:

| URL | What it is |
| --- | --- |
| [`/docs`](https://qsprpred.cloud.vhp4safety.nl/docs) | Browsable API reference |
| [`/openapi.yaml`](https://qsprpred.cloud.vhp4safety.nl/openapi.yaml) | The machine-readable contract (also `/openapi.json`) |
| [`/llms.txt`](https://qsprpred.cloud.vhp4safety.nl/llms.txt) | Short orientation for LLM agents, including how to read the output |

#### Which models can be used how

Not all models work on all surfaces. `GET /api/v1/models` reports this per model, so a
client never has to guess:

| Models | JSON API | Web UI | Applicability domain |
| --- | --- | --- | --- |
| `TRalpha`, `TRbeta` (thyroid) | yes | yes | yes |
| `4or2_final`, `8t6j_b_final` (Parkinson's) | **no** — `supports_batch_api: false` | yes, ~45 s/molecule | yes |
| the five archived models | yes | **no** — no training set is shipped | **no** |

#### Listing the models

```sh
curl https://qsprpred.cloud.vhp4safety.nl/api/v1/models
```

#### Predicting

```sh
curl -X POST https://qsprpred.cloud.vhp4safety.nl/api/v1/predict \
  -H "Content-Type: application/json" \
  -d '{"smiles": ["CC(=O)OC1=CC=CC=C1C(=O)O"], "models": ["TRalpha"]}'
```

```json
{
  "api_version": "1.0.0",
  "service": {"name": "qsprpred-docker", "version": "1.1.0", "qsprpred_version": "3.1.1"},
  "request": {"n_input": 1, "n_predicted": 1, "models": ["TRalpha"]},
  "results": [
    {
      "input_smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
      "standardized_smiles": "CC(=O)Oc1ccccc1C(=O)O",
      "valid": true,
      "predictions": [
        {
          "model": "TRalpha", "task": "REGRESSION",
          "value": 5.9989, "value_type": "pchembl", "unit": "-log10(mol/L)", "label": null,
          "applicability_domain": {"supported": true, "within": true}
        }
      ]
    }
  ],
  "warnings": []
}
```

Predictions are made on the **standardized** structure, which can differ from the SMILES
submitted; both are returned. `applicability_domain.supported: false` means the model has
no applicability domain, so no check was made — which is not the same as passing it. An
unparseable SMILES is reported in `warnings` and does not fail the request. Errors use a
single envelope with a stable `error.code` (`unknown_model`,
`model_not_available_via_api`, `too_many_smiles`, ...); see `/docs`.

#### Deprecated: `POST /api` and `GET /get_models`

The unversioned endpoints still work and are unchanged, but new clients should use
`/api/v1/*`. `POST /api` builds its response keys from model names
(`"prediction (TRalpha)"`), returns predicted values as **strings**, and crashes the
worker with a 502 if given a model that requires docking. `/api/v1/predict` returns a
clean `400 model_not_available_via_api` instead.

<details>
<summary>Legacy <code>POST /api</code> examples (unchanged)</summary>

JSON output:
```sh
curl -X POST localhost:5000/api -H "Content-Type: application/json" -d '{
        "smiles": ["C1=CC=CC=C1C(=O)NC2=CC=CC=C2", "CC(=O)OC1=CC=CC=C1C(=O)O"],
        "models": ["TRalpha", "TRbeta"],
        "format": "json"
    }'
```
CSV output:
```sh
curl -X POST localhost:5000/api -H "Content-Type: application/json" -d '{
        "smiles": ["C1=CC=CC=C1C(=O)NC2=CC=CC=C2", "CC(=O)OC1=CC=CC=C1C(=O)O"],
        "models": ["TRalpha", "TRbeta"],
        "format": "csv"
    }' -o predictions.csv
```
</details>

### Local Development Setup
To run the application locally without Docker:
```sh
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python app.py
```
Then access the application at http://localhost:5000.

### Running the tests
```sh
pip install -r requirements-dev.txt
pytest
```
The suite covers the `/api/v1` request/response contract and guards `openapi.yaml`
against drifting from the code: every `/api/v1` route must be documented, every
documented path must be served, every example in the document must validate against its
own schema, and the raised error codes must match the documented ones.

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
- **Archived models work only through the API.** They ship without their training set, which the web UI needs for its nearest-neighbour panel, so selecting one in the UI fails. `GET /api/v1/models` reports this as `training_data_available: false`.
- **No nephrotoxicity model is installed.** The case study is described above and the UI has a tab for it, but no model currently ships for it.
