"""Functional tests for the /api/v1 surface.

``api_v1`` receives its dependencies through :func:`init_api_v1`, so the entire
request path - validation, error codes, response shape, number formatting - is
exercised here with stubs, without rdkit, qsprpred, vina or the real model
files. The real models are covered by the end-to-end curls in the README/plan.
"""

import pathlib

import pytest
from flask import Flask

from api_v1 import init_api_v1

REPO = pathlib.Path(__file__).resolve().parent.parent

ASPIRIN = "CC(=O)OC1=CC=CC=C1C(=O)O"
ASPIRIN_STD = "CC(=O)Oc1ccccc1C(=O)O"


def _meta(name, **overrides):
    """Mimic one entry of app.extract_model_info()."""
    base = {
        "name": name,
        "pref_name": name,
        "case_study": "thyroid",
        "target_property_name": "pchembl_value_Mean",
        "target_property_task": "REGRESSION",
        "feature_calculator": "MorganFP",
        "algorithm": "RandomForestRegressor",
        "currDir": f"models/{name}",
        "has_qmrf": True,
        "has_applicability_domain": True,
    }
    base.update(overrides)
    return base


MODEL_META = [
    _meta("TRalpha"),
    _meta(
        "4or2_final",
        case_study="pd",
        target_property_name="active",
        target_property_task="SINGLECLASS",
        feature_calculator="DataFrameDescriptorSet",
        algorithm="XGBClassifier",
    ),
    _meta(
        "P31644_RF_Model",
        case_study="archived",
        pref_name="GABA",
        has_applicability_domain=False,
        has_qmrf=False,
    ),
    _meta(
        "StubClassifier",
        target_property_name="active",
        target_property_task="SINGLECLASS",
        algorithm="RandomForestClassifier",
    ),
]


class StubModel:
    """Stands in for qsprpred's SklearnModel."""

    def __init__(self, values, within=True, raises=False):
        self.values = values
        self.within = within
        self.raises = raises

    def predictMols(self, smiles, use_applicability_domain=False):
        if self.raises:
            raise RuntimeError("stub failure")
        predictions = [[self.values[index % len(self.values)]] for index in range(len(smiles))]
        if use_applicability_domain:
            return predictions, [[self.within] for _ in smiles]
        return predictions


def _standardize(smiles):
    """Stub standardizer: anything containing 'bad' is unparseable."""
    if "bad" in smiles:
        return None
    return ASPIRIN_STD if smiles == ASPIRIN else smiles.lower()


@pytest.fixture
def models():
    return {
        "TRalpha": StubModel([5.99887]),
        "P31644_RF_Model": StubModel([6.968]),
        "StubClassifier": StubModel([1]),
    }


@pytest.fixture
def client(models, tmp_path):
    # Only TRalpha gets a training-set directory, so training_data_available is
    # exercised in both states.
    (tmp_path / "TRalpha_Data").mkdir()
    (tmp_path / "TRalpha_Data" / "train_full_model_TRalpha.csv").write_text("SMILES\n")

    app = Flask(__name__)
    app.config["TESTING"] = True
    init_api_v1(
        app,
        extract_model_info=lambda directory: MODEL_META,
        standardize_smiles=_standardize,
        get_timing_stats=lambda: {"ifp_time_per_mol": 45.0, "standard_time_per_mol": 2.0},
        load_model=lambda name: models[name],
        # Real models dir so the curated targets.json is loaded for real.
        models_dir=str(REPO / "models"),
        data_dir=str(tmp_path),
        max_smiles=1000,
        service_version="1.1.0-test",
        qsprpred_version="3.1.1",
    )
    return app.test_client()


def _predict(client, **body):
    return client.post("/api/v1/predict", json=body)


# --- GET /api/v1/models ----------------------------------------------------


def test_models_lists_everything(client):
    payload = client.get("/api/v1/models").get_json()
    assert payload["api_version"] == "1.0.0"
    assert payload["service"]["qsprpred_version"] == "3.1.1"
    assert {model["name"] for model in payload["models"]} == {
        meta["name"] for meta in MODEL_META
    }


def test_models_never_leaks_internal_paths(client):
    payload = client.get("/api/v1/models").get_json()
    assert all("currDir" not in model for model in payload["models"])


def test_only_docking_models_are_excluded_from_the_batch_api(client):
    models = {m["name"]: m for m in client.get("/api/v1/models").get_json()["models"]}
    assert models["4or2_final"]["requires_docking"] is True
    assert models["4or2_final"]["supports_batch_api"] is False
    assert models["4or2_final"]["estimated_seconds_per_molecule"] == 45.0
    for name in ("TRalpha", "P31644_RF_Model", "StubClassifier"):
        assert models[name]["supports_batch_api"] is True
        assert models[name]["estimated_seconds_per_molecule"] == 2.0


def test_models_report_applicability_domain_and_training_data(client):
    models = {m["name"]: m for m in client.get("/api/v1/models").get_json()["models"]}
    assert models["TRalpha"]["applicability_domain_supported"] is True
    assert models["TRalpha"]["training_data_available"] is True
    assert models["P31644_RF_Model"]["applicability_domain_supported"] is False
    # Archived models ship without a training set: fine over the API, broken in the UI.
    assert models["P31644_RF_Model"]["training_data_available"] is False


def test_models_describe_value_semantics(client):
    models = {m["name"]: m for m in client.get("/api/v1/models").get_json()["models"]}
    assert models["TRalpha"]["value_type"] == "pchembl"
    assert models["TRalpha"]["unit"] == "-log10(mol/L)"
    assert models["4or2_final"]["value_type"] == "class_label"
    assert models["4or2_final"]["unit"] is None


def test_models_carry_curated_target_annotations(client):
    models = {m["name"]: m for m in client.get("/api/v1/models").get_json()["models"]}
    assert models["TRalpha"]["target"]["uniprot"] == "P10827"
    # Not in targets.json, so the field must be absent rather than guessed.
    assert "target" not in models["StubClassifier"]


# --- POST /api/v1/predict, happy paths -------------------------------------


def test_predict_returns_numbers_not_strings(client):
    payload = _predict(client, smiles=[ASPIRIN], models=["TRalpha"]).get_json()
    prediction = payload["results"][0]["predictions"][0]
    assert prediction["value"] == 5.9989  # rounded to 4 dp, as a JSON number
    assert isinstance(prediction["value"], float)
    assert prediction["unit"] == "-log10(mol/L)"
    assert prediction["label"] is None


def test_predict_reports_both_input_and_standardized_smiles(client):
    payload = _predict(client, smiles=[ASPIRIN], models=["TRalpha"]).get_json()
    result = payload["results"][0]
    assert result["input_smiles"] == ASPIRIN
    assert result["standardized_smiles"] == ASPIRIN_STD
    assert result["valid"] is True


def test_predict_distinguishes_unsupported_ad_from_outside_ad(client):
    supported = _predict(client, smiles=[ASPIRIN], models=["TRalpha"]).get_json()
    domain = supported["results"][0]["predictions"][0]["applicability_domain"]
    assert domain == {"supported": True, "within": True}

    unsupported = _predict(
        client, smiles=[ASPIRIN], models=["P31644_RF_Model"]
    ).get_json()
    domain = unsupported["results"][0]["predictions"][0]["applicability_domain"]
    assert domain == {"supported": False, "within": None}


def test_predict_marks_molecules_outside_the_domain(models, client):
    models["TRalpha"].within = False
    payload = _predict(client, smiles=[ASPIRIN], models=["TRalpha"]).get_json()
    assert payload["results"][0]["predictions"][0]["applicability_domain"] == {
        "supported": True,
        "within": False,
    }


def test_predict_labels_classification_output(client):
    payload = _predict(client, smiles=[ASPIRIN], models=["StubClassifier"]).get_json()
    prediction = payload["results"][0]["predictions"][0]
    assert prediction["value"] == 1
    assert prediction["label"] == "Active"
    assert prediction["task"] == "SINGLECLASS"


def test_predict_applies_every_model_to_every_molecule(client):
    payload = _predict(
        client,
        smiles=[ASPIRIN, "CCO"],
        models=["TRalpha", "P31644_RF_Model"],
    ).get_json()
    assert payload["request"] == {
        "n_input": 2,
        "n_predicted": 2,
        "models": ["TRalpha", "P31644_RF_Model"],
    }
    for result in payload["results"]:
        assert [p["model"] for p in result["predictions"]] == [
            "TRalpha",
            "P31644_RF_Model",
        ]


def test_predict_skips_unparseable_smiles_without_failing(client):
    payload = _predict(
        client, smiles=[ASPIRIN, "bad-molecule"], models=["TRalpha"]
    ).get_json()
    assert payload["request"]["n_input"] == 2
    assert payload["request"]["n_predicted"] == 1
    # Order is preserved and the invalid input is still represented.
    assert payload["results"][0]["valid"] is True
    assert payload["results"][1]["valid"] is False
    assert payload["results"][1]["standardized_smiles"] is None
    assert payload["results"][1]["predictions"] == []
    assert payload["warnings"] == [
        {
            "code": "invalid_smiles",
            "smiles": "bad-molecule",
            "message": "SMILES could not be parsed and was skipped.",
        }
    ]


def test_predict_keeps_order_when_the_first_input_is_invalid(client):
    payload = _predict(
        client, smiles=["bad", ASPIRIN], models=["TRalpha"]
    ).get_json()
    assert [r["valid"] for r in payload["results"]] == [False, True]
    assert payload["results"][1]["standardized_smiles"] == ASPIRIN_STD


def test_predict_deduplicates_repeated_model_names(client):
    payload = _predict(
        client, smiles=[ASPIRIN], models=["TRalpha", "TRalpha"]
    ).get_json()
    assert payload["request"]["models"] == ["TRalpha"]
    assert len(payload["results"][0]["predictions"]) == 1


# --- POST /api/v1/predict, guardrails --------------------------------------


def test_docking_models_are_rejected_cleanly(client):
    """The regression that matters: a clean 400, never a crashed worker."""
    response = _predict(client, smiles=[ASPIRIN], models=["4or2_final"])
    assert response.status_code == 400
    error = response.get_json()["error"]
    assert error["code"] == "model_not_available_via_api"
    assert error["details"]["models"] == ["4or2_final"]
    assert "/predict" in error["message"]


def test_unknown_model_is_rejected_with_the_available_list(client):
    response = _predict(client, smiles=[ASPIRIN], models=["TRgamma"])
    assert response.status_code == 400
    error = response.get_json()["error"]
    assert error["code"] == "unknown_model"
    assert error["details"]["unknown_models"] == ["TRgamma"]
    assert "TRalpha" in error["details"]["available_models"]


def test_batch_size_limit_is_enforced(client):
    response = _predict(client, smiles=[ASPIRIN] * 1001, models=["TRalpha"])
    assert response.status_code == 413
    error = response.get_json()["error"]
    assert error["code"] == "too_many_smiles"
    assert error["details"] == {"received": 1001, "limit": 1000}


def test_all_smiles_invalid_is_an_error(client):
    response = _predict(client, smiles=["bad1", "bad2"], models=["TRalpha"])
    assert response.status_code == 400
    assert response.get_json()["error"]["code"] == "no_valid_smiles"


@pytest.mark.parametrize(
    "body",
    [
        {},
        {"smiles": [ASPIRIN]},
        {"models": ["TRalpha"]},
        {"smiles": ASPIRIN, "models": ["TRalpha"]},
        {"smiles": [], "models": ["TRalpha"]},
        {"smiles": [ASPIRIN], "models": []},
        {"smiles": [1, 2], "models": ["TRalpha"]},
        {"smiles": [ASPIRIN], "models": [""]},
        {"smiles": [ASPIRIN], "models": "TRalpha"},
    ],
)
def test_malformed_bodies_are_rejected(client, body):
    response = _predict(client, **body)
    assert response.status_code == 400
    assert response.get_json()["error"]["code"] == "invalid_request"


def test_non_object_body_is_rejected(client):
    response = client.post("/api/v1/predict", json=[ASPIRIN])
    assert response.status_code == 400
    assert response.get_json()["error"]["code"] == "invalid_request"


def test_a_failing_model_returns_a_structured_error(models, client):
    models["TRalpha"].raises = True
    response = _predict(client, smiles=[ASPIRIN], models=["TRalpha"])
    assert response.status_code == 500
    error = response.get_json()["error"]
    assert error["code"] == "prediction_failed"
    assert error["details"]["model"] == "TRalpha"
