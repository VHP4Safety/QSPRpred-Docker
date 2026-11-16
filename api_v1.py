"""Versioned JSON API (``/api/v1``) for QSPRpred-Docker.

This module holds the *formal*, schematised API contract described by
``openapi.yaml``. It is deliberately separate from the legacy form/HTML routes
in ``app.py``, and it receives its dependencies through :func:`init_api_v1`
instead of importing them, for two reasons:

* it avoids a circular import (``app.py`` imports this module), and
* ``app.py`` must import ``ifp`` (which loads vina) *before* rdkit or the
  process segfaults, so this module must not drag rdkit/qsprpred in at import
  time ahead of that. See the import-order note at the top of ``app.py``.

Design notes for anyone extending this:

* Every response is an object (never a bare array) so new top-level fields can
  be added without breaking clients.
* Numbers are returned as JSON numbers, not strings.
* ``applicability_domain`` is an object so "this model has no AD" is
  distinguishable from "this molecule is outside the AD".
* Errors always use the same envelope: ``{"error": {"code", "message", ...}}``.
"""

import json
import logging
import os

from flask import Blueprint, current_app, jsonify, request

logger = logging.getLogger(__name__)

#: Version of the API contract itself (independent of the service version).
API_VERSION = "1.0.0"

#: Key under ``app.extensions`` holding the injected dependencies.
EXT_KEY = "qsprpred_api_v1"

#: Feature calculator used by the docking/interaction-fingerprint models. Those
#: models cannot be served from a plain SMILES list: their descriptors are
#: looked up from precomputed docking poses, so ``predictMols`` finds no row for
#: an unseen molecule. ``/predict`` (the web UI) handles them via ``ifp.predict``;
#: this API rejects them with a clear 400 instead.
DOCKING_FEATURE_CALCULATOR = "DataFrameDescriptorSet"

#: Human-readable labels for the ``case_study`` identifiers used in model metadata.
CASE_STUDY_LABELS = {
    "thyroid": "Thyroid",
    "pd": "Parkinson's disease",
    "nephrotox": "Nephrotoxicity",
    "archived": "Archived",
}

bp = Blueprint("api_v1", __name__, url_prefix="/api/v1")


def init_api_v1(
    app,
    *,
    extract_model_info,
    standardize_smiles,
    get_timing_stats,
    load_model,
    models_dir,
    data_dir,
    max_smiles,
    service_version,
    qsprpred_version,
):
    """Register the ``/api/v1`` blueprint on ``app`` with its dependencies.

    Must be called *after* ``app.py`` has imported ``ifp`` (see module docstring).
    """
    app.extensions[EXT_KEY] = {
        "extract_model_info": extract_model_info,
        "standardize_smiles": standardize_smiles,
        "get_timing_stats": get_timing_stats,
        "load_model": load_model,
        "models_dir": models_dir,
        "data_dir": data_dir,
        "max_smiles": max_smiles,
        "service_version": service_version,
        "qsprpred_version": qsprpred_version,
        "targets": _load_targets(models_dir),
    }
    app.register_blueprint(bp)
    return bp


def _load_targets(models_dir):
    """Load the hand-curated model -> target-identifier map, if present."""
    path = os.path.join(models_dir, "targets.json")
    try:
        with open(path) as handle:
            return json.load(handle)
    except FileNotFoundError:
        logger.warning("No %s found; model target annotations will be omitted.", path)
        return {}
    except (OSError, ValueError) as exc:
        logger.error("Could not read %s: %s", path, exc)
        return {}


def _deps():
    return current_app.extensions[EXT_KEY]


def _error(code, message, status, **details):
    """Build the single error envelope used by every ``/api/v1`` failure."""
    payload = {"error": {"code": code, "message": message}}
    if details:
        payload["error"]["details"] = details
    return jsonify(payload), status


def _service_block():
    deps = _deps()
    return {
        "name": "qsprpred-docker",
        "version": deps["service_version"],
        "qsprpred_version": deps["qsprpred_version"],
    }


def _value_semantics(task, target_property_name):
    """Describe what a prediction value means, so clients need not guess."""
    if task == "REGRESSION":
        if target_property_name.lower().startswith("pchembl"):
            return "pchembl", "-log10(mol/L)"
        return "numeric", None
    return "class_label", None


def describe_models():
    """Return the public capability description for every installed model.

    Built on top of ``app.extract_model_info`` and deliberately drops its
    ``currDir`` field, which is an internal filesystem path.
    """
    deps = _deps()
    timing = deps["get_timing_stats"]()
    described = []

    for info in deps["extract_model_info"](deps["models_dir"]):
        name = info["name"]
        requires_docking = info["feature_calculator"] == DOCKING_FEATURE_CALCULATOR
        task = info["target_property_task"]
        value_type, unit = _value_semantics(task, info["target_property_name"])

        # The nearest-neighbour step of the web UI reads the training set from
        # data/<accession>_Data/. Models without it predict fine over this API
        # but fail in the web UI, so say so rather than letting a client assume
        # the two surfaces are equivalent.
        accession = name.split("_")[0]
        train_csv = os.path.join(
            deps["data_dir"], f"{accession}_Data", f"train_full_model_{accession}.csv"
        )

        seconds = (
            timing["ifp_time_per_mol"] if requires_docking else timing["standard_time_per_mol"]
        )

        entry = {
            "name": name,
            "pref_name": info["pref_name"],
            "case_study": {
                "id": info["case_study"],
                "label": CASE_STUDY_LABELS.get(info["case_study"], info["case_study"]),
            },
            "task": task,
            "target_property_name": info["target_property_name"],
            "value_type": value_type,
            "unit": unit,
            "algorithm": info["algorithm"],
            "feature_calculator": info["feature_calculator"],
            "requires_docking": requires_docking,
            "supports_batch_api": not requires_docking,
            "applicability_domain_supported": bool(info.get("has_applicability_domain")),
            "training_data_available": os.path.isfile(train_csv),
            "estimated_seconds_per_molecule": round(float(seconds), 1),
            "max_smiles": deps["max_smiles"],
            "has_qmrf": info["has_qmrf"],
        }

        target = deps["targets"].get(name)
        if target:
            entry["target"] = target

        described.append(entry)

    described.sort(key=lambda m: (m["case_study"]["id"], m["name"]))
    return described


@bp.get("/models")
def models_v1():
    """List the available models and what each one can actually do."""
    return jsonify(
        {
            "api_version": API_VERSION,
            "service": _service_block(),
            "models": describe_models(),
        }
    )


def _validate_request(body, max_smiles):
    """Validate the request body. Returns an error response, or None if valid."""
    if not isinstance(body, dict):
        return _error(
            "invalid_request",
            "Request body must be a JSON object with 'smiles' and 'models' arrays.",
            400,
        )

    for field in ("smiles", "models"):
        value = body.get(field)
        if not isinstance(value, list) or not value:
            return _error(
                "invalid_request",
                f"'{field}' must be a non-empty array of strings.",
                400,
                field=field,
            )
        if not all(isinstance(item, str) and item.strip() for item in value):
            return _error(
                "invalid_request",
                f"'{field}' must contain only non-empty strings.",
                400,
                field=field,
            )

    received = len(body["smiles"])
    if received > max_smiles:
        return _error(
            "too_many_smiles",
            f"Too many SMILES submitted: {received} > {max_smiles}. "
            "Please split your input into smaller batches.",
            413,
            received=received,
            limit=max_smiles,
        )
    return None


def _format_prediction(name, meta, predictions, domain, position):
    """Turn one raw model output into the public prediction object."""
    raw = predictions[position]
    value = raw[0] if hasattr(raw, "__len__") else raw

    prediction = {
        "model": name,
        "task": meta["task"],
        "value_type": meta["value_type"],
        "unit": meta["unit"],
    }

    if meta["task"] == "REGRESSION":
        prediction["value"] = round(float(value), 4)
        prediction["label"] = None
    else:
        as_int = int(value)
        prediction["value"] = as_int
        prediction["label"] = "Active" if as_int == 1 else "Inactive"

    within = None
    if meta["applicability_domain_supported"] and domain is not None:
        try:
            flag = domain[position]
            within = bool(flag[0] if hasattr(flag, "__len__") else flag)
        except (TypeError, IndexError, KeyError):
            logger.warning("Could not read applicability domain for model %s", name)

    prediction["applicability_domain"] = {
        "supported": meta["applicability_domain_supported"],
        "within": within,
    }
    return prediction


@bp.post("/predict")
def predict_v1():
    """Predict the target property of one or more molecules with one or more models."""
    deps = _deps()

    body = request.get_json(silent=True)
    invalid = _validate_request(body, deps["max_smiles"])
    if invalid is not None:
        return invalid

    smiles_input = body["smiles"]
    catalog = {model["name"]: model for model in describe_models()}

    unknown = [name for name in body["models"] if name not in catalog]
    if unknown:
        return _error(
            "unknown_model",
            f"Unknown model(s): {', '.join(unknown)}.",
            400,
            unknown_models=unknown,
            available_models=sorted(catalog),
        )

    docking = [name for name in body["models"] if catalog[name]["requires_docking"]]
    if docking:
        return _error(
            "model_not_available_via_api",
            f"Model(s) {', '.join(docking)} require molecular docking and are not "
            "available over this API. Use the web UI at /predict instead.",
            400,
            models=docking,
        )

    model_names = list(dict.fromkeys(body["models"]))

    # Standardize once; a None result means rdkit could not parse the input.
    standardized = [deps["standardize_smiles"](smile) for smile in smiles_input]
    warnings = [
        {
            "code": "invalid_smiles",
            "smiles": smiles_input[index],
            "message": "SMILES could not be parsed and was skipped.",
        }
        for index, value in enumerate(standardized)
        if value is None
    ]

    valid_positions = {}
    valid_smiles = []
    for index, value in enumerate(standardized):
        if value is not None:
            valid_positions[index] = len(valid_smiles)
            valid_smiles.append(value)

    if not valid_smiles:
        return _error(
            "no_valid_smiles",
            "None of the submitted SMILES could be parsed.",
            400,
            invalid_smiles=smiles_input,
        )

    outputs = {}
    for name in model_names:
        try:
            model = deps["load_model"](name)
            if catalog[name]["applicability_domain_supported"]:
                predictions, domain = model.predictMols(
                    valid_smiles, use_applicability_domain=True
                )
            else:
                predictions, domain = model.predictMols(valid_smiles), None
        except Exception:
            logger.exception("Prediction failed for model %s", name)
            return _error(
                "prediction_failed",
                f"Prediction failed for model {name}.",
                500,
                model=name,
            )
        outputs[name] = (predictions, domain)

    results = []
    for index, raw_smiles in enumerate(smiles_input):
        value = standardized[index]
        entry = {
            "input_smiles": raw_smiles,
            "standardized_smiles": value,
            "valid": value is not None,
            "predictions": [],
        }
        if value is not None:
            position = valid_positions[index]
            for name in model_names:
                predictions, domain = outputs[name]
                entry["predictions"].append(
                    _format_prediction(name, catalog[name], predictions, domain, position)
                )
        results.append(entry)

    return jsonify(
        {
            "api_version": API_VERSION,
            "service": _service_block(),
            "request": {
                "n_input": len(smiles_input),
                "n_predicted": len(valid_smiles),
                "models": model_names,
            },
            "results": results,
            "warnings": warnings,
        }
    )
