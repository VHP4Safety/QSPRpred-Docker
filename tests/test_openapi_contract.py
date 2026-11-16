"""Guards against openapi.yaml and the code drifting apart.

These tests deliberately avoid importing ``app.py``, which pulls in rdkit,
qsprpred and vina. The ``/api/v1`` rules come from registering the blueprint on
a throwaway Flask app, and the legacy routes are read statically out of
``app.py``'s source. That keeps the contract test runnable anywhere.
"""

import json
import pathlib
import re

import jsonschema
import pytest
import yaml
from flask import Flask

REPO = pathlib.Path(__file__).resolve().parent.parent
SPEC_PATH = REPO / "openapi.yaml"

# Routes that are intentionally NOT described in openapi.yaml. They build
# filesystem paths out of raw query parameters, so publishing them as part of a
# machine-readable contract would advertise that surface to automated clients.
# Documenting them is blocked on hardening them first.
UNDOCUMENTED_BY_DESIGN = {"/download", "/downloadqprf", "/downloadposes"}


@pytest.fixture(scope="module")
def spec():
    return yaml.safe_load(SPEC_PATH.read_text())


@pytest.fixture(scope="module")
def v1_rules():
    """The real /api/v1 rules, via a throwaway app (no heavy imports)."""
    import sys

    sys.path.insert(0, str(REPO))
    from api_v1 import bp

    probe = Flask(__name__)
    probe.register_blueprint(bp)
    return {
        rule.rule
        for rule in probe.url_map.iter_rules()
        if rule.rule.startswith("/api/v1")
    }


@pytest.fixture(scope="module")
def app_routes():
    """Route paths declared in app.py, read from source rather than imported."""
    source = (REPO / "app.py").read_text()
    return set(re.findall(r"@app\.route\(\s*'([^']+)'", source))


def _spec_paths(spec, prefix=None):
    paths = set(spec["paths"])
    if prefix is None:
        return paths
    return {path for path in paths if path.startswith(prefix)}


def test_every_v1_route_is_documented(spec, v1_rules):
    """No undocumented endpoint may exist under /api/v1."""
    assert v1_rules == _spec_paths(spec, "/api/v1")


def test_every_documented_path_exists(spec, v1_rules, app_routes):
    """The document may not describe endpoints the service does not serve."""
    served = v1_rules | app_routes
    missing = _spec_paths(spec) - served
    assert not missing, f"documented but not served: {sorted(missing)}"


def test_file_download_routes_stay_undocumented(spec, app_routes):
    """Pin the deliberate omission, so it is a decision and not an oversight."""
    assert UNDOCUMENTED_BY_DESIGN <= app_routes, "expected routes have moved"
    assert not (UNDOCUMENTED_BY_DESIGN & _spec_paths(spec))


def test_spec_version_matches_api_version(spec):
    import sys

    sys.path.insert(0, str(REPO))
    import api_v1

    assert spec["info"]["version"] == api_v1.API_VERSION


def _as_2020_12(schema, spec):
    """Inline components.schemas as $defs so jsonschema can resolve the $refs."""
    document = {"$defs": spec["components"]["schemas"], **schema}
    rewritten = json.dumps(document).replace(
        "#/components/schemas/", "#/$defs/"
    )
    return json.loads(rewritten)


def _iter_examples(spec):
    for path, operations in spec["paths"].items():
        for method, operation in operations.items():
            if not isinstance(operation, dict):
                continue
            body = operation.get("requestBody")
            if body:
                for media in body.get("content", {}).values():
                    for name, example in (media.get("examples") or {}).items():
                        yield (
                            f"{method.upper()} {path} request[{name}]",
                            media.get("schema"),
                            example["value"],
                        )
            for status, response in (operation.get("responses") or {}).items():
                for media in (response.get("content") or {}).values():
                    for name, example in (media.get("examples") or {}).items():
                        yield (
                            f"{method.upper()} {path} {status}[{name}]",
                            media.get("schema"),
                            example["value"],
                        )


def test_spec_is_a_valid_openapi_document():
    """Validate against the OpenAPI 3.1 meta-schema, not just our own checks."""
    from openapi_spec_validator import validate
    from openapi_spec_validator.readers import read_from_filename

    document, _ = read_from_filename(str(SPEC_PATH))
    validate(document)


def test_spec_has_examples(spec):
    assert len(list(_iter_examples(spec))) >= 8


def test_every_example_validates(spec):
    """Every example in the document must satisfy its own schema.

    This is what stops the examples - the part an LLM imitates most directly -
    from describing a response the service would never send.
    """
    failures = []
    for label, schema, value in _iter_examples(spec):
        if not schema:
            continue
        try:
            jsonschema.Draft202012Validator(_as_2020_12(schema, spec)).validate(value)
        except jsonschema.ValidationError as exc:
            failures.append(f"{label}: {exc.message} at {list(exc.absolute_path)}")
    assert not failures, "invalid examples:\n" + "\n".join(failures)


def test_error_codes_match_implementation(spec):
    """The documented error codes and the raised ones must be the same set."""
    source = (REPO / "api_v1.py").read_text()
    raised = set(re.findall(r'_error\(\s*"([a-z_]+)"', source))
    documented = set(
        spec["components"]["schemas"]["Error"]["properties"]["error"]["properties"][
            "code"
        ]["enum"]
    )
    assert raised == documented, (
        f"only raised: {sorted(raised - documented)}; "
        f"only documented: {sorted(documented - raised)}"
    )


def test_llms_txt_mentions_every_installed_model():
    """llms.txt names a static model list; keep it honest."""
    llms = (REPO / "llms.txt").read_text()
    installed = {
        path.name
        for path in (REPO / "models").iterdir()
        if path.is_dir() and (path / f"{path.name}_meta.json").is_file()
    }
    missing = {name for name in installed if name not in llms}
    assert not missing, f"models absent from llms.txt: {sorted(missing)}"


def test_targets_file_is_loadable_and_annotates_installed_models():
    targets = json.loads((REPO / "models" / "targets.json").read_text())
    annotated = {key for key in targets if not key.startswith("_")}
    installed = {
        path.name
        for path in (REPO / "models").iterdir()
        if path.is_dir() and (path / f"{path.name}_meta.json").is_file()
    }
    assert annotated <= installed, f"annotates unknown models: {sorted(annotated - installed)}"
    for name in annotated:
        assert "uniprot" in targets[name], f"{name} has no uniprot accession"
        assert "source" in targets[name], f"{name} does not record its provenance"
