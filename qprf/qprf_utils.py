import os
from datetime import datetime
from importlib.metadata import metadata, packages_distributions, version

import pystow
import qsprpred
from docxtpl import DocxTemplate
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

DATA_DIR = pystow.join("QSPRpred", "tutorials")
MODELS_DIR = pystow.join("QSPRpred", "tutorials", "tutorial_output", "models")


enumerator = rdMolStandardize.TautomerEnumerator()

def general_info():
    context = {
        "general": {},
        "substance": {},
        "model": {},
        "prediction": {},
        "input": {},
        "ad": {},
        "reliability": {},
        "analogues": {},
        "purpose": {}
    }

    context["general"]["date_QPRF"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    context["general"]["authors"] = ["Linde Schoenmaker (l.schoenmaker@lacdr.leidenuniv.nl)"]
    context["model"]["availability"] = "Model is non-proprietary: training and test sets are available in model repository X"
    context["prediction"]["property"] = "Protein-binding"
    context["prediction"]["dependent_variable"] = "For modelling purposes the bioactivity values were transformed to logarithmic units. The dependent variable is defined as: -Log(molar IC50, XC50, EC50, AC50, Ki, Kd or Potency)."
    context["input"]["descriptors"] = "Descriptors are same as those used for model development and validation"
    context["ad"]["description"] = "k-Nearest Neighbors (KNN) is used for AD evaluation. The distance of the predicted compound to the nearest neighbors in the training set is compared to a threshold. The applicability domain threshold is based on the 95% percentile of the training set. Attached in the supporting information is a figure of the residuals plotted against the KNN distance."
    context["reliability"]["reproducibility"] = "The software implementing the model is publicly available and the model is clearly described."
    context["input"]["model"] = "No custom settings were used to generate the prediction."
    return context


def substance_info(context, smiles):
    context["substance"]["SMILES"] = smiles
    mol = Chem.MolFromSmiles(smiles)
    context["substance"]["stereochemical"] = f"Number of stereocenters {len(Chem.FindPotentialStereo(mol))}. Assesed by FindPotentialStereo from rdkit version {version('rdkit')}"

    context["input"]["structure"] = f"The kind of input used to represent the input structure is SMILES. The associated value is {smiles}"
    context["input"]["stereochemistry"] = f"Number of stereocenters {len(Chem.FindPotentialStereo(mol))}. Assesed by FindPotentialStereo from rdkit version {version('rdkit')}. The model was trained on molecules with stereochemical features removed. If a molecule has stereoisomers the reliability of predictions is lower."
    tauts = enumerator.Enumerate(mol)
    context["input"]["tautomerism"] = f"Number of tautomers: {len(tauts)}. Assesed by TautomerEnumerator from rdkit version {version('rdkit')}"


    return context

def model_info(context, model):
    context["model"]["identifier"] = model.name
    software_package = packages_distributions()[qsprpred.__name__.split(".")[0]][0]
    context["model"]["software"] = software_package
    context["model"]["version"] = version(software_package)
    name = qsprpred.__package__
    context["model"]["reference"] = metadata(name).get_all("Project-URL")

    return context

def prediction_info(context, prediction, model):
    context["prediction"]["value"] = prediction
    context["prediction"]["unit"] = "pChEMBL value"
    if model.task.isClassification():
        context["prediction"]["value"] = f"Cut-off values: {', '.join(model.targetProperties[0].th)}"
        context["prediction"]["unit"] = "The prediction is a classification and therefore has no unit. Original training values had pChEMBL units."
    return context

def ad_info(context, value, model):
    if model.applicabilityDomain.threshold:
        assessment = f"Input is {'within' if value <= model.applicabilityDomain.threshold else 'not within'} AD."
        context["ad"]["value"] = f"Input is has distance of {value}. Inputs with distance {model.applicabilityDomain.direction} {model.applicabilityDomain.threshold} are within AD"
    else:
        assessment = f"Input is {'within' if value else 'not within'} AD."
    context["ad"]["assessment"] = assessment
    context["reliability"]["descriptor"] = assessment
    return context

def nn_info(context, nearest_neighbor):
    context["analogues"] = nearest_neighbor
    context["analogues"]["source"] = "Nearest neighbor was picked out of the dataset used for model training."
    if isinstance(nearest_neighbor["predicted_value"], float): 
        context["analogues"]["accuracy"] = f'Difference between experimental value and predicted value is {abs(nearest_neighbor["predicted_value"]-nearest_neighbor["value"]):.2f}. It should be noted that the nearest neighbor was part of the model training set so the accuracy is expected to be high.'

    return context

def render_qprf(smile, model, prediction, ad, nearest_neighbor):
    context = general_info()
    context = substance_info(context, smile)
    context = prediction_info(context, prediction, model)
    context = ad_info(context, ad, model)
    context = model_info(context, model)
    context = nn_info(context, nearest_neighbor)

    tpl = DocxTemplate(
        "/usr/src/app/qprf/qsar-assessment-framework-annex-2-qsar-prediction-reporting-format.docx"
    )
    tpl.render(context)

    directory = f"qprf/output/{model.name}"
    if not os.path.exists(directory):
        os.makedirs(directory)
    tpl.save(f"qprf/output/{model.name}/{smile}.docx")