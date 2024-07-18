from qsprpred.models import SklearnModel
from sys import argv

model_name = 'models/P10828_RF_Model/P10828_RF_Model_meta.json'

# smiles = ['Cc1c(Cc2cc(I)c(OCC(=O)O)c(I)c2)c2c(cccc2)o1', 'O=c1cnn(-c2cc(Cl)c(Oc3ccc(O)c(S(=O)(=O)N4CCc5ccccc54)c3)c(Cl)c2)c(=O)[nH]1']

# Collect all SMILES strings from the command line arguments
smiles = argv[1:]  # This takes all arguments after the script name

# Ensure at least one SMILES is provided
if not smiles:
    print("Error: No SMILES strings provided.")
    exit(1)

model = SklearnModel.fromFile(model_name)
predictions = model.predictMols(smiles)

# Print each SMILES string with its corresponding prediction
for sm, pred in zip(smiles, predictions):
    print(f"SMILES: {sm} -> Prediction: {pred[0]:.4f}")

