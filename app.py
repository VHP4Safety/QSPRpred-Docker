from flask import Flask, request, render_template, jsonify
from qsprpred.models import SklearnModel
import os
import json
import traceback

app = Flask(__name__)

# Define the models directory
MODELS_DIR = 'models'

@app.route('/')
def home():
    # List available models
    available_models = [d for d in os.listdir(MODELS_DIR) if os.path.isdir(os.path.join(MODELS_DIR, d))]
    return render_template('index.html', models=available_models)

@app.route('/predict', methods=['POST'])
def predict():
    try:
        smiles = request.form.get('smiles')
        model_name = request.form.get('model')
        
        if not smiles:
            return render_template('index.html', error="No SMILES strings provided.")
        if not model_name:
            return render_template('index.html', error="No model selected.")
        
        model_path = os.path.join(MODELS_DIR, model_name, f"{model_name}_meta.json")
        
        if not os.path.exists(model_path):
            return render_template('index.html', error="Selected model does not exist.")
        
        model = SklearnModel.fromFile(model_path)
        smiles_list = [smile.strip() for smile in smiles.split(',')]
        predictions = model.predictMols(smiles_list)
        predictions_formatted = [f"{pred[0]:.4f}" for pred in predictions]
        
        # Prepare data as a list of tuples
        data = list(zip(smiles_list, predictions_formatted))
        
        available_models = [d for d in os.listdir(MODELS_DIR) if os.path.isdir(os.path.join(MODELS_DIR, d))]
        return render_template('index.html', models=available_models, data=data)
    except Exception as e:
        # Print the exception to the Flask server logs
        print(f"Exception occurred: {e}")
        traceback.print_exc()
        available_models = [d for d in os.listdir(MODELS_DIR) if os.path.isdir(os.path.join(MODELS_DIR, d))]
        return render_template('index.html', models=available_models, error="An error occurred while processing the request.")

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)  # Enable debug mode

