from flask import Flask, request, render_template, jsonify, Response
from qsprpred.models import SklearnModel
import os
import json
import traceback
import pandas as pd

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
        smiles_input = request.form.get('smiles')
        uploaded_file = request.files.get('file')
        model_names = request.form.getlist('model')
        
        if not model_names:
            return render_template('index.html', error="No model selected.")
        
        if smiles_input:
            smiles_list = []
            smiles_list = [smile.strip() for smile in smiles_input.split(',')]
        
        if uploaded_file:
            smiles_list = []
            uploaded_df = pd.read_csv(uploaded_file)
            if 'SMILES' in uploaded_df.columns:
                smiles_list.extend(uploaded_df['SMILES'].tolist())
        
        if not smiles_list:
            return render_template('index.html', error="No SMILES strings provided.")
        
        all_predictions = {}
        
        for model_name in model_names:
            model_path = os.path.join(MODELS_DIR, model_name, f"{model_name}_meta.json")
            if not os.path.exists(model_path):
                return render_template('index.html', error=f"Model {model_name} does not exist.")
            
            model = SklearnModel.fromFile(model_path)
            predictions = model.predictMols(smiles_list)
            predictions_formatted = [f"{pred[0]:.4f}" for pred in predictions]
            all_predictions[model_name] = predictions_formatted
        
        # Prepare data for table
        table_data = []
        for i, smile in enumerate(smiles_list):
            row = [smile]
            for model_name in model_names:
                row.append(all_predictions[model_name][i])
            table_data.append(row)
        
        # Prepare header
        headers = ['SMILES'] + [f'Prediction ({model})' for model in model_names]
        
        available_models = [d for d in os.listdir(MODELS_DIR) if os.path.isdir(os.path.join(MODELS_DIR, d))]
        return render_template('index.html', models=available_models, headers=headers, data=table_data)
    except Exception as e:
        print(f"Exception occurred: {e}")
        traceback.print_exc()
        available_models = [d for d in os.listdir(MODELS_DIR) if os.path.isdir(os.path.join(MODELS_DIR, d))]
        return render_template('index.html', models=available_models, error="An error occurred while processing the request.")

@app.route('/api', methods=['POST'])
def apipredict():
    data = request.json
    smiles = data.get('smiles', [])
    output_format = data.get('format', 'json')  # Default to JSON format

    if not smiles or not isinstance(smiles, list):
        return jsonify({'error': 'Invalid input: please provide a list of SMILES strings.'}), 400

    predictions = model.predictMols(smiles)
    result = [{'smiles': sm, 'prediction': pred[0]} for sm, pred in zip(smiles, predictions)]

    if output_format == 'text':
        result_text = '\n'.join([f"SMILES: {sm} -> Prediction: {pred[0]:.4f}" for sm, pred in zip(smiles, predictions)])
        return Response(result_text, mimetype='text/plain')

    return jsonify(result)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)

