from flask import Flask, request, render_template, jsonify, Response
from qsprpred.models import SklearnModel
import os
import io
import csv
import json
import traceback
import pandas as pd

app = Flask(__name__)

# Define the models directory
MODELS_DIR = 'models'

def extract_model_info(directory):
    models_info = []
    for d in os.listdir(directory):
        meta_path = os.path.join(directory, d, f"{d}_meta.json")
        if os.path.isfile(meta_path):
            with open(meta_path, 'r') as meta_file:
                meta_data = json.load(meta_file)
                state = meta_data['py/state']
                model_info = {
                    'name': state['name'],
                    'target_property_name': state['targetProperties'][0]['py/state']['name'],
                    'target_property_task': state['targetProperties'][0]['py/state']['task']['py/reduce'][1]['py/tuple'][0],
                    'feature_calculator': state['featureCalculators'][0]['py/object'].split('.')[-1],
                    'radius': state['featureCalculators'][0]['py/state']['radius'],
                    'nBits': state['featureCalculators'][0]['py/state']['nBits'],
                    'algorithm': state['alg'].split('.')[-1]
                }
                models_info.append(model_info)
    return models_info

@app.route('/')
def home():
    available_models = extract_model_info(MODELS_DIR)
    return render_template('index.html', models=available_models)

@app.route('/predict', methods=['POST'])
def predict():
    available_models = extract_model_info(MODELS_DIR)
    try:
        smiles_input = request.form.get('smiles')
        uploaded_file = request.files.get('file')
        model_names = request.form.getlist('model')
        
        if not model_names:
            return render_template('index.html', models=available_models, error="No model selected.")
        
        smiles_list = [smile.strip() for smile in smiles_input.split(',')] if smiles_input else []
        if uploaded_file:
            uploaded_df = pd.read_csv(uploaded_file)
            smiles_list.extend(uploaded_df['SMILES'].tolist())
        
        if not smiles_list:
            return render_template('index.html', models=available_models, error="No SMILES strings provided.")
        
        all_predictions = {}
        for model_name in model_names:
            model_path = os.path.join(MODELS_DIR, model_name, f"{model_name}_meta.json")
            model = SklearnModel.fromFile(model_path)
            predictions = model.predictMols(smiles_list)
            all_predictions[model_name] = [f"{pred[0]:.4f}" for pred in predictions]
        
        table_data = [[smile] + [all_predictions[model][i] for model in model_names] for i, smile in enumerate(smiles_list)]
        headers = ['SMILES'] + [f'Prediction ({model})' for model in model_names]
        return render_template('index.html', models=available_models, headers=headers, data=table_data)
    except Exception as e:
        print(f"Exception occurred: {e}")
        traceback.print_exc()
        return render_template('index.html', models=available_models, error="An error occurred while processing the request.")


@app.route('/api', methods=['POST'])
def apipredict():
    data = request.json
    smiles = data.get('smiles', [])
    models = data.get('models', [])
    output_format = data.get('format', 'json')  # Default to JSON format

    if not smiles or not isinstance(smiles, list):
        return jsonify({'error': 'Invalid input: please provide a list of SMILES strings.'}), 400

    if not models or not isinstance(models, list):
        return jsonify({'error': 'Invalid input: please provide a list of model names.'}), 400

    all_predictions = {}

    for model_name in models:
        model_path = os.path.join(MODELS_DIR, model_name, f"{model_name}_meta.json")
        if not os.path.exists(model_path):
            return jsonify({'error': f"Model {model_name} does not exist."}), 400

        model = SklearnModel.fromFile(model_path)
        predictions = model.predictMols(smiles)
        predictions_formatted = [f"{pred[0]:.4f}" for pred in predictions]
        all_predictions[model_name] = predictions_formatted

    # Format the result
    result = []
    for i, smile in enumerate(smiles):
        result_entry = {'smiles': smile}
        for model_name in models:
            result_entry[f'prediction ({model_name})'] = all_predictions[model_name][i]
        result.append(result_entry)

    if output_format == 'text':
        result_text = '\n'.join([f"SMILES: {entry['smiles']} -> " + ", ".join([f"{model}: {pred}" for model, pred in entry.items() if model != 'smiles']) for entry in result])
        return Response(result_text, mimetype='text/plain')

    if output_format == 'csv':
        # Create an in-memory output file for CSV data
        output = io.StringIO()
        writer = csv.DictWriter(output, fieldnames=result[0].keys())
        writer.writeheader()
        writer.writerows(result)
        output.seek(0)
        return Response(output.getvalue(), mimetype='text/csv', headers={"Content-Disposition": "attachment;filename=predictions.csv"})


    return jsonify(result)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)

