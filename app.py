from flask import Flask, request, render_template, jsonify, Response
from qsprpred.models import SklearnModel
import os
import json
import traceback
import pandas as pd

app = Flask(__name__)

# Define the models directory
MODELS_DIR = 'models'

# Load the default model when the application starts
model_name = 'models/P10828_RF_Model/P10828_RF_Model_meta.json'
model = SklearnModel.fromFile(model_name)

@app.route('/')
def home():
    # List available models
    available_models = [d for d in os.listdir(MODELS_DIR) if os.path.isdir(os.path.join(MODELS_DIR, d))]
    return render_template('index.html', models=available_models)

@app.route('/predict', methods=['POST'])
def predict():
    try:
        model_name = request.form.get('model')
        smiles = request.form.get('smiles')
        file = request.files.get('file')

        if not model_name:
            return render_template('index.html', error="No model selected.")

        model_path = os.path.join(MODELS_DIR, model_name, f"{model_name}_meta.json")
        
        if not os.path.exists(model_path):
            return render_template('index.html', error="Selected model does not exist.")
        
        model = SklearnModel.fromFile(model_path)

        if file:
            # Handle file upload
            if file.filename.endswith('.csv'):
                df = pd.read_csv(file)
                smiles_list = df.iloc[:, 0].dropna().tolist()  # Assuming SMILES are in the first column
            else:
                return render_template('index.html', error="Invalid file type. Please upload a CSV file.")
        else:
            # Handle SMILES input
            if not smiles:
                return render_template('index.html', error="No SMILES strings provided.")
            smiles_list = [smile.strip() for smile in smiles.split(',')]

        # Perform prediction
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

@app.route('/api', methods=['POST'])
def apipredict():
    try:
        data = request.json
        smiles = data.get('smiles', [])
        output_format = data.get('format', 'json')  # Default to JSON format

        if not smiles or not isinstance(smiles, list):
            return jsonify({'error': 'Invalid input: please provide a list of SMILES strings.'}), 400

        # Perform prediction
        predictions = model.predictMols(smiles)

        # Create result
        result = [{'smiles': sm, 'prediction': round(pred[0], 4)} for sm, pred in zip(smiles, predictions)]

        if output_format == 'text':
            # Format result as plain text
            result_text = '\n'.join([f"{sm},{pred}" for sm, pred in zip(smiles, predictions)])
            return Response(result_text, mimetype='text/plain')

        if output_format == 'csv':
        # Create CSV output
            output = io.StringIO()
            writer = csv.writer(output)
            writer.writerow(['SMILES', 'Predicted pChEMBL value'])
            for item in result:
                writer.writerow([item['smiles'], f"{item['prediction']:.4f}"])
            output.seek(0)
            return Response(output.getvalue(), mimetype='text/csv', headers={'Content-Disposition': 'attachment;filename=predictions.csv'})


        # Default to JSON format
        return jsonify(result)
    except Exception as e:
        print(f"Exception occurred: {e}")
        traceback.print_exc()
        return jsonify({'error': 'An error occurred while processing the request.'}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)  # Enable debug mode

