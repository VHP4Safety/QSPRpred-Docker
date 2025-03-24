import base64
import csv
import io
import json
import logging
import os

import numpy as np
import pandas as pd
from flask import Flask, Response, jsonify, render_template, request, send_file
from flask_cors import CORS
from qsprpred.models import SklearnModel
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

from qprf.qprf_utils import render_qprf

app = Flask(__name__)
app.jinja_env.filters['zip'] = zip
CORS(app)  # Allow all origins by default

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Define the models directory
MODELS_DIR = 'models'

# define RDKit image implementer
def smiles_to_image(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles) # attempt conversion to RDKit molecule
        if mol is None: # return None if not possible
            return None
        img = Draw.MolToImage(mol) # Image generation
        buffer = io.BytesIO()
        img.save(buffer, format="PNG")
        buffer.seek(0)
        img_base64 = base64.b64encode(buffer.read()).decode('utf-8')
        buffer.close()
        return f"data:image/png;base64,{img_base64}"
    except Exception as e:
        logging.error(f"Error generating image for SMILES {smiles}: {e}") # Log the error message if any exception occurs during the process.
        return None

# define invalid SMILES scrubber
def validate_smiles(smiles_list):
    """
    Validates a list of SMILES strings. Returns a list of invalid SMILES.
    """
    invalid_smiles = []
    for smile in smiles_list:
        if Chem.MolFromSmiles(smile) is None:
            invalid_smiles.append(smile)
    return invalid_smiles

class SimilaritySearcher():
    def __init__(self):
        self.descgen = AllChem.GetMorganGenerator(radius=3)
        self.scorer = DataStructs.BulkTanimotoSimilarity
        
    def get_nearest_neighbors(self, smile, ms):
        """_summary_

        Args:
            smile (str): smiles string
            ms (list): list of rdkit molecules from reference set

        Returns:
            id for most similar molecule in reference set
        """
        m1 = Chem.MolFromSmiles(smile)
        query_fp = self.descgen.GetSparseCountFingerprint(m1)

        target_fingerprints = [self.descgen.GetSparseCountFingerprint(x) for x in ms]
        scores = self.scorer(query_fp, target_fingerprints)

        id_top = np.argmax(np.array(scores))
        
        return id_top
    
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
                    'pref_name': state['pref_name'],
                    'case_study': state['case_study'],
                    'target_property_name': state['targetProperties'][0]['py/state']['name'],
                    'target_property_task': state['targetProperties'][0]['py/state']['task']['py/reduce'][1]['py/tuple'][0],
                    'feature_calculator': state['featureCalculators'][0]['py/object'].split('.')[-1],
                    'radius': state['featureCalculators'][0]['py/state']['radius'],
                    'nBits': state['featureCalculators'][0]['py/state']['nBits'],
                    'algorithm': state['alg'].split('.')[-1],
                    'currDir': os.path.join(directory, d),
                }
                models_info.append(model_info)
                logging.info(f"Loaded model metadata: {model_info}")
    return models_info

@app.route('/download')
def download_qmrf():
    path = request.args.get('path')
    return send_file(path + '/qmrf.docx', as_attachment=True)

@app.route('/downloadqprf')
def download_qprf():
    model = request.args.get('model')
    smile = request.args.get('smile')
    return send_file(f"qprf/output/{model}/{smile}.docx", as_attachment=True)

@app.route('/')
@app.route('/predict')
def home():
    available_models = extract_model_info(MODELS_DIR)
    return render_template('index.html', models=available_models)

@app.route('/predict', methods=['POST'])
def predict():
    logging.info("Handling prediction request.")
    available_models = extract_model_info(MODELS_DIR)
    try:
        smiles_input = request.form.get('smiles')
        uploaded_file = request.files.get('file')
        model_names = request.form.getlist('model')
        file_name = request.form.get('uploaded_file_name')
        
        logging.debug(f"Received SMILES input: {smiles_input}")
        logging.debug(f"Uploaded file: {uploaded_file}")
        logging.debug(f"Selected models: {model_names}")
        logging.debug(f"Previous uploaded file name: {file_name}")
        
        if not model_names:
            logging.error("No model selected.")
            return render_template('index.html', models=available_models, error="No model selected.")
        
        smiles_list = []
        invalid_smiles = []

        # Handle SMILES string input
        if smiles_input:
            input_smiles = [smile.strip() for smile in smiles_input.split(',')]
            
            # Check if only one SMILES string is entered
            if len(input_smiles) == 1:
                if Chem.MolFromSmiles(input_smiles[0]) is None:  # Check for invalid single SMILES
                    logging.error(f"Invalid SMILES string: {input_smiles[0]}")  # Log the invalid SMILES
                    return render_template('index.html', models=available_models, error="Invalid SMILES string")  # Display error for single invalid SMILES
                else:
                    smiles_list.extend(input_smiles)  # Add valid SMILES to processing list
            else:
                invalid_smiles.extend([smile for smile in input_smiles if Chem.MolFromSmiles(smile) is None])  # Collect invalid SMILES
                smiles_list.extend([smile for smile in input_smiles if Chem.MolFromSmiles(smile) is not None])  # Collect valid SMILES
        
        # Handle uploaded file
        if uploaded_file and uploaded_file.filename != '':
            file_name = uploaded_file.filename
            logging.debug("Processing uploaded file.")
            uploaded_df = pd.read_csv(uploaded_file)
            logging.debug(f"Uploaded file contents: {uploaded_df.head()}")
            if 'SMILES' in uploaded_df.columns:
                file_smiles = uploaded_df['SMILES'].tolist()
                invalid_smiles.extend([smile for smile in file_smiles if Chem.MolFromSmiles(smile) is None])  # Collect invalid SMILES from file
                smiles_list.extend([smile for smile in file_smiles if Chem.MolFromSmiles(smile) is not None])  # Collect valid SMILES from file
        elif file_name:
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], file_name)
            if os.path.exists(file_path):
                logging.debug("Reprocessing previous uploaded file.")
                uploaded_df = pd.read_csv(file_path)
                if 'SMILES' in uploaded_df.columns:
                    file_smiles = uploaded_df['SMILES'].tolist()
                    invalid_smiles.extend([smile for smile in file_smiles if Chem.MolFromSmiles(smile) is None])  # Collect invalid SMILES from previous file
                    smiles_list.extend([smile for smile in file_smiles if Chem.MolFromSmiles(smile) is not None])  # Collect valid SMILES from previous file
        
        logging.debug(f"Final SMILES list: {smiles_list}")
        logging.debug(f"Invalid SMILES detected: {invalid_smiles}")  # Log invalid SMILES
        
        if not smiles_list and not invalid_smiles:
            error_message = "No SMILES strings provided"
            logging.error(error_message)
            return render_template('index.html', models=available_models, error=error_message)
        
        all_predictions = {}
        all_ads = {}
        model_info_list = []
        for model_name in model_names:
            logging.debug(f"Processing model: {model_name}")
            model_path = os.path.join(MODELS_DIR, model_name, f"{model_name}_meta.json")
            model = SklearnModel.fromFile(model_path)
            ad = []
            if getattr(model, 'applicabilityDomain', None):
                predictions, ad = model.predictMols(smiles_list, use_applicability_domain=True)
                ad = list(ad)
            else:
                predictions = model.predictMols(smiles_list)
            
            # Check if the model is regression or classification
            if model.task.isRegression():
                formatted_predictions = [f"{pred[0]:.2f}" for pred in predictions]
            else:
                # Format classification output as Active/Inactive
                formatted_predictions = ["Active" if pred[0] == 1 else "Inactive" for pred in predictions]

            all_predictions[model_name] = formatted_predictions
            all_ads[model_name] = ad

            
            # Extract model info for report
            with open(model_path, 'r') as meta_file:
                meta_data = json.load(meta_file)
                state = meta_data['py/state']
                model_info = {
                    'name': state['name'],
                    'pref_name': state['pref_name'],
                    'case_study': state['case_study'],
                    'target_property_name': state['targetProperties'][0]['py/state']['name'],
                    'target_property_task': state['targetProperties'][0]['py/state']['task']['py/reduce'][1]['py/tuple'][0],
                    'feature_calculator': state['featureCalculators'][0]['py/object'].split('.')[-1],
                    'radius': state['featureCalculators'][0]['py/state']['radius'],
                    'nBits': state['featureCalculators'][0]['py/state']['nBits'],
                    'algorithm': state['alg'].split('.')[-1]
                }
                model_info_list.append(model_info)
        
        table_data = []
        
        for i, smile in enumerate(smiles_list): 
            image_data = smiles_to_image(smile)
            if getattr(model, 'applicabilityDomain', None):
                row = [image_data] + [smile] + [all_predictions[model][i] + f' ({str(all_ads[model][i])})' for model in model_names]
            else:
                row = [image_data] + [smile] + [all_predictions[model][i] for model in model_names]

            table_data.append(row)
                        
        # Update headers
        table_data_extensive = []
        headers = ['Structure', 'SMILES']
        tooltips = ['2D depiction of input molecule', 'Line representation of input molecule']
        headers_extensive = ['Model', 'Structure', 'SMILES', 'Nearest Neighbor', 'Source', 'Predicted pChEMBL Value', 'Within Applicability Domain']
        tooltips_extensive = [
            'Unique identifier of the model that made the prediction', 
            '2D depiction of input molecule', 
            'Line representation of input molecule', 
            '2D depiction of nearest neighbor of input molecule in model training set. More information available in QPRF', 
            'Document(s) containing experimental data for nearest neighbor',
            'Model prediction for input molecule. pChEMBL is defined as -log(response). More information available in QMRF & QPRF',
            'AD is based on descriptors of training set. An input molecule is within AD if the distance to the training set is lower than a set threshold. More information available in QMRF & QPRF',
            ]
        searcher = SimilaritySearcher()
        for model_name in model_names:
            accession = model_name.split("_")[0]
            train_df = pd.read_csv(f'data/{accession}_Data/train_full_model_{accession}.csv').reset_index()
            train_smiles = train_df['SMILES'].to_list()
            ms = [Chem.MolFromSmiles(x) for x in train_smiles]
            model_path = os.path.join(MODELS_DIR, model_name, f"{model_name}_meta.json")               
            model = SklearnModel.fromFile(model_path)
            
            if getattr(model, 'applicabilityDomain', None):
                if model.task.isRegression():
                    # Format regression table header
                    headers.append(f'Predicted pChEMBL Value ({model_name})')
                    tooltips.append('Model prediction for input molecule. pChEMBL is defined as -log(response). The value in brackets shows if input molecule is within AD')
                else:
                    # Format classification table header
                    headers.append(f'Predicted class label ({model_name})')
                    tooltips.append('Model prediction for input molecule. The value in brackets shows if input molecule is within AD')
            else:
                if model.task.isRegression():
                    # Format regression table header
                    headers.append(f'Predicted pChEMBL Value ({model_name})')
                    tooltips.append('Model prediction for input molecule. pChEMBL is defined as -log(response)')
                else:
                    # Format classification table header
                    headers.append(f'Predicted class label ({model_name})')
                    tooltips.append('Model prediction for input molecule')

            
            for i, smile in enumerate(smiles_list): 
                image_data = smiles_to_image(smile)
                id_top = searcher.get_nearest_neighbors(smile, ms)
                nearest_neighbor = {}
                nn_smiles = train_df.iloc[id_top]['SMILES']
                nearest_neighbor["smiles"] = nn_smiles
                doi_nn = train_df.iloc[id_top]['doi']
                if doi_nn:
                    doi_nn = 'https://doi.org/' + doi_nn
                else:
                    doi_nn = train_df.iloc[id_top]['all_doc_ids']
                nearest_neighbor["reference"] = doi_nn
                nearest_neighbor["value"] = train_df.iloc[id_top]['pchembl_value']
                nearest_neighbor["predicted_value"] = model.predictMols([nn_smiles])[0][0]
                nearest_neighbor["similarity"] = f"Nearest neighbor was found using {searcher.scorer.__name__} based on {searcher.descgen.__class__.__name__}"
                image_data_nn = smiles_to_image(nn_smiles)
                if getattr(model, 'applicabilityDomain', None):
                    row = [model_name] + [image_data] + [smile] + [image_data_nn] + [nn_smiles] + [doi_nn] + [all_predictions[model_name][i]] + [all_ads[model_name][i]]
                else:
                    row = [model_name] + [image_data] + [smile] + [image_data_nn] + [nn_smiles] + [doi_nn] + [all_predictions[model_name][i]]

                table_data_extensive.append(row)

                render_qprf(smile, model, predictions[i], ad[i], nearest_neighbor)

        error_message = None
        if invalid_smiles:
            error_message = f"Invalid SMILES, could not be processed: {', '.join(invalid_smiles)}"  # Mention invalid SMILES in error message
        
        return render_template('index.html', models=available_models, headers=headers, tooltips=tooltips, data=table_data, headers_extensive=headers_extensive, tooltips_extensive = tooltips_extensive, data_extensive=table_data_extensive, smiles_input=smiles_input, model_names=model_names, file_name=file_name, error=error_message)
    except Exception:
        logging.exception("An error occurred while processing the request.")
        return render_template('index.html', models=available_models, error="An error occurred while processing the request.")

def create_report(model_info_list, headers, table_data):
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter)
    elements = []

    # Add title
    styles = getSampleStyleSheet()
    title_style = styles['Title']
    elements.append(Paragraph("Prediction Report", title_style))
    elements.append(Spacer(1, 12))

    # Add model metadata
    for model_info in model_info_list:
        elements.append(Paragraph(f"Model: {model_info['name']}", styles['Heading2']))
        elements.append(Paragraph(f"Model Name: {model_info['pref_name']}", styles['Heading2']))
        elements.append(Paragraph(f"Target Property Name: {model_info['target_property_name']}", styles['Normal']))
        elements.append(Paragraph(f"Target Property Task: {model_info['target_property_task']}", styles['Normal']))
        elements.append(Paragraph(f"Feature Calculator: {model_info['feature_calculator']}", styles['Normal']))
        elements.append(Paragraph(f"Radius: {model_info['radius']}", styles['Normal']))
        elements.append(Paragraph(f"nBits: {model_info['nBits']}", styles['Normal']))
        elements.append(Paragraph(f"Algorithm: {model_info['algorithm']}", styles['Normal']))
        elements.append(Spacer(1, 12))

    # Convert headers to Paragraphs for wrapping
    header_paragraphs = [Paragraph(header, styles['Normal']) for header in headers]

    # Convert SMILES strings to Paragraphs for wrapping
    for i in range(len(table_data)):
        table_data[i][0] = Paragraph(table_data[i][0], styles['Normal'])

    # Add prediction table
    data = [header_paragraphs] + table_data
    table = Table(data, repeatRows=1)

    # Apply style to table
    table.setStyle(TableStyle([
#        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
#        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 12),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
#        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
    ]))

    # Ensure SMILES strings can wrap
    table._argW[0] = 2.5 * inch  # width of SMILES column
    elements.append(table)

    doc.build(elements)
    buffer.seek(0)
    return buffer

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
