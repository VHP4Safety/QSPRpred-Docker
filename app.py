from flask import Flask, request, render_template
from qsprpred.models import SklearnModel
import traceback

app = Flask(__name__)

model_name = 'models/P10828_RF_Model/P10828_RF_Model_meta.json'
model = SklearnModel.fromFile(model_name)

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/predict', methods=['POST'])
def predict():
    try:
        smiles = request.form.get('smiles')
        if not smiles:
            return render_template('index.html', error="No SMILES strings provided.")
        
        smiles_list = [smile.strip() for smile in smiles.split(',')]
        predictions = model.predictMols(smiles_list)
        predictions_formatted = [f"{pred[0]:.4f}" for pred in predictions]
        
        # Prepare data as a list of tuples
        data = list(zip(smiles_list, predictions_formatted))
        
        return render_template('index.html', data=data)
    except Exception as e:
        # Print the exception to the Flask server logs
        print(f"Exception occurred: {e}")
        traceback.print_exc()
        return render_template('index.html', error="An error occurred while processing the request.")

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)  # Enable debug mode

