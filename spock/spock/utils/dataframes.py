import os.path
import warnings

import pandas as pd


def load_dataframe(path: str, **kwargs):
    """
    Load a library file into a pandas dataframe.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f'File not found: {path}')
    extension = path.split('.')[-1]
    if extension == 'xlsx':
        df = pd.read_excel(path, **kwargs)
    elif extension == 'csv':
        df = pd.read_csv(path, sep=',', **kwargs)
    elif extension == 'tsv':
        df = pd.read_csv(path, sep='\t', **kwargs)
    else:
        raise ValueError(f'Unsupported file format: .{extension}')
    return df