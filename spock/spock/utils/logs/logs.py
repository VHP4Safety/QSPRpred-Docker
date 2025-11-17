import os
from spock.utils import encode, decode

def write_to_log(log_path, smiles, target, scores, pdbqt_path=None):
    """
    Writes a log file
    Arguments:
        log_path (str)            : path to log file
        smiles (str)              : SMILES of ligand
        target (str)              : target name
        scores (list)             : list of scores
        pdbqt_path (str)          : path to pdbqt file
    """

    if not os.path.isfile(log_path):
        with open(os.path.join(log_path), 'w') as f:
            header = '\t'.join(['smiles', 'target', 'scores', 'pdbqt'])
            f.write(header + '\n')

    if pdbqt_path is not None:
        with open(pdbqt_path, 'r') as f:
            pdbqt = f.read()
        pdbqt = encode(pdbqt)
    else:
        pdbqt = ''
    
    if not isinstance(scores, list):
        scores = [scores]
    
    z = [str(score) for score in  scores]
    if len(z) == 1:
        scores = z[0]
    else:
        scores = ';'.join(z)

    with open(log_path, 'a') as f:
        f.write('\t'.join([smiles, target, scores, pdbqt])+'\n')
    

def read_log(log_path):
    """
    Reads a log file
    Arguments:
        log_path (str)            : path to log file
    Returns:
        log (list)                : list of log entries
    """
    log = []
    with open(log_path, 'r') as f:
        lines = f.readlines()[1:]
        for line in lines:
            smiles, target, scores, pdbqt = line.strip().split('\t')
            scores = [float(score) for score in scores.split(';')]
            pdbqt = decode(pdbqt)
            log += [(smiles, target, scores, pdbqt)]
    return log