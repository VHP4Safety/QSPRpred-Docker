from spock.prep import Ligprep
from multiprocessing import Pool, current_process, Queue
import os
from spock.utils.logs.logs import read_log

############################################
## Deprecated, use spock.prep.ligprep instead
############################################
def ligprep_job(smiles_path: str):
    """
    Run ligprep on a single SMILES string.
    """
    ident = current_process().ident
    device_id = queue.get()
        
    if verbosity:
        print(f'{ident}: starting ligprep process on {device_id}')

    try:
        ligprep_kwargs['smiles'] = smiles_path
        _ = runners[device_id].run(**ligprep_kwargs)
        print('{}: finished'.format(ident))
    except Exception as e:
        print(e)
    finally:
        queue.put(device_id)

def parallel_ligprep(smiles=[], smi_paths=[], output_subfolder='',
                     num_workers=1):
    global verbosity
    global runners
    global queue
    global ligprep_kwargs
    verbosity = locals['verbose']
    runners = [Ligprep() for _ in range(num_workers)]
    queue = Queue()
    ligprep_kwargs = locals()

    # initialize the queue with the GPU ids
    for process_id in range(num_workers):
        queue.put(process_id)

    # Start the worker pool
    pool = Pool(processes=num_workers)
    for _ in pool.imap_unordered(ligprep_job, smiles):
        pass
    pool.close()
    pool.join()

    ## Read generated scores from the log file
    path = os.path.join('output', output_subfolder, 'log.tsv')
    log = read_log(path)
    scores = []
    processed_smiles = [entry[0] for entry in log]
    for ligand in smiles:
        if ligand in processed_smiles:
            idx = processed_smiles.index(ligand)
            best_score = log[idx][2][0]
            scores.append(best_score)
        else:
            scores.append(100.0) # Arbitrarily high score for failed docking

    return scores