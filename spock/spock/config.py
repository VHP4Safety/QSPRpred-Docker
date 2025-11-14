from dataclasses import dataclass
from typing import List
import json

#### Not used for now, but will be in the future
#### One stop shop for all configuration parameters

@dataclass
class BaseConfig:

    def to_json(self):
        return json.dumps(self.__dict__, indent=4, sort_keys=True)

    def from_json(self, json_path):
        with open(json_path, 'r') as f:
            config = json.load(f)
        for key, value in config.items():
            setattr(self, key, value)

    def as_dict(self):
        """
        Return the configuration as a dictionary
        """
        return self.__dict__
    

@dataclass
class PathConfig(BaseConfig):
    
    ## Backends
    ligprep: str = 'rdkit' # rdkit, meeko, schrodinger, dimorphite
    target_prep: str = 'ADFR' # rdkit, schrodinger
    docking: str = 'vina'  # vina, vina-gpu, autodock, diffdock
    standardizer: str = 'papyrus' # papyrus, rdkit, schrodinger
    database: str = 'sqlite' # sqlite, postgresql

    ## Paths
    schrodinger_path: str = '/data/bernataviciusa/schrodinger2023'
    vina_gpu_path: str = '/vina-gpu-dockerized/vina'

    ## Vina-GPU arguments
    executable: str = 'Vina-GPU'
    help_flag: str = '--help'
    gpu_ids: int = 1
    workers_per_gpu: int = 1
    search_depth: int = 10
    threads: int = 1024
    threads_per_call: int = 1024


@dataclass
class LigprepParams(BaseConfig):
    schrodinger_path: str = PathConfig().schrodinger_path
    max_atoms: int = 500
    force_field: int = 16
    epik: str = 'epik'
    optimize: bool = True
    ph:  float = 6.5
    th:   float = 1.0
    i: int = 2
    delete: bool = False
    num_stereoisomers: int = 32
    cpus: int = 32
    workers: int = 1
    path: str = None
    executable: str = f'{schrodinger_path}/ligprep'
    help_flag: str = '--help'


@dataclass
class DefaultConfig(BaseConfig):
    """
    Main configuration class, contains all the configuration parameters
    Can be loaded/saved to restore state
    """
    paths: PathConfig = PathConfig()
    ligprep: LigprepParams = LigprepParams()
    # ligprep = SchrodingerLigprepConfig()
    # docking = VinaConfig()
    # standardizer = PapyrusConfig()

    def as_dict(self):
        """
        Return the configuration as a dictionary
        """
        dictionary = {}
        for key, value in self.__dict__.items():
            print(key, value)
            dictionary[key] = value.as_dict()
        return dictionary

    @classmethod
    def from_file(cls, path):
        """
        Load a configuration from a file
        """
        config = cls()
        config.from_json(path)
        return config


def main():
    config = DefaultConfig()
    print(config.as_dict())

if __name__ == '__main__':
    main()