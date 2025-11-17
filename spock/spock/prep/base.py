from abc import ABC, abstractmethod
import os
import sys
from dataclasses import dataclass
from spock.storage import TabularStorage


class LigandPreparation(ABC):
    def __init__(self, store: TabularStorage=None):
        self.store = store

    @abstractmethod
    def run(self, sdf_path: str, output_folder: str) -> None:
        pass

    @abstractmethod
    def process_store(self, store: TabularStorage, output_folder: str) -> None:
        pass


