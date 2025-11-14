import os
import tempfile

from spock.prep import Dimorphite
from spock.storage import TabularStorage


class TestStorageMixIn:

    def initStorage(self):
        # create random storage with random name
        self.chunk_size = 2
        self.n_jobs = 2
        self.storage = self.create_storage()
        self.prepare_storage()

    def create_storage(self):
        example_lib = os.path.join(
            os.path.dirname(__file__), "test_files/example_smiles.csv"
        )
        temp_path = tempfile.mkdtemp()
        store = TabularStorage(path=temp_path)
        store.add_library(
            example_lib,
            parallel=True,
            chunk_size=self.chunk_size,
            n_jobs=self.n_jobs,
            save=True,
        )
        store.summary()
        return store

    def prepare_storage(self):
        ligprep = Dimorphite()
        ligprep.process_store(self.storage, self.n_jobs)
        self.storage.summary()
