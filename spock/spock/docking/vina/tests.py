import os
import unittest

from spock.storage.tests.utils import TestStorageMixIn
from spock.docking.vina.cpu_local import VinaDockingCPULocal, VinaProtein


class ProteinMixIn:

    def initProtein(self):
        self.pdb_file = os.path.join(
            os.path.dirname(__file__),
            "test_files/5T1A_clean_mutations_reversed_withHs.pdb",
        )
        self.pdbqt_file = os.path.join(
            os.path.dirname(__file__),
            "test_files/5T1A_clean_mutations_reversed_withHs.pdbqt",
        )
        self.protein = VinaProtein(
            "CCR2",  # just a name
            self.pdb_file,
            self.pdbqt_file,
        )
        self.box_spec = {"center": [5.1, 28.0, 187.6], "box_size": [16.2, 17.8, 17.4]}


class DockingTestCase(TestStorageMixIn, ProteinMixIn, unittest.TestCase):

    def setUp(self):
        super().setUp()
        self.initProtein()
        self.initStorage()
        self.docking_jobs = 2
        self.docking_chunk_size = 1

    def check_poses(self, n_expected=9):
        self.assertEqual(len(self.storage.targets), 1)
        for mol in self.storage:
            poses = self.storage.get_poses(mol.id, self.storage.targets[0])
            self.assertTrue(len(poses) == n_expected)

    def test_basic_local_docking(self):
        # dock storage
        docking = VinaDockingCPULocal(
            protein=self.protein,
            n_cpus=self.docking_jobs,
            box_spec=self.box_spec,
            embed_mols=True,
            exhaustiveness=8,
            seed=42,
        )
        docking.dock_storage(
            self.storage, chunk_size=self.docking_chunk_size, overwrite=True, save=True
        )
        self.check_poses(n_expected=9)
        # add poses
        docking.n_cpus = 1
        docking.dock_storage(
            self.storage, chunk_size=len(self.storage), overwrite=False, save=True
        )
        self.storage.summary()
        self.check_poses(n_expected=18)
        # use timeout
        docking.timeout = 1
        docking.n_cpus = 2
        docking.dock_storage(
            self.storage, chunk_size=self.chunk_size, overwrite=False, save=True
        )
        self.check_poses(n_expected=18)
