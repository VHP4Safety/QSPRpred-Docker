import pymol
from pymol import cmd

class Pymol:
    def __init__(self):
        self.is_running = False

    def start(self):
        if not self.is_running:
            pymol.finish_launching(['pymol', '-qc'])
            self.is_running = True

    def finish(self):
        self.is_running = False
        pymol.cmd.quit()

    def load(self, path, clear=False):
        if clear:
            cmd.reinitialize()
        cmd.load(path)