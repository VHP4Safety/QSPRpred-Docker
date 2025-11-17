import subprocess as sp
import re


def run_executable(cmd, shell=True, **kwargs) -> tuple[bytes, bytes]:
    """
    Run executable command and return output from stdout and stderr
    """
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=shell, **kwargs)
    stdout, stderr = proc.communicate()
    return stdout, stderr
