## Helper functions for merging/splitting sdf files\
import re

def sdf_to_lines(sdf_lines):
    previous_lines = []
    for line in sdf_lines:
        strip_number = re.compile(r"^(>  <\S+>)(\s+\(\d+\)\s*)$", re.MULTILINE)
        line = strip_number.sub(r"\1", line)
        if '$$$$' == line.strip():
            previous_lines.append(line)
            yield previous_lines
            previous_lines = []
        else:
            previous_lines.append(line)
    yield previous_lines  # Yield any remaining lines after the last "$$$$" line

def sdf_iterator(sdf_file):
    """
    Iterates over an sdf file, yielding each molecule as a list of lines
    """
    with open(sdf_file, 'r') as file:
        return sdf_to_lines(file.readlines())

def parse_sdf(raw_sdf):
    """
    Parses a raw sdf file into a dictionary of properties

    input: yield from sdf_iterator(sdf_file)
    output: {key: value}
    """
    sdf_dict = {}
    key = None
    for line in raw_sdf:
        if line.startswith('>'):
            key = line[3:].strip('<').strip('>').strip().strip('>')
        else:
            if key not in [None, ''] and line.strip() not in ['', '$$$$']:
                sdf_dict[key] = line.strip()
    return sdf_dict