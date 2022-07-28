import re
from itertools import zip_longest
import numpy as np


def read_jaspar(f: str):
    """"
    Read in JASPAR style matrix files (header following by 4 rows)
    """
    with open(f, 'r') as fh:
        for rows in grouper(5, fh):
            header = rows[0].rstrip()[1:] # strip leading '>'
            matrix = np.array([extract_ints(row) for row in rows[1:]])
            yield (header, matrix)


def extract_ints(s: str):
    """Extract and convert integers in a string."""
    return [int(x) for x in re.findall("\d+", s)]


def grouper(n, iterable):
    "grouper(3, 'ABCDEFG') --> ABC DEF"
    args = [iter(iterable)] * n
    return zip_longest(*args)
