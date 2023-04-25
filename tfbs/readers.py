import re
from itertools import zip_longest
import numpy as np
from pathlib import Path
import json

from .motif import Motif


def load_pwms(jsn: str, pvalue: float, pfm="PFM", name="name", **kwargs) -> list[Motif]:
    pfms = json.load(open(jsn, "r"))
    return [Motif.from_pfm(p[pfm], p[name], **kwargs) for p in pfms]


def read_moods(d: str, suffix="*.pfm"):
    """
    Read in a directory of .pfm files, as required by MOODS (if not using this package).
    """
    d = Path(d)

    if not d.is_dir: raise Exception("read_moods() must be supplied a directory.")

    for pfm in d.glob(suffix):
        yield {"name": pfm.stem, "PFM": read_pfm_from_file(pfm)}


def read_pfm_from_file(f: str) -> np.ndarray:
    rows = open(f, 'r').readlines()
    return np.array([extract_ints(row) for row in rows])


def read_jaspar(f: str):
    """"
    Read in JASPAR style matrix files (header following by 4 rows).
    """
    with open(f, 'r') as fh:
        for rows in grouper(5, fh):
            header = rows[0].rstrip()[1:] # strip leading '>'
            matrix = np.array([extract_ints(row) for row in rows[1:]])
            yield {"name": header, "PFM": matrix}


def extract_ints(s: str):
    """Extract and convert integers in a string."""
    return [int(x) for x in re.findall(r"\d+", s)]


def grouper(n, iterable):
    "grouper(3, 'ABCDEFG') --> ABC DEF"
    args = [iter(iterable)] * n
    return zip_longest(*args)


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
