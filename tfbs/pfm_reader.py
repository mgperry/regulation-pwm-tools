import re
import pydantic
from typing import List
from itertools import zip_longest

class PFM(pydantic.BaseModel):

    id: str
    metadata = {}
    PFM: List[List[int]]

    @pydantic.validator("PFM")
    @classmethod
    def validate_PFM(cls, rows):
        
        if len(rows) != 4:
            raise ValueError("PFM should have 4 rows")

        if len(set(len(row) for row in rows)) > 1:
            raise ValueError("All PFM rows should be the same length")
        
        return rows


def pfm_reader(f, delimiter="\s+"):
    for lines in grouper(5, f):
        header = lines[0].rstrip()
        if header.startswith(">"):
            header = header[1:]
        info = re.split(delimiter, header)

        matrix = [strip(re.split(delimiter, l)) for l in lines[1:]]

        yield info, matrix


def strip(xs):
    r = re.compile(r"\d+")
    return [int(x) for x in xs if r.match(x)]


def grouper(n, iterable):
    "grouper(3, 'ABCDEFG') --> ABC DEF"
    args = [iter(iterable)] * n
    return zip_longest(*args)
