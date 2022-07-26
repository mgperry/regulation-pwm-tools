import re
import pydantic
from typing import Iterator, TextIO
from itertools import zip_longest

class PFM(pydantic.BaseModel):

    id: str
    metadata = {}
    PFM: list[list[int]]

    @pydantic.validator("PFM")
    @classmethod
    def validate_PFM(cls, rows):
        
        if len(rows) != 4:
            raise ValueError("PFM should have 4 rows")

        if len(set(len(row) for row in rows)) > 1:
            raise ValueError("All PFM rows should be the same length")
        
        return rows

def parse_integer_matrix(rows: list[str]) -> list[list[int]]:
    """
    Parse a set of matrix rows stripping out non-numeric values.
    """

    if len(rows) != 4:
        raise ValueError("PFM should have 4 rows")

    int_rows = [extract_ints(row) for row in rows]

    if len({len(row) for row in int_rows}) > 1:
        raise ValueError("All PFM rows should be the same length")
        
    return int_rows


def pfm_reader(f: TextIO, header_start=">") -> Iterator[dict]:
    """"
    Return an iterator to read in JASPAR style matrix files (header following by 4 rows)
    """
    return (parse_jaspar_entry(lines, header_start) for lines in grouper(5, f))


def parse_jaspar_entry(lines: list[str], header_start=">"):
    header = lines[0].rstrip()
        
    if header.startswith(header_start):
        header = header[1:]

    matrix = parse_integer_matrix(lines[1:])

    return {
        "header": header,
        "PFM": matrix
    }

def extract_ints(s: str):
    """Extract and convert integers in a string."""
    return [int(x) for x in re.findall("\d+", s)]


def grouper(n, iterable):
    "grouper(3, 'ABCDEFG') --> ABC DEF"
    args = [iter(iterable)] * n
    return zip_longest(*args)
