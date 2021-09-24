import sys
import json
import re
from itertools import zip_longest


def pfm_reader(f, delimiter="\s+"):
    for lines in grouper(5, f):
        header = lines[0].rstrip()
        if header[0] == ">":
            header = header[1:]
        info = re.split(delimiter, header)

        matrix = [strip(re.split(delimiter, l)) for l in lines[1:]]

        yield {"info": info, "PFM": matrix}


def strip(xs):
    r = re.compile(r"\d+")
    return [int(x) for x in xs if r.match(x)]


def grouper(n, iterable):
    "grouper(3, 'ABCDEFG') --> ABC DEF"
    args = [iter(iterable)] * n
    return zip_longest(*args)
