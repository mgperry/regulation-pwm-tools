#%%
import json
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import polars as pl

from ncls import NCLS
from pyfaidx import Fasta
from tfbs import Scanner, Motif

#%%
@dataclass
class Range:
    chr: str
    start: int
    end: int
    name: str = "."
    score: float = 0
    strand: str = "."


def read_bed(f):
    return pl.read_csv(
        f,
        has_header=False,
        separator="\t",
        new_columns = ["chr", "start", "end", "name", "score", "strand"],
        dtypes = {
            "start": pl.Int32,
            "end": pl.Int32
        }
    )


#%%
hg38 = Fasta(Path.home() / "Data/Annotation/ENCODE_GRCh38.fasta")

peak_df = read_bed("test_data/k562_atac_peaks_chr11.bed")

def get_seq(r: Range):
    return hg38[r.chr][r.start:r.end]


#%%
with open("test_data/motifs.json") as jsn:
    motifs = [
        Motif.with_p(tf["id"], tf["matrix"], p=0.0001)
        for tf in json.load(jsn)
    ]


# %%
scanner = Scanner(motifs)

motif_hits = {}

peaks = (Range(*peak[:6]) for peak in peak_df.rows())

for peak in peaks:
    seq = get_seq(peak)
    id = f"{peak.chr}:{peak.start}:{peak.end}"
    
    hits = {}

    for h in scanner.scan(seq):
        if not h.name in hits:
            hits[h.name] = []
        hits[h.name].append(h)

    motif_hits[id] = hits

