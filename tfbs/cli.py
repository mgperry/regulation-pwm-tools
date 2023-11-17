import click
import sys

from pyfaidx import Fasta
from pybedtools import BedTool

from .scanner import Scanner
from .readers import load_pwms


@click.command()
@click.option("--pfms", required=True)
@click.option("--fasta", required=True)
@click.option("--bed", help="Restrict calling to given ranges.", type=str)
@click.option("--pvalue", type=float, default=0.0001)
@click.option("--cores", default=1, type=int)
def scan(pfms: str, fasta: str, bed: str, pvalue: float, cores: int):
    if cores > 1:
        print("multicore not yet available", file=sys.stderr)

    pwms = load_pwms(pfms, pvalue)

    fa = Fasta(fasta)

    if bed:
        seqs = (fa[i.chrom][(i.start + 1):i.end] for i in BedTool(bed))
    else:
        seqs = (record[:] for record in fa)

    s = Scanner(pwms)

    hits_by_seq = (s.scan(seq) for seq in seqs)

    for hits in hits_by_seq:
        for hit in hits:
            print(hit.to_bed())
