import json
import click
import sys

from pyfaidx import Fasta
from pybedtools import BedTool

from .scanner import Scanner
from .pwm import PWM
from .io import read_jaspar, read_moods, NumpyEncoder

@click.group()
def cli():
    pass

@click.command()
@click.option("--pfms", required=True)
@click.option("--fasta", required=True)
@click.option("--bed", help="Restrict calling to given ranges.", type=str)
@click.option("--pvalue", type=float, default=0.0001)
@click.option("--cores", default=1, type=int)
def scan(pfms: str, fasta: str, bed: str, pvalue: float, cores: int):
    if cores > 1:
        print("multicore not yet available", file=sys.stderr)

    pfms = json.load(open(pfms, "r"))
    pwms = [PWM(p["PFM"], p["name"], pvalue=pvalue) for p in pfms]

    fa = Fasta(fasta)

    if bed:
        seqs = (fa[i.chrom][(i.start + 1):i.end] for i in BedTool(bed))
    else:
        seqs = (record[:] for record in fa)

    s = Scanner(pwms)

    hits_by_seq = (s.scan(seq) for seq in seqs)

    for hits in hits_by_seq:
        for hit in hits:
            bed = hit._replace(start = hit[1] - 1)
            print("%s\t%d\t%d\t%s\t%.2f\t%s" % bed)


readers = {
    "jaspar": read_jaspar,
    "moods": read_moods
}


@click.command()
@click.argument('format', type=click.Choice(['jaspar', 'moods'], case_sensitive=False))
@click.argument('input')
def parse(input: str, format: str):
    pfms = list(readers[format](input))
    print(json.dumps(pfms, indent=4, cls=NumpyEncoder))


cli.add_command(scan)
cli.add_command(parse)
