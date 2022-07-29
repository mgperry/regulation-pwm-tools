import json
from multiprocessing import Pool
import click

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
@click.option("--cores", default=4, type=int)
def moods(pfms: str, fasta: str, bed: str, pvalue: float, cores: int):
    pfms = json.load(open(pfms, "r"))
    pwms = [PWM(p["PFM"], p["name"], pvalue=pvalue) for p in pfms]

    fa = Fasta(fasta)

    if bed:
        seqs = (fa[i.chrom][(i.start + 1):i.end] for i in BedTool(bed))
    else:
        seqs = (record[:] for record in fa)

    s = Scanner(pwms)

    def scan(fa: tuple[str, str]):
        return s.scan(fa)

    results = list(map(scan, seqs))

    print(json.dumps(results, indent=4)) # TODO better than JSON results?

@click.command()
@click.argument('input')
def jaspar_to_json(input: str):
    pfms = list(read_jaspar(input))
    print(json.dumps(pfms, indent=4, cls=NumpyEncoder))
    return 0


@click.command()
@click.argument('input')
def moods_to_json(input: str):
    pfms = list(read_moods(input))
    print(json.dumps(pfms, indent=4, cls=NumpyEncoder))
    return 0


cli.add_command(moods)
cli.add_command(jaspar_to_json)
cli.add_command(moods_to_json)
