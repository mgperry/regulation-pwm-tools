import json
from multiprocessing import Pool
import click

from pyfaidx import Fasta
from pybedtools import BedTool

from .scanner import Scanner
from .pwm import PWM

@click.command()
@click.option("--pfms", required=True)
@click.option("--fasta", required=True)
@click.option("--bed", help="Restrict calling to given ranges.", type=str)
@click.option("--pvalue", type=float, default=0.0001)
@click.option("--cores", default=4, type=int)
def cli(pfms: str, fasta: str, bed: str, pvalue: float, cores: int):
    pfms = json.load(open(pfms, "r"))
    pwms = [PWM(p["PFM"], p["name"], pvalue=pvalue) for p in pfms]

    fa = Fasta(fasta)

    if bed:
        seqs = (fa[i.chr][i.start + 1][i.end] for i in BedTool(bed))
    else:
        seqs = (record[:].seq for record in fa)

    s = Scanner(pwms)

    def scan(fa: tuple[str, str]):
        return s.scan(fa)

    with Pool(cores) as p:
        results = list(p.map(scan, seqs))

    print(json.dumps(results, indent=4)) # TODO better than JSON results?

