import json
from multiprocessing import Pool
import click

from pyfaidx import Fasta

from .scanner import Scanner
from .pwm import PWM

@click.command()
@click.argument("--pfms", required=True)
@click.argument("--fasta", required=True)
@click.argument("--pvalue", type=float, default=0.0001)
@click.argument("--cores", default=4, type=int)
def cli(pfms: str, fasta: str, pvalue: float, cores: int):
    pfms = json.load(open(pfms, "r"))
    pwms = [PWM(p["PFM"], p["label"], pvalue=pvalue) for p in pfms]

    seqs = (record[:].seq for record in Fasta(fasta))

    s = Scanner(pwms)

    def scan(fa: tuple[str, str]):
        return s.scan(fa)

    with Pool(cores) as p:
        results = list(p.map(scan, seqs))

    print(json.dumps(results, indent=4)) # TODO better than JSON results?

