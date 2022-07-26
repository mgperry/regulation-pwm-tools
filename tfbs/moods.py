import argparse
import json
import MOODS.parsers
import MOODS.scan
import MOODS.tools
import numpy as np
from multiprocessing import Pool
from collections import namedtuple

from .pwm import PWM
from .utils import iter_fasta


Hit = namedtuple("Hit", ["seqname", "TF", "start", "end", "score", "strand"])


class Scanner:
    def __init__(self, pwms, background=(0.25, 0.25, 0.25, 0.25)):

        matrices = []
        thresholds = []

        for pwm in pwms:
            matrices.extend([pwm.PWM, np.flip(pwm.PWM, axis=[0, 1])])
            thresholds.extend([pwm.threshold] * 2)

        scanner = MOODS.scan.Scanner(7)  # why 7 lol
        scanner.set_motifs(matrices, background, thresholds)

        self.pwms = pwms
        self.scanner = scanner
        self.background = background

    def scan(self, fa_record: tuple[str, str]):
        header, seq = fa_record

        raw = self.scanner.scan(seq)

        results = []

        for i, rs in enumerate(raw):
            strand = "+" if i % 2 == 1 else "-"  # odd pwms are reverse complement
            pwm_index = i // 2  # divide index by 2 due to rc pwms
            id = self.pwms[pwm_index].id
            width = self.pwms[pwm_index].width
            results.extend(
                [Hit(header, id, r.pos, r.pos + width, r.score, strand) for r in rs]
            )

        return results


parser = argparse.ArgumentParser()

parser.add_argument("--pfms", required=True)
parser.add_argument("--fasta", required=True)
parser.add_argument("--pvalue", type=float, default=0.0001)
parser.add_argument("--cores", default=4, type=int)


if __name__ == "__main__":

    args = parser.parse_args()

    pfms = json.load(open(args.pfms, "r"))
    pwms = [PWM(p["PFM"], p["label"], pvalue=args.pvalue) for p in pfms]

    seqs = iter_fasta(args.fasta)

    s = Scanner(pwms)

    def scan(fa: tuple[str, str]):
        return s.scan(fa)

    with Pool(args.cores) as p:
        results = list(p.map(scan, seqs))

    print(json.dumps(results, indent=4))
