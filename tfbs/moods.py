import argparse
import json
import numpy as np
import MOODS.parsers
import MOODS.scan
import MOODS.tools
import itertools
from multiprocessing import Pool
from itertools import groupby, chain, islice

# from dataclasses import dataclass
from collections import namedtuple


class PWM:
    def __init__(
        self,
        PFM,
        id,
        background=(0.25, 0.25, 0.25, 0.25),
        pseudocount=0.8,
        threshold=("pvalue", 0.001),
    ):
        self.PFM = np.array(PFM)
        self.id = id
        self.background = background
        self.pseudocount = pseudocount
        self._bg_col_vec = np.array([background]).T

        self.PWM = self.calculate_pwm()
        self.threshold = self.calculate_threshold(threshold)
        self.width = np.size(self.PWM, 1)

    def calculate_pwm(self):
        return np.log(self.PPM()) - np.log(self._bg_col_vec)

    def PPM(self):
        pfm_adj = self.PFM + (self._bg_col_vec * self.pseudocount)
        col_sums = np.sum(pfm_adj, axis=0)
        return pfm_adj / col_sums

    def tuples(self):
        return tuple(map(tuple, self.PWM))

    def tuples_rc(self):
        return tuple(map(tuple, np.flip(self.PWM, axis=[0, 1])))

    def calculate_threshold(self, threshold):
        if threshold is None:
            th = None
        elif type(threshold) == "float":
            th = threshold
        elif threshold[0] == "pvalue":
            # magic numbers copied from moods python script
            th = MOODS.tools.threshold_from_p_with_precision(
                self.PWM, self.background, threshold[1], 2000.0, 4
            )
        elif threshold[0] == "score":
            # manually calculate percentage of maximum score
            range = self.PWM.max(axis=0) - self.PWM.min(axis=0)
            th = np.sum(self.PWM.min(axis=0) + (threshold[1] * range))
        else:
            raise Exception("threshold not formatted correctly")

        return th


def iter_fasta(filename):
    with open(filename, "r") as f:
        iterator = groupby(f, lambda line: line[0] == ">")
        for is_header, group in iterator:
            if is_header:
                header = next(group)[1:].strip()
            else:
                yield header, "".join(s.strip() for s in group)


Hit = namedtuple("Hit", "TF start end score strand")


class Scanner:
    def __init__(self, pwms, background=(0.25, 0.25, 0.25, 0.25)):

        matrices = []
        thresholds = []

        for pwm in pwms:
            matrices.extend([pwm.tuples(), pwm.tuples_rc()])
            thresholds.extend([pwm.threshold] * 2)

        scanner = MOODS.scan.Scanner(7)  # why 7 lol
        scanner.set_motifs(matrices, background, thresholds)

        self.pwms = pwms
        self.scanner = scanner
        self.background = background

    def scan(self, seq: tuple[str, str]):
        results = {
            "header": seq[0],
            "hits": [],
        }

        raw = self.scanner.scan(seq[1])

        for i, rs in enumerate(raw):
            strand = "+" if i % 2 == 1 else "-"  # odd pwms are reverse complement
            pwm_index = i // 2 # divide index by 2 due to rc pwms
            id = self.pwms[pwm_index].id
            width = self.pwms[pwm_index].width
            results["hits"].extend([Hit(id, r.pos, r.pos + width, r.score, strand) for r in rs])

        return results


parser = argparse.ArgumentParser()

parser.add_argument("--pfms", required=True)
parser.add_argument("--fasta", required=True)
parser.add_argument("--pvalue", type=float)
parser.add_argument("--min-score", type=float)
parser.add_argument("--limit", default=0, type=int)
parser.add_argument("--cores", default=4, type=int)


if __name__ == "__main__":

    args = parser.parse_args()

    if args.pvalue and args.min_score:
        raise argparse.ArgumentError("cannot set pvalue and min_score at the same time")

    if args.pvalue:
        threshold = "pvalue"
        value = args.pvalue
    elif args.min_score:
        threshold = "min_score"
        value = args.min_score
    else:
        threshold = "pvalue"
        value = 0.01

    pfms = json.load(open(args.pfms, "r"))
    pwms = [PWM(p["PFM"], p["label"], threshold=(threshold, value)) for p in pfms]

    if args.limit != 0:
        seqs = list(islice(iter_fasta(args.fasta), 0, args.limit))
    else:
        seqs = iter_fasta(args.fasta)

    s = Scanner(pwms)

    def scan(fa: tuple[str, str]):
        return s.scan(fa)

    # could potentially cause memory issues with many sequences or low threshold
    with Pool(args.cores) as p:
        results = list(p.map(scan, seqs))

    print(json.dumps(results, indent=4))
