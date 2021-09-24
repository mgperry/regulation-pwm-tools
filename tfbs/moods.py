import argparse
import json
import numpy as np
import multiprocessing
import MOODS.parsers
import MOODS.scan
import MOODS.tools
from itertools import groupby, chain, islice

# from dataclasses import dataclass
from collections import namedtuple


class PWM:
    def __init__(self, PFM, id, background=(0.25, 0.25, 0.25, 0.25), pseudocount=0.8):
        self.PFM = np.array(PFM)
        self.id = id
        self.background = background
        self.pseudocount = pseudocount
        self._bg_col_vec = np.array([background]).T

        self.PWM = self.calculate_pwm()

    def get_pvalue_threshold(self, pvalue):
        return MOODS.tools.threshold_from_p_with_precision(
            self.PWM, self.background, pvalue, 2000.0, 4
        )

    def PPM(self):
        pfm_adj = self.PFM + (self._bg_col_vec * self.pseudocount)
        col_sums = np.sum(pfm_adj, axis=0)
        return pfm_adj / col_sums

    def PWM_rc(self):
        return tuple(map(tuple, np.flip(self.PWM, axis=[0, 1])))

    def calculate_pwm(self):
        return np.log(self.PPM()) - np.log(self._bg_col_vec)

    def get_percent_score(self, pc):
        range = self.PWM.max(axis=0) - self.PWM.min(axis=0)
        return np.sum(self.PWM.min(axis=0) + (pc * range))

    def PWM_to_tuples(self):
        return tuple(map(tuple, self.PWM))


# creates tuple of tuples, make this into a test of pfm_to_log_odds function
# pfm_file = "/home/malcolm/Data/Regulation/SELEX/matrices/ALX3_AE_TGCAAG20NGA_NTAATYNRATTAN_m1_c2_Cell2013.pfm"
# x = MOODS.parsers.pfm_to_log_odds(pfm_file, [0.25,0.25,0.25,0.25], 0.8)


def iter_fasta(filename):
    with open(filename, "r") as f:
        iterator = groupby(f, lambda line: line[0] == ">")
        for is_header, group in iterator:
            if is_header:
                header = next(group)[1:].strip()
            else:
                yield header, "".join(s.strip() for s in group)


# @dataclass
# class Hit:
#     '''Class for keeping track of PWM hits.'''
#     TF: str
#     pos: float
#     score: float
#     strand: str

Hit = namedtuple("Hit", "TF start score strand")


class Scanner:
    def __init__(
        self, pwms, background=(0.25, 0.25, 0.25, 0.25), threshold="pvalue", value=0.001
    ):
        # update TFs with PWMs and thresholds
        matrices = [p.PWM_to_tuples() for p in pwms]
        matrices_rc = [p.PWM_rc() for p in pwms]

        if threshold == "pvalue":
            thresholds = [p.get_pvalue_threshold(value) for p in pwms]
        elif threshold == "min_score":
            thresholds = [p.get_percent_score(value) for p in pwms]

        scanner = MOODS.scan.Scanner(7)  # why 7 lol
        scanner.set_motifs(matrices + matrices_rc, background, thresholds + thresholds)

        self.pwms = pwms
        self.scanner = scanner
        self.background = background
        self.threshold = {"method": threshold, "value": value}

    def scan(self, seq):
        raw = self.scanner.scan(seq)
        results = []
        n_pwms = len(self.pwms)
        for i, rs in enumerate(raw):
            if i >= n_pwms:
                i = i % len(self.pwms)
                strand = "-"
            else:
                strand = "+"
            tf_name = self.pwms[i].id
            results.extend([Hit(tf_name, r.pos, r.score, strand) for r in rs])
            # results.extend([{"TF": tf_name, "pos": r.pos, "score": r.score, "strand": "+"} for r in rs])

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

    pfms = json.load(open(args.pfms, "r"))
    pwms = [PWM(p["PFM"], p["label"]) for p in pfms]

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

    if args.limit != 0:
        seqs = list(islice(iter_fasta(args.fasta), 0, args.limit))
    else:
        seqs = iter_fasta(args.fasta)

    s = Scanner(pwms, threshold=threshold, value=value)

    def scan(fa_record):
        return {
            "header": fa_record[0],
            "seq": fa_record[1],
            "hits": s.scan(fa_record[1]),
        }

    p = multiprocessing.Pool(args.cores)

    print("ID,TF,start,score,strand")
    for result in p.imap(scan, seqs):
        id = result["header"].split(":")[-1]
        for hit in result["hits"]:
            print(id + "," + ",".join(str(x) for x in hit))

    # print(json.dumps(results, indent=4))

    # print("ID,TF,start,score,strand")

    # for r in results:
    #     id = r["header"].split(":")[-1]
    #     for hit in r["hits"]:
    #         print(id + ',' + ','.join(str(x) for x in hit))
