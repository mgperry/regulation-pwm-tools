from collections import namedtuple

import numpy as np

import MOODS.parsers
import MOODS.scan
import MOODS.tools
import pyfaidx


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

    def scan(self, seq: pyfaidx.FastaRecord):
        header = seq.name
        seq = seq.seq

        raw = self.scanner.scan(seq)

        results = []

        for i, rs in enumerate(raw):
            strand = "+" if i % 2 else "-"  # odd pwms are reverse complement
            pwm_index = i // 2  # divide index by 2 due to rc pwms
            id = self.pwms[pwm_index].id
            width = self.pwms[pwm_index].width
            results.extend(
                [Hit(header, id, r.pos, r.pos + width, r.score, strand) for r in rs]
            )

        return results

