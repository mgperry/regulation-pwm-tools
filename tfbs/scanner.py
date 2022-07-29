from collections import namedtuple

import numpy as np

import MOODS.parsers
import MOODS.scan
import MOODS.tools
import pyfaidx

Hit = namedtuple("Hit", ["seqname", "start", "end", "name", "score", "strand"])

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

    def scan(self, seq: pyfaidx.Sequence):

        hits_by_tf = self.scanner.scan(seq.seq)

        results = []

        for i, hits in enumerate(hits_by_tf):
            strand = "+" if i % 2 else "-"  # odd pwms are reverse complement
            pwm_index = i // 2  # divide index by 2 due to rc pwms
            id = self.pwms[pwm_index].id
            width = self.pwms[pwm_index].width
            results.extend(
                [Hit(seq.name, seq.start + h.pos - 1, seq.start + h.pos + width - 1, id, h.score, strand) for h in hits]
            )

        return results

