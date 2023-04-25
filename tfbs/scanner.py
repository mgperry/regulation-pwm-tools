from dataclasses import dataclass, replace

import numpy as np

import MOODS.parsers
import MOODS.scan
import MOODS.tools
import pyfaidx

from .motif import Motif

@dataclass(frozen=True)
class TFBS:
    """
    Class representing a TFBS as an interval, same fields as bed/GRanges.

    NB this uses "GRanges" style inclusive ranges, not BED-style insanity.
    """

    __slots__ = ["chr", "start", "end", "name", "score", "strand"]

    chr: str
    start: int
    end: int
    name: str
    score: float
    strand: str

    def to_bed(self) -> str:
        bed = replace(self, start = self.start - 1)
        return "%s\t%d\t%d\t%s\t%.2f\t%s" % bed


class Scanner:
    def __init__(self, motifs: list[Motif], background: tuple[int] =(0.25, 0.25, 0.25, 0.25)):

        matrices = []
        thresholds = []

        for motif in motifs:
            matrices.extend([motif.matrix, motif.complement])
            thresholds.extend([motif.threshold] * 2)

        scanner = MOODS.scan.Scanner(7)  # why 7 lol
        scanner.set_motifs(matrices, background, thresholds)

        self.motifs = motifs
        self.scanner = scanner
        self.background = background

    def scan(self, seq: pyfaidx.Sequence, simplify=True) -> list[TFBS]:

        hits_by_tf = self.scanner.scan(seq.seq)

        results = []

        for i, hits in enumerate(hits_by_tf):
            # check for empty
            if not hits:
                continue

            strand = "+" if i % 2 else "-"  # odd pwms are reverse complement
            motif_ix = i // 2  # divide index by 2 due to rc pwms
            id = self.motifs[motif_ix].id
            width = self.motifs[motif_ix].width
            start = seq.start - 1
            end = seq.start + width - 1
            rs = [TFBS(seq.name, start + h.pos, end + h.pos, id, h.score, strand) for h in hits]
            if simplify:
                results.extend(rs)
            else:
                results.append(rs)


        return results

