from dataclasses import dataclass, replace

import MOODS.parsers
import MOODS.scan
import MOODS.tools

from pyfaidx import Sequence

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

    def scan(self, seq: Sequence) -> list[TFBS]:

        results = []

        for i, tf_matches in enumerate(self.scanner.scan(seq.seq)):
            # check for empty
            if not tf_matches:
                continue

            strand = "+" if i % 2 else "-"  # odd pwms are reverse complement
            motif_ix = i // 2  # divide index by 2 due to rc pwms

            name = self.motifs[motif_ix].id
            width = self.motifs[motif_ix].width

            start = seq.start - 1
            end = start + width

            rs = [TFBS(seq.name, start + m.pos, end + m.pos, name, m.score, strand) for m in tf_matches]

            results.extend(rs)

        return results

