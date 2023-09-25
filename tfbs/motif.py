import numpy as np
import MOODS.tools


class Motif:
    def __init__(
        self,
        id: str,
        matrix: np.ndarray,
        threshold: float = None,
    ):

        self.id = id
        self.matrix = np.array(matrix)
        self.threshold = threshold

        self.width = np.size(self.matrix, 1)
        self.max_score = np.sum(self.matrix.max(axis=0))
        self.consensus = "".join("ACGT"[i] for i in self.matrix.argmax(0))

    @property
    def complement(self):
        return np.flip(self.matrix, axis=[0, 1])
        
    @classmethod
    def with_p(cls, id, matrix, p, bg=None):
        motif = cls(id, matrix)
        motif.threshold = motif.threshold_from_p(p, bg)
        return motif

    def threshold_from_p(self, p: float, bg: tuple[float] = None):
        if bg is None:
            bg = (0.25,) * 4
        return MOODS.tools.threshold_from_p(self.matrix, bg, p)

    @classmethod
    def from_pfm(cls, *args, **kwargs):
        return make_pwm_motif_from_pfm(*args, **kwargs)
        

def as_col_vec(xs):
    return np.array([xs]).T


def calculate_PPM(pfm, bg=(0.25, 0.25, 0.25, 0.25), p=1):
    pfm_adj = pfm + (as_col_vec(bg) * p)
    return pfm_adj / pfm_adj.sum(axis=0)


def calculate_PWM(ppm, bg=(0.25, 0.25, 0.25, 0.25)):
    return np.log2(ppm / as_col_vec(bg))


# def calculate_ICM(ppm, background=(0.25, 0.25, 0.25, 0.25)):
#     return ppm * np.log2(ppm / as_col_vec(background))


def make_pwm_motif_from_pfm(
    id,
    PFM,
    pvalue=None,
    background=(0.25, 0.25, 0.25, 0.25),
    pseudocount=0.8
):

    PFM = np.array(PFM)
    PPM = calculate_PPM(PFM, background, pseudocount)
    PWM = calculate_PWM(PPM, background)
    
    matrix = Motif(id, PWM)
    if pvalue is not None:
        matrix.threshold = matrix.threshold_from_p(background, pvalue)

    matrix.background = background
    matrix.pvalue = pvalue

    return matrix

