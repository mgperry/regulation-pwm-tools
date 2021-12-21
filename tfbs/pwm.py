import numpy as np
import MOODS.tools


def as_col_vec(xs):
    return np.array([xs]).T


def calculate_PPM(pfm, background=(0.25, 0.25, 0.25, 0.25), p=1):
    pfm_adj = pfm + (as_col_vec(background) * p)
    return pfm_adj / pfm_adj.sum(axis=0)


def calculate_PWM(ppm, background=(0.25, 0.25, 0.25, 0.25)):
    return np.log2(ppm / as_col_vec(background))


def calculate_ICM(ppm, background=(0.25, 0.25, 0.25, 0.25)):
    return ppm * np.log2(ppm / as_col_vec(background))


class PWM:
    def __init__(
        self,
        PFM,
        id,
        background=(0.25, 0.25, 0.25, 0.25),
        pseudocount=0.8,
        pvalue=0.001,
    ):
        self.id = id
        self.background = background
        self.pseudocount = pseudocount

        self.PFM = np.array(PFM)
        self.PPM = calculate_PPM(self.PFM, background, pseudocount)
        self.PWM = calculate_PWM(self.PPM, background)
        self.ICM = calculate_ICM(self.PPM, background)
        self.width = np.size(self.PFM, 1)

        self.pvalue = pvalue
        self.threshold = MOODS.tools.threshold_from_p(self.PWM, background, pvalue)

        self.max_score = self.calculate_max_score()
        self.scaling_factor = self.calculate_lambda()

    def calculate_max_score(self):
        return np.sum(self.PWM.max(axis=0))

    def calculate_lambda(self):
        delta_S = self.max_score - self.threshold
        total_IC = np.sum(self.ICM)
        mismatch_bits = total_IC * 6 / 13.2
        return delta_S / mismatch_bits
