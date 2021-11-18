import numpy as np
import MOODS.tools


def as_col_vec(xs):
    return np.array([xs]).T


def calculate_PPM(pfm, background=(0.25, 0.25, 0.25, 0.25), p=1):
    pfm_adj = pfm + (as_col_vec(background) * p)
    return pfm_adj / pfm_adj.sum(axis=0)


def calculate_PWM(ppm, background=(0.25, 0.25, 0.25, 0.25)):
    return np.log2(ppm / as_col_vec(background))

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
        self.pvalue = pvalue

        self.PFM = np.array(PFM)
        self.PPM = calculate_PPM(self.PFM, background, pseudocount)
        self.PWM = calculate_PWM(self.PPM, background)

        self.threshold = MOODS.tools.threshold_from_p(self.PWM, background, pvalue)
        self.width = np.size(self.PFM, 1)

    def pwm_to_tuples(self):
        return tuple(map(tuple, self.PWM))

    def rev_pwm_to_tuples(self):
        return tuple(map(tuple, np.flip(self.PWM, axis=[0, 1])))
