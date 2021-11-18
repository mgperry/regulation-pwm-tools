#%%
import numpy as np
from tfbs.pwm import PWM, calculate_PPM, calculate_PWM, as_col_vec
from MOODS.tools import threshold_from_p

myc = {
    "id": "MYCN_HUMAN.H11MO.0.A",
    "PFM": [
        [146, 132, 80, 12, 426, 3, 40, 44, 1, 32, 109, 97],
        [85, 24, 282, 485, 3, 448, 6, 38, 0, 142, 60, 204],
        [227, 314, 78, 0, 41, 2, 448, 1, 486, 267, 203, 163],
        [41, 29, 59, 2, 29, 46, 5, 416, 12, 58, 127, 35],
    ],
}

# %%
params = {
    "mismatch_tolerance": 6,
    "max_mismatch": 13.2,  # bits
    "top_sites": 1e-3,  # top 0.1% of sites
}

# %%
def calculate_ICM(ppm, background=[0.25, 0.25, 0.25, 0.25]):
    return ppm * np.log2(ppm / as_col_vec(background))

max_score = 0

def calculate_lambda(pwm: PWM):
    max_score = np.sum(pwm.PWM.max(axis=0))
    delta_S = max_score - threshold_from_p(pwm.PWM, pwm.background, params["top_sites"])
    mismatch_bits = np.sum(calculate_ICM(pwm.PPM), background=pwm.background) * 6/13.2
    return delta_S / mismatch_bits


# %%
