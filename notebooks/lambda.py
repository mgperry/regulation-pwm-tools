#%%
import numpy as np
from tfbs.pwm import PWM, as_col_vec
from MOODS.tools import threshold_from_p

test_pfm = {
    "id": "test_pfm",
    "PFM": [
        [6, 4, 0, 5, 5, 4],
        [0, 0, 2, 0, 0, 0],
        [0, 0, 3, 0, 0, 0],
        [0, 2, 1, 1, 1, 2],
    ],
}

# %%
params = {
    # "mismatch_tolerance": 6,
    # "max_mismatch": 13.2,  # bits
    "top_sites": 1e-3,  # top 0.1% of sites
}

pwm = PWM(test_pfm["PFM"], test_pfm["id"], pvalue=params["top_sites"])

pwm.calculate_lambda()

print(f"scaling factor for pwm {pwm.id} is {pwm.scaling_factor:.3f}")

# %%
