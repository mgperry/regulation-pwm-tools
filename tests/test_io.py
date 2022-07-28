from tfbs.io import read_jaspar, read_moods
import numpy as np

def test_jaspar():
    pfms = list(read_jaspar("example_data/jaspar.txt"))
    assert len(pfms) == 2
    assert pfms[0]["name"] == "CTCF"
    assert type(pfms[1]["PFM"]) == np.ndarray
    assert pfms[1]["PFM"].shape == (4, 14)

def test_moods():
    pfms = list(read_moods("example_data/moods/"))
    assert len(pfms) == 2
    assert pfms[0]["name"] == "CTCF"
    assert type(pfms[1]["PFM"]) == np.ndarray
    assert pfms[1]["PFM"].shape == (4, 19)
