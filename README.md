# MOODS-tools

This package provides a more pythonic wrapper around the MOODS Cpp // swig motif
scanning tool.

There is also a notebook and some test data showing an example workflow. The deps
for this are included as poetry development dependencies.

## Package Structure

There are two main classes in the package, `Scanner` and `Motif`. `Motif` provides a
representation of a scanning matrix (for example a PWM, but not limited to just this) in
`numpy` format, alongside an ID field and a threshold. There are convenience functions
for loading PWMs from PFMs, and for calculating p-values (this just wraps MOODS).

The `Scanner` class takes a list of motifs, and builds a scanner, which will efficiently
scan a nucleotide sequence in one pass. This is a bit fiddly, since the scanner is not
strand-aware, and returns hits by motif order. I've taken a similar approach to the
`moods-py` script we use in regulation, but this is wrapped up so that the `.scan()`
method takes a `pyfaidx.Sequence` and outputs a list of `TFBS`s (small dataclass, effectively
just a `BED` row, tuples might be faster).

There is also a cli `scan` command, which should be very similar to the `moods-py` script,
but taking PWMs in `json` input, and some functions in `readers.py` for parsing MOODS
or JASPAR PWMs.

## Multicore

Currently the package doesn't have an easy way to use multiple cores: because `Scanner`
uses Swig objects, they can't be pickled. The most straightforward way would probably
be to create `N` processes to do the scanning, and just build a new `Scanner` in each
one.
