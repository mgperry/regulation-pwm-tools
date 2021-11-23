from itertools import groupby


def iter_fasta(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    fa_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for group in fa_iter:
        # drop the ">"
        header = next(group)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in next(fa_iter))
        yield header, seq
