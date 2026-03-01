#!/usr/bin/env python3
# Takes in sequences in tab file format, reverse complements, and prints:
# sequence name \n sequence

import sys

class Gene:
    # Supports uppercase + lowercase DNA bases
    _COMPLEMENT_TABLE = str.maketrans("ATGCatgc", "TACGtacg")

    def __init__(self, creationid, creationseq):
        self.id = creationid
        self.seq = creationseq

    def reverse_complement(self):
        return self.seq[::-1].translate(self._COMPLEMENT_TABLE)

for line in sys.stdin:
    seq_list = line.strip().split("\t")
    if len(seq_list) < 2:
        continue  # skip malformed/blank lines

    seq_id = seq_list[0]
    seq_nu = seq_list[1]
    new_seq = Gene(seq_id, seq_nu)

    print(seq_id)
    print(new_seq.reverse_complement())