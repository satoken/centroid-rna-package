#!/usr/bin/env python

import CentroidFold

cf = CentroidFold.CentroidFold(CentroidFold.CentroidFold.CONTRAFOLD)
cf.calculate_posterior("GGGCCCAUAGCUCAGUGGUAGAGUGCCUCCUUUGCAAGGAGGAUGCCCUGGGUUCGAAUCCCAGUGGGUCCA")
s,ea=cf.decode_structure(4)
print s,ea

aln = [
    "-----GCUA-AUAUCGCUGUGGAAACACCUGGAACCAUCCCGAACCCAGC-AGUUAAGCACAGUGGAGCUAAAU--GUA--G--G-UAGUAAUACUG----AG-AAUA",
    "UCCGGUGACUUUACGCGUGAGGAAACACUCGUUCCCAUUCCGAACACGAC-AGUUAAGCUCCCG-CGGCCGAUGA--UAGUGCC--CA-CCA----GCGUGAA-AGUA"
    ]
cf.calculate_posterior(aln)
s,ea=cf.decode_structure(4)
print s,ea
