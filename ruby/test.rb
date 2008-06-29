#!/usr/bin/env ruby

require 'CentroidFold'
include CentroidFold

cf = CentroidFold::CentroidFold.new(CentroidFold::CentroidFold::CONTRAFOLD)
cf.calculate_posterior("GGGCCCAUAGCUCAGUGGUAGAGUGCCUCCUUUGCAAGGAGGAUGCCCUGGGUUCGAAUCCCAGUGGGUCCA")
[0.5, 1.0, 2.0, 4.0].each do |g|
  s,ea=cf.decode_structure(g)
  p [s,ea]
end

aln=[
     "-----GCUA-AUAUCGCUGUGGAAACACCUGGAACCAUCCCGAACCCAGC-AGUUAAGCACAGUGGAGCUAAAU--GUA--G--G-UAGUAAUACUG----AG-AAUA",
     "UCCGGUGACUUUACGCGUGAGGAAACACUCGUUCCCAUUCCGAACACGAC-AGUUAAGCUCCCG-CGGCCGAUGA--UAGUGCC--CA-CCA----GCGUGAA-AGUA"
    ]
cf.calculate_posterior(aln)
[0.5, 1.0, 2.0, 4.0].each do |g|
  s,ea=cf.decode_structure(g)
  p [s,ea]
end
