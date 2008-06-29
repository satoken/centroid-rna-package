#!/usr/bin/perl
use CentroidFold;

$cf = new CentroidFold::CentroidFold($CentroidFold::CentroidFold::CONTRAFOLD);
$cf->calculate_posterior("GGGCCCAUAGCUCAGUGGUAGAGUGCCUCCUUUGCAAGGAGGAUGCCCUGGGUUCGAAUCCCAGUGGGUCCA");
($a, $s)=$cf->decode_structure(4);
print $a, " ", $s, "\n";

$aln = [
    "-----GCUA-AUAUCGCUGUGGAAACACCUGGAACCAUCCCGAACCCAGC-AGUUAAGCACAGUGGAGCUAAAU--GUA--G--G-UAGUAAUACUG----AG-AAUA",
    "UCCGGUGACUUUACGCGUGAGGAAACACUCGUUCCCAUUCCGAACACGAC-AGUUAAGCUCCCG-CGGCCGAUGA--UAGUGCC--CA-CCA----GCGUGAA-AGUA"
    ];

$cf->calculate_posterior($aln);
($a,$s)=$cf->decode_structure(4);
print $a, " ", $s, "\n";
