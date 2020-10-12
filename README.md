CentroidFold for predicting RNA secondary structures
====================================================

CentroidFold is one of the most accurate tools for predicting RNA
secondary structures. It is based on the gamma-centroid estimator
(Hamada et.al., 2009) for high-dimensional discrete spaces. Generally,
a gamma-centroid estimator is *slightly* more accurate than an MEA
estimator (Do et.al., 2005) under the same probability distribution. 

CentroidAlifold can predict common secondary structures for multiple
alignments of RNA sequences by using an averaged gamma-centroid
estimator (Hamada et al., 2010).

CentroidFold can employ various probability distributions. Currently,
* the CONTRAfold model,
* the McCaskill model implemented in the Vienna RNA package,
* the RNAalifold model implemented in the Vienna RNA package, and
* the pfold model
are supported. 

According to our benchmark, CentroidFold with the McCaskill model
using Boltzmann likelihood parameters (Andronescue et al., 2010)
will predict the most accurate RNA secondary structures among
currently available prediction tools at this time.


Requirements
------------

* [Boost C++ Library](http://www.boost.org/) (>= 1.40.0)
* [Vienna RNA package](http://www.tbi.univie.ac.at/~ivo/RNA/) (>= 1.8)
* [pfold](http://www.daimi.au.dk/~compbio/pfold/) (optional)


Install
-------

### Using configure

	./configure && make && make install

If the Vienna RNA package has been installed at the non-standard
directory, specify the location like:

	./configure --with-vienna-rna=/path/to/vienna-rna

### Using cmake

	mkdir build && cd build
	cmake .. && make


Usage
-----

`centroid_fold` can take the FASTA format as an input sequence,
then predict its secondary structure.
`centroid_alifold` can take the CLUSTAL format as an input
aligment, the predict its common secondary structure.

### Secondary structure prediction for single sequences

`centroid_fold` is the main program of this package. It calculates
a base-pairing probability matrix for a given sequence with one of
several algorithms including the CONTRAfold model and the McCaskill
algorithm using `pf_fold` routine in Vienna RNA package. Then,
`centroid_fold` estimates a gamma-centroid estimator for the
base-pairing probability matrix. 

	centroid_fold [options] seq

	Options:
	  -e [ --engine ] arg     specify the inference engine (default: "McCaskill")
	  -g [ --gamma ] arg      weight of base pairs
	  --noncanonical          allow non-canonical base-pairs
	  -C [ --constraints ]    use structure constraints
	  --postscript arg        draw predicted secondary structures with the 
                              postscript (PS) format

By `-e` or `--engine` option, you can select the inference engine for
calculating base-pairing probabilities from the CONTRAfold model
("CONTRAfold"), the McCaskill model ("McCaskill" if available) and
the pfold model ("pfold" if available). If you use the pfold model,
several environment variables should be set: `PFOLD_BIN_DIR` to the
directory which contains pfold binaries, and `AWK_BIN` and `SED_BIN`
to awk and sed program available on your system, respectively.

If a negative value is given for the option `--gamma`,
`centroid_fold` calculates secondary structures for several values
of gamma at the same time ({2^k | -5 <= k <= 10} and 6). 

For long sequences, you can use `-d` options to restrict the maximum
distance of base-pairs. This can reduce the computation time and the
memory requirement: O(L^3) to O(LW^2), and O(L^2) to O(LW),
respectively, where L is length of a sequence and W is the specified
maximum distance of base-pairs. This option is available only for the
CONTRAfold model.

If a part of secondary structure for a given sequence is known, you
can specify it by a modified FASTA format like:

	> RF00008_B
	CAAAAGUCUGGGCUAAGCCCACUGAUGAGCCGCUGAAAUGCGGCGAAACUUUUG
	.(((((........???????????????????????...........))))).

and run `centroid_fold` with `-C` options. The positions at '('
and ')' are restricted to base-pairs, the position at '.' is
restricted to unpaired bases, and the position at '?' is
unrestricted.

If `--sampling` option is given, `centroid_fold` uses the
stochastic traceback algorithm instead of the McCaskill's base-pairing
probability matrix. Like Sfold (Ding et al., 2005), build clusters of
secondary structures, and then compute their centroids. The number of
clusters can be specified by `--max-clusters` option.

Example:

	% centroid_fold -g -1 RF00008_B.fa
	> RF00008_B
	CAAAAGUCUGGGCUAAGCCCACUGAUGAGCCGCUGAAAUGCGGCGAAACUUUUG
	.(((((...(((.....)))........(((((......)))))....))))). (g=0.03125,th=0.969697)
	.(((((...(((.....)))........(((((......)))))....))))). (g=0.0625,th=0.941176)
	((((((...(((.....)))........(((((......)))))....)))))) (g=0.125,th=0.888889)
	((((((...((((...))))........(((((......)))))....)))))) (g=0.25,th=0.8)
	(((((((..((((...))))........(((((......)))))...))))))) (g=0.5,th=0.666667)
	(((((((.(((((...))))).......(((((......)))))...))))))) (g=1,th=0.5)
	(((((((.(((((...))))).......(((((......)))))...))))))) (g=2,th=0.333333)
	(((((((.(((((...))))).......(((((......)))))...))))))) (g=4,th=0.2)
	(((((((.(((((...))))).......(((((......)))))...))))))) (g=6,th=0.142857)
	(((((((.(((((...))))).......(((((......)))))...))))))) (g=8,th=0.111111)
	(((((((((((((...))))).......(((((......))))))..))))))) (g=16,th=0.0588235)
	(((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=32,th=0.030303)
	(((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=64,th=0.0153846)
	(((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=128,th=0.00775194)
	(((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=256,th=0.00389105)
	(((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=512,th=0.00194932)
	(((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=1024,th=0.00097561)


### Common secondary structure prediction for multiple alignments

For the CLUSTAL format, `centroid_alifold` predicts common
secondary structures for the given multiple alignments.

	centroid_alifold [options] seq

	Options:
	  -h [ --help ]           show this message
	  -e [ --engine ] arg     specify the inference engine (default: "McCaskill & 
                              Alifold")
	  -w [ --mixture ] arg    mixture weights of inference engines
      -g [ --gamma ] arg      weight of base pairs
	  --noncanonical          allow non-canonical base-pairs
	  -C [ --constraints ]    use structure constraints
	  --postscript arg        draw predicted secondary structures with the 
                              postscript (PS) format

By `-e` or `--engine` option, you can select the inference engine for
calculating base-pairing probabilities from the CONTRAfold model
("CONTRAfold"), the McCaskill model ("McCaskill" if available),
the RNAalifold model ("Alifold" if available) and the pfold model
("pfold" if avilable).

If you specify the inference engines multiply, `centroid_alifold`
employs a mixtured baes-pairing probability matrix. The mixture
weight can be set by `-w` or `--mixture` option. The default
setting of `centroid_alifold` is `-e McCaskill -w 1.0 -e Alifold
-w 1.0`. See more detail in (Hamada et al., 2010).

Example:

	% centroid_alifold -g -1 RF00436.aln
	>AB029447-1/1210-1265
	--BCAHuUGYAVgUCGCUUUGGAYAaaAG--CGUCUGCUAAAUGM-VURwrukKAAAUDu-
	............................................................. (g=0.03125,th=0.969697)
	............................................................. (g=0.0625,th=0.941176)
	............................................................. (g=0.125,th=0.888889)
	............................................................. (g=0.25,th=0.8)
	...............(((.........))..)............................. (g=0.5,th=0.666667)
	..............(((((.......)))..))............................ (g=1,th=0.5)
	............(.((((((.....))))..)).).......................... (g=2,th=0.333333)
	........((..(.((((((.....))))..)).).))....................... (g=4,th=0.2)
	...((((.((.((.((((((.....))))..)).))))..))))................. (g=6,th=0.142857)
	...(((((((.((.((((((.....))))..)).)))).)))))................. (g=8,th=0.111111)
	..((((((((.((.((((((.....))))..)).)))).))))))................ (g=16,th=0.0588235)
	..((((((((.((.((((((.....))))..)).)))).))))))((.....))....... (g=32,th=0.030303)
	..((((((((.((.((((((.....))))..)).)))).))))))((.....))....... (g=64,th=0.0153846)
	..((((((((.((.((((((.....))))..)).)))).))))))((.....))....... (g=128,th=0.00775194)
	..((((((((.((.((((((.....))))..)).)))).))))))(((...)))....... (g=256,th=0.00389105)
	.(((((((((.((.((((((.....))))..)).)))).)))))).).((((....)))). (g=512,th=0.00194932)
	.(((((((((.((.((((((.....))))..)).)))).)))))).)(((((....))))) (g=1024,th=0.00097561)

The first line of the result is the description of the first sequence
in the given alignment. The second is the "most informative sequence"
(Freyhult et al., 2005), which is similar to IUPAC ambiguity
characters, produced by a library routine of the Vienna Package.


### Secondary structure prediction using homologous sequence information 

If homologous sequences to the target sequence are available,
`centroid_homfold` can predict secondary structures for the target
sequence more accurately than `centroid_fold` using homologous
sequence information with the probabilistic consistency transformation
for base-pairing probabilities (Hamada et al, 2009).

	centroid_homfold [options] seq

	Options:
	  -H [ --homologous ] arg fasta file containing homologous sequences (REQUIRED)
	  --engine_s arg          specify the inference engine for secondary structures
                              (default: "McCaskill")
	  --engine_a arg          specify the inference engine for pairwise alignments 
                              (default: "CONTRAlign")
	  -g [ --gamma ] arg      weight of base pairs
	  --postscript arg        draw predicted secondary structures with the 
                              postscript (PS) format

You should give homologous sequences to the target sequence in the
FASTA format by `-H` option.

Example:

	% centroid_homfold -g -1 -H RF00005.fa seq.fa 
	>X12857.1/421-494
	GCGGAUGUAGCCAAGUGGAUCAAGGCAGUGGAUUGUGAAUCCACCAUGCGCGGGUUCAAUUCCCGUCAUUCGCC
	.......................................................................... (g=0.03125,th=0.969697,e=0)
	.......................................................................... (g=0.0625,th=0.941176,e=0)
	.......................................................................... (g=0.125,th=0.888889,e=0)
	.......................................................................... (g=0.25,th=0.8,e=0)
	.......................................................................... (g=0.5,th=0.666667,e=0)
	.......................................................................... (g=1,th=0.5,e=0)
	.(((((.............................................(((.......)))...))))).. (g=2,th=0.333333,e=-1.21)
	(((((((...........................................((((.......)))).))))))). (g=4,th=0.2,e=-10.8)
	(((((((.....................((((.......))))......(((((.......)))))))))))). (g=6,th=0.142857,e=-17.9)
	(((((((...(.............)..(((((.......))))).....(((((.......)))))))))))). (g=8,th=0.111111,e=-17.22)
	(((((((..(((...........))).(((((.......))))).....(((((.......)))))))))))). (g=16,th=0.0588235,e=-23.8)
	(((((((..((((.........))))(((((((.....))))))..)..((((((....).)))))))))))). (g=32,th=0.030303,e=-13.3)
	(((((((..((((.(....)..))))((((((((...)))))))..)..((((((....).)))))))))))). (g=64,th=0.0153846,e=-8.9)
	(((((((..((((((....)).))))((((((((...)))))))..)..((((((....).)))))))))))). (g=128,th=0.00775194,e=-10.6)
	(((((((..((((((....)).))))((((((((...)))))))..)..((((((....).)))))))))))). (g=256,th=0.00389105,e=-10.6)
	(((((((.(((((((....)).))))((((((((...)))))))..).)((((((....).)))))))))))). (g=512,th=0.00194932,e=-4.9)
	(((((((.(((((((....)).))))((((((((...)))))))..).)((((((....).)))))))))))). (g=1024,th=0.00097561,e=-4.9)


References
----------

* Centroid estimators
  * Carvalho, L.E. and Lawrence, C.E.: Centroid estimation in
    discrete high-dimensional spaces with applications in
    biology. *Proc Natl Sci USA*, 105:3209-3214, 2008.
  * Ding, Y., Chan, C. Y., and Lawrence, C.E.: RNA secondary
    structure prediction by centroids in a Boltzmann weighted
    ensemble, *RNA*, 11:1157-1166, 2005
  * Hamada, M., Kiryu, H., Sato, K., Mituyama, T. and Asai, K.:
    Predictions of RNA secondary structure using generalized centroid
    estimators, *Bioinformatics*, 25:465-473, 2009
  * Hamada, M., Sato, K., Kiryu, H., Mituyama, T. and Asai, K.:
    Predictions of RNA secondary structures by combining homologous
    sequence information, *Bioinformatics*, 23:i330-338, 2009.
  * Hamada, M., Sato, K., and Asai, K.: Improving the accuracy of
    predicting secondary structure for aligned RNA sequences, 
	*Nucleic Acids Res*, 39:393–402, 2011.
* The CONTRAfold model and MEA estimators
  * Do, C.B., Woods, D.A. and Batzoglou, S.: CONTRAfold: RNA
    secondary structure prediction without physics-based
    models. *Bioinformatics*, 22:e90-e98, 2006.
* The McCaskill model
  * McCaskill, J.S.: The equilibrium partition function and base pair
    binding probabilities for RNA secondary structure. *Biopolymers*,
    29, 1105-1119, 1990.
  * Hofacker, I.L.: Vienna RNA secondary structure server. 
    *Nucleic Acids Res*, 31:3429-3431, 2003.
* The RNAalifold model
  * Bernahart, S., Hofacker, I.L., Will, S., Gruber, A.R., and
    Stadler, P.F.: RNAalifold: improved consensus structure prediction
    for RNA alignments, *BMC Bioinformatics*, 9:474, 2008.
* The pfold model
  * Knudsen, B., and Hein, J.: Using stochastic context free grammars
    and molecular evolution to predict RNA secondary
    structure. *Bioinformatics*, 15:446-454, 1999.
  * Knudsen, B., and Hein, J.: Pfold: RNA secondary structure
    prediction using stochastic context-free grammars. 
	*Nucleic Acids Res*, 31:3423-3428, 2003.
* Boltzmann likelihood parameters
  * Andronescu, M., Condon, A., Hoos, H.H., Mathews, D.H., and Murphy,
    K.P.: Computational approaches for RNA energy parameter
    estimation. *RNA*, 16:2304-18, 2010
* Others
  * Freyhult, E., Moulton, V., and Gardner, PP.: Predicting RNA
    structure using mutual information. *Appl Bioinformatics*. 4:53-59,
    2004.
