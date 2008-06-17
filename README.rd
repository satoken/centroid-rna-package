= CentroidFold for predicting RNA secondary structures

CentroidFold is one of the most accurate tools for predicting RNA
secondary structures. It is based on the gamma-centroid estimator
(Hamada et.al., 2008) for high-dimensional discrete spaces. Generally,
a gamma-centroid estimator is ((*slightly*)) more accurate than an MEA
estimator (Do et.al., 2005) under the same probability distribution. 

CentroidFold can take various probability distributions. Currently,
* McCaskill models implemented in Vienna RNA package,
* CONTRAfold models, and
* CG optimized models implemented in MultiRNAFold package
are supported. 

According to our benchmark, CentroidFold with CONTRAfold models will
predict the most accurate RNA secondary structures among currently
available prediction tools at this time.


== Requirements

* Boost C++ Library (>= 1.34.1)
  ((<URL:http://www.boost.org/>))
* Vienna RNA package (>= 1.6)
  ((<URL:http://www.tbi.univie.ac.at/~ivo/RNA/>))
* ruby (>= 1.8) (optional)
  ((<URL:http://www.ruby-lang.org>))
* CONTRAfold (>= 2.0) (optional)
  ((<URL:http://contra.stanford.edu/contrafold/>))
* MultiRNAFold package (>= 1.10) (optional)
  ((<URL:http://www.rnasoft.ca/>))


== Install

=== simple build

 ./configure && make && make install

=== build with CONTRAfold

Currently, a patch for CONTRAfold v2.0.0 is provided, which makes it
possible to link the routine for calculating base-pairing
probabilities in CONTRAfold with CentroidFold.

After applying this patch, you can build (({libcontrafold.a})).

 cd /somewhere/contrafold/src
 patch -p2 < /somewhere/centroid_fold-x.x.x/contrafold.patch
 make lib

Then, put (({libcontrafold.a})) into the CentroidFold build tree.

 cp libcontrafold.a /somewhere/centroid_fold-x.x.x/src/

Finally, build CentroidFold.

 cd /somewhere/centroid_fold-x.x.x/
 ./configure --with-contrafold && make && make install

If you build CentroidFold with (({--with-contrafold})),
CONTRAfold models are used as the default posterior probability
distribution rather than McCaskill models.

== Usage

Our tools accept the FASTA format as an input sequence.

=== CentroidFold with McCaskill models

(({centroid_fold})) is the main program of this package. It calculates
a base-pairing probability matrix for a given sequence with McCaskill
algorithm using (({pf_fold})) routine in Vienna RNA package, and then
estimates a gamma-centroid estimator for the base-pairing probability
matrix. 

 centroid_fold [options] seq

 Options:
   -g [ --gamma ] arg    weight of base pairs (default: 1.0 for Centroid, 6.0 for MEA)
   --mea                 run as an MEA estimator
   --aux                 use auxiliary base-pairing probabilities

If a negative value is given for the option '--gamma',
(({centroid_fold})) calculates secondary structures for several values
of gamma at the same time ({2^k | -5 <= k <= 10} and 6). 

Using the option '--aux', (({centroid_fold})) can take an auxiliary
base-pairing probability matrix instead of McCaskill models.

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


=== CentroidFold with CONTRAfold models

Even (({centroid_fold})) without (({libcontrafold.a})) can employ
CONTRAfold models using a wrapper script, (({contrafold.rb})),
which executes (({contrafold})) to calculate a base-pairing
probability matrix for a given sequence, and then executes
(({centroid_fold})) with '--aux' to estimate the most accurate
secondary structure from the base-pairing probability matrix.

 Usage: contrafold.rb [options] seq
    -g, --gamma gamma                weight of base-pairs
        --mea                        MEA estimators
        --viterbi                    viterbi estimators
    -t, --threshold th               threshold of posteriors
        --centroid_fold path         exec path of centroid_fold
        --contrafold path            exec path of contrafold
        --params params              specify a parameter file for contrafold

Example:
 % contrafold.rb -g -1 RF00008_B.fa
 > RF00008_B
 CAAAAGUCUGGGCUAAGCCCACUGAUGAGCCGCUGAAAUGCGGCGAAACUUUUG
 ...................................................... (g=0.03125,th=0.969697)
 .............................(((........)))........... (g=0.0625,th=0.941176)
 .........(((.....)))........((((........)))).......... (g=0.125,th=0.888889)
 .........(((.....)))........(((((......))))).......... (g=0.25,th=0.8)
 .(((((..((((.....)))).......(((((......)))))....))))). (g=0.5,th=0.666667)
 ((((((..(((((...))))).......(((((......)))))....)))))) (g=1,th=0.5)
 (((((((.(((((...))))).......(((((......)))))...))))))) (g=2,th=0.333333)
 (((((((.(((((...))))).......(((((......)))))...))))))) (g=4,th=0.2)
 (((((((((((((...)))))..)....(((((......)))))...))))))) (g=6,th=0.142857)
 (((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=8,th=0.111111)
 (((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=16,th=0.0588235)
 (((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=32,th=0.030303)
 (((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=64,th=0.0153846)
 (((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=128,th=0.00775194)
 (((((((((((((...)))))..).(..((((((...).)))))..)))))))) (g=256,th=0.00389105)
 (((((((((((((...)))))..).(..((((((...).)))))..)))))))) (g=512,th=0.00194932)
 (((((((((((((...)))))..).(..((((((...).)))))..)))))))) (g=1024,th=0.00097561)


=== CentroidFold with CG optimized models

(({simfold.rb})) is also a wrapper script which executes
(({simfold_pf})) instead of (({contrafold})).

 Usage: simfold.rb [options] seq
    -g, --gamma gamma                weight of base-pairs
        --mea                        MEA estimators
        --viterbi                    viterbi estimators
    -t, --threshold th               threshold of posteriors
        --centroid_fold path         exec path of centroid_fold
        --simfold path               exec path of simfold
        --simfold_pf path            exec path of simfold_pf
        --params params              specify a parameter file for simfold

Example:
 simfold.rb -g -1 --simfold_pf ~/MultiRNAFold-1.10/simfold_pf --params ~/MultiRNAFold-1.10/params/CG_best_parameters_ISMB2007.txt RF00008_B.fa
 > RF00008_B
 CAAAAGUCUGGGCUAAGCCCACUGAUGAGCCGCUGAAAUGCGGCGAAACUUUUG
 .........(((.....))).........((((......))))........... (g=0.03125,th=0.969697)
 .........(((.....)))........(((((......))))).......... (g=0.0625,th=0.941176)
 .........(((.....)))........(((((......))))).......... (g=0.125,th=0.888889)
 .........((((...))))........(((((......))))).......... (g=0.25,th=0.8)
 ....(...(((((...))))).......(((((......))))).....).... (g=0.5,th=0.666667)
 ((((((..(((((...))))).......(((((......)))))....)))))) (g=1,th=0.5)
 (((((((.(((((...))))).......(((((......)))))...))))))) (g=2,th=0.333333)
 (((((((.(((((...))))).......(((((......)))))...))))))) (g=4,th=0.2)
 (((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=6,th=0.142857)
 (((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=8,th=0.111111)
 (((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=16,th=0.0588235)
 (((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=32,th=0.030303)
 (((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=64,th=0.0153846)
 (((((((((((((...)))))..).(..(((((......)))))..)))))))) (g=128,th=0.00775194)
 (((((((((((((...)))))..).(..((((((...).)))))..)))))))) (g=256,th=0.00389105)
 (((((((((((((...)))))..).(..((((((...).)))))..)))))))) (g=512,th=0.00194932)
 (((((((((((((...)))))..).(..((((((...).)))))..)))))))) (g=1024,th=0.00097561)


== References

* Centroid estimators
  * Carvalho, L. E. and Lawrence, C. E.: Centroid estimation in
    discrete high-dimensional spaces with applications in
    biology. Proc Natl Sci USA, 105, 3209-3214, 2008.
  * Hamada, M., Kiryu, H., Sato, K., Mituyama, T. and Asai, K.:
    Estimation of secondary structure of an RNA sequence using
    generalized centroid estimator, submitted, 2008.
* CONTRAfold models and MEA estimators
  * Do, C. B., Woods, D. A. and Batzoglou, S.: CONTRAfold: RNA
    secondary structure prediction without physics-based
    models. Bioinformatics, 22, e90-e98, 2006.
* McCaskill models
  * McCaskill, J. S.: The equilibrium partition function and base pair
    binding probabilities for RNA secondary structure. Biopolymers,
    29, 1105-1119, 1990.
  * Hofacker, I. L.: Vienna RNA secondary structure server. Nucleic
    Acids Res, 31, 3429-3431, 2003
* CG optimized models
  * Andronescu, M., Condon, A., Hoos, H. H., Mathews, D. H. and
    Murphy, K. P.: Efficient parameter estimation for RNA secondary 
    structure prediction. Bioinformatics, 23, i19-i28, 2007
