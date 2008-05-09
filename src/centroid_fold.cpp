// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <iostream>
#include <string>
#include <cassert>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/range.hpp>
#include <boost/algorithm/string.hpp>
#include <deque>
//#include "rule.h"
#include "iupac.h"
#include "mea.h"
#include "centroid.h"
#include "fa.h"
#include "aln.h"

#ifdef HAVE_LIBRNA
namespace Vienna {
extern "C" {
#include <ViennaRNA/PS_dot.h>
#include <ViennaRNA/aln_util.h>
};
};
#endif

namespace po = boost::program_options;

typedef SCFG::BP::Table<double> BPTable;
typedef boost::shared_ptr<BPTable> BPTablePtr;

#ifdef HAVE_LIBRNA
bool ps_out;
bool svg_out;
#endif

template < class BPTable >
void
output(const std::string& name, const std::string& seq, const BPTable& bp,
       bool centroid, const std::vector<double>& gamma)
{
  std::cout << ">" << name << std::endl
	    << seq << std::endl;
  std::string paren(seq.size(), '.');
  if (centroid) {
    std::vector<double>::const_iterator g;
    for (g=gamma.begin(); g!=gamma.end(); ++g) {
      std::fill(paren.begin(), paren.end(), '.');
      //double p =
      SCFG::Centroid::execute(bp, paren, 1.0, *g);
      //std::cout << paren << " (g=" << *g << ",EA=" << p << ")" << std::endl;
      std::cout << paren << " (g=" << *g << ",th=" << (1.0/(1.0+*g)) << ")" << std::endl;
    }
  } else {
    std::vector<double>::const_iterator g;
    for (g=gamma.begin(); g!=gamma.end(); ++g) {
      std::fill(paren.begin(), paren.end(), '.');
      double p = SCFG::MEA::execute(bp, *g, paren);
      std::cout << paren << " (g=" << *g << ",EA=" << p << ")" << std::endl;
    }
  }
#ifdef HAVE_LIBRNA
  {
    std::string fbase(name);
    boost::algorithm::trim(fbase);
    std::replace(boost::begin(fbase), boost::end(fbase), ' ', '_');
    std::replace(boost::begin(fbase), boost::end(fbase), '/', '_');
    if (ps_out) {
      char fname[100];
      sscanf(fbase.c_str(), "%12s", fname);
      //Vienna::rna_plot_type=0;
      strcat(fname, "_ss.ps");
      Vienna::PS_rna_plot(const_cast<char*>(seq.c_str()),
			  const_cast<char*>(paren.c_str()),
			  fname);
    }
    if (svg_out) {
      char fname[100];
      sscanf(name.c_str(), "%12s", fname);
      //Vienna::rna_plot_type=0;
      strcat(fname, "_ss.svg");
      Vienna::svg_rna_plot(const_cast<char*>(seq.c_str()),
			   const_cast<char*>(paren.c_str()), fname);
    }
  }
#endif
}

int
main(int argc, char* argv[])
{
  std::vector<double> gamma;
  std::string input;
  std::vector<std::string> model;
  bool centroid=false;
  float p_th=0.0;
  
  // parse command line options
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "show this message")
    ("gamma,g",
     po::value<std::vector<double> >(&gamma),
     "weight of base pairs (default: 1.0 for Centroid, 6.0 for MEA)")
    ("mea", "run as an MEA estimator")
#ifdef HAVE_LIBRNA
    ("alipf_fold", "use alipf_fold base-pairing probabilities")
#endif
    ("aux", "use auxiliary base-pairing probabilities")
    ("posteriors", po::value<float>(&p_th),
     "output base-pairing probability matrices which contain base-pairing probabilities more than the given value.")
#ifdef HAVE_LIBRNA
    ("ps", "draw secondary structures into a PS file")
    ("svg", "draw secondary structures into a SVG file")
#endif
    ;
  po::options_description opts("Options");
  opts.add_options()
    ("seq-file", po::value<std::string>(&input), "training sequence filename")
    ("model-file", po::value<std::vector<std::string> >(&model), "model filename");

  opts.add(desc);
  po::positional_options_description pd;
  pd.add("seq-file", 1).add("model-file", -1);
  po::parsed_options parsed =
    po::command_line_parser(argc, argv).options(opts).positional(pd).run();
  po::variables_map vm;
  po::store(parsed, vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("seq-file") ||
      (vm.count("aux") && model.empty())) {
    std::cout << "Centroid fold for structural RNAs" << std::endl
	      << "Usage:" << std::endl
	      << " " << argv[0]
	      << " [options] seq [model ...]\n\n"
	      << desc << std::endl;
    return 1;
  }

  centroid = !vm.count("mea");
  if (gamma.size()==1 && gamma[0]==-1.0) {
    double g[] = { 0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0,
		   8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0 };
    gamma.resize(boost::size(g));
    std::copy(boost::begin(g), boost::end(g), gamma.begin());
  }
  if (gamma.empty()) gamma.push_back(centroid ? 1.0 : 6.0);
#ifdef HAVE_LIBRNA
  ps_out = vm.count("ps");
  svg_out = vm.count("svg");
#endif

  boost::spirit::file_iterator<> fi(input.c_str());
  if (!fi) {
    perror(input.c_str());
    return 1;
  }
  while (1) {
    Fasta fa;
    Aln aln;
    BPTable bp;
    if (fa.load(fi)) {		// for single sequences
      if (vm.count("aux")) {
	bp.load(model[0].c_str());
      } else {
	bp.pf_fold(fa.seq());
      }

      if (!vm.count("posteriors")) {
	output(fa.name(), fa.seq(), bp, centroid, gamma);
      } else {
	bp.save(std::cout, fa.seq(), p_th);
      }

    } else if (aln.load(fi)) {	// for multiple alignments
      if (vm.count("alipf_fold")) {
	BPTable bp;
	bp.alipf_fold(aln.seq());
	if (!vm.count("posteriors")) {
	  output(aln.name().front(), aln.seq().front(), bp, centroid, gamma);
	} else {
	  bp.save(std::cout, aln.seq().front(), p_th);
	}

      } else {
	if (vm.count("aux")) {
	  if (aln.seq().size()!=model.size()) {
	    std::cerr << "Given BP matrices and alignments were inconsistent.\n";
	    break;
	  }
	}
	std::list<BPTablePtr> bps;
	std::list<std::string> seqs;
	std::list<std::vector<uint> > idxmaps;
	BPTable::convert_to_raw_sequences(aln.seq(), seqs, idxmaps);
	uint i=0;
	std::list<std::string>::const_iterator x;
	for (x=seqs.begin(), i=0; x!=seqs.end(); ++x) {
	  BPTablePtr bp(new BPTable);
	  if (vm.count("aux")) {
	    bp->load(model[i++].c_str());
	  } else {
	    bp->pf_fold(*x);
	  }
	  bps.push_back(bp);
	}
	bp.average(bps, idxmaps);
	if (!vm.count("posteriors")) {
	  output(aln.name().front(), aln.seq().front(), bp, centroid, gamma);
	} else {
	  bp.save(std::cout, aln.seq().front(), p_th);
	}
      }
    } else {
      break;
    }
  }

  return 0;
}
