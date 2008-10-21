/*
 * $Id$
 *
 * CentroidFold: A generalized centroid estimator for predicting RNA
 *               secondary structures
 *
 * Copyright (C) 2008 Kengo Sato
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <iostream>
#include <string>
#include <cassert>
#include <boost/program_options.hpp>
#include <boost/range.hpp>
#include "centroid_fold.h"
#include "iupac.h"
#include "mea.h"
#include "centroid.h"
#include "fa.h"
#include "aln.h"

namespace po = boost::program_options;

int
main(int argc, char* argv[])
{

  std::vector<double> gamma, gamma_ali;
  std::string input;
  std::vector<std::string> model;
  float p_th=0.0;
  uint max_bp_dist;
  std::string param;
  
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
#ifdef HAVE_LIBCONTRAFOLD
    ("pf_fold", "use pf_fold base-pairing probabilities")
    ("params", po::value<std::string>(&param),
     "use the parameter file")
#endif
#endif
    ("aux", "use auxiliary base-pairing probabilities")
    ("max-dist", po::value<uint>(&max_bp_dist)->default_value(0),
      "the maximum distance of base-pairs")
    ("posteriors", po::value<float>(&p_th),
     "output base-pairing probability matrices which contain base-pairing probabilities more than the given value.")
#ifdef HAVE_LIBRNA
    ("ps", "draw predicted secondary structures with the postscript (PS) format")
    ("svg", "draw predicted secondary structures with the scalable vector graphics (SVG) format")
#endif
    ;
  po::options_description opts("Options");
  opts.add_options()
    ("seq-file", po::value<std::string>(&input), "training sequence filename")
    ("model-file", po::value<std::vector<std::string> >(&model), "model filename");

  opts.add(desc);
  po::positional_options_description pd;
  pd.add("seq-file", 1); pd.add("model-file", -1);
  po::parsed_options parsed =
    po::command_line_parser(argc, argv).options(opts).positional(pd).run();
  po::variables_map vm;
  po::store(parsed, vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("seq-file") ||
      (vm.count("aux") && model.empty()))
  {
    std::string features("aux files");
#ifdef HAVE_LIBRNA
    features += ", McCaskill model";
#endif
#ifdef HAVE_LIBCONTRAFOLD
    features += ", CONTRAfold model";
#endif
    std::cout << "CentroidFold v" << VERSION 
	      << " for predicting RNA secondary structures" << std::endl
	      << "  (enabled features: " << features << ")" << std::endl
	      << "Usage:" << std::endl
	      << " " << argv[0]
	      << " [options] seq [bp_matrix ...]\n\n"
	      << desc << std::endl;
    return 1;
  }

#ifdef HAVE_LIBCONTRAFOLD
  unsigned int engine = CentroidFold::CONTRAFOLD;
#else
  unsigned int engine = CentroidFold::PFFOLD;
#endif
  if (vm.count("aux")) engine = CentroidFold::AUX;
  if (vm.count("pf_fold")) engine = CentroidFold::PFFOLD;
  if (vm.count("alipf_fold")) engine = CentroidFold::ALIPFFOLD;

  if (gamma.size()==1 && gamma[0]<0.0)
  {
    double g[] = { 0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0,
		   8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0 };
    gamma.resize(boost::size(g));
    std::copy(boost::begin(g), boost::end(g), gamma.begin());
    gamma_ali.resize(boost::size(g));
    std::copy(boost::begin(g), boost::end(g), gamma_ali.begin());
  }
  if (gamma.empty())
  {
    switch (engine)
    {
#ifdef HAVE_LIBCONTRAFOLD
    case CentroidFold::CONTRAFOLD:
      gamma.push_back(vm.count("mea") ? 6.0 : 2.0);
      gamma_ali.push_back(vm.count("mea") ? 6.0 : 4.0);
      break;
#endif
    case CentroidFold::PFFOLD:
    case CentroidFold::ALIPFFOLD:
      gamma.push_back(vm.count("mea") ? 6.0 : 1.0);
      gamma_ali.push_back(vm.count("mea") ? 6.0 : 2.0);
      break;
    default:
      gamma.push_back(vm.count("mea") ? 6.0 : 1.0);
      gamma.push_back(vm.count("mea") ? 6.0 : 1.0);
      break;
    }
  }

  boost::spirit::file_iterator<> fi(input.c_str());
  if (!fi)
  {
    perror(input.c_str());
    return 1;
  }

  CentroidFold cf(engine, vm.count("mea"));
#ifdef HAVE_LIBCONTRAFOLD
  switch (engine) {
  case CentroidFold::CONTRAFOLD:
    cf.set_options(param, max_bp_dist);
    break;
  }
#endif
  
  while (1)
  {
    Fasta fa;
    Aln aln;
    if (fa.load(fi))
    {
      if (!vm.count("aux"))
        cf.calculate_posterior(fa.seq());
      else
        cf.calculate_posterior(fa.seq(), model[0]);
      if (!vm.count("posteriors"))
	cf.print(std::cout, fa.name(), fa.seq(), gamma);
      else
	cf.print_posterior(std::cout, fa.seq(), p_th);
      if (vm.count("ps")) cf.ps_plot(fa.name(), fa.seq(), gamma[0]);
      if (vm.count("svg")) cf.svg_plot(fa.name(), fa.seq(), gamma[0]);
    }
    else if (aln.load(fi))
    {
      if (vm.count("aux"))
      {
	if (aln.seq().size()!=model.size())
	{
	  std::cerr << "Given BP matrices and alignments were inconsistent.\n";
	  break;
	}
        cf.calculate_posterior(aln.seq(), model);
      } else {
        cf.calculate_posterior(aln.seq());
      }
      if (!vm.count("posteriors"))
	cf.print(std::cout, aln.name().front(), aln.consensus(), gamma_ali);
      else
	cf.print_posterior(std::cout, aln.consensus(), p_th);
      if (vm.count("ps")) cf.ps_plot(aln.name().front(), aln.consensus(), gamma_ali[0]);
      if (vm.count("svg")) cf.svg_plot(aln.name().front(), aln.consensus(), gamma_ali[0]);
    }
    else
    {
      break;
    }
  }

  return 0;
}
