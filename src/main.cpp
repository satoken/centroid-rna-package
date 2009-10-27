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
#include <fstream>
#include <string>
#include <cassert>
#include <boost/program_options.hpp>
#include <boost/range.hpp>
#include <boost/algorithm/string.hpp>
#include "centroid_fold.h"
#include "iupac.h"
#include "mea.h"
#include "centroid.h"
#include "fa.h"
#include "aln.h"

namespace po = boost::program_options;

int
centroid_fold_main(int argc, char* argv[])
{
  std::vector<float> gamma;
  std::vector<float> th;
  std::string input;
  std::vector<std::string> model;
  float p_th=0.0;
  std::string p_outname;
  std::string outname;
  std::string ps_outname;
  uint max_bp_dist;
  std::string param;
  uint max_clusters;
  uint num_samples;
  uint seed;
  //
  int num_ea_samples = -1;
  int max_mcc = -1;

  // parse command line options
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "show this message")
    ("gamma,g", po::value<std::vector<float> >(&gamma), "weight of base pairs")
    ("threshold,t", po::value<std::vector<float> >(&th),
     "thereshold of base pairs (this option overwrites 'gamma')")
    //
    ("ea", po::value<int>(&num_ea_samples), 
     "compute (pseudo-)expected accuracy (pseudo if arg==0, sampling if arg>0; arg: # of sampling)")
    ("max-mcc", po::value<int>(&max_mcc), 
     "predict secondary structure by maximizing pseudo-expected MCC (arg: # of sampling)")
    // added by M. Hamada
    ("mea", "run as an MEA estimator")
    ("noncanonical", "allow non-canonical base-pairs")
    ("aux", "use auxiliary base-pairing probabilities")
    ("constraints,C", "use structure constraints")
    ("output,o", po::value<std::string>(&outname),
     "specify filename to output predicted secondary structures. If empty, use the standard output.")
    ("posteriors", po::value<float>(&p_th),
     "output base-pairing probability matrices which contain base-pairing probabilities more than the given value.")
    ("posteriors-output", po::value<std::string>(&p_outname),
     "specify filename to output base-pairing probability matrices. If empty, use the standard output.")
    ("postscript", po::value<std::string>(&ps_outname),
     "draw predicted secondary structures with the postscript (PS) format")
    ("monochrome", "draw the postscript with monochrome");

  po::options_description opts_contrafold("Options for CONTRAfold model");
  opts_contrafold.add_options()
#ifdef HAVE_LIBRNA
    ("pf_fold", "use pf_fold base-pairing probabilities rather than those of CONTRAfold model")
#endif
    ("params", po::value<std::string>(&param), "use the parameter file")
    ("max-dist,d", po::value<uint>(&max_bp_dist)->default_value(0),
     "the maximum distance of base-pairs");

  po::options_description opts_sampling("Options for sampling");
  opts_sampling.add_options()
    ("sampling", 
     po::value<uint>(&num_samples),
     "specify the number of samples to be generated for each sequence")
    ("max-clusters,c",
     po::value<uint>(&max_clusters)->default_value(10),
     "the maximum number of clusters for the stochastic sampling algorithm")
    ("seed",
     po::value<uint>(&seed)->default_value(0),
     "specify the seed for the random number generator (set this automatically if seed=0)");
  
  po::options_description opts("Options");
  opts.add_options()
    ("seq-file", po::value<std::string>(&input), "training sequence filename")
    ("model-file", po::value<std::vector<std::string> >(&model), "model filename");

  opts.add(desc);
  opts.add(opts_contrafold);
  opts.add(opts_sampling);
  po::positional_options_description pd;
  pd.add("seq-file", 1); pd.add("model-file", -1);
  po::variables_map vm;
  bool usage=false;
  try {
    po::parsed_options parsed =
      po::command_line_parser(argc, argv).options(opts).positional(pd).run();
    po::store(parsed, vm);
    po::notify(vm);
  } catch (...) {
    usage=true;
  }

  if (usage || vm.count("help") || !vm.count("seq-file") ||
      (vm.count("aux") && model.empty()))
  {
    std::string features("aux files");
    features += ", CONTRAfold model";
#ifdef HAVE_LIBRNA
    features += ", McCaskill model";
#endif
    std::cout << "CentroidFold v" << VERSION 
	      << " for predicting RNA secondary structures" << std::endl
	      << "  (enabled features: " << features << ")" << std::endl
	      << "Usage:" << std::endl
	      << " " << argv[0]
	      << " [options] seq [bp_matrix ...]\n\n"
	      << desc << std::endl
              << opts_contrafold << std::endl
              << opts_sampling << std::endl;
    return 1;
  }

  unsigned int engine = CentroidFold::CONTRAFOLD;
  if (vm.count("aux")) engine = CentroidFold::AUX;
  if (vm.count("pf_fold")) engine = CentroidFold::PFFOLD;
  if (vm.count("alipf_fold")) engine = CentroidFold::ALIPFFOLD;

  if (th.size()>0)
  {
    gamma.resize(th.size());
    for (uint i=0; i!=th.size(); ++i)
      gamma[i] = 1.0/th[i]-1.0;
  }

  if (gamma.size()==1 && gamma[0]<0.0)
  {
    float g[] = { 0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0,
                  8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0 };
    gamma.resize(boost::size(g));
    std::copy(boost::begin(g), boost::end(g), gamma.begin());
  }
  if (gamma.empty())
  {
    switch (engine)
    {
    case CentroidFold::CONTRAFOLD:
      gamma.push_back(vm.count("mea") ? 6.0 : 2.0);
      break;
    case CentroidFold::PFFOLD:
    case CentroidFold::ALIPFFOLD:
      gamma.push_back(vm.count("mea") ? 6.0 : 1.0);
      break;
    default:
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

  CentroidFold cf(engine, vm.count("mea"), num_ea_samples, 0, seed);
  switch (engine) {
#ifdef HAVE_LIBRNA
  case CentroidFold::PFFOLD:
  case CentroidFold::ALIPFFOLD:
    cf.set_options_for_pf_fold(!vm.count("noncanonical"), max_bp_dist);
    break;
#endif
  case CentroidFold::CONTRAFOLD:
    cf.set_options_for_contrafold(param, !vm.count("noncanonical"), max_bp_dist);
    break;
  }

  std::ostream* out = &std::cout;
  if (vm.count("output"))
  {
    out = new std::ofstream(outname.c_str());
    if (out->fail())
    {
      perror(outname.c_str());
      delete out;
      return 1;
    }
  }

  std::ostream* p_out = &std::cout;
  if (vm.count("posteriors") && vm.count("posteriors-output") && !vm.count("sampling"))
  {
    p_out = new std::ofstream(p_outname.c_str());
    if (p_out->fail())
    {
      perror(p_outname.c_str());
      delete p_out;
      return 1;
    }
  }

  Fasta fa;
  uint n=0;
  while (fa.load(fi))
  {
    n++;
    if (max_mcc > 0) {
      cf.max_mcc_fold (fa.name(), fa.seq(), max_mcc, *out);
      continue;
    }    
    if (vm.count("sampling"))
    {
      cf.stochastic_fold(fa.name(), fa.seq(), num_samples, max_clusters,
                         gamma, *out, p_outname, p_th);
      continue;
    }
    if (!vm.count("aux"))
    {
      if (!vm.count("constraints"))
      {
        cf.calculate_posterior(fa.seq());
      }
      else if (fa.str().empty())
      {
        std::cout << "Input constraints:" << std::endl;
        std::string str;
        std::getline(std::cin, str);
        cf.calculate_posterior(fa.seq(), str);
      }
      else
      {
        cf.calculate_posterior(fa.seq(), fa.str());
      }
    }
    else
    {
      CentroidFold::BPTable bp;
      bp.load(model[0].c_str());
      cf.calculate_posterior(fa.seq(), bp);
    }
    cf.print(*out, fa.name(), fa.seq(), gamma);
    if (vm.count("posteriors")) cf.print_posterior(*p_out, fa.seq(), p_th);
    if (!ps_outname.empty())
    {
      char buf[PATH_MAX];
      if (n==1)
        strncpy(buf, ps_outname.c_str(), sizeof(buf));
      else
        snprintf(buf, sizeof(buf), "%s-%d", ps_outname.c_str(), n-1);
      cf.ps_plot(std::string(buf), fa.seq(), gamma[0], !vm.count("monochrome"));
    }
  }

  if (out != &std::cout) delete out;
  if (p_out != &std::cout) delete p_out;

  return 0;
}

int
centroid_alifold_main(int argc, char* argv[])
{
  std::vector<float> gamma;
  std::vector<float> th;
  std::string input;
  std::vector<std::string> model;
  float p_th=0.0;
  std::string p_outname;
  std::string outname;
  std::string ps_outname;
  uint max_bp_dist;
  std::string param;
  uint max_clusters;
  uint num_samples;
  uint seed;
  //
  int num_ea_samples = -1;
  int max_mcc = -1;
  
  // parse command line options
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "show this message")
    ("gamma,g", po::value<std::vector<float> >(&gamma), "weight of base pairs")
    ("threshold,t", po::value<std::vector<float> >(&th),
     "thereshold of base pairs (this option overwrites 'gamma')")
    //
    ("ea", po::value<int>(&num_ea_samples), 
     "compute (pseudo-)expected accuracy (pseudo if arg==0, sampling if arg>0; arg: # of sampling)")
    ("max-mcc", po::value<int>(&max_mcc), 
     "predict secondary structure by maximizing pseudo-expected MCC (arg: # of sampling)")
    // added by M. Hamada
    ("mea", "run as an MEA estimator")
    ("noncanonical", "allow non-canonical base-pairs")
    ("aux", "use auxiliary base-pairing probabilities")
    ("constraints,C", "use structure constraints")
    ("output,o", po::value<std::string>(&outname),
     "specify filename to output predicted secondary structures. If empty, use the standard output.")
    ("posteriors", po::value<float>(&p_th),
     "output base-pairing probability matrices which contain base-pairing probabilities more than the given value.")
    ("posteriors-output", po::value<std::string>(&p_outname),
     "specify filename to output base-pairing probability matrices. If empty, use the standard output.")
    ("postscript", po::value<std::string>(&ps_outname),
     "draw predicted secondary structures with the postscript (PS) format")
    ("monochrome", "draw the postscript with monochrome");

  po::options_description opts_contrafold("Options for CONTRAfold model");
  opts_contrafold.add_options()
#ifdef HAVE_LIBRNA
    ("alipf_fold", "use alipf_fold base-pairing probabilities rather than those of CONTRAfold model")
    ("pf_fold", "use pf_fold base-pairing probabilities rather than those of CONTRAfold model")
#endif
    ("params", po::value<std::string>(&param), "use the parameter file")
    ("max-dist,d", po::value<uint>(&max_bp_dist)->default_value(0),
     "the maximum distance of base-pairs");

  po::options_description opts_sampling("Options for sampling");
  opts_sampling.add_options()
    ("sampling", 
     po::value<uint>(&num_samples),
     "specify the number of samples to be generated for each sequence")
    ("max-clusters,c",
     po::value<uint>(&max_clusters)->default_value(10),
     "the maximum number of clusters for the stochastic sampling algorithm")
    ("seed",
     po::value<uint>(&seed)->default_value(0),
     "specify the seed for the random number generator (set this automatically if seed=0)");

  po::options_description opts("Options");
  opts.add_options()
    ("seq-file", po::value<std::string>(&input), "training sequence filename")
    ("model-file", po::value<std::vector<std::string> >(&model), "model filename");

  opts.add(desc);
  opts.add(opts_contrafold);
  opts.add(opts_sampling);
  po::positional_options_description pd;
  pd.add("seq-file", 1); pd.add("model-file", -1);
  po::variables_map vm;
  bool usage=false;
  try {
    po::parsed_options parsed =
      po::command_line_parser(argc, argv).options(opts).positional(pd).run();
    po::store(parsed, vm);
    po::notify(vm);
  } catch (...) {
    usage=true;
  }

  if (usage || vm.count("help") || !vm.count("seq-file") ||
      (vm.count("aux") && model.empty()))
  {
    std::string features("aux files");
    features += ", CONTRAfold model";
#ifdef HAVE_LIBRNA
    features += ", McCaskill model";
#endif
    std::cout << "CentroidAlifold v" << VERSION 
	      << " for predicting RNA secondary structures" << std::endl
	      << "  (enabled features: " << features << ")" << std::endl
	      << "Usage:" << std::endl
	      << " " << argv[0]
	      << " [options] seq [bp_matrix ...]\n\n"
	      << desc << std::endl
              << opts_contrafold << std::endl
              << opts_sampling << std::endl;
    return 1;
  }

  unsigned int engine = CentroidFold::CONTRAFOLD;
  if (vm.count("aux")) engine = CentroidFold::AUX;
  if (vm.count("pf_fold")) engine = CentroidFold::PFFOLD;
  if (vm.count("alipf_fold")) engine = CentroidFold::ALIPFFOLD;

  if (th.size()>0)
  {
    gamma.resize(th.size());
    for (uint i=0; i!=th.size(); ++i)
      gamma[i] = 1.0/th[i]-1.0;
  }

  if (gamma.size()==1 && gamma[0]<0.0)
  {
    float g[] = { 0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0,
                  8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0 };
    gamma.resize(boost::size(g));
    std::copy(boost::begin(g), boost::end(g), gamma.begin());
  }
  if (gamma.empty())
  {
    switch (engine)
    {
    case CentroidFold::CONTRAFOLD:
      gamma.push_back(vm.count("mea") ? 6.0 : 4.0);
      break;
    case CentroidFold::PFFOLD:
    case CentroidFold::ALIPFFOLD:
      gamma.push_back(vm.count("mea") ? 6.0 : 2.0);
      break;
    default:
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

  CentroidFold cf(engine, vm.count("mea"), num_ea_samples, seed);
  switch (engine) {
#ifdef HAVE_LIBRNA
  case CentroidFold::PFFOLD:
  case CentroidFold::ALIPFFOLD:
    cf.set_options_for_pf_fold(!vm.count("noncanonical"), max_bp_dist);
    break;
#endif
  case CentroidFold::CONTRAFOLD:
    cf.set_options_for_contrafold(param, !vm.count("noncanonical"), max_bp_dist);
    break;
  }

  std::ostream* out = &std::cout;
  if (vm.count("output"))
  {
    out = new std::ofstream(outname.c_str());
    if (out->fail())
    {
      perror(outname.c_str());
      delete out;
      return 1;
    }
  }

  std::ostream* p_out = &std::cout;
  if (vm.count("posteriors") && vm.count("posteriors-output") && !vm.count("sampling"))
  {
    p_out = new std::ofstream(p_outname.c_str());
    if (p_out->fail())
    {
      perror(p_outname.c_str());
      delete p_out;
      return 1;
    }
  }

  Aln aln;
  uint n=0;
  while (aln.load(fi))
  {
    n++;

    if (max_mcc > 0) {
      cf.max_mcc_fold (aln, max_mcc, *out);
      continue;
    }    
    if (vm.count("sampling"))
    {
      cf.stochastic_fold(aln, num_samples, max_clusters,
                         gamma, *out, p_outname, p_th);
      continue;
    }
    if (vm.count("aux"))
    {
      if (aln.seq().size()!=model.size())
      {
        std::cerr << "Given BP matrices and alignments were inconsistent.\n";
        break;
      }
      std::list<boost::shared_ptr<CentroidFold::BPTable> > bps;
      std::vector<std::string>::const_iterator m;
      for (m=model.begin(); m!=model.end(); ++m)
      {
        boost::shared_ptr<CentroidFold::BPTable> v(new CentroidFold::BPTable);
        v->load(m->c_str());
        bps.push_back(v);
      }
      cf.calculate_posterior(aln.seq(), bps);
    }
    else
    {
      if (!vm.count("constraints"))
      {
        cf.calculate_posterior(aln.seq());
      }
      else
      {
        std::cout << "Input constraints:" << std::endl;
        std::string str;
        std::getline(std::cin, str);
        cf.calculate_posterior(aln.seq(), str);
      }
    }
      
    cf.print(*out, aln, gamma);
    if (vm.count("posteriors")) cf.print_posterior(*p_out, aln.consensus(), p_th);
    if (!ps_outname.empty())
    {
      char buf[PATH_MAX];
      if (n==1)
        strncpy(buf, ps_outname.c_str(), sizeof(buf));
      else
        snprintf(buf, sizeof(buf), "%s-%d", ps_outname.c_str(), n-1);
      cf.ps_plot(std::string(buf), aln.consensus(), gamma[0], !vm.count("monochrome"));
    }
  }

  if (out != &std::cout) delete out;
  if (p_out != &std::cout) delete p_out;

  return 0;
}

int
main(int argc, char* argv[])
{
  if (strstr(argv[0], "alifold")==NULL)
    return centroid_fold_main(argc, argv);
  else
    return centroid_alifold_main(argc, argv);
}
