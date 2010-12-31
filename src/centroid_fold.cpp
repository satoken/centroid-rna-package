/*
 * $Id$
 *
 * CentroidFold: A generalized centroid estimator for predicting RNA
 *               secondary structures
 *
 * Copyright (C) 2008-2010 Kengo Sato
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

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cassert>
#include <boost/program_options.hpp>
#include <boost/range.hpp>
#include <boost/algorithm/string.hpp>
#include "folding_engine.h"
#include "mea.h"
#include "centroid.h"
#include "fa.h"
#include "aln.h"

#include "engine/contrafold.h"
#include "engine/mccaskill.h"
#include "engine/pfold.h"
#include "engine/mixture.h"
#include "engine/aux.h"

namespace po = boost::program_options;

int
centroid_fold(int argc, char* argv[])
{
  std::vector<float> gamma;
  std::vector<float> th;
  std::vector<std::string> engine;
  std::vector<float> mix_w;
  std::string input;
  std::vector<std::string> model;
  float p_th=0.0;
  std::string p_outname;
  std::string outname;
  std::string ps_outname;
  uint max_bp_dist;
  std::string param;
  uint max_clusters;
  uint num_samples=0;
  uint seed;
  //
  int num_ea_samples = -1;
  int max_mcc = -1;

  // parse command line options
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "show this message")
#ifdef HAVE_LIBRNA
    ("engine,e", po::value<std::vector<std::string> >(&engine),
     "specify the inference engine (default: \"McCaskill\")")
#else
    ("engine,e", po::value<std::vector<std::string> >(&engine),
     "specify the inference engine (default: \"CONTRAfold\")")
#endif
    ("mixture,w", po::value<std::vector<float> >(&mix_w), "mixture weights of inference engines")
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
    ("constraints,C", "use structure constraints")
    ("output,o", po::value<std::string>(&outname),
     "specify filename to output predicted secondary structures. If empty, use the standard output.")
    ("posteriors", po::value<float>(&p_th),
     "output base-pairing probability matrices which contain base-pairing probabilities more than the given value.")
    ("oposteriors", po::value<std::string>(&p_outname),
     "specify filename to output base-pairing probability matrices. If empty, use the standard output.")
    ("postscript", po::value<std::string>(&ps_outname),
     "draw predicted secondary structures with the postscript (PS) format")
    /*("monochrome", "draw the postscript with monochrome")*/
    ("params", po::value<std::string>(&param), "use the parameter file");

  po::options_description opts_contrafold("Options for CONTRAfold model");
  opts_contrafold.add_options()
#if 0 // HAVE_LIBRNA // move to hidden options
    ("pf_fold", "use pf_fold base-pairing probabilities rather than those of CONTRAfold model")
#endif
    ("max-dist,d", po::value<uint>(&max_bp_dist)->default_value(0),
     "the maximum distance of base-pairs");

  po::options_description opts_sampling("Options for sampling");
  opts_sampling.add_options()
    ("sampling,s", 
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
#ifdef HAVE_LIBRNA
    ("pf_fold", "use pf_fold base-pairing probabilities rather than those of CONTRAfold model")
#endif
    ("aux", "use auxiliary base-pairing probabilities")
    ("monochrome", "draw the postscript with monochrome")
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
    std::string features("CONTRAfold");
#ifdef HAVE_LIBRNA
    features += ", McCaskill";
#endif
    features += ", pfold";
    features += ", AUX";

    std::cout << "CentroidFold v" << VERSION 
	      << " for predicting RNA secondary structures" << std::endl
	      << "  (available engines: " << features << ")" << std::endl
	      << "Usage:" << std::endl
	      << " " << argv[0]
	      << " [options] seq [bp_matrix ...]\n\n"
	      << desc << std::endl
              << opts_contrafold << std::endl
              << opts_sampling << std::endl;
    return 1;
  }

  BOOST_SPIRIT_CLASSIC_NS::file_iterator<> fi(input.c_str());
  if (!fi)
  {
    perror(input.c_str());
    return 1;
  }

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

#ifdef HAVE_LIBRNA
  if (engine.empty()) engine.push_back("McCaskill");
#else
  if (engine.empty()) engine.push_back("CONTRAfold");  
#endif
  if (vm.count("pf_fold")) { engine.resize(1); engine[0]="McCaskill"; }
  if (vm.count("aux")) { engine.resize(1); engine[0]="AUX"; }

  FoldingEngine<std::string>* cf=NULL;
  std::vector<FoldingEngine<std::string>*> cf_list(engine.size(), NULL);
  for (uint i=0; i!=engine.size(); ++i)
  {
    if (engine[i]=="CONTRAfold")
    {
      if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 2.0);
      cf_list[i] = new CONTRAfoldModel(param, !vm.count("noncanonical"), max_bp_dist, seed, vm.count("mea"));
    }
#ifdef HAVE_LIBRNA
    else if (engine[i]=="McCaskill")
    {
      if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 1.0);
      cf_list[i] = new McCaskillModel(!vm.count("noncanonical"), max_bp_dist,
                                      param.empty() ? NULL : param.c_str(),
                                      seed, vm.count("mea"));
    }
#endif
    else if (engine[i]=="pfold")
    {
      if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 1.0);
      std::string pfold_bin_dir(getenv("PFOLD_BIN_DIR") ? getenv("PFOLD_BIN_DIR") : ".");
      std::string awk_bin(getenv("AWK_BIN") ? getenv("AWK_BIN") : "mawk");
      std::string sed_bin(getenv("SED_BIN") ? getenv("SED_BIN") : "sed");
      cf_list[i] = new PfoldModel<std::string>(pfold_bin_dir, awk_bin, sed_bin, vm.count("mea"));
    }
    else if (engine[i]=="AUX")
    {
      if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 1.0);
      cf_list[i] = new AuxModel(model, vm.count("mea"));
    }
    else
    {
      throw std::logic_error("unsupported inference engine");
    }
  }

  if (engine.size()==1)
    cf=cf_list[0];
  else
  {
    if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 1.0);
    std::vector<std::pair<FoldingEngine<std::string>*,float> > models;
    for (uint i=0; i!=engine.size(); ++i)
    {
      if (engine.size()!=mix_w.size())
        models.push_back(std::make_pair(cf_list[i], 1.0));
      else
        models.push_back(std::make_pair(cf_list[i], mix_w[i]));
    }
    cf = new MixtureModel<std::string>(models, vm.count("mea"));
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
  if (vm.count("posteriors") && vm.count("oposteriors") && !vm.count("sampling"))
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

    std::replace(fa.seq().begin(), fa.seq().end(), 't', 'u');
    std::replace(fa.seq().begin(), fa.seq().end(), 'T', 'U');

    if (vm.count("constraints"))
    {
      if (!fa.str().empty())
        cf->set_constraint(fa.str());
      else
      {
        std::cout << "Input constraints:" << std::endl;
        std::string str;
        std::getline(std::cin, str);
        cf->set_constraint(str);
      }
    }

    if (num_samples>0)
    {
      if (max_clusters>0)
        cf->stochastic_fold(fa.name(), fa.seq(), num_samples, gamma, max_clusters, *out, p_outname, p_th);
      else
        cf->stochastic_fold(fa.seq(), num_samples, *out);
      continue;
    }

    if (max_mcc>0)
    {
      cf->max_mcc_fold (fa.name(), fa.seq(), *out, max_mcc);
      continue;
    }

    if (num_ea_samples>=0)
      cf->centroid_fold(fa.name(), fa.seq(), gamma, *out, num_ea_samples);
    else
      cf->centroid_fold(fa.name(), fa.seq(), gamma, *out);

    if (vm.count("posteriors")) cf->get_bp().save(*p_out, fa.seq(), p_th);

    if (!ps_outname.empty())
    {
      char buf[PATH_MAX];
      if (n==1)
        strncpy(buf, ps_outname.c_str(), sizeof(buf));
      else
        snprintf(buf, sizeof(buf), "%s-%d", ps_outname.c_str(), n-1);
      cf->ps_plot(std::string(buf), fa.seq(), gamma[0], !vm.count("monochrome"));
    }
  }

  if (engine.size()!=1 && cf) delete cf;
  for (uint i=0; i!=cf_list.size(); ++i) if (cf_list[i]) delete cf_list[i];
  if (out!=&std::cout) delete out;
  if (p_out!=&std::cout) delete p_out;

  return 0;
}

int
main(int argc, char* argv[])
{
  try
  {
    return centroid_fold(argc, argv);
  }
  catch (std::logic_error e)
  {
    std::cerr << e.what() << std::endl;
  }
  catch (std::runtime_error e)
  {
    std::cerr << e.what() << std::endl;
  }
  return 1;
}
