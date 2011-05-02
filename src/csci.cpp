/*
 * $Id: centroid_alifold.cpp 110 2011-04-26 11:31:27Z sato-kengo $
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
#include "engine/contrafoldm.h"
#include "engine/mccaskill.h"
#include "engine/alifold.h"
#include "engine/pfold.h"
#include "engine/averaged.h"
#include "engine/mixture.h"
#include "engine/aux.h"

namespace po = boost::program_options;

int
csci(int argc, char* argv[])
{
  std::vector<float> gamma_a, gamma_s;
  std::vector<std::string> engine;
  std::vector<float> mix_w;
  std::string input;
  std::vector<std::string> model;
  uint max_bp_dist;
  std::string param;
  
  // parse command line options
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "show this message")
#ifdef HAVE_LIBRNA
    ("engine,e", po::value<std::vector<std::string> >(&engine),
     "specify the inference engine (default: \"McCaskill & Alifold\")")
#else
    ("engine,e", po::value<std::vector<std::string> >(&engine),
     "specify the inference engine (default: \"CONTRAfold\")")
#endif
    ("mixture,w", po::value<std::vector<float> >(&mix_w), "mixture weights of inference engines")
    ("gamma_a", po::value<std::vector<float> >(&gamma_a), "weight of base pairs for alignments")
    ("gamma_s", po::value<std::vector<float> >(&gamma_s), "weight of base pairs for individual sequences")
    //
    ("mea", "run as an MEA estimator")
    ("noncanonical", "allow non-canonical base-pairs")
    ("params", po::value<std::string>(&param), "use the parameter file");

  po::options_description opts_contrafold("Options for CONTRAfold model");
  opts_contrafold.add_options()
#if 0 // HAVE_LIBRNA // move to hidden options
    ("alipf_fold", "use alipf_fold base-pairing probabilities rather than those of CONTRAfold model")
    ("pf_fold", "use pf_fold base-pairing probabilities rather than those of CONTRAfold model")
#endif
    ("max-dist,d", po::value<uint>(&max_bp_dist)->default_value(0),
     "the maximum distance of base-pairs");

  po::options_description opts("Options");
  opts.add_options()
#ifdef HAVE_LIBRNA
    ("alipf_fold", "use alipf_fold base-pairing probabilities rather than those of CONTRAfold model")
    ("pf_fold", "use pf_fold base-pairing probabilities rather than those of CONTRAfold model")
#endif
    ("aux", "use auxiliary base-pairing probabilities")
    ("seq-file", po::value<std::string>(&input), "training sequence filename")
    ("model-file", po::value<std::vector<std::string> >(&model), "model filename");

  opts.add(desc);
  opts.add(opts_contrafold);
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
    features += ", Alifold";
#endif
    features += ", pfold";
    features += ", AUX";

    std::cout << "CSCI v" << VERSION << std::endl
	      << "  (available engines: " << features << ")" << std::endl
	      << "Usage:" << std::endl
	      << " " << argv[0]
	      << " [options] seq [bp_matrix ...]\n\n"
	      << desc << std::endl
              << opts_contrafold << std::endl;
    return 1;
  }

  BOOST_SPIRIT_CLASSIC_NS::file_iterator<> fi(input.c_str());
  if (!fi)
  {
    perror(input.c_str());
    return 1;
  }

  if (engine.empty())
  {
#ifdef HAVE_LIBRNA
    engine.push_back("McCaskill");
    engine.push_back("Alifold");
#else
    engine.push_back("CONTRAfold");
#endif
  }
  if (vm.count("pf_fold")) { engine.resize(1); engine[0]="McCaskill"; }
  if (vm.count("alipf_fold")) { engine.resize(1); engine[0]="Alifold"; }
  if (vm.count("aux")) { engine.resize(1); engine[0]="AUX"; }

  FoldingEngine<Aln>* cf=NULL;
  std::vector<FoldingEngine<Aln>*> cf_list(engine.size(), NULL);
  std::vector<FoldingEngine<std::string>*> src_list(engine.size(), NULL);
  for (uint i=0; i!=engine.size(); ++i)
  {
    if (engine[i]=="CONTRAfold")
    {
      src_list[i] = new CONTRAfoldModel(param, !vm.count("noncanonical"), max_bp_dist, 0, vm.count("mea"));
      cf_list[i] = new AveragedModel(src_list[i], max_bp_dist, vm.count("mea"));
    }
    else if (engine[i]=="CONTRAfoldM")
    {
      cf_list[i] = new CONTRAfoldMultiModel(param, !vm.count("noncanonical"), max_bp_dist, 0, vm.count("mea"));
    }
#ifdef HAVE_LIBRNA
    else if (engine[i]=="McCaskill" || engine[i]=="McCaskillVienna")
    {
      src_list[i] = new McCaskillModel(!vm.count("noncanonical"), max_bp_dist,
                                       engine[i]!="McCaskillVienna",
                                       param.empty() ? NULL : param.c_str(),
                                       0, vm.count("mea"));
      cf_list[i] = new AveragedModel(src_list[i], max_bp_dist, vm.count("mea"));
    }
    else if (engine[i]=="Alifold" || engine[i]=="AlifoldVienna")
    {
      cf_list[i] = new AliFoldModel(!vm.count("noncanonical"), max_bp_dist,
                                    engine[i]!="AlifoldVienna",
                                    param.empty() ? NULL : param.c_str(),
                                    0, vm.count("mea"));
    }
#endif
    else if (engine[i]=="pfold")
    {
      std::string pfold_bin_dir(getenv("PFOLD_BIN_DIR") ? getenv("PFOLD_BIN_DIR") : ".");
      std::string awk_bin(getenv("AWK_BIN") ? getenv("AWK_BIN") : "mawk");
      std::string sed_bin(getenv("SED_BIN") ? getenv("SED_BIN") : "sed");
      cf_list[i] = new PfoldModel<Aln>(pfold_bin_dir, awk_bin, sed_bin, vm.count("mea"));
    }
    else if (engine[i]=="AUX")
    {
      src_list[i] = new AuxModel(model, vm.count("mea"));
      cf_list[i] = new AveragedModel(src_list[i], 0, vm.count("mea"));
    }
    else
    {
      throw std::logic_error("unsupported engine");
    }
  }

  if (engine.size()==1)
    cf=cf_list[0];
  else
  {
    std::vector<std::pair<FoldingEngine<Aln>*,float> > models;
    for (uint i=0; i!=engine.size(); ++i)
    {
      if (engine.size()!=mix_w.size())
        models.push_back(std::make_pair(cf_list[i], 1.0));
      else
        models.push_back(std::make_pair(cf_list[i], mix_w[i]));
    }
    cf = new MixtureModel<Aln>(models, vm.count("mea"));
  }

  if (gamma_a.empty()) gamma_a.push_back(1.0);
  if (gamma_s.empty()) gamma_s.push_back(1.0);
  
  Aln aln;
  uint n=0;
  uint bytes=0;
  uint total_bytes=0;
  while ((bytes=aln.load(fi))>0)
  {
    n++;
    total_bytes+=bytes;

    std::list<std::string>::iterator seq;
    for (seq=aln.seq().begin(); seq!=aln.seq().end(); ++seq)
    {
      std::replace(seq->begin(), seq->end(), 't', 'u');
      std::replace(seq->begin(), seq->end(), 'T', 'U');
    }

    std::vector<float> r = cf->csci(aln, gamma_a, gamma_s);
    std::copy(r.begin(), r.end(), std::ostream_iterator<float>(std::cout, ","));
    std::cout << std::endl;
  }
  if (fi!=fi.make_end())
    std::cout << "parse error after " << total_bytes << " bytes were loaded" << std::endl;

  if (engine.size()!=1 && cf) delete cf;
  for (uint i=0; i!=cf_list.size(); ++i) if (cf_list[i]) delete cf_list[i];
  for (uint i=0; i!=src_list.size(); ++i) if (src_list[i]) delete src_list[i];

  return 0;
}

int
main(int argc, char* argv[])
{
  try
  {
    return csci(argc, argv);
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
