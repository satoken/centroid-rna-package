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

#include "engine/contrafoldhom.h"
#include "engine/mccaskillhom.h"

namespace po = boost::program_options;

int
centroid_homfold(int argc, char* argv[])
{
  std::vector<float> gamma;
  std::vector<float> th;
  std::vector<std::string> engine; // for secondary structure
  std::vector<std::string> engine_a; // for alignment
  std::vector<float> mix_w;
  std::string input;
  std::string hom_seqs = "";
  std::vector<std::string> model;
  float p_th=0.0;
  std::string p_outname;
  std::string outname;
  std::string ps_outname;
  uint max_bp_dist=0;
  std::string param;
  uint seed;

  // parse command line options
  po::options_description desc("Options");
  desc.add_options()
    ("homologous,H", po::value<std::string>(&hom_seqs), "fasta file containing homologous sequences (REQUIRED)")
    ("help,h", "show this message")
#ifdef HAVE_LIBRNA
    ("engine_s", po::value<std::vector<std::string> >(&engine),
     "specify the inference engine for secondary structures (default: \"McCaskill\")")
#else
    ("engine_s", po::value<std::vector<std::string> >(&engine),
     "specify the inference engine for secondary structures (default: \"CONTRAfold\")")
#endif
    ("engine_a", po::value<std::vector<std::string> >(&engine_a),
     "specify the inference engine for pairwise alignments (default: \"CONTRAlign\")")
    //("mixture,w", po::value<std::vector<float> >(&mix_w), "mixture weights of inference engines")
    ("gamma,g", po::value<std::vector<float> >(&gamma), "weight of base pairs")
    ("threshold,t", po::value<std::vector<float> >(&th),
     "thereshold of base pairs (this option overwrites 'gamma')")
    //
    ("ea", "compute (pseudo-)expected accuracy")
    //("max-mcc", po::value<int>(&max_mcc), 
    // "predict secondary structure by maximizing pseudo-expected MCC (arg: # of sampling)")
    //("mea", "run as an MEA estimator")
    //("noncanonical", "allow non-canonical base-pairs")
    //("constraints,C", "use structure constraints")
    ("output,o", po::value<std::string>(&outname),
     "specify filename to output predicted secondary structures. If empty, use the standard output.")
    ("posteriors", po::value<float>(&p_th),
     "output base-pairing probability matrices which contain base-pairing probabilities more than the given value.")
    ("oposteriors", po::value<std::string>(&p_outname),
     "specify filename to output base-pairing probability matrices. If empty, use the standard output.")
    ("postscript", po::value<std::string>(&ps_outname),
     "draw predicted secondary structures with the postscript (PS) format")
    ("params", po::value<std::string>(&param), "use the parameter file");

  po::options_description opts("Options");
  opts.add_options()
#ifdef HAVE_LIBRNA
    ("pf_fold", "use pf_fold base-pairing probabilities rather than those of CONTRAfold model")
#endif
    ("monochrome", "draw the postscript with monochrome")
    ("seq-file", po::value<std::string>(&input), "training sequence filename")
    ("model-file", po::value<std::vector<std::string> >(&model), "model filename");

  opts.add(desc);
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

  if (usage || vm.count("help") || !vm.count("seq-file") || hom_seqs == "" ||
      (vm.count("aux") && model.empty()))
  {
    std::string features("CONTRAfold");
#ifdef HAVE_LIBRNA
    features += ", McCaskill";
#endif

    std::cout << "CentroidHomfold v" << VERSION 
	      << " for predicting RNA secondary structures" << std::endl
	      << "  (available engines: " << features << ")" << std::endl
	      << "Usage:" << std::endl
	      << " " << argv[0]
	      << " [options] seq [bp_matrix ...]\n\n"
	      << desc << std::endl;
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
                  8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0};
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

  if (engine_a.empty()) engine_a.push_back("CONTRAlign"); 

  std::vector<std::string> homs;
  if (hom_seqs != "") {
    BOOST_SPIRIT_CLASSIC_NS::file_iterator<> fi2(hom_seqs.c_str());
    if (!fi2) {
      perror(hom_seqs.c_str());
      return 1;
    }
    while (1) {
      Fasta fa;
      if (fa.load(fi2)) {
	homs.push_back (fa.seq());
      } else break;
    }
  }  

  FoldingEngine<TH>* cf=NULL;
  std::vector<FoldingEngine<TH>*> cf_list(engine.size(), NULL);
  for (uint i=0; i!=engine.size(); ++i)
  {
    if (engine[i]=="CONTRAfold")
    {
      if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 2.0);
      cf_list[i] = new CONTRAfoldHomModel(param, engine_a[0], !vm.count("noncanonical"), max_bp_dist, seed, vm.count("mea"));
    }
#ifdef HAVE_LIBRNA
    else if (engine[i]=="McCaskill")
    {
      if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 1.0);
      cf_list[i] = new McCaskillHomModel(engine_a[0], !vm.count("noncanonical"), max_bp_dist,
					 param.empty() ? NULL : param.c_str(),
                                         seed, vm.count("mea"));
    }
#endif
    else
    {
      std::cerr << engine[i] << std::endl;
      throw std::logic_error("unsupported inference engine");
    }
  }

  if (engine.size()==1)
    cf=cf_list[0];

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

    if (vm.count("ea"))
      cf->centroid_fold(fa.name(), TH (fa.seq(), homs), gamma, *out, 0);
    else
      cf->centroid_fold(fa.name(), TH (fa.seq(), homs), gamma, *out);

    if (vm.count("posteriors")) cf->get_bp().save(*p_out, fa.seq(), p_th);

    if (!ps_outname.empty())
    {
      char buf[PATH_MAX];
      if (n==1)
        strncpy(buf, ps_outname.c_str(), sizeof(buf));
      else
        snprintf(buf, sizeof(buf), "%s-%d", ps_outname.c_str(), n-1);
      cf->ps_plot(std::string(buf), TH(fa.seq(), homs), gamma[0], !vm.count("monochrome"));
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
    return centroid_homfold(argc, argv);
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
