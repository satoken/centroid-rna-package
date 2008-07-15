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

#include <string>
#include <cassert>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include "centroid_fold.h"
#include "mea.h"
#include "centroid.h"

#ifdef HAVE_LIBRNA
namespace Vienna {
extern "C" {
#include <ViennaRNA/PS_dot.h>
#include <ViennaRNA/aln_util.h>
};
};
#endif

CentroidFold::
CentroidFold(unsigned int engine, bool run_as_mea,
	     unsigned int reserved_size)
  : engine_(engine),
    mea_(run_as_mea),
    bp_(reserved_size)
{
}

CentroidFold::
~CentroidFold()
{
}

#if 0
void
CentroidFold::
execute(const std::string& seq, const std::vector<float>& gamma,
	std::vector<std::string>& paren, std::vector<float>& ea,
	const std::string& model)
{
  calculate_posterior(seq, model);
  paren.resize(gamma.size());
  ea.resize(gamma.size());
  for (uint i=0; i!=gamma.size(); ++i) {
    ea[i] = decode_structure(gamma[i], paren[i]);
  }
}

void
CentroidFold::
execute(const std::list<std::string>& seq, const std::vector<float>& gamma,
	std::vector<std::string>& paren, std::vector<float>& ea,
	const std::vector<std::string>& model)
{
  calculate_posterior(seq, model);
  paren.resize(gamma.size());
  ea.resize(gamma.size());
  for (uint i=0; i!=gamma.size(); ++i) {
    ea[i] = decode_structure(gamma[i], paren[i]);
  }
}

void
CentroidFold::
execute(const std::vector<std::string>& seq, const std::vector<float>& gamma,
	std::vector<std::string>& paren, std::vector<float>& ea,
	const std::vector<std::string>& model)
{
  calculate_posterior(seq, model);
  paren.resize(gamma.size());
  ea.resize(gamma.size());
  for (uint i=0; i!=gamma.size(); ++i) {
    ea[i] = decode_structure(gamma[i], paren[i]);
  }
}

void
CentroidFold::
execute(const std::list<std::string>& seq, const std::vector<float>& gamma,
	std::vector<std::string>& paren, std::vector<float>& ea)
{
  std::vector<std::string> model;
  execute(seq, gamma, paren, ea, model);
}

void
CentroidFold::
execute(const std::vector<std::string>& seq, const std::vector<float>& gamma,
	std::vector<std::string>& paren, std::vector<float>& ea)
{
  std::vector<std::string> model;
  execute(seq, gamma, paren, ea, model);
}
#endif

void
CentroidFold::
calculate_posterior(const std::string& seq, const std::string& model)
{
  std::string seq2(seq);
  std::replace(seq2.begin(), seq2.end(), 't', 'u');
  std::replace(seq2.begin(), seq2.end(), 'T', 'U');
  switch (engine_) {
  case AUX:
    bp_.load(model.c_str());
    break;
#ifdef HAVE_LIBRNA
  case PFFOLD:
    bp_.pf_fold(seq2);
    break;
#endif
#ifdef HAVE_LIBCONTRAFOLD
  case CONTRAFOLD:
    bp_.contra_fold(seq2);
    break;
#endif
  default:
    assert(!"never come here");
    break;
  }
}

void
CentroidFold::
calculate_posterior(const std::list<std::string>& seq,
		    const std::vector<std::string>& model)
{
  std::list<std::string> seq2;
  std::list<std::string>::const_iterator x;
  for (x=seq.begin(); x!=seq.end(); ++x) {
    seq2.push_back(*x);
    std::replace(seq2.back().begin(), seq2.back().end(), 't', 'u');
    std::replace(seq2.back().begin(), seq2.back().end(), 'T', 'U');
  }
#ifdef HAVE_LIBRNA
  if (engine_==ALIPFFOLD)  {
    bp_.alipf_fold(seq2);
  } else {
#endif
    typedef boost::shared_ptr<BPTable> BPTablePtr;
    std::list<BPTablePtr> bps;
    std::list<std::string> seqs;
    std::list<std::vector<uint> > idxmaps;
    BPTable::convert_to_raw_sequences(seq2, seqs, idxmaps);
    uint i;
    std::list<std::string>::const_iterator x;
    for (x=seqs.begin(), i=0; x!=seqs.end(); ++x, ++i) {
      BPTablePtr bpi(new BPTable);
      switch (engine_) {
      case AUX:
	bpi->load(model[i].c_str());
	break;
#ifdef HAVE_LIBRNA
      case PFFOLD:
	bpi->pf_fold(*x);
	break;
#endif
#ifdef HAVE_LIBCONTRAFOLD
      case CONTRAFOLD:
	bpi->contra_fold(*x);
	break;
#endif
      default:
	assert(!"never come here");
	break;
      }
      bps.push_back(bpi);
    }
    bp_.average(bps, idxmaps);
#ifdef HAVE_LIBRNA
  }
#endif
}

void
CentroidFold::
calculate_posterior(const std::vector<std::string>& seq,
		    const std::vector<std::string>& model)
{
  std::list<std::string> s(seq.size());
  std::copy(seq.begin(), seq.end(), s.begin());
  calculate_posterior(s, model);
}

void
CentroidFold::
calculate_posterior(const std::list<std::string>& seq)
{
  const std::vector<std::string> model;
  calculate_posterior(seq, model);
}

void
CentroidFold::
calculate_posterior(const std::vector<std::string>& seq)
{
  const std::vector<std::string> model;
  calculate_posterior(seq, model);
}

float
CentroidFold::
decode_structure(float gamma, std::string& paren) const
{
  float p=0.0;
  paren.resize(bp_.size());
  std::fill(paren.begin(), paren.end(), '.');
  if (!mea_) {
    p = SCFG::Centroid::execute(bp_, paren, gamma);
  } else {
    p = SCFG::MEA::execute(bp_, paren, gamma);
  }
  return p;
}

std::pair<std::string,float>
CentroidFold::
decode_structure(float gamma) const
{
  std::string paren;
  float p = decode_structure(gamma, paren);
  return std::make_pair(paren, p);
}

void
CentroidFold::
print(std::ostream& out, const std::string& name, const std::string& seq,
      const std::vector<double>& gamma) const
{
  out << ">" << name << std::endl << seq << std::endl;
  std::vector<double>::const_iterator g;
  for (g=gamma.begin(); g!=gamma.end(); ++g) {
    std::string paren;
    float p = decode_structure(*g, paren);
    out << paren;
    if (!mea_)
      out << " (g=" << *g << ",th=" << (1.0/(1.0+*g)) << ")";
    else
      out << " (g=" << *g << ",EA=" << p << ")";
    out << std::endl;
  }
}

void
CentroidFold::
print_posterior(std::ostream& out, const std::string& seq, float th) const
{
  bp_.save(out, seq, th);
}

#ifdef HAVE_LIBRNA
void
CentroidFold::
ps_plot(const std::string& name, const std::string& seq, float g) const
{
  std::string paren;
  decode_structure(g, paren);
  std::string fbase(name);
  boost::algorithm::trim(fbase);
  std::replace(boost::begin(fbase), boost::end(fbase), ' ', '_');
  std::replace(boost::begin(fbase), boost::end(fbase), '/', '_');
  if (fbase.empty()) fbase="rna";
  char fname[100];
  sscanf(fbase.c_str(), "%12s", fname);
  //Vienna::rna_plot_type=0;
  strcat(fname, "_ss.ps");
  Vienna::PS_rna_plot(const_cast<char*>(seq.c_str()),
		      const_cast<char*>(paren.c_str()), fname);
}

void
CentroidFold::
svg_plot(const std::string& name, const std::string& seq, float g) const
{
  std::string paren;
  decode_structure(g, paren);
  std::string fbase(name);
  boost::algorithm::trim(fbase);
  std::replace(boost::begin(fbase), boost::end(fbase), ' ', '_');
  std::replace(boost::begin(fbase), boost::end(fbase), '/', '_');
  if (fbase.empty()) fbase="rna";
  char fname[100];
  sscanf(fbase.c_str(), "%12s", fname);
  //Vienna::rna_plot_type=0;
  strcat(fname, "_ss.ps");
  Vienna::svg_rna_plot(const_cast<char*>(seq.c_str()),
		       const_cast<char*>(paren.c_str()), fname);
}
#endif
