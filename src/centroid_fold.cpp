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

#include <cmath>
#include <string>
#include <stack>
#include <cassert>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include "centroid_fold.h"
#include "mea.h"
#include "centroid.h"

#ifdef HAVE_LIBRNA
namespace Vienna {
extern "C" {
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/alifold.h>
#include <ViennaRNA/PS_dot.h>
#include <ViennaRNA/aln_util.h>
#if 0
  extern int pfl_fold(char *sequence, int winSize, int pairdist,
		      float cutoff, struct plist **pl);
  extern void init_pf_foldLP(int length);
  extern void free_pf_arraysLP(void);
#endif
  extern char* pbacktrack(char *sequence);
  extern int   st_back;
};
};
#endif

#ifdef HAVE_LIBCONTRAFOLD
#include <contrafold.h>
#endif

// folding routines
#ifdef HAVE_LIBRNA
template < class T >
static
void
pf_fold(T& bp, const std::string& seq, const std::string& str="")
{
  bp.resize(seq.size());
  Vienna::pf_scale = -1;
  Vienna::init_pf_fold(seq.size());
  if (str.empty()) {
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  } else {
    char *str2 = new char[str.size()+1];
    strcpy(str2, str.c_str());
    int bk = Vienna::fold_constrained;
    Vienna::fold_constrained = 1;
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), str2);
    delete[] str2;
    Vienna::fold_constrained = bk;
  }
  for (uint j=2; j!=bp.size()+1; ++j) {
    for (uint i=j-1; ; --i) {
      bp.update(i-1, j-1, Vienna::pr[Vienna::iindx[i]-j]);
      if (i==1) break;
    }
  }
  Vienna::free_pf_arrays();
}

template < class T >
static
void
alipf_fold(T& bp, const std::list<std::string>& ma, const std::string& str="")
{
  bp.resize(ma.front().size());
  // prepare an alignment
  uint length = ma.front().size();
  char** seqs = new char*[ma.size()+1];
  seqs[ma.size()]=NULL;
  std::list<std::string>::const_iterator x;
  uint i=0;
  for (x=ma.begin(); x!=ma.end(); ++x) {
    assert(x->size()==length);
    seqs[i] = new char[length+1];
    strcpy(seqs[i++], x->c_str());
  }

  char* str2 = new char[length+1];
  int bk = Vienna::fold_constrained;
  if (!str.empty()) Vienna::fold_constrained=1;
  {
    // scaling parameters to avoid overflow
    strcpy(str2, str.c_str());
    double min_en = Vienna::alifold(seqs, str2);
    Vienna::free_alifold_arrays();
    double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
    Vienna::pf_scale = exp(-(1.07*min_en)/kT/length);
  }
  // build a base pair probablity matrix
  Vienna::pair_info* pi;
  strcpy(str2, str.c_str());
  Vienna::alipf_fold(seqs, str2, &pi);
  for (uint k=0; pi[k].i!=0; ++k)
    bp.update(pi[k].i-1, pi[k].j-1, pi[k].p);
  free(pi);
  for (uint i=0; seqs[i]!=NULL; ++i) delete[] seqs[i];

  Vienna::fold_constrained = bk;
  delete[] seqs;
  delete[] str2;
}
#endif

#ifdef HAVE_LIBCONTRAFOLD
template < class T, class U >
static
void
contra_fold(T& bp, CONTRAfold<U>& cf, const std::string& seq, const std::string& str="")
{
  bp.resize(seq.size(), cf.max_bp_dist());
  std::vector<float> posterior;
  cf.SetConstraint(str);
  cf.ComputePosterior(seq, posterior);

  if (cf.max_bp_dist()==0) {
    uint k=0;
    for (uint i=0; i!=seq.size()+1; ++i) {
      for (uint j=i; j!=seq.size()+1; ++j) {
        if (i!=0) bp.update(i-1, j-1, posterior[k]);
        ++k;
      }
    }
  } else {
    uint k=0;
    for (uint i=0; i!=seq.size()+1; ++i) {
      for (uint j=i; j!=i+cf.max_bp_dist(); ++j) {
        if (i!=0 && j<=seq.size()) bp.update(i-1, j-1, posterior[k]);
        ++k;
      }
    }
  }
}
#endif


CentroidFold::
CentroidFold(unsigned int engine, bool run_as_mea,
	     unsigned int reserved_size)
  : engine_(engine),
    mea_(run_as_mea),
    bp_(reserved_size)
#ifdef HAVE_LIBCONTRAFOLD
  ,
    contrafold_(),
    model_(),
    canonical_only_(true),
    max_bp_dist_(0)
#endif
{
}

CentroidFold::
~CentroidFold()
{
}

void
CentroidFold::
calculate_posterior(const std::string& seq)
{
  std::string seq2(seq);
  std::replace(seq2.begin(), seq2.end(), 't', 'u');
  std::replace(seq2.begin(), seq2.end(), 'T', 'U');
  switch (engine_) {
  case AUX:
    assert(!"AUX should be given a bp matrix.");
    break;
#ifdef HAVE_LIBRNA
  case PFFOLD:
    assert(max_bp_dist_==0);
    pf_fold(bp_, seq2);
    break;
#endif
#ifdef HAVE_LIBCONTRAFOLD
  case CONTRAFOLD:
    if (contrafold_==NULL) {
      contrafold_ =
        boost::shared_ptr<CONTRAfold<float> >(new CONTRAfold<float>(canonical_only_, max_bp_dist_));
      if (!model_.empty()) contrafold_->SetParameters(model_);
    }
    contra_fold(bp_, *contrafold_, seq2);
    break;
#endif
  default:
    assert(!"never come here");
    break;
  }
}

void
CentroidFold::
calculate_posterior(const std::string& seq, const std::string& str)
{
  std::string seq2(seq);
  std::replace(seq2.begin(), seq2.end(), 't', 'u');
  std::replace(seq2.begin(), seq2.end(), 'T', 'U');
  switch (engine_) {
  case AUX:
    assert(!"AUX should be given a bp matrix.");
    break;
#ifdef HAVE_LIBRNA
  case PFFOLD:
    assert(max_bp_dist_==0);
    pf_fold(bp_, seq2, str);
    break;
#endif
#ifdef HAVE_LIBCONTRAFOLD
  case CONTRAFOLD:
    if (contrafold_==NULL) {
      contrafold_ =
        boost::shared_ptr<CONTRAfold<float> >(new CONTRAfold<float>(canonical_only_, max_bp_dist_));
      if (!model_.empty()) contrafold_->SetParameters(model_);
    }
    contra_fold(bp_, *contrafold_, seq2, str);
    break;
#endif
  default:
    assert(!"never come here");
    break;
  }
}

void
CentroidFold::
calculate_posterior(const std::string& seq, const BPTable& bp)
{
  std::string seq2(seq);
  std::replace(seq2.begin(), seq2.end(), 't', 'u');
  std::replace(seq2.begin(), seq2.end(), 'T', 'U');
  switch (engine_) {
  case AUX:
    assert(max_bp_dist_==0);
    bp_ = bp;
    break;
  default:
    assert(!"unreachable");
    break;
  }
}

void
CentroidFold::
calculate_posterior(const std::list<std::string>& seq)
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
    assert(max_bp_dist_==0);
    alipf_fold(bp_, seq2);
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
        assert(!"AUX should be given bp matrices.");
	break;
#ifdef HAVE_LIBRNA
      case PFFOLD:
        assert(max_bp_dist_==0);
	pf_fold(*bpi, *x);
	break;
#endif
#ifdef HAVE_LIBCONTRAFOLD
      case CONTRAFOLD:
        if (contrafold_==NULL) {
          contrafold_ =
            boost::shared_ptr<CONTRAfold<float> >(new CONTRAfold<float>(canonical_only_, max_bp_dist_));
          if (!model_.empty()) contrafold_->SetParameters(model_);
        }
        contra_fold(*bpi, *contrafold_, *x);
	break;
#endif
      default:
	assert(!"never come here");
	break;
      }
      bps.push_back(bpi);
    }
    bp_.average(bps, idxmaps, max_bp_dist_);
#ifdef HAVE_LIBRNA
  }
#endif
}

void
CentroidFold::
calculate_posterior(const std::list<std::string>& seq, const std::string& str)
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
    assert(max_bp_dist_==0);
    alipf_fold(bp_, seq2, str);
  } else {
#endif
    // make an alignment-to-sequence map
    typedef boost::shared_ptr<BPTable> BPTablePtr;
    std::list<BPTablePtr> bps;
    std::list<std::string> seqs;
    std::list<std::vector<uint> > idxmaps;
    BPTable::convert_to_raw_sequences(seq2, seqs, idxmaps);

    // make a base-pair map
    std::vector<uint> bpmap(idxmaps.front().size(), static_cast<uint>(-1));
    std::stack<uint> stk;
    for (uint i=0; i!=str.size(); ++i) {
      if (str[i]=='(') {
        stk.push(i);
      } else if (str[i]==')') {
        bpmap[i]=stk.top();
        bpmap[stk.top()]=i;
        stk.pop();
      }
    }

    std::list<std::string>::const_iterator x;
    std::list<std::vector<uint> >::const_iterator y;
    for (x=seqs.begin(), y=idxmaps.begin(); x!=seqs.end(); ++x, ++y) {
      // map the consensus structure for constraints
      char unknow_str = engine_==CONTRAFOLD ? '?' : '.';
      BPTablePtr bpi(new BPTable);
      std::string str2(x->size(), ' ');
      for (uint i=0, j=0; i!=y->size(); ++i) {
        if ((*y)[i]!=static_cast<uint>(-1)) {
          if ((*y)[bpmap[i]]!=static_cast<uint>(-1))
            str2[j++] = str[i];
          else
            str2[j++] = unknow_str;
        }
      }
      
      switch (engine_) {
      case AUX:
        assert(!"AUX should be given bp matrices.");
	break;
#ifdef HAVE_LIBRNA
      case PFFOLD:
        assert(max_bp_dist_==0);
	pf_fold(*bpi, *x, str2);
	break;
#endif
#ifdef HAVE_LIBCONTRAFOLD
      case CONTRAFOLD:
        if (contrafold_==NULL) {
          contrafold_ =
            boost::shared_ptr<CONTRAfold<float> >(new CONTRAfold<float>(canonical_only_, max_bp_dist_));
          if (!model_.empty()) contrafold_->SetParameters(model_);
        }
        contra_fold(*bpi, *contrafold_, *x, str2);
	break;
#endif
      default:
	assert(!"never come here");
	break;
      }
      bps.push_back(bpi);
    }
    bp_.average(bps, idxmaps, max_bp_dist_);
#ifdef HAVE_LIBRNA
  }
#endif
}

void
CentroidFold::
calculate_posterior(const std::list<std::string>& seq,
		    const std::list<boost::shared_ptr<BPTable> >& bps)
{
  std::list<std::string> seq2;
  std::list<std::string>::const_iterator x;
  for (x=seq.begin(); x!=seq.end(); ++x) {
    seq2.push_back(*x);
    std::replace(seq2.back().begin(), seq2.back().end(), 't', 'u');
    std::replace(seq2.back().begin(), seq2.back().end(), 'T', 'U');
  }
  typedef boost::shared_ptr<BPTable> BPTablePtr;
  std::list<std::string> seqs;
  std::list<std::vector<uint> > idxmaps;
  BPTable::convert_to_raw_sequences(seq2, seqs, idxmaps);
  switch (engine_) {
  case AUX:
    assert(max_bp_dist_==0);
    bp_.average(bps, idxmaps);
    break;
  default:
    assert(!"unreachable");
    break;
  }
}

float
CentroidFold::
decode_structure(float gamma, std::string& paren) const
{
  float p=0.0;
  paren.resize(bp_.size());
  std::fill(paren.begin(), paren.end(), '.');
  if (!mea_) {
    if (max_bp_dist_==0)
      p = SCFG::Centroid::execute(bp_, paren, gamma);
    else
      p = SCFG::Centroid::execute(bp_, paren, max_bp_dist_, gamma);
  } else {
    if (max_bp_dist_==0)
      p = SCFG::MEA::execute(bp_, paren, gamma);
    else
      p = SCFG::MEA::execute(bp_, paren, max_bp_dist_, gamma);
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
