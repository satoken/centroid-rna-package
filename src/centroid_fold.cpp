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
#include <sstream>
#include <cmath>
#include <string>
#include <stack>
#include <cstdio>
#include <cassert>
#include <sys/time.h>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/multi_array.hpp>
#include <boost/bind.hpp>
#include "centroid_fold.h"
#include "mea.h"
#include "centroid.h"
#include "diana.h"
#include "ps_plot.h"

#ifdef HAVE_LIBRNA
namespace Vienna {
extern "C" {
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/alifold.h>
#include <ViennaRNA/LPfold.h>
#include <ViennaRNA/PS_dot.h>
#include <ViennaRNA/aln_util.h>
#include <ViennaRNA/utils.h>
  extern int   st_back;
};
};
#endif

#include "contrafold/contrafold.h"

typedef boost::dynamic_bitset<> BPvec;
typedef boost::shared_ptr<BPvec> BPvecPtr;

inline
uint
hamming_distance(const BPvec& x, const BPvec& y)
{
  return (x^y).count();
}

struct EncodeBP
{
  EncodeBP()
    : t_(NULL), vec_(), p_(0)
  {}

  BPvecPtr execute(const SCFG::CYKTable<uint>& t)
  {
    t_ = &t;
    p_ = 0;
    vec_ = BPvecPtr(new BPvec(t_->table_size(), false));
    SCFG::inside_traverse(0, t_->size()-1, *this);
    t_ = NULL;
    return vec_;
  }

  void operator()(uint i, uint j)
  {
    if (i!=j && (*t_)(i,j)>0) (*vec_)[p_]=true;
    p_++;
  }

private:
  const SCFG::CYKTable<uint>* t_;
  BPvecPtr vec_;
  uint p_;
};

template < class L >
struct CountBP
{
  CountBP(const L& bp, uint sz)
    : bp_(bp), table(sz, 0), p(0)
  {}

  void operator()(uint i, uint j)
  {
    typename L::const_iterator x;
    uint c=0;
    for (x=bp_.begin(); x!=bp_.end(); ++x) {
      if (i!=j && (**x)[p]) c++;
    }
    table.update(i, j, static_cast<double>(c)/bp_.size());
    p++;
  }
  
  const L& bp_;
  CentroidFold::BPTable table;
  uint p;
};


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
pfl_fold(T& bp, const std::string& seq, uint max_dist, float th=1e-5)
{
  assert(seq.size()>max_dist);
  bp.resize(seq.size(), max_dist);
  Vienna::pf_scale = -1;
  Vienna::plist *pl;
#ifdef HAVE_VIENNA18
  pl = Vienna::pfl_fold(const_cast<char*>(seq.c_str()), max_dist, max_dist, th, NULL, NULL, NULL, NULL);  
#else
  pl = Vienna::pfl_fold(const_cast<char*>(seq.c_str()), max_dist, max_dist, th, NULL);
#endif
  for (uint k=0; pl[k].i!=0; ++k)
    bp.update(pl[k].i, pl[k].j, pl[k].p);
  free(pl);
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
  strcpy(str2, str.c_str());
#ifdef HAVE_VIENNA18
      Vienna::plist* pi;
#else
      Vienna::pair_info* pi;
#endif
  Vienna::alipf_fold(seqs, str2, &pi);
  for (uint k=0; pi[k].i!=0; ++k)
    bp.update(pi[k].i-1, pi[k].j-1, pi[k].p);
  free(pi);
  for (uint i=0; seqs[i]!=NULL; ++i) delete[] seqs[i];
  Vienna::free_alipf_arrays();
  Vienna::fold_constrained = bk;
  delete[] seqs;
  delete[] str2;
}
#endif  // HAVE_LIBRNA

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

#ifdef HAVE_LIBRNA
static
uint
pf_fold_st(uint num_samples, const std::string& seq, std::list<BPvecPtr>& bp)
{
  int bk_st_back=Vienna::st_back;
  Vienna::st_back=1;

  Vienna::pf_scale = -1;
  Vienna::init_pf_fold(seq.size());
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  EncodeBP encoder;

  // stochastic sampling
  for (uint n=0; n!=num_samples; ++n) {
    SCFG::CYKTable<uint> bp_pos(seq.size(), 0);
    char *str = Vienna::pbacktrack(const_cast<char*>(seq.c_str()));
    assert(seq.size()==strlen(str));
    std::stack<uint> st;
    for (uint i=0; i!=seq.size(); ++i) {
      if (str[i]=='(') {
	st.push(i);
      } else if (str[i]==')') {
	bp_pos(st.top(), i)++;
	st.pop();
      }
    }
    bp.push_back(encoder.execute(bp_pos));
    free(str);
  }

  Vienna::st_back=bk_st_back;
  Vienna::free_pf_arrays();

  return num_samples;
}

static
uint
pf_fold_st(uint num_samples, const std::string& seq, std::list<BPvecPtr>& bp,
           const std::vector<uint>& idxmap, uint sz)
{
  int bk_st_back=Vienna::st_back;
  Vienna::st_back=1;

  Vienna::pf_scale = -1;
  Vienna::init_pf_fold(seq.size());
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  EncodeBP encoder;

  // stochastic sampling
  for (uint n=0; n!=num_samples; ++n) {
    SCFG::CYKTable<uint> bp_pos(sz, 0);
    char *str = Vienna::pbacktrack(const_cast<char*>(seq.c_str()));
    assert(seq.size()==strlen(str));
    std::stack<uint> st;
    for (uint i=0; i!=seq.size(); ++i) {
      if (str[i]=='(') {
	st.push(i);
      } else if (str[i]==')') {
	bp_pos(idxmap[st.top()], idxmap[i])++;
	st.pop();
      }
    }
    bp.push_back(encoder.execute(bp_pos));
    free(str);
  }

  Vienna::st_back=bk_st_back;
  Vienna::free_pf_arrays();

  return num_samples;
}

static
uint
pf_fold_st(uint num_samples, const std::string& seq, std::ostream& out)
{
  int bk_st_back=Vienna::st_back;
  Vienna::st_back=1;

  Vienna::pf_scale = -1;
  Vienna::init_pf_fold(seq.size());
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);

  // stochastic sampling
  for (uint n=0; n!=num_samples; ++n) {
    char *str = Vienna::pbacktrack(const_cast<char*>(seq.c_str()));
    assert(seq.size()==strlen(str));
    out << str << std::endl;
    free(str);
  }

  Vienna::st_back=bk_st_back;
  Vienna::free_pf_arrays();

  return num_samples;
}

static
uint
pf_fold_st(uint num_samples, const std::string& seq, std::ostream& out,
           const std::vector<uint>& idxmap, uint sz)
{
  int bk_st_back=Vienna::st_back;
  Vienna::st_back=1;

  Vienna::pf_scale = -1;
  Vienna::init_pf_fold(seq.size());
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);

  // stochastic sampling
  for (uint n=0; n!=num_samples; ++n) {
    SCFG::CYKTable<uint> bp_pos(sz, 0);
    std::string p(sz, '.');
    char *str = Vienna::pbacktrack(const_cast<char*>(seq.c_str()));
    assert(seq.size()==strlen(str));
    std::stack<uint> st;
    for (uint i=0; i!=seq.size(); ++i) {
      if (str[i]=='(') {
	st.push(i);
      } else if (str[i]==')') {
	p[idxmap[st.top()]] = '(';
        p[idxmap[i]] = ')';
	st.pop();
      }
    }
    out << p << std::endl;
    free(str);
  }

  Vienna::st_back=bk_st_back;
  Vienna::free_pf_arrays();

  return num_samples;
}

#ifdef HAVE_VIENNA18
static
uint
alipf_fold_st(uint num_samples, const std::vector<std::string>& ma, std::list<BPvecPtr>& bp)
{
  EncodeBP encoder;
  // prepare an alignment
  uint length = ma.front().size();
  char** seqs = new char*[ma.size()+1];
  seqs[ma.size()]=NULL;
  std::vector<std::string>::const_iterator x;
  uint i=0;
  for (x=ma.begin(); x!=ma.end(); ++x) {
    assert(x->size()==length);
    seqs[i] = new char[length+1];
    strcpy(seqs[i++], x->c_str());
  }

  char* str2 = new char[length+1];
  //int bk = Vienna::fold_constrained;
  //if (!str.empty()) Vienna::fold_constrained=1;
  // scaling parameters to avoid overflow
  //strcpy(str2, str.c_str());
  double min_en = Vienna::alifold(seqs, str2);
  Vienna::free_alifold_arrays();
  double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(1.07*min_en)/kT/length);
  // build a base pair probablity matrix
  //strcpy(str2, str.c_str());
  Vienna::plist* pi;
  Vienna::alipf_fold(seqs, str2, &pi);

  // stochastic sampling
  for (uint n=0; n!=num_samples; ++n) {
    SCFG::CYKTable<uint> bp_pos(length, 0);
    double prob=1;
    char *str = Vienna::alipbacktrack(&prob);
    //std::cout << s << std::endl << str << std::endl;
    assert(length==strlen(str));
    std::stack<uint> st;
    for (uint i=0; i!=length; ++i) {
      if (str[i]=='(') {
	st.push(i);
      } else if (str[i]==')') {
	bp_pos(st.top(), i)++;
	st.pop();
      }
    }
    bp.push_back(encoder.execute(bp_pos));
    free(str);
  }

  Vienna::free_alipf_arrays();
  free(pi);
  for (uint i=0; seqs[i]!=NULL; ++i) delete[] seqs[i];
  //Vienna::fold_constrained = bk;
  delete[] seqs;
  delete[] str2;

  return num_samples;
}

static
uint
alipf_fold_st(uint num_samples, const std::vector<std::string>& ma, std::ostream& out)
{
  // prepare an alignment
  uint length = ma.front().size();
  char** seqs = new char*[ma.size()+1];
  seqs[ma.size()]=NULL;
  std::vector<std::string>::const_iterator x;
  uint i=0;
  for (x=ma.begin(); x!=ma.end(); ++x) {
    assert(x->size()==length);
    seqs[i] = new char[length+1];
    strcpy(seqs[i++], x->c_str());
  }

  char* str2 = new char[length+1];
  //int bk = Vienna::fold_constrained;
  //if (!str.empty()) Vienna::fold_constrained=1;
  // scaling parameters to avoid overflow
  //strcpy(str2, str.c_str());
  double min_en = Vienna::alifold(seqs, str2);
  Vienna::free_alifold_arrays();
  double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(1.07*min_en)/kT/length);
  // build a base pair probablity matrix
  //strcpy(str2, str.c_str());
  Vienna::plist* pi;
  Vienna::alipf_fold(seqs, str2, &pi);

  // stochastic sampling
  for (uint n=0; n!=num_samples; ++n) {
    SCFG::CYKTable<uint> bp_pos(length, 0);
    double prob=1;
    char *str = Vienna::alipbacktrack(&prob);
    //std::cout << s << std::endl << str << std::endl;
    assert(length==strlen(str));
    out << str << std::endl;
    free(str);
  }

  Vienna::free_alipf_arrays();
  free(pi);
  for (uint i=0; seqs[i]!=NULL; ++i) delete[] seqs[i];
  //Vienna::fold_constrained = bk;
  delete[] seqs;
  delete[] str2;

  return num_samples;
}
#endif  // HAVE_VIENNA18
#endif  // HAVE_LIBRNA

template < class U >
static
uint
contra_fold_st(uint num_samples, CONTRAfold<U>& cf, const std::string& seq, std::list<BPvecPtr>& bp)
{
  EncodeBP encoder;
  cf.PrepareStochasticTraceback(seq);
  for (uint n=0; n!=num_samples; ++n) {
    SCFG::CYKTable<uint> bp_pos(seq.size(), 0);
    std::vector<int> paren = cf.StochasticTraceback();
    for (uint i=0; i!=paren.size(); ++i) {
      if (paren[i]>int(i)) {
        bp_pos(i-1, paren[i]-1)++;
      }
    }
    bp.push_back(encoder.execute(bp_pos));
  }

  return num_samples;
}

template < class U >
static
uint
contra_fold_st(uint num_samples, CONTRAfoldM<U>& cf, const std::vector<std::string>& aln, std::list<BPvecPtr>& bp)
{
  EncodeBP encoder;
  cf.PrepareStochasticTraceback(aln);
  for (uint n=0; n!=num_samples; ++n) {
    SCFG::CYKTable<uint> bp_pos(aln.front().size(), 0);
    std::vector<int> paren = cf.StochasticTraceback();
    for (uint i=0; i!=paren.size(); ++i) {
      if (paren[i]>int(i)) {
        bp_pos(i-1, paren[i]-1)++;
      }
    }
    bp.push_back(encoder.execute(bp_pos));
  }

  return num_samples;
}

template < class U >
static
uint
contra_fold_st(uint num_samples, CONTRAfold<U>& cf, const std::string& seq, std::ostream& out)
{
  cf.PrepareStochasticTraceback(seq);
  for (uint n=0; n!=num_samples; ++n) {
    std::string p(seq.size(), '.');
    std::vector<int> paren = cf.StochasticTraceback();
    for (uint i=0; i!=paren.size(); ++i) {
      if (paren[i]>int(i)) {
        p[i-1] = '(';
        p[paren[i]-1] = ')';
      }
    }
    out << p << std::endl;
  }

  return num_samples;
}

template < class U >
static
uint
contra_fold_st(uint num_samples, CONTRAfoldM<U>& cf, const std::vector<std::string>& aln, std::ostream& out)
{
  EncodeBP encoder;
  cf.PrepareStochasticTraceback(aln);
  for (uint n=0; n!=num_samples; ++n) {
    std::string p(aln.front().size(), '.');
    std::vector<int> paren = cf.StochasticTraceback();
    for (uint i=0; i!=paren.size(); ++i) {
      if (paren[i]>int(i)) {
        p[i-1] = '(';
        p[paren[i]-1] = ')';
      }
    }
    out << p << std::endl;
  }

  return num_samples;
}

CentroidFold::
CentroidFold(uint engine,
             bool run_as_mea,
             uint reserved_size,
             uint seed)
  : engine_(engine),
    mea_(run_as_mea),
    bp_(reserved_size),
    canonical_only_(true),
    contrafold_(NULL),
    model_(),
    max_bp_dist_(0),
    seed_(seed)
{
  if (seed_==0)
  {
#if 0
    time_t t;
    time(&t);
    seed_ = (uint)t;
#else
    timeval t;
    gettimeofday(&t, NULL);
    seed_ = t.tv_usec * t.tv_sec;
#endif
  }
#ifdef HAVE_LIBRNA
  using namespace Vienna;
  // copy from ViennaRNA-x.x.x/lib/utils.c
  xsubi[0] = xsubi[1] = xsubi[2] = (unsigned short) seed_;  /* lower 16 bit */
  xsubi[1] += (unsigned short) ((unsigned)seed_ >> 6);
  xsubi[2] += (unsigned short) ((unsigned)seed_ >> 12);
#endif  // HAVE_LIBRNA
  srand((uint)seed_);
}

CentroidFold::
~CentroidFold()
{
  if (contrafold_) delete contrafold_;
}


#ifdef HAVE_LIBRNA
void
CentroidFold::
set_options_for_pf_fold(bool canonical_only, uint max_dist)
{
  canonical_only_ = canonical_only;
  if (!canonical_only_)
    Vienna::nonstandards = const_cast<char*>("AAACAGCACCCUGAGGUCUU");
  max_bp_dist_ = max_dist;
}
#endif

void
CentroidFold::
set_options_for_contrafold(const std::string& model, bool canonical_only, uint max_bp_dist)
{
  model_ = model;
  canonical_only_ = canonical_only;
  max_bp_dist_ = max_bp_dist;
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
    if (max_bp_dist_==0 || max_bp_dist_>=seq2.size())
      pf_fold(bp_, seq2);
    else
      pfl_fold(bp_, seq2, max_bp_dist_);
    break;
#endif
  case CONTRAFOLD:
    if (contrafold_==NULL) {
      contrafold_ = new CONTRAfold<float>(canonical_only_, max_bp_dist_);
      if (!model_.empty()) contrafold_->SetParameters(model_);
    }
    contra_fold(bp_, *contrafold_, seq2);
    break;
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
  case PFFOLD: {
    std::string str2(str);
    std::replace(str2.begin(), str2.end(), '.', 'x');
    std::replace(str2.begin(), str2.end(), '?', '.');
    if (max_bp_dist_==0 || max_bp_dist_>seq2.size())
      pf_fold(bp_, seq2, str);
    else
      pfl_fold(bp_, seq2, max_bp_dist_);
    break;
  }
#endif
  case CONTRAFOLD:
    if (contrafold_==NULL) {
      contrafold_ = new CONTRAfold<float>(canonical_only_, max_bp_dist_);
      if (!model_.empty()) contrafold_->SetParameters(model_);
    }
    contra_fold(bp_, *contrafold_, seq2, str);
    break;
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
    //std::list<std::string>::const_iterator x;
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
      case CONTRAFOLD:
        if (contrafold_==NULL) {
          contrafold_ = new CONTRAfold<float>(canonical_only_, max_bp_dist_);
          if (!model_.empty()) contrafold_->SetParameters(model_);
        }
        contra_fold(*bpi, *contrafold_, *x);
	break;
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
calculate_posterior(const std::vector<std::string>& seq)
{
  std::list<std::string> seq2(seq.size());
  std::copy(seq.begin(), seq.end(), seq2.begin());
  calculate_posterior(seq2);
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
    std::string str2(str);
    std::replace(str2.begin(), str2.end(), '.', 'x');
    std::replace(str2.begin(), str2.end(), '?', '.');
    alipf_fold(bp_, seq2, str2);
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

    //std::list<std::string>::const_iterator x;
    std::list<std::vector<uint> >::const_iterator y;
    for (x=seqs.begin(), y=idxmaps.begin(); x!=seqs.end(); ++x, ++y) {
      // map the consensus structure for constraints
      BPTablePtr bpi(new BPTable);
      std::string str2(x->size(), ' ');
      for (uint i=0, j=0; i!=y->size(); ++i) {
        if ((*y)[i]!=static_cast<uint>(-1)) {
          if ((*y)[bpmap[i]]!=static_cast<uint>(-1))
            str2[j++] = str[i];
          else
            str2[j++] = '?';
        }
      }
      
      switch (engine_) {
      case AUX:
        assert(!"AUX should be given bp matrices.");
	break;
#ifdef HAVE_LIBRNA
      case PFFOLD: {
        assert(max_bp_dist_==0);
        std::string str3(str);
        std::replace(str3.begin(), str3.end(), '.', 'x');
        std::replace(str3.begin(), str3.end(), '?', '.');
	pf_fold(*bpi, *x, str3);
	break;
      }
#endif
      case CONTRAFOLD:
        if (contrafold_==NULL) {
          contrafold_ = new CONTRAfold<float>(canonical_only_, max_bp_dist_);
          if (!model_.empty()) contrafold_->SetParameters(model_);
        }
        contra_fold(*bpi, *contrafold_, *x, str2);
	break;
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
calculate_posterior(const std::vector<std::string>& seq, const std::string& str)
{
  std::list<std::string> seq2(seq.size());
  std::copy(seq.begin(), seq.end(), seq2.begin());
  calculate_posterior(seq2, str);
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
      const std::vector<float>& gamma) const
{
  out << ">" << name << std::endl << seq << std::endl;
  std::vector<float>::const_iterator g;
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

std::string
CentroidFold::
posterior(const std::string& seq, float th) const
{
  std::ostringstream out;
  bp_.save(out, seq, th);
  return out.str();
}


static
void
stochastic_fold_helper(const std::string& name, const std::string& seq,
                       const std::vector<BPvecPtr>& bpv, uint max_clusters,
                       const std::vector<float>& gamma, std::ostream& out,
                       const std::string& p_outname, float th)
{    
  if (max_clusters>=2) {
    // calculate the distance matrix
    typedef boost::multi_array<double, 2> DMatrix;
    DMatrix dmatrix(boost::extents[bpv.size()][bpv.size()]);
    for (uint i=0; i!=bpv.size(); ++i) {
      for (uint j=i; j!=bpv.size(); ++j) {
	dmatrix[i][j]=dmatrix[j][i]=hamming_distance(*bpv[i], *bpv[j]);
      }
    }

    // hierarchical clustering by DIANA
    HCLUST::Diana<DMatrix> diana(dmatrix);
    diana.build(max_clusters);
    uint opt_n = diana.optimal_size();
    std::vector<uint> res;
    std::vector<uint> num;
    diana.get_clusters(opt_n, res, num);

    // centroid estimation for each cluster
    uint p=0;
    for (uint i=0; i!=num.size(); ++i) {
      std::vector<BPvecPtr> v(num[i]);
      for (uint j=0; j!=v.size(); ++j) v[j] = bpv[res[p++]];
      CountBP< std::vector<BPvecPtr> > count_bp(v, seq.size());
      SCFG::inside_traverse(0, seq.size()-1, count_bp);
      char buf[100];
      sprintf(buf, "%3.1f%%", static_cast<float>(num[i])/bpv.size()*100);
      out << ">" << name << " (" << i+1 << " of " << num.size() << ", size="
          << buf << ")" << std::endl
          << seq << std::endl;

      std::vector<float>::const_iterator g;
      for (g=gamma.begin(); g!=gamma.end(); ++g) {
        std::string paren(seq.size(), '.');
        SCFG::Centroid::execute(count_bp.table, paren, *g);
        out << paren << " (g=" << *g << ",th=" << (1.0/(1.0+*g)) << ")" << std::endl;
      }

      if (!p_outname.empty())
      {
        char buf[1024];
        snprintf(buf, sizeof(buf), "%s-%d", p_outname.c_str(), i);
        count_bp.table.save(buf, seq, th);
      }
    }
  } else {
    CountBP< std::vector<BPvecPtr> > count_bp(bpv, seq.size());
    SCFG::inside_traverse(0, seq.size()-1, count_bp);
    std::string paren(seq.size(), '.');
    std::vector<float>::const_iterator g;
    std::cout << ">" << name << std::endl
	      << seq << std::endl;
    for (g=gamma.begin(); g!=gamma.end(); ++g) {
      std::fill(paren.begin(), paren.end(), '.');
      SCFG::Centroid::execute(count_bp.table, paren, *g);
      std::cout << paren << " (g=" << *g << ",th=" << (1.0/(1.0+*g)) << ")" << std::endl;
    }
  }
}

void
CentroidFold::
stochastic_fold(const std::string& name, const std::string& seq,
                uint num_samples, uint max_clusters,
                const std::vector<float>& gamma, std::ostream& out,
                const std::string& p_outname, float th)
{
  std::list<BPvecPtr> bpvl;

  std::string seq2(seq);
  std::replace(seq2.begin(), seq2.end(), 't', 'u');
  std::replace(seq2.begin(), seq2.end(), 'T', 'U');
  switch (engine_) {
#ifdef HAVE_LIBRNA
  case PFFOLD:
    if (max_clusters>0)
      pf_fold_st(num_samples, seq2, bpvl);
    else
      pf_fold_st(num_samples, seq2, out);
    break;
#endif
  case CONTRAFOLD:
    if (contrafold_==NULL) {
      contrafold_ = new CONTRAfold<float>(canonical_only_, max_bp_dist_);
      if (!model_.empty()) contrafold_->SetParameters(model_);
      contrafold_->init_rand(seed_);
    }
    if (max_clusters>0)
      contra_fold_st(num_samples, *contrafold_, seq2, bpvl);
    else
      contra_fold_st(num_samples, *contrafold_, seq2, out);
    break;
  default:
    throw "not supported yet";
    break;
  }

  if (max_clusters>0) {
    std::vector<BPvecPtr> bpv(bpvl.size());
    std::copy(bpvl.begin(), bpvl.end(), bpv.begin());
    bpvl.clear();
    stochastic_fold_helper(name, seq, bpv, max_clusters, gamma, out, p_outname, th);
  }
}

void
CentroidFold::
stochastic_fold(const std::string& name, const std::string& consensus,
                const std::vector<std::string>& seq,
                uint num_samples, uint max_clusters,
                const std::vector<float>& gamma, std::ostream& out,
                const std::string& p_outname, float th)
{
  std::list<std::string> seq2(seq.size());
  std::copy(seq.begin(), seq.end(), seq2.begin());
  stochastic_fold(name, consensus, seq2, num_samples, max_clusters, gamma, out, p_outname, th);
}

void
CentroidFold::
stochastic_fold(const std::string& name, const std::string& consensus,
                const std::list<std::string>& seq,
                uint num_samples, uint max_clusters,
                const std::vector<float>& gamma, std::ostream& out,
                const std::string& p_outname, float th)
{
  std::list<BPvecPtr> bpvl;

  switch (engine_)
  {
#ifdef HAVE_LIBRNA
    case PFFOLD:
    {
      uint sz = seq.front().size();
      std::list<std::string>::const_iterator x;
      for (x=seq.begin(); x!=seq.end(); ++x)
      {
        std::string s;
        std::vector<uint> idxmap(x->size(), static_cast<uint>(-1));
        uint j=0;
        for (uint i=0; i!=x->size(); ++i)
        {
          if ((*x)[i]!='-')
          {
            switch ((*x)[i])
            {
              case 't': s += 'u'; break;
              case 'T': s += 'U'; break;
              default:  s += (*x)[i];
            }
            idxmap[j++] = i;
          }
        }
        idxmap.resize(j);
        if (max_clusters>0)
          pf_fold_st(num_samples, s, bpvl, idxmap, sz);
        else
          pf_fold_st(num_samples, s, out, idxmap, sz);
      }
      break;
    }
#ifdef HAVE_VIENNA18
    case ALIPFFOLD:
    {
      std::vector<std::string> aln(seq.size());
      std::copy(seq.begin(), seq.end(), aln.begin());
      for (uint i=0; i!=aln.size(); ++i)
      {
        std::replace(aln[i].begin(), aln[i].end(), 't', 'u');
        std::replace(aln[i].begin(), aln[i].end(), 'T', 'U');
      }
      if (max_clusters>0)
        alipf_fold_st(num_samples, aln, bpvl);
      else
        alipf_fold_st(num_samples, aln, out);
      break;
    }
#endif
#endif
    case CONTRAFOLD:
    {
      CONTRAfoldM<float>* cf = new CONTRAfoldM<float>(canonical_only_, max_bp_dist_);
      if (!model_.empty()) cf->SetParameters(model_);
      cf->init_rand(seed_);
      std::vector<std::string> aln(seq.size());
      std::copy(seq.begin(), seq.end(), aln.begin());
      for (uint i=0; i!=aln.size(); ++i)
      {
        std::replace(aln[i].begin(), aln[i].end(), 't', 'u');
        std::replace(aln[i].begin(), aln[i].end(), 'T', 'U');
      }
      if (max_clusters>0)
        contra_fold_st(num_samples, *cf, aln, bpvl);
      else
        contra_fold_st(num_samples, *cf, aln, out);
      delete cf;
      break;
    }

    default:
      throw "not supported yet";
      break;
  }

  if (max_clusters>0) {
    std::vector<BPvecPtr> bpv(bpvl.size());
    std::copy(bpvl.begin(), bpvl.end(), bpv.begin());
    bpvl.clear();
    stochastic_fold_helper(name, consensus, bpv, max_clusters, gamma, out, p_outname, th);
  }
}

void
CentroidFold::
ps_plot(const std::string& name, const std::string& seq, float g, bool color) const
{
  std::string paren;
  decode_structure(g, paren);
  if (color)
    ::ps_color_plot(seq.c_str(), paren.c_str(), bp_, name.c_str());
  else
    ::ps_plot(seq.c_str(), paren.c_str(), name.c_str());
}

#ifdef HAVE_LIBRNA
void
CentroidFold::
svg_plot(const std::string& name, const std::string& seq, float g) const
{
  std::string paren;
  decode_structure(g, paren);
  Vienna::svg_rna_plot(const_cast<char*>(seq.c_str()),
		       const_cast<char*>(paren.c_str()),
                       const_cast<char*>(name.c_str()));
}
#endif
