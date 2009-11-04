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
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/alifold.h>
#include <ViennaRNA/LPfold.h>
#include <ViennaRNA/PS_dot.h>
#include <ViennaRNA/aln_util.h>
#include <ViennaRNA/utils.h>
  extern int eos_debug;
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
contra_fold_st(uint num_samples, CONTRAfold<U>& cf, const std::string& seq, std::list<BPvecPtr>& bp,
               const std::vector<uint>& idxmap, uint sz)
{
  EncodeBP encoder;
  cf.PrepareStochasticTraceback(seq);
  for (uint n=0; n!=num_samples; ++n) {
    SCFG::CYKTable<uint> bp_pos(sz, 0);
    std::vector<int> paren = cf.StochasticTraceback();
    for (uint i=0; i!=paren.size(); ++i) {
      if (paren[i]>int(i)) {
        bp_pos(idxmap[i-1], idxmap[paren[i]-1])++;
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

template < class U >
static
uint
contra_fold_st(uint num_samples, CONTRAfold<U>& cf, const std::string& seq, std::ostream& out,
               const std::vector<uint>& idxmap, uint sz)
{
  EncodeBP encoder;
  cf.PrepareStochasticTraceback(seq);
  for (uint n=0; n!=num_samples; ++n) {
    std::string p(sz, '.');
    std::vector<int> paren = cf.StochasticTraceback();
    for (uint i=0; i!=paren.size(); ++i) {
      if (paren[i]>int(i)) {
        p[idxmap[i-1]] = '(';
        p[idxmap[paren[i]-1]] = ')';
      }
    }
    out << p << std::endl;
  }

  return num_samples;
}

CentroidFold::
CentroidFold(uint engine,
             bool run_as_mea,
	     int num_ea_samples,
             uint reserved_size,
             uint seed)
  : engine_(engine),
    mea_(run_as_mea),
    num_ea_samples_(num_ea_samples),
    bp_(reserved_size),
    canonical_only_(true),
    contrafold_(NULL),
    model_(),
    max_bp_dist_(0),
    seed_(seed),
    dist_type_(0)
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
  srand((uint)seed_);

  switch (engine_)
  {
#ifdef HAVE_LIBRNA
    case PFFOLD:
    case ALIPFFOLD:
    {
      using namespace Vienna;
      // copy from ViennaRNA-x.x.x/lib/utils.c
      xsubi[0] = xsubi[1] = xsubi[2] = (unsigned short) seed_;  /* lower 16 bit */
      xsubi[1] += (unsigned short) ((unsigned)seed_ >> 6);
      xsubi[2] += (unsigned short) ((unsigned)seed_ >> 12);
      break;
    }
#endif  // HAVE_LIBRNA
    case CONTRAFOLD:
      contrafold_ = new CONTRAfold<float>(canonical_only_, max_bp_dist_);
      contrafold_->init_rand(seed_);
      break;
    default:
      break;
  }
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
set_options_for_contrafold(const std::string& model, bool canonical_only, uint max_bp_dist, uint dist_type/*=0*/)
{
  model_ = model;
  canonical_only_ = canonical_only;
  max_bp_dist_ = max_bp_dist;
  dist_type_ = dist_type;
  if (!model_.empty()) contrafold_->SetParameters(model_);
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

struct Result
{
  Result(float g_, float sen_, float ppv_, float mcc_, float ea_, float e_,
         const std::string& paren_)
    : g(g_), sen(sen_), ppv(ppv_), mcc(mcc_), ea(ea_), e(e_), paren(paren_)
  { }

  float g;
  float sen;
  float ppv;
  float mcc;
  float ea;
  float e;
  std::string paren;
};

struct cmp_by_mcc
{
  bool operator()(const Result& x, const Result& y) const
  {
    return x.paren!=y.paren ? x.mcc>y.mcc : x.g<y.g;
  }
};

void
CentroidFold::
print(std::ostream& out, const std::string& name, const std::string& seq,
      const std::vector<float>& gamma) const
{
#ifdef HAVE_LIBRNA
  Vienna::eos_debug = -1;
#endif
  if (num_ea_samples_>=0)
  {
    std::vector<Result> res;
    std::vector<float>::const_iterator g;
    for (g=gamma.begin(); g!=gamma.end(); ++g) {
      std::string paren;
      float p = decode_structure(*g, paren);
      double sen = 0.0, ppv = 0.0, mcc = 0.0;
      if (num_ea_samples_==0) 
        compute_expected_accuracy (paren, bp_, sen, ppv, mcc);
      else if (num_ea_samples_>0) 
        compute_expected_accuracy_sampling (paren, seq, num_ea_samples_, sen, ppv, mcc);
      float e=0;
#ifdef HAVE_LIBRNA
      e = Vienna::energy_of_struct(seq.c_str(), paren.c_str());
#endif
      res.push_back(Result(*g, sen, ppv, mcc, p, e, paren));
    }
    std::sort(res.begin(), res.end(), cmp_by_mcc());

    std::string prev_paren;
    out << ">" << name << std::endl << seq << std::endl;
    std::vector<Result>::const_iterator r;
    for (r=res.begin(); r!=res.end(); ++r) {
      if (prev_paren!=r->paren) {
        out << r->paren;
        if (!mea_) {
          out << " (g=" << r->g << ",th=" << (1.0/(1.0+r->g));
#ifdef HAVE_LIBRNA
          out << ",e=" << r->e;
#endif
          out << ")";
        } else {
          out << " (g=" << r->g << ",EA=" << r->ea;
#ifdef HAVE_LIBRNA
          out << ",e=" << r->e;
#endif
          out << ")";
        }
        if (num_ea_samples_>=0) {
          out << " [" << "SEN=" << r->sen << ",PPV=" << r->ppv <<  ",MCC=" << r->mcc << "]";
        }
        out << std::endl;
        prev_paren = r->paren;
      }
    }
  }
  else
  {
    out << ">" << name << std::endl << seq << std::endl;
    std::vector<float>::const_iterator g;
    for (g=gamma.begin(); g!=gamma.end(); ++g) {
      std::string paren;
      float p = decode_structure(*g, paren);
      out << paren;
#ifdef HAVE_LIBRNA
      float e = Vienna::energy_of_struct(seq.c_str(), paren.c_str());
#endif
      if (!mea_) {
        out << " (g=" << *g << ",th=" << (1.0/(1.0+*g));
#ifdef HAVE_LIBRNA
        out << ",e=" << e;
#endif
        out  << ")";
      } else {
        out << " (g=" << *g << ",EA=" << p;
#ifdef HAVE_LIBRNA
        out << ",e=" << e;
#endif
        out << ")";
      }
      out << std::endl;
    }
  }
}

void
CentroidFold::
print(std::ostream& out, const Aln& aln, const std::vector<float>& gamma) const
{
#ifdef HAVE_LIBRNA
  Vienna::eos_debug = -1;
#endif
  if (num_ea_samples_>=0)
  {
    std::vector<Result> res;
    std::vector<float>::const_iterator g;
    for (g=gamma.begin(); g!=gamma.end(); ++g) {
      std::string paren;
      float p = decode_structure(*g, paren);
      double sen = 0.0, ppv = 0.0, mcc = 0.0;
      if (num_ea_samples_==0) 
        compute_expected_accuracy (paren, bp_, sen, ppv, mcc);
      else if (num_ea_samples_>0)
        compute_expected_accuracy_sampling (paren, aln, num_ea_samples_, sen, ppv, mcc);
      float e=0;
#ifdef HAVE_LIBRNA
      e = aln.energy_of_struct(paren);
#endif
      res.push_back(Result(*g, sen, ppv, mcc, p, e, paren));
    }
    std::sort(res.begin(), res.end(), cmp_by_mcc());

    std::string prev_paren;
    out << ">" << aln.name().front() << std::endl << aln.consensus() << std::endl;
    std::vector<Result>::const_iterator r;
    for (r=res.begin(); r!=res.end(); ++r) {
      if (prev_paren!=r->paren) {
        out << r->paren;
        if (!mea_) {
          out << " (g=" << r->g << ",th=" << (1.0/(1.0+r->g));
#ifdef HAVE_LIBRNA
          out << ",e=" << r->e;
#endif
          out << ")";
        } else {
          out << " (g=" << r->g << ",EA=" << r->ea;
#ifdef HAVE_LIBRNA
          out << ",e=" << r->e;
#endif
          out << ")";
        }
        if (num_ea_samples_>=0) {
          out << " [" << "SEN=" << r->sen << ",PPV=" << r->ppv <<  ",MCC=" << r->mcc << "]";
        }
        out << std::endl;
        prev_paren = r->paren;
      }
    }
  }
  else
  {
    out << ">" << aln.name().front() << std::endl << aln.consensus() << std::endl;
    std::vector<float>::const_iterator g;
    for (g=gamma.begin(); g!=gamma.end(); ++g) {
      std::string paren;
      float p = decode_structure(*g, paren);
      out << paren;
#ifdef HAVE_LIBRNA
      float e = aln.energy_of_struct(paren);
#endif
      if (!mea_) {
        out << " (g=" << *g << ",th=" << (1.0/(1.0+*g));
#ifdef HAVE_LIBRNA
        out << ",e=" << e;
#endif
        out  << ")";
      } else {
        out << " (g=" << *g << ",EA=" << p;
#ifdef HAVE_LIBRNA
        out << ",e=" << e;
#endif
        out << ")";
      }
      out << std::endl;
    }
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

struct cmp_by_size
{
  cmp_by_size(const std::vector<uint>& num) : num_(num) { }

  bool operator()(uint i, uint j) const
  {
    return num_[i]>num_[j];
  }

  const std::vector<uint>& num_;  
};

const std::string& seq(const std::string& t) { return t; }
uint length(const std::string& t) { return t.size(); }
#ifdef HAVE_LIBRNA
float energy(const std::string& t, const std::string& s)
{
  Vienna::eos_debug = -1;
  return Vienna::energy_of_struct(t.c_str(), s.c_str());
}
#endif

std::string seq(const Aln& aln) { return aln.consensus(); }
uint length(const Aln& aln) { return aln.seq().front().size(); }
#ifdef HAVE_LIBRNA
float energy(const Aln& aln, const std::string& s)
{
  Vienna::eos_debug = -1;
  return aln.energy_of_struct(s);
}
#endif

template < class T >
static
void
stochastic_fold_helper(const std::string& name, const T& t,
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

    // sort by size of clusters
    std::vector<uint> s(num.size(), 0);
    std::vector<uint> idx(num.size());
    for (uint i=0; i!=num.size(); ++i) {
      idx[i]=i;
      if (i!=0) s[i]=s[i-1]+num[i-1];
    }
    std::sort(idx.begin(), idx.end(), cmp_by_size(num));

    // centroid estimation for each cluster
    for (uint i=0; i!=idx.size(); ++i) {
      std::vector<BPvecPtr> v(num[idx[i]]);
      for (uint j=0; j!=v.size(); ++j) v[j] = bpv[res[s[idx[i]]+j]];
      CountBP< std::vector<BPvecPtr> > count_bp(v, length(t));
      SCFG::inside_traverse(0, length(t)-1, count_bp);
      char buf[100];
      sprintf(buf, "%3.1f%%", static_cast<float>(num[idx[i]])/bpv.size()*100);
      out << ">" << name << " (" << i+1 << " of " << num.size() << ", size="
          << buf << ")" << std::endl
          << seq(t) << std::endl;

      std::vector<float>::const_iterator g;
      for (g=gamma.begin(); g!=gamma.end(); ++g) {
        std::string paren(length(t), '.');
        SCFG::Centroid::execute(count_bp.table, paren, *g);
        out << paren << " (g=" << *g << ",th=" << (1.0/(1.0+*g));
#ifdef HAVE_LIBRNA
        out << ",e=" << energy(t, paren);
#endif
        out << ")" << std::endl;
      }

      if (!p_outname.empty())
      {
        char buf[1024];
        snprintf(buf, sizeof(buf), "%s-%d", p_outname.c_str(), i);
        count_bp.table.save(buf, seq(t), th);
      }
    }
  } else {
    CountBP< std::vector<BPvecPtr> > count_bp(bpv, length(t));
    SCFG::inside_traverse(0, length(t)-1, count_bp);
    std::string paren(length(t), '.');
    std::vector<float>::const_iterator g;
    out << ">" << name << std::endl
        << seq(t) << std::endl;
    for (g=gamma.begin(); g!=gamma.end(); ++g) {
      std::fill(paren.begin(), paren.end(), '.');
      SCFG::Centroid::execute(count_bp.table, paren, *g);
      out << paren << " (g=" << *g << ",th=" << (1.0/(1.0+*g));
#ifdef HAVE_LIBRNA
        out << ",e=" << energy(t, paren);
#endif
      out << ")" << std::endl;
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
stochastic_fold(const Aln& aln,
                uint num_samples, uint max_clusters,
                const std::vector<float>& gamma, std::ostream& out,
                const std::string& p_outname, float th)
{
  std::list<BPvecPtr> bpvl;
  const std::list<std::string>& seq = aln.seq();

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
      if (dist_type_==0)
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
      }
      else
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
            contra_fold_st(num_samples, *contrafold_, s, bpvl, idxmap, sz);
          else
            contra_fold_st(num_samples, *contrafold_, s, out, idxmap, sz);
        }
      }
      break;

    default:
      throw "not supported yet";
      break;
  }

  if (max_clusters>0) {
    std::vector<BPvecPtr> bpv(bpvl.size());
    std::copy(bpvl.begin(), bpvl.end(), bpv.begin());
    bpvl.clear();
    stochastic_fold_helper(aln.name().front(), aln, bpv, max_clusters, gamma, out, p_outname, th);
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

// added by M. Hamada

inline
void
get_TP_TN_FP_FN (const BPvec& ref, const BPvec& pre, 
		 double& TP, double& TN, double& FP, double& FN)
{
  TP = (ref&pre).count();
  FP = pre.count() - TP;
  FN = ref.count() - TP;
  TN = ref.size() - ref.count() - FP; 
}

// Added by M. Hamada (2009/09/18)
void CentroidFold::compute_expected_accuracy (const std::string& paren,
					      const BPTable& bp, 
					      double& sen, double& ppv, double& mcc)
{
  std::vector< std::pair<int, int> > BP;
  std::vector<int> st;
  for (uint n=0; n<paren.size(); ++n) {
    const char c = paren[n];
    if (c == '(') {
      st.push_back (n);
    }
    else if (c == ')') {
      int left = st.back();
      st.pop_back();
      BP.push_back ( std::make_pair<int, int>(left, n) );
    }
  }
  assert (st.size()==0);

  int L  = paren.size ();
  int L2 = L*(L-1)/2;
  int N  = BP.size();

  double sump = 0.0;
  double etp  = 0.0;

  if (bp.size() != paren.size()) 
    std::cerr << bp.size() << "!=" << paren.size() << std::endl << paren << std::endl;; 
  assert (bp.size()==paren.size());

  for (uint i=0; i<bp.size(); ++i) {
    for (uint j = i+1; j<bp.size(); ++j) {
      double p = bp(i,j);
      sump += p;
    }
  }

  for (std::vector<std::pair<int, int> >::const_iterator it=BP.begin();
       it != BP.end(); ++it) {
    etp += bp(it->first, it->second);
  }

  double etn = L2 - N - sump + etp;
  double efp = N - etp;
  double efn = sump - etp;

  sen = ppv = mcc = 0;
  if (etp+efn!=0) sen = etp / (etp + efn);
  if (etp+efp!=0) ppv = etp / (etp + efp);
  if (etp+efp!=0 && etp+efn!=0 && etn+efp!=0 && etn+efn!=0)
    mcc = (etp*etn-efp*efn) / std::sqrt((etp+efp)*(etp+efn)*(etn+efp)*(etn+efn));
}

#ifdef HAVE_LIBRNA
static
void
pf_fold_ea (const BPvecPtr& bp, const std::string& seq, int num_samples,
	    double& sen, double& ppv, double& mcc)
{
  int bk_st_back=Vienna::st_back;
  Vienna::st_back=1;

  Vienna::pf_scale = -1;
  Vienna::init_pf_fold(seq.size());
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  EncodeBP encoder;

  sen = ppv = mcc = 0.0;
  int N = 0;
  // stochastic sampling
  for (int n=0; n!=num_samples; ++n) {
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
    free(str);    
    double TP, TN, FP, FN;
    get_TP_TN_FP_FN (*encoder.execute(bp_pos), *bp, TP, TN, FP, FN);
    if (TP+FN!=0) sen += (double) TP / (TP+FN);
    if (TP+FP!=0) ppv += (double) TP / (TP+FP);
    if (TP+FP!=0 && TP+FN!=0 && TN+FP!=0 && TN+FN!=0)
      mcc += (TP*TN-FP*FN)/std::sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    ++N;
  }

  sen /= N;
  ppv /= N;
  mcc /= N;

  Vienna::st_back=bk_st_back;
  Vienna::free_pf_arrays();
}

static
void
pf_fold_ea(const BPvecPtr& bp, const std::vector<std::string>& ma, uint num_samples,
           double& sen, double& ppv, double& mcc)
{
  int bk_st_back=Vienna::st_back;
  Vienna::st_back=1;

  sen = ppv = mcc = 0.0;
  int N = 0;

  uint sz = ma.front().size();
  std::vector<std::string>::const_iterator x;
  for (x=ma.begin(); x!=ma.end(); ++x)
  {
    std::string seq;
    std::vector<uint> idxmap(x->size(), static_cast<uint>(-1));
    uint j=0;
    for (uint i=0; i!=x->size(); ++i)
    {
      if ((*x)[i]!='-')
      {
        switch ((*x)[i])
        {
          case 't': seq += 'u'; break;
          case 'T': seq += 'U'; break;
          default:  seq += (*x)[i];
        }
        idxmap[j++] = i;
      }
    }

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
      free(str);

      double TP, TN, FP, FN;
      get_TP_TN_FP_FN (*encoder.execute(bp_pos), *bp, TP, TN, FP, FN);
      if (TP+FN!=0) sen += (double) TP / (TP+FN);
      if (TP+FP!=0) ppv += (double) TP / (TP+FP);
      if (TP+FP!=0 && TP+FN!=0 && TN+FP!=0 && TN+FN!=0)
        mcc += (TP*TN-FP*FN)/std::sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
      ++N;
    }

    Vienna::free_pf_arrays();
  }

  sen /= N;
  ppv /= N;
  mcc /= N;

  Vienna::st_back=bk_st_back;
}

#ifdef HAVE_VIENNA18
static
void
alipf_fold_ea(const BPvecPtr& bp, const std::vector<std::string>& ma, uint num_samples,
              double& sen, double& ppv, double& mcc)
{
  EncodeBP encoder;
  sen = ppv = mcc = 0.0;
  int N = 0;

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
    free(str);

    double TP, TN, FP, FN;
    get_TP_TN_FP_FN (*encoder.execute(bp_pos), *bp, TP, TN, FP, FN);
    if (TP+FN!=0) sen += (double) TP / (TP+FN);
    if (TP+FP!=0) ppv += (double) TP / (TP+FP);
    if (TP+FP!=0 && TP+FN!=0 && TN+FP!=0 && TN+FN!=0)
      mcc += (TP*TN-FP*FN)/std::sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    ++N;
  }

  sen /= N;
  ppv /= N;
  mcc /= N;

  Vienna::free_alipf_arrays();
  free(pi);
  for (uint i=0; seqs[i]!=NULL; ++i) delete[] seqs[i];
  //Vienna::fold_constrained = bk;
  delete[] seqs;
  delete[] str2;
}
#endif
#endif

template <class U>
static
void
contra_fold_ea (CONTRAfold<U>& cf, const BPvecPtr& bp, const std::string& seq, int num_samples,
		double& sen, double& ppv, double& mcc)
{
  sen = ppv = mcc = 0.0;
  int N = 0;

  EncodeBP encoder;
  cf.PrepareStochasticTraceback(seq);
  for (int n=0; n!=num_samples; ++n) {
    SCFG::CYKTable<uint> bp_pos(seq.size(), 0);
    std::vector<int> paren = cf.StochasticTraceback();
    for (uint i=0; i!=paren.size(); ++i) {
      if (paren[i]>int(i)) {
        bp_pos(i-1, paren[i]-1)++;
      }
    }
    double TP, TN, FP, FN;
    get_TP_TN_FP_FN (*encoder.execute(bp_pos), *bp, TP, TN, FP, FN);
    if (TP+FN!=0) sen += (double) TP / (TP+FN);
    if (TP+FP!=0) ppv += (double) TP / (TP+FP);
    if (TP+FP!=0 && TP+FN!=0 && TN+FP!=0 && TN+FN!=0)
      mcc += (TP*TN-FP*FN)/std::sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    ++N;
  }  

  sen /= N;
  ppv /= N;
  mcc /= N;
}

template < class U >
static
void
contra_fold_ea(CONTRAfoldM<U>& cf, const BPvecPtr& bp, const std::vector<std::string>& aln, uint num_samples,
               double& sen, double& ppv, double& mcc)
{
  sen = ppv = mcc = 0.0;
  int N = 0;

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
    double TP, TN, FP, FN;
    get_TP_TN_FP_FN (*encoder.execute(bp_pos), *bp, TP, TN, FP, FN);
    if (TP+FN!=0) sen += (double) TP / (TP+FN);
    if (TP+FP!=0) ppv += (double) TP / (TP+FP);
    if (TP+FP!=0 && TP+FN!=0 && TN+FP!=0 && TN+FN!=0)
      mcc += (TP*TN-FP*FN)/std::sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    ++N;
  }

  sen /= N;
  ppv /= N;
  mcc /= N;
}

template < class U >
static
void
contra_fold_ea(CONTRAfold<U>& cf, const BPvecPtr& bp, const std::vector<std::string>& ma, uint num_samples,
               double& sen, double& ppv, double& mcc)
{
  sen = ppv = mcc = 0.0;
  int N = 0;

  uint sz = ma.front().size();
  std::vector<std::string>::const_iterator x;
  for (x=ma.begin(); x!=ma.end(); ++x)
  {
    std::string seq;
    std::vector<uint> idxmap(x->size(), static_cast<uint>(-1));
    uint j=0;
    for (uint i=0; i!=x->size(); ++i)
    {
      if ((*x)[i]!='-')
      {
        switch ((*x)[i])
        {
          case 't': seq += 'u'; break;
          case 'T': seq += 'U'; break;
          default:  seq += (*x)[i];
        }
        idxmap[j++] = i;
      }
    }
    idxmap.resize(j);

    EncodeBP encoder;
    cf.PrepareStochasticTraceback(seq);
    // stochastic sampling
    for (uint n=0; n!=num_samples; ++n) {
      SCFG::CYKTable<uint> bp_pos(sz, 0);
      std::vector<int> paren = cf.StochasticTraceback();
      for (uint i=0; i!=paren.size(); ++i) {
        if (paren[i]>int(i)) {
          bp_pos(idxmap[i-1], idxmap[paren[i]-1])++;
        }
      }
      double TP, TN, FP, FN;
      get_TP_TN_FP_FN (*encoder.execute(bp_pos), *bp, TP, TN, FP, FN);
      if (TP+FN!=0) sen += (double) TP / (TP+FN);
      if (TP+FP!=0) ppv += (double) TP / (TP+FP);
      if (TP+FP!=0 && TP+FN!=0 && TN+FP!=0 && TN+FN!=0)
        mcc += (TP*TN-FP*FN)/std::sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
      ++N;
    }
  }

  sen /= N;
  ppv /= N;
  mcc /= N;
}

void CentroidFold::compute_expected_accuracy_sampling (const std::string& paren,  
						       const std::string& seq, 
						       uint num_samples, 
						       double& sen, double& ppv, double& mcc) const
{
  SCFG::CYKTable<uint> bp_pos(paren.size(), 0);
  std::stack<uint> st;
  for (uint i=0; i<paren.size(); ++i) {
    if (paren[i]=='(') {
      st.push(i);
    } else if (paren[i]==')') {
      bp_pos(st.top(), i)++;
      st.pop();
    }
  }
  EncodeBP encoder;
  BPvecPtr bp = encoder.execute(bp_pos);

  // stochastic sampling
  std::string seq2(seq);
  std::replace(seq2.begin(), seq2.end(), 't', 'u');
  std::replace(seq2.begin(), seq2.end(), 'T', 'U');

  switch (engine_) {
#ifdef HAVE_LIBRNA
  case PFFOLD:
    pf_fold_ea (bp, seq2, num_samples, sen , ppv, mcc);
    break;
#endif
  case CONTRAFOLD:
    contra_fold_ea (*contrafold_, bp, seq2, num_samples, sen , ppv, mcc);
    break;
  default:
    throw "not supported yet";
    break;
  }    
}

void CentroidFold::compute_expected_accuracy_sampling (const std::string& paren,  
						       const Aln& aln, uint num_samples, 
						       double& sen, double& ppv, double& mcc) const
{
  SCFG::CYKTable<uint> bp_pos(paren.size(), 0);
  std::stack<uint> st;
  for (uint i=0; i<paren.size(); ++i) {
    if (paren[i]=='(') {
      st.push(i);
    } else if (paren[i]==')') {
      bp_pos(st.top(), i)++;
      st.pop();
    }
  }
  EncodeBP encoder;
  BPvecPtr bp = encoder.execute(bp_pos);

  // stochastic sampling
  const std::list<std::string>& seq = aln.seq();

  switch (engine_)
  {
#ifdef HAVE_LIBRNA
    case PFFOLD:
    {
      std::vector<std::string> aln(seq.size());
      std::copy(seq.begin(), seq.end(), aln.begin());
      for (uint i=0; i!=aln.size(); ++i)
      {
        std::replace(aln[i].begin(), aln[i].end(), 't', 'u');
        std::replace(aln[i].begin(), aln[i].end(), 'T', 'U');
      }
      pf_fold_ea(bp, aln, num_samples, sen, ppv, mcc);
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
      alipf_fold_ea(bp, aln, num_samples, sen, ppv, mcc);
      break;
    }
#endif
#endif

    case CONTRAFOLD:
      if (dist_type_==0)
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
        contra_fold_ea(*cf, bp, aln, num_samples, sen, ppv, mcc);
        delete cf;
      }
      else
      {
        std::vector<std::string> aln(seq.size());
        std::copy(seq.begin(), seq.end(), aln.begin());
        for (uint i=0; i!=aln.size(); ++i)
        {
          std::replace(aln[i].begin(), aln[i].end(), 't', 'u');
          std::replace(aln[i].begin(), aln[i].end(), 'T', 'U');
        }
        contra_fold_ea(*contrafold_, bp, aln, num_samples, sen, ppv, mcc);
      }
      break;

    default:
      throw "not supported yet";
      break;
  }
}

#ifdef HAVE_LIBRNA
static
uint
pf_fold_max_mcc(uint num_samples, 
		const std::string& name, const std::string& seq, 
		const CentroidFold::BPTable& bp,
		std::ostream& out)
{
  int bk_st_back=Vienna::st_back;
  Vienna::st_back=1;

  Vienna::pf_scale = -1;
  Vienna::init_pf_fold(seq.size());
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);

  double MCC = -10000;
  double SEN, PPV;
  std::string PAREN;

  // stochastic sampling
  for (uint n=0; n!=num_samples; ++n) {
    SCFG::CYKTable<uint> bp_pos(seq.size(), 0);
    char *str = Vienna::pbacktrack(const_cast<char*>(seq.c_str()));
    assert(seq.size()==strlen(str));
    std::string paren (str);
    free(str);

    double sen, ppv, mcc;
    CentroidFold::compute_expected_accuracy (paren, bp, sen, ppv, mcc);
    if (MCC<mcc) {
      PAREN = paren;
      MCC = mcc; SEN = sen; PPV = ppv; 
    }
  }

  Vienna::st_back=bk_st_back;
  Vienna::free_pf_arrays();
  out << ">" << name << "\n" << seq << "\n" << PAREN << " [MCC=" << MCC << "]" << std::endl;
  return num_samples;
}

static
uint
pf_fold_max_mcc(uint num_samples, const std::string& name, const std::string& seq,
                const std::vector<std::string>& ma,
                const CentroidFold::BPTable& bp, std::ostream& out)
{
  int bk_st_back=Vienna::st_back;
  Vienna::st_back=1;

  double MCC = -10000;
  double SEN, PPV;
  std::string PAREN;

  uint sz = seq.size();
  std::vector<std::string>::const_iterator x;
  for (x=ma.begin(); x!=ma.end(); ++x)
  {
    std::string seq;
    std::vector<uint> idxmap(x->size(), static_cast<uint>(-1));
    uint j=0;
    for (uint i=0; i!=x->size(); ++i)
    {
      if ((*x)[i]!='-')
      {
        switch ((*x)[i])
        {
          case 't': seq += 'u'; break;
          case 'T': seq += 'U'; break;
          default:  seq += (*x)[i];
        }
        idxmap[j++] = i;
      }
    }

    Vienna::pf_scale = -1;
    Vienna::init_pf_fold(seq.size());
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
    EncodeBP encoder;

    // stochastic sampling
    for (uint n=0; n!=num_samples; ++n) {
      SCFG::CYKTable<uint> bp_pos(sz, 0);
      char *str = Vienna::pbacktrack(const_cast<char*>(seq.c_str()));
      assert(seq.size()==strlen(str));
      std::string paren(sz, '.');
      for (uint i=0; i!=seq.size(); ++i) paren[idxmap[i]] = str[i];
      free(str);

      double sen, ppv, mcc;
      CentroidFold::compute_expected_accuracy (paren, bp, sen, ppv, mcc);
      if (MCC<mcc) {
        PAREN = paren;
        MCC = mcc; SEN = sen; PPV = ppv; 
      }
    }

    Vienna::free_pf_arrays();
  }

  Vienna::st_back=bk_st_back;
  out << ">" << name << std::endl << seq << std::endl << PAREN << " [MCC=" << MCC << "]" << std::endl;

  return num_samples;
}

#ifdef HAVE_VIENNA18
static
uint
alipf_fold_max_mcc(uint num_samples, const std::string& name, const std::string& seq,
                   const std::vector<std::string>& ma,
                   const CentroidFold::BPTable& bp, std::ostream& out)
{
  double MCC = -10000;
  double SEN, PPV;
  std::string PAREN;

  // prepare an alignment
  uint length = seq.size();
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
    std::string paren(str);
    free(str);

    double sen, ppv, mcc;
    CentroidFold::compute_expected_accuracy (paren, bp, sen, ppv, mcc);
    if (MCC<mcc) {
      PAREN = paren;
      MCC = mcc; SEN = sen; PPV = ppv; 
    }
  }

  Vienna::free_alipf_arrays();
  free(pi);
  for (uint i=0; seqs[i]!=NULL; ++i) delete[] seqs[i];
  //Vienna::fold_constrained = bk;
  delete[] seqs;
  delete[] str2;

  out << ">" << name << std::endl << seq << std::endl << PAREN << " [MCC=" << MCC << "]" << std::endl;

  return num_samples;
}
#endif
#endif

template < class U >
static
uint
contra_fold_max_mcc(uint num_samples, CONTRAfold<U>& cf, 
		    const std::string& name, const std::string& seq, 
		    const CentroidFold::BPTable& bp, 
		    std::ostream& out)
{
  cf.PrepareStochasticTraceback(seq);
  double MCC = -10000;
  double SEN, PPV;
  std::string PAREN;
  for (uint n=0; n!=num_samples; ++n) {
    std::vector<int> paren_i = cf.StochasticTraceback();
    std::string paren (paren_i.size()-1, '.');

    for (uint i=0; i!=paren_i.size(); ++i) {
      if (paren_i[i]>int(i)) {
	paren[i-1] = '(';
	paren[paren_i[i]-1] = ')';
      }
    }
    double sen, ppv, mcc;
    CentroidFold::compute_expected_accuracy (paren, bp, sen, ppv, mcc);
    if (MCC<mcc) {
      PAREN = paren;
      MCC = mcc; SEN = sen; PPV = ppv; 
    }
  }
  out << ">" << name << "\n" << seq << "\n" << PAREN << " [MCC=" << MCC << "]" << std::endl;
  return num_samples;
}

template < class U >
static
uint
contra_fold_max_mcc(uint num_samples, CONTRAfoldM<U>& cf,
                    const std::string& name, const std::string& seq,
                    const std::vector<std::string>& ma,
                    const CentroidFold::BPTable& bp, std::ostream& out)
{
  cf.PrepareStochasticTraceback(ma);

  double MCC = -10000;
  double SEN, PPV;
  std::string PAREN;
  for (uint n=0; n!=num_samples; ++n) {
    std::string p(seq.size(), '.');
    std::vector<int> paren = cf.StochasticTraceback();
    for (uint i=0; i!=paren.size(); ++i) {
      if (paren[i]>int(i)) {
        p[i-1] = '(';
        p[paren[i]-1] = ')';
      }
    }
    double sen, ppv, mcc;
    CentroidFold::compute_expected_accuracy (p, bp, sen, ppv, mcc);
    if (MCC<mcc) {
      PAREN = p;
      MCC = mcc; SEN = sen; PPV = ppv; 
    }
  }

  out << ">" << name << "\n" << seq << "\n" << PAREN << " [MCC=" << MCC << "]" << std::endl;
  return num_samples;
}

template < class U >
static
uint
contra_fold_max_mcc(uint num_samples, CONTRAfold<U>& cf, const std::string& name, const std::string& seq,
                    const std::vector<std::string>& ma,
                    const CentroidFold::BPTable& bp, std::ostream& out)
{
  double MCC = -10000;
  double SEN, PPV;
  std::string PAREN;

  uint sz = seq.size();
  std::vector<std::string>::const_iterator x;
  for (x=ma.begin(); x!=ma.end(); ++x)
  {
    std::string seq;
    std::vector<uint> idxmap(x->size(), static_cast<uint>(-1));
    uint j=0;
    for (uint i=0; i!=x->size(); ++i)
    {
      if ((*x)[i]!='-')
      {
        switch ((*x)[i])
        {
          case 't': seq += 'u'; break;
          case 'T': seq += 'U'; break;
          default:  seq += (*x)[i];
        }
        idxmap[j++] = i;
      }
    }

    // stochastic sampling
    cf.PrepareStochasticTraceback(seq);
    for (uint n=0; n!=num_samples; ++n) {
      SCFG::CYKTable<uint> bp_pos(sz, 0);
      std::vector<int> paren_i = cf.StochasticTraceback();
      std::string paren (sz, '.');
      
      for (uint i=0; i!=paren_i.size(); ++i) {
        if (paren_i[i]>int(i)) {
          paren[idxmap[i-1]] = '(';
          paren[idxmap[paren_i[i]-1]] = ')';
        }
      }

      double sen, ppv, mcc;
      CentroidFold::compute_expected_accuracy (paren, bp, sen, ppv, mcc);
      if (MCC<mcc) {
        PAREN = paren;
        MCC = mcc; SEN = sen; PPV = ppv; 
      }
    }
  }

  out << ">" << name << std::endl << seq << std::endl << PAREN << " [MCC=" << MCC << "]" << std::endl;

  return num_samples;
}

void 
CentroidFold::
max_mcc_fold (const std::string& name, const std::string& seq, uint num_samples, std::ostream& out)
{
  std::string seq2(seq);
  std::replace(seq2.begin(), seq2.end(), 't', 'u');
  std::replace(seq2.begin(), seq2.end(), 'T', 'U');

  // calculate base-pairing probabililty matrix
  calculate_posterior(seq2);

  switch (engine_) {
  case CONTRAFOLD:
    contra_fold_max_mcc (num_samples, *contrafold_, name, seq, bp_, out); 
    break;
#ifdef HAVE_LIBRNA
  case PFFOLD:
    pf_fold_max_mcc (num_samples, name, seq, bp_, out);
    break;
#endif
  default:
    assert (false);
  }
}

void 
CentroidFold::
max_mcc_fold (const Aln& aln, uint num_samples, std::ostream& out)
{
  std::vector<std::string> ma(aln.seq().size());
  std::copy(aln.seq().begin(), aln.seq().end(), ma.begin());
  for (uint i=0; i!=ma.size(); ++i)
  {
    std::replace(ma[i].begin(), ma[i].end(), 't', 'u');
    std::replace(ma[i].begin(), ma[i].end(), 'T', 'U');
  }

  // calculate base-pairing probabililty matrix
  calculate_posterior(ma);

  switch (engine_) {
    case CONTRAFOLD:
      if (dist_type_==0)
      {
        CONTRAfoldM<float>* cf = new CONTRAfoldM<float>(canonical_only_, max_bp_dist_);
        if (!model_.empty()) cf->SetParameters(model_);
        cf->init_rand(seed_);
        contra_fold_max_mcc (num_samples, *cf, aln.name().front(), aln.consensus(), ma, bp_, out); 
        delete cf;
      }
      else
      {
        contra_fold_max_mcc (num_samples, *contrafold_, aln.name().front(), aln.consensus(), ma, bp_, out);
      }
      break;
#ifdef HAVE_LIBRNA
    case PFFOLD:
      pf_fold_max_mcc (num_samples, aln.name().front(), aln.consensus(), ma, bp_, out);
      break;
#ifdef HAVE_VIENNA18
    case ALIPFFOLD:
      alipf_fold_max_mcc(num_samples, aln.name().front(), aln.consensus(), ma, bp_, out);
      break;
#endif    
#endif
  default:
    assert (false);
  }
}

