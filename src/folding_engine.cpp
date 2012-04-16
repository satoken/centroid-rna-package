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
#include <iostream>
#include <cmath>
#include <string>
#include <stack>
#include <cstdio>
#include <cassert>
#include <stdexcept>
#include <boost/multi_array.hpp>
#include "folding_engine.h"
#include "centroid.h"
#include "mea.h"
#include "diana.h"
#include "ps_plot.h"

#ifdef HAVE_LIBRNA
#ifndef __INC_LIBRNA_H
#define __INC_LIBRNA_H
namespace Vienna {
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/PS_dot.h>
  extern int eos_debug;
};
};
#endif
#endif

static
const std::string& get_seq(const std::string& t) { return t; }

#ifdef HAVE_LIBRNA
static
float calc_energy(const std::string& t, const std::string& s)
{
#ifdef HAVE_VIENNA20
  return Vienna::energy_of_structure(t.c_str(), s.c_str(), -1);
#else
  Vienna::eos_debug = -1;
  return Vienna::energy_of_struct(t.c_str(), s.c_str());
#endif
}
#endif

static
std::string get_seq(const Aln& aln) { return aln.consensus(); }

static
std::string get_seq (const TH& th) { return th.first; }

#ifdef HAVE_LIBRNA
static
float calc_energy(const Aln& aln, const std::string& s)
{
  Vienna::eos_debug = -1;
  return aln.energy_of_struct(s);
}
#endif

#ifdef HAVE_LIBRNA
static
float calc_energy(const TH& t, const std::string& s)
{
#ifdef HAVE_VIENNA20
  return Vienna::energy_of_structure(t.first.c_str(), s.c_str(), -1);
#else
  Vienna::eos_debug = -1;
  return Vienna::energy_of_struct(t.first.c_str(), s.c_str());
#endif
}
#endif


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
  BPTable table;
  uint p;
};

//////////////////////////////////////////////////////////////////////////////////

template < class SEQ >
void
FoldingEngine<SEQ>::
centroid_fold(const std::string& name, const SEQ& seq,
              const std::vector<float>& gamma, std::ostream& out)
{
#ifdef HAVE_LIBRNA
  Vienna::eos_debug = -1;
#endif
  calculate_posterior(seq);
  out << ">" << name << std::endl << get_seq(seq) << std::endl;
  std::vector<float>::const_iterator g;
  for (g=gamma.begin(); g!=gamma.end(); ++g) {
    std::string paren;
    float p = decode_structure(*g, paren);
    out << paren;
#ifdef HAVE_LIBRNA
    float e = calc_energy(seq, paren.c_str());
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

template < class SEQ >
void
FoldingEngine<SEQ>::
centroid_fold(const std::string& name, const SEQ& seq,
              const std::vector<float>& gamma, std::ostream& out,
              uint num_ea_samples)
{
#ifdef HAVE_LIBRNA
  Vienna::eos_debug = -1;
#endif
  calculate_posterior(seq);

  std::vector<Result> res;
  std::vector<float>::const_iterator g;
  for (g=gamma.begin(); g!=gamma.end(); ++g)
  {
    std::string paren;
    float p = decode_structure(*g, paren);
    double sen = 0.0, ppv = 0.0, mcc = 0.0;
    if (num_ea_samples>0) 
      compute_expected_accuracy_sampling(paren, seq, num_ea_samples, sen, ppv, mcc);
    else
      compute_expected_accuracy(paren, sen, ppv, mcc);
    float e=0;
#ifdef HAVE_LIBRNA
    e = calc_energy(seq, paren.c_str());
#endif
    res.push_back(Result(*g, sen, ppv, mcc, p, e, paren));
  }
  std::sort(res.begin(), res.end(), cmp_by_mcc());

  std::string prev_paren;
  out << ">" << name << std::endl << get_seq(seq) << std::endl;
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
      if (num_ea_samples>=0) {
        out << " [" << "SEN=" << r->sen << ",PPV=" << r->ppv <<  ",MCC=" << r->mcc << "]";
      }
      out << std::endl;
      prev_paren = r->paren;
    }
  }
}

template < class SEQ >
void
FoldingEngine<SEQ>::
stochastic_fold(const SEQ& seq, uint num_samples, std::ostream& out)
{
  prepare_stochastic_traceback(seq);
  for (uint n=0; n!=num_samples; ++n)
  {
    std::string p(seq.size(), '.');
    std::vector<int> paren = stochastic_traceback(seq);
    for (uint i=0; i!=paren.size(); ++i)
    {
      if (paren[i]>int(i))
      {
        p[i-1] = '(';
        p[paren[i]-1] = ')';
      }
    }
    out << p << std::endl;
  }
  clean_stochastic_traceback(seq);
}

template < class SEQ >
void
FoldingEngine<SEQ>::
stochastic_fold(const SEQ& seq, uint num_samples, std::vector<BPvecPtr>& bpv)
{
  EncodeBP encoder;
  prepare_stochastic_traceback(seq);
  for (uint n=0; n!=num_samples; ++n)
  {
    CYKTable<uint> bp_pos(seq.size(), 0);
    std::vector<int> paren = stochastic_traceback(seq);
    for (uint i=0; i!=paren.size(); ++i)
    {
      if (paren[i]>int(i)) bp_pos(i-1, paren[i]-1)++;
    }
    bpv.push_back(encoder.execute(bp_pos));
  }
  clean_stochastic_traceback(seq);
}

inline
uint
hamming_distance(const BPvec& x, const BPvec& y)
{
  return (x^y).count();
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

template < class SEQ >
void
FoldingEngine<SEQ>::
stochastic_fold(const std::string& name, const SEQ& seq, uint num_samples, 
                const std::vector<float>& gamma, uint max_clusters, std::ostream& out, 
                const std::string& p_outname, float th)
{
  std::vector<BPvecPtr> bpv;
  stochastic_fold(seq, num_samples, bpv);
  
  if (max_clusters>=2)
  {
    // calculate the distance matrix
    typedef boost::multi_array<double, 2> DMatrix;
    DMatrix dmatrix(boost::extents[bpv.size()][bpv.size()]);
    for (uint i=0; i!=bpv.size(); ++i)
      for (uint j=i; j!=bpv.size(); ++j)
        dmatrix[i][j]=dmatrix[j][i]=hamming_distance(*bpv[i], *bpv[j]);

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
    for (uint i=0; i!=num.size(); ++i)
    {
      idx[i]=i;
      if (i!=0) s[i]=s[i-1]+num[i-1];
    }
    std::sort(idx.begin(), idx.end(), cmp_by_size(num));

    // centroid estimation for each cluster
    for (uint i=0; i!=idx.size(); ++i)
    {
      std::vector<BPvecPtr> v(num[idx[i]]);
      for (uint j=0; j!=v.size(); ++j) v[j] = bpv[res[s[idx[i]]+j]];
      CountBP< std::vector<BPvecPtr> > count_bp(v, seq.size());
      inside_traverse(0, seq.size()-1, count_bp);
      char buf[100];
      sprintf(buf, "%3.1f%%", static_cast<float>(num[idx[i]])/bpv.size()*100);
      out << ">" << name << " (" << i+1 << " of " << num.size() << ", size="
          << buf << ")" << std::endl
          << get_seq(seq) << std::endl;

      std::vector<float>::const_iterator g;
      for (g=gamma.begin(); g!=gamma.end(); ++g) {
        std::string paren(seq.size(), '.');
        Centroid::execute(count_bp.table, paren, *g);
        out << paren << " (g=" << *g << ",th=" << (1.0/(1.0+*g));
#ifdef HAVE_LIBRNA
        out << ",e=" << calc_energy(seq, paren);
#endif
        out << ")" << std::endl;
      }

      if (!p_outname.empty())
      {
        char buf[1024];
        snprintf(buf, sizeof(buf), "%s-%d", p_outname.c_str(), i);
        count_bp.table.save(buf, get_seq(seq), th);
      }
    }
  }
  else
  {
    CountBP< std::vector<BPvecPtr> > count_bp(bpv, seq.size());
    inside_traverse(0, seq.size()-1, count_bp);
    std::string paren(seq.size(), '.');
    std::vector<float>::const_iterator g;
    out << ">" << name << std::endl
        << get_seq(seq) << std::endl;
    for (g=gamma.begin(); g!=gamma.end(); ++g)
    {
      std::fill(paren.begin(), paren.end(), '.');
      Centroid::execute(count_bp.table, paren, *g);
      out << paren << " (g=" << *g << ",th=" << (1.0/(1.0+*g));
#ifdef HAVE_LIBRNA
      out << ",e=" << calc_energy(seq, paren);
#endif
      out << ")" << std::endl;
    }

    if (!p_outname.empty())
    {
      char buf[1024];
      snprintf(buf, sizeof(buf), "%s-1", p_outname.c_str());
      count_bp.table.save(buf, get_seq(seq), th);
    }
  }
}

template <class SEQ>
void
FoldingEngine<SEQ>::
max_mcc_fold(const std::string& name, const SEQ& seq, std::ostream& out, uint num_samples)
{
  double MCC = -1.0;
  double SEN, PPV;
  std::string PAREN;

  calculate_posterior(seq);
  prepare_stochastic_traceback(seq);
  for (uint n=0; n!=num_samples; ++n)
  {
    std::string p(seq.size(), '.');
    std::vector<int> paren = stochastic_traceback(seq);
    for (uint i=0; i!=paren.size(); ++i)
    {
      if (paren[i]>int(i))
      {
        p[i-1] = '(';
        p[paren[i]-1] = ')';
      }
    }
    double sen, ppv, mcc;
    compute_expected_accuracy(p, sen, ppv, mcc);
    if (MCC<mcc) {
      PAREN = p;
      MCC = mcc; SEN = sen; PPV = ppv; 
    }
  }
  clean_stochastic_traceback(seq);

  out << ">" << name << "\n" << get_seq(seq) << "\n" << PAREN << " [MCC=" << MCC << "]" << std::endl;
}

template <class SEQ>
void
FoldingEngine<SEQ>::
compute_expected_accuracy(const std::string& paren,
                          double& sen, double& ppv, double& mcc) const
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

  if (bp_.size() != paren.size()) 
    std::cerr << bp_.size() << "!=" << paren.size() << std::endl << paren << std::endl;; 
  assert (bp_.size()==paren.size());

  for (uint i=0; i<bp_.size(); ++i) {
    for (uint j = i+1; j<bp_.size(); ++j) {
      double p = bp_(i,j);
      sump += p;
    }
  }

  for (std::vector<std::pair<int, int> >::const_iterator it=BP.begin();
       it != BP.end(); ++it) {
    etp += bp_(it->first, it->second);
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

template <class SEQ>
//static
void
FoldingEngine<SEQ>::
get_TP_TN_FP_FN (const BPvec& ref, const BPvec& pre, 
		 double& TP, double& TN, double& FP, double& FN)
{
  TP = (ref&pre).count();
  FP = pre.count() - TP;
  FN = ref.count() - TP;
  TN = ref.size() - ref.count() - FP; 
}

template <class SEQ>
void
FoldingEngine<SEQ>::
compute_expected_accuracy_sampling(const std::string& paren, const SEQ& seq,
                                   uint num_ea_samples, double& sen, double& ppv, double& mcc)
{
  CYKTable<uint> bp_pos(paren.size(), 0);
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

  sen = ppv = mcc = 0.0;
  int N = 0;

  prepare_stochastic_traceback(seq);
  for (int n=0; n!=int(num_ea_samples); ++n) {
    CYKTable<uint> bp_pos(seq.size(), 0);
    std::vector<int> paren = stochastic_traceback(seq);
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
  clean_stochastic_traceback(seq);

  sen /= N;
  ppv /= N;
  mcc /= N;
}

template <class SEQ>
float
FoldingEngine<SEQ>::
decode_structure(float gamma, std::string& paren) const
{
  float p=0.0;
  paren.resize(bp_.size());
  std::fill(paren.begin(), paren.end(), '.');
  if (!mea_) {
    if (max_bp_dist_==0)
      p = Centroid::execute(bp_, paren, gamma);
    else
      p = Centroid::execute(bp_, paren, max_bp_dist_, gamma);
  } else {
    if (max_bp_dist_==0)
      p = MEA::execute(bp_, paren, gamma);
    else
      p = MEA::execute(bp_, paren, max_bp_dist_, gamma);
  }
  return p;
}

template <class SEQ>
void
FoldingEngine<SEQ>::
ps_plot(const std::string& name, const SEQ& seq, float g, bool color) const
{
  std::string paren;
  decode_structure(g, paren);
  if (color)
    ::ps_color_plot(get_seq(seq).c_str(), paren.c_str(), bp_, name.c_str());
  else
    ::ps_plot(get_seq(seq).c_str(), paren.c_str(), name.c_str());
}

#ifdef HAVE_LIBRNA
template <class SEQ>
void
FoldingEngine<SEQ>::
svg_plot(const std::string& name, const SEQ& seq, float g) const
{
  std::string paren;
  decode_structure(g, paren);
  Vienna::svg_rna_plot(const_cast<char*>(get_seq(seq).c_str()),
		       const_cast<char*>(paren.c_str()),
                       const_cast<char*>(name.c_str()));
}
#endif

template < class SEQ >
void
FoldingEngine<SEQ>::
set_constraint(const std::string& str)
{
  throw std::logic_error("unsupported method: set_constraint()");
}
  
template < class SEQ >
void
FoldingEngine<SEQ>::
prepare_stochastic_traceback(const SEQ& seq)
{
  throw std::logic_error("unsupported method: prepare_stochastic_traceback()");
}

template < class SEQ >
std::vector<int>
FoldingEngine<SEQ>::
stochastic_traceback(const SEQ& seq)
{
  throw std::logic_error("unsupported method: stochastic_traceback()");
  return std::vector<int>();
}

template < class SEQ >
void
FoldingEngine<SEQ>::
clean_stochastic_traceback(const SEQ& seq)
{
  throw std::logic_error("unsupported method: clean_stochastic_traceback()");
}

// instantiation
template class FoldingEngine<std::string>;
template class FoldingEngine<Aln>;
template class FoldingEngine<TH>;
