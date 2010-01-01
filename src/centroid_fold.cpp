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
#include <stdexcept>
#include <sys/time.h>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/multi_array.hpp>
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

const std::string& get_seq(const std::string& t) { return t; }
#ifdef HAVE_LIBRNA
float calc_energy(const std::string& t, const std::string& s)
{
  Vienna::eos_debug = -1;
  return Vienna::energy_of_struct(t.c_str(), s.c_str());
}
#endif

std::string get_seq(const Aln& aln) { return aln.consensus(); }
#ifdef HAVE_LIBRNA
float calc_energy(const Aln& aln, const std::string& s)
{
  Vienna::eos_debug = -1;
  return aln.energy_of_struct(s);
}
#endif

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
  BPTable table;
  uint p;
};

//////////////////////////////////////////////////////////////////////////////////

template < class SEQ >
void
CentroidFold<SEQ>::
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
CentroidFold<SEQ>::
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
CentroidFold<SEQ>::
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
CentroidFold<SEQ>::
stochastic_fold(const SEQ& seq, uint num_samples, std::vector<BPvecPtr>& bpv)
{
  EncodeBP encoder;
  prepare_stochastic_traceback(seq);
  for (uint n=0; n!=num_samples; ++n)
  {
    SCFG::CYKTable<uint> bp_pos(seq.size(), 0);
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
CentroidFold<SEQ>::
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
      SCFG::inside_traverse(0, seq.size()-1, count_bp);
      char buf[100];
      sprintf(buf, "%3.1f%%", static_cast<float>(num[idx[i]])/bpv.size()*100);
      out << ">" << name << " (" << i+1 << " of " << num.size() << ", size="
          << buf << ")" << std::endl
          << get_seq(seq) << std::endl;

      std::vector<float>::const_iterator g;
      for (g=gamma.begin(); g!=gamma.end(); ++g) {
        std::string paren(seq.size(), '.');
        SCFG::Centroid::execute(count_bp.table, paren, *g);
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
    SCFG::inside_traverse(0, seq.size()-1, count_bp);
    std::string paren(seq.size(), '.');
    std::vector<float>::const_iterator g;
    out << ">" << name << std::endl
        << get_seq(seq) << std::endl;
    for (g=gamma.begin(); g!=gamma.end(); ++g)
    {
      std::fill(paren.begin(), paren.end(), '.');
      SCFG::Centroid::execute(count_bp.table, paren, *g);
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
CentroidFold<SEQ>::
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
CentroidFold<SEQ>::
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

inline
static
void
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
CentroidFold<SEQ>::
compute_expected_accuracy_sampling(const std::string& paren, const SEQ& seq,
                                   uint num_ea_samples, double& sen, double& ppv, double& mcc)
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

  sen = ppv = mcc = 0.0;
  int N = 0;

  prepare_stochastic_traceback(seq);
  for (int n=0; n!=int(num_ea_samples); ++n) {
    SCFG::CYKTable<uint> bp_pos(seq.size(), 0);
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
CentroidFold<SEQ>::
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

template <class SEQ>
void
CentroidFold<SEQ>::
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
CentroidFold<SEQ>::
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
CentroidFold<SEQ>::
set_constraint(const std::string& str)
{
  throw std::logic_error("unsupported method: set_constraint()");
}
  
template < class SEQ >
void
CentroidFold<SEQ>::
prepare_stochastic_traceback(const SEQ& seq)
{
  throw std::logic_error("unsupported method: prepare_stochastic_traceback()");
}

template < class SEQ >
std::vector<int>
CentroidFold<SEQ>::
stochastic_traceback(const SEQ& seq)
{
  throw std::logic_error("unsupported method: stochastic_traceback()");
  return std::vector<int>();
}

template < class SEQ >
void
CentroidFold<SEQ>::
clean_stochastic_traceback(const SEQ& seq)
{
  throw std::logic_error("unsupported method: clean_stochastic_traceback()");
}

// instantiation
template class CentroidFold<std::string>;
template class CentroidFold<Aln>;

//////////////////////////////////////////////////////////////////////////////////

CONTRAfoldModel::
CONTRAfoldModel(const std::string& model, bool canonical_only, uint max_bp_dist,
                uint seed /*=0*/, bool run_as_mea /*=false*/)
  : CentroidFold<std::string>(run_as_mea, max_bp_dist),
    contrafold_(NULL)
{
  contrafold_ = new CONTRAfold<float>(canonical_only, max_bp_dist);
  if (seed==0)
  {
    timeval t;
    gettimeofday(&t, NULL);
    seed = t.tv_usec * t.tv_sec;
  }
  contrafold_->init_rand(seed);
  if (!model.empty()) contrafold_->SetParameters(model);
}

CONTRAfoldModel::
~CONTRAfoldModel()
{
  if (contrafold_) delete contrafold_;
}

void
CONTRAfoldModel::
set_constraint(const std::string& str)
{
  contrafold_->SetConstraint(str);
}

inline
void
CONTRAfoldModel::
calculate_posterior(const std::string& seq)
{
  bp_.resize(seq.size(), contrafold_->max_bp_dist());
  std::vector<float> posterior;
  contrafold_->ComputePosterior(seq, posterior);

  if (contrafold_->max_bp_dist()==0) {
    uint k=0;
    for (uint i=0; i!=seq.size()+1; ++i) {
      for (uint j=i; j!=seq.size()+1; ++j) {
        if (i!=0) bp_.update(i-1, j-1, posterior[k]);
        ++k;
      }
    }
  } else {
    uint k=0;
    for (uint i=0; i!=seq.size()+1; ++i) {
      for (uint j=i; j!=i+contrafold_->max_bp_dist(); ++j) {
        if (i!=0 && j<=seq.size()) bp_.update(i-1, j-1, posterior[k]);
        ++k;
      }
    }
  }
}

void
CONTRAfoldModel::
prepare_stochastic_traceback(const std::string& seq)
{
  contrafold_->PrepareStochasticTraceback(seq);
}

std::vector<int>
CONTRAfoldModel::
stochastic_traceback(const std::string& seq)
{
  return contrafold_->StochasticTraceback();
}

//////////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_LIBRNA

McCaskillModel::
McCaskillModel(bool canonical_only, uint max_bp_dist,
               uint seed /*=0*/, bool run_as_mea /*=false*/)
  : CentroidFold<std::string>(run_as_mea, max_bp_dist),
    canonical_only_(canonical_only), bk_st_back_(Vienna::st_back), str_()
{
  if (!canonical_only_)
    Vienna::nonstandards = const_cast<char*>("AAACAGCACCCUGAGGUCUU");

  if (seed==0)
  {
    timeval t;
    gettimeofday(&t, NULL);
    seed = t.tv_usec * t.tv_sec;
  }
  srand((uint)seed);

  {
    using namespace Vienna;
    // copy from ViennaRNA-x.x.x/lib/utils.c
    xsubi[0] = xsubi[1] = xsubi[2] = (unsigned short) seed;  /* lower 16 bit */
    xsubi[1] += (unsigned short) ((unsigned)seed >> 6);
    xsubi[2] += (unsigned short) ((unsigned)seed >> 12);
  }
}

void
McCaskillModel::
set_constraint(const std::string& str)
{
  str_=str;
  std::replace(str_.begin(), str_.end(), '.', 'x');
  std::replace(str_.begin(), str_.end(), '?', '.');
}

void
McCaskillModel::
calculate_posterior(const std::string& seq)
{
  bp_.resize(seq.size());
  Vienna::pf_scale = -1;
  Vienna::init_pf_fold(seq.size());
  if (str_.empty()) {
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  } else {
    char *str2 = new char[str_.size()+1];
    strcpy(str2, str_.c_str());
    int bk = Vienna::fold_constrained;
    Vienna::fold_constrained = 1;
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), str2);
    delete[] str2;
    Vienna::fold_constrained = bk;
  }
  for (uint j=2; j!=bp_.size()+1; ++j) {
    for (uint i=j-1; ; --i) {
      bp_.update(i-1, j-1, Vienna::pr[Vienna::iindx[i]-j]);
      if (i==1) break;
    }
  }
  Vienna::free_pf_arrays();
}  

void
McCaskillModel::
prepare_stochastic_traceback(const std::string& seq)
{
  bk_st_back_=Vienna::st_back;
  Vienna::st_back=1;
  Vienna::pf_scale = -1;
  Vienna::init_pf_fold(seq.size());
  if (str_.empty()) {
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  } else {
    char *str2 = new char[str_.size()+1];
    strcpy(str2, str_.c_str());
    int bk = Vienna::fold_constrained;
    Vienna::fold_constrained = 1;
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), str2);
    delete[] str2;
    Vienna::fold_constrained = bk;
  }
}  

std::vector<int>
McCaskillModel::
stochastic_traceback(const std::string& seq)
{
  char *str = Vienna::pbacktrack(const_cast<char*>(seq.c_str()));
  std::vector<int> paren(seq.size()+1, 0);
  std::stack<uint> st;
  for (uint i=0; i!=seq.size(); ++i) {
    if (str[i]=='(') {
      st.push(i);
    } else if (str[i]==')') {
      paren[st.top()+1]=i+1;
      paren[i+1]=st.top()+1;
      st.pop();
    }
  }
  free(str);
  return paren;
}

void
McCaskillModel::
clean_stochastic_traceback(const std::string& seq)
{
  Vienna::st_back=bk_st_back_;
  Vienna::free_pf_arrays();
}
#endif

//////////////////////////////////////////////////////////////////////////////////

static
void
remove_gaps(const std::string& seq_in, std::string& seq_out, std::vector<uint>& idx)
{
  idx.resize(seq_in.size());
  std::fill(idx.begin(), idx.end(), static_cast<uint>(-1));
  seq_out.clear();
  for (uint i=0, j=0; i!=seq_in.size(); ++i)
  {
    if (seq_in[i]!='-')
    {
      seq_out += seq_in[i];
      idx[i] = j++;
    }
  }
}

static
std::string
remove_gaps(const std::string& seq_in, std::string& seq_out, std::vector<uint>& idx,
            const std::vector<uint>& paren)
{
  remove_gaps(seq_in, seq_out, idx);

  std::string str(seq_out.size(), '.');
  for (uint i=0; i!=paren.size(); ++i)
  {
    if (paren[i]!=static_cast<uint>(-1) && int(i)<int(paren[i]) &&
        idx[i]!=static_cast<uint>(-1) && idx[paren[i]]!=static_cast<uint>(-1))
    {
      str[idx[i]]='(';
      str[idx[paren[i]]]=')';
    }
  }
  return str;
}

struct AverageBP
{
  AverageBP(BPTable& bp,
	    const std::list<BPTablePtr>& bps,
	    const std::list<std::vector<uint> >& idxmap,
            uint max_dist)
    : bp_(bp), bps_(bps), idxmap_(idxmap), max_dist_(max_dist)
  {
  }

  void operator()(uint i, uint j)
  {
    double v=0.0;
    uint c=0;
    std::list<BPTablePtr>::const_iterator b = bps_.begin();
    std::list<std::vector<uint> >::const_iterator idx = idxmap_.begin();
    while (b!=bps_.end() && idx!=idxmap_.end()) {
      uint ii=(*idx)[i];
      uint jj=(*idx)[j];
      if (ii!=static_cast<uint>(-1) && jj!=static_cast<uint>(-1)) v += (**b)(ii,jj);
      ++c; ++b; ++idx;
    }
    bp_.update(i, j, v/c);
  }

  void make()
  {
    if (max_dist_==0)
      SCFG::inside_traverse(0, bp_.size()-1, *this);
    else
      SCFG::inside_traverse(0, bp_.size()-1, max_dist_, *this);
  }

private:
  BPTable& bp_;
  const std::list<BPTablePtr>& bps_;
  const std::list<std::vector<uint> >& idxmap_;
  uint max_dist_;
};

AveragedModel::
AveragedModel(CentroidFold<std::string>* cf, uint max_bp_dist, bool run_as_mea /*=false*/)
  : CentroidFold<Aln>(run_as_mea, max_bp_dist), cf_(cf), paren_()
{
}

void
AveragedModel::
stochastic_fold(const Aln& aln, uint num_samples, std::ostream& out)
{
  std::list<std::string>::const_iterator s;
  for (s=aln.seq().begin(); s!=aln.seq().end(); ++s)
  {
    std::string seq;
    std::vector<uint> idx;
    if (paren_.empty())
      remove_gaps(*s, seq, idx);
    else
      cf_->set_constraint(remove_gaps(*s, seq, idx, paren_));
    std::vector<uint> rev(seq.size());
    for (uint i=0; i!=idx.size(); ++i)
      if (idx[i]!=static_cast<uint>(-1)) rev[idx[i]]=i;

    cf_->prepare_stochastic_traceback(seq);
    for (uint n=0; n!=num_samples/aln.num_aln(); ++n)
    {
      std::string p(aln.size(), '.');
      std::vector<int> paren = cf_->stochastic_traceback(seq);
      for (uint i=0; i!=paren.size(); ++i)
      {
        if (paren[i]>int(i))
        {
          p[rev[i-1]] = '(';
          p[rev[paren[i]-1]] = ')';
        }
      }
      out << p << std::endl;
    }
    cf_->clean_stochastic_traceback(seq);
  }
}

void
AveragedModel::
stochastic_fold(const Aln& aln, uint num_samples, std::vector<BPvecPtr>& bpv)
{
  EncodeBP encoder;
  std::list<std::string>::const_iterator s;
  for (s=aln.seq().begin(); s!=aln.seq().end(); ++s)
  {
    std::string seq;
    std::vector<uint> idx;
    if (paren_.empty())
      remove_gaps(*s, seq, idx);
    else
      cf_->set_constraint(remove_gaps(*s, seq, idx, paren_));
    std::vector<uint> rev(seq.size());
    for (uint i=0; i!=idx.size(); ++i)
      if (idx[i]!=static_cast<uint>(-1)) rev[idx[i]]=i;

    cf_->prepare_stochastic_traceback(seq);
    for (uint n=0; n!=num_samples/aln.num_aln(); ++n)
    {
      SCFG::CYKTable<uint> bp_pos(aln.size(), 0);
      std::vector<int> paren = cf_->stochastic_traceback(seq);
      for (uint i=0; i!=paren.size(); ++i)
      {
        if (paren[i]>int(i))
          bp_pos(rev[i-1], rev[paren[i]-1])++;
      }
      bpv.push_back(encoder.execute(bp_pos));
    }
    cf_->clean_stochastic_traceback(seq);
  }
}

void
AveragedModel::
max_mcc_fold(const std::string& name, const Aln& aln, std::ostream& out, uint num_samples)
{
  double MCC = -1.0;
  double SEN, PPV;
  std::string PAREN;

  calculate_posterior(aln);
  std::list<std::string>::const_iterator s;
  for (s=aln.seq().begin(); s!=aln.seq().end(); ++s)
  {
    std::string seq;
    std::vector<uint> idx;
    if (paren_.empty())
      remove_gaps(*s, seq, idx);
    else
      cf_->set_constraint(remove_gaps(*s, seq, idx, paren_));
    std::vector<uint> rev(seq.size());
    for (uint i=0; i!=idx.size(); ++i)
      if (idx[i]!=static_cast<uint>(-1)) rev[idx[i]]=i;

    cf_->prepare_stochastic_traceback(seq);
    for (uint n=0; n!=num_samples/aln.num_aln(); ++n)
    {
      std::string p(aln.size(), '.');
      std::vector<int> paren = cf_->stochastic_traceback(seq);
      for (uint i=0; i!=paren.size(); ++i)
      {
        if (paren[i]>int(i))
        {
          p[rev[i-1]] = '(';
          p[rev[paren[i]-1]] = ')';
        }
      }
      double sen, ppv, mcc;
      compute_expected_accuracy(p, sen, ppv, mcc);
      if (MCC<mcc) {
        PAREN = p;
        MCC = mcc; SEN = sen; PPV = ppv; 
      }
    }
    cf_->clean_stochastic_traceback(seq);
  }
  out << ">" << name << std::endl
      << get_seq(aln) << std::endl
      << PAREN << " [MCC=" << MCC << "]" << std::endl;
}
  
void
AveragedModel::
compute_expected_accuracy_sampling(const std::string& paren, const Aln& aln,
                                   uint num_ea_samples, double& sen, double& ppv, double& mcc)
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

  sen = ppv = mcc = 0.0;
  int N = 0;

  std::list<std::string>::const_iterator s;
  for (s=aln.seq().begin(); s!=aln.seq().end(); ++s)
  {
    std::string seq;
    std::vector<uint> idx;
    if (paren_.empty())
      remove_gaps(*s, seq, idx);
    else
      cf_->set_constraint(remove_gaps(*s, seq, idx, paren_));
    std::vector<uint> rev(seq.size());
    for (uint i=0; i!=idx.size(); ++i)
      if (idx[i]!=static_cast<uint>(-1)) rev[idx[i]]=i;

    cf_->prepare_stochastic_traceback(seq);
    for (int n=0; n!=int(num_ea_samples/aln.num_aln()); ++n) {
      SCFG::CYKTable<uint> bp_pos(aln.size(), 0);
      std::vector<int> paren = cf_->stochastic_traceback(seq);
      for (uint i=0; i!=paren.size(); ++i) {
        if (paren[i]>int(i))
          bp_pos(rev[i-1], rev[paren[i]-1])++;
      }
      double TP, TN, FP, FN;
      get_TP_TN_FP_FN (*encoder.execute(bp_pos), *bp, TP, TN, FP, FN);
      if (TP+FN!=0) sen += (double) TP / (TP+FN);
      if (TP+FP!=0) ppv += (double) TP / (TP+FP);
      if (TP+FP!=0 && TP+FN!=0 && TN+FP!=0 && TN+FN!=0)
        mcc += (TP*TN-FP*FN)/std::sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
      ++N;
    }
    cf_->clean_stochastic_traceback(seq);
  }

  sen /= N;
  ppv /= N;
  mcc /= N;
}

void
AveragedModel::
set_constraint(const std::string& str)
{
  if (str.empty())
    paren_.clear();
  else
  {
    paren_.resize(str.size());
    std::fill(paren_.begin(), paren_.end(), static_cast<uint>(-1));
    std::stack<uint> st;
    for (uint i=0; i<str.size(); ++i) {
      if (str[i]=='(') {
        st.push(i);
      } else if (str[i]==')') {
        paren_[i]=st.top();
        paren_[st.top()]=i;
        st.pop();
      }
    }
  }
}

void
AveragedModel::
calculate_posterior(const Aln& aln)
{
  std::list<BPTablePtr> bps;
  std::list<std::vector<uint> > idxmaps;
  std::list<std::string>::const_iterator s;
  for (s=aln.seq().begin(); s!=aln.seq().end(); ++s)
  {
    std::string seq;
    std::vector<uint> idx;
    if (paren_.empty())
      remove_gaps(*s, seq, idx);
    else
      cf_->set_constraint(remove_gaps(*s, seq, idx, paren_));
    cf_->calculate_posterior(seq);
    BPTablePtr bpi(new BPTable(cf_->get_bp()));
    bps.push_back(bpi);
    idxmaps.push_back(idx);
  }
  bp_.resize(idxmaps.front().size(), max_bp_dist_);
  AverageBP avg(bp_, bps, idxmaps, max_bp_dist_);
  avg.make();
}

//////////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_LIBRNA

AliFoldModel::
AliFoldModel(bool canonical_only, uint max_bp_dist, uint seed /*=0*/, bool run_as_mea /*=false*/)
  : CentroidFold<Aln>(run_as_mea, max_bp_dist), canonical_only_(canonical_only)
{
  if (!canonical_only_)
    Vienna::nonstandards = const_cast<char*>("AAACAGCACCCUGAGGUCUU");

  if (seed==0)
  {
    timeval t;
    gettimeofday(&t, NULL);
    seed = t.tv_usec * t.tv_sec;
  }
  srand((uint)seed);

  {
    using namespace Vienna;
    // copy from ViennaRNA-x.x.x/lib/utils.c
    xsubi[0] = xsubi[1] = xsubi[2] = (unsigned short) seed;  /* lower 16 bit */
    xsubi[1] += (unsigned short) ((unsigned)seed >> 6);
    xsubi[2] += (unsigned short) ((unsigned)seed >> 12);
  }
}

void
AliFoldModel::
set_constraint(const std::string& str)
{
  str_ = str;
  std::replace(str_.begin(), str_.end(), '.', 'x');
  std::replace(str_.begin(), str_.end(), '?', '.');
}

char**
alloc_aln(const Aln& aln)
{
  char** seqs = new char*[aln.num_aln()+1];
  seqs[aln.num_aln()]=NULL;
  std::list<std::string>::const_iterator x;
  uint i=0;
  for (x=aln.seq().begin(); x!=aln.seq().end(); ++x)
  {
    seqs[i] = new char[aln.size()+1];
    strcpy(seqs[i++], x->c_str());
  }
  return seqs;
}

void
free_aln(char** seqs)
{
  for (uint i=0; seqs[i]!=NULL; ++i) delete[] seqs[i];
  delete[] seqs;
}

void
AliFoldModel::
calculate_posterior(const Aln& aln)
{
  bp_.resize(aln.size());
  char** seqs=alloc_aln(aln);
  char* str2=new char[aln.size()+1];
  int bk=Vienna::fold_constrained;
  
  if (!str_.empty())
  {
    Vienna::fold_constrained=1;    
    strcpy(str2, str_.c_str());
  }
  // scaling parameters to avoid overflow
  double min_en = Vienna::alifold(seqs, str2);
  double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(1.07*min_en)/kT/aln.size());
  Vienna::free_alifold_arrays();

#ifdef HAVE_VIENNA18
  Vienna::plist* pi;
#else
  Vienna::pair_info* pi;
#endif
  if (!str_.empty()) strcpy(str2, str_.c_str());
  Vienna::alipf_fold(seqs, str2, &pi);
  for (uint k=0; pi[k].i!=0; ++k)
    bp_.update(pi[k].i-1, pi[k].j-1, pi[k].p);
  free(pi);

  Vienna::free_alipf_arrays();
  if (str2) delete[] str2;
  free_aln(seqs);
  Vienna::fold_constrained = bk;
}

void
AliFoldModel::
prepare_stochastic_traceback(const Aln& aln)
{
  char** seqs=alloc_aln(aln);
  char* str2=new char[aln.size()+1];
  int bk=Vienna::fold_constrained;
  
  if (!str_.empty())
  {
    Vienna::fold_constrained=1;    
    strcpy(str2, str_.c_str());
  }
  // scaling parameters to avoid overflow
  double min_en = Vienna::alifold(seqs, str2);
  double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(1.07*min_en)/kT/aln.size());
  Vienna::free_alifold_arrays();

#ifdef HAVE_VIENNA18
  Vienna::plist* pi;
#else
  Vienna::pair_info* pi;
#endif
  if (!str_.empty()) strcpy(str2, str_.c_str());
  Vienna::alipf_fold(seqs, str2, &pi);

  if (str2) delete[] str2;
  free_aln(seqs);
  Vienna::fold_constrained = bk;
}

std::vector<int>
AliFoldModel::
stochastic_traceback(const Aln& aln)
{
  double prob=1;
  char *str = Vienna::alipbacktrack(&prob);
  std::vector<int> paren(aln.size()+1, 0);
  std::stack<uint> st;
  for (uint i=0; i!=aln.size(); ++i) {
    if (str[i]=='(') {
      st.push(i);
    } else if (str[i]==')') {
      paren[st.top()+1]=i+1;
      paren[i+1]=st.top()+1;
      st.pop();
    }
  }
  free(str);
  return paren;
}

void
AliFoldModel::
clean_stochastic_traceback(const Aln& aln)
{
  Vienna::st_back=bk_st_back_;
  Vienna::free_alipf_arrays();
}

#endif

//////////////////////////////////////////////////////////////////////////////////

CONTRAfoldMultiModel::
CONTRAfoldMultiModel(const std::string& model, bool canonical_only, uint max_bp_dist,
                     uint seed /*=0*/, bool run_as_mea /*=false*/)
  : CentroidFold<Aln>(run_as_mea, max_bp_dist),
    contrafold_(NULL), contrafoldm_(NULL), avg_(NULL)
{
  contrafold_ = new CONTRAfoldModel(model, canonical_only, max_bp_dist, seed, run_as_mea);

  contrafoldm_ = new CONTRAfoldM<float>(canonical_only, max_bp_dist);
  if (seed==0)
  {
    timeval t;
    gettimeofday(&t, NULL);
    seed = t.tv_usec * t.tv_sec;
  }
  contrafoldm_->init_rand(seed);
  if (!model.empty()) contrafoldm_->SetParameters(model);

  avg_ = new AveragedModel(contrafold_, max_bp_dist, run_as_mea);
}

CONTRAfoldMultiModel::
~CONTRAfoldMultiModel()
{
  if (contrafold_) delete contrafold_;
  if (contrafoldm_) delete contrafoldm_;
  if (avg_) delete avg_;
}

void
CONTRAfoldMultiModel::
calculate_posterior(const Aln& aln)
{
  avg_->calculate_posterior(aln);
  bp_=avg_->get_bp();
}

void
CONTRAfoldMultiModel::
prepare_stochastic_traceback(const Aln& aln)
{
  std::vector<std::string> seqs(aln.seq().size());
  std::copy(aln.seq().begin(), aln.seq().end(), seqs.begin());
  contrafoldm_->PrepareStochasticTraceback(seqs);
}

std::vector<int>
CONTRAfoldMultiModel::
stochastic_traceback(const Aln& aln)
{
  return contrafoldm_->StochasticTraceback();
}

//////////////////////////////////////////////////////////////////////////////////

AuxModel::
AuxModel(const std::vector<std::string>& bpfiles, bool run_as_mea /*=false*/)
  : CentroidFold<std::string>(run_as_mea, 0), bpfiles_(bpfiles), pos_(0)
{
}

void
AuxModel::
calculate_posterior(const std::string& seq)
{
  if (bpfiles_.size()<=pos_) throw "More BP files are needed.";
  bp_.load(bpfiles_[pos_++].c_str());
}
