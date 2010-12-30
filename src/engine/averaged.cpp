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

#include "engine/averaged.h"
#include <stack>
#include <cmath>

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
      inside_traverse(0, bp_.size()-1, *this);
    else
      inside_traverse(0, bp_.size()-1, max_dist_, *this);
  }

private:
  BPTable& bp_;
  const std::list<BPTablePtr>& bps_;
  const std::list<std::vector<uint> >& idxmap_;
  uint max_dist_;
};

AveragedModel::
AveragedModel(FoldingEngine<std::string>* cf, uint max_bp_dist, bool run_as_mea /*=false*/)
  : FoldingEngine<Aln>(run_as_mea, max_bp_dist), cf_(cf), paren_()
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
      CYKTable<uint> bp_pos(aln.size(), 0);
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
      << aln.consensus() << std::endl
      << PAREN << " [MCC=" << MCC << "]" << std::endl;
}
  
void
AveragedModel::
compute_expected_accuracy_sampling(const std::string& paren, const Aln& aln,
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
      CYKTable<uint> bp_pos(aln.size(), 0);
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
