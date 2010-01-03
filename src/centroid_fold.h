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
#ifndef __INC_CENTROID_FOLD_H__
#define __INC_CENTROID_FOLD_H__

#include <string>
#include <list>
#include <vector>
#include <iosfwd>
#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include "bp.h"
#include "aln.h"

typedef unsigned int uint;
typedef SCFG::BP::Table<float> BPTable;
typedef boost::shared_ptr<BPTable> BPTablePtr;
typedef boost::dynamic_bitset<> BPvec;
typedef boost::shared_ptr<BPvec> BPvecPtr;

template < class SEQ >
class CentroidFold
{
public:
  CentroidFold(bool run_as_mea=false, uint max_bp_dist=0)
    : mea_(run_as_mea), bp_(), max_bp_dist_(max_bp_dist) { }
  virtual ~CentroidFold() { }

  // folding routines
  void centroid_fold(const std::string& name, const SEQ& seq,
                     const std::vector<float>& gamma, std::ostream& out);
  void centroid_fold(const std::string& name, const SEQ& seq,
                     const std::vector<float>& gamma, std::ostream& out, uint num_ea_samples);
  virtual void stochastic_fold(const SEQ& seq, uint num_samples, std::ostream& out);
  virtual void stochastic_fold(const SEQ& seq, uint num_samples, std::vector<BPvecPtr>& bpv);
  void stochastic_fold(const std::string& name, const SEQ& seq, uint num_samples,
                       const std::vector<float>& gamma, uint max_clusters, std::ostream& out,
                       const std::string& p_outname, float th);
  virtual void max_mcc_fold(const std::string& name, const SEQ& seq, std::ostream& out, uint num_samples);
  
  void compute_expected_accuracy(const std::string& paren, double& sen, double& ppv, double& mcc) const;
  virtual void compute_expected_accuracy_sampling(const std::string& paren, const SEQ& seq,
                                                  uint num_ea_samples, double& sen, double& ppv, double& mcc);

  float decode_structure(float gamma, std::string& paren) const;
  
  const BPTable& get_bp() const { return bp_; }

  void ps_plot(const std::string& name, const SEQ& seq, float g, bool color=true) const;
#ifdef HAVE_LIBRNA
  void svg_plot(const std::string& name, const SEQ& seq, float g) const;
#endif

  // interfaces
  virtual void set_constraint(const std::string& str);
  virtual void calculate_posterior(const SEQ& seq) = 0;
  virtual void prepare_stochastic_traceback(const SEQ& seq);
  virtual std::vector<int> stochastic_traceback(const SEQ& seq);
  virtual void clean_stochastic_traceback(const SEQ& seq);

protected:
  bool mea_;
  BPTable bp_;
  uint max_bp_dist_;

protected:
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

  static
  void get_TP_TN_FP_FN (const BPvec& ref, const BPvec& pre, 
                        double& TP, double& TN, double& FP, double& FN);
};
#endif	// #ifndef __INC_CENTROID_FOLD_H__
