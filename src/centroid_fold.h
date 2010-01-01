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
#if 0
#include <contrafold/contrafold.h>
#else
template < class RealT > class CONTRAfold;
template < class RealT > class CONTRAfoldM;
#endif

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
};

class CONTRAfoldModel : public CentroidFold<std::string>
{
public:
  CONTRAfoldModel(const std::string& model, bool canonical_only, uint max_bp_dist,
                  uint seed=0, bool run_as_mea=false);
  virtual ~CONTRAfoldModel();

  // interface implementations
  virtual void set_constraint(const std::string& str);
  virtual void calculate_posterior(const std::string& seq);
  virtual void prepare_stochastic_traceback(const std::string& seq);
  virtual std::vector<int> stochastic_traceback(const std::string& seq);
  virtual void clean_stochastic_traceback(const std::string& seq) { }

private:
  CONTRAfold<float>* contrafold_;
};

#ifdef HAVE_LIBRNA
class McCaskillModel : public CentroidFold<std::string>
{
public:
  McCaskillModel(bool canonical_only, uint max_bp_dist,
                 uint seed=0, bool run_as_mea=false);
  virtual ~McCaskillModel() { }

  // interface implementations
  virtual void set_constraint(const std::string& str);
  virtual void calculate_posterior(const std::string& seq);
  virtual void prepare_stochastic_traceback(const std::string& seq);
  virtual std::vector<int> stochastic_traceback(const std::string& seq);
  virtual void clean_stochastic_traceback(const std::string& seq);

private:
  bool canonical_only_;
  int bk_st_back_;
  std::string str_;
};
#endif

class AveragedModel : public CentroidFold<Aln>
{
public:
  AveragedModel(CentroidFold<std::string>* cf, uint max_bp_dist, bool run_as_mea=false);
  virtual ~AveragedModel() {}

  virtual void stochastic_fold(const Aln& aln, uint num_samples, std::ostream& out);
  virtual void stochastic_fold(const Aln& aln, uint num_samples, std::vector<BPvecPtr>& bpv);
  virtual void max_mcc_fold(const std::string& name, const Aln& aln, std::ostream& out, uint num_samples);
  
  virtual void compute_expected_accuracy_sampling(const std::string& paren, const Aln& aln,
                                                  uint num_ea_samples, double& sen, double& ppv, double& mcc);

  // interface implementations
  virtual void set_constraint(const std::string& str);
  virtual void calculate_posterior(const Aln& aln);

private:
  CentroidFold<std::string>* cf_;
  std::vector<uint> paren_;
};

#ifdef HAVE_LIBRNA
class AliFoldModel : public CentroidFold<Aln>
{
public:
  AliFoldModel(bool canonical_only, uint max_bp_dist, uint seed=0, bool run_as_mea=false);
  virtual ~AliFoldModel() { }

  // interface implementations
  virtual void set_constraint(const std::string& str);
  virtual void calculate_posterior(const Aln& aln);
  virtual void prepare_stochastic_traceback(const Aln& aln);
  virtual std::vector<int> stochastic_traceback(const Aln& aln);
  virtual void clean_stochastic_traceback(const Aln& aln);

private:
  bool canonical_only_;
  int bk_st_back_;
  std::string str_;
};
#endif

class CONTRAfoldMultiModel : public CentroidFold<Aln>
{
public:
  CONTRAfoldMultiModel(const std::string& model, bool canonical_only, uint max_bp_dist,
                       uint seed=0, bool run_as_mea=false);
  virtual ~CONTRAfoldMultiModel();

  // interface implementations
  //void set_constraint(const std::string& str);
  virtual void calculate_posterior(const Aln& aln);
  virtual void prepare_stochastic_traceback(const Aln& aln);
  virtual std::vector<int> stochastic_traceback(const Aln& aln);
  virtual void clean_stochastic_traceback(const Aln& aln) { }

private:
  CONTRAfoldModel* contrafold_;
  CONTRAfoldM<float>* contrafoldm_;
  AveragedModel* avg_;
};

class AuxModel : public CentroidFold<std::string>
{
public:
  AuxModel(const std::vector<std::string>& bpfiles, bool run_as_mea=false);
  virtual ~AuxModel() { }

  // interface implementations
  virtual void calculate_posterior(const std::string& seq);

private:
  std::vector<std::string> bpfiles_;
  uint pos_;
};

#endif	// #ifndef __INC_CENTROID_FOLD_H__
