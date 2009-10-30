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
#include "bp.h"
#include "aln.h"
#if 0
#include <contrafold/contrafold.h>
#else
template < class RealT >
class CONTRAfold;
#endif

typedef unsigned int uint;

class CentroidFold
{
public:
  typedef SCFG::BP::Table<float> BPTable;

  enum {
    AUX,
    PFFOLD,
    CONTRAFOLD,
    ALIPFFOLD
  };
  
public:
  CentroidFold(uint engine,
               bool run_as_mea=false,
	       int num_ea_samples=-1,
	       uint reserved_size=0,
               uint seed=0);
  ~CentroidFold();

  void set_options_for_contrafold(const std::string& model, bool canonical_only, uint max_bp_dist,
                                  uint dist_type=0);

#ifdef HAVE_LIBRNA
  void set_options_for_pf_fold(bool canonical_only, uint max_dist);
#endif

  void calculate_posterior(const std::string& seq);
  void calculate_posterior(const std::string& seq, const std::string& str);
  void calculate_posterior(const std::string& seq, const BPTable& bp);
  void calculate_posterior(const std::list<std::string>& seq);
  void calculate_posterior(const std::vector<std::string>& seq);
  void calculate_posterior(const std::list<std::string>& seq, const std::string& str);
  void calculate_posterior(const std::vector<std::string>& seq, const std::string& str);
  void calculate_posterior(const std::list<std::string>& seq,
                           const std::list<boost::shared_ptr<BPTable> >& bps);

  float decode_structure(float gamma, std::string& paren) const;
  std::pair<std::string,float> decode_structure(float gamma) const;

  void print(std::ostream& out, const std::string& name, const std::string& seq,
	     const std::vector<float>& gamma) const;
  void print(std::ostream& out, const Aln& aln, const std::vector<float>& gamma) const;
  void print_posterior(std::ostream& out, const std::string& seq, float th) const;
  std::string posterior(const std::string& seq, float th) const;

  void stochastic_fold(const std::string& name, const std::string& seq,
                       uint num_samples, uint max_clusters,
                       const std::vector<float>& gamma, std::ostream& out,
                       const std::string& p_outname, float th);
  void stochastic_fold(const Aln& aln,
                       uint num_samples, uint max_clusters,
                       const std::vector<float>& gamma, std::ostream& out,
                       const std::string& p_outname, float th);

  const BPTable& get_bp() const { return bp_; }

  // added by M. Hamada
  void max_mcc_fold(const std::string& name, const std::string& seq, uint num_samples, std::ostream& out);
  void max_mcc_fold(const Aln& aln, uint num_samples, std::ostream& out);
  static void compute_expected_accuracy (const std::string& paren, 
					 const BPTable& bp,
					 double& sen, double& ppv, double& mcc);
  void compute_expected_accuracy_sampling (const std::string& paren, const std::string& seq,
					   uint num_samples,
					   double& sen, double& ppv, double& mcc) const;
  void compute_expected_accuracy_sampling (const std::string& paren, const Aln& aln,
					   uint num_samples,
					   double& sen, double& ppv, double& mcc) const;

  void ps_plot(const std::string& name, const std::string& seq, float g, bool color=true) const;
#ifdef HAVE_LIBRNA
  void svg_plot(const std::string& name, const std::string& seq, float g) const;
#endif

private:
  uint engine_;
  bool mea_;
  int num_ea_samples_; // for computing expected accuracies
  BPTable bp_;
  bool canonical_only_;
  mutable CONTRAfold<float>* contrafold_;
  std::string model_;
  uint max_bp_dist_;
  uint seed_;
  uint dist_type_;
};

#endif	// #ifndef __INC_CENTROID_FOLD_H__
