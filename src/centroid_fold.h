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
#ifdef HAVE_LIBCONTRAFOLD
#if 0
#include <contrafold.h>
#else
template < class RealT >
class CONTRAfold;
#endif
#endif

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
  CentroidFold(unsigned int engine,
               bool run_as_mea=false,
	       unsigned int reserved_size=0,
               unsigned int seed=0);
  ~CentroidFold();

#ifdef HAVE_LIBCONTRAFOLD
  void set_options_for_contrafold(const std::string& model, bool canonical_only, uint max_bp_dist);
#endif

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
  void print_posterior(std::ostream& out, const std::string& seq, float th) const;
  std::string posterior(const std::string& seq, float th) const;

  void stochastic_fold(const std::string& name, const std::string& seq,
                       uint num_samples, uint max_clusters,
                       const std::vector<float>& gamma, std::ostream& out);
  void stochastic_fold(const std::string& name, const std::string& consensus,
                       const std::vector<std::string>& seq,
                       uint num_samples, uint max_clusters,
                       const std::vector<float>& gamma, std::ostream& out);
  void stochastic_fold(const std::string& name, const std::string& consensus,
                       const std::list<std::string>& seq,
                       uint num_samples, uint max_clusters,
                       const std::vector<float>& gamma, std::ostream& out);


#ifdef HAVE_LIBRNA
  void ps_plot(const std::string& name, const std::string& seq, float g) const;
  void svg_plot(const std::string& name, const std::string& seq, float g) const;
#endif

private:
  unsigned int engine_;
  bool mea_;
  BPTable bp_;
  bool canonical_only_;
#ifdef HAVE_LIBCONTRAFOLD
  CONTRAfold<float>* contrafold_;
  std::string model_;
#endif
  uint max_bp_dist_;
};

#endif	// #ifndef __INC_CENTROID_FOLD_H__
