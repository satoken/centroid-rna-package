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

#ifndef __INC_ENGINE_AVERAGED_H__
#define __INC_ENGINE_AVERAGED_H__

#include "../folding_engine.h"

class AveragedModel : public FoldingEngine<Aln>
{
public:
  AveragedModel(FoldingEngine<std::string>* cf, uint max_bp_dist, bool run_as_mea=false);
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
  FoldingEngine<std::string>* cf_;
  std::vector<uint> paren_;
};

#endif  //  __INC_ENGINE_AVERAGED_H__
