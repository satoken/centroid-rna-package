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

#ifndef __INC_ENGINE_CONTRAFOLDM_H__
#define __INC_ENGINE_CONTRAFOLDM_H__

#include "../folding_engine.h"
#include "engine/contrafold.h"
#include "engine/averaged.h"

class CONTRAfoldMultiModel : public FoldingEngine<Aln>
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

#endif  //  __INC_ENGINE_CONTRAFOLDM_H__
