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

#include <sys/time.h>
#include "engine/contrafold.h"
#include "../contrafold/contrafold.h"

CONTRAfoldModel::
CONTRAfoldModel(const std::string& model, bool canonical_only, uint max_bp_dist,
                uint seed /*=0*/, bool run_as_mea /*=false*/)
  : FoldingEngine<std::string>(run_as_mea, max_bp_dist),
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

