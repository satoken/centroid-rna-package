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

#include <vector>
#include <sys/time.h>
#include "../contrafold/contrafold.h"
#include "engine/contrafoldm.h"

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

