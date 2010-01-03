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

#include "engine/mixture.h"

template <class SEQ>
MixtureModel<SEQ>::
MixtureModel(const std::vector<std::pair<CentroidFold<SEQ>*,float> >& models,
             bool run_as_mea /*=false*/)
  : CentroidFold<SEQ>(run_as_mea, 0), models_(models)
{
}

struct MulAdd
{
  MulAdd(BPTable& bp, float w, const BPTable& t) : bp_(bp), w_(w), t_(t) {}
  void operator()(uint i, uint j) { bp_.update(i, j, bp_(i,j)+w_*t_(i,j)); }
  
  BPTable& bp_;
  float w_;
  const BPTable& t_;
};

struct Div
{
  Div(BPTable& bp, float d) : bp_(bp), d_(d) {}
  void operator()(uint i, uint j) { bp_.update(i, j, bp_(i,j)/d_); }

  BPTable& bp_;
  float d_;
};

template <class SEQ>
void
MixtureModel<SEQ>::
calculate_posterior(const SEQ& seq)
{
  float sum_w=0.0;
  bp_.resize(seq.size());
  typename std::vector<std::pair<CentroidFold<SEQ>*,float> >::iterator x;
  for (x=models_.begin(); x!=models_.end(); ++x)
  {
    MulAdd ma(bp_, x->second, x->first->get_bp());
    SCFG::inside_traverse(0, bp_.size(), ma);
    sum_w+=x->second;
  }
  Div div(bp_, sum_w);
  SCFG::inside_traverse(0, bp_.size(), div);
}

// instantiation
template class MixtureModel<Aln>;
