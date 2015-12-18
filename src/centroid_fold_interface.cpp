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

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include "centroid_fold_interface.h"
#include "folding_engine.h"
#include "engine/contrafold.h"
#include "engine/mccaskill.h"
#include "engine/pfold.h"
#include "engine/mixture.h"
#include "engine/aux.h"
//#include "engine/contrafoldm.h"
#include "engine/alifold.h"
#include "engine/averaged.h"
#include "centroid.h"
#include "mea.h"

CentroidFold::
CentroidFold(int engine,
             bool run_as_mea /*=false*/,
	     int num_ea_samples /*=-1*/,
             int reserved_size /*=0*/,
             int seed /*=0*/)
  : engine_(engine),
    mea_(run_as_mea),
    num_ea_samples_(num_ea_samples),
    canonical_only_(true),
    model_(),
    max_bp_dist_(0),
    seed_(seed),
    dist_type_(0),
    cf_single_(NULL),
    cf_multiple_(NULL),
    cf_list_(),
    bp_(NULL)
{
  switch (engine_)
  {
    default:
#ifdef HAVE_LIBRNA
    case PFFOLD:
      cf_single_ = new McCaskillModel(canonical_only_, max_bp_dist_, NULL, seed_, mea_);
      cf_list_.push_back(new AveragedModel(cf_single_, max_bp_dist_, mea_));
      break;
    case ALIPFFOLD:
      cf_single_ = new McCaskillModel(canonical_only_, max_bp_dist_, NULL, seed_, mea_);
      cf_list_.push_back(new AliFoldModel(canonical_only_, max_bp_dist_, NULL, seed_, mea_));
      break;
#endif  // HAVE_LIBRNA
    case CONTRAFOLD:
      cf_single_ = new CONTRAfoldModel(canonical_only_, max_bp_dist_, "", seed_, mea_);
      cf_list_.push_back(new AveragedModel(cf_single_, max_bp_dist_, mea_));
      break;
  }

  if (cf_list_.size()==1)
  {
    cf_multiple_ = cf_list_[0];
    cf_list_.resize(0);
  }
  else 
  {
    std::vector<std::pair<FoldingEngine<Aln>*,float> > models;
    for (uint i=0; i!=cf_list_.size(); ++i)
      models.push_back(std::make_pair(cf_list_[i], 1.0));
    cf_multiple_ = new MixtureModel<Aln>(models, mea_);
  }
}

CentroidFold::
~CentroidFold()
{
  if (cf_single_) delete cf_single_;
  if (cf_multiple_) delete cf_multiple_;
  for (uint i=0; i!=cf_list_.size(); ++i)
    if (cf_list_[i]) delete cf_list_[i];
}


#if 0
#ifdef HAVE_LIBRNA
void
CentroidFold::
set_options_for_pf_fold(bool canonical_only, uint max_dist)
{
  canonical_only_ = canonical_only;
  if (!canonical_only_)
    Vienna::nonstandards = const_cast<char*>("AAACAGCACCCUGAGGUCUU");
  max_bp_dist_ = max_dist;
}
#endif

void
CentroidFold::
set_options_for_contrafold(const std::string& model, bool canonical_only, uint max_bp_dist, uint dist_type/*=0*/)
{
  model_ = model;
  canonical_only_ = canonical_only;
  max_bp_dist_ = max_bp_dist;
  dist_type_ = dist_type;
  if (!model_.empty()) contrafold_->SetParameters(model_);
}
#endif

void
CentroidFold::
calculate_posterior(const std::string& seq)
{
  std::string seq2(seq);
  std::replace(seq2.begin(), seq2.end(), 't', 'u');
  std::replace(seq2.begin(), seq2.end(), 'T', 'U');
  cf_single_->calculate_posterior(seq2);
  bp_=&cf_single_->get_bp();
}

void
CentroidFold::
calculate_posterior(const std::string& seq, const std::string& str)
{
  std::string seq2(seq);
  std::replace(seq2.begin(), seq2.end(), 't', 'u');
  std::replace(seq2.begin(), seq2.end(), 'T', 'U');
  cf_single_->set_constraint(str);
  cf_single_->calculate_posterior(seq2);
  bp_=&cf_single_->get_bp();
}

void
CentroidFold::
calculate_posterior(const std::list<std::string>& seq)
{
  std::list<std::string> name;
  std::list<std::string> seq2;
  std::list<std::string>::const_iterator x;
  for (x=seq.begin(); x!=seq.end(); ++x) {
    name.push_back("seq");
    seq2.push_back(*x);
    std::replace(seq2.back().begin(), seq2.back().end(), 't', 'u');
    std::replace(seq2.back().begin(), seq2.back().end(), 'T', 'U');
  }
  Aln aln(name, seq2);
  cf_multiple_->calculate_posterior(aln);
  bp_=&cf_multiple_->get_bp();
}

void
CentroidFold::
calculate_posterior(const std::vector<std::string>& seq)
{
  std::list<std::string> seq2(seq.size());
  std::copy(seq.begin(), seq.end(), seq2.begin());
  calculate_posterior(seq2);
}

void
CentroidFold::
calculate_posterior(const std::list<std::string>& seq, const std::string& str)
{
  std::list<std::string> name;
  std::list<std::string> seq2;
  std::list<std::string>::const_iterator x;
  for (x=seq.begin(); x!=seq.end(); ++x) {
    name.push_back("seq");
    seq2.push_back(*x);
    std::replace(seq2.back().begin(), seq2.back().end(), 't', 'u');
    std::replace(seq2.back().begin(), seq2.back().end(), 'T', 'U');
  }
  cf_multiple_->set_constraint(str);
  Aln aln(name, seq2);
  cf_multiple_->calculate_posterior(aln);
  bp_=&cf_multiple_->get_bp();
}

void
CentroidFold::
calculate_posterior(const std::vector<std::string>& seq, const std::string& str)
{
  std::list<std::string> seq2(seq.size());
  std::copy(seq.begin(), seq.end(), seq2.begin());
  calculate_posterior(seq2, str);
}

float
CentroidFold::
decode_structure(float gamma, std::string& paren) const
{
  float p=0.0;
  paren.resize(bp_->size());
  std::fill(paren.begin(), paren.end(), '.');
  if (!mea_) {
    if (max_bp_dist_==0)
      p = Centroid::execute(*bp_, paren, gamma);
    else
      p = Centroid::execute(*bp_, paren, max_bp_dist_, gamma);
  } else {
    if (max_bp_dist_==0)
      p = MEA::execute(*bp_, paren, gamma);
    else
      p = MEA::execute(*bp_, paren, max_bp_dist_, gamma);
  }
  return p;
}

std::pair<std::string,float>
CentroidFold::
decode_structure(float gamma) const
{
  std::string paren;
  float p = decode_structure(gamma, paren);
  return std::make_pair(paren, p);
}

void
CentroidFold::
ps_plot(const std::string& name, const std::string& seq, float g, bool color) const
{
  cf_single_->ps_plot(name, seq, g, color);
}

#ifdef HAVE_LIBRNA
void
CentroidFold::
svg_plot(const std::string& name, const std::string& seq, float g) const
{
  cf_single_->svg_plot(name, seq, g);
}
#endif

