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

#ifndef __INC_ENGINE_AUX_H__
#define __INC_ENGINE_AUX_H__

#include "../folding_engine.h"

class AuxModel : public FoldingEngine<std::string>
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

#endif  //  __INC_ENGINE_AUX_H__
