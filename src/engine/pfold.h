/*
 * $Id: mccaskill.h 91 2010-01-03 16:46:51Z sato-kengo $
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

#ifndef __INC_ENGINE_PFOLD_H__
#define __INC_ENGINE_PFOLD_H__

#include "../folding_engine.h"

template <class SEQ>
class PfoldModel : public FoldingEngine<SEQ>
{
public:
  PfoldModel(const std::string& pfold_bin_dir,
             const std::string& awk_bin, const std::string& sed_bin, bool run_as_mea=false)
    : FoldingEngine<SEQ>(0, run_as_mea),
      pfold_bin_dir_(pfold_bin_dir), awk_bin_(awk_bin), sed_bin_(sed_bin)
  { }
  virtual ~PfoldModel() { }

  // interface implementations
  virtual void calculate_posterior(const SEQ& seq);

private:
  using FoldingEngine<SEQ>::bp_;
  std::string pfold_bin_dir_;
  std::string awk_bin_;
  std::string sed_bin_;
};

#endif  //  __INC_ENGINE_PFOLD_H__

// Local Variables:
// mode: C++
// End:
