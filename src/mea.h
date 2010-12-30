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

#ifndef __INC_MEA_H__
#define __INC_MEA_H__

#include <string>
//#include "dptable.h"
//#include "inside.h"
//#include "outside.h"
#include "bp.h"

namespace MEA
{
  template < class V >
  V execute(const BPTableTmpl<V>& bp, std::string& paren, float gamma);

  template < class V >
  V execute(const BPTableTmpl<V>& bp, std::string& paren, uint max_dist, float gamma);

  template < class V, class Seq, class RuleSet >
  V execute(const Seq& seq, const RuleSet& rules, std::string& paren, float gamma)
  {
    BPTableTmpl<double> bp;
    bp.parse(seq, rules);
    return MEA::execute(bp, paren, gamma);
  }
}

#endif	// __INC_CYK_H__

// Local Variables:
// mode: C++
// End:
