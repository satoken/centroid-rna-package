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

#include <iostream>
#include <cassert>
#include "rule.h"
#include "centroid.h"
#include "cyktable.h"
//#include "dptable.h"

namespace Centroid
{
  template < class T >
  struct Estimator
  {
    Estimator(const T& t, std::string& paren, float gamma)
      : t_(t), paren_(paren), th_(1.0/(gamma+1.0)), gamma_(gamma), ea_(0.0) { }

    void operator()(uint i, uint j)
    {
      if (i!=j && t_(i,j-1)>th_) {
        paren_[i] = '(';
        paren_[j-1] = ')';
        ea_ += (gamma_+1.0)*2 * t_(i,j-1) - 2;
      }
    }

    float expected_accuracy() const { return ea_; }

    const T& t_;
    std::string& paren_;
    float th_;
    float gamma_;
    float ea_;
  };

  template < class V >
  struct Cell
  {
    typedef V value_type;

    value_type val;
    Rule::type_t type;
    uint br_pos;

    Cell()
      : val(static_cast<V>(0.0)),
        type(static_cast<Rule::type_t>(-1)),
        br_pos(static_cast<uint>(-1))
    {}

    bool
    update(value_type val, Rule::type_t type,
           uint br_pos=static_cast<uint>(-1))
    {
      if (val>this->val) {
        this->val = val;
        this->type = type;
        this->br_pos = br_pos;
        return true;
      }
      return false;
    }
  };

  template < class BPTable, class DPTable >
  struct Updater
  {
    typedef typename BPTable::value_type value_type;

    Updater(const BPTable& bp, DPTable& dp, float gamma, uint min_loop=3)
      : bp_(bp), dp_(dp), gamma_(gamma), min_loop_(min_loop)
    {
    }

    void operator()(uint i, uint j)
    {
      dp_(i,j).val = static_cast<value_type>(-1e100);
      // rule type E: X -> epsilon
      if (i==j) {
        dp_(i,j).update(0.0, Rule::E);
      } else {
        // rule type L: X -> a Y
        if (i+1<=j) {
          value_type v = dp_(i+1,j).val;
          dp_(i,j).update(v, Rule::L);
        }
        // rule type R: X -> Y a
        if (i<=j-1) {
          value_type v = dp_(i,j-1).val;
          dp_(i,j).update(v, Rule::R);
        }
        // rule type P: X -> a Y b
        if (min_loop_+i+1<=j-1) {
          value_type v = dp_(i+1,j-1).val + (gamma_+1.0)*2 * bp_(i,j-1) - 2;
          dp_(i,j).update(v, Rule::P);
        }
        // rule type B: X -> Y Z
        for (uint k=i+1; k<j; ++k) {
          value_type v = dp_(i,k).val + dp_(k,j).val;
          dp_(i,j).update(v, Rule::B, k);
        }
      }
      assert(dp_(i,j).type!=static_cast<Rule::type_t>(-1));
    }

    const BPTable& bp_;
    DPTable& dp_;
    float gamma_;
    uint min_loop_;
  };

  template < class DPTable >
  struct TraceBack
  {
    TraceBack(std::string& paren, const DPTable& dp)
      : paren_(paren), dp_(dp)
    {}

    void
    operator()(uint i, uint j)
    {
      assert(i!=static_cast<uint>(-1));
      assert(j!=static_cast<uint>(-1));

      switch (dp_(i,j).type) {
	case Rule::E:
	  assert(i==j);
	  return;
	  break;
	case Rule::L:
	  (*this)(i+1, j);
	  break;
	case Rule::R:
	  (*this)(i, j-1);
	  break;
	case Rule::P:
	  paren_[i]='('; paren_[j-1]=')';
	  (*this)(i+1, j-1);
	  break;
	case Rule::B:
	  assert(dp_(i,j).br_pos!=static_cast<uint>(-1));
	  (*this)(i, dp_(i,j).br_pos);
	  (*this)(dp_(i,j).br_pos, j);
	  break;
	default:
	  assert(!"unreachable");
	  break;
      }
    }

    std::string& paren_;
    const DPTable& dp_;
  };

  template < class T >
  float
  execute(const T& table, std::string& paren, float gamma)
  {
    if (gamma<1.0) {
      Estimator<T> est(table, paren, gamma);
      inside_traverse(0, table.size(), est);
      return est.expected_accuracy();
    } else {
      typedef typename T::value_type value_type;
      typedef CYKTable< Cell<value_type> > DPTable;

      DPTable dp(table.size()+1);
      Updater<T,DPTable> update(table, dp, gamma);
      inside_traverse(0, table.size(), update);

      TraceBack<DPTable> traceback(paren, dp);
      traceback(0, dp.size()-1);

      return dp(0, dp.size()-1).val;
    }
  }

  template < class T >
  float
  execute(const T& table, std::string& paren, uint max_dist, float gamma)
  {
    if (gamma<1.0) {
      Estimator<T> est(table, paren, gamma);
      inside_traverse(0, table.size(), max_dist, est);
      return est.expected_accuracy();
    } else {
      typedef typename T::value_type value_type;
      typedef CYKTable< Cell<value_type> > DPTable;

      DPTable dp(table.size()+1, max_dist);
      std::vector<Cell<value_type> > outer(table.size()+1);
      Updater<T,DPTable> update(table, dp, gamma);
      inside_traverse(0, table.size(), max_dist, update);
      outer[0].val = -1e100;
      outer[0].update(0.0, Rule::E);
      for (uint j=1; j!=outer.size(); ++j) {
        outer[j].val = -1e100;
        outer[j].update(outer[j-1].val, Rule::R);
        uint l = j>max_dist ? j-max_dist : 0;
        for (uint k=l; k<j; ++k) {
          if (dp(k,j).type==Rule::P)
            outer[j].update(outer[k].val+dp(k,j).val, Rule::B, k);
        }
      }

      int j=outer.size()-1;
      while (j>0) {
        if (outer[j].type==Rule::B) {
          TraceBack<DPTable> traceback(paren, dp);
          traceback(outer[j].br_pos, j);
          j = outer[j].br_pos;
        } else {
          j--;
        }
      }

      return outer.back().val;
    }
  }
};

// instantiation
#include "cyktable.h"
#include "bp.h"

template
float
Centroid::
execute(const BPTableTmpl<float>& table, std::string& paren, float gamma);

template
float
Centroid::
execute(const CYKTable<float>& table, std::string& paren, float gamma);

template
float
Centroid::
execute(const BPTableTmpl<float>& table, std::string& paren, uint max_dist, float gamma);

