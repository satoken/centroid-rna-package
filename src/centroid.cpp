// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <iostream>
#include <cassert>
#include "rule.h"
#include "centroid.h"
#include "cyktable.h"
//#include "dptable.h"

namespace SCFG
{
  namespace Centroid
  {
    template < class T >
    struct Estimator
    {
      Estimator(const T& t, std::string& paren, double n, double gamma)
	: t_(t), paren_(paren), th_(n/(gamma+1.0)), gamma_(gamma), n_(n), ea_(0.0) { }

      void operator()(uint i, uint j)
      {
	if (i!=j && t_(i,j-1)>th_) {
	  paren_[i] = '(';
	  paren_[j-1] = ')';
	  ea_ += (gamma_+1.0)/n_*2 * t_(i,j-1) - 2;
	}
      }

      double expected_accuracy() const { return ea_; }

      const T& t_;
      std::string& paren_;
      double th_;
      double gamma_;
      double n_;
      double ea_;
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

      Updater(const BPTable& bp, DPTable& dp, double n, double gamma,
	      uint min_loop=3)
	: bp_(bp), dp_(dp), n_(n), gamma_(gamma), min_loop_(min_loop)
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
	    value_type v = dp_(i+1,j-1).val + (gamma_+1.0)/n_*2 * bp_(i,j-1) - 2;
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
      double n_;
      double gamma_;
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
	  assert(!"never come here.");
	  break;
	}
      }

      std::string& paren_;
      const DPTable& dp_;
    };

    template < class T >
    double
    execute(const T& table, std::string& paren, double n /*=1.0*/, double gamma /*=1.0*/)
    {
      if (gamma<=1.0) {
	Estimator<T> est(table, paren, n, gamma);
	SCFG::inside_traverse(0, table.size(), est);
	return est.expected_accuracy();
      } else {
	typedef typename T::value_type value_type;
	typedef CYKTable< Cell<value_type> > DPTable;

	DPTable dp(table.size()+1);
	Updater<T,DPTable> update(table, dp, n, gamma);
	SCFG::inside_traverse(0, table.size(), update);

	TraceBack<DPTable> traceback(paren, dp);
	traceback(0, dp.size()-1);

	return dp(0, dp.size()-1).val;
      }
    }
  };
};

// instantiation
#include "cyktable.h"
#include "bp.h"

template
double
SCFG::Centroid::
execute(const CYKTable<uint>& table, std::string& paren, double n, double gamma);

template
double
SCFG::Centroid::
execute(const SCFG::BP::Table<double>& table, std::string& paren, double n, double gamma);

