// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <cassert>
#include "mea.h"
#include "rule.h"
//#include "dptable.h"
#include "cyktable.h"

namespace SCFG
{
  namespace MEA
  {
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

      Updater(const BPTable& bp, DPTable& dp,
	      value_type gamma=1.0, uint min_loop=3)
	: bp_(bp), dp_(dp), gamma_(gamma), min_loop_(min_loop)
      {
      }

      void operator()(uint i, uint j)
      {
	dp_(i,j).val = -1e100;
	// rule type E: X -> epsilon
	if (i==j) {
	  dp_(i,j).update(0.0, Rule::E);
	} else {
	  // rule type L: X -> a Y
	  if (i+1<=j) {
	    value_type v = dp_(i+1,j).val + bp_[i];
	    dp_(i,j).update(v, Rule::L);
	  }
	  // rule type R: X -> Y a
	  if (i<=j-1) {
	    value_type v = dp_(i,j-1).val + bp_[j-1];
	    dp_(i,j).update(v, Rule::R);
	  }
	  // rule type P: X -> a Y b
	  if (min_loop_+i+1<=j-1) {
	    value_type v = dp_(i+1,j-1).val + 2*gamma_*bp_(i,j-1);
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
      value_type gamma_;
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

    template < class V >
    V
    execute(const SCFG::BP::Table<V>& bp, V gamma, std::string& paren)
    {
      typedef V value_type;
      typedef SCFG::BP::Table<V> BPTable;
      typedef CYKTable< Cell<value_type> > DPTable;

      DPTable dp(bp.size()+1);
      Updater<BPTable,DPTable> update(bp, dp, gamma);
      SCFG::inside_traverse(0, bp.size(), update);

      paren.resize(dp.size()-1);
      std::fill(paren.begin(), paren.end(), '.');
      TraceBack<DPTable> traceback(paren, dp);
      traceback(0, dp.size()-1);

      return dp(0, dp.size()-1).val;
    }
  };
};

//#include "rna.h"
//#include "log_value.h"
#include "bp.h"

// instantiate
template
double
SCFG::MEA::execute(const SCFG::BP::Table<double>& bp, double gamma,
		   std::string& paren);
