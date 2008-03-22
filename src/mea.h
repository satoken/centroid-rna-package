// $Id$

#ifndef __INC_MEA_H__
#define __INC_MEA_H__

#include <string>
//#include "dptable.h"
//#include "inside.h"
//#include "outside.h"
#include "bp.h"

namespace SCFG
{
  namespace MEA
  {
    template < class V >
    V execute(const SCFG::BP::Table<V>& bp, V gamma, std::string& paren);

    template < class V, class Seq, class RuleSet >
    V execute(const Seq& seq, const RuleSet& rules, V gamma, std::string& paren)
    {
      SCFG::BP::Table<double> bp;
      bp.parse(seq, rules);
      return SCFG::MEA::execute(bp, gamma, paren);
    }
  }
};

#endif	// __INC_CYK_H__

// Local Variables:
// mode: C++
// End:
