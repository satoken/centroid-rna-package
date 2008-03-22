// $Id$

#ifndef __INC_RULE_H__
#define __INC_RULE_H__

namespace Rule
{
  // rule type identifier
  enum type_t {
    E = 0,			// X -> epsilon
    N = 1,			// X -> Y
    L = 2,			// X -> a Y
    R = 3,			// X -> Y a
    P = 4,			// X -> a Y b
    B = 5,			// X -> Y Z
    N_TYPE = 6
  };

  // parameter type identifier
  enum { TYPE, BINARY, UNARY_N, UNARY_L, UNARY_R, UNARY_P, EMIT_S, EMIT_P };

  typedef unsigned short symbol_t;
  const symbol_t START = 0;
};

#endif	// __INC_RULE_H__

// Local Variables:
// mode: C++
// End:

