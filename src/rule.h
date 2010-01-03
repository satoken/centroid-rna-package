/*
 * $Id$
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

