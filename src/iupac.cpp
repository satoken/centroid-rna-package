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

#include "iupac.h"
#include <cassert>
//#include <iostream>

IUPACsymbol::
IUPACsymbol(char c /*='-'*/) : iupac_(GAP)
{
  switch (c) {
  case '-': iupac_ = GAP; break;
  case 'a': case 'A': iupac_ = A; break;
  case 'c': case 'C': iupac_ = C; break;
  case 'g': case 'G': iupac_ = G; break;
  case 't': case 'T': iupac_ = T; break;
  case 'u': case 'U': iupac_ = U; break;
  case 'r': case 'R': iupac_ = R; break;
  case 'y': case 'Y': iupac_ = Y; break;
  case 'm': case 'M': iupac_ = M; break;
  case 'k': case 'K': iupac_ = K; break;
  case 's': case 'S': iupac_ = S; break;
  case 'w': case 'W': iupac_ = W; break;
  case 'b': case 'B': iupac_ = B; break;
  case 'd': case 'D': iupac_ = D; break;
  case 'h': case 'H': iupac_ = H; break;
  case 'v': case 'V': iupac_ = V; break;
  case 'n': case 'N': iupac_ = N; break;
  default:
    //std::cout << c << std::endl;
    assert(!"never come here.");
    break;
  }
}

char
IUPACsymbol::
to_char() const
{
  static char char_table[] = "-acgurymkswbdhvn";
  return char_table[index()];
}

// static
void
IUPACsymbol::
transform(IUPACsequence& r, const std::string& s)
{
  r.resize(s.size());
  for (unsigned int i=0; i!=s.size(); ++i)
    r[i] = IUPACsymbol(s[i]);
}
