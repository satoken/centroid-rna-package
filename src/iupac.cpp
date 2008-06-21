// $Id$

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
