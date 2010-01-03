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

#ifndef __INC_IUPAC_H__
#define __INC_IUPAC_H__

#include <vector>
#include <string>

class IUPACsymbol;

typedef std::vector<IUPACsymbol> IUPACsequence;

enum {
  RNA_A = 0,
  RNA_C = 1,
  RNA_G = 2,
  RNA_T = 3,
  RNA_U = RNA_T,
  N_RNA = 4,
  N_IUPAC = 16
};

class IUPACsymbol
{
  enum {
    GAP = 0,
    A = 1<<RNA_A,
    C = 1<<RNA_C,
    G = 1<<RNA_G,
    T = 1<<RNA_T,
    U = T,
    R = A | G,
    Y = C | T,
    M = A | C,
    K = G | T,
    S = C | G,
    W = A | T,
    B = C | G | T,
    D = A | G | T,
    H = A | C | T,
    V = A | C | G,
    N = A | C | G | T
  };

public:
  explicit IUPACsymbol(char c = '-');
  IUPACsymbol(const IUPACsymbol& s) : iupac_(s.iupac_) {}
  IUPACsymbol& operator=(const IUPACsymbol& s)
  {
    if (this!=&s) iupac_ = s.iupac_;
    return *this;
  }
  unsigned int index() const { return static_cast<unsigned int>(iupac_); }
  char to_char() const;

  static void transform(IUPACsequence& r, const std::string& s);

private:
  unsigned char iupac_;
};

#endif	// __INC_IUPAC_H__

// Local Variables:
// mode: C++
// End:
