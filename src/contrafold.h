/*
 * wrapper routines for CONTRAfold
 *
 * Copyright (C) 2008 Kengo Sato
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

#ifndef __INC_CONTRAFOLD_H__
#define __INC_CONTRAFOLD_H__

#include <string>
#include <vector>

namespace CONTRAfold {

  template < class T >
  struct WSImpl;

  template < class T >
  struct WS {
    WSImpl<T>* impl;

    WS(unsigned int size);
    ~WS();
  };
  
  template < class T >
  const T* ComputePosterior(const std::string& seq, bool canonical_only = true);

  template < class T >
  const T* ComputePosterior(const std::string& seq, WS<T>& ws, bool canonical_only = true);
};

#endif	// __INC_CONTRAFOLD_H__

// Local Variables:
// mode: C++
// End:
