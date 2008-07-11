/*
 * $Id$
 *
 * Copyright (C) 2007,2008 Kengo Sato
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


#ifndef __INC_CYKTABLE_H__
#define __INC_CYKTABLE_H__

#include <cassert>
#include <vector>
#include <utility>

typedef unsigned int uint;

namespace SCFG
{
  typedef std::pair<uint,uint> Pos;

  /**
   * @brief A templete class of CYK tables
   *        used for Cocke-Younger-Kasami (CYK) algorithm
   */
  template <class T>
  class CYKTable
  {
  public:
    typedef T value_type;

  private:
    std::vector<T> table_;
    std::vector<T*> ptr_;
    uint i_from_;
    uint i_to_;
    uint j_from_;
    uint j_to_;
  
  public:
    CYKTable()
      : table_(), ptr_(), j_to_(0)
    {
    }

    CYKTable(uint i_from, uint i_to, uint j_from, uint j_to)
      : table_(), ptr_(), j_to_(0)
    {
      setup(i_from, i_to, j_from, j_to);
    }

    CYKTable(uint size)
      : table_(), ptr_(), j_to_(0)
    {
      setup(0, size, 0, size);
    }

    CYKTable(uint size, const T& val)
      : table_(), ptr_(), j_to_(0)
    {
      setup(0, size, 0, size);
      fill(val);
    }
    
    uint estimate_size(uint i_from, uint i_to, uint j_from, uint j_to)
    {
      uint ret=0;
      for (uint j=j_from; j!=j_to; ++j) {
	if (j==i_from) {
	  ret++;
	  continue;
	}
#if 0
	for (uint i=std::min(j,i_to-1); ; --i) {
	  ret++;
	  if (i==i_from) break;
	}
#else
	ret += std::min(j,i_to-1)-i_from+1;
#endif
      }
      return ret;
    }

    void setup(uint i_from, uint i_to, uint j_from, uint j_to)
    {
      assert(i_from<=i_to);
      assert(j_from<=j_to);
      assert(i_from<=j_from);
      assert(i_to<=j_to);

      uint new_size = estimate_size(i_from, i_to, j_from, j_to);
      table_.resize(new_size);
      ptr_.resize(j_to);
      std::fill(ptr_.begin(), ptr_.end(), static_cast<T*>(NULL));

      uint x=0;
      for (uint j=j_from; j!=j_to; ++j) {
	ptr_[j] = &table_[x - i_from];
	if (j==i_from) {
	  x++;
	  continue;
	}
	for (uint i=std::min(j,i_to-1); ; --i) {
	  x++;
	  if (i==i_from) break;
	}
      }

      i_from_ = i_from;
      i_to_ = i_to;
      j_from_ = j_from;
      j_to_ = j_to;
    }

    uint size() const
    {
      return std::max(i_to_-i_from_, j_to_-j_from_);
    }

    void resize(uint i_from, uint i_to, uint j_from, uint j_to)
    {
      setup(i_from, i_to, j_from, j_to);
    }

    void resize(uint size)
    {
      resize(0, size, 0, size);
    }

    void fill(const T& val)
    {
      std::fill(table_.begin(), table_.end(), val);
    }

    uint table_size() const { return table_.size(); }

    inline
    void put(uint i, uint j, const T& val)
    {
      assert(i>=i_from_);
      assert(i<i_to_);
      assert(j>=j_from_);
      assert(j<j_to_);
      assert(i<=j);
      assert(ptr_[j]!=NULL);
      ptr_[j][i]=val;
    }

    inline
    void put(const Pos& pos, const T& val)
    {
      put(pos.first, pos.second, val);
    }

    inline
    const T& get(uint i, uint j) const
    {
      assert(i>=i_from_);
      assert(i<=i_to_);
      assert(j>=j_from_);
      assert(j<=j_to_);
      assert(i<=j);
      assert(ptr_[j]!=NULL);
      return ptr_[j][i];
    }

    inline
    const T& get(const Pos& pos) const
    {
      return get(pos.first, pos.second);
    }

    inline
    T& get(uint i, uint j)
    {
      assert(i>=i_from_);
      assert(i<=i_to_);
      assert(j>=j_from_);
      assert(j<=j_to_);
      assert(i<=j);
      assert(ptr_[j]!=NULL);
      return ptr_[j][i];
    }

    inline
    T& get(const Pos& pos)
    {
      return get(pos.first, pos.second);
    }


    inline
    const T& operator()(uint i, uint j) const
    {
      return get(i,j);
    }

    inline
    const T& operator()(const Pos& pos) const
    {
      return get(pos);    
    }

    inline
    T& operator()(uint i, uint j)
    {
      return get(i,j);
    }

    inline
    T& operator()(const Pos& pos)
    {
      return get(pos);    
    }
  };

  template <class Update>
  void
  inside_traverse(uint from, uint to, Update& update)
  {
    // for each position
    for (uint j=from; j!=to+1; ++j) {
      update(j, j);
      if (j==from) continue;
    
      // for each substring
      for (uint i=j-1; ; --i) {
	update(i, j);
	if (i==from) break;
      }
    }
  }

  template <class Update>
  void
  outside_traverse(uint from, uint to, Update& update)
  {
    // for each position
    for (uint j=to; ; --j) {
      // for each substring
      for (uint i=from; i!=j; ++i) {
	update(i, j);
      }
      update(j, j);
      if (j==from) break;
    }
  }
};

#endif // __INC_CYKTABLE_H__

// Local Variables:
// mode: C++
// End:
