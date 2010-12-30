/*
 * $Id$
 *
 * Copyright (C) 2007-2010 Kengo Sato
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
  uint size_;
  uint max_dist_;
  
public:
  CYKTable()
    : table_(), ptr_(), size_(0), max_dist_(0)
  {
  }

  CYKTable(uint size, uint max_dist=0)
    : table_(), ptr_(), size_(size), max_dist_(std::min(max_dist,size))
  {
    setup();
  }

  CYKTable(const CYKTable& x)
    : table_(), ptr_(), size_(x.size_), max_dist_(x.max_dist_)
  {
    setup();
    std::copy(x.table_.begin(), x.table_.end(), table_.begin());
  }
    
  uint estimate_size(uint size, uint max_dist)
  {
    if (max_dist==0)
      return size_*(size_+1)/2;
    else
      return size_*max_dist;
  }

  void setup()
  {
    uint new_size = estimate_size(size_, max_dist_);
    table_.resize(new_size);
    ptr_.resize(size_);
    std::fill(ptr_.begin(), ptr_.end(), static_cast<T*>(NULL));

    if (max_dist_==0) {
      for (uint i=0; i!=size_; ++i)
        ptr_[i] = &table_[size_*(size_+1)/2 - (size_-i)*(size_-i+1)/2 - i];
    } else {
      for (uint i=0; i!=size_; ++i)
        ptr_[i] = &table_[i*max_dist_-i];
    }
  }

  uint size() const { return size_; }

  uint table_size() const { return table_.size(); }

  uint max_dist() const { return max_dist_; }

  void resize(uint size, uint max_dist=0)
  {
    size_ = size;
    max_dist_ = std::min(max_dist, size);
    setup();
  }

  void fill(const T& val)
  {
    std::fill(table_.begin(), table_.end(), val);
  }


  inline
  void put(uint i, uint j, const T& val)
  {
    assert(i<=j);
    assert(max_dist_==0 || j-i<=max_dist_);
    assert(ptr_[j]!=NULL);
    ptr_[i][j]=val;
  }

  inline
  void put(const Pos& pos, const T& val)
  {
    put(pos.first, pos.second, val);
  }

  inline
  const T& get(uint i, uint j) const
  {
    assert(i<=j);
    assert(max_dist_==0 || j-i<=max_dist_);
    assert(ptr_[i]!=NULL);
    return ptr_[i][j];
  }

  inline
  const T& get(const Pos& pos) const
  {
    return get(pos.first, pos.second);
  }

  inline
  T& get(uint i, uint j)
  {
    assert(i<=j);
    assert(ptr_[i]!=NULL);
    return ptr_[i][j];
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
inside_traverse(uint from, uint to, uint width, Update& update)
{
  // for each position
  for (uint j=from; j!=to+1; ++j) {
    update(j, j);
    if (j==from) continue;
    
    // for each substring
    for (uint i=j-1; ; --i) {
      update(i, j);
      if (i==from || j-i>=width-1) break;
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

template <class Update>
void
outside_traverse(uint from, uint to, uint width, Update& update)
{
  // for each position
  for (uint j=to; ; --j) {
    // for each substring
    uint l=std::max(from, j>width ? j-width-1 : 0);
    for (uint i=l; i!=j; ++i) {
      update(i, j);
    }
    update(j, j);
    if (j==from) break;
  }
}

#endif // __INC_CYKTABLE_H__

// Local Variables:
// mode: C++
// End:
