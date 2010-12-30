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

#ifndef __INC_BP_H__
#define __INC_BP_H__

#include <iosfwd>
#include <list>
#include "cyktable.h"

template < class T >
class BPTableTmpl
{
public:
  typedef T value_type;
      
public:
  BPTableTmpl() : bp_(), size_(0), reserved_size_(0) { }
      
  BPTableTmpl(uint sz, uint max_dist=0) : bp_(sz, max_dist), size_(sz), reserved_size_(sz)
  {
    bp_.fill(0.0);
  }

  BPTableTmpl(const BPTableTmpl& x) : bp_(x.bp_), size_(x.size_), reserved_size_(x.reserved_size_) { } 

  void reserve(uint sz, uint max_dist=0)
  {
    bp_.resize(sz, max_dist);
    bp_.fill(0);
    reserved_size_=sz;
  }

  uint reserved_size() const { return reserved_size_; }

  void resize(uint size, uint max_dist=0)
  {
    if (size>reserved_size() || max_dist!=this->max_dist()) {
      reserve(size, max_dist);
    }
    size_ = size;
    bp_.fill(0);
  }

  uint size() const { return size_; }

  uint max_dist() const { return bp_.max_dist(); };

  void update(uint i, uint j, T v) { bp_(i,j)=v; }

  void add(uint i, uint j, T v) { bp_(i,j)+=v; }

  T operator()(uint i, uint j) const { return bp_(i,j); }

  template < class Seq, class RuleSet >
  void parse(const Seq& seq, const RuleSet& rules);

  bool parse(const std::string& str, bool ignore_alone_pair=false, uint minloop=3);

  bool load(const char* filename);
  bool save(const char* filename, const std::string& seq, float th) const;
  bool save(std::ostream& out, const std::string& seq, float th) const;

  std::vector<T> calc_nonbp_prob() const;

private:
  CYKTable<T> bp_;
  uint size_;
  uint reserved_size_;
};

#endif	// __INC_BP_H__

// Local Variables:
// mode: C++
// End:
