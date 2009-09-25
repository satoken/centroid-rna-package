/*
 * $Id$
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

#ifndef __INC_BP_H__
#define __INC_BP_H__

#include <iosfwd>
#include <list>
#include "cyktable.h"

namespace SCFG
{
  namespace BP
  {
    template < class T >
    class Table
    {
    public:
      typedef T value_type;
      
    public:
      Table() : bp_(), q_(), size_(0) { }
      
      Table(uint sz, uint max_dist=0) : bp_(sz), q_(sz), size_(sz)
      {
	reserve(sz, max_dist);
      }

      Table(const Table& x) : bp_(x.bp_), q_(x.q_), size_(x.size_) { } 

      void reserve(uint sz, uint max_dist=0)
      {
	bp_.resize(sz, max_dist);
	bp_.fill(0);
	q_.resize(sz);
	std::fill(q_.begin(), q_.end(), 1);
      }

      uint reserved_size() const { return q_.size(); }

      void resize(uint size, uint max_dist=0)
      {
	if (size>reserved_size() || max_dist!=this->max_dist()) {
	  reserve(size, max_dist);
	}
	size_ = size;
	bp_.fill(0);
	std::fill(q_.begin(), q_.end(), 1);
      }

      uint size() const { return size_; }

      uint max_dist() const { return bp_.max_dist(); };

      void update(uint i, uint j, T v)
      {
	value_type d = v-bp_(i,j);
	q_[i] -= d;
	q_[j] -= d;
	bp_(i,j) += d;
      }

      void add(uint i, uint j, T v)
      {
	q_[i] -= v;
	q_[j] -= v;
	bp_(i,j) += v;
      }

      T operator()(uint i, uint j) const { return bp_(i,j); }
      T operator[](uint i) const { return q_[i]; }

      template < class Seq, class RuleSet >
      void parse(const Seq& seq, const RuleSet& rules);

      bool parse(const std::string& str, bool ignore_alone_pair=false, uint minloop=3);

      static
      void
      convert_to_raw_sequences(const std::list<std::string>& ma,
			       std::list<std::string>& seqs,
			       std::list<std::vector<uint> >& idxmaps);
      
      template <class BPTablePtr>
      void
      average(const std::list<BPTablePtr>& bps,
	      const std::list<std::vector<uint> >& idxmaps,
              uint max_dist=0);

      bool load(const char* filename);
      bool save(const char* filename, const std::string& seq, float th) const;
      bool save(std::ostream& out, const std::string& seq, float th) const;

    private:
      CYKTable<T> bp_;
      std::vector<T> q_;
      uint size_;
    };
  }
};

#endif	// __INC_EM_H__

// Local Variables:
// mode: C++
// End:
