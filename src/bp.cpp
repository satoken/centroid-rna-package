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

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <cassert>
#include <cmath>
#include <stack>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include "bp.h"

template < class Table >
struct IgnoreAlonePair
{
  IgnoreAlonePair(Table& bp) : bp_(bp) { }

  void operator()(uint i, uint j)
  {
    if (i<j && bp_(i,j)>0) {
      uint c=0;
      if (i+1<j-1 && bp_(i+1,j-1)>0) c++;
      if (i>0 && j+1<bp_.size() && bp_(i-1,j+1)>0) c++;
      if (c==0)	bp_.update(i, j, 0);
    }
  }

  Table& bp_;
};

template < class V >
bool
BPTableTmpl<V>::
parse(const std::string& str, bool ignore_alone_pair/*=false*/, uint min_loop /*=3*/)
{
  resize(str.size());
  std::stack<uint> st;
  for (uint i=0; i!=str.size(); ++i) {
    switch (str[i]) {
      case '(':
        st.push(i);
        break;
      case ')':
        if (st.empty()) return false;
        if (min_loop+st.top()<i)
          update(st.top(), i, 1);
        st.pop();
        break;
      default:
        break;
    }
  }

  if (ignore_alone_pair) {
    IgnoreAlonePair<BPTableTmpl<V> > f(*this);
    inside_traverse(0, size()-1, f);
  }
      
  return st.empty();
}

template < class V >
bool
BPTableTmpl<V>::
load(const char* filename)
{
  std::string l;
  uint len=0;
  std::ifstream in(filename);
  if (!in) return false;
  while (std::getline(in, l)) ++len;
  resize(len);
  in.clear();
  in.seekg(0, std::ios::beg);
  while (std::getline(in, l)) {
    std::vector<std::string> v;
    boost::algorithm::split(v, l, boost::is_space(), boost::algorithm::token_compress_on);
    uint up = atoi(v[0].c_str());
    for (uint i=2; i!=v.size(); ++i) {
      uint down;
      float p;
      if (sscanf(v[i].c_str(), "%u:%f", &down, &p)==2) {
        update(up-1, down-1, p);
      }
    }
  }
  return true;
}

template < class V >
bool
BPTableTmpl<V>::
save(const char* filename, const std::string& seq, float th) const
{
  std::ofstream out(filename);
  if (!out) return false;
  return save(out, seq, th);
}

template < class V >
bool
BPTableTmpl<V>::
save(std::ostream& out, const std::string& seq, float th) const
{
  for (uint i=0; i!=seq.size(); ++i) {
    out << (i+1) << ' ' << seq[i] << ' ';
    uint l = max_dist()==0 ? seq.size() : std::min(size_t(i+max_dist()), seq.size());
    for (uint j=i; j!=l; ++j) {
      if ((*this)(i,j)>=th && (*this)(i,j)>0.0)
        out << (j+1) << ':' << (*this)(i,j) << ' ';
    }
    out << std::endl;
  }
  return true;
}

template <class V>
struct Updater
{
  Updater(const CYKTable<V>& bp, uint sz) : bp_(bp), q(sz, 1.0) { }

  void operator()(uint i, uint j)
  {
    if (i<j && bp_(i,j)>0.0)
    {
      q[i]-=bp_(i,j);
      q[j]-=bp_(i,j);
    }
  }
  CYKTable<V> bp_;
  std::vector<V> q;
};

template < class V >
std::vector<V>
BPTableTmpl<V>::
calc_nonbp_prob() const
{
  Updater<V> updater(bp_, size_);
  inside_traverse(0, size(), updater);
  return updater.q;
}


// instantiate
template
class BPTableTmpl<float>;

template
bool
BPTableTmpl<uint>::
parse(const std::string& str, bool ignore_alone_pair, uint min_loop);
