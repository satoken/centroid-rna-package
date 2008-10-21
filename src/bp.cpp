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

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <cassert>
#include <cmath>
#include <stack>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
//#include "inside.h"
//#include "outside.h"
#include "bp.h"
#include "rule.h"
//#include "dptable.h"

#ifdef HAVE_LIBRNA
namespace Vienna {
extern "C" {
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/alifold.h>
#include <ViennaRNA/PS_dot.h>
#if 0
  extern int pfl_fold(char *sequence, int winSize, int pairdist,
		      float cutoff, struct plist **pl);
  extern void init_pf_foldLP(int length);
  extern void free_pf_arraysLP(void);
#endif
  extern char* pbacktrack(char *sequence);
  extern int   st_back;
};
};
#endif

#ifdef HAVE_LIBCONTRAFOLD
#include <contrafold.h>
#endif

template < class BPTable >
struct IgnoreAlonePair
{
  IgnoreAlonePair(BPTable& bp) : bp_(bp) { }

  void operator()(uint i, uint j)
  {
    if (i<j && bp_(i,j)>0) {
      uint c=0;
      if (i+1<j-1 && bp_(i+1,j-1)>0) c++;
      if (i>0 && j+1<bp_.size() && bp_(i-1,j+1)>0) c++;
      if (c==0)	bp_.update(i, j, 0);
    }
  }

  BPTable& bp_;
};

template < class BPTable, class BPTablePtr >
struct AverageBP
{
  AverageBP(BPTable& bp,
	    const std::list<BPTablePtr>& bps,
	    const std::list<std::vector<uint> >& idxmap)
    : bp_(bp), bps_(bps), idxmap_(idxmap)
  {
  }

  void operator()(uint i, uint j)
  {
    double v=0.0;
    uint c=0;
    typename std::list<BPTablePtr>::const_iterator b = bps_.begin();
    std::list<std::vector<uint> >::const_iterator idx = idxmap_.begin();
    while (b!=bps_.end() && idx!=idxmap_.end()) {
      uint ii=(*idx)[i];
      uint jj=(*idx)[j];
      if (ii!=static_cast<uint>(-1) &&
	  jj!=static_cast<uint>(-1)) {
	v += (**b)(ii,jj);
      }
      ++c; ++b; ++idx;
    }
    bp_.update(i, j, v/c);
  }

  void make()
  {
    SCFG::inside_traverse(0, bp_.size()-1, *this);
  }

private:
  BPTable& bp_;
  const std::list<BPTablePtr>& bps_;
  const std::list<std::vector<uint> >& idxmap_;
};

namespace SCFG
{
  namespace BP
  {
    template < class V >
    bool
    Table<V>::
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
	IgnoreAlonePair<Table<V> > f(*this);
	SCFG::inside_traverse(0, size()-1, f);
      }
      
      return st.empty();
    }

#ifdef HAVE_LIBRNA
    template < class V >
    void
    Table<V>::
    pf_fold(const std::string& seq)
    {
      resize(seq.size());
      Vienna::pf_scale = -1;
      Vienna::init_pf_fold(seq.size());
      Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
      for (uint j=2; j!=size()+1; ++j) {
	for (uint i=j-1; ; --i) {
	  update(i-1, j-1, Vienna::pr[Vienna::iindx[i]-j]);
	  if (i==1) break;
	}
      }
      Vienna::free_pf_arrays();
    }

    template < class V >
    void
    Table<V>::
    alipf_fold(const std::list<std::string>& ma)
    {
      resize(ma.front().size());
      // prepare an alignment
      uint length = ma.front().size();
      char** seqs = new char*[ma.size()+1];
      seqs[ma.size()]=NULL;
      std::list<std::string>::const_iterator x;
      uint i=0;
      for (x=ma.begin(); x!=ma.end(); ++x) {
	assert(x->size()==length);
	seqs[i] = new char[length+1];
	strcpy(seqs[i++], x->c_str());
      }
      {
	// scaling parameters to avoid overflow
	char* str = new char[length+1];
	double min_en = Vienna::alifold(seqs, str);
	delete[] str;
	Vienna::free_alifold_arrays();
	double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
	Vienna::pf_scale = exp(-(1.07*min_en)/kT/length);
      }
      // build a base pair probablity matrix
      Vienna::pair_info* pi;
      Vienna::alipf_fold(seqs, NULL, &pi);
      for (uint k=0; pi[k].i!=0; ++k)
	update(pi[k].i-1, pi[k].j-1, pi[k].p);
      free(pi);
      for (uint i=0; seqs[i]!=NULL; ++i) delete[] seqs[i];
      delete[] seqs;
      
    }
#endif

#ifdef HAVE_LIBCONTRAFOLD
    template < class V >
    void
    Table<V>::
    contra_fold(const std::string& seq)
    {
      resize(seq.size());
      if (contrafold_==NULL)
	contrafold_ = boost::shared_ptr<CONTRAfold<float> >(new CONTRAfold<float>(true));
      std::vector<float> posterior;
      contrafold_->ComputePosterior(seq, posterior);
      uint k=0;
      for (uint i=0; i!=seq.size()+1; ++i) {
	for (uint j=i; j!=seq.size()+1; ++j) {
	  if (i!=0) update(i-1, j-1, posterior[k]);
	  ++k;
	}
      }
    }

    template < class V >
    void
    Table<V>::
    contra_fold(const std::string& seq, const std::string& model)
    {
      resize(seq.size());
      if (contrafold_==NULL)
	contrafold_ = boost::shared_ptr<CONTRAfold<float> >(new CONTRAfold<float>(true));
      contrafold_->SetParameters(model);
      std::vector<float> posterior;
      contrafold_->ComputePosterior(seq, posterior);
      uint k=0;
      for (uint i=0; i!=seq.size()+1; ++i) {
	for (uint j=i; j!=seq.size()+1; ++j) {
	  if (i!=0) update(i-1, j-1, posterior[k]);
	  ++k;
	}
      }
    }
#endif

    template < class V >
    void
    Table<V>::
    convert_to_raw_sequences(const std::list<std::string>& ma,
			     std::list<std::string>& seqs,
			     std::list<std::vector<uint> >& idxmaps)
    {
      std::list<std::string>::const_iterator x;
      for (x=ma.begin(); x!=ma.end(); ++x) {
	std::string s;
	std::vector<uint> idxmap(x->size(), static_cast<uint>(-1));
	for (uint i=0, j=0; i!=x->size(); ++i) {
	  if ((*x)[i]!='-') {
	    s += (*x)[i];
	    idxmap[i] = j++;
	  }
	}
	seqs.push_back(s);
	idxmaps.push_back(idxmap);
      }
    }

    template < class V >
    template <class BPTablePtr>
    void
    Table<V>::
    average(const std::list<BPTablePtr>& bps,
	    const std::list<std::vector<uint> >& idxmaps)
    {
      resize(idxmaps.front().size());
      AverageBP<Table<V>, BPTablePtr> avg(*this, bps, idxmaps);
      avg.make();
    }

    template < class V >
    bool
    Table<V>::
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
    Table<V>::
    save(const char* filename, const std::string& seq, float th) const
    {
      std::ofstream out(filename);
      if (!out) return false;
      return save(out, seq, th);
    }

    template < class V >
    bool
    Table<V>::
    save(std::ostream& out, const std::string& seq, float th) const
    {
      for (uint i=0; i!=seq.size(); ++i) {
	out << (i+1) << ' ' << seq[i] << ' ';
	for (uint j=i; j!=seq.size(); ++j) {
	  if ((*this)(i,j)>=th && (*this)(i,j)>0.0)
	    out << (j+1) << ':' << (*this)(i,j) << ' ';
	}
	out << std::endl;
      }
      return true;
    }
  };
};


// instantiate
#include <boost/shared_ptr.hpp>
//#include "rna.h"
#include "iupac.h"
//#include "log_value.h"

template
class SCFG::BP::Table<double>;

#if 0
template
void
SCFG::BP::Table<double>::
parse(const RNASequence& seq,
      const Rule::Set<LogValue<double> >& rules);
#endif

#if 0
template
void
SCFG::BP::Table<double>::
parse(const IUPACsequence& seq,
      const Rule::Set<LogValue<double> >& rules);
#endif

template
bool
SCFG::BP::Table<uint>::
parse(const std::string& str, bool ignore_alone_pair, uint min_loop);

template
void
SCFG::BP::Table<double>::
average(const std::list<boost::shared_ptr<SCFG::BP::Table<double> > >& bps,
	const std::list<std::vector<uint> >& idxmaps);
