/*
 * $Id$
 *
 * CentroidFold: A generalized centroid estimator for predicting RNA
 *               secondary structures
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

#include <cstdio>
#include <cstring>
#include <cmath>
#include <stack>
#include <sys/time.h>
#include "engine/alifold.h"

#ifdef HAVE_LIBRNA
namespace Vienna {
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/alifold.h>
#include <ViennaRNA/LPfold.h>
#include <ViennaRNA/PS_dot.h>
#include <ViennaRNA/aln_util.h>
#include <ViennaRNA/utils.h>
  extern int eos_debug;
  extern int   st_back;
  extern void read_parameter_file(const char fname[]);
};
};

extern "C" {
#include "engine/boltzmann_param.h"
};

AliFoldModel::
AliFoldModel(bool canonical_only, uint max_bp_dist,
             const char* param /*=NULL*/, uint seed /*=0*/, bool run_as_mea /*=false*/)
  : FoldingEngine<Aln>(run_as_mea, max_bp_dist), canonical_only_(canonical_only)
{
  if (!canonical_only_)
    Vienna::nonstandards = const_cast<char*>("AAACAGCACCCUGAGGUCUU");

  if (seed==0)
  {
    timeval t;
    gettimeofday(&t, NULL);
    seed = t.tv_usec * t.tv_sec;
  }
  srand((uint)seed);

  {
    using namespace Vienna;
    // copy from ViennaRNA-x.x.x/lib/utils.c
    xsubi[0] = xsubi[1] = xsubi[2] = (unsigned short) seed;  /* lower 16 bit */
    xsubi[1] += (unsigned short) ((unsigned)seed >> 6);
    xsubi[2] += (unsigned short) ((unsigned)seed >> 12);
  }

  copy_boltzmann_parameters();
  if (param) Vienna::read_parameter_file(param);
}

void
AliFoldModel::
set_constraint(const std::string& str)
{
  str_ = str;
  std::replace(str_.begin(), str_.end(), '.', 'x');
  std::replace(str_.begin(), str_.end(), '?', '.');
}

char**
alloc_aln(const Aln& aln)
{
  char** seqs = new char*[aln.num_aln()+1];
  seqs[aln.num_aln()]=NULL;
  std::list<std::string>::const_iterator x;
  uint i=0;
  for (x=aln.seq().begin(); x!=aln.seq().end(); ++x)
  {
    seqs[i] = new char[aln.size()+1];
    strcpy(seqs[i++], x->c_str());
  }
  return seqs;
}

void
free_aln(char** seqs)
{
  for (uint i=0; seqs[i]!=NULL; ++i) delete[] seqs[i];
  delete[] seqs;
}

void
AliFoldModel::
calculate_posterior(const Aln& aln)
{
  bp_.resize(aln.size());
  char** seqs=alloc_aln(aln);
  char* str2=new char[aln.size()+1];
  int bk=Vienna::fold_constrained;
  
  if (!str_.empty())
  {
    Vienna::fold_constrained=1;    
    strcpy(str2, str_.c_str());
  }
  // scaling parameters to avoid overflow
#ifdef HAVE_VIENNA20
  double min_en = Vienna::alifold((const char**)seqs, str2);
#else
  double min_en = Vienna::alifold(seqs, str2);
#endif
  double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(1.07*min_en)/kT/aln.size());
  Vienna::free_alifold_arrays();

#ifdef HAVE_VIENNA18
  Vienna::plist* pi;
#else
  Vienna::pair_info* pi;
#endif
  if (!str_.empty()) strcpy(str2, str_.c_str());
#ifdef HAVE_VIENNA20
  Vienna::alipf_fold((const char**)seqs, str2, &pi);
#else
  Vienna::alipf_fold(seqs, str2, &pi);
#endif
  for (uint k=0; pi[k].i!=0; ++k)
    bp_.update(pi[k].i-1, pi[k].j-1, pi[k].p);
  free(pi);

  Vienna::free_alipf_arrays();
  if (str2) delete[] str2;
  free_aln(seqs);
  Vienna::fold_constrained = bk;
}

void
AliFoldModel::
prepare_stochastic_traceback(const Aln& aln)
{
  char** seqs=alloc_aln(aln);
  char* str2=new char[aln.size()+1];
  int bk=Vienna::fold_constrained;
  
  if (!str_.empty())
  {
    Vienna::fold_constrained=1;    
    strcpy(str2, str_.c_str());
  }
  // scaling parameters to avoid overflow
#ifdef HAVE_VIENNA20
  double min_en = Vienna::alifold((const char**)seqs, str2);
#else
  double min_en = Vienna::alifold(seqs, str2);
#endif
  double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(1.07*min_en)/kT/aln.size());
  Vienna::free_alifold_arrays();

#ifdef HAVE_VIENNA18
  Vienna::plist* pi;
#else
  Vienna::pair_info* pi;
#endif
  if (!str_.empty()) strcpy(str2, str_.c_str());
#ifdef HAVE_VIENNA20
  Vienna::alipf_fold((const char**)seqs, str2, &pi);
#else
  Vienna::alipf_fold(seqs, str2, &pi);
#endif

  if (str2) delete[] str2;
  free_aln(seqs);
  Vienna::fold_constrained = bk;
}

std::vector<int>
AliFoldModel::
stochastic_traceback(const Aln& aln)
{
  double prob=1;
  char *str = Vienna::alipbacktrack(&prob);
  std::vector<int> paren(aln.size()+1, 0);
  std::stack<uint> st;
  for (uint i=0; i!=aln.size(); ++i) {
    if (str[i]=='(') {
      st.push(i);
    } else if (str[i]==')') {
      paren[st.top()+1]=i+1;
      paren[i+1]=st.top()+1;
      st.pop();
    }
  }
  free(str);
  return paren;
}

void
AliFoldModel::
clean_stochastic_traceback(const Aln& aln)
{
  Vienna::st_back=bk_st_back_;
  Vienna::free_alipf_arrays();
}

#endif


