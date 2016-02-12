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

#include <stack>
#include <cstdio>
#include <cstring>
#include <sys/time.h>
#include "engine/mccaskill.h"

#ifdef HAVE_LIBRNA
#ifndef __INC_LIBRNA_H
#define __INC_LIBRNA_H
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
#ifndef FLT_OR_DBL
typedef Vienna::FLT_OR_DBL FLT_OR_DBL;
#endif
#endif

extern "C" {
#include "engine/boltzmann_param.h"
};

McCaskillModel::
McCaskillModel(bool canonical_only, uint max_bp_dist,
               const char* param /*=NULL*/, uint seed /*=0*/, bool run_as_mea /*=false*/)
  : FoldingEngine<std::string>(run_as_mea, max_bp_dist),
    canonical_only_(canonical_only), bk_st_back_(Vienna::st_back), str_()
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

  if (param)
  {
    if (strcmp(param, "default")!=0)
      Vienna::read_parameter_file(param);
  }
  else
    copy_boltzmann_parameters();
}

void
McCaskillModel::
set_constraint(const std::string& str)
{
  str_=str;
  std::replace(str_.begin(), str_.end(), '.', 'x');
  std::replace(str_.begin(), str_.end(), '?', '.');
}

void
McCaskillModel::
calculate_posterior(const std::string& seq)
{
  bp_.resize(seq.size());
  Vienna::pf_scale = -1;
  char *str2 = new char[seq.size()+1];
  int bk = Vienna::fold_constrained;
  if (!str_.empty())
  {
    Vienna::fold_constrained = 1;
    strcpy(str2, str_.c_str());
  }
  if (seq.size()>1600)
  {
    double min_en = Vienna::fold(seq.c_str(), str2);
    double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
    Vienna::pf_scale = exp(-(1.07*min_en)/kT/seq.size());
  }
#ifndef HAVE_VIENNA20
  Vienna::init_pf_fold(seq.size());
#endif
  if (str_.empty()) {
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  } else {
    strcpy(str2, str_.c_str());
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), str2);
  }
#ifdef HAVE_VIENNA20
  FLT_OR_DBL* pr = Vienna::export_bppm();
  int* iindx = Vienna::get_iindx(seq.size());
#else
  FLT_OR_DBL* pr = Vienna::pr;
  int* iindx = Vienna::iindx;
#endif
  for (uint j=2; j!=bp_.size()+1; ++j) {
    for (uint i=j-1; ; --i) {
      bp_.update(i-1, j-1, pr[iindx[i]-j]);
      if (i==1) break;
    }
  }
  Vienna::free_pf_arrays();
  Vienna::fold_constrained = bk;
#ifdef HAVE_VIENNA20
  free(iindx);
#endif
  if (str2) delete[] str2;
}  

void
McCaskillModel::
prepare_stochastic_traceback(const std::string& seq)
{
  bk_st_back_=Vienna::st_back;
  Vienna::st_back=1;
  Vienna::pf_scale = -1;
  char *str2 = new char[seq.size()+1];
  int bk = Vienna::fold_constrained;
  if (!str_.empty())
  {
    Vienna::fold_constrained = 1;
    strcpy(str2, str_.c_str());
  }
  if (seq.size()>1600)
  {
    double min_en = Vienna::fold(seq.c_str(), str2);
    double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
    Vienna::pf_scale = exp(-(1.07*min_en)/kT/seq.size());
  }
#ifndef HAVE_VIENNA20
  Vienna::init_pf_fold(seq.size());
#endif
  if (str_.empty()) {
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  } else {
    strcpy(str2, str_.c_str());
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), str2);
  }
  Vienna::fold_constrained = bk;
  if (str2) delete[] str2;
}  

std::vector<int>
McCaskillModel::
stochastic_traceback(const std::string& seq)
{
  char *str = Vienna::pbacktrack(const_cast<char*>(seq.c_str()));
  std::vector<int> paren(seq.size()+1, 0);
  std::stack<uint> st;
  for (uint i=0; i!=seq.size(); ++i) {
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
McCaskillModel::
clean_stochastic_traceback(const std::string& seq)
{
  Vienna::st_back=bk_st_back_;
  Vienna::free_pf_arrays();
}
#endif

