// Copyright 2009, 2010 Michiaki Hamada

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <stack>
#include <cstdio>
#include <cstring>
#include <sys/time.h>
#include "engine/mccaskillhom.h"
#include "../util.h"

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
#endif

McCaskillHomModel::
McCaskillHomModel(const std::string& engine_a, bool canonical_only, uint max_bp_dist,
               const char* param /*=NULL*/, uint seed /*=0*/, bool run_as_mea /*=false*/)
  : FoldingEngine<TH>(run_as_mea, max_bp_dist),
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
  
  if (param && std::string(param) != "") Vienna::read_parameter_file(param);
  param_ = param;

  engine_a_ = engine_a;
}

void
McCaskillHomModel::
set_constraint(const std::string& str)
{
  str_=str;
  std::replace(str_.begin(), str_.end(), '.', 'x');
  std::replace(str_.begin(), str_.end(), '?', '.');
}

void
McCaskillHomModel::
calculate_posterior(const TH& th)
{
  const std::string& seq = th.first;
  const std::vector<std::string>& hom = th.second;

  typedef const CONTRALIGN::SparseMatrixEntry<float>* cIter;
  typedef CONTRALIGN::SparseMatrixEntry<float>* Iter;

  //PROBCONS::Probcons pc;
  //CONTRALIGN::CONTRAlign<float> ca;
  //CONTRAfold<float> cf2;
  CONTRALIGN::CONTRAlign<float>* ca = NULL;
  PROBCONS::Probcons* pc = NULL;
  if (engine_a_ == "CONTRAlign") ca =new CONTRALIGN::CONTRAlign<float> ();
  else if (engine_a_ == "ProbCons") pc = new PROBCONS::Probcons ();
  else {
    std::cerr << "No supported engine: " << engine_a_ << std::endl;
    assert (false);
  }

  std::vector< std::vector<float> > tmp;
  tmp.resize (seq.size()+1);
  for (uint k=0; k<tmp.size(); ++k) tmp[k].resize (seq.size()+1, 0.0);

  for (uint n=0; n<hom.size(); ++n) {
    CONTRALIGN::SparseMatrix<float>* ap = NULL;
    if (ca != NULL) contra_align (ap, *ca, seq, hom[n], 0.0001);
    else if (pc != NULL) probcons (ap, *pc, seq, hom[n], 0.0001);

    //ap->PrintSparse (std::cout);

    CONTRALIGN::SparseMatrix<float>* bp2 = NULL;
    pf_fold (bp2, hom[n], "", 0.01, param_);

    //bp2->PrintSparse (std::cout);

    const CONTRALIGN::SparseMatrix<float> apt (*ap, CONTRALIGN::SparseMatrix<float>::TRANSPOSE);
    for (uint i=1; i<=seq.size(); ++i) {
      for (cIter itr = ap->GetRowBegin (i); itr != ap->GetRowEnd (i); ++itr) {
	const int k = itr->column;
	const float pik = itr->value;
	for (cIter itr2 = bp2->GetRowBegin (k); itr2 != bp2->GetRowEnd(k); ++itr2) {
	  const int l = itr2->column;
	  const float pkl = itr2->value;
	  for (cIter itr3 = apt.GetRowBegin (l); itr3 != apt.GetRowEnd(l); ++itr3) {
	    const int j= itr3->column;
	    const float plj = itr3->value;
	    if ((int)i<j) tmp[i][j] += pik*pkl*plj; 
	  }
	}
      } // i
    }// j
    delete ap; delete bp2;
  } // n

  double alpha = 1.0 / (hom.size() + 1.0);

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
  for (uint j=2; j!=bp_.size()+1; ++j) {
    for (uint i=j-1; ; --i) {
#ifdef HAVE_VIENNA20
      FLT_OR_DBL* pr = Vienna::export_bppm();
      int* iindx = Vienna::get_iindx(seq.size());
#else
      FLT_OR_DBL* pr = Vienna::pr;
      int* iindx = Vienna::iindx;
#endif
      bp_.update(i-1, j-1, alpha * pr[iindx[i]-j] + (1-alpha) * tmp[i][j] / (float)hom.size());
      if (i==1) break;
    }
  }
  Vienna::free_pf_arrays();
  Vienna::fold_constrained = bk;
  if (str2) delete[] str2;
  delete pc; delete ca;
}  

/*
void
McCaskillModel::
prepare_stochastic_traceback(const std::string& seq)
{
  bk_st_back_=Vienna::st_back;
  Vienna::st_back=1;
  Vienna::pf_scale = -1;
  Vienna::init_pf_fold(seq.size());
  if (str_.empty()) {
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  } else {
    char *str2 = new char[str_.size()+1];
    strcpy(str2, str_.c_str());
    int bk = Vienna::fold_constrained;
    Vienna::fold_constrained = 1;
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), str2);
    delete[] str2;
    Vienna::fold_constrained = bk;
  }
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
*/
#endif

