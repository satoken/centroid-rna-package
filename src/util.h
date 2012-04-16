#ifndef __H__UTIL__MH
#define __H__UTIL__MH

#include "probconsRNA/probcons.h"
#include "contralign/contralign.h"
#include "contralign/SparseMatrix.hpp"
#include "contrafold/contrafold.h"


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
#endif //__INC_LIBRNA_H
#endif

template < class T >
static 
void
contra_align (std::vector<T>& ap, CONTRALIGN::CONTRAlign<T>& ca, 
	      const std::string& seq1, const std::string& seq2, double th=0.0)
{
  ca.ComputePosterior (seq1, seq2, ap, th);
}

template < class T >
static
void
contra_align (CONTRALIGN::SparseMatrix<T>*& ap, CONTRALIGN::CONTRAlign<T>& ca,
	      const std::string& seq1, const std::string& seq2, double th)
{
  if (ap!=NULL) delete ap;
  std::vector<T> table;
  contra_align (table, ca, seq1, seq2, th);
  ap = new CONTRALIGN::SparseMatrix<T> (&table[0], seq1.size()+1, seq2.size()+1, 0.0);
}

template< class T >
static
void
probcons (std::vector<T>& ap, PROBCONS::Probcons& pc,
	  const std::string& seq1, const std::string& seq2, double th=0.0)
{
  pc.ComputePosterior (seq1, seq2, ap, th);
}

template < class T >
static
void
probcons (CONTRALIGN::SparseMatrix<T>*& ap, PROBCONS::Probcons& pc,
	  const std::string& seq1, const std::string& seq2, double th)
{
  if (ap!=NULL) delete ap;
  std::vector<T> table;
  probcons (table, pc, seq1, seq2, th);
  ap = new CONTRALIGN::SparseMatrix<T> (&table[0], seq1.size()+1, seq2.size()+1, 0.0);
}


template < class T, class U >
static
void
contra_fold(CONTRALIGN::SparseMatrix<T>*& bp, CONTRAfold<U>& cf, 
	    const std::string& seq, const std::string& str="", double th=0.0)
{
  if (bp != NULL) delete bp;
  std::vector<float> posterior;
  cf.SetConstraint(str);
  cf.ComputePosterior(seq, posterior);

  std::map<std::pair<int, int>, T> elems;
  uint k=0;
  for (uint i=0; i!=seq.size()+1; ++i) {
    for (uint j=i; j!=seq.size()+1; ++j) {
      if (i!=0 && posterior[k] > th) elems[std::make_pair(i,j)]=posterior[k]; 
      ++k;
    }
  }
  bp = new CONTRALIGN::SparseMatrix<T>(elems, seq.size()+1, seq.size()+1, 0);
  assert (bp != NULL);
}

template < class T >
static
void
pf_fold(CONTRALIGN::SparseMatrix<T>*& bp, const std::string& seq, const std::string& str="", double th=0.0, const char* param=NULL)
{
  //if (param) Vienna::read_parameter_file (param);
  if (bp!=NULL) delete bp;
  Vienna::pf_scale = -1;
#ifndef HAVE_VIENNA20
  Vienna::init_pf_fold(seq.size());
#endif
  if (str.empty()) {
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  } else {
    char *str2 = new char[str.size()+1];
    strcpy(str2, str.c_str());
    int bk = Vienna::fold_constrained;
    Vienna::fold_constrained = 1;
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), str2);
    delete[] str2;
    Vienna::fold_constrained = bk;
  }
  std::map<std::pair<int, int>, T> elems;
  for (uint j=2; j!=seq.size()+1; ++j) {
    for (uint i=j-1; ; --i) {
      double prob = Vienna::pr[Vienna::iindx[i]-j] >= th ? Vienna::pr[Vienna::iindx[i]-j] : 0.0; 
      if (prob != 0) elems[std::make_pair(i,j)] = prob; 
      if (i==1) break;
    }
  }
  Vienna::free_pf_arrays();
  bp = new CONTRALIGN::SparseMatrix<T>(elems, seq.size()+1, seq.size()+1, 0);
}

#endif
