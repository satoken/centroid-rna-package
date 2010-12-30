// Copyright 2009, 2010 Michiaki Hamada

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <sys/time.h>
#include "engine/contrafoldhom.h"
#include "../contrafold/contrafold.h"

#include "../util.h"

CONTRAfoldHomModel::
CONTRAfoldHomModel(const std::string& model, const std::string& engine_a, bool canonical_only, uint max_bp_dist,
		   uint seed /*=0*/, bool run_as_mea /*=false*/)
  : FoldingEngine<TH>(run_as_mea, max_bp_dist),
    contrafold_(NULL)
{
  contrafold_ = new CONTRAfold<float>(canonical_only, max_bp_dist);
  if (seed==0)
  {
    timeval t;
    gettimeofday(&t, NULL);
    seed = t.tv_usec * t.tv_sec;
  }
  contrafold_->init_rand(seed);
  if (!model.empty()) contrafold_->SetParameters(model);
  engine_a_ = engine_a;
}

CONTRAfoldHomModel::
~CONTRAfoldHomModel()
{
  if (contrafold_) delete contrafold_;
}

void
CONTRAfoldHomModel::
set_constraint(const std::string& str)
{
  contrafold_->SetConstraint(str);
}

inline
void
CONTRAfoldHomModel::
calculate_posterior(const TH& th)
{
  const std::string& seq = th.first;
  const std::vector<std::string>& hom = th.second;

  typedef const CONTRALIGN::SparseMatrixEntry<float>* cIter;
  typedef CONTRALIGN::SparseMatrixEntry<float>* Iter;

  //PROBCONS::Probcons pc;
  CONTRALIGN::CONTRAlign<float>* ca = NULL;
  PROBCONS::Probcons* pc = NULL;
  if (engine_a_ == "CONTRAlign") ca =new CONTRALIGN::CONTRAlign<float> ();
  else if (engine_a_ == "ProbCons") pc = new PROBCONS::Probcons ();
  else {
    std::cerr << "No supported engine: " << engine_a_ << std::endl;
    assert (false);
  }

  CONTRAfold<float> cf2;
  std::vector< std::vector<float> > tmp;
  tmp.resize (seq.size()+1);
  for (uint k=0; k<tmp.size(); ++k) tmp[k].resize (seq.size()+1, 0.0);

  for (uint n=0; n<hom.size(); ++n) {
    CONTRALIGN::SparseMatrix<float>* ap = NULL;
    //probcons (ap, pc, seq, hom[n], 0.0001);
    if (ca != NULL) contra_align (ap, *ca, seq, hom[n], 0.0001);
    else if (pc != NULL) probcons (ap, *pc, seq, hom[n], 0.0001);
    //ap->PrintSparse (std::cout);

    CONTRALIGN::SparseMatrix<float>* bp2 = NULL;
    contra_fold (bp2, cf2, hom[n], "", 0.01);

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

  bp_.resize(seq.size(), contrafold_->max_bp_dist());
  std::vector<float> posterior;
  contrafold_->ComputePosterior(seq, posterior);

  double alpha = 1.0 / (hom.size() + 1.0);
  //std::cout << "alpha=" << alpha << std::endl;
  if (contrafold_->max_bp_dist()==0) {
    uint k=0;
    for (uint i=0; i!=seq.size()+1; ++i) {
      for (uint j=i; j!=seq.size()+1; ++j) {
        if (i!=0) bp_.update(i-1, j-1, alpha * posterior[k] + (1-alpha) * tmp[i][j] / (float)hom.size());
        ++k;
      }
    }
  } else {
    uint k=0;
    for (uint i=0; i!=seq.size()+1; ++i) {
      for (uint j=i; j!=i+contrafold_->max_bp_dist(); ++j) {
        if (i!=0 && j<=seq.size()) bp_.update(i-1, j-1, alpha * posterior[k] + (1-alpha) * tmp[i][j] / (float)hom.size());
        ++k;
      }
    }
  }
  delete ca;  delete pc;
  //bp_.save (std::cout, seq, 0.01);
}

// void
// CONTRAfoldModel::
// prepare_stochastic_traceback(const std::string& seq)
// {
//   contrafold_->PrepareStochasticTraceback(seq);
// }

// std::vector<int>
// CONTRAfoldModel::
// stochastic_traceback(const std::string& seq)
// {
//   return contrafold_->StochasticTraceback();
// }

