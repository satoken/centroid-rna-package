// Copyright 2010 Michiaki Hamada

#ifndef __INC_ENGINE_CONTRAFOLDHOM_H__
#define __INC_ENGINE_CONTRAFOLDHOM_H__

#include "../folding_engine.h"

template < class RealT > class CONTRAfold;

// CentroidHomfold based on CONTRAfold model
class CONTRAfoldHomModel : public FoldingEngine<TH>
{
public:
  CONTRAfoldHomModel(const std::string& model, const std::string& engine_a, bool canonical_only, uint max_bp_dist,
		     uint seed=0, bool run_as_mea=false);
  virtual ~CONTRAfoldHomModel();

  // interface implementations
  virtual void set_constraint(const std::string& str);
  virtual void calculate_posterior(const TH& th);
  //virtual void prepare_stochastic_traceback(const std::string& seq);
  //virtual std::vector<int> stochastic_traceback(const std::string& seq);
  //virtual void clean_stochastic_traceback(const std::string& seq) { }

private:
  CONTRAfold<float>* contrafold_;
  std::string engine_a_;
};

#endif  //  __INC_ENGINE_CONTRAFOLDHOM_H__
