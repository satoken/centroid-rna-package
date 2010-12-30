// Copyright 2009, 2010 Michiaki Hamada

#ifndef __INC_ENGINE_MCCASKILLHOM_H__
#define __INC_ENGINE_MCCASKILLHOM_H__

#include "../folding_engine.h"

#ifdef HAVE_LIBRNA
class McCaskillHomModel : public FoldingEngine<TH>
{
public:
  McCaskillHomModel(const std::string& engine_a, bool canonical_only, uint max_bp_dist,
		    const char* param=NULL, uint seed=0, bool run_as_mea=false);
  virtual ~McCaskillHomModel() { }

  // interface implementations
  virtual void set_constraint(const std::string& str);
  virtual void calculate_posterior(const TH& seq);
/*   virtual void prepare_stochastic_traceback(const std::string& seq); */
/*   virtual std::vector<int> stochastic_traceback(const std::string& seq); */
/*   virtual void clean_stochastic_traceback(const std::string& seq); */

private:
  bool canonical_only_;
  int bk_st_back_;
  std::string str_;
  const char* param_;
  std::string engine_a_;
};
#endif

#endif  //  __INC_ENGINE_MCCASKILLHOM_H__
