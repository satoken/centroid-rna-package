// 
#include <vector>
#include <string>
#include "InferenceEngine.hpp"
#include "ParameterManager.hpp"
#include "SStruct.hpp"
#include "Defaults.ipp"
#include "contrafold.h"

template < class T > 
struct CONTRAfold<T>::Impl
{
  Impl(bool canonical_only, int max_bp_dist);
  ~Impl() {};

  void SetParameters(const std::string& params);
  void SetConstraint(const std::string& paren);
  
  const T* ComputePosterior(const std::string& seq);
  const T* ComputePosterior(const std::string& seq, std::vector<T>& p);
  T ComputeInside(const std::string& seq);
  T ComputeLogPartitionCoefficient() const;
  T ComputeViterbi(const std::string& seq);

  void PrepareStochasticTraceback(const std::string& seq);
  std::vector<int> StochasticTraceback() const;

  ParameterManager<T> pm_;
  InferenceEngine<T> engine_;
  std::vector<T> w_;
  std::string paren_;
};

template < class T >
CONTRAfold<T>::
CONTRAfold(bool canonical_only, int max_bp_dist)
  : impl_(new Impl(canonical_only, max_bp_dist)), max_bp_dist_(max_bp_dist)
{
}

template < class T >
CONTRAfold<T>::
~CONTRAfold()
{
  delete impl_;
}

template < class T >
void
CONTRAfold<T>::
SetParameters(const std::string& params)
{
  impl_->SetParameters(params);
}

template < class T >
void
CONTRAfold<T>::
SetConstraint(const std::string& paren)
{
  impl_->SetConstraint(paren);
}

template < class T >
const T* 
CONTRAfold<T>::
ComputePosterior(const std::string& seq)
{
  return impl_->ComputePosterior(seq);
}

template < class T >
const T* 
CONTRAfold<T>::
ComputePosterior(const std::string& seq, std::vector<T>& p)
{
  return impl_->ComputePosterior(seq, p);
}

template < class T >
T
CONTRAfold<T>::
ComputeInside(const std::string& seq)
{
  return impl_->ComputeInside(seq);
}

template < class T >
T
CONTRAfold<T>::
ComputeLogPartitionCoefficient() const
{
  return impl_->ComputeLogPartitionCoefficient();
}

template < class T >
T
CONTRAfold<T>::
ComputeViterbi(const std::string& seq)
{
  return impl_->ComputeViterbi(seq);
}

template < class T >
void
CONTRAfold<T>::
PrepareStochasticTraceback(const std::string& seq)
{
  return impl_->PrepareStochasticTraceback(seq);
}

template < class T >
std::vector<int>
CONTRAfold<T>::
StochasticTraceback() const
{
  return impl_->StochasticTraceback();
}

template < class T > 
CONTRAfold<T>::Impl::
Impl(bool canonical_only, int max_bp_dist)
  : pm_(), engine_(!canonical_only, max_bp_dist), paren_()
{
  engine_.RegisterParameters(pm_);
  if (canonical_only)
    w_ = GetDefaultComplementaryValues<float>();
  else
    w_ = GetDefaultNoncomplementaryValues<float>();
}

template < class T >
void
CONTRAfold<T>::Impl::
SetParameters(const std::string& params)
{
  pm_.ReadFromFile(params, w_);
}

template < class T >
void
CONTRAfold<T>::Impl::
SetConstraint(const std::string& paren)
{
  paren_ = paren;
}

template < class T > 
const T*
CONTRAfold<T>::Impl::
ComputePosterior(const std::string& seq)
{
  SStruct* s=NULL;
  if (paren_.empty())
    s = new SStruct("unknown", seq);
  else
    s = new SStruct("unknown", seq, paren_);
  engine_.LoadSequence(*s);
  if (!paren_.empty())
    engine_.UseConstraints(s->GetMapping());
  engine_.LoadValues(w_);
  engine_.ComputeInside();
  engine_.ComputeOutside();
  engine_.ComputePosterior();
  delete s;
  return engine_.GetPosterior(0.0);
}

template < class T > 
const T*
CONTRAfold<T>::Impl::
ComputePosterior(const std::string& seq, std::vector<T>& p)
{
  SStruct* s=NULL;
  if (paren_.empty())
    s = new SStruct("unknown", seq);
  else
    s = new SStruct("unknown", seq, paren_);
  engine_.LoadSequence(*s);
  if (!paren_.empty())
    engine_.UseConstraints(s->GetMapping());
  engine_.LoadValues(w_);
  engine_.ComputeInside();
  engine_.ComputeOutside();
  engine_.ComputePosterior();
  delete s;
  return engine_.GetPosterior(0.0, p);
}

template < class T > 
T
CONTRAfold<T>::Impl::
ComputeInside(const std::string& seq)
{
  SStruct* s=NULL;
  if (paren_.empty())
    s = new SStruct("unknown", seq);
  else
    s = new SStruct("unknown", seq, paren_);
  engine_.LoadSequence(*s);
  if (!paren_.empty())
    engine_.UseConstraints(s->GetMapping());
  engine_.LoadValues(w_);
  engine_.ComputeInside();
  delete s;
  return engine_.ComputeLogPartitionCoefficient();
}

template < class T > 
T
CONTRAfold<T>::Impl::
ComputeLogPartitionCoefficient() const
{
  return engine_.ComputeLogPartitionCoefficient();
}

template < class T > 
T
CONTRAfold<T>::Impl::
ComputeViterbi(const std::string& seq)
{
  SStruct* s=NULL;
  if (paren_.empty())
    s = new SStruct("unknown", seq);
  else
    s = new SStruct("unknown", seq, paren_);
  engine_.LoadSequence(*s);
  if (!paren_.empty())
    engine_.UseConstraints(s->GetMapping());
  engine_.LoadValues(w_);
  engine_.ComputeViterbi();

  return engine_.GetViterbiScore();
}

template < class T > 
void
CONTRAfold<T>::Impl::
PrepareStochasticTraceback(const std::string& seq)
{
  SStruct* s=NULL;
  if (paren_.empty())
    s = new SStruct("unknown", seq);
  else
    s = new SStruct("unknown", seq, paren_);
  engine_.LoadSequence(*s);
  if (!paren_.empty())
    engine_.UseConstraints(s->GetMapping());
  engine_.LoadValues(w_);
  engine_.ComputeInside();
  delete s;
}

template < class T > 
std::vector<int>
CONTRAfold<T>::Impl::
StochasticTraceback() const
{
  return engine_.PredictPairingsStochasticTraceback();
}

template
class CONTRAfold<float>;
