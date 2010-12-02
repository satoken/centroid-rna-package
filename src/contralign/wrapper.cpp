// Copyright 2009, 2010 Michiaki Hamada

#include "Defaults.ipp"
#include "InferenceEngine.hpp"
#include "ParameterManager.hpp"
#include "contralign.h" // wrapper's header

namespace CONTRALIGN
{
  template<class T>
  struct CONTRAlign<T>::Impl
  {
    Impl ();
    ~Impl () {}

    const T* ComputePosterior (const std::string& seq1, const std::string& seq2, float th=0.0);
    const T* ComputePosterior (const std::string& seq1, const std::string& seq2, 
			       std::vector<T>& p, float th=0.0);

    InferenceEngine<T> inference_engine;
  };

  template < class T >
  CONTRAlign<T>::
  CONTRAlign () : impl_ (new Impl())
  {
  }

  template < class T >
  CONTRAlign<T>::
  ~CONTRAlign () 
  {
    delete impl_;
  }

  template < class T >
  const T*
  CONTRAlign<T>::
  ComputePosterior (const std::string& seq1, const std::string& seq2, float th)
  {
    return impl_->ComputePosterior (seq1, seq2, th);
  }

  template < class T >
  const T*
  CONTRAlign<T>::
  ComputePosterior (const std::string& seq1, const std::string& seq2, std::vector<T>& p, float th)
  {
    return impl_->ComputePosterior (seq1, seq2, p, th);
  }


  template < class T >
  CONTRAlign<T>::Impl::
  Impl () 
  {
    ParameterManager<T> parameter_manager;
    std::vector<T> w = GetDefaultRNAValues<T>();
    inference_engine.RegisterParameters (parameter_manager);
    inference_engine.LoadValues(w);    
  }

  template < class T >
  const T* 
  CONTRAlign<T>::Impl::
  ComputePosterior (const std::string& seq1, const std::string& seq2, float th)
  {
    MultiSequence seqs;
    seqs.AddSequence (new Sequence ('@'+seq1, "unkown", Sequence::UNKNOWN)); 
    seqs.AddSequence (new Sequence ('@'+seq2, "unkown", Sequence::UNKNOWN));     

    inference_engine.LoadSequences (seqs);

    inference_engine.ComputeForward();
    inference_engine.ComputeBackward();
    inference_engine.ComputePosterior();

    return inference_engine.GetPosterior (th);
  }

  template < class T >
  const T* 
  CONTRAlign<T>::Impl::
  ComputePosterior (const std::string& seq1, const std::string& seq2, std::vector<T>& p, float th)
  {
    MultiSequence seqs;
    seqs.AddSequence (new Sequence ('@'+seq1, "unkown", Sequence::UNKNOWN)); 
    seqs.AddSequence (new Sequence ('@'+seq2, "unkown", Sequence::UNKNOWN));     

    inference_engine.LoadSequences (seqs);

    inference_engine.ComputeForward();
    inference_engine.ComputeBackward();
    inference_engine.ComputePosterior();
    
    return inference_engine.GetPosterior (th, p);
  }

  template 
  class CONTRAlign<float>;
}
