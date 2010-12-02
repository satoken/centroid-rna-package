// Copyright 2009, 2010 Michiaki Hamada

#ifndef __H_INCLUDE_PROBCONS_MH
#define __H_INCLUDE_PROBCONS_MH

namespace PROBCONS {
  class Probcons
  {
  private:
    struct Impl;

  public:
    Probcons ();
    ~Probcons ();
    
    const float* ComputePosterior (const std::string& seq1, const std::string& seq2, float th=0.0);
    const float* ComputePosterior (const std::string& seq1, const std::string& seq2, 
				   std::vector<float>& p, float th=0.0);
    void PrintParameters (const char* message) const;
  private:
    Impl* impl_;
  };
}
#endif // __H_INCLUDE_PROBCONS_MH
