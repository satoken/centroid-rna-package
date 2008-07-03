// wrapper routines for CONTRAfold
#ifndef __INC_CONTRAFOLD_H__
#define __INC_CONTRAFOLD_H__

#include <string>
#include <vector>

namespace CONTRAfold {

  template < class T >
  struct WSImpl;

  template < class T >
  struct WS {
    WSImpl<T>* impl;

    WS(unsigned int size);
    ~WS();
  };
  
  template < class T >
  const T* ComputePosterior(const std::string& seq, bool canonical_only = true);

  template < class T >
  const T* ComputePosterior(const std::string& seq, WS<T>& ws, bool canonical_only = true);
};

#endif	// __INC_CONTRAFOLD_H__
