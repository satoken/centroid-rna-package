// wrapper routines for CONTRAfold
#ifndef __INC_CONTRAFOLD_H__
#define __INC_CONTRAFOLD_H__

#include <string>

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
  T* ComputePosterior(const std::string& seq);

  template < class T >
  T* ComputePosterior(const std::string& seq, WS<T>& ws);
};

#endif	// __INC_CONTRAFOLD_H__
