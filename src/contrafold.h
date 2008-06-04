// wrapper routines for CONTRAfold
#ifndef __INC_CONTRAFOLD_H__
#define __INC_CONTRAFOLD_H__

#include <string>

namespace CONTRAfold {
  template < class T >
  T* ComputePosterior(const std::string& seq);
};

#endif	// __INC_CONTRAFOLD_H__
