// $Id$

#ifndef __INC_CENTROID_H__
#define __INC_CENTROID_H__

#include <string>

namespace SCFG
{
  namespace Centroid
  {
    template < class T >
    double
    execute(const T& table, std::string& paren, double n=1.0, double th=1.0);
  };
};

#endif	//  __INC_CENTROID_H__

// Local Variables:
// mode: C++
// End:
