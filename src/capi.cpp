// $Id$
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <cstdlib>
#include <cstring>
#include <new>
#include <stdexcept>
#include "capi.h"
#include "centroid_fold.h"

class CFWrapper
{
public:
  CFWrapper()
  {
    impl_ = new CentroidFold(CentroidFold::CONTRAFOLD);
  }

  ~CFWrapper()
  {
    delete impl_;
  }

  int
  centroid_fold(const char* seq, char* str, float g)
  {
    try {
      impl_->calculate_posterior(seq);
      unsigned int l=strlen(seq);
      std::string paren(l, '.');
      impl_->decode_structure(g, paren);
      strcpy(str, paren.c_str());
    } catch (std::bad_alloc& e) {
      return CENTROID_BAD_ALLOC;
    } catch (std::length_error& e) {
      return CENTROID_LENGTH_ERR;
    } catch (std::exception& e) {
      return CENTROID_UNKNOWN_ERR;
    }
    return CENTROID_SUCCESS;
  }

private:
  CentroidFold* impl_;
};

static CFWrapper cf;

int
centroidfold(const char *seq, char *str, float g)
{
  return cf.centroid_fold(seq, str, g);
}
