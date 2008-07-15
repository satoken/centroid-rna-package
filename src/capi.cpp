// $Id$
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <cstdlib>
#include <cstring>
#include <new>
#include <stdexcept>
#include "capi.h"
#include "centroid.h"
#include "bp.h"

typedef SCFG::BP::Table<double> BPTable;

static BPTable bp;

int
centroidfold(const char *seq, char *str, float g)
{
  try {
    if (bp.reserved_size()<RESERVED_LENGTH)
      bp.reserve(RESERVED_LENGTH);
    bp.contra_fold(seq);
    unsigned int l=strlen(seq);
    std::string paren(l, '.');
    SCFG::Centroid::execute(bp, paren, g);
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
