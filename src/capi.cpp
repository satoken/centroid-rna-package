// $Id$
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <cstdlib>
#include <cstring>
#include "capi.h"
#include "centroid.h"
#include "bp.h"

typedef SCFG::BP::Table<double> BPTable;

int
centroidfold(const char *seq, char *str)
{
  BPTable bp(MAX_SEQ_LEN);
  bp.contra_fold(seq);
  unsigned int l=strlen(seq);
  std::string paren(l, '.');
  SCFG::Centroid::execute(bp, paren, 1.0, 1.0);
  strcpy(str, paren.c_str());

  return 0;
}
