#include <stdio.h>
#include "capi.h"

int
main()
{
  char seq[]="augcaugcaugcaugcaugcaugc";
  char str[]="........................";
  centroidfold(seq, str);
  printf("%s\n%s\n", seq, str);
  return 0;
}
