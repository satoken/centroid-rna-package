#include <stdio.h>
#include "capi.h"

int
main()
{
  char seq[]="augcaugcaunnnngcaugcaugcaugc";
  char str[]="............................";
  int x=centroidfold(seq, str);
  switch (x) {
  case CENTROID_SUCCESS:
    printf("%s\n%s\n", seq, str);
    break;
  default:
    printf("err: %d\n", x);
    break;
  }
  return 0;
}
