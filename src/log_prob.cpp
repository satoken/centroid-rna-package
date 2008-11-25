// example: calculate the log-probability
// compile this program by:
//   g++ log_prob.cpp libcontrafold.a

#include <iostream>
#include <vector>
#include <string>
#include "contrafold.h"

int
main(int argc, char *argv[])
{
  bool canonical_only = true;

  while (argc>1)
  {
    if (argv[1][0]=='-')
    {
      switch (argv[1][1])
      {
      case 'n':
        // allow non-canonical base-pairs
        canonical_only = false;
        break;
      }
      argc--; argv++;
    } else {
      break;
    }
  }
  
  CONTRAfold<float> cf(canonical_only);

  std::string seq;
  std::string paren;
  
  while (std::cin >> seq >> paren)
  {
    float pf = cf.ComputeInside(seq);
    cf.SetConstraint(paren);
    float sc = cf.ComputeViterbi(seq);
    //float sc = cf.ComputeInside(seq);
    std::cout << sc - pf << std::endl;
  }
  
  return 0;
}
