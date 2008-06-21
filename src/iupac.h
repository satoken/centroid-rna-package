// $Id$
#ifndef __INC_IUPAC_H__
#define __INC_IUPAC_H__

#include <vector>
#include <string>

class IUPACsymbol;

typedef std::vector<IUPACsymbol> IUPACsequence;

enum {
  RNA_A = 0,
  RNA_C = 1,
  RNA_G = 2,
  RNA_T = 3,
  RNA_U = RNA_T,
  N_RNA = 4,
  N_IUPAC = 16
};

class IUPACsymbol
{
  enum {
    GAP = 0,
    A = 1<<RNA_A,
    C = 1<<RNA_C,
    G = 1<<RNA_G,
    T = 1<<RNA_T,
    U = T,
    R = A | G,
    Y = C | T,
    M = A | C,
    K = G | T,
    S = C | G,
    W = A | T,
    B = C | G | T,
    D = A | G | T,
    H = A | C | T,
    V = A | C | G,
    N = A | C | G | T
  };

public:
  explicit IUPACsymbol(char c = '-');
  IUPACsymbol(const IUPACsymbol& s) : iupac_(s.iupac_) {}
  IUPACsymbol& operator=(const IUPACsymbol& s)
  {
    if (this!=&s) iupac_ = s.iupac_;
    return *this;
  }
  unsigned int index() const { return static_cast<unsigned int>(iupac_); }
  char to_char() const;

  static void transform(IUPACsequence& r, const std::string& s);

private:
  unsigned char iupac_;
};

#endif	// __INC_IUPAC_H__

// Local Variables:
// mode: C++
// End:
