// $Id$

#ifndef __INC_FA_H__
#define __INC_FA_H__

#include <boost/spirit/iterator/file_iterator.hpp>
#include <string>
#include <list>
//#include "rna.h"

class Fasta
{
public:
  Fasta() : name_(), seq_(), str_() { }

  Fasta(const std::string& name,
	const std::string& seq,
	const std::string& str="")
    : name_(name), seq_(seq), str_(str)
  { }

  Fasta(const Fasta& fa)
    : name_(fa.name_), seq_(fa.seq_), str_(fa.str_)
  { }
    
  Fasta&
  operator=(const Fasta& fa)
  {
    if (this != &fa) {
      name_ = fa.name_;
      seq_ = fa.seq_;
      str_ = fa.str_;
    }
    return *this;
  }

  const std::string& name() const { return name_; }
  const std::string& seq() const { return seq_; }
  const std::string& str() const { return str_; }
  uint size() const { return seq_.size(); }

  static
  uint
  load(std::list<Fasta>& data, const char* file);
  
  bool
  load(boost::spirit::file_iterator<>& fi);

private:
  std::string name_;
  std::string seq_;
  std::string str_;
};

#endif //  __INC_FA_H__

// Local Variables:
// mode: C++
// End:
