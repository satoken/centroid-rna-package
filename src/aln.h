// $Id$

#ifndef __INC_ALN_H__
#define __INC_ALN_H__

#include <boost/spirit/iterator/file_iterator.hpp>
#include <string>
#include <list>

class Aln
{
public:
  Aln() : name_(), seq_()
  {}

  Aln(const std::list<std::string>& name,
      const std::list<std::string>& seq)
    : name_(name), seq_(seq)
  {}

  Aln(const Aln& aln)
    : name_(aln.name_), seq_(aln.seq_)
  {}

  Aln&
  operator=(const Aln& aln)
  {
    if (this != &aln) {
      name_ = aln.name_;
      seq_  = aln.seq_;
    }
    return *this;
  }
    
  const std::list<std::string>& name() const { return name_; }
  const std::list<std::string>& seq() const { return seq_; }
  unsigned int size() const { return seq_.begin()->size(); }
  unsigned int num_aln() const { return name_.size(); }
  std::string consensus() const;
  
  static
  unsigned int
  load(std::list<Aln>& data, const char* file );

  bool 
  load(boost::spirit::file_iterator<>& fi);

private:
  std::list<std::string> name_;
  std::list<std::string> seq_;
};

#endif	// __INC_ALN_H__

// Local Variables:
// mode: C++
// End:
