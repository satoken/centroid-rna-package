// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include "fa.h"
#include <iostream>
#include <sstream>
#include <string>
#include <boost/spirit/core.hpp>
#include <boost/spirit/utility/chset.hpp>
#include <boost/spirit/utility/confix.hpp>
#include <boost/spirit/iterator/file_iterator.hpp>

using namespace boost::spirit;

struct fa_parser : public grammar< fa_parser >
{
  fa_parser(std::string& n, std::string& s, std::string& st)
    : name(n), seq(s), str(st) { }

  std::string& name;
  std::string& seq;
  std::string& str;

  struct append_seq
  {
    std::string& seq;

    append_seq(std::string& s) : seq(s) { }

    template <class Ite>
    void operator()(Ite i1, Ite i2) const
    {
      std::string s(i1, i2);
      seq += s;
    }
  };


  template <class ScannerT>
  struct definition
  {
    typedef rule<ScannerT> rule_t;
    rule_t fa;
    rule_t head;
    rule_t seq_l;
    rule_t seq;
    rule_t str;
    rule_t str_l;

    definition(const fa_parser& self)
    {
      fa = head >> +(seq | str);
      head = ch_p('>') >> (*(blank_p | graph_p))[assign_a(self.name)] >> eol_p;
      seq_l = *(alpha_p | ch_p('-'));
      str_l = *(chset_p("().") | blank_p);
      seq = +(seq_l[append_seq(self.seq)] >> eol_p);
      str = +(str_l[append_seq(self.str)] >> eol_p);
    }

    const rule_t& start() const { return fa; }
  };
};

bool
Fasta::
load(file_iterator<>& fi)
{
  file_iterator<> s = fi;
  seq_.clear();
  fa_parser parser(name_, seq_, str_);
  parse_info<file_iterator<> > info =  parse(fi, fi.make_end(), parser);
#if 0
  std::cout << "name: " << name() << std::endl
	    << " seq: " << seq() << std::endl
	    << " str: " << str() << std::endl;
#endif
  if (!info.hit) {
    fi = s;
    return false;
  } else {
    fi = info.stop;
    return true;
  }
}

// static
unsigned int
Fasta::
load(std::list<Fasta>& data, const char* filename)
{
  unsigned int n=0;
  file_iterator<> fi(filename);
  if (!fi) {
    std::ostringstream os;
    os << filename << ": no such file";
    throw os.str().c_str();
    //return false;
  }
  while (1) {
    Fasta fa;
    if (fa.load(fi)) {
      n++;
      data.push_back(fa);
    } else {
      break;
    }
  }
  return n;
}
