// $Id$

#ifndef __INC_BP_H__
#define __INC_BP_H__

#include <iosfwd>
#include <list>
#include "cyktable.h"

namespace SCFG
{
  namespace BP
  {
    template < class T >
    class Table
    {
    public:
      typedef T value_type;
      
    public:
      Table() : bp_(), q_() { }
      
      Table(uint size) : bp_(size), q_(size)
      {
	bp_.fill(0);
	std::fill(q_.begin(), q_.end(), 1);
      }

      void resize(uint size)
      {
	bp_.resize(size);
	bp_.fill(0);
	q_.resize(size);
	std::fill(q_.begin(), q_.end(), 1);
      }

      uint size() const { return q_.size(); }

      void update(uint i, uint j, T v)
      {
	value_type d = v-bp_(i,j);
	q_[i] -= d;
	q_[j] -= d;
	bp_(i,j) += d;
      }

      void add(uint i, uint j, T v)
      {
	q_[i] -= v;
	q_[j] -= v;
	bp_(i,j) += v;
      }

      T operator()(uint i, uint j) const { return bp_(i,j); }
      T operator[](uint i) const { return q_[i]; }

      template < class Seq, class RuleSet >
      void parse(const Seq& seq, const RuleSet& rules);

      bool parse(const std::string& str, bool ignore_alone_pair=false, uint minloop=3);

#ifdef HAVE_LIBRNA
      void pf_fold(const std::string& seq);

      void alipf_fold(const std::list<std::string>& ma);
#endif

#ifdef HAVE_LIBCONTRAFOLD
      void contra_fold(const std::string& seq);
#endif

      static
      void
      convert_to_raw_sequences(const std::list<std::string>& ma,
			       std::list<std::string>& seqs,
			       std::list<std::vector<uint> >& idxmaps);
      
      template <class BPTablePtr>
      void
      average(const std::list<BPTablePtr>& bps,
	      const std::list<std::vector<uint> >& idxmaps);

      bool load(const char* filename);
      bool save(const char* filename, const std::string& seq, float th) const;
      bool save(std::ostream& out, const std::string& seq, float th) const;

    private:
      CYKTable<T> bp_;
      std::vector<T> q_;
    };
  }
};

#endif	// __INC_EM_H__

// Local Variables:
// mode: C++
// End:
