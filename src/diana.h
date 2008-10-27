// $Id$

#ifndef __INC_DIANA_H__
#define __INC_DIANA_H__

#include <vector>
#include <list>
#include <boost/range.hpp>

typedef unsigned int uint;

namespace HCLUST
{
  // Divisive Analysis: a hierarchical clustering by the top-down order
  template <class DMatrix>
  class Diana
  {
  public:
    typedef std::pair<uint*,uint*> range_t;
    typedef typename DMatrix::element value_type;

  private:
    struct clust_t
    {
      std::pair<uint,uint> ch;
      range_t range;
      value_type diameter;
    };
  
  public:
    Diana(const DMatrix& dmatrix);

    void build(uint max_clusters=static_cast<uint>(-1));

    double get_clusters(uint n_clusters, std::vector<uint>& objects,
			std::vector<uint>& num) const;

    uint optimal_size() const;

  private:
    value_type diameter(const clust_t& c) const;

    boost::range_const_iterator<range_t>::type find_seed_splinter(const clust_t& c) const;

    value_type splinter_value(uint x, const clust_t& c,
			      boost::range_const_iterator<range_t>::type dp) const;

    boost::range_const_iterator<range_t>::type divide(const clust_t& c);

    double cluster_index(uint k, double& out_c, double& in_c) const;

  private:
    const DMatrix& dmatrix_;
    std::vector<uint> objects_;
    std::vector<value_type> in_sum_;
    std::vector<clust_t> clusters_;
    std::list<uint> active_;
    std::vector<uint> hist_;
    std::vector<double> ch_index_;
    value_type in_c_;
    value_type out_c_;
  };
};

#endif //  __INC_DIANA_H__

// Local Variables:
// mode: C++
// End:
