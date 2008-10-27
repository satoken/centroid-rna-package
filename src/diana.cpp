// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <cmath>

#include <iostream>
#include <vector>
#include <list>
#include <utility>
#include <boost/range.hpp>
#include "diana.h"

using namespace HCLUST;
  
template <class DMatrix>
Diana<DMatrix>::
Diana(const DMatrix& dmatrix)
  : dmatrix_(dmatrix),
    objects_(dmatrix_.size()),
    in_sum_(dmatrix_.size(), 0.0),
    clusters_(dmatrix_.size()*2-1),
    active_(),
    hist_(dmatrix_.size()-1),
    ch_index_(dmatrix_.size()-1, 0.0),
    in_c_(0.0),
    out_c_(0.0)
{
}

template <class DMatrix>
void
Diana<DMatrix>::
build(uint max_clusters /*=static_cast<uint>(-1)*/)
{
  if (max_clusters<2) max_clusters=static_cast<uint>(-1);
  if (max_clusters!=static_cast<uint>(-1)) {
    hist_.resize(max_clusters-2);
    ch_index_.resize(max_clusters-2);
  }
  // initialize the object table
  for (uint i=0; i!=objects_.size(); ++i) {
    objects_[i] = i;
    in_sum_[i] = 0.0;
    for (uint j=0; j!=objects_.size(); ++j) {
      in_sum_[i] += dmatrix_[i][j];
      in_c_ += dmatrix_[i][j]*dmatrix_[i][j];
    }
  }

  // make the initial cluster
  uint last=0;
  clusters_[last].range = std::make_pair(&objects_[0], &objects_[0]+objects_.size());
  clusters_[last].diameter = diameter(clusters_[last]);
  active_.push_back(last);
  last++;

  uint h=0;
  while (!active_.empty() && h!=max_clusters-2) {
    // find maximum diamter cluster
    std::list<uint>::iterator max_x=active_.begin();
    value_type max_d=clusters_[*max_x].diameter;
    std::list<uint>::iterator x;
    for (x=active_.begin(); x!=active_.end(); ++x) {
      if (max_d<clusters_[*x].diameter) {
	max_d = clusters_[*x].diameter;
	max_x = x;
      }
    }
    uint cur = *max_x;
    active_.erase(max_x);
      
    if (boost::size(clusters_[cur].range)>1) {
      // divide the cluster 
      uint* new_dp = divide(clusters_[cur]);

      // make new clusters
      clusters_[cur].ch = std::make_pair(last, last+1);
      clusters_[last].range = std::make_pair(boost::begin(clusters_[cur].range), new_dp);
      clusters_[last].diameter = diameter(clusters_[last]);
      clusters_[last+1].range = std::make_pair(new_dp, boost::end(clusters_[cur].range));
      clusters_[last+1].diameter = diameter(clusters_[last+1]);

      // push them into queue
      active_.push_back(last);
      active_.push_back(last+1);
      last += 2;

      // calculate cluster index
#if 0
      cluster_index(h+2, out_c_, in_c_);
#else
      for (uint* x=boost::begin(clusters_[cur].range); x!=new_dp; ++x) {
	for (uint* y=new_dp; y!=boost::end(clusters_[cur].range); ++y) {
	  value_type v = dmatrix_[*x][*y];
	  v = v*v*2;
	  in_c_ -= v;
	  out_c_ += v;
	}
      }
#endif
      if (in_c_!=0.0 && objects_.size()-(h+2)!=0)
	ch_index_[h] = (out_c_/((h+2)-1))/(in_c_/(objects_.size()-(h+2)));
      else
	ch_index_[h] = 0.0;

      // record the divided cluster
      hist_[h++] = cur;
    }
  }
}

template <class DMatrix>
double
Diana<DMatrix>::
get_clusters(uint n_clusters, std::vector<uint>& objects, std::vector<uint>& num) const
{
  std::list<uint> a;
  a.push_back(hist_[0]);
  for (uint i=0; i!=n_clusters-1; ++i) {
    a.remove(hist_[i]);
    a.push_back(clusters_[hist_[i]].ch.first);
    a.push_back(clusters_[hist_[i]].ch.second);
  }

  objects.resize(objects_.size());
  num.resize(n_clusters);
  std::list<uint>::const_iterator x;
  uint c;
  uint n=0;
  for (x=a.begin(), c=0; x!=a.end(); ++x, ++c) {
    num[c] = boost::size(clusters_[*x].range);
    std::copy(boost::begin(clusters_[*x].range), boost::end(clusters_[*x].range),
	      &objects[n]);
    n += num[c];
  }

  return n_clusters>1 ? ch_index_[n_clusters-2] : 0.0;
}

template <class DMatrix>
uint
Diana<DMatrix>::
optimal_size() const
{
  std::vector<double>::const_iterator x=std::max_element(ch_index_.begin(), ch_index_.end());
  return std::distance(ch_index_.begin(), x)+2;
}

template <class DMatrix>
typename Diana<DMatrix>::value_type
Diana<DMatrix>::
diameter(const clust_t& c) const
{
  value_type d=0.0;
  boost::range_const_iterator<range_t>::type x, y;
  for (x=boost::begin(c.range); x!=boost::end(c.range); ++x) {
    for (y=x; y!=boost::end(c.range); ++y) {
      d=std::max(d, dmatrix_[*x][*y]);
    }
  }
  return d;
}

template <class DMatrix>
uint*
Diana<DMatrix>::
find_seed_splinter(const clust_t& c) const
{
  double max_d = 0.0;
  boost::range_const_iterator<range_t>::type max_obj = boost::begin(c.range);
  boost::range_const_iterator<range_t>::type x;
  for (x=boost::begin(c.range); x!=boost::end(c.range); ++x) {
#if 0
    boost::range_const_iterator<range_t>::type y;
    value_type d = 0.0;
    for (y=boost::begin(c.range); y!=boost::end(c.range); ++y)
      if (*x!=*y) d += dmatrix_[*x][*y];
    assert(in_sum_[*x]==d);
    d /= boost::size(c.range)-1;
#else
    value_type d = in_sum_[*x] / (boost::size(c.range)-1);
#endif
    if (d>max_d) {
      max_obj = x;
      max_d = d;
    }
  }
  return max_obj;
}

template <class DMatrix>
typename Diana<DMatrix>::value_type 
Diana<DMatrix>::
splinter_value(uint x, const clust_t& c,
	       boost::range_const_iterator<range_t>::type dp) const
{
  value_type new_dist=0.0, old_dist=0.0;
  uint new_cnt;
  boost::range_const_iterator<range_t>::type y;

  for (y=boost::begin(c.range), new_cnt=0; y!=dp; ++y, ++new_cnt)
    new_dist += dmatrix_[x][*y];
  new_dist /= new_cnt;

  uint old_cnt;
  for (y=dp, old_cnt=0; y!=boost::end(c.range); ++y, ++old_cnt)
    old_dist += dmatrix_[x][*y];
  assert(in_sum_[x]==old_dist);
  old_dist /= old_cnt;

  return old_dist - new_dist;
}

template <class DMatrix>
uint*
Diana<DMatrix>::
divide(const clust_t& c)
{
  // initialize
  std::vector<double> cur_in_sum(boost::size(c.range));
  boost::range_iterator<range_t>::type p;
  for (p=boost::begin(c.range); p!=boost::end(c.range); ++p)
    cur_in_sum[std::distance(boost::begin(c.range),p)] = in_sum_[*p];

  // find a seed
  boost::range_iterator<range_t>::type x = find_seed_splinter(c);
  boost::range_iterator<range_t>::type dp = boost::begin(c.range);
  uint m = *x;
  std::swap(*x, *dp);
  std::swap(cur_in_sum[std::distance(boost::begin(c.range), x)],
	    cur_in_sum[std::distance(boost::begin(c.range), dp)]);
  ++dp;
  in_sum_[m] = dmatrix_[m][m];
  for (x=dp; x!=boost::end(c.range); ++x) {
    in_sum_[*x] -= dmatrix_[m][*x];
  }

  // find the most far object from the old cluster
  while (1) {
    value_type max_v=0.0;
    boost::range_iterator<range_t>::type max_x = boost::end(c.range);
    for (x=dp; x!=boost::end(c.range); ++x) {
#if 0
      value_type v = splinter_value(*x, c, dp);
#else
      value_type old_dist = in_sum_[*x] / std::distance(dp, boost::end(c.range));
      value_type new_dist = cur_in_sum[std::distance(boost::begin(c.range), x)] - in_sum_[*x];
      new_dist /= std::distance(boost::begin(c.range), dp);
      value_type v = old_dist-new_dist;
#endif
      if (v>max_v) {
	max_v = v;
	max_x = x;
      }
    }
    if (max_x != boost::end(c.range)) {
      m = *max_x;
      std::swap(*max_x, *dp);
      std::swap(cur_in_sum[std::distance(boost::begin(c.range), max_x)],
		cur_in_sum[std::distance(boost::begin(c.range), dp)]);
      in_sum_[m] = dmatrix_[m][m];
      for (x=boost::begin(c.range); x!=dp; ++x) {
	in_sum_[m]  += dmatrix_[m][*x];
	in_sum_[*x] += dmatrix_[m][*x];
      }
      ++dp;
      for (x=dp; x!=boost::end(c.range); ++x) {
	in_sum_[*x] -= dmatrix_[m][*x];
      }
    } else {
      // no positive object requires any updates
      break;
    }
  }
  return dp;
}

template <class DMatrix>
double
Diana<DMatrix>::
cluster_index(uint k, double& out_c, double& in_c) const
{
  out_c = 0.0;
  in_c = 0.0;
  std::list<uint>::const_iterator x;
  for (x=active_.begin(); x!=active_.end(); ++x) {
    boost::range_const_iterator<range_t>::type y;
    for (y=boost::begin(clusters_[*x].range);
	 y!=boost::end(clusters_[*x].range); ++y) {
      boost::range_const_iterator<range_t>::type z;
      for (z=boost::begin(clusters_[0].range);
	   z!=boost::begin(clusters_[*x].range); z++) {
	out_c += dmatrix_[*y][*z]*dmatrix_[*y][*z];
      }
      for (z=boost::begin(clusters_[*x].range);
	   z!=boost::end(clusters_[*x].range); z++) {
	in_c += dmatrix_[*y][*z]*dmatrix_[*y][*z];
      }
      for (z=boost::end(clusters_[*x].range);
	   z!=boost::end(clusters_[0].range); z++) {
	out_c += dmatrix_[*y][*z]*dmatrix_[*y][*z];
      }
    }
  }
  return (out_c/(k-1))/(in_c/(objects_.size()-k));
}

// instantiation
#include <boost/multi_array.hpp>

typedef boost::multi_array<double,2> DMatrix;

template class Diana<DMatrix>;

#if 0
#include <iostream>
#include <cmath>

int
main(int argc, char* argv[])
{
  double data[] = { 1, 2, 5, 4, 10, 20, -1, 3, 2, 6, 7, 15, -10, 9, 12, -10, 2, 3, 4, 6 };
  uint sz = boost::size(data);

  typedef boost::multi_array<double,2> DMatrix;
  DMatrix dmatrix(boost::extents[sz][sz]);
  for (uint i=0; i!=sz; ++i) {
    for (uint j=0; j!=sz; ++j) {
      dmatrix[i][j]=std::fabs(data[i]-data[j]);
    }
  }
  
  Diana<DMatrix> diana(dmatrix);
  diana.build();
  std::vector<uint> res;
  std::vector<uint> num;
  for (uint n=2; n!=sz+1; ++n) {
    double ch=diana.get_clusters(n, res, num);
    std::cout << ch << ": ";
    for (uint i=0, k=0; i!=num.size(); ++i) {
      for (uint j=0; j!=num[i]; ++j, ++k) {
	std::cout << data[res[k]] << " ";
      }
      std::cout << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << diana.optimal_size() << std::endl;
  {
    uint n=diana.optimal_size();
    double ch=diana.get_clusters(n, res, num);
    std::cout << ch << ": ";
    for (uint i=0, k=0; i!=num.size(); ++i) {
      for (uint j=0; j!=num[i]; ++j, ++k) {
	std::cout << data[res[k]] << " ";
      }
      std::cout << ", ";
    }
    std::cout << std::endl;
  }
  return 0;
}
#endif
