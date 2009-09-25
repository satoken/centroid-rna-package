// $Id:$

#ifndef __INC_RAND_H__
#define __INC_RAND_H__

#include <boost/random.hpp>

class Die
{
  typedef boost::mt19937 RNG;
  typedef boost::uniform_real<> DST;
  typedef boost::variate_generator<RNG&, DST> DIE;

public:
  Die(unsigned int seed)
    : rng_(seed), dst_(0, 1), die_(rng_, dst_)
  { }

  double operator()() { return die_(); }

private:
  RNG rng_;
  DST dst_;
  DIE die_;
};

#endif  //  __INC_RAND_H__
