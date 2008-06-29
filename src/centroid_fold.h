// $Id:$

#ifndef __INC_CENTROID_FOLD_H__
#define __INC_CENTROID_FOLD_H__

#include <string>
#include <list>
#include <vector>
#include <iosfwd>
#include "bp.h"

class CentroidFold
{
public:
  typedef SCFG::BP::Table<double> BPTable;

  enum {
    AUX,
    PFFOLD,
    CONTRAFOLD,
    ALIPFFOLD
  };
  
public:
  CentroidFold(unsigned int engine, bool run_as_mea=false,
	       unsigned int reserved_size=0);
  ~CentroidFold();

#if 0
  void execute(const std::string& seq, const std::vector<float>& gamma,
	       std::vector<std::string>& paren, std::vector<float>& ea,
	       const std::string& model="");
  void execute(const std::list<std::string>& seq, const std::vector<float>& gamma,
	       std::vector<std::string>& paren, std::vector<float>& ea,
	       const std::vector<std::string>& model);
  void execute(const std::vector<std::string>& seq, const std::vector<float>& gamma,
	       std::vector<std::string>& paren, std::vector<float>& ea,
	       const std::vector<std::string>& model);
  void execute(const std::list<std::string>& seq, const std::vector<float>& gamma,
	       std::vector<std::string>& paren, std::vector<float>& ea);
  void execute(const std::vector<std::string>& seq, const std::vector<float>& gamma,
	       std::vector<std::string>& paren, std::vector<float>& ea);
#endif

  void calculate_posterior(const std::string& seq, const std::string& model="");
  void calculate_posterior(const std::list<std::string>& seq,
			   const std::vector<std::string>& model);
  void calculate_posterior(const std::list<std::string>& seq);
  void calculate_posterior(const std::vector<std::string>& seq,
			   const std::vector<std::string>& model);
  void calculate_posterior(const std::vector<std::string>& seq);
  float decode_structure(float gamma, std::string& paren) const;
  std::pair<std::string,float> decode_structure(float gamma) const;
  void print(std::ostream& out, const std::string& name, const std::string& seq,
	     const std::vector<double>& gamma) const;
  void print_posterior(std::ostream& out, const std::string& seq, float th) const;
#ifdef HAVE_LIBRNA
  void ps_plot(const std::string& name, const std::string& seq, float g) const;
  void svg_plot(const std::string& name, const std::string& seq, float g) const;
#endif

private:
  unsigned int engine_;
  bool mea_;
  BPTable bp_;
};

#endif	// #ifndef __INC_CENTROID_FOLD_H__
