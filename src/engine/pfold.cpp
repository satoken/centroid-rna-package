/*
 * $Id:$
 *
 * CentroidFold: A generalized centroid estimator for predicting RNA
 *               secondary structures
 *
 * Copyright (C) 2008-2010 Kengo Sato
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include "engine/pfold.h"

static
std::string
mk_tempname(const char* base)
{
  char* b=strdup(base);
  int fd=mkstemp(b);
  std::string tempname(b);
  close(fd);
  return tempname;
}

static
void
mk_fa(const std::string& fname, const std::string& fa)
{
  std::ofstream os(fname.c_str());
  os << ">" << "noname" << std::endl
     << fa << std::endl;
}

static
void
mk_fa(const std::string& fname, const Aln& aln)
{
  std::ofstream os(fname.c_str());
  std::list<std::string>::const_iterator n;
  std::list<std::string>::const_iterator s;
  for (n=aln.name().begin(), s=aln.seq().begin(); n!=aln.name().end(); ++n, ++s)
  {
    os << ">" << *n << std::endl
       << *s << std::endl;
  }
}

template <class SEQ>
void
PfoldModel<SEQ>::
calculate_posterior(const SEQ& seq)
{
  const char* tempname="centroid_fold.XXXXXX";

  std::string fa, col, nj_col, ml_col, res_col, pp;
  // assume seq=ALN
  try {
    fa=mk_tempname(tempname);
    mk_fa(fa, seq);

    col=mk_tempname(tempname);
    {
      char* buf=(char*)malloc(awk_bin_.size()+pfold_bin_dir_.size()+fa.size()+sed_bin_.size()+col.size()+60);
      sprintf(buf, "%s -f %s/fasta2col %s | %s 's/arbitrary/RNA/g' > %s",
              awk_bin_.c_str(), pfold_bin_dir_.c_str(), fa.c_str(), sed_bin_.c_str(), col.c_str());
      if (system(buf)) throw std::runtime_error("cannot execute pfold");
      free(buf);
    }

    nj_col=mk_tempname(tempname);
    {
      char* buf=(char*)malloc(pfold_bin_dir_.size()*2+col.size()+nj_col.size()+40);
      sprintf(buf, "%s/findphyl %s/scfg.rate %s > %s",
              pfold_bin_dir_.c_str(), pfold_bin_dir_.c_str(), col.c_str(), nj_col.c_str());
      if (system(buf)) throw std::runtime_error("cannot execute pfold");
      free(buf);
    }

    ml_col=mk_tempname(tempname);
    {
      char* buf=(char*)malloc(pfold_bin_dir_.size()*2+nj_col.size()+ml_col.size()+40);
      sprintf(buf, "%s/mltree %s/scfg.rate %s > %s",
              pfold_bin_dir_.c_str(), pfold_bin_dir_.c_str(), nj_col.c_str(), ml_col.c_str());
      if (system(buf)) throw std::runtime_error("cannot execute pfold");
      free(buf);
    }

    res_col=mk_tempname(tempname);
    pp=mk_tempname(tempname);
    {
      char* buf=(char*)malloc(pfold_bin_dir_.size()*2+pp.size()+res_col.size()+60);
      sprintf(buf, "%s/scfg --treeinfile --ppfile %s %s/article.grm %s > %s",
              pfold_bin_dir_.c_str(), pp.c_str(), pfold_bin_dir_.c_str(), ml_col.c_str(), res_col.c_str());
      if (system(buf)) throw std::runtime_error("cannot execute pfold");
      free(buf);
    }
  } catch (...) {
#if 0
    std::cout << "fa: " << fa << std::endl
              << "col: " << col << std::endl
              << "nj_col: " << nj_col << std::endl
              << "ml_col: " << ml_col << std::endl
              << "res_col: " << res_col << std::endl
              << "pp: " << pp << std::endl;
#endif
#if 1
    unlink(fa.c_str());
    unlink(col.c_str());
    unlink(nj_col.c_str());
    unlink(ml_col.c_str());
    unlink(res_col.c_str());
    unlink(pp.c_str());
#endif
    throw;
  }
  
  std::ifstream is(pp.c_str());
  if (!is) throw "pfold failed";

  uint sz;
  float p;
  is >> sz;
  bp_.resize(sz);
  for (uint i=0; i!=sz; ++i)
  {
    for (uint j=0; j!=sz; ++j)
    {
      is >> p;
      if (i<j) bp_.update(i,j,p);
    }
  }

  unlink(fa.c_str());
  unlink(col.c_str());
  unlink(nj_col.c_str());
  unlink(ml_col.c_str());
  unlink(res_col.c_str());
  unlink(pp.c_str());
}

// instantiation
template class PfoldModel<std::string>;
template class PfoldModel<Aln>;
