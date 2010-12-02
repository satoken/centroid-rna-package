// Copyright 2009, 2010 Michiaki Hamada

#include "Defaults.h"
#include "ScoreType.h"
#include "ProbabilisticModel.h"
//
#include "probcons.h"

#include <iomanip>

namespace PROBCONS {

  VF initDistrib (NumMatrixTypes);
  VF gapOpen (2*NumInsertStates);
  VF gapExtend (2*NumInsertStates);
  VVF emitPairs (256, VF (256, 1e-10));
  VF emitSingle (256, 1e-5);
  string alphabet = alphabetDefault;

  struct Probcons::Impl
  {
    Impl ();
    ~Impl () {}

    const float* ComputePosterior (const std::string& seq1, const std::string& seq2, float th=0.0);
    const float* ComputePosterior (const std::string& seq1, const std::string& seq2, 
				   std::vector<float>& p, float th=0.0);
    void ReadParameters ();
    void PrintParameters (const char* message) const;
  };

  Probcons::Probcons () : impl_(new Impl ())
  {
  }

  Probcons::~Probcons () 
  {
    delete impl_;
  }

  const float* Probcons::ComputePosterior (const std::string& seq1, const std::string& seq2, float th)
  {
    return impl_->ComputePosterior (seq1, seq2, th);
  }

  const float* Probcons::ComputePosterior (const std::string& seq1, const std::string& seq2, std::vector<float>& p, 
					   float th)
  {
    return impl_->ComputePosterior (seq1, seq2, p, th);
  }

  /*
  const float* Probcons::ComputePosterior (const std::string& seq1, const std::string& seq2, 
					   std::vector<float>& p, float th)
  {
    return impl_->ComputePosterior (seq1, seq2, p, th);
  }
  */
  void Probcons::PrintParameters (const char* message) const
  {
    impl_->PrintParameters (message);
  }
  
  ///
  Probcons::Impl::Impl()
  {
    ReadParameters ();
  }

  const float* Probcons::Impl::ComputePosterior (const std::string& seq1str, const std::string& seq2str,
						 float th)
  {
    SafeVector<char>* data1 = new SafeVector<char> ();
    data1->resize (seq1str.size()+1);
    (*data1)[0] = '@';
    copy (seq1str.begin(), seq1str.end(), data1->begin()+1);
    Sequence seq1 (data1, "", seq1str.size(), 0, 0);

    SafeVector<char>* data2 = new SafeVector<char> ();
    data2->resize (seq2str.size()+1);
    (*data2)[0]='@';
    copy (seq2str.begin(), seq2str.end(), data2->begin()+1);
    Sequence seq2 (data2, "", seq2str.size(), 0, 0);

    ProbabilisticModel model (initDistrib, gapOpen, gapExtend, emitPairs, emitSingle);
    VF* forward  = model.ComputeForwardMatrix  (&seq1, &seq2); assert (forward);
    VF* backward = model.ComputeBackwardMatrix (&seq1, &seq2); assert (backward);
    VF* posterior = model.ComputePosteriorMatrix (&seq1, &seq2, *forward, *backward);
    delete forward;
    delete backward;
    unsigned s = posterior->size();
    float* ret = new float[s];
    for (unsigned i=0; i<s; i++) {
      ret[i] = ( (*posterior)[i] >= th ? (*posterior)[i] : float(0) );
      assert (ret[i]<=1.0);
    }
    delete posterior;
    return ret;
  }

  const float* Probcons::Impl::ComputePosterior (const std::string& seq1str, const std::string& seq2str,
						 std::vector<float>& p, float th)
  {
    SafeVector<char>* data1 = new SafeVector<char> ();
    data1->resize (seq1str.size()+1);
    (*data1)[0] = '@';
    copy (seq1str.begin(), seq1str.end(), data1->begin()+1);
    Sequence seq1 (data1, "", seq1str.size(), 0, 0);

    SafeVector<char>* data2 = new SafeVector<char> ();
    data2->resize (seq2str.size()+1);
    (*data2)[0]='@';
    copy (seq2str.begin(), seq2str.end(), data2->begin()+1);
    Sequence seq2 (data2, "", seq2str.size(), 0, 0);

    ProbabilisticModel model (initDistrib, gapOpen, gapExtend, emitPairs, emitSingle);
    VF* forward  = model.ComputeForwardMatrix  (&seq1, &seq2); assert (forward);
    VF* backward = model.ComputeBackwardMatrix (&seq1, &seq2); assert (backward);
    VF* posterior = model.ComputePosteriorMatrix (&seq1, &seq2, *forward, *backward);
    delete forward;
    delete backward;
    unsigned s = posterior->size();
    p.resize (s);
    //float* ret = new float[s];
    for (unsigned i=0; i<s; i++) {
      p[i] = ( (*posterior)[i] >= th ? (*posterior)[i] : float(0) );
      assert (p[i]<=1);
    }
    delete posterior;
    return &p[0];
  }
  

  void Probcons::Impl::ReadParameters ()
  {
    emitPairs = VVF (256, VF (256, 1e-10));
    emitSingle = VF (256, 1e-5);

    if (NumInsertStates == 1){
      for (int i = 0; i < NumMatrixTypes; i++) initDistrib[i] = initDistrib1Default[i];
      for (int i = 0; i < 2*NumInsertStates; i++) gapOpen[i] = gapOpen1Default[i];
      for (int i = 0; i < 2*NumInsertStates; i++) gapExtend[i] = gapExtend1Default[i];
    }
    else if (NumInsertStates == 2){
      for (int i = 0; i < NumMatrixTypes; i++) initDistrib[i] = initDistrib2Default[i];
      for (int i = 0; i < 2*NumInsertStates; i++) gapOpen[i] = gapOpen2Default[i];
      for (int i = 0; i < 2*NumInsertStates; i++) gapExtend[i] = gapExtend2Default[i];
    }
    else {
      cerr << "ERROR: No default initial distribution/parameter settings exist" << endl
           << "       for " << NumInsertStates << " pairs of insert states.  Use --paramfile." << endl;
      exit (1);
    }

    alphabet = alphabetDefault;

    for (int i = 0; i < (int) alphabet.length(); i++){
      emitSingle[(unsigned char) tolower(alphabet[i])] = emitSingleDefault[i];
      emitSingle[(unsigned char) toupper(alphabet[i])] = emitSingleDefault[i];
      for (int j = 0; j <= i; j++){
        emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) tolower(alphabet[j])] = emitPairsDefault[i][j];
        emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) toupper(alphabet[j])] = emitPairsDefault[i][j];
        emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) tolower(alphabet[j])] = emitPairsDefault[i][j];
        emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) toupper(alphabet[j])] = emitPairsDefault[i][j];
        emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) tolower(alphabet[i])] = emitPairsDefault[i][j];
        emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) toupper(alphabet[i])] = emitPairsDefault[i][j];
        emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) tolower(alphabet[i])] = emitPairsDefault[i][j];
        emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) toupper(alphabet[i])] = emitPairsDefault[i][j];
      }
    }
  }

  void Probcons::Impl::PrintParameters (const char* message) const
  {
   // print parameters to the screen
    cerr << message << endl
         << "    initDistrib[] = { ";
    for (int i = 0; i < NumMatrixTypes; i++) cerr << setprecision (10) << initDistrib[i] << " ";
    cerr << "}" << endl
         << "        gapOpen[] = { ";
    for (int i = 0; i < NumInsertStates*2; i++) cerr << setprecision (10) << gapOpen[i] << " ";
    cerr << "}" << endl
         << "      gapExtend[] = { ";
    for (int i = 0; i < NumInsertStates*2; i++) cerr << setprecision (10) << gapExtend[i] << " ";
    cerr << "}" << endl
         << endl;
  }
}
