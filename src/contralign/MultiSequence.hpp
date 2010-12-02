/////////////////////////////////////////////////////////////////
// MultiSequence.hpp
/////////////////////////////////////////////////////////////////

#ifndef MULTISEQUENCE_HPP
#define MULTISEQUENCE_HPP

#include <cstdio>
#include <cctype>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include "Config.hpp"
#include "Utilities.hpp"
#include "Sequence.hpp"

/////////////////////////////////////////////////////////////////
// class MultiSequence
/////////////////////////////////////////////////////////////////
namespace CONTRALIGN {
class MultiSequence {
    std::vector<Sequence *> sequences;
    std::vector<double> weights;
    
    // compute position-based sequence weights
    std::vector<double> ComputePositionBasedSequenceWeights() const;

public:

    // default constructor
    MultiSequence();

    // constructor for reading either alignments or unaligned sequences
    MultiSequence(const std::string &filename, bool ignore_gaps);

    // copy constructor
    MultiSequence(const MultiSequence &rhs);

    // copy constructor, with projection to a subset of sequences
    MultiSequence(const MultiSequence &rhs, const std::vector<int> &subset);

    // destructor
    ~MultiSequence();

    // assignment operator
    const MultiSequence &operator=(const MultiSequence &rhs);

    // equality comparator
    bool operator==(const MultiSequence &rhs) const;

    // add a pre-allocated sequence to the multisequence
    void AddSequence(Sequence *seq);

    // sort sequences according to their saved indices
    void Sort();

    // load alignment or unaligned sequences from file
    void LoadMFA(const std::string &filename, bool ignore_gaps);

    // load alignment or unaligned sequences from file
    void LoadMFA(std::istream &infile, bool ignore_gaps);

    // write alignment or unaligned sequences to a file
    void WriteMFA(std::ostream &outfile, const int num_columns = 1000000000) const;

    // return a vector<int> indicating which columns are occupied by some sequence in subset
    std::vector<int> ComputeOccupiedColumns(const std::vector<int> &subset) const;

    // return a vector<int> indicating which columns are core blocks in subset
    std::vector<int> ComputeCoreBlockColumns(const std::vector<int> &subset) const;

    // retrieve consensus sequence for aligned sequences
    std::string GetConsensus() const;

    // compute edit string of (mxyMXY)* for the alignment
    std::string ComputeEditString(const std::vector<int> &subset1, const std::vector<int> &subset2) const;

    // retrieve CLUSTALW annotation character for a column of characters
    char GetAnnotationChar(const std::vector<char> &column) const;

    // write alignment in CLUSTALW format
    void WriteCLUSTALW(std::ostream &outfile, const int num_columns = 60) const;

    // return true if alignment
    bool IsAlignment() const;

    // return number of sequences
    int GetNumSequences() const;

    // return length of alignment (all sequences must have equal length)
    int GetLength() const;

    // return length of alignment of some subset (all sequences in subset must have equal length)
    int GetLength(const std::vector<int> &subset) const;

    // return a particular sequence
    const Sequence &GetSequence(int index) const;

    // return the position to which each letter is aligned (only for pairwise alignments)
    std::pair<std::vector<int>, std::vector<int> > GetAlignedTo() const;

    // return weight for a pair of sequences
    double GetPairWeight(int i, int j) const;
};
}
#endif
