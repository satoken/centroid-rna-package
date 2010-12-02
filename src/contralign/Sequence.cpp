/////////////////////////////////////////////////////////////////
// Sequence.cpp
/////////////////////////////////////////////////////////////////

#include "Sequence.hpp"

namespace CONTRALIGN {
const int Sequence::UNKNOWN = -1;
const int Sequence::UNALIGNED = 0;

/////////////////////////////////////////////////////////////////
// Sequence::Sequence()
//
// Constructor.
/////////////////////////////////////////////////////////////////

Sequence::Sequence(const std::string &data, const std::string &name, int id) :
    data(data), name(name), id(id)
{}

/////////////////////////////////////////////////////////////////
// Sequence::Sequence()
//
// Copy constructor.
/////////////////////////////////////////////////////////////////

Sequence::Sequence(const Sequence &rhs) : 
    data(rhs.data), name(rhs.name), id(rhs.id)
{}

/////////////////////////////////////////////////////////////////
// Sequence::operator=()
//
// Assignment operator.
/////////////////////////////////////////////////////////////////

const Sequence& Sequence::operator=(const Sequence &rhs)
{
    if (this != &rhs)
    {
        data = rhs.data;
        name = rhs.name;
        id = rhs.id;
    }
    return *this;
}

/////////////////////////////////////////////////////////////////
// Sequence::operator==()
//
// Equality operator.
/////////////////////////////////////////////////////////////////

bool Sequence::operator==(const Sequence &rhs) const 
{
    if (this == &rhs) return true;
    
    if (data != rhs.data) return false;
    if (name != rhs.name) return false;
    if (id != rhs.id) return false;
    
    return true;
}

/////////////////////////////////////////////////////////////////
// Sequence::Sequence()
//
// Retrieve a substring of the current sequence, from 
// positions "start" to "end", inclusive.
/////////////////////////////////////////////////////////////////

Sequence::Sequence(const Sequence &rhs, ExtractPositions, int start, int end) :
    data("@"), name(rhs.name), id(rhs.id)
{
    Assert(1 <= start && start <= end && end <= rhs.GetLength(), "Invalid positions.");
    
    const int new_length = std::max(end - start + 1, 0);
    data += rhs.data.substr(start, new_length);
}

/////////////////////////////////////////////////////////////////
// Sequence::Sequence()
//
// Return a sequence without gaps.
/////////////////////////////////////////////////////////////////

Sequence::Sequence(const Sequence &rhs, CompressGaps) :
  data("@"), name(rhs.name), id(rhs.id) 
{
    for (size_t i = 1; i < rhs.data.length(); i++)
        if (rhs.data[i] != '-') data.push_back(rhs.data[i]);
}

/////////////////////////////////////////////////////////////////
// Sequence::Sequence()
//
// Insert gaps into a sequence.  Takes in an edit string,
// such as @XXXYYXYXYYMMMM, and returns a sequence with the
// appropriate positions gapped.  Leave existing gaps in
// the sequence.
/////////////////////////////////////////////////////////////////

Sequence::Sequence(const Sequence &rhs, InsertGaps, const std::string &edit_string, char ch) :
    data("@"), name(rhs.name), id(rhs.id)
{
    size_t j = 1;
    ch = toupper(ch);
    for (size_t i = 1; i < edit_string.length(); i++)
    {
        char this_ch = toupper(edit_string[i]);
        if (this_ch == ch || this_ch == 'M')
        {
            Assert(j < rhs.data.length(), "Sequence does not match edit string.");
            data.push_back (rhs.data[j++]);
        }
        else 
            data.push_back ('-');
    }
}

/////////////////////////////////////////////////////////////////
// Sequence::ComputeMapping()
//
// Return a vector called "mapping" such that
// (1) mapping[0] = 0
// (2) mapping[i] = the index of the ith non-gapped letter 
//                  in the sequence
/////////////////////////////////////////////////////////////////

std::vector<int> Sequence::ComputeMapping() const
{
    std::vector<int> ret;
    for (int i = 0; i < int(data.length()); i++)
        if (i == 0 || data[i] != '-')
            ret.push_back(i);
    return ret;
}

/////////////////////////////////////////////////////////////////
// Sequence::GetCompressedLength()
//
// Return length of sequence without gaps.
/////////////////////////////////////////////////////////////////

int Sequence::GetCompressedLength() const
{
    int ret = 0;
    for (int i = 1; i < int(data.length()); i++) 
        if (isalpha(data[i])) ret++;
    return ret;
}
}
