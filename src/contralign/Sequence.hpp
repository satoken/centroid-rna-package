/////////////////////////////////////////////////////////////////
// Sequence.hpp
/////////////////////////////////////////////////////////////////

#ifndef SEQUENCE_HPP
#define SEQUENCE_HPP

#include <cstdlib>
#include <cstring>
#include "Utilities.hpp"

/////////////////////////////////////////////////////////////////
// class Sequence
/////////////////////////////////////////////////////////////////
namespace CONTRALIGN {
class Sequence {
    std::string data;
    std::string name;
    int id;

public:  
    enum ExtractPositions { EXTRACT_POSITIONS };
    enum CompressGaps { COMPRESS_GAPS };
    enum InsertGaps { INSERT_GAPS };

    static const int UNKNOWN;
    static const int UNALIGNED;
    
    Sequence(const std::string &data, const std::string &name, int id);
    Sequence(const Sequence &rhs);
    const Sequence &operator=(const Sequence &rhs);
    bool operator==(const Sequence &rhs) const;
    
    Sequence(const Sequence &rhs, ExtractPositions, int start, int end);
    Sequence(const Sequence &rhs, CompressGaps);
    Sequence(const Sequence &rhs, InsertGaps, const std::string &edit_string, char ch);
    
    std::vector<int> ComputeMapping() const;
    
    // getters
    const std::string &GetData() const { return data; }
    const std::string &GetName() const { return name; }
    int GetLength() const { return data.size()-1; }
    int GetCompressedLength() const;
    int GetID() const { return id; } 
};
}
#endif
