/////////////////////////////////////////////////////////////////
// MultiSequence.cpp
/////////////////////////////////////////////////////////////////

#include "MultiSequence.hpp"

namespace CONTRALIGN {
/////////////////////////////////////////////////////////////////
// MultiSequence::MultiSequence()
//
// Default constructor.
/////////////////////////////////////////////////////////////////

MultiSequence::MultiSequence() {}

/////////////////////////////////////////////////////////////////
// MultiSequence::MultiSequence()
//
// Constructor.  Also load an MFA file, specified by
// "filename".  If "ignore_gaps" is true, then all gaps are
// removed from the sequences when they are read.
/////////////////////////////////////////////////////////////////

MultiSequence::MultiSequence(const std::string &filename, bool ignore_gaps)
{
    LoadMFA(filename, ignore_gaps);
}

/////////////////////////////////////////////////////////////////
// MultiSequence::MultiSequence()
//
// Copy constructor.
/////////////////////////////////////////////////////////////////

MultiSequence::MultiSequence(const MultiSequence &rhs)
{
    for (size_t i = 0; i < rhs.sequences.size(); i++)
        sequences.push_back(new Sequence(*rhs.sequences[i]));
    weights = rhs.weights;
}

/////////////////////////////////////////////////////////////////
// MultiSequence::MultiSequence()
//
// Copy constructor, which projects onto a subset of the 
// sequences.  In particular, suppose "rhs" contains a 
// multiple alignment of K sequences.  Then
// "subset" contains a list of indices, each between 0 and
// K-1.  This function returns an alignment which has
// been projected down to those particular sequences,
// by removing any columns with nothing but gaps.
/////////////////////////////////////////////////////////////////

MultiSequence::MultiSequence(const MultiSequence &rhs, const std::vector<int> &subset)
{
    Assert(subset.size() > 0, "Cannot project onto empty subset.");
    std::vector<int> column_occupied = rhs.ComputeOccupiedColumns(subset);
    
    // now perform the projections
    for (size_t i = 0; i < subset.size(); i++)
    {
        const Sequence &seq = *rhs.sequences[subset[i]];
        const std::string &seq_data = seq.GetData();
        
        // save only columns that are not completely gapped
        std::string projected_data = "@";
        projected_data.reserve(seq_data.size());
        for (size_t j = 1; j < column_occupied.size(); j++)
            if (column_occupied[j]) projected_data.push_back(seq_data[j]);
        
        // add sequence to return value
        sequences.push_back(new Sequence(projected_data, seq.GetName(), seq.GetID()));
    }

    weights = ComputePositionBasedSequenceWeights();
}

/////////////////////////////////////////////////////////////////
// MultiSequence::~MultiSequence()
//
// Destructor.
/////////////////////////////////////////////////////////////////

MultiSequence::~MultiSequence() 
{
    for (size_t i = 0; i < sequences.size(); i++)
        delete sequences[i];
}

/////////////////////////////////////////////////////////////////
// MultiSequence::operator=()
//
// Assignment operator.
/////////////////////////////////////////////////////////////////

const MultiSequence& MultiSequence::operator=(const MultiSequence &rhs) {
    if (this != &rhs)
    {
        for (size_t i = 0; i < sequences.size(); i++)
            delete sequences[i];
        sequences.clear();
        for (size_t i = 0; i < rhs.sequences.size(); i++)
            sequences.push_back(new Sequence(*rhs.sequences[i]));
        weights = rhs.weights;
    }
    return *this;
}

/////////////////////////////////////////////////////////////////
// MultiSequence::operator==()
//
// Equality operator.
/////////////////////////////////////////////////////////////////

bool MultiSequence::operator==(const MultiSequence &rhs) const
{
    if (this == &rhs) return true;
    
    if (sequences.size() != rhs.sequences.size()) return false;
    for (size_t i = 0; i < sequences.size(); i++)
        if (!(*sequences[i] == *rhs.sequences[i])) return false;
    if (weights != rhs.weights) return false;
    
    return true;
}

/////////////////////////////////////////////////////////////////
// MultiSequence::AddSequence()
//
// Add sequence to MultiSequence object.
/////////////////////////////////////////////////////////////////

void MultiSequence::AddSequence(Sequence *seq)
{
    sequences.push_back(seq);
    weights = ComputePositionBasedSequenceWeights();
}

/////////////////////////////////////////////////////////////////
// MultiSequence::Sort()
//
// Sort sequences by ID number.  Using a slow O(N^2)
// sorting algorithm here.
/////////////////////////////////////////////////////////////////

void MultiSequence::Sort()
{
    for (size_t i = 0; i < sequences.size(); i++)
    {
        for (size_t j = i+1; j < sequences.size(); j++)
        {
            if (sequences[i]->GetID() > sequences[j]->GetID())
            {
                std::swap(sequences[i], sequences[j]);
                std::swap(weights[i], weights[j]);
            }
        }
    }
}

/////////////////////////////////////////////////////////////////
// MultiSequence::LoadMFA()
//
// Load data from MFA file.
/////////////////////////////////////////////////////////////////

void MultiSequence::LoadMFA(const std::string &filename, bool ignore_gaps)
{
    std::ifstream infile(filename.c_str());
    if (infile.fail()) Error(("Unable to open input file: " + filename).c_str());
    LoadMFA(infile, ignore_gaps);
    infile.close();
}

/////////////////////////////////////////////////////////////////
// MultiSequence::LoadMFA()
//
// Load data from MFA file.
/////////////////////////////////////////////////////////////////

void MultiSequence::LoadMFA(std::istream &infile, bool ignore_gaps)
{ 
    std::string s, name, data;
    
    while (getline(infile, s))
    {
        if (s.length() == 0) continue;
        
        if (s[0] == '>')
        {
            if (data.length() != 0) sequences.push_back(new Sequence(data, name, int(sequences.size())));
            name = s.substr(1);
            data = "@";
        }
        else
        {
            if (data.length() == 0) Error("MFA sequence header should begin with '>'.");
            for (size_t i = 0; i < s.length(); i++)
            {
                char ch = s[i];
                if (isspace(ch)) continue;
                if (ch == '.') ch = '-';
#if RNA
                if (ch == 't') ch = 'u';
                if (ch == 'T') ch = 'U';
#endif
                if (ignore_gaps && ch == '-') continue;
                data.push_back (ch);
            }
        }
    }
    
    if (data.length() == 0)
        Error("MFA input file contained no sequences.");
    else
        sequences.push_back(new Sequence(data, name, int(sequences.size())));

    weights = ComputePositionBasedSequenceWeights();
}

/////////////////////////////////////////////////////////////////
// MultiSequence::WriteMFA()
//
// Write MFA to file.
/////////////////////////////////////////////////////////////////

void MultiSequence::WriteMFA(std::ostream &outfile, const int num_columns) const
{
    for (size_t i = 0; i < sequences.size(); i++)
    {
        outfile << ">" << sequences[i]->GetName() << std::endl;
        int length = sequences[i]->GetLength();
        const std::string &data = sequences[i]->GetData();
        
        for (int j = 1; j <= length; j++)
        {
            outfile << char(data[j]);
            if (j % num_columns == 0) outfile << std::endl;
        }
        if (length % num_columns != 0) outfile << std::endl;
    }
}

/////////////////////////////////////////////////////////////////
// MultiSequence::ComputeOccupiedColumns()
//
// Compute occupied columns in MultiSequence.
/////////////////////////////////////////////////////////////////

std::vector<int> MultiSequence::ComputeOccupiedColumns(const std::vector<int> &subset) const
{
    // column_occupied[i] is true whenever there is at least one 
    // letter in that column for sequences in subset 
    
    std::vector<int> column_occupied(GetLength(subset) + 1, 0);
    
    for (size_t i = 0; i < subset.size(); i++)
    {
        Assert(subset[i] >= 0 && subset[i] < int(sequences.size()), 
               "Reference to invalid sequence in MultiSequence object!");
        
        const std::string &seq_data = sequences[subset[i]]->GetData();
        for (size_t j = 1; j < column_occupied.size(); j++)
            column_occupied[j] = column_occupied[j] || (seq_data[j] != '-');
    }
    
    return column_occupied;
}

/////////////////////////////////////////////////////////////////
// MultiSequence::ComputeCoreBlockColumns()
//
// Identify columns which have no lowercase letters.
/////////////////////////////////////////////////////////////////

std::vector<int> MultiSequence::ComputeCoreBlockColumns(const std::vector<int> &subset) const 
{
    // is_core_block[i] is false whenever there is at least one 
    // lowercase letter in that column for sequences in subset 
    
    std::vector<int> is_core_block (GetLength(subset) + 1, 1);
    
    for (size_t i = 0; i < subset.size(); i++)
    {
        Assert(subset[i] >= 0 && subset[i] < int(sequences.size()), 
               "Reference to invalid sequence in MultiSequence object!");
        
        const std::string &seq_data = sequences[subset[i]]->GetData();
        for (size_t j = 1; j < is_core_block.size(); j++)
            if (islower(seq_data[j])) is_core_block[j] = 0;
    }
    
    return is_core_block;
}

/////////////////////////////////////////////////////////////////
// MultiSequence::GetConsensus()
//
// Convert multiple alignment into a consensus sequence.  The
// scheme used here is that 'X' is used whenever more than one
// type of amino acid appears in a column. 
/////////////////////////////////////////////////////////////////

std::string MultiSequence::GetConsensus() const
{
    std::string consensus = sequences[0]->GetData();
    for (size_t i = 1; i < sequences.size(); i++)
    {
        const std::string &seq_data = sequences[i]->GetData();
        for (size_t j = 1; j < seq_data.length(); j++)
        {
            if (consensus[j] == '-') 
                consensus[j] = toupper(seq_data[j]);
            else if (seq_data[j] == '-')
                continue;
            else if (seq_data[j] != consensus[j])
                consensus[j] = 'X';
        }
    }
    
    return consensus;
}

/////////////////////////////////////////////////////////////////
// MultiSequence::ComputeEditString()
//
// Project alignment onto two subsets and finds the edit string.
// relating them.
/////////////////////////////////////////////////////////////////

std::string MultiSequence::ComputeEditString(const std::vector<int> &subset1, const std::vector<int> &subset2) const
{
    Assert(subset1.size() > 0, "Cannot project onto empty subset.");
    Assert(subset2.size() > 0, "Cannot project onto empty subset.");
    
    const int length = GetLength();
    
    std::vector<int> column_occupied1 = ComputeOccupiedColumns(subset1);
    std::vector<int> column_occupied2 = ComputeOccupiedColumns(subset2);
    std::vector<int> is_core_block1 = ComputeCoreBlockColumns(subset1);
    std::vector<int> is_core_block2 = ComputeCoreBlockColumns(subset2);
    
    std::string ret = "@";
    ret.reserve (length+1);
    for (int i = 1; i <= length; i++)
    {
        char ch = '~';
        if (column_occupied1[i] && column_occupied2[i]) ch = 'm';
        else if (column_occupied1[i]) ch = 'x';
        else if (column_occupied2[i]) ch = 'y';
        
        if (ch != '~')
        {
            if (is_core_block1[i] && is_core_block2[i]) ch = toupper(ch);
            ret.push_back (ch);
        }    
    }

    return ret;
}

/////////////////////////////////////////////////////////////////
// MultiSequence::GetAnnotationChar()
//
// Return CLUSTALW annotation for column.
/////////////////////////////////////////////////////////////////

char MultiSequence::GetAnnotationChar(const std::vector<char> &column) const
{
    std::vector<int> counts(256, 0);
    char annot[] = ":::::::::...........";
    const char *groups[20] = { "STA", "NEQK", "NHQK", "NDEQ", "QHRK",
                               "MILV", "MILF", "HY", "FYW", "CSA", 
                               "ATV", "SAG", "STNK", "STPA", "SGND", 
                               "SNDEQK", "NDEQHK", "NEHQRK", "FVLIM", "HFY" };
    
    for (size_t i = 0; i < column.size(); i++)
        counts[BYTE(toupper(column[i]))]++;
    
    for (int i = 0; i < 256; i++)
        if (char(i) != '-' && counts[i] == int(column.size()))
            return '*';
    
    for (int i = 0; i < 20; i++)
    {
        int total_counts = 0;
        for (size_t j = 0; j < strlen(groups[i]); j++) 
            total_counts += counts[BYTE(groups[i][j])];
        if (total_counts == int(column.size()))
            return annot[i];
    }

    return ' ';
}

/////////////////////////////////////////////////////////////////
// MultiSequence::WriteCLUSTALW()
//
// Write MultiSequence in CLUSTALW format.
/////////////////////////////////////////////////////////////////

void MultiSequence::WriteCLUSTALW(std::ostream &outfile, const int num_columns) const
{  
    // print header
    
    outfile << "CLUSTALW (CentroidAlign) multiple sequence alignment" << std::endl;
    
    // compute length of longest name
    
    int name_len = 0;
    for (size_t i = 0; i < sequences.size(); i++)
        name_len = std::max(name_len, int(sequences[i]->GetName().length()));
    name_len += 4;
    
    int written_chars = 0;
    bool all_done = false;
    
    // continue writing characters, until finished
    
    while (!all_done)
    {
        outfile << std::endl;
        all_done = true;
        
        // write each sequence
        
        for (size_t i = 0; i < sequences.size(); i++)
        {
            // write sequence name
            
            outfile << sequences[i]->GetName();
            for (size_t j = name_len - sequences[i]->GetName().length(); j > 0; j--)
                outfile << ' ';
            
            const int len = sequences[i]->GetLength();
            const std::string &data = sequences[i]->GetData();
            
            // write characters in the sequence (as long as there are
            // characters left to be written)
            
            if (written_chars < len)
            {
                for (int j = 1; j <= num_columns; j++)
                {
                    if (written_chars + j <= len)
                        outfile << data[written_chars + j];
                    else
                        break;
                }
            }
            outfile << std::endl;
            
            if (written_chars + num_columns < len) all_done = false;
        }
        
        // write annotation
        
        for (int j = 0; j < name_len; j++) outfile << ' ';
        for (int j = 1; j <= num_columns; j++)
        {
            std::vector<char> column;
            for (size_t i = 0; i < sequences.size(); i++)
                if (written_chars + j <= sequences[i]->GetLength())
                    column.push_back(sequences[i]->GetData()[written_chars + j]);
            if (column.size() > 0) outfile << GetAnnotationChar (column);     
        }
        outfile << std::endl;

        written_chars += num_columns;
    }
}


/////////////////////////////////////////////////////////////////
// MultiSequence::IsAlignment()
//
// Check if all sequences have the same length.
/////////////////////////////////////////////////////////////////

bool MultiSequence::IsAlignment() const
{
    for (size_t i = 1; i < sequences.size(); i++) 
        if (sequences[i]->GetLength() != sequences[i-1]->GetLength())
            return false;
    return true;
}

/////////////////////////////////////////////////////////////////
// MultiSequence::GetAlignedTo()
//
// Compute pair of mappings (aligned_to_x, aligned_to_y) such
// that aligned_to_x[i] states the index of the position in y to
// which x[i] is aligned (and vice versa), or 0 if unaligned.
// This routine is only valid for pairwise alignments.
/////////////////////////////////////////////////////////////////

std::pair<std::vector<int>, std::vector<int> > MultiSequence::GetAlignedTo() const
{
    Assert(sequences.size() == 2, "Invalid number of sequences.");
    
    const std::string &x = sequences[0]->GetData();
    const std::string &y = sequences[1]->GetData();
    Assert(x.length() == y.length(), "Length mismatch.");

#if !FORCE_UNIQUE_PARSES
    std::vector<int> aligned_to_x(sequences[0]->GetCompressedLength() + 1, Sequence::UNKNOWN);
    std::vector<int> aligned_to_y(sequences[1]->GetCompressedLength() + 1, Sequence::UNKNOWN);
#else
    std::vector<int> aligned_to_x(sequences[0]->GetCompressedLength() + 1, Sequence::UNALIGNED);
    std::vector<int> aligned_to_y(sequences[1]->GetCompressedLength() + 1, Sequence::UNALIGNED);
#endif
    
    int i = 0;
    int j = 0;
    for (size_t k = 1; k < x.length(); k++)
    {
        if (isalpha(x[k])) i++;
        if (isalpha(y[k])) j++;
        if (isalpha(x[k]) && isalpha(y[k]) &&
            isupper(x[k]) && isupper(y[k]))
        {
            Assert(i >= 1 && i < int(aligned_to_x.size()), "Array access out-of-bounds.");
            Assert(j >= 1 && j < int(aligned_to_y.size()), "Array access out-of-bounds.");
            aligned_to_x[i] = j;
            aligned_to_y[j] = i;
        }
    }
    
    return std::make_pair(aligned_to_x, aligned_to_y);
}

/////////////////////////////////////////////////////////////////
// MultiSequence::GetNumSequences()
//
// Retrieve number of sequences.
/////////////////////////////////////////////////////////////////

int MultiSequence::GetNumSequences() const
{ 
    return sequences.size();
}
  
/////////////////////////////////////////////////////////////////
// MultiSequence::GetLength()
//
// Retrieve alignment length.
/////////////////////////////////////////////////////////////////

int MultiSequence::GetLength() const
{
    Assert(sequences.size() > 0, "MultiSequence must have at least one sequence to retrieve length.");
    const int length = sequences[0]->GetLength();
#ifndef NDEBUG
    for (int i = 1; i < (int) sequences.size(); i++)
        Assert(sequences[i]->GetLength() == length, "Not all sequences have the same length.");
#endif
    return length;
}

/////////////////////////////////////////////////////////////////
// MultiSequence::GetLength()
//
// Retrieve length of sequences in "subset".
/////////////////////////////////////////////////////////////////

int MultiSequence::GetLength(const std::vector<int> &subset) const
{
    Assert(sequences.size() > 0, "MultiSequence must have at least one sequence to retrieve length.");
    const int length = sequences[subset[0]]->GetLength();
#ifndef NDEBUG
    for (int i = 1; i < (int) subset.size(); i++)
        Assert(sequences[subset[i]]->GetLength() == length, "Not all sequences have the same length.");
#endif
    return length;
}

/////////////////////////////////////////////////////////////////
// MultiSequence::GetSequence()
//
// Retrieve sequence.
/////////////////////////////////////////////////////////////////

const Sequence &MultiSequence::GetSequence(int index) const
{
    Assert(0 <= index && index < (int) sequences.size(), "Invalid sequence index.");
    return *sequences[index];
}

//////////////////////////////////////////////////////////////////////
// MultiSequence::ComputePositionBasedSequenceWeights()
//
// Compute sequence weights according to:
//    Henikoff, S., and Henikoff J.  1994.  Position-based
//    sequence weights.  J Mol Biol 243(4):574-578.
//////////////////////////////////////////////////////////////////////

std::vector<double> MultiSequence::ComputePositionBasedSequenceWeights() const
{
    std::vector<double> weights(sequences.size(), 0.0);
    std::vector<int> counts(256);

    for (int i = 1; i <= sequences[0]->GetLength(); i++)
    {
        int diversity = 0;
        std::fill(counts.begin(), counts.end(), 0);
        
        for (size_t j = 0; j < sequences.size(); j++)
        {
            if (counts[BYTE(sequences[j]->GetData()[i])] == 0) diversity++;
            ++(counts[BYTE(sequences[j]->GetData()[i])]);
        }
        
        for (size_t j = 0; j < sequences.size(); j++)
            weights[j] += 1.0 / (diversity * counts[BYTE(sequences[j]->GetData()[i])]);
    }

    weights /= Sum(weights);
    return weights;
}

//////////////////////////////////////////////////////////////////////
// MultiSequence::GetPairWeight()
//
// Get weight for a pair of sequences.
//////////////////////////////////////////////////////////////////////

double MultiSequence::GetPairWeight(int i, int j) const
{
    Assert(0 <= i && i < (int) sequences.size(), "Invalid sequence index.");
    Assert(0 <= j && j < (int) sequences.size(), "Invalid sequence index.");

    double total_weight = (1.0 - Sum(weights * weights)) / 2.0;
    return weights[i] * weights[j] / total_weight;
}
}
