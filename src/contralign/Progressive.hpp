/////////////////////////////////////////////////////////////////
// Progressive.hpp
/////////////////////////////////////////////////////////////////

#ifndef PROGRESSIVE_HPP
#define PROGRESSIVE_HPP

#include "MultiSequence.hpp"
#include "SparseMatrix.hpp"
#include "Utilities.hpp"
#include "Tree.hpp"

/////////////////////////////////////////////////////////////////
// class Progressive
/////////////////////////////////////////////////////////////////
namespace CONTRALIGN{
template<class RealT>
class Progressive
{
    const MultiSequence &multi_seqs;
    const std::vector<SparseMatrix<RealT> *> &posteriors;
    const bool toggle_verbose;

    typename Tree<RealT>::TreeNode *root;
    MultiSequence *alignment;

    void BuildTree();
    MultiSequence *AlignSubtree(typename Tree<RealT>::TreeNode *node);
    MultiSequence *AlignGroups(const MultiSequence &left, const MultiSequence &right);
    RealT *ComputeScoringMatrix(const MultiSequence &left, const MultiSequence &right);
    std::string ComputeAlignment(const RealT *scoring, int LX, int LY);

    
public:
    Progressive(const MultiSequence &multi_seqs,
                const std::vector<SparseMatrix<RealT> *> &posteriors,
                const bool toggle_verbose);
    ~Progressive();

    void DoAlignment();
    
    const MultiSequence &GetAlignment();
};

#include "Progressive.ipp"
}
#endif
