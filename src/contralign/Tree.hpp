/////////////////////////////////////////////////////////////////
// Tree.hpp
//
// Routines for constructing UPGMA trees.
/////////////////////////////////////////////////////////////////

#ifndef TREE_HPP
#define TREE_HPP

#include "Utilities.hpp"

/////////////////////////////////////////////////////////////////
// class Tree
/////////////////////////////////////////////////////////////////
namespace CONTRALIGN{
template <class RealT>
class Tree
{

public:
    
    /////////////////////////////////////////////////////////////////
    // struct TreeNode
    /////////////////////////////////////////////////////////////////
    
    struct TreeNode
    {
        std::vector<int> ids;
        
        RealT left_dist;
        RealT right_dist;
        TreeNode *left_child;
        TreeNode *right_child;
        
        TreeNode();
        TreeNode(const TreeNode &rhs);
        TreeNode &operator=(const TreeNode &rhs);
        virtual ~TreeNode();
        void Print(std::ostream &outfile, bool root = true) const;
    };
    
    static TreeNode *UPGMA(std::vector<std::vector<RealT> > distances);
};
}
#include "Tree.ipp"
    
#endif
