/////////////////////////////////////////////////////////////////
// Tree.cpp
//
// Routines for constructing UPGMA trees.
/////////////////////////////////////////////////////////////////

#include "Tree.hpp"

/////////////////////////////////////////////////////////////////
// Tree::TreeNode::Print()
//
// Print tree recursively.
/////////////////////////////////////////////////////////////////
namespace CONTRALIGN {
void Tree::TreeNode::Print(std::ostream &outfile, bool root) const
{
    if (ids.size() == 1)
    {
        outfile << ids[0];
    }
    else
    {
        outfile << "(";
        left_child->Print(outfile, false);
        if (left_child->ids.size() == 1 || right_child->ids.size() == 1) outfile << " ";
        right_child->Print(outfile, false);
        outfile << ")";
    }
    if (root) outfile << std::endl;    
}

/////////////////////////////////////////////////////////////////
// Tree::TreeNode::TreeNode()
//
// Default constructor.
/////////////////////////////////////////////////////////////////

Tree::TreeNode::TreeNode() :
    ids(0), left_dist(0), right_dist(0), left_child(NULL), right_child(NULL)
{}

/////////////////////////////////////////////////////////////////
// Tree::TreeNode::TreeNode()
//
// Copy constructor.
/////////////////////////////////////////////////////////////////

Tree::TreeNode::TreeNode(const TreeNode &rhs) :
    ids(rhs.ids), left_dist(rhs.left_dist), right_dist(rhs.right_dist),
    left_child(rhs.left_child), right_child(rhs.right_child)
{}

/////////////////////////////////////////////////////////////////
// Tree::TreeNode::operator=()
//
// Assignment operator.
/////////////////////////////////////////////////////////////////

Tree::TreeNode &Tree::TreeNode::operator=(const TreeNode &rhs)
{
    if (this != &rhs)
    {
        ids = rhs.ids;
        left_dist = rhs.left_dist;
        right_dist = rhs.right_dist;
        left_child = rhs.left_child;
        right_child = rhs.right_child;
    }
    return *this;
}

/////////////////////////////////////////////////////////////////
// Tree::TreeNode::TreeNode()
//
// Default constructor.
/////////////////////////////////////////////////////////////////

Tree::TreeNode::~TreeNode()
{
    delete left_child;
    delete right_child;
}

/////////////////////////////////////////////////////////////////
// Tree::UPGMA()
//
// Build tree via UPGMA.
/////////////////////////////////////////////////////////////////

Tree::TreeNode *Tree::UPGMA(std::vector<std::vector<Real> > distances)
{
    
    // set diagonal to infinite distance
    
    int n = int(distances.size());
    for (int i = 0; i < n; i++)
        distances[i][i] = Real(1e10);
    
    // construct initial forest
    
    std::vector<Tree::TreeNode *> trees (n);
    for (int i = 0; i < n; i++)
    {
        trees[i] = new Tree::TreeNode();
        trees[i]->ids = std::vector<int>(1, i);
    }
    
    // perform n-1 merges
    
    for (int k = 0; k < n-1; k++)
    {
        // find nearest neighbors (currently O(N^3) implementation)
        
        int bi = 0, bj = 0;
        for (int i = 0; i < n; i++)
        {
            for (int j = i+1; j < n; j++)
            {
                if (distances[i][j] < distances[bi][bj])
                {
                    bi = i;
                    bj = j;
                }
            }
        }

        // merge trees at slot bi
        
        Tree::TreeNode *merged = new Tree::TreeNode();
        merged->ids = Concatenate (trees[bi]->ids, trees[bj]->ids);
        merged->left_dist = distances[bi][bj] / Real(2.0);
        merged->right_dist = distances[bi][bj] / Real(2.0);
        merged->left_child = trees[bi];
        merged->right_child = trees[bj];
        
        trees[bi] = merged;
        trees[bj] = NULL;
        
        // update distance matrix
        
        distances[bi][bj] = Real(1e10);
        for (int m = 0; m < n; m++) if (m != bi && m != bj)
        {
            distances[m][bi] = distances[bi][m] =
                Real(0.9) * std::min(distances[m][bi], distances[m][bj]) +
                Real(0.1) * (distances[m][bi] + distances[m][bj]) / Real(2.0);
            distances[m][bj] = distances[bj][m] = Real(1e10);
        }
    }
    
    return trees[0];
}

}
