/////////////////////////////////////////////////////////////////
// Progressive.ipp
/////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
// Progressive::Progressive()
//
// Constructor.
/////////////////////////////////////////////////////////////////

template<class RealT>
Progressive<RealT>::Progressive(const MultiSequence &multi_seqs,
                                const std::vector<SparseMatrix<RealT> *> &posteriors,
                                const bool toggle_verbose) :
    multi_seqs(multi_seqs), posteriors(posteriors), toggle_verbose(toggle_verbose),
    root(NULL), alignment(NULL)
{}

/////////////////////////////////////////////////////////////////
// Progressive::~Progressive()
//
// Destructor.
/////////////////////////////////////////////////////////////////

template<class RealT>
Progressive<RealT>::~Progressive()
{
    delete root;
    delete alignment;
}

/////////////////////////////////////////////////////////////////
// Progressive::BuildTree()
//
// Compute guide tree.
/////////////////////////////////////////////////////////////////

template<class RealT>
void Progressive<RealT>::BuildTree()
{
    // first, compute pairwise distances

    const int m = multi_seqs.GetNumSequences();
    std::vector<std::vector<RealT> > distances(m, std::vector<RealT>(m));

    for (int i = 0; i < m; i++)
    {
        distances[i][i] = 0;
        for (int j = i+1; j < m; j++)
        {
            distances[i][j] = posteriors[i*m+j]->GetSum();
        }
    }

    // second, construct tree
    
    root = Tree<RealT>::UPGMA(distances);
}	

/////////////////////////////////////////////////////////////////
// Progressive::DoAlignment()
//
// Perform progressive alignment starting from root.
/////////////////////////////////////////////////////////////////

template<class RealT>
void Progressive<RealT>::DoAlignment()
{
    BuildTree();
    if (toggle_verbose)
    {
        root->Print(std::cerr);
    }
    alignment = AlignSubtree(root);
    alignment->Sort();
}

/////////////////////////////////////////////////////////////////
// Progressive::AlignSubtree()
//
// Perform progressive alignment for a subtree.
/////////////////////////////////////////////////////////////////

template<class RealT>
MultiSequence *Progressive<RealT>::AlignSubtree(typename Tree<RealT>::TreeNode *node)
{
    if (!node->left_child && !node->right_child)
        return new MultiSequence(multi_seqs, node->ids);

    MultiSequence *left = AlignSubtree(node->left_child);
    MultiSequence *right = AlignSubtree(node->right_child);
    MultiSequence *ret = AlignGroups(*left, *right);    
    delete left;
    delete right;    
    return ret;
}

/////////////////////////////////////////////////////////////////
// Progressive::AlignGroups()
//
// Align two groups of sequences.
/////////////////////////////////////////////////////////////////

template<class RealT>
MultiSequence *Progressive<RealT>::AlignGroups(const MultiSequence &left, const MultiSequence &right)
{
    RealT *scoring = ComputeScoringMatrix(left, right);
    const std::string edit_string = ComputeAlignment(scoring, left.GetLength(), right.GetLength());
    delete [] scoring;

    MultiSequence *ret = new MultiSequence();
    for (int i = 0; i < left.GetNumSequences(); i++)
        ret->AddSequence(new Sequence(left.GetSequence(i), Sequence::INSERT_GAPS, edit_string, 'X'));
    for (int i = 0; i < right.GetNumSequences(); i++)
        ret->AddSequence(new Sequence(right.GetSequence(i), Sequence::INSERT_GAPS, edit_string, 'Y'));

    return ret;    
}

/////////////////////////////////////////////////////////////////
// Progressive::ComputeScoringMatrix()
//
// Compute scoring matrix by stacking posterior matrices.
/////////////////////////////////////////////////////////////////

template<class RealT>
RealT *Progressive<RealT>::ComputeScoringMatrix(const MultiSequence &left, const MultiSequence &right)
{
    const int m = multi_seqs.GetNumSequences();
    const int LX = left.GetLength();
    const int LY = right.GetLength();
    const int SIZE = (LX+1)*(LY+1);

    RealT *ret = new RealT[SIZE];
    std::fill(ret, ret + SIZE, RealT(0));

    for (int lc = 0; lc < left.GetNumSequences(); lc++)
    {
        for (int rc = 0; rc < right.GetNumSequences(); rc++)
        {
            const int l = left.GetSequence(lc).GetID();
            const int r = right.GetSequence(rc).GetID();
            const std::vector<int> l_positions = GetSequencePositions(left.GetSequence(lc).GetData());
            const std::vector<int> r_positions = GetSequencePositions(right.GetSequence(rc).GetData());

            if (l < r)
            {
                const SparseMatrix<RealT> &sparse = *posteriors[l*m+r];
                for (int i = 1; i < sparse.GetNumRows(); i++)
                {
                    for (const SparseMatrixEntry<RealT> *iter = sparse.GetRowBegin(i); iter != sparse.GetRowEnd(i); ++iter)
                        ret[l_positions[i] * (LY+1) + r_positions[iter->column]] += iter->value;
                }
            }
            else
            {
                const SparseMatrix<RealT> &sparse = *posteriors[r*m+l];
                for (int j = 1; j < sparse.GetNumRows(); j++)
                {
                    for (const SparseMatrixEntry<RealT> *iter = sparse.GetRowBegin(j); iter != sparse.GetRowEnd(j); ++iter)
                        ret[l_positions[iter->column] * (LY+1) + r_positions[j]] += iter->value;
                }
            }
        }
    }

    return ret;
}

/////////////////////////////////////////////////////////////////
// Progressive::ComputeAlignment()
//
// Compute alignment given a scoring matrix.
/////////////////////////////////////////////////////////////////

template<class RealT>
std::string Progressive<RealT>::ComputeAlignment(const RealT *scoring, int LX, int LY)
{
    const int SIZE = (LX+1)*(LY+1);
    RealT *best = new RealT[SIZE];
    int *traceback = new int[SIZE];

    // dynamic programming

    for (int i = 0; i <= LX; i++)
    {
        for (int j = 0; j <= LY; j++)
        {
            int ij = i*(LY+1)+j;
            int i1j = ij-LY-1;
            int ij1 = ij-1;
            int i1j1 = ij-LY-2;

            best[ij] = (ij == 0 ? RealT(0) : RealT(NEG_INF));
            traceback[ij] = -1;

            if (i > 0 && j > 0 && best[ij] < best[i1j1] + scoring[ij])
            {
                best[ij] = best[i1j1] + scoring[ij];
                traceback[ij] = 0;
            }

            if (i > 0 && best[ij] < best[i1j])
            {
                best[ij] = best[i1j];
                traceback[ij] = 1;
            }

            if (j > 0 && best[ij] < best[ij1])
            {
                best[ij] = best[ij1];
                traceback[ij] = 2;
            }
        }
    }

    delete [] best;

    // traceback

    std::string edit_string;
    int ci = LX;
    int cj = LY;
    
    while (ci > 0 || cj > 0)
    {
        switch (traceback[ci*(LY+1)+cj])
        {
            case 0: edit_string.push_back('M'); ci--; cj--; break;
            case 1: edit_string.push_back('X'); ci--; break;
            case 2: edit_string.push_back('Y'); cj--; break;
            default: Error("Should not get here."); break;
        }
    }

    delete [] traceback;

    edit_string.push_back('@');
    std::reverse(edit_string.begin(), edit_string.end());

    return edit_string;
}

/////////////////////////////////////////////////////////////////
// Progressive::GetAlignment()
//
// Get final alignment.
/////////////////////////////////////////////////////////////////

template<class RealT>
const MultiSequence &Progressive<RealT>::GetAlignment()
{
    return *alignment;
}
