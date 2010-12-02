/////////////////////////////////////////////////////////////////
// InferenceEngine.hpp
/////////////////////////////////////////////////////////////////

#ifndef INFERENCEENGINE_HPP
#define INFERENCEENGINE_HPP

#include <queue>
#include <vector>
#include <string>
#include "Config.hpp"
#include "MultiSequence.hpp"
#include "ParameterManager.hpp"
#include "Utilities.hpp"
#include "LogSpace.hpp"

/////////////////////////////////////////////////////////////////
// class InferenceEngine
/////////////////////////////////////////////////////////////////
namespace CONTRALIGN {
template<class RealT>
class InferenceEngine
{
    enum STATE_TYPE
    {
        MATCH,
        INS_X,
        INS_Y,
#if PARAMS_DOUBLE_AFFINE
        INS2_X,
        INS2_Y,
#endif
        K
    };

    ParameterManager<RealT> *parameter_manager;
    unsigned char char_mapping[256];
#if PARAMS_HYDROPATHY
    bool is_hydrophilic[256];
#endif
    RealT per_position_loss;

    // dimensions   
    int LX, LY, SIZE;

    // sequence data
    std::vector<int> x, y;
    std::vector<int> aligned_to_x, aligned_to_y;
    std::vector<int> true_aligned_to_x, true_aligned_to_y;
    
#if PARAMS_HYDROPATHY
    std::vector<int> hydrophilic_count_x, hydrophilic_count_y;
#endif
    
    // dynamic programming matrices
    std::vector<int> Ft[K];                  // traceback
    std::vector<RealT> Fv[K];                // Viterbi
    std::vector<RealT> Ff[K];                // forward
    std::vector<RealT> Fb[K];                // backward
    std::vector<RealT> posterior;
    
    // parameters
    
#if PARAMS_MATCH
    std::pair<RealT, RealT> score_match[M+1][M+1];
#if PARAMS_COMPRESSED
    std::pair<RealT, RealT> score_compressed_match[M+1][M+1];
#endif
#endif
#if PARAMS_INSERT
    std::pair<RealT, RealT> score_insert[M+1];
#if PARAMS_COMPRESSED
    std::pair<RealT, RealT> score_compressed_insert[M+1];
#endif
#endif
#if PARAMS_SINGLE
    std::pair<RealT, RealT> score_single[K];
#endif
#if PARAMS_PAIR
    std::pair<RealT, RealT> score_pair[K][K];
#endif
#if PARAMS_HYDROPATHY
    std::pair<RealT, RealT> score_hydropathy[K][K][H+1];
#endif
#if PARAMS_TERMINAL_INSERTS
    std::pair<RealT, RealT> score_terminal_insert_boundary;
    std::pair<RealT, RealT> score_terminal_insert;
#endif
    
    void FillScores(typename std::vector<std::pair<RealT, RealT> >::iterator begin, typename std::vector<std::pair<RealT, RealT> >::iterator end, RealT value);
    void FillCounts(typename std::vector<std::pair<RealT, RealT> >::iterator begin, typename std::vector<std::pair<RealT, RealT> >::iterator end, RealT value);

    RealT ScoreMatch(int i, int j, int s) const;
    RealT ScoreInsertX(int i, int j, int s) const;
    RealT ScoreInsertY(int i, int j, int s) const;
    RealT ScoreInsert2X(int i, int j, int s) const;
    RealT ScoreInsert2Y(int i, int j, int s) const;
    
    void CountMatch(int i, int j, int s, RealT value);
    void CountInsertX(int i, int j, int s, RealT value);
    void CountInsertY(int i, int j, int s, RealT value);
    void CountInsert2X(int i, int j, int s, RealT value);
    void CountInsert2Y(int i, int j, int s, RealT value);

    void ClearCounts();
    std::vector<RealT> GetCounts();

    /*
    void InitializeParameters();
    void InitializeCounts();
    
    void TransformParameters();
    void TransformCounts();
    void WriteCounts(std::vector<RealT> &counts);
    */
    
public:

    // constructor and destructor
    InferenceEngine();
    ~InferenceEngine();		   

    // register params with the parameter manager
    void RegisterParameters(ParameterManager<RealT> &parameter_manager);
                            
    // load sequence
    void LoadSequences(const MultiSequence &seqs);
    
    // load parameter values                        
    void LoadValues(const std::vector<RealT> &values);
    
    // load loss function
    void UseLoss(const std::pair<std::vector<int>, std::vector<int> > &true_aligned_to, RealT example_loss);

    // use constraints
    void UseConstraints(const std::pair<std::vector<int>, std::vector<int> > &true_mapping);

    // Viterbi inference
    void ComputeViterbi();
    RealT GetViterbiScore() const;
    std::string PredictAlignmentViterbi();
    std::vector<RealT> ComputeViterbiFeatureCounts();

    // MEA inference
    void ComputeForward();
    RealT ComputeLogPartitionCoefficient() const;
    RealT ComputeForwardLogPartitionCoefficient() const;
    void ComputeBackward();
    RealT ComputeBackwardLogPartitionCoefficient() const;
    std::vector<RealT> ComputeFeatureCountExpectations();
    void ComputePosterior();
    std::string PredictAlignmentPosterior(const RealT gamma);
    std::string PredictAlignmentPosterior(const RealT gamma, const RealT *posterior);
    RealT *GetPosterior(const RealT posterior_cutoff) const;
    RealT *GetPosterior(const RealT posterior_cutoff, std::vector<RealT>& p) const;
};
}
#include "InferenceEngine.ipp"

#endif
