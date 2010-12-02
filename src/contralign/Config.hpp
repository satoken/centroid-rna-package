//////////////////////////////////////////////////////////////////////
// Config.hpp
//
// Global configuration file.
//////////////////////////////////////////////////////////////////////

#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>

#define COMMENT                                    0
namespace CONTRALIGN {
//////////////////////////////////////////////////////////////////////
// Miscellaneous options
//////////////////////////////////////////////////////////////////////

// upper bound on the number of logical parameters in the model; the
// program will fail to operate properly if this value is set too low
const int SHARED_PARAMETER_SIZE = 5000;

//////////////////////////////////////////////////////////////////////
// Options related to general inference
//////////////////////////////////////////////////////////////////////

// showing timings for inference routines
#define SHOW_TIMINGS                               0

//////////////////////////////////////////////////////////////////////
// Options related to training mode configuration
//////////////////////////////////////////////////////////////////////

#define STOCHASTIC_GRADIENT                        0
#define FORCE_UNIQUE_PARSES                        0

//////////////////////////////////////////////////////////////////////
// (A) Options related to max-margin training
//////////////////////////////////////////////////////////////////////

// the maximum loss DELTA(y,y') allocated to each training example; if
// this symbol is undefined, then the DELTA(y,y') loss function is not
// included
//    -- for a straight CRF, this value should be undefined
//    -- for a max-margin model, this value should be set to 1
// #define HAMMING_LOSS                            1

// multiplier used in the iterative convex-concave procedure (CCCP) for
// improving the solution of a max-margin model
//    -- for a regular max-margin model, this value should be set to 0
//    -- for a standard nonconvex model, this value should be set to 1
const double NONCONVEX_MULTIPLIER = 0.0;

// number of steps of CCCP; there is no need to change this to 1 in the
// case that NONCONVEX_MULTIPLIER == 0, as the code will detect this
// and abort after the first CCCP iteration by default
const int NUM_CCCP_STEPS = 5;

// use smooth approximation of max-margin algorithm for inference
// during training
#define SMOOTH_MAX_MARGIN                          0

//////////////////////////////////////////////////////////////////////
// (B) Regularization type
//////////////////////////////////////////////////////////////////////

#define SINGLE_HYPERPARAMETER                      0
#define MULTIPLE_HYPERPARAMETERS                   1
#define ARD_HYPERPARAMETERS                        0

//////////////////////////////////////////////////////////////////////
// (C) Options related to regularization hyperparameter estimation
//////////////////////////////////////////////////////////////////////

// Three possible modes:
//    -- holdout cross-validation via grid search
//    -- holdout cross-validation via gradient-based optimization
//    -- majorization-minimization

#define HYPERPARAMETER_GRID_SEARCH                 0
#define HYPERPARAMETER_GRADIENT_OPTIMIZATION       1
#define HYPERPARAMETER_MAJORIZATION_MINIMIZATION   0

//////////////////////////////////////////////////////////////////////
// (C1) Grid-search options
//////////////////////////////////////////////////////////////////////

// use logloss instead of regular holdout loss for holdout cross-validation
#define CROSS_VALIDATE_USING_LOGLOSS               1

//////////////////////////////////////////////////////////////////////
// (C2) Gradient-based optimization options
//////////////////////////////////////////////////////////////////////

// starting regularization parameter
const double INITIAL_LOG_C = 5.0;

//////////////////////////////////////////////////////////////////////
// (C3) majorization-minimization-options
//////////////////////////////////////////////////////////////////////

// number of iterative relinearization steps if using
// majorization-minmimization algorithm
const int NUM_ITERATIVE_RELINEARIZATION_STEPS = 5;

// smoothing used for majorization-minimization algorithm
const double MM_SMOOTHING = 1.0;

//////////////////////////////////////////////////////////////////////
// (D) Input type
//////////////////////////////////////////////////////////////////////

#ifndef RNA
#error Compilation failed!  RNA variable should be defined.
#endif

//////////////////////////////////////////////////////////////////////
// (E) Used parameter groups
//////////////////////////////////////////////////////////////////////

#if RNA

#define PARAMS_MATCH                               1
#define PARAMS_INSERT                              1
#define PARAMS_SINGLE                              1
#define PARAMS_PAIR                                1
#define PARAMS_HYDROPATHY                          0
#define PARAMS_COMPRESSED                          0
#if FORCE_UNIQUE_PARSES
#define PARAMS_DOUBLE_AFFINE                       0
#else
#define PARAMS_DOUBLE_AFFINE                       1
#endif
#define PARAMS_TERMINAL_INSERTS                    0

const std::string alphabet = "ACGU";                       // allowed symbols -- all other letters ignored
const int M = 4;                                           // number of alphabet symbols

#else

#define PARAMS_MATCH                               1       // 23x23 substitution matrix
#define PARAMS_INSERT                              1       // 23x1 emission matrix
#define PARAMS_SINGLE                              1       
#define PARAMS_PAIR                                1       
#define PARAMS_HYDROPATHY                          1
#define PARAMS_COMPRESSED                          1
#if FORCE_UNIQUE_PARSES
#define PARAMS_DOUBLE_AFFINE                       0
#else
#define PARAMS_DOUBLE_AFFINE                       1
#endif
#define PARAMS_TERMINAL_INSERTS                    0

const std::string alphabet = "ARNDCQEGHILKMFPSTWYVBZX";    // allowed symbols -- all other letters ignored
const int M = 23;                                          // number of alphabet symbols

const std::string hydrophilic_alphabet = "DEGKNQPRS";      // hydrophilic symbols
const int H = 6;                                           // hydropathy window size

const int COMPRESSED_M = 6;
const std::string compressed_alphabet[COMPRESSED_M] = {"AGPST", "C", "DENQ", "FWY", "HKR", "ILMV"};

#endif

/////////////////////////////////////////////////////////////////////
// (G) BMRM stuff
//////////////////////////////////////////////////////////////////////

#define BMRM_AVAILABLE                              0
// #define DAIFLETCHER
}
#endif
