/////////////////////////////////////////////////////////////////
// Contralign.cpp
/////////////////////////////////////////////////////////////////

// include files
#ifdef MULTI
#include <mpi.h>
#endif
#include "Config.hpp"
#include "Options.hpp"
#include "Utilities.hpp"
#include "ComputationWrapper.hpp"
#include "FileDescription.hpp"
#include "InferenceEngine.hpp"
#include "ParameterManager.hpp"
#include "OptimizationWrapper.hpp"

namespace CONTRALIGN{
// constants
const double GAMMA_DEFAULT = 100000;
const double REGULARIZATION_DEFAULT = 0;
const int PC_ITERS_DEFAULT = 2;
const int SC_ITERS_DEFAULT = 0;

// function prototypes
void Usage(const Options &options);
void Version();
void ParseArguments(int argc, char **argv, Options &options, std::vector<std::string> &filenames);
void MakeFileDescriptions(const Options &options, const std::vector<std::string> &filenames, std::vector<FileDescription> &descriptions);

template<class RealT>
void RunGradientSanityCheck(const Options &options, const std::vector<FileDescription> &descriptions);

template<class RealT>
void RunTrainingMode(const Options &options, const std::vector<FileDescription> &descriptions);

template<class RealT>
void RunPredictionMode(const Options &options, const std::vector<FileDescription> &descriptions);
}
// default parameters
#include "Defaults.ipp"

/////////////////////////////////////////////////////////////////
// main()
//
// Main program.
/////////////////////////////////////////////////////////////////
namespace CONTRALIGN{
int _main(int argc, char **argv)
{
#ifdef MULTI
    MPI_Init(&argc, &argv);
#endif

    // first, parse arguments
    Options options;
    std::vector<std::string> filenames;
    ParseArguments(argc, argv, options, filenames);

    // second, read input files
    std::vector<FileDescription> descriptions;
    MakeFileDescriptions(options, filenames, descriptions);
    
    // perform required task
    if (options.GetBoolValue("gradient_sanity_check"))
    {
        RunGradientSanityCheck<double>(options, descriptions);
    }
    else if (options.GetBoolValue("training_mode"))
    {
        RunTrainingMode<double>(options, descriptions);
    }
    else
    {
        RunPredictionMode<float>(options, descriptions);
    }
    
#ifdef MULTI
    MPI_Finalize();
#endif
    return 1;
}

/////////////////////////////////////////////////////////////////
// Usage()
//
// Display program usage.
/////////////////////////////////////////////////////////////////

void Usage(const Options &options)
{
    std::cerr << std::endl
              << "Usage: contralign [predict|train] [OPTION]... INFILE(s)" << std::endl 
              << std::endl
              << "       where [OPTION]...   is a list of zero or more optional arguments" << std::endl
              << "             INFILE(s)     is the name of the input FASTA file(s)" << std::endl 
              << std::endl
              << "Miscellaneous arguments:" << std::endl
              << "  --version                display program version information" << std::endl
              << "  --verbose                show detailed console output" << std::endl
              << "  --logbase LOG_BASE       set base of log-sum-exp" << std::endl
              << "  --viterbi                use Viterbi instead of posterior decoding for prediction, " << std::endl
              << "                           or max-margin instead of log-likelihood for training" << std::endl
              << std::endl 
              << "Additional arguments for 'predict' mode:" << std::endl
              << "  --params FILENAME        use particular model parameters" << std::endl
              << "  --constraints            use existing constraints (requires input to contain pre-aligned sequences)" << std::endl
              << "  --gamma GAMMA            set sensivity/specificity tradeoff parameter (default: GAMMA=" << options.GetRealValue("gamma") << ")" << std::endl
              << "                             if GAMMA > 1, emphasize sensitivity" << std::endl
              << "                             if 0 <= GAMMA <= 1, emphasize specificity" << std::endl
              << "                             if GAMMA < 0, try tradeoff parameters of 2^k for k = -5,...,10" << std::endl
              << "  --conflate               process all input filenames as sequences in a single input" << std::endl
              << "  --pairwise               produce all-pairs pairwise alignments" << std::endl
              << "  --annealing              use sequence annealing algorithm" << std::endl
              << "  --pc_iters NUM           perform NUM probabilistic consistency iterations (default: NUM=" << options.GetIntValue("probabilistic_consistency_iterations") << ")" << std::endl
              << "  --sc_iters NUM           perform NUM spectral consistency iterations (default: NUM=" << options.GetIntValue("spectral_consistency_iterations") << ")" << std::endl
              << std::endl
              << "  --mfa OUTFILEORDIR       write MFA output to file or directory" << std::endl
              << "  --clustalw OUTFILEORDIR  write CLUSTALW output to file or directory" << std::endl
              << "  --posteriors CUTOFF OUTFILEORDIR" << std::endl
              << "                           write posterior match probabilities to file or directory" << std::endl
              << "  --partition              compute the partition function or Viterbi score only" << std::endl
              << std::endl
              << "Additional arguments for training (many input files may be specified):" << std::endl
              << "  --sanity                 perform gradient sanity check" << std::endl
              << "  --holdout F              use fraction F of training data for holdout cross-validation" << std::endl
              << "  --regularize C           perform BFGS training, using a single regularization coefficient C" << std::endl
              << std::endl;
    exit(0);
}

/////////////////////////////////////////////////////////////////
// Version()
//
// Display program version.
/////////////////////////////////////////////////////////////////

void Version()
{
#if RNA
    std::cerr << "CONTRAlign version 2.00 - RNA multiple sequence alignment" << std::endl << std::endl
#else
    std::cerr << "CONTRAlign version 2.00 - Protein multiple sequence alignment" << std::endl << std::endl
#endif
              << "Written by Chuong B. Do" << std::endl;
    exit(0);
}

/////////////////////////////////////////////////////////////////
// ParseArguments()
//
// Parse command line parameters.
/////////////////////////////////////////////////////////////////

void ParseArguments(int argc,
                    char **argv,
                    Options &options,
                    std::vector<std::string> &filenames)
{
    // register default options
    options.SetBoolValue("training_mode", false);

    options.SetBoolValue("verbose_output", false);
    options.SetRealValue("log_base", 1.0);
    options.SetBoolValue("viterbi_parsing", false);

    options.SetStringValue("parameter_filename", "");
    options.SetBoolValue("use_constraints", false);
    options.SetRealValue("gamma", GAMMA_DEFAULT);
    options.SetBoolValue("conflate", false);
    options.SetBoolValue("pairwise", false);
    options.SetBoolValue("annealing", false);
    options.SetIntValue("probabilistic_consistency_iterations", PC_ITERS_DEFAULT);
    options.SetIntValue("spectral_consistency_iterations", SC_ITERS_DEFAULT);
    options.SetStringValue("output_mfa_destination", "");
    options.SetStringValue("output_clustalw_destination", "");
    options.SetRealValue("output_posteriors_cutoff", 0);
    options.SetStringValue("output_posteriors_destination", "");
    options.SetBoolValue("partition_function_only", false);

    options.SetBoolValue("gradient_sanity_check", false);
    options.SetRealValue("holdout_ratio", 0);
    options.SetRealValue("regularization_coefficient", REGULARIZATION_DEFAULT);

    // check for sufficient arguments
    if (argc < 2) Usage(options);
    filenames.clear();

    // check for prediction or training mode    
    if (!strcmp(argv[1], "train"))
        options.SetBoolValue("training_mode", true);
    else
        if (strcmp(argv[1], "predict"))
            Error("CONTRAlign must be run in either 'predict' or 'train' mode.");

    // go through remaining arguments
    for (int argno = 2; argno < argc; argno++)
    {
        // parse optional arguments
        if (argv[argno][0] == '-')
        {
            // miscellaneous options
            if (!strcmp(argv[argno], "--version"))
            {
                Version();
            }
            else if (!strcmp(argv[argno], "--verbose"))
            {
                options.SetBoolValue("verbose_output", true);
            }
            else if (!strcmp(argv[argno], "--logbase"))
            {
                if (argno == argc - 1) Error("Must specify log base LOG_BASE after --logbase.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse log base.");
                if (value <= 0)
                    Error("Log base must be positive.");
                options.SetRealValue("log_base", value);
            }
            else if (!strcmp(argv[argno], "--viterbi"))
            {
                options.SetBoolValue("viterbi_parsing", true);
            }
            
            // prediction options
            else if (!strcmp(argv[argno], "--params"))
            {
                if (argno == argc - 1) Error("Must specify FILENAME after --params.");
                options.SetStringValue("parameter_filename", argv[++argno]);
            }
            else if (!strcmp(argv[argno], "--constraints"))
            {
                options.SetBoolValue("use_constraints", true);
            }
            else if (!strcmp(argv[argno], "--gamma"))
            {
                if (argno == argc - 1) Error("Must specify trade-off parameter GAMMA after --gamma.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse value after --gamma.");
                options.SetRealValue("gamma", value);
            }
            else if (!strcmp(argv[argno], "--conflate"))
            {
                options.SetBoolValue("conflate", true);
            }
            else if (!strcmp(argv[argno], "--pairwise"))
            {
                options.SetBoolValue("pairwise", true);
            }
            else if (!strcmp(argv[argno], "--annealing"))
            {
                options.SetBoolValue("annealing", true);
            }
            else if (!strcmp(argv[argno], "--pc_iters"))
            {
                if (argno == argc - 1) Error("Must specify parameter NUM after --pc_iters.");
                int value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse value after --pc_iters.");
                options.SetIntValue("probabilistic_consistency_iterations", value);
            }
            else if (!strcmp(argv[argno], "--sc_iters"))
            {
                if (argno == argc - 1) Error("Must specify parameter NUM after --sc_iters.");
                int value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse value after --sc_iters.");
                options.SetIntValue("spectral_consistency_iterations", value);
            }
            else if (!strcmp(argv[argno], "--mfa"))
            {
                if (argno == argc - 1) Error("Must specify output file or directory name after --mfa.");
                options.SetStringValue("output_mfa_destination", argv[++argno]);
            }
            else if (!strcmp(argv[argno], "--clustalw"))
            {
                if (argno == argc - 1) Error("Must specify output file or directory name after --clustalw.");
                options.SetStringValue("output_clustalw_destination", argv[++argno]);
            }
            else if (!strcmp(argv[argno], "--posteriors"))
            {
                if (argno == argc - 1) Error("Must specify posterior probability threshold CUTOFF after --posteriors.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse cutoff value after --posteriors.");
                options.SetRealValue("output_posteriors_cutoff", value);
                if (argno == argc - 1) Error("Must specify output file or directory for --posteriors.");
                options.SetStringValue("output_posteriors_destination", argv[++argno]);
            }
            else if (!strcmp(argv[argno], "--partition"))
            {
                options.SetBoolValue("partition_function_only", true);
            }
            
            // training options
            else if (!strcmp(argv[argno], "--sanity"))
            {
                options.SetBoolValue("gradient_sanity_check", true);
            }
            else if (!strcmp(argv[argno], "--holdout"))
            {
                if (argno == argc - 1) Error("Must specify holdout ratio F after --holdout.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse holdout ratio.");
                if (value < 0 || value > 1)
                    Error("Holdout ratio must be between 0 and 1.");
                options.SetRealValue("holdout_ratio", value);
            }
            else if (!strcmp(argv[argno], "--regularize"))
            {
                if (argno == argc - 1) Error("Must specify regularization parameter C after --regularize.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse regularization parameter after --regularize.");
                if (value < 0)
                    Error("Regularization parameter should not be negative.");
                options.SetRealValue("regularization_coefficient", value);
            }
            else
            {
                Error("Unknown option \"%s\" specified.  Run program without any arguments to see command-line options.", argv[argno]);
            }
            
        }
        else
        {
            filenames.push_back(argv[argno]);
        }
    }

    // ensure that at least one input file specified
    if (filenames.size() == 0)
        Error("No filenames specified.");

    // check to make sure that arguments make sense
    if (options.GetBoolValue("training_mode"))
    {
        if (options.GetStringValue("parameter_filename") != "")
            Error("Should not specify parameter file for training mode.");
        if (options.GetBoolValue("use_constraints"))
            Error("The --constraints flag has no effect in training mode.");
        if (options.GetRealValue("gamma") != GAMMA_DEFAULT)
            Error("Gamma parameter should not be specified in training mode.");
        if (options.GetBoolValue("conflate"))
            Error("The --conflate flags should not be used in training mode.");
        if (options.GetBoolValue("pairwise"))
            Error("The --pairwise flag has no effect in training mode.");
        if (options.GetBoolValue("annealing"))
            Error("The --annealing flag has no effect in training mode.");
        if (options.GetIntValue("probabilistic_consistency_iterations") != PC_ITERS_DEFAULT)
            Error("The --pc_iters flag has no effect in training mode.");
        if (options.GetIntValue("spectral_consistency_iterations") != SC_ITERS_DEFAULT)
            Error("The --sc_iters flag has no effect in training mode.");
        if (options.GetStringValue("output_mfa_destination") != "")
            Error("The --mfa option cannot be used in training mode.");
        if (options.GetStringValue("output_clustalw_destination") != "")
            Error("The --clustalw option cannot be used in training mode.");
        if (options.GetStringValue("output_posteriors_destination") != "" ||
            options.GetRealValue("output_posteriors_cutoff") != 0)
            Error("The --posteriors option cannot be used in training mode.");
        if (options.GetBoolValue("partition_function_only"))
            Error("The --partition flag cannot be used in training mode.");
        if (options.GetRealValue("regularization_coefficient") != REGULARIZATION_DEFAULT &&
            options.GetRealValue("holdout_ratio") > 0)
            Error("The --holdout and --regularize options cannot be specified simultaneously.");
    }
    else
    {
        if (options.GetRealValue("gamma") < 0 &&
            options.GetStringValue("output_mfa_destination") == "" &&
            options.GetStringValue("output_clustalw_destination") == "" &&
            options.GetStringValue("output_posteriors_destination") == "")
            Error("Output directory must be specified when using GAMMA < 0.");

#ifdef MULTI
        if (filenames.size() > 1 &&
            options.GetStringValue("output_mfa_destination") == "" &&
            options.GetStringValue("output_clustalw_destination") == "" &&
            options.GetStringValue("output_posteriors_destination") == "")
            Error("Output directory must be specified when performing predictions for multiple input files.");
#endif
        if (options.GetBoolValue("viterbi_parsing") &&
            options.GetStringValue("output_posteriors_destination") != "")
            Error("The --posteriors option cannot be used with Viterbi parsing.");
    }
}

/////////////////////////////////////////////////////////////////
// MakeFileDescriptions()
//
// Build file descriptions
/////////////////////////////////////////////////////////////////

void MakeFileDescriptions(const Options &options,
                          const std::vector<std::string> &filenames,
                          std::vector<FileDescription> &descriptions)
{
    descriptions.clear();
    for (size_t i = 0; i < filenames.size(); i++)
    {
        descriptions.push_back(FileDescription(filenames[i]));
    }

    if (options.GetBoolValue("conflate"))
    {
        for (size_t i = 1; i < descriptions.size(); i++)
        {
            for (int j = 0; j < descriptions[i].seqs.GetNumSequences(); j++)
            {
                descriptions[0].seqs.AddSequence(new Sequence(descriptions[i].seqs.GetSequence(j).GetData(),
                                                              descriptions[i].seqs.GetSequence(j).GetName(),
                                                              descriptions[0].seqs.GetNumSequences()));
            }
        }
        descriptions.resize(1);
    }
    
    std::sort(descriptions.begin(), descriptions.end());
}

/////////////////////////////////////////////////////////////////
// RunGradientSanityCheck()
//
// Compute gradient sanity check.
/////////////////////////////////////////////////////////////////

template<class RealT>
void RunGradientSanityCheck(const Options &options,
                            const std::vector<FileDescription> &descriptions)
{
    // The architecture of the code is somewhat complicated here, so
    // here's a quick explanation:
    // 
    //    ParameterManager: associates each parameter of the model
    //                      with a name and manages hyperparameter
    //                      groups
    //                     
    //    InferenceEngine: performs application-specific
    //                     (loss-augmented) inference
    //
    //    ComputationEngine: makes all necessary calls to dynamic
    //                       programming routines for processing
    //                       individual sequences and interfaces with
    //                       distributed computation module
    //
    //    ComputationWrapper: provides a high-level interface for
    //                        performing computations on groups of
    //                        sequences
    //
    //    OuterOptimizationWrapper / InnerOptimizationWrapper:
    //                        interface between computation routines
    //                        and optimization routines
    
    ParameterManager<RealT> parameter_manager;
    InferenceEngine<RealT> inference_engine;
    inference_engine.RegisterParameters(parameter_manager);
    ComputationEngine<RealT> computation_engine(options, descriptions, inference_engine, parameter_manager);
    ComputationWrapper<RealT> computation_wrapper(computation_engine);

    // decide whether I'm a compute node or master node
    if (computation_engine.IsComputeNode())
    {
        computation_engine.RunAsComputeNode();
        return;
    }

    std::vector<RealT> w(parameter_manager.GetNumLogicalParameters(), RealT(0));
    computation_wrapper.SanityCheckGradient(computation_wrapper.GetAllUnits(), w);
    computation_engine.StopComputeNodes();
}

/////////////////////////////////////////////////////////////////
// RunTrainingMode()
//
// Run CONTRAlign in training mode.
/////////////////////////////////////////////////////////////////

template<class RealT>
void RunTrainingMode(const Options &options,
                     const std::vector<FileDescription> &descriptions)
{
    ParameterManager<RealT> parameter_manager;
    InferenceEngine<RealT> inference_engine;
    inference_engine.RegisterParameters(parameter_manager);
    ComputationEngine<RealT> computation_engine(options, descriptions, inference_engine, parameter_manager);
    ComputationWrapper<RealT> computation_wrapper(computation_engine);

    // decide whether I'm a compute node or master node
    if (computation_engine.IsComputeNode())
    {
        computation_engine.RunAsComputeNode();
        return;
    }

    std::vector<RealT> w(parameter_manager.GetNumLogicalParameters(), RealT(0));

    for (size_t i = 0; i < w.size(); i++)
    {
        w[i] = RealT(drand48() * 1e-4 - 5e-5);
    }
    
    std::vector<int> units = computation_wrapper.FilterNonparsable(computation_wrapper.GetAllUnits());
    OptimizationWrapper<RealT> optimization_wrapper(computation_wrapper);
    
    // decide between using a fixed regularization parameter or
    // using cross-validation to determine regularization parameters
    if (options.GetRealValue("holdout_ratio") <= 0)
    {
        std::vector<RealT> regularization_coefficients(parameter_manager.GetNumParameterGroups(), options.GetRealValue("regularization_coefficient"));
        optimization_wrapper.Train(units, w, regularization_coefficients);
    }
    else
    {
        optimization_wrapper.LearnHyperparameters(units, w);
    }
    
    parameter_manager.WriteToFile("optimize.params.final", w);
    computation_engine.StopComputeNodes();
}

/////////////////////////////////////////////////////////////////
// RunPredictionMode()
//
// Run CONTRAlign in prediction mode.
/////////////////////////////////////////////////////////////////

template<class RealT>
void RunPredictionMode(const Options &options,
                       const std::vector<FileDescription> &descriptions)
{
    ParameterManager<RealT> parameter_manager;
    InferenceEngine<RealT> inference_engine;
    inference_engine.RegisterParameters(parameter_manager);
    ComputationEngine<RealT> computation_engine(options, descriptions, inference_engine, parameter_manager);
    ComputationWrapper<RealT> computation_wrapper(computation_engine);
    
    // decide whether I'm a compute node or master node
    if (computation_engine.IsComputeNode())
    {
        computation_engine.RunAsComputeNode();
        return;
    }

    const std::string output_mfa_destination = options.GetStringValue("output_mfa_destination");
    const std::string output_clustalw_destination = options.GetStringValue("output_clustalw_destination");
    const std::string output_posteriors_destination = options.GetStringValue("output_posteriors_destination");

    // load parameters
    std::vector<RealT> w;

    if (options.GetStringValue("parameter_filename") != "")
    {
        parameter_manager.ReadFromFile(options.GetStringValue("parameter_filename"), w);
    }
    else
    {
#if RNA
        w = GetDefaultRNAValues<RealT>();
#else
        w = GetDefaultProteinValues<RealT>();
#endif
    }

    if (options.GetRealValue("gamma") < 0)
    {
        // create directories for storing each run
        if (output_mfa_destination != "") MakeDirectory(output_mfa_destination);
        if (output_clustalw_destination != "") MakeDirectory(output_clustalw_destination);
        if (output_posteriors_destination != "") MakeDirectory(output_posteriors_destination);
        
        // try different values of gamma
        for (int k = -5; k <= 10; k++)
        {
            // create output subdirectories, if needed
            const double gamma = Pow(2.0, double(k));

            if (descriptions.size() > 1)
            {
                if (output_mfa_destination != "")
                    MakeDirectory(SPrintF("%s%c%s.gamma=%lf",
                                          output_mfa_destination.c_str(),
                                          DIR_SEPARATOR_CHAR,
                                          GetBaseName(output_mfa_destination).c_str(), gamma));
                if (output_clustalw_destination != "")
                    MakeDirectory(SPrintF("%s%c%s.gamma=%lf",
                                          output_clustalw_destination.c_str(),
                                          DIR_SEPARATOR_CHAR,
                                          GetBaseName(output_clustalw_destination).c_str(), gamma));
                if (output_posteriors_destination != "")
                    MakeDirectory(SPrintF("%s%c%s.gamma=%lf",
                                          output_posteriors_destination.c_str(),
                                          DIR_SEPARATOR_CHAR,
                                          GetBaseName(output_posteriors_destination).c_str(), gamma));
            }
            
            // perform predictions
            computation_wrapper.Predict(computation_wrapper.GetAllUnits(), w, gamma, options.GetRealValue("log_base"));
        }
    }
    else
    {
        // create output directories for output files, if needed
        if (descriptions.size() > 1)
        {
            if (output_mfa_destination != "") MakeDirectory(output_mfa_destination);
            if (output_clustalw_destination != "") MakeDirectory(output_clustalw_destination);
            if (output_posteriors_destination != "") MakeDirectory(output_posteriors_destination);
        }
        
        computation_wrapper.Predict(computation_wrapper.GetAllUnits(), w, options.GetRealValue("gamma"), options.GetRealValue("log_base"));
    }
    computation_engine.StopComputeNodes();
}
}

int main(int argc, char **argv) {
  return CONTRALIGN::_main(argc, argv);
}
