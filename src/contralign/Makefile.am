MODEL_TYPE = -DRNA=1

AM_CXXFLAGS = $(MODEL_TYPE) -Wall -Wundef -O3 -DNDEBUG -fomit-frame-pointer \
	-ffast-math -funroll-all-loops -funsafe-math-optimizations \
	-fpeel-loops -Winline --param large-function-growth=100000 \
	--param max-inline-insns-single=100000 \
	--param inline-unit-growth=100000

noinst_LIBRARIES = libcontralign.a

libcontralign_a_SOURCES =  wrapper.cpp MultiSequence.cpp Sequence.cpp Utilities.cpp \
	Defaults.ipp InferenceEngine.hpp ParameterManager.hpp \
	Config.hpp MultiSequence.hpp Utilities.hpp ParameterManager.ipp InferenceEngine.ipp \
	Sequence.hpp Utilities.ipp LogSpace.hpp SparseMatrix.hpp SparseMatrix.ipp
