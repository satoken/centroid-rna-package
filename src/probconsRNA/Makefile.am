OTHERFLAGS = -DNumInsertStates=1 -DVERSION="1.10"
AM_CXXFLAGS = -O3 -W -Wall -pedantic -DNDEBUG $(OTHERFLAGS) -funroll-loops 

noinst_LIBRARIES = libprobcons.a

libprobcons_a_SOURCES =  wrapper.cpp \
	Defaults.h ScoreType.h ProbabilisticModel.h probcons.h \
	SafeVector.h SparseMatrix.h MultiSequence.h Sequence.h FileBuffer.h 