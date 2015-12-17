AUTOMAKE_OPTIONS = subdir-objects

SUBDIRS = contrafold probconsRNA contralign
bin_PROGRAMS = centroid_fold centroid_alifold centroid_homfold
noinst_LIBRARIES = libcentroid.a
dist_bin_SCRIPTS = contrafold.rb simfold.rb

centroid_fold_SOURCES = centroid_fold.cpp 
centroid_fold_LDADD = libcentroid.a contrafold/libcontrafold.a

centroid_alifold_SOURCES = centroid_alifold.cpp
centroid_alifold_LDADD = libcentroid.a contrafold/libcontrafold.a

centroid_homfold_SOURCES = centroid_homfold.cpp
centroid_homfold_LDADD = libcentroid.a contrafold/libcontrafold.a contralign/libcontralign.a probconsRNA/libprobcons.a

libcentroid_a_SOURCES = \
	folding_engine.cpp folding_engine.h centroid.cpp centroid.h \
	mea.cpp mea.h bp.cpp bp.h cyktable.h rule.h util.h\
	fa.cpp fa.h aln.cpp aln.h \
	engine/aux.cpp engine/aux.h \
	engine/alifold.cpp engine/alifold.h \
	engine/contrafold.cpp engine/contrafold.h \
	engine/mccaskill.cpp engine/mccaskill.h \
	engine/averaged.cpp engine/averaged.h \
	engine/contrafoldm.cpp engine/contrafoldm.h \
	engine/mixture.cpp engine/mixture.h \
	engine/pfold.cpp engine/pfold.h \
	engine/contrafoldhom.cpp engine/contrafoldhom.h \
	engine/mccaskillhom.cpp engine/mccaskillhom.h \
	contrafold/contrafold.h \
        probconsRNA/probcons.h \
	contralign/contralign.h \
        diana.cpp diana.h ps_plot.cpp ps_plot.h naview.c \
	engine/boltzmann_param.h engine/boltzmann_param.c
