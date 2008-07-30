bin_PROGRAMS = centroid_fold
noinst_LIBRARIES = libcentroid.a
dist_bin_SCRIPTS = contrafold.rb simfold.rb 

centroid_fold_SOURCES = \
	main.cpp fa.cpp aln.cpp iupac.cpp fa.h aln.h 
centroid_fold_LDADD = libcentroid.a ${LIBCONTRAFOLD} ${BOOST_LIBS}

libcentroid_a_SOURCES = \
	centroid_fold.cpp centroid.cpp mea.cpp bp.cpp iupac.cpp \
	centroid_fold.h centroid.h mea.h bp.h iupac.h cyktable.h rule.h contrafold.h
