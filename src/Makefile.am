bin_PROGRAMS = centroid_fold
dist_bin_SCRIPTS = contrafold.rb contrafold_m.rb simfold.rb simfold_m.rb

centroid_fold_SOURCES = \
	centroid_fold.cpp centroid.cpp mea.cpp bp.cpp fa.cpp aln.cpp iupac.cpp \
	centroid.h mea.h bp.h fa.h aln.h iupac.h cyktable.h rule.h contrafold.h
centroid_fold_LDADD = ${LIBCONTRAFOLD} ${BOOST_LIBS}
