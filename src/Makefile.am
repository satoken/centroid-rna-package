SUBDIRS = contrafold
bin_PROGRAMS = centroid_fold
noinst_LIBRARIES = libcentroid.a
dist_bin_SCRIPTS = contrafold.rb simfold.rb
CLEANFILES = centroid_alifold$(EXEEXT)
all: centroid_alifold$(EXEEXT)

centroid_fold_SOURCES = \
	main.cpp fa.cpp aln.cpp iupac.cpp fa.h aln.h
centroid_fold_LDADD = libcentroid.a contrafold/libcontrafold.a ${BOOST_LIBS}

libcentroid_a_SOURCES = \
	centroid_fold.cpp centroid.cpp mea.cpp bp.cpp iupac.cpp \
	centroid_fold.h centroid.h mea.h bp.h iupac.h cyktable.h rule.h contrafold/contrafold.h \
        diana.cpp diana.h

install-exec-hook:
	$(LN_S) $(DESTDIR)$(bindir)/centroid_fold$(EXEEXT) \
		$(DESTDIR)$(bindir)/centroid_alifold$(EXEEXT)

centroid_alifold$(EXEEXT): centroid_fold$(EXEEXT)
	$(LN_S) centroid_fold$(EXEEXT) centroid_alifold$(EXEEXT)
