SUBDIRS = contrafold
bin_PROGRAMS = centroid_fold
noinst_LIBRARIES = libcentroid.a
dist_bin_SCRIPTS = contrafold.rb simfold.rb
CLEANFILES = centroid_alifold$(EXEEXT)
all: centroid_alifold$(EXEEXT)

centroid_fold_SOURCES = \
	main.cpp fa.cpp aln.cpp fa.h aln.h
centroid_fold_LDADD = libcentroid.a contrafold/libcontrafold.a ${BOOST_LIBS}

libcentroid_a_SOURCES = \
	centroid_fold.cpp centroid.cpp mea.cpp bp.cpp \
	engine/aux.cpp engine/alifold.cpp engine/contrafold.cpp engine/mccaskill.cpp \
	engine/averaged.cpp engine/contrafoldm.cpp engine/mixture.cpp \
	centroid_fold.h centroid.h mea.h bp.h cyktable.h rule.h contrafold/contrafold.h \
	engine/aux.h engine/alifold.h engine/contrafold.h engine/mccaskill.h \
	engine/averaged.h engine/contrafoldm.h engine/mixture.h \
        diana.cpp diana.h ps_plot.cpp ps_plot.h naview.c 

install-exec-hook:
	$(LN_S) $(DESTDIR)$(bindir)/centroid_fold$(EXEEXT) \
		$(DESTDIR)$(bindir)/centroid_alifold$(EXEEXT)

centroid_alifold$(EXEEXT): centroid_fold$(EXEEXT)
	$(LN_S) centroid_fold$(EXEEXT) centroid_alifold$(EXEEXT)
