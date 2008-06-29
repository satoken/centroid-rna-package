SUBDIRS = src
EXTRA_DIST = README.rd README.html contrafold.patch \
	ruby/centroid_fold_wrap.cxx ruby/centroid_fold.i ruby/extconf.rb ruby/test.rb \
	perl/centroid_fold_wrap.cxx perl/centroid_fold.i perl/Makefile.PL perl/test.pl \
	python/centroid_fold_wrap.cxx python/centroid_fold.i python/setup.py python/test.py

README.html: README.rd
	rd2 README.rd > README.html

ruby/centroid_fold_wrap.cxx: ruby/centroid_fold.i
	( cd ruby; swig -c++ -ruby centroid_fold.i )

perl/centroid_fold_wrap.cxx: perl/centroid_fold.i
	( cd perl; swig -c++ -perl centroid_fold.i )

python/centroid_fold_wrap.cxx: python/centroid_fold.i
	( cd python; swig -c++ -python centroid_fold.i )
