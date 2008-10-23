SUBDIRS = src
EXTRA_DIST = README.rd README.html contrafold_v2_02.patch \
	ruby/README.rd ruby/centroid_fold_wrap.cxx ruby/centroid_fold.i ruby/extconf.rb ruby/test.rb \
	perl/README.rd perl/centroid_fold_wrap.cxx perl/centroid_fold.i perl/Makefile.PL perl/test.pl \
	python/README.rd python/centroid_fold_wrap.cxx python/centroid_fold.i python/setup.py python/test.py

README.html: README.rd
	rd2 README.rd > README.html

ruby/centroid_fold_wrap.cxx: ruby/centroid_fold.i
	( cd ruby; swig -c++ -ruby centroid_fold.i )

perl/centroid_fold_wrap.cxx: perl/centroid_fold.i
	( cd perl; swig -c++ -perl centroid_fold.i )

python/centroid_fold_wrap.cxx: python/centroid_fold.i
	( cd python; swig -c++ -python centroid_fold.i )
