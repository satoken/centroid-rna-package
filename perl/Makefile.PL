# File : Makefile.PL
use ExtUtils::MakeMaker;
WriteMakefile(
	'NAME'    => 'CentroidFold',                  # Name of package
        'INC'  => '-I../src -I..',
        'DEFINE' => ' -DHAVE_CONFIG_H',
        'CC'  => 'g++',
        'LD'  => 'g++',
	'LIBS'    => ['-lRNA -lm'],                    # Name of custom libraries
	'OBJECT'  => 'centroid_fold_wrap.o ../src/libcentroid.a ../src/libcontrafold.a '  # Object files
);
