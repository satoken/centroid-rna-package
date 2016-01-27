#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(centroid_fold, 0.0.12, no-define)
AM_INIT_AUTOMAKE([-Wall foreign])
AM_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR([src/centroid_fold.cpp])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_LANG(C++)
AC_CXX_HAVE_STL
AC_PROG_RANLIB
AC_PROG_LN_S
AC_PROG_AWK
AC_PROG_SED
AM_PROG_AR

dnl # Memory Resevation
dnl AC_ARG_WITH([reserve], AS_HELP_STRING([--with-reserve],[the length of seqeunces for reserved memory (default: 500)]),
dnl     [
dnl       AC_DEFINE_UNQUOTED([RESERVED_LENGTH], $withval, [the length of seqeunces for reserved memory])
dnl     ],
dnl     [
dnl       AC_DEFINE_UNQUOTED([RESERVED_LENGTH], 500, [the length of seqeunces for reserved memory])
dnl     ]
dnl )

# debug option
AC_ARG_ENABLE([debug], AS_HELP_STRING([--enable-debug], [turn on debugging (default: off)]))
if test "$enable_debug" == "yes"; then
  CFLAGS="-DDEBUG $CFLAGS"
  CXXFLAGS="-DDEBUG $CXXFLAGS"
else
  CFLAGS="-DNDEBUG $CFLAGS"
  CXXFLAGS="-DNDEBUG $CXXFLAGS"
fi

# random number generator
AC_ARG_ENABLE([mersenne-twister], AS_HELP_STRING([--enable-mersenne-twister], [use Mersenne Twister as a random number generator (default: on)]))
if test "$enable_mersenne_twister" == "no"; then
   AC_DEFINE([STD_RAND], [1], [use the standard random number generator instead of Mersenne Twister])
fi

# boost libraries
AX_BOOST_BASE(1.34, , AC_MSG_ERROR("the boost library is unavailable"))
if test "x$BOOST_CPPFLAGS" != x; then
  CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
fi
if test "x$BOOST_LDFLAGS" != x; then
  LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
fi
dnl AC_CHECK_HEADER(boost/lambda/lambda.hpp, ,
dnl                 AC_MSG_ERROR("boost/lambda/lambda.hpp is unavailable"))
AC_CHECK_HEADER(boost/multi_array.hpp, ,
                AC_MSG_ERROR("boost/multi_array.hpp is unavailable"))
AC_CHECK_HEADER(boost/random.hpp, ,
                AC_MSG_ERROR("boost/random.hpp is unavailable"))

AX_BOOST_PROGRAM_OPTIONS
if test "x$BOOST_PROGRAM_OPTIONS_LIB" != x; then
  LIBS="$LIBS $BOOST_PROGRAM_OPTIONS_LIB"
else
  AC_MSG_ERROR("the Boost::Program_Options library is unavailable")  
fi

AX_BOOST_SYSTEM
if test "x$BOOST_SYSTEM_LIB" != x; then
  LIBS="$LIBS $BOOST_SYSTEM_LIB"
else
  AC_MSG_ERROR("the Boost::Program_Options library is unavailable")  
fi



# Checks for libraries.
AC_ARG_WITH([vienna-rna], AS_HELP_STRING([--with-vienna-rna], [support Vienna RNA package]),
  if test "$with_vienna_rna" != "yes" && test "$with_vienna_rna" != "no" ; then
     CPPFLAGS="-I${with_vienna_rna}/include $CPPFLAGS"
     LDFLAGS="-L${with_vienna_rna}/lib $LDFLAGS"
  fi
)
AC_CHECK_LIB([RNA], [pf_fold], [],
[
  AC_MSG_NOTICE([try '-fopenmp'])
  CXXFLAGS_BAK=${CXXFLAGS}
  CXXFLAGS="-fopenmp ${CXXFLAGS}"
  unset ac_cv_lib_RNA_pf_fold
  AC_CHECK_LIB([RNA], [pf_fold], [], [ CXXFLAGS=${CXXFLAGS_BAK} ])
])

AC_CHECK_FUNC(energy_of_structure,
	      [AC_DEFINE([HAVE_VIENNA20], [], [have the Vienna RNA package version 1.8 or above])])

# Checks for header files.
AC_CHECK_HEADER([ViennaRNA/ribo.h],
		[AC_DEFINE([HAVE_VIENNA18], [], [have the Vienna RNA package version 2.0 or above])])

AC_CHECK_FUNC([mkstemp],
	      [AC_DEFINE([HAVE_MKSTEMP], [], [have mkstemp() in the standard library])])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE

# Checks for library functions.
AC_HEADER_STDC

AC_CONFIG_FILES([Makefile src/Makefile src/contrafold/Makefile src/probconsRNA/Makefile src/contralign/Makefile])
AC_OUTPUT
