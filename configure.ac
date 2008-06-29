#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(centroid_fold, 0.0.4, no-define)
AM_INIT_AUTOMAKE([-Wall foreign])
AM_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR([src/centroid_fold.cpp])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_LANG(C++)
AC_CXX_HAVE_STL
AC_PROG_RANLIB

# CONTRAfold
AC_ARG_WITH([contrafold], AS_HELP_STRING([--with-contrafold],[link the CONTRAfold library]),
    [
      AC_DEFINE([HAVE_LIBCONTRAFOLD], [], [link the CONTRAfold library])
      LIBCONTRAFOLD="libcontrafold.a"
    ]
)
AC_SUBST(LIBCONTRAFOLD)

# Memory Resevation
AC_ARG_WITH([reserve], AS_HELP_STRING([--with-reserve],[the length of seqeunces for reserved memory (default: 500)]),
    [
      AC_DEFINE_UNQUOTED([RESERVED_LENGTH], $withval, [the length of seqeunces for reserved memory])
    ],
    [
      AC_DEFINE_UNQUOTED([RESERVED_LENGTH], 500, [the length of seqeunces for reserved memory])
    ]
)

# boost libraries
dnl AC_CHECK_HEADER(boost/lambda/lambda.hpp, ,
dnl                 AC_MSG_ERROR("boost/lambda/lambda.hpp is unavailable"))
AC_CHECK_HEADER(boost/multi_array.hpp, ,
                AC_MSG_ERROR("boost/multi_array.hpp is unavailable"))

AX_BOOST_PROGRAM_OPTIONS
if test "x$BOOST_PROGRAM_OPTIONS_LIB" != x; then
  BOOST_LIBS="$BOOST_LIBS $BOOST_PROGRAM_OPTIONS_LIB"
else
  AC_MSG_ERROR("the Boost::Program_Options library is unavailable")  
fi
dnl AX_BOOST_THREAD
dnl if test "x$BOOST_THREAD_LIB" != "x"; then
dnl   BOOST_LIBS="$BOOST_LIBS $BOOST_THREAD_LIB"
dnl fi
dnl AX_BOOST_IOSTREAMS
dnl if test "x$BOOST_IOSTREAMS_LIB" != "x"; then
dnl   BOOST_LIBS="$BOOST_LIBS $BOOST_IOSTREAMS_LIB"
dnl fi
AC_SUBST(BOOST_LIBS)

# Checks for libraries.
AC_CHECK_LIB(RNA, init_pf_fold)

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE

# Checks for library functions.
AC_HEADER_STDC

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
