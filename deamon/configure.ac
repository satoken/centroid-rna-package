#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(centroidfold_deamon, 0.0.4, no-define)
AC_CONFIG_SRCDIR([centroidfold_deamon.cpp])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX

# Syslog facility
AC_ARG_WITH([syslog],AS_HELP_STRING([--with-syslog],[Syslog facility (local0-7)]),
  [ AC_DEFINE_UNQUOTED([SYSLOG],"$withval",[Syslog facility]) ],
  [ AC_DEFINE_UNQUOTED([SYSLOG],"none",[Syslog facility]) ]
)

# Timeout period for receiving a request
AC_ARG_WITH([timeout],AS_HELP_STRING([--with-timeout],[Timeout period for receiving a request (second, default: 10)]),
  [ AC_DEFINE_UNQUOTED([TIMEOUT],$withval,[Timeout period for receiving a request]) ],
  [ AC_DEFINE_UNQUOTED([TIMEOUT],10,[Timeout period for receiving a request]) ]
)

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([arpa/inet.h fcntl.h netdb.h netinet/in.h stdlib.h string.h sys/socket.h sys/time.h syslog.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_TIME
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_SELECT_ARGTYPES
AC_CHECK_FUNCS([bzero inet_ntoa select socket strcasecmp strchr strpbrk strstr])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
