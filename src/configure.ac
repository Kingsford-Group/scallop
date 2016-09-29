#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(scallop, 1.0, mingfu.shao@cs.cmu.edu)
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/main.cc])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

### require environment variables BOOST_HOME GUROBI_HOME GUROBI_LICENCE_FILE HTSLIB
AS_IF([test "s$HTSLIB" = "s"],[AC_MSG_ERROR([*** Please set environment variable HTSLIB ***])]) 
AS_IF([test "s$BOOST_HOME" = "s"],[AC_MSG_ERROR([*** Please set environment variable BOOST_HOME ***])]) 
AS_IF([test "s$GUROBI_HOME" = "s"],[AC_MSG_ERROR([*** Please set environment variable GUROBI_HOME ***])]) 
AS_IF([test "s$GRB_LICENSE_FILE" = "s"],[AC_MSG_ERROR([*** Please set environment variable GRB_LICENSE_FILE ***])]) 

#AC_ARG_WITH(samtools, [--with-samtools prefix for samtools], SAMTOOLS=$withval, SAMTOOLS=)
#AC_SUBST(SAMTOOLS)
#AS_IF([test "x$SAMTOOLS" = "x"],[AC_MSG_ERROR([*** Please use the --with-samtools to setup the prefix path for samtools ***])]) 


# Checks for libraries.
# FIXME: Replace `main' with a function in `-lbam':
#AC_CHECK_LIB([bam], [main])
# FIXME: Replace `main' with a function in `-lhts':
#AC_CHECK_LIB([hts], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdint.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT32_T

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
