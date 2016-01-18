#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(scallop, 1.0, mingfu.shao@cs.cmu.edu)
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/boundary.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

### require environment variables BOOST_HOME GUROBI_HOME GUROBI_LICENCE_FILE HTSLIB

## only htslib is enough
#AC_ARG_WITH(samtools, [--with-samtools prefix for samtools], SAMTOOLS=$withval, SAMTOOLS=)
#AC_SUBST(SAMTOOLS)
#AS_IF([test "x$SAMTOOLS" = "x"],[AC_MSG_ERROR([*** Please use the --with-samtools to setup the prefix path for samtools ***])]) 

#AC_ARG_WITH(htslib, [--with-htslib prefix for htslib], HTSLIB=$withval, HTSLIB=)
#AC_SUBST(HTSLIB)
#AS_IF([test "x$HTSLIB" = "x"],[AC_MSG_ERROR([*** Please use the --with-htslib to setup the prefix path for htslib ***])]) 

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
