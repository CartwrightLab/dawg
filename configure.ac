m4_define([dawg_rev], [m4_esyscmd([svnversion -n .])])
AC_INIT([dawg], [1.1.0.dawg_rev()], [http://scit.us/dawg/])
AC_CONFIG_SRCDIR([src/dawg.h])
AC_CONFIG_AUX_DIR([config])
AC_PREREQ(2.59)
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_YACC
AM_PROG_LEX
AC_CHECK_PROG(DIFF, diff, diff)

AC_LANG_PUSH(C++)

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([time.h float.h stdio.h unistd.h math.h assert.h stdarg.h stddef.h stdtypes.h process.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_C_BIGENDIAN
AC_CXX_MEMBER_TEMPLATES_OUTSIDE_CLASS
AC_CXX_MEMBER_TEMPLATES

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([getpid _getpid copysign _copysign])
dnl AC_CHECK_FUNCS([floor pow sqrt])

AC_LANG_POP(C++)

AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile tests/Makefile doc/Makefile])
AC_OUTPUT
