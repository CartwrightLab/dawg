AC_PREREQ(2.53)
AC_INIT(src/dawg.cpp)
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([src/dawg.cc])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE(dawg, 1.0.0)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AM_PROG_LEX
AC_PROG_YACC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([time.h float.h stdio.h unistd.h math.h assert.h stdarg.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_C_BIGENDIAN

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([sqrt])

AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile tests/Makefile doc/Makefile])
AC_OUTPUT
