AC_INIT(dawg, 1.0.0-rc1, rac@uga.edu)
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

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([time.h float.h stdio.h unistd.h math.h assert.h stdarg.h stddef.h process.h])

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
AC_CHECK_FUNCS([getpid _getpid])

AC_ARG_ENABLE(mmx, AC_HELP_STRING([--enable-mmx], [use mmx if available (default=yes)]), [enable_mmx=$enableval], [enable_mmx=yes])

if test x$enable_mmx = xyes; then
AX_CFLAGS_GCC_OPTION([-mmmx])
AX_CXXFLAGS_GCC_OPTION([-mmmx])

AC_MSG_CHECKING(for mmintrin.h)
AC_LANG_PUSH(C++)
AC_COMPILE_IFELSE(AC_LANG_PROGRAM([#include <mmintrin.h>],[_mm_empty()]),
[AC_MSG_RESULT(yes)
 AC_DEFINE(USE_MMX, 1, [Define to 1 to use MMX optimizations.])
 AC_DEFINE(HAVE_MMINTRIN_H, 1, [Define to 1 if you have the <mmintrin.h> header file.])],
 [AC_MSG_RESULT(no)
  enable_mmx=no])
AC_LANG_POP(C++)
fi

AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile tests/Makefile doc/Makefile])
AC_OUTPUT
