INCLUDES =  -I$(top_builddir) -I$(top_srcdir)
AM_YFLAGS = -d
AM_LFLAGS = -CF
BUILT_SOURCES = parser.h

bin_PROGRAMS = dawg
dawg_SOURCES = dawg.cc eigen.cc  indel.cc matrix.cc output.cc rand.cc tree.cc var.cc dawg.h indel.h tree.h rand.h var.h matrix.h dawgtxt.h parser.yy lexer.ll
dawg_LDADD = -lstdc++ @LEXLIB@
CLEANFILES = parser.cc parser.h lexer.cc
