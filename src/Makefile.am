INCLUDES =  -I$(top_builddir) -I$(top_srcdir)
AM_YFLAGS = -d
AM_LFLAGS = -CF
EXTRA_DIST = $(top_builddir)/examples/*.fud

bin_PROGRAMS = dawg
dawg_SOURCES = dawg.cpp eigen.cpp  indel.cpp matrix.cpp output.cpp rand.cpp tree.cpp var.cpp dawg.h indel.h tree.h rand.h var.h matrix.h parser.yy lexer.ll
dawg_LDADD = -lstdc++ @LEXLIB@

lexer.ll: parser.cc
