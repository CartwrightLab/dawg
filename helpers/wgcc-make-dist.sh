#!/bin/sh

export CC=wgcc
export CXX=wgcc
export LDFLAGS="$LFLAGS -L/dev/fs/C/Programs/GnuWin32/lib -lgetopt"
DISTDIR=`echo ~/dawg-dist-dir/`

if [ -e ${DISTDIR} ]
then
	rm -rf ${DISTDIR}
fi

gmake distclean

./autogen.sh && \
./configure --prefix=${DISTDIR} && \
gmake check && \
gmake install
mv ${DISTDIR}/bin/dawg ${DISTDIR}/bin/dawg.exe

PACKAGE=`grep "^PACKAGE = " Makefile | sed "s/PACKAGE = //"`
VERSION=`grep "^VERSION = " Makefile | sed "s/VERSION = //"`
DIRNAME=$PACKAGE-$VERSION
ZIPNAME=`echo ${DIRNAME} | tr "A-Z" "a-z"`
rm -f ${ZIPNAME}.zip
mv ${DISTDIR} ${DIRNAME} && zip -rq9 ${ZIPNAME}.zip ${DIRNAME}
rm -rf ${DIRNAME}
