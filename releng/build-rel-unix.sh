#!/bin/sh

PROJ=dawg
PROJ_DISTS=dawg-2*
MAKE=make
CMAKE=cmake
REPOS=`svn info | grep URL: | perl -pe 's!^URL: (.+)/releng$!$1!'`

build_mingw32=
build_m32=
for build_option; do
	case $build_option in
	--mingw32)
		build_mingw32=yes ;;
	--m32)
		build_m32=yes ;;
	esac
done

echo 
echo "Building distributions for $REPOS ..."


RELENG_DIR=`mktemp -q -d -t ${PROJ}-releng.XXX`
if [ $? -ne 0 ]; then
        echo "$0: Can't create temp directory, exiting..."
        exit 1
fi

echo Using temp directory $RELENG_DIR ...
echo

DEST_DIR=`pwd`
SOURCE_DIR="${RELENG_DIR}/source"
BUILD_DIR="${RELENG_DIR}/build"

svn co -q $REPOS $SOURCE_DIR || exit 1

mkdir $BUILD_DIR || exit 1
cd $BUILD_DIR || exit 1

CMAKE_ARGS="\
	-DCMAKE_BUILD_TYPE=Release
	-DBoost_USE_STATIC_LIBS=ON
	-DGSL_USE_STATIC_LIBS=ON
"

if test $build_mingw32; then
	$CMAKE $SOURCE_DIR ${CMAKE_ARGS} \
		-DCMAKE_TOOLCHAIN_FILE="${SOURCE_DIR}/releng/mingw32.cmake"
elif test $build_m32; then
	$CMAKE $SOURCE_DIR ${CMAKE_ARGS} \
		-DCMAKE_CXX_FLAGS=-m32 \
		-DCMAKE_C_FLAGS=-m32
else
	$CMAKE $SOURCE_DIR ${CMAKE_ARGS}
fi

$MAKE
$MAKE package
$MAKE package_source

echo
echo Moving distribution packages ...

mv -v $PROJ_DISTS $DEST_DIR

echo
echo Cleaning up ...

cd $DEST_DIR
rm -rf $RELENG_DIR
