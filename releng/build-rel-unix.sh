#!/bin/sh

MAKE=make
CMAKE=cmake
REPOS=svn://scit.us/dawg/stable

RELENG_DIR=`mktemp -q -d -t dawg-releng`
if [ $? -ne 0 ]; then
        echo "$0: Can't create temp directory, exiting..."
        exit 1
fi

DEST_DIR=`pwd`
SOURCE_DIR="${RELENG_DIR}/source"
BUILD_DIR="${RELENG_DIR}/build"

svn co -q $REPOS $SOURCE_DIR || exit 1

mkdir $BUILD_DIR || exit 1
cd $BUILD_DIR || exit 1

$CMAKE $SOURCE_DIR -DCMAKE_BUILD_TYPE=Release
$MAKE
$MAKE package
$MAKE package_source

mv dawg-1* $DEST_DIR

cd $DEST_DIR
rm -rf $RELENG_DIR

