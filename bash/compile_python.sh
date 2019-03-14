#!/bin/bash
#
# Purpose:
#
# This script will build and install Python and pip and put the new
# Python in the linker's path via ld.so.cond.d.  This script targets
# RedHat 6.
#
# Usage:
#
# 1. Become root
# 2. Set PYTHON_DEST_DIR below.
# 3. Set Python MAJVER suffix to match whatever python is calling itself: i.e. python$MAJVER
# 4. Download and untar Python
# 5. Run this script in the top level Python source directory
#
# Tips:
# Useful dependencies
#
#          BZip         libz        Open SSL       NCurses                    LZMA compression  
# RedHat   bzip-devel   libz-devel  openssl-devel  ncurses-devel  gdbm-devel  xz-devel          sqlite-devel  tk-devel  readline-devel
# Ubuntu   ?            ?
#

PYTHON_DEST_DIR=/usr/local/python-3.5.3
PYTHON_LIB_DIR=$PYTHON_DEST_DIR/lib
PYTHON_BIN_DIR=$PYTHON_DEST_DIR/bin
LD_CONF_FILE=/etc/ld.so.conf.d/python353.conf
MAJVER=3

if [ -d $PYTHON_DEST_DIR ]; then
    echo "WARNING: Python destination directory already exists"
else 
    echo "Making Python destination directory $PYTHON_DEST_DIR"
    mkdir $PYTHON_DEST_DIR
fi

if [ -d $PYTHON_LIB_DIR ]; then
    echo "WARNING: Python lib directory already exists"
else
    echo "Making lib directory $PYTHON_LIB_DIR"
    mkdir $PYTHON_LIB_DIR
fi

if [ -f $LD_CONF_FILE ]; then
    echo "WARNING: $LD_CONF_FILE already exists"
else
    echo "Making ld conf file $LD_CONF_FILE"
    echo $PYTHON_LIB_DIR > $LD_CONF_FILE
    echo "Running ldconfig"
    ldconfig
fi

if [ -f $PYTHON_BIN_DIR/python$MAJVER ]; then
    echo "WARNING: python already exists"
else
    echo "Building python and installing pip"
    ./configure --prefix=$PYTHON_DEST_DIR --enable-shared LDFLAGS="-Wl,-rpath $PYTHON_LIB_DIR"
    make
    make install

    cd $PYTHON_BIN_DIR

    ./python$MAJVER -m ensurepip --upgrade

    ./pip$MAJVER install -U pip
fi
