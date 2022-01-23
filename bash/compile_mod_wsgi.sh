#!/bin/bash

PYTHON=/usr/local/python-3.6.1/bin/python3
WSGIPREFIX=/usr/local/wsgi-4.6.5-py-3.6.1

./configure --prefix=$WSGIPREFIX --with-python=$PYTHON
make

echo "You might want to back up the old mod_wsgi.so before typing make install."
