#!/bin/bash
BASE_DIR="`dirname \"$0\"`"
BASE_DIR="`( cd \"$BASE_DIR/..\" && pwd )`"

set -e

# build dependencies:
## blitz
cd $BASE_DIR/libraries
tar -xzf blitz-1.0.1.tar.gz
cd blitz-1.0.1
./configure --with-pic --prefix=`pwd`/install
make && make install

## cvodes
cd $BASE_DIR/libraries
tar -xzf cvodes-2.7.0.tar.gz
cd cvodes-2.7.0
./configure --with-pic  --prefix=`pwd`/install
make && make install

# run example
#$BASE_DIR/models/BirthDeath.py
cd $BASE_DIR/models
python3 ./BirthDeath.py
