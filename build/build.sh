#!/bin/bash
pushd DEBUG
gcc -g -O0 -c ../../models/Wang2010/C/*.c
gcc -g -O0 -c ../../src/computeLinearNoise.cpp ../../src/main.cpp -I ../../include/ -I ../../libraries/blitz-1.0.1 -I$HOME/include  -I ../../models/Wang2010/C
gcc *.o -L$HOME/lib -L../../libraries/install/blitz-1.0.1/lib/ -lsundials_cvodes -lblitz -lsundials_nvecserial -lstdc++ -o Wang2010_DEBUG
popd

pushd RUN
gcc -c ../../models/Wang2010/C/*.c
gcc -c ../../src/computeLinearNoise.cpp ../../src/main.cpp -I ../../include/ -I ../../libraries/blitz-1.0.1 -I$HOME/include  -I ../../models/Wang2010/C
gcc *.o -L$HOME/lib -L../../libraries/install/blitz-1.0.1/lib/ -lsundials_cvodes -lblitz -lsundials_nvecserial -lstdc++ -o Wang2010_RUN
popd

