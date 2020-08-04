#!/bin/bash

LANGUAGE=$1

PATH=$PATH:/opt/local/bin
export PATH

OPTIONS=""

if [ "$LANGUAGE" == "C++98" ]; then
  echo "Building with C++98";
elif [ "$LANGUAGE" == "C++11" ]; then
  echo "Building with C++11";
  OPTIONS="$OPTIONS -DCMAKE_CXX_FLAGS='-std=c++11' "
fi


#########################################
# Build release version
#########################################

if [ ! -d build-release ]; then
  mkdir build-release
fi

cd build-release

cmake -DCMAKE_BUILD_TYPE=Release -DOPENMESH_BUILD_UNIT_TESTS=TRUE $OPTIONS ../

#build it
make

#build the unit tests
make unittests


#########################################
# Run Release Unittests
#########################################
cd Unittests

#execute tests
./unittests --gtest_color=yes --gtest_output=xml

cd ..
cd ..


#########################################
# Build Debug version and Unittests
#########################################

if [ ! -d build-debug ]; then
  mkdir build-debug
fi

cd build-debug

cmake -DCMAKE_BUILD_TYPE=Debug -DOPENMESH_BUILD_UNIT_TESTS=TRUE $OPTIONS ../

#build the unit tests
make unittests


#########################################
# Run Debug Unittests
#########################################

cd Unittests

# Run the unittests
./unittests --gtest_color=yes --gtest_output=xml


