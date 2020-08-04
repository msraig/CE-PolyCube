#!/bin/bash

COMPILER=$1
LANGUAGE=$2

OPTIONS=""

if [ "$COMPILER" == "gcc" ]; then
  echo "Building with GCC";
elif [ "$COMPILER" == "clang" ]; then

  OPTIONS="$OPTIONS -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang"
  echo "Building with CLANG";
fi  

if [ "$LANGUAGE" == "C++98" ]; then
  echo "Building with C++98";
elif [ "$LANGUAGE" == "C++11" ]; then
  echo "Building with C++11";
  OPTIONS="$OPTIONS -DCMAKE_CXX_FLAGS='-std=c++11' "
fi  

#########################################

# Make release build folder
if [ ! -d build-release ]; then
  mkdir build-release
fi

cd build-release

cmake -DCMAKE_BUILD_TYPE=Release -DOPENMESH_BUILD_UNIT_TESTS=TRUE -DSTL_VECTOR_CHECKS=ON -DOPENMESH_BUILD_PYTHON_UNIT_TESTS=ON $OPTIONS ../

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

# Execute Python unittests
cd Python-Unittests
python -m unittest discover -v

cd ..
cd ..


#########################################
# Build Debug version and Unittests
#########################################

if [ ! -d build-debug ]; then
  mkdir build-debug
fi

cd build-debug

cmake -DCMAKE_BUILD_TYPE=Debug -DOPENMESH_BUILD_UNIT_TESTS=TRUE -DSTL_VECTOR_CHECKS=ON -DOPENMESH_BUILD_PYTHON_UNIT_TESTS=ON $OPTIONS ../

#build the unit tests
make unittests


#########################################
# Run Debug Unittests
#########################################

cd Unittests

# Run the unittests
./unittests --gtest_color=yes --gtest_output=xml


