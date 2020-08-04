#include <gtest/gtest.h>

#include "unittests_common.hh"
#include "unittests_basics.hh"
#include "unittests_iterators.hh"
#include "unittests_properties.hh"
#include "unittests_files.hh"

int main(int _argc, char** _argv) {

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();
}
