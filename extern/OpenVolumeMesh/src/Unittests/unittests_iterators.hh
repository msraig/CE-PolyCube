#ifndef INCLUDE_UNITTESTS_ITERATORS_HH
#define INCLUDE_UNITTESTS_ITERATORS_HH

#include <gtest/gtest.h>
#include <Unittests/unittests_common.hh>

TEST_F(HexahedralMeshBase, HexVertexIterTest) {

    generateHexahedralMesh(mesh_);

    HexVertexIter hv_it = mesh_.hv_iter(CellHandle(0));

    EXPECT_TRUE(hv_it.valid());

    EXPECT_EQ(VertexHandle(0), *hv_it); ++hv_it;
    EXPECT_EQ(VertexHandle(1), *hv_it); ++hv_it;
    EXPECT_EQ(VertexHandle(2), *hv_it); ++hv_it;
    EXPECT_EQ(VertexHandle(3), *hv_it); ++hv_it;
    EXPECT_EQ(VertexHandle(4), *hv_it); ++hv_it;
    EXPECT_EQ(VertexHandle(7), *hv_it); ++hv_it;
    EXPECT_EQ(VertexHandle(6), *hv_it); ++hv_it;
    EXPECT_EQ(VertexHandle(5), *hv_it);
}

#endif // INCLUDE GUARD
