#ifndef INCLUDE_UNITTESTS_FILES_HH
#define INCLUDE_UNITTESTS_FILES_HH

#include <gtest/gtest.h>
#include <Unittests/unittests_common.hh>

#include <OpenVolumeMesh/FileManager/FileManager.hh>

TEST_F(PolyhedralMeshBase, LoadFile) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("Cylinder.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());
}

TEST_F(PolyhedralMeshBase, LoadNonManifoldMesh) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("NonManifold.ovm", mesh_));

  EXPECT_EQ(12u, mesh_.n_vertices());
  EXPECT_EQ(20u, mesh_.n_edges());
  EXPECT_EQ(11u, mesh_.n_faces());
  EXPECT_EQ(2u, mesh_.n_cells());
}

TEST_F(HexahedralMeshBase, LoadFile) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.isHexahedralMesh("Cylinder.ovm"));
  ASSERT_TRUE(fileManager.readFile("Cylinder.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());
}

TEST_F(PolyhedralMeshBase, SaveFile) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("Cylinder.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());

  // Write file
  ASSERT_TRUE(fileManager.writeFile("Cylinder.copy.ovm", mesh_));

  mesh_.clear();

  ASSERT_TRUE(fileManager.readFile("Cylinder.copy.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());
}

TEST_F(PolyhedralMeshBase, SaveFileWithProps) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("Cylinder.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());

  // Attach non-persistent properties
  HalfFacePropertyT<float> hfprop = mesh_.request_halfface_property<float>("MyHalfFaceProp");
  VertexPropertyT<unsigned int> vprop = mesh_.request_vertex_property<unsigned int>("MyVertexProp");

  for(unsigned int i = 0; i < mesh_.n_halffaces(); ++i) {
      hfprop[i] = (float)i/2.0f;
  }
  for(unsigned int i = 0; i < mesh_.n_vertices(); ++i) {
      vprop[i] = i;
  }

  mesh_.set_persistent(hfprop);
  mesh_.set_persistent(vprop);

  // Write file
  ASSERT_TRUE(fileManager.writeFile("Cylinder.copy.ovm", mesh_));

  mesh_.clear();

  ASSERT_TRUE(fileManager.readFile("Cylinder.copy.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());

  EXPECT_EQ(1u, mesh_.n_halfface_props());
  EXPECT_EQ(1u, mesh_.n_vertex_props());

  HalfFacePropertyT<float> hfprop2 = mesh_.request_halfface_property<float>("MyHalfFaceProp");
  VertexPropertyT<unsigned int> vprop2 = mesh_.request_vertex_property<unsigned int>("MyVertexProp");

  for(unsigned int i = 0; i < mesh_.n_halffaces(); ++i) {
      EXPECT_FLOAT_EQ((float)i/2.0f, hfprop2[i]);
  }
  for(unsigned int i = 0; i < mesh_.n_vertices(); ++i) {
      EXPECT_EQ(i, vprop2[i]);
  }
}

TEST_F(PolyhedralMeshBase, LoadFileWithProps) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("Cube_with_props.ovm", mesh_));

  EXPECT_EQ(8u, mesh_.n_vertices());
  EXPECT_EQ(12u, mesh_.n_edges());
  EXPECT_EQ(6u, mesh_.n_faces());
  EXPECT_EQ(1u, mesh_.n_cells());

  EXPECT_EQ(1u, mesh_.n_vertex_props());
  EXPECT_EQ(1u, mesh_.n_edge_props());
  EXPECT_EQ(0u, mesh_.n_halfedge_props());
  EXPECT_EQ(1u, mesh_.n_face_props());
  EXPECT_EQ(1u, mesh_.n_halfface_props());
  EXPECT_EQ(0u, mesh_.n_cell_props());
}

TEST_F(PolyhedralMeshBase, SaveFileWithProps2) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("Cube_with_props.ovm", mesh_));

  EXPECT_EQ(8u, mesh_.n_vertices());
  EXPECT_EQ(12u, mesh_.n_edges());
  EXPECT_EQ(6u, mesh_.n_faces());
  EXPECT_EQ(1u, mesh_.n_cells());

  EXPECT_EQ(1u, mesh_.n_vertex_props());
  EXPECT_EQ(1u, mesh_.n_edge_props());
  EXPECT_EQ(0u, mesh_.n_halfedge_props());
  EXPECT_EQ(1u, mesh_.n_face_props());
  EXPECT_EQ(1u, mesh_.n_halfface_props());
  EXPECT_EQ(0u, mesh_.n_cell_props());

  ASSERT_TRUE(fileManager.writeFile("Cube_with_props.copy.ovm", mesh_));

  mesh_.clear();

  ASSERT_TRUE(fileManager.readFile("Cube_with_props.copy.ovm", mesh_));

  EXPECT_EQ(8u, mesh_.n_vertices());
  EXPECT_EQ(12u, mesh_.n_edges());
  EXPECT_EQ(6u, mesh_.n_faces());
  EXPECT_EQ(1u, mesh_.n_cells());

  EXPECT_EQ(1u, mesh_.n_vertex_props());
  EXPECT_EQ(1u, mesh_.n_edge_props());
  EXPECT_EQ(0u, mesh_.n_halfedge_props());
  EXPECT_EQ(1u, mesh_.n_face_props());
  EXPECT_EQ(1u, mesh_.n_halfface_props());
  EXPECT_EQ(0u, mesh_.n_cell_props());
}

#endif // INCLUDE GUARD
