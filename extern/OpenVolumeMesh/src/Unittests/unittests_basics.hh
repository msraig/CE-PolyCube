#ifndef INCLUDE_UNITTESTS_BASICS_HH
#define INCLUDE_UNITTESTS_BASICS_HH

#include <gtest/gtest.h>
#include <Unittests/unittests_common.hh>

#include <OpenVolumeMesh/Attribs/StatusAttrib.hh>
#include <OpenVolumeMesh/Attribs/NormalAttrib.hh>
#include <OpenVolumeMesh/Attribs/ColorAttrib.hh>

using namespace OpenVolumeMesh;

/*
 * ====================================================================
 * Define tests below
 * ====================================================================
 */

TEST_F(PolyhedralMeshBase, CreateSimpleMesh) {

  /*
   * Add vertices
   */

  VertexHandle v0 = mesh_.add_vertex(Vec3d(-1.0, -1.0, -1.0));
  VertexHandle v1 = mesh_.add_vertex(Vec3d( 1.0, -1.0, -1.0));
  VertexHandle v2 = mesh_.add_vertex(Vec3d( 1.0,  1.0, -1.0));
  VertexHandle v3 = mesh_.add_vertex(Vec3d(-1.0,  1.0, -1.0));
  VertexHandle v4 = mesh_.add_vertex(Vec3d(-1.0, -1.0,  1.0));
  VertexHandle v5 = mesh_.add_vertex(Vec3d( 1.0, -1.0,  1.0));
  VertexHandle v6 = mesh_.add_vertex(Vec3d( 1.0,  1.0,  1.0));
  VertexHandle v7 = mesh_.add_vertex(Vec3d(-1.0,  1.0,  1.0));

  EXPECT_EQ(8u, mesh_.n_vertices()) << "The number of vertices is not correct!";

  /*
   * Add faces
   */

  std::vector<VertexHandle> fvertices;

  fvertices.push_back(v3);
  fvertices.push_back(v2);
  fvertices.push_back(v1);
  fvertices.push_back(v0);

  FaceHandle fh0 = mesh_.add_face(fvertices);

  fvertices.clear();

  fvertices.push_back(v4);
  fvertices.push_back(v5);
  fvertices.push_back(v6);
  fvertices.push_back(v7);

  FaceHandle fh1 = mesh_.add_face(fvertices);

  fvertices.clear();

  fvertices.push_back(v0);
  fvertices.push_back(v4);
  fvertices.push_back(v7);
  fvertices.push_back(v3);

  FaceHandle fh2 = mesh_.add_face(fvertices);

  fvertices.clear();

  fvertices.push_back(v1);
  fvertices.push_back(v2);
  fvertices.push_back(v6);
  fvertices.push_back(v5);

  FaceHandle fh3 = mesh_.add_face(fvertices);

  fvertices.clear();

  fvertices.push_back(v7);
  fvertices.push_back(v6);
  fvertices.push_back(v2);
  fvertices.push_back(v3);

  FaceHandle fh4 = mesh_.add_face(fvertices);

  fvertices.clear();

  fvertices.push_back(v0);
  fvertices.push_back(v1);
  fvertices.push_back(v5);
  fvertices.push_back(v4);

  FaceHandle fh5 = mesh_.add_face(fvertices);

  EXPECT_EQ(12u, mesh_.n_edges()) << "The number of edges is not correct!";
  EXPECT_EQ(6u, mesh_.n_faces())  << "The number of faces is not correct!";

  /*
   * Add cell
   */

  std::vector<HalfFaceHandle> chfaces;

  chfaces.push_back(mesh_.halfface_handle(fh0, 0));
  chfaces.push_back(mesh_.halfface_handle(fh1, 0));
  chfaces.push_back(mesh_.halfface_handle(fh2, 0));
  chfaces.push_back(mesh_.halfface_handle(fh3, 0));
  chfaces.push_back(mesh_.halfface_handle(fh4, 0));
  chfaces.push_back(mesh_.halfface_handle(fh5, 0));

  mesh_.add_cell(chfaces);

  EXPECT_EQ(1u, mesh_.n_cells())  << "The number of cells is not correct!";
}

//===========================================================================

TEST_F(PolyhedralMeshBase, CreateSimpleMeshWithoutCells) {

    Vec3d p1(0.0, 0.0, 0.0);
    Vec3d p2(1.0, 0.0, 0.0);
    Vec3d p3(1.0, 1.0, 0.0);
    Vec3d p4(0.0, 1.0, 0.0);

    Vec3d p5(0.0, 0.0, 1.0);
    Vec3d p6(1.0, 0.0, 1.0);
    Vec3d p7(1.0, 1.0, 1.0);
    Vec3d p8(0.0, 1.0, 1.0);

    VertexHandle v1 = mesh_.add_vertex(p1);
    VertexHandle v2 = mesh_.add_vertex(p2);
    VertexHandle v3 = mesh_.add_vertex(p3);
    VertexHandle v4 = mesh_.add_vertex(p4);

    VertexHandle v5 = mesh_.add_vertex(p5);
    VertexHandle v6 = mesh_.add_vertex(p6);
    VertexHandle v7 = mesh_.add_vertex(p7);
    VertexHandle v8 = mesh_.add_vertex(p8);

    EXPECT_EQ(0, v1);
    EXPECT_EQ(1, v2);
    EXPECT_EQ(2, v3);
    EXPECT_EQ(3, v4);
    EXPECT_EQ(4, v5);
    EXPECT_EQ(5, v6);
    EXPECT_EQ(6, v7);
    EXPECT_EQ(7, v8);

    EdgeHandle e1 = mesh_.add_edge(v1, v2);
    EdgeHandle e2 = mesh_.add_edge(v2, v3);
    EdgeHandle e3 = mesh_.add_edge(v3, v4);
    EdgeHandle e4 = mesh_.add_edge(v4, v1);

    EdgeHandle e5 = mesh_.add_edge(v5, v6);
    EdgeHandle e6 = mesh_.add_edge(v6, v7);
    EdgeHandle e7 = mesh_.add_edge(v7, v8);
    EdgeHandle e8 = mesh_.add_edge(v8, v5);

    EXPECT_EQ(0, e1);
    EXPECT_EQ(1, e2);
    EXPECT_EQ(2, e3);
    EXPECT_EQ(3, e4);
    EXPECT_EQ(4, e5);
    EXPECT_EQ(5, e6);
    EXPECT_EQ(6, e7);
    EXPECT_EQ(7, e8);

    // Get halfedges
    HalfEdgeHandle h1 = mesh_.halfedge_handle(e1, 0u);
    HalfEdgeHandle h2 = mesh_.halfedge_handle(e2, 0u);
    HalfEdgeHandle h3 = mesh_.halfedge_handle(e3, 0u);
    HalfEdgeHandle h4 = mesh_.halfedge_handle(e4, 0u);

    HalfEdgeHandle h5 = mesh_.halfedge_handle(e5, 0u);
    HalfEdgeHandle h6 = mesh_.halfedge_handle(e6, 0u);
    HalfEdgeHandle h7 = mesh_.halfedge_handle(e7, 0u);
    HalfEdgeHandle h8 = mesh_.halfedge_handle(e8, 0u);

    EXPECT_EQ(v1, mesh_.halfedge(h1).from_vertex());
    EXPECT_EQ(v2, mesh_.halfedge(h1).to_vertex());
    EXPECT_EQ(v2, mesh_.halfedge(h2).from_vertex());
    EXPECT_EQ(v3, mesh_.halfedge(h2).to_vertex());
    EXPECT_EQ(v3, mesh_.halfedge(h3).from_vertex());
    EXPECT_EQ(v4, mesh_.halfedge(h3).to_vertex());
    EXPECT_EQ(v4, mesh_.halfedge(h4).from_vertex());
    EXPECT_EQ(v1, mesh_.halfedge(h4).to_vertex());

    EXPECT_EQ(v5, mesh_.halfedge(h5).from_vertex());
    EXPECT_EQ(v6, mesh_.halfedge(h5).to_vertex());
    EXPECT_EQ(v6, mesh_.halfedge(h6).from_vertex());
    EXPECT_EQ(v7, mesh_.halfedge(h6).to_vertex());
    EXPECT_EQ(v7, mesh_.halfedge(h7).from_vertex());
    EXPECT_EQ(v8, mesh_.halfedge(h7).to_vertex());
    EXPECT_EQ(v8, mesh_.halfedge(h8).from_vertex());
    EXPECT_EQ(v5, mesh_.halfedge(h8).to_vertex());

    // Check opposite halfedges
    EXPECT_EQ(v2, mesh_.opposite_halfedge(h1).from_vertex());
    EXPECT_EQ(v1, mesh_.opposite_halfedge(h1).to_vertex());
    EXPECT_EQ(v3, mesh_.opposite_halfedge(h2).from_vertex());
    EXPECT_EQ(v2, mesh_.opposite_halfedge(h2).to_vertex());
    EXPECT_EQ(v4, mesh_.opposite_halfedge(h3).from_vertex());
    EXPECT_EQ(v3, mesh_.opposite_halfedge(h3).to_vertex());
    EXPECT_EQ(v1, mesh_.opposite_halfedge(h4).from_vertex());
    EXPECT_EQ(v4, mesh_.opposite_halfedge(h4).to_vertex());

    EXPECT_EQ(v6, mesh_.opposite_halfedge(h5).from_vertex());
    EXPECT_EQ(v5, mesh_.opposite_halfedge(h5).to_vertex());
    EXPECT_EQ(v7, mesh_.opposite_halfedge(h6).from_vertex());
    EXPECT_EQ(v6, mesh_.opposite_halfedge(h6).to_vertex());
    EXPECT_EQ(v8, mesh_.opposite_halfedge(h7).from_vertex());
    EXPECT_EQ(v7, mesh_.opposite_halfedge(h7).to_vertex());
    EXPECT_EQ(v5, mesh_.opposite_halfedge(h8).from_vertex());
    EXPECT_EQ(v8, mesh_.opposite_halfedge(h8).to_vertex());

    // Add a face via vertices
    std::vector<VertexHandle> vertices;
    vertices.push_back(v2); vertices.push_back(v6);
    vertices.push_back(v7); vertices.push_back(v3);
    FaceHandle f1 = mesh_.add_face(vertices);

    EXPECT_EQ(0, f1);

    // Get halfedges of face
    std::vector<HalfEdgeHandle> halfedges = mesh_.face(f1).halfedges();

    std::vector<HalfEdgeHandle>::iterator it = halfedges.begin();

    EXPECT_EQ(8, mesh_.edge_handle(*it)); ++it;
    EXPECT_EQ(5, mesh_.edge_handle(*it)); ++it;
    EXPECT_EQ(9, mesh_.edge_handle(*it)); ++it;
    EXPECT_EQ(1, mesh_.edge_handle(*it));

    // Add invalid face
    halfedges.clear();
    halfedges.push_back(mesh_.halfedge_handle(e1, 0)); halfedges.push_back(mesh_.halfedge_handle(e2, 0));
    halfedges.push_back(mesh_.halfedge_handle(e7, 0)); halfedges.push_back(mesh_.halfedge_handle(e4, 0));

    FaceHandle fI = mesh_.add_face(halfedges, true);

    EXPECT_EQ(PolyhedralMesh::InvalidFaceHandle, fI);

    // Now add valid face via edges
    halfedges.clear();
    halfedges.push_back(mesh_.halfedge_handle(e1, 0)); halfedges.push_back(mesh_.halfedge_handle(e2, 0));
    halfedges.push_back(mesh_.halfedge_handle(e3, 0)); halfedges.push_back(mesh_.halfedge_handle(e4, 0));

    FaceHandle f2 = mesh_.add_face(halfedges);

    EXPECT_EQ(1, f2);

    // Get halfedges of face
    halfedges = mesh_.face(f2).halfedges();
    int handle = 0;
    for(it = halfedges.begin(); it != halfedges.end(); ++it) {
        EXPECT_EQ(handle, mesh_.edge_handle(*it)); handle++;
    }
}

TEST_F(PolyhedralMeshBase, TopologyCheckPass) {

    // Add eight vertices
    VertexHandle v0 = mesh_.add_vertex(Vec3d(-1.0, 0.0, 0.0));
    VertexHandle v1 = mesh_.add_vertex(Vec3d( 0.0, 0.0, 1.0));
    VertexHandle v2 = mesh_.add_vertex(Vec3d( 1.0, 0.0, 0.0));
    VertexHandle v3 = mesh_.add_vertex(Vec3d( 0.0, 1.0, 0.0));

    std::vector<VertexHandle> vertices;

    // Add faces
    vertices.push_back(v0); vertices.push_back(v1);vertices.push_back(v3);
    FaceHandle f0 = mesh_.add_face(vertices);

    vertices.clear();
    vertices.push_back(v1); vertices.push_back(v2);vertices.push_back(v3);
    FaceHandle f1 = mesh_.add_face(vertices);

    vertices.clear();
    vertices.push_back(v0); vertices.push_back(v1);vertices.push_back(v2);
    FaceHandle f2 = mesh_.add_face(vertices);

    vertices.clear();
    vertices.push_back(v0); vertices.push_back(v3);vertices.push_back(v2);
    FaceHandle f3 = mesh_.add_face(vertices);

    std::vector<HalfFaceHandle> halffaces;

    // Add first tetrahedron
    halffaces.push_back(mesh_.halfface_handle(f0, 1));
    halffaces.push_back(mesh_.halfface_handle(f1, 1));
    halffaces.push_back(mesh_.halfface_handle(f2, 0));
    halffaces.push_back(mesh_.halfface_handle(f3, 1));
    EXPECT_NE(PolyhedralMesh::InvalidCellHandle, mesh_.add_cell(halffaces));
}

TEST_F(PolyhedralMeshBase, TopologyCheckFail) {

    // Add eight vertices
    VertexHandle v0 = mesh_.add_vertex(Vec3d(-1.0, 0.0, 0.0));
    VertexHandle v1 = mesh_.add_vertex(Vec3d( 0.0, 0.0, 1.0));
    VertexHandle v2 = mesh_.add_vertex(Vec3d( 1.0, 0.0, 0.0));
    VertexHandle v3 = mesh_.add_vertex(Vec3d( 0.0, 1.0, 0.0));

    std::vector<VertexHandle> vertices;

    // Add faces
    vertices.push_back(v0); vertices.push_back(v1);vertices.push_back(v3);
    FaceHandle f0 = mesh_.add_face(vertices);

    vertices.clear();
    vertices.push_back(v1); vertices.push_back(v2);vertices.push_back(v3);
    FaceHandle f1 = mesh_.add_face(vertices);

    vertices.clear();
    vertices.push_back(v0); vertices.push_back(v1);vertices.push_back(v2);
    FaceHandle f2 = mesh_.add_face(vertices);

    vertices.clear();
    vertices.push_back(v0); vertices.push_back(v3);vertices.push_back(v2);
    FaceHandle f3 = mesh_.add_face(vertices);

    std::vector<HalfFaceHandle> halffaces;

    // Add first tetrahedron
    halffaces.push_back(mesh_.halfface_handle(f0, 1));
    halffaces.push_back(mesh_.halfface_handle(f1, 1));
    halffaces.push_back(mesh_.halfface_handle(f2, 0));
    halffaces.push_back(mesh_.halfface_handle(f3, 0));
    EXPECT_EQ(PolyhedralMesh::InvalidCellHandle, mesh_.add_cell(halffaces, true));
}

TEST_F(HexahedralMeshBase, TopologyCheckPass) {

    VertexHandle v0 = mesh_.add_vertex(Vec3d(-1.0, -1.0, -1.0));
    VertexHandle v1 = mesh_.add_vertex(Vec3d( 1.0, -1.0, -1.0));
    VertexHandle v2 = mesh_.add_vertex(Vec3d( 1.0,  1.0, -1.0));
    VertexHandle v3 = mesh_.add_vertex(Vec3d(-1.0,  1.0, -1.0));
    VertexHandle v4 = mesh_.add_vertex(Vec3d(-1.0, -1.0,  1.0));
    VertexHandle v5 = mesh_.add_vertex(Vec3d( 1.0, -1.0,  1.0));
    VertexHandle v6 = mesh_.add_vertex(Vec3d( 1.0,  1.0,  1.0));
    VertexHandle v7 = mesh_.add_vertex(Vec3d(-1.0,  1.0,  1.0));

    std::vector<VertexHandle> fvertices;

    fvertices.push_back(v3);
    fvertices.push_back(v2);
    fvertices.push_back(v1);
    fvertices.push_back(v0);

    FaceHandle fh0 = mesh_.add_face(fvertices);

    fvertices.clear();

    fvertices.push_back(v4);
    fvertices.push_back(v5);
    fvertices.push_back(v6);
    fvertices.push_back(v7);

    FaceHandle fh1 = mesh_.add_face(fvertices);

    fvertices.clear();

    fvertices.push_back(v0);
    fvertices.push_back(v4);
    fvertices.push_back(v7);
    fvertices.push_back(v3);

    FaceHandle fh2 = mesh_.add_face(fvertices);

    fvertices.clear();

    fvertices.push_back(v1);
    fvertices.push_back(v2);
    fvertices.push_back(v6);
    fvertices.push_back(v5);

    FaceHandle fh3 = mesh_.add_face(fvertices);

    fvertices.clear();

    fvertices.push_back(v7);
    fvertices.push_back(v6);
    fvertices.push_back(v2);
    fvertices.push_back(v3);

    FaceHandle fh4 = mesh_.add_face(fvertices);

    fvertices.clear();

    fvertices.push_back(v0);
    fvertices.push_back(v1);
    fvertices.push_back(v5);
    fvertices.push_back(v4);

    FaceHandle fh5 = mesh_.add_face(fvertices);

    std::vector<HalfFaceHandle> chfaces;

    chfaces.push_back(mesh_.halfface_handle(fh0, 0));
    chfaces.push_back(mesh_.halfface_handle(fh1, 0));
    chfaces.push_back(mesh_.halfface_handle(fh2, 0));
    chfaces.push_back(mesh_.halfface_handle(fh3, 0));
    chfaces.push_back(mesh_.halfface_handle(fh4, 0));
    chfaces.push_back(mesh_.halfface_handle(fh5, 0));

    EXPECT_NE(HexahedralMesh::InvalidCellHandle, mesh_.add_cell(chfaces, true));
}

TEST_F(HexahedralMeshBase, TopologyCheckFail) {

    VertexHandle v0 = mesh_.add_vertex(Vec3d(-1.0, -1.0, -1.0));
    VertexHandle v1 = mesh_.add_vertex(Vec3d( 1.0, -1.0, -1.0));
    VertexHandle v2 = mesh_.add_vertex(Vec3d( 1.0,  1.0, -1.0));
    VertexHandle v3 = mesh_.add_vertex(Vec3d(-1.0,  1.0, -1.0));
    VertexHandle v4 = mesh_.add_vertex(Vec3d(-1.0, -1.0,  1.0));
    VertexHandle v5 = mesh_.add_vertex(Vec3d( 1.0, -1.0,  1.0));
    VertexHandle v6 = mesh_.add_vertex(Vec3d( 1.0,  1.0,  1.0));
    VertexHandle v7 = mesh_.add_vertex(Vec3d(-1.0,  1.0,  1.0));

    std::vector<VertexHandle> fvertices;

    fvertices.push_back(v3);
    fvertices.push_back(v2);
    fvertices.push_back(v1);
    fvertices.push_back(v0);

    FaceHandle fh0 = mesh_.add_face(fvertices);

    fvertices.clear();

    fvertices.push_back(v4);
    fvertices.push_back(v5);
    fvertices.push_back(v6);
    fvertices.push_back(v7);

    FaceHandle fh1 = mesh_.add_face(fvertices);

    fvertices.clear();

    fvertices.push_back(v0);
    fvertices.push_back(v4);
    fvertices.push_back(v7);
    fvertices.push_back(v3);

    FaceHandle fh2 = mesh_.add_face(fvertices);

    fvertices.clear();

    fvertices.push_back(v1);
    fvertices.push_back(v2);
    fvertices.push_back(v6);
    fvertices.push_back(v5);

    FaceHandle fh3 = mesh_.add_face(fvertices);

    fvertices.clear();

    fvertices.push_back(v7);
    fvertices.push_back(v6);
    fvertices.push_back(v2);
    fvertices.push_back(v3);

    FaceHandle fh4 = mesh_.add_face(fvertices);

    fvertices.clear();

    fvertices.push_back(v0);
    fvertices.push_back(v1);
    fvertices.push_back(v5);
    fvertices.push_back(v4);

    FaceHandle fh5 = mesh_.add_face(fvertices);

    std::vector<HalfFaceHandle> chfaces;

    chfaces.push_back(mesh_.halfface_handle(fh0, 0));
    chfaces.push_back(mesh_.halfface_handle(fh1, 1));
    chfaces.push_back(mesh_.halfface_handle(fh2, 0));
    chfaces.push_back(mesh_.halfface_handle(fh3, 0));
    chfaces.push_back(mesh_.halfface_handle(fh4, 0));
    chfaces.push_back(mesh_.halfface_handle(fh5, 0));

    EXPECT_EQ(HexahedralMesh::InvalidCellHandle, mesh_.add_cell(chfaces, true));
}

TEST_F(PolyhedralMeshBase, VolumeMeshConnectivity) {

    generatePolyhedralMesh(mesh_);

    // Add invalid cell
    std::vector<HalfFaceHandle> hfaces;
    hfaces.push_back(1); hfaces.push_back(5);
    hfaces.push_back(7); hfaces.push_back(9);
    hfaces.push_back(10); hfaces.push_back(21);
    CellHandle i_cell = mesh_.add_cell(hfaces, true);

    EXPECT_EQ(PolyhedralMesh::InvalidCellHandle, i_cell);

    EXPECT_EQ(CellHandle(0), mesh_.incident_cell(1));
    EXPECT_EQ(CellHandle(0), mesh_.incident_cell(2));
    EXPECT_EQ(CellHandle(0), mesh_.incident_cell(5));
    EXPECT_EQ(CellHandle(0), mesh_.incident_cell(7));
    EXPECT_EQ(CellHandle(0), mesh_.incident_cell(9));
    EXPECT_EQ(CellHandle(0), mesh_.incident_cell(10));

    EXPECT_EQ(CellHandle(1), mesh_.incident_cell(3));
    EXPECT_EQ(CellHandle(1), mesh_.incident_cell(12));
    EXPECT_EQ(CellHandle(1), mesh_.incident_cell(15));
    EXPECT_EQ(CellHandle(1), mesh_.incident_cell(17));
    EXPECT_EQ(CellHandle(1), mesh_.incident_cell(19));
    EXPECT_EQ(CellHandle(1), mesh_.incident_cell(20));

    // Test adjacency function
    HalfFaceHandle ad_hf1 = mesh_.adjacent_halfface_in_cell(1, 3);
    // Should be halfface 5
    EXPECT_EQ(HalfFaceHandle(5), ad_hf1);

    HalfFaceHandle ad_hf2 = mesh_.adjacent_halfface_in_cell(1, 7);
    // Should be halfface 7
    EXPECT_EQ(HalfFaceHandle(7), ad_hf2);

    HalfFaceHandle ad_hf3 = mesh_.adjacent_halfface_in_cell(5, 24);
    // Should be invalid
    EXPECT_EQ(PolyhedralMesh::InvalidHalfFaceHandle, ad_hf3);

    HalfFaceHandle ad_hf4 = mesh_.adjacent_halfface_in_cell(12, 24);
    // Should be invalid
    EXPECT_EQ(HalfFaceHandle(20), ad_hf4);

    HalfFaceHandle ad_hf5 = mesh_.adjacent_halfface_in_cell(0, 0);
    // Should be invalid
    EXPECT_EQ(PolyhedralMesh::InvalidHalfFaceHandle, ad_hf5);

    EXPECT_EQ(12u, mesh_.n_vertices());
    EXPECT_EQ(20u, mesh_.n_edges());
    EXPECT_EQ(2u,  mesh_.n_cells());
    EXPECT_EQ(11u, mesh_.n_faces());
}

TEST_F(PolyhedralMeshBase, VolumeMeshNormals) {

    generatePolyhedralMesh(mesh_);

    NormalAttrib<GeometricPolyhedralMeshV3d> normals(mesh_);

    Vec3d n_x(1.0, 0.0, 0.0);
    Vec3d n_y(0.0, 1.0, 0.0);
    Vec3d n_z(0.0, 0.0, 1.0);

    normals.update_face_normals();

    // Should be positive z-axis
    Vec3d& n = normals[FaceHandle(0)];
    EXPECT_DOUBLE_EQ(n_z[0], n[0]);
    EXPECT_DOUBLE_EQ(n_z[1], n[1]);
    EXPECT_DOUBLE_EQ(n_z[2], n[2]);

    n = normals[HalfFaceHandle(1)];
    EXPECT_DOUBLE_EQ(-n_z[0], n[0]);
    EXPECT_DOUBLE_EQ(-n_z[1], n[1]);
    EXPECT_DOUBLE_EQ(-n_z[2], n[2]);

    // Should be negative x-axis
    n = normals[FaceHandle(2)];
    EXPECT_DOUBLE_EQ(-n_x[0], n[0]);
    EXPECT_DOUBLE_EQ(-n_x[1], n[1]);
    EXPECT_DOUBLE_EQ(-n_x[2], n[2]);

    n = normals[HalfFaceHandle(4)];
    EXPECT_DOUBLE_EQ(-n_x[0], n[0]);
    EXPECT_DOUBLE_EQ(-n_x[1], n[1]);
    EXPECT_DOUBLE_EQ(-n_x[2], n[2]);

    n = normals[HalfFaceHandle(5)];
    EXPECT_DOUBLE_EQ(n_x[0], n[0]);
    EXPECT_DOUBLE_EQ(n_x[1], n[1]);
    EXPECT_DOUBLE_EQ(n_x[2], n[2]);

    // Should be negative y-axis
    n = normals[FaceHandle(4)];
    EXPECT_DOUBLE_EQ(-n_y[0], n[0]);
    EXPECT_DOUBLE_EQ(-n_y[1], n[1]);
    EXPECT_DOUBLE_EQ(-n_y[2], n[2]);

    n = normals[HalfFaceHandle(9)];
    EXPECT_DOUBLE_EQ(n_y[0], n[0]);
    EXPECT_DOUBLE_EQ(n_y[1], n[1]);
    EXPECT_DOUBLE_EQ(n_y[2], n[2]);

    // Should be positive y-axis
    n = normals[FaceHandle(5)];
    EXPECT_DOUBLE_EQ(-n_y[0], n[0]);
    EXPECT_DOUBLE_EQ(-n_y[1], n[1]);
    EXPECT_DOUBLE_EQ(-n_y[2], n[2]);
}

TEST_F(PolyhedralMeshBase, PolyhedralMeshStatusTest) {

    generatePolyhedralMesh(mesh_);

    // Request status
    StatusAttrib status(mesh_);

    // Select a few faces
    status[FaceHandle(1)].set_tagged(true);
    status[FaceHandle(4)].set_tagged(true);

    status[HalfFaceHandle(21)].set_deleted(true);
    status[HalfFaceHandle(0)].set_deleted(true);

    status[VertexHandle(3)].set_selected(true);
    status[VertexHandle(8)].set_selected(true);

    EXPECT_TRUE(status[FaceHandle(1)].tagged());
    EXPECT_TRUE(status[FaceHandle(4)].tagged());
    EXPECT_FALSE(status[FaceHandle(7)].tagged());
    EXPECT_FALSE(status[FaceHandle(2)].tagged());

    EXPECT_TRUE(status[HalfFaceHandle(21)].deleted());
    EXPECT_TRUE(status[HalfFaceHandle(0)].deleted());
    EXPECT_FALSE(status[HalfFaceHandle(13)].deleted());
    EXPECT_FALSE(status[HalfFaceHandle(20)].deleted());

    EXPECT_TRUE(status[VertexHandle(3)].selected());
    EXPECT_TRUE(status[VertexHandle(8)].selected());
    EXPECT_FALSE(status[VertexHandle(1)].selected());
    EXPECT_FALSE(status[VertexHandle(9)].selected());
}

TEST_F(PolyhedralMeshBase, PolyhedralMeshColorTest) {

    generatePolyhedralMesh(mesh_);

    typedef OpenVolumeMesh::Geometry::Vec4f Vec4f;

    // Request colors
    ColorAttrib<Vec4f> colors(mesh_);

    colors[VertexHandle(7)] = Vec4f(1.0f, 1.0f, 0.0f, 1.0f);
    colors[EdgeHandle(6)] = Vec4f(1.0f, 1.0f, 0.0f, 1.0f);
    colors[HalfEdgeHandle(5)] = Vec4f(1.0f, 1.0f, 0.0f, 1.0f);
    colors[FaceHandle(4)] = Vec4f(1.0f, 1.0f, 0.0f, 1.0f);
    colors[HalfFaceHandle(3)] = Vec4f(1.0f, 1.0f, 0.0f, 1.0f);
    colors[CellHandle(1)] = Vec4f(1.0f, 1.0f, 0.0f, 1.0f);

    EXPECT_FLOAT_EQ(1.0f, colors[VertexHandle(7)][0]);
    EXPECT_FLOAT_EQ(1.0f, colors[VertexHandle(7)][1]);
    EXPECT_FLOAT_EQ(0.0f, colors[VertexHandle(7)][2]);
    EXPECT_FLOAT_EQ(1.0f, colors[VertexHandle(7)][3]);
    EXPECT_FLOAT_EQ(1.0f, colors[EdgeHandle(6)][0]);
    EXPECT_FLOAT_EQ(1.0f, colors[EdgeHandle(6)][1]);
    EXPECT_FLOAT_EQ(0.0f, colors[EdgeHandle(6)][2]);
    EXPECT_FLOAT_EQ(1.0f, colors[EdgeHandle(6)][3]);
    EXPECT_FLOAT_EQ(1.0f, colors[HalfEdgeHandle(5)][0]);
    EXPECT_FLOAT_EQ(1.0f, colors[HalfEdgeHandle(5)][1]);
    EXPECT_FLOAT_EQ(0.0f, colors[HalfEdgeHandle(5)][2]);
    EXPECT_FLOAT_EQ(1.0f, colors[HalfEdgeHandle(5)][3]);
    EXPECT_FLOAT_EQ(1.0f, colors[FaceHandle(4)][0]);
    EXPECT_FLOAT_EQ(1.0f, colors[FaceHandle(4)][1]);
    EXPECT_FLOAT_EQ(0.0f, colors[FaceHandle(4)][2]);
    EXPECT_FLOAT_EQ(1.0f, colors[FaceHandle(4)][3]);
    EXPECT_FLOAT_EQ(1.0f, colors[HalfFaceHandle(3)][0]);
    EXPECT_FLOAT_EQ(1.0f, colors[HalfFaceHandle(3)][1]);
    EXPECT_FLOAT_EQ(0.0f, colors[HalfFaceHandle(3)][2]);
    EXPECT_FLOAT_EQ(1.0f, colors[HalfFaceHandle(3)][3]);
    EXPECT_FLOAT_EQ(1.0f, colors[CellHandle(1)][0]);
    EXPECT_FLOAT_EQ(1.0f, colors[CellHandle(1)][1]);
    EXPECT_FLOAT_EQ(0.0f, colors[CellHandle(1)][2]);
    EXPECT_FLOAT_EQ(1.0f, colors[CellHandle(1)][3]);
}

TEST_F(PolyhedralMeshBase, PolyhedralMeshProperties) {

    generatePolyhedralMesh(mesh_);

    VertexPropertyT<Vec3d> vp = mesh_.request_vertex_property<Vec3d>("VProp");

    EXPECT_TRUE(mesh_.vertex_property_exists<Vec3d>("VProp"));

    for(VertexIter v_it = mesh_.v_iter(); v_it.valid(); ++v_it) {
        vp[v_it->idx()] = Vec3d(1.0, 0.0, 0.0);
    }

    for(VertexIter v_it = mesh_.v_iter(); v_it.valid(); ++v_it) {
        Vec3d t;
        t = vp[v_it->idx()];
        EXPECT_DOUBLE_EQ(1.0, t[0]);
        EXPECT_DOUBLE_EQ(0.0, t[1]);
        EXPECT_DOUBLE_EQ(0.0, t[2]);
    }

    VertexHandle vh = mesh_.add_vertex(Vec3d(3.0,3.0,3.0));
    vp[vh.idx()] = Vec3d(0.0);
    Vec3d p = vp[vh.idx()];
    EXPECT_DOUBLE_EQ(0.0, p[0]);
    EXPECT_DOUBLE_EQ(0.0, p[1]);
    EXPECT_DOUBLE_EQ(0.0, p[2]);

    EdgePropertyT<unsigned int> ep = mesh_.request_edge_property<unsigned int>("EProp");

    EXPECT_TRUE(mesh_.edge_property_exists<unsigned int>("EProp"));

    unsigned int i = 0;
    for(EdgeIter e_it = mesh_.e_iter(); e_it.valid(); ++e_it) {
        ep[e_it->idx()] = i++;
    }

    i = 0;
    for(EdgeIter e_it = mesh_.e_iter(); e_it.valid(); ++e_it) {
        EXPECT_EQ(i++, ep[e_it->idx()]);
    }

    HalfFacePropertyT<bool> hfp = mesh_.request_halfface_property<bool>("HFProp");

    EXPECT_TRUE(mesh_.halfface_property_exists<bool>("HFProp"));

    bool b = false;
    for(HalfFaceIter hf_it = mesh_.hf_iter(); hf_it.valid(); ++hf_it) {
        hfp[hf_it->idx()] = b;
        b = !b;
    }

    b = false;
    for(HalfFaceIter hf_it = mesh_.hf_iter(); hf_it.valid(); ++hf_it) {
        EXPECT_EQ(b, hfp[hf_it->idx()]);
        b = !b;
    }

    // Request halfface properties
    CellPropertyT<std::string> cp = mesh_.request_cell_property<std::string>("CProp");

    EXPECT_TRUE(mesh_.cell_property_exists<std::string>("CProp"));

    for(CellIter c_it = mesh_.c_iter(); c_it.valid(); ++c_it) {
        cp[c_it->idx()] = std::string("MyTestString");
    }

    for(CellIter c_it = mesh_.c_iter(); c_it.valid(); ++c_it) {
        EXPECT_EQ(std::string("MyTestString"), cp[c_it->idx()]);
    }

    EXPECT_FALSE(mesh_.halfedge_property_exists<unsigned char>("HEProp"));
    EXPECT_FALSE(mesh_.vertex_property_exists<size_t>(""));
}

TEST_F(PolyhedralMeshBase, STLCompliance) {

    generatePolyhedralMesh(mesh_);

    Print p;
    p.mute(true);
    //std::cerr << "Vertices:" << std::endl;
    std::for_each(mesh_.vertices_begin(), mesh_.vertices_end(), p);
    //std::cerr << "Edges:" << std::endl;
    std::for_each(mesh_.edges_begin(), mesh_.edges_end(), p);
    //std::cerr << "HalfEdges:" << std::endl;
    std::for_each(mesh_.halfedges_begin(), mesh_.halfedges_end(), p);
    //std::cerr << "Faces:" << std::endl;
    std::for_each(mesh_.faces_begin(), mesh_.faces_end(), p);
    //std::cerr << "HalfFaces:" << std::endl;
    std::for_each(mesh_.halffaces_begin(), mesh_.halffaces_end(), p);
    //std::cerr << "Cells:" << std::endl;
    std::for_each(mesh_.cells_begin(), mesh_.cells_end(), p);
}

/*
 * Hexahedral mesh tests
 */

TEST_F(HexahedralMeshBase, SimpleHexMeshNavigation) {

    generateHexahedralMesh(mesh_);

    EXPECT_EQ(12u, mesh_.n_vertices());
    EXPECT_EQ(20u, mesh_.n_edges());
    EXPECT_EQ(11u, mesh_.n_faces());
    EXPECT_EQ(2u,  mesh_.n_cells());

    EXPECT_EQ(1,  mesh_.xfront_halfface(CellHandle(0)));
    EXPECT_EQ(2, mesh_.xback_halfface(CellHandle(0)));
    EXPECT_EQ(5,  mesh_.yfront_halfface(CellHandle(0)));
    EXPECT_EQ(6,  mesh_.yback_halfface(CellHandle(0)));
    EXPECT_EQ(8,  mesh_.zfront_halfface(CellHandle(0)));
    EXPECT_EQ(11,  mesh_.zback_halfface(CellHandle(0)));

    EXPECT_EQ(12, mesh_.opposite_halfface_handle_in_cell(
            HalfFaceHandle(3), CellHandle(1)));

    EXPECT_EQ(HalfFaceHandle(20), mesh_.adjacent_halfface_on_sheet(
            HalfFaceHandle(9), HalfEdgeHandle(12)));
    EXPECT_EQ(HalfFaceHandle(21), mesh_.adjacent_halfface_on_sheet(
            HalfFaceHandle(8), HalfEdgeHandle(12)));

    HexahedralMesh::CellSheetCellIter csc_it = mesh_.csc_iter(0, HexahedralMesh::YF);
    EXPECT_EQ(CellHandle(1), *csc_it);

    HexahedralMesh::HalfFaceSheetHalfFaceIter hfshf_it = mesh_.hfshf_iter(5);
    EXPECT_EQ(HalfFaceHandle(15), *hfshf_it);
    hfshf_it = mesh_.hfshf_iter(6);
    EXPECT_EQ(HalfFaceHandle(16), *hfshf_it);
}

TEST_F(HexahedralMeshBase, BottomUpIncidenceUpdate1) {

    generateHexahedralMesh(mesh_);

    EXPECT_EQ(12u, mesh_.n_vertices());
    EXPECT_EQ(20u, mesh_.n_edges());
    EXPECT_EQ(11u, mesh_.n_faces());
    EXPECT_EQ(2u, mesh_.n_cells());

    mesh_.delete_vertex(VertexHandle(0));

    EXPECT_EQ(11u, mesh_.n_vertices());
    EXPECT_EQ(17u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_faces());
    EXPECT_EQ(1u, mesh_.n_cells());

    HexVertexIter hv_it = mesh_.hv_iter(CellHandle(0));

    EXPECT_EQ(3, hv_it->idx()); ++hv_it;
    EXPECT_EQ(4, hv_it->idx()); ++hv_it;
    EXPECT_EQ(5, hv_it->idx()); ++hv_it;
    EXPECT_EQ(6, hv_it->idx()); ++hv_it;
    EXPECT_EQ(7, hv_it->idx()); ++hv_it;
    EXPECT_EQ(10, hv_it->idx()); ++hv_it;
    EXPECT_EQ(9, hv_it->idx()); ++hv_it;
    EXPECT_EQ(8, hv_it->idx());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTest1) {

    generateHexahedralMesh(mesh_);

    EXPECT_EQ(12u, mesh_.n_vertices());
    EXPECT_EQ(20u, mesh_.n_edges());
    EXPECT_EQ(11u, mesh_.n_faces());
    EXPECT_EQ(2u, mesh_.n_cells());

    StatusAttrib status(mesh_);

    HexVertexIter hv_it = mesh_.hv_iter(CellHandle(1));

    EXPECT_EQ(4, hv_it->idx()); ++hv_it;
    EXPECT_EQ(5, hv_it->idx()); ++hv_it;
    EXPECT_EQ(6, hv_it->idx()); ++hv_it;
    EXPECT_EQ(7, hv_it->idx()); ++hv_it;
    EXPECT_EQ(8, hv_it->idx()); ++hv_it;
    EXPECT_EQ(11, hv_it->idx()); ++hv_it;
    EXPECT_EQ(10, hv_it->idx()); ++hv_it;
    EXPECT_EQ(9, hv_it->idx());

    status[VertexHandle(0)].set_deleted(true);

    status.garbage_collection(false);

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(11u, mesh_.n_vertices());
    EXPECT_EQ(17u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_faces());

    hv_it = mesh_.hv_iter(CellHandle(0));

    EXPECT_EQ(3, hv_it->idx()); ++hv_it;
    EXPECT_EQ(4, hv_it->idx()); ++hv_it;
    EXPECT_EQ(5, hv_it->idx()); ++hv_it;
    EXPECT_EQ(6, hv_it->idx()); ++hv_it;
    EXPECT_EQ(7, hv_it->idx()); ++hv_it;
    EXPECT_EQ(10, hv_it->idx()); ++hv_it;
    EXPECT_EQ(9, hv_it->idx()); ++hv_it;
    EXPECT_EQ(8, hv_it->idx());

    status.garbage_collection(true);

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(8u, mesh_.n_vertices());
    EXPECT_EQ(12u, mesh_.n_edges());
    EXPECT_EQ(6u, mesh_.n_faces());

    hv_it = mesh_.hv_iter(CellHandle(0));

    EXPECT_EQ(0, hv_it->idx()); ++hv_it;
    EXPECT_EQ(1, hv_it->idx()); ++hv_it;
    EXPECT_EQ(2, hv_it->idx()); ++hv_it;
    EXPECT_EQ(3, hv_it->idx()); ++hv_it;
    EXPECT_EQ(4, hv_it->idx()); ++hv_it;
    EXPECT_EQ(7, hv_it->idx()); ++hv_it;
    EXPECT_EQ(6, hv_it->idx()); ++hv_it;
    EXPECT_EQ(5, hv_it->idx());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTest2) {

    generateHexahedralMesh(mesh_);

    StatusAttrib status(mesh_);

    status[VertexHandle(0)].set_deleted(true);

    status.garbage_collection(false);

    status.garbage_collection(false);

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(11u, mesh_.n_vertices());
    EXPECT_EQ(17u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_faces());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTest3) {

    generateHexahedralMesh(mesh_);

    StatusAttrib status(mesh_);

    status[EdgeHandle(0)].set_deleted(true);

    status.garbage_collection(false);

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(12u, mesh_.n_vertices());
    EXPECT_EQ(19u, mesh_.n_edges());
    EXPECT_EQ(9u, mesh_.n_faces());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTest4) {

    generateHexahedralMesh(mesh_);

    StatusAttrib status(mesh_);

    status[EdgeHandle(5)].set_deleted(true);

    status.garbage_collection(false);

    EXPECT_EQ(0u, mesh_.n_cells());
    EXPECT_EQ(12u, mesh_.n_vertices());
    EXPECT_EQ(19u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_faces());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTest5) {

    generateHexahedralMesh(mesh_);

    StatusAttrib status(mesh_);

    status[EdgeHandle(0)].set_deleted(true);
    status[EdgeHandle(1)].set_deleted(true);
    status[EdgeHandle(2)].set_deleted(true);

    status.garbage_collection(true);

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(8u, mesh_.n_vertices());
    EXPECT_EQ(12u, mesh_.n_edges());
    EXPECT_EQ(6u, mesh_.n_faces());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTestManifoldness1) {

    generateHexahedralMesh(mesh_);

    StatusAttrib status(mesh_);

    status[VertexHandle(0)].set_deleted(true);

    status.garbage_collection(true);

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(6u, mesh_.n_faces());
    EXPECT_EQ(12u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_vertices());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTestManifoldness2) {

    generateHexahedralMesh(mesh_);

    StatusAttrib status(mesh_);

    status[EdgeHandle(0)].set_deleted(true);

    status.garbage_collection(true);

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(6u, mesh_.n_faces());
    EXPECT_EQ(12u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_vertices());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTestManifoldness3) {

    generateHexahedralMesh(mesh_);

    StatusAttrib status(mesh_);

    status[FaceHandle(0)].set_deleted(true);

    status.garbage_collection(true);

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(6u, mesh_.n_faces());
    EXPECT_EQ(12u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_vertices());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTestManifoldness4) {

    generateHexahedralMesh(mesh_);

    StatusAttrib status(mesh_);

    status[CellHandle(0)].set_deleted(true);

    status.garbage_collection(true);

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(6u, mesh_.n_faces());
    EXPECT_EQ(12u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_vertices());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTestManifoldness5) {

    generateHexahedralMesh(mesh_);

    StatusAttrib status(mesh_);

    status[EdgeHandle(5)].set_deleted(true);

    status.garbage_collection(false);

    EXPECT_EQ(0u, mesh_.n_cells());
    EXPECT_EQ(8u, mesh_.n_faces());
    EXPECT_EQ(19u, mesh_.n_edges());
    EXPECT_EQ(12u, mesh_.n_vertices());

    status.garbage_collection(true);

    EXPECT_EQ(0u, mesh_.n_cells());
    EXPECT_EQ(0u, mesh_.n_faces());
    EXPECT_EQ(0u, mesh_.n_edges());
    EXPECT_EQ(0u, mesh_.n_vertices());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTestProps1) {

    generateHexahedralMesh(mesh_);

    StatusAttrib status(mesh_);

    FacePropertyT<int> fprop = mesh_.request_face_property<int>("FProp");

    FaceHandle fh0(0);
    FaceHandle fh1(1);
    FaceHandle fh2(2);
    FaceHandle fh3(3);
    FaceHandle fh4(4);
    FaceHandle fh5(5);
    FaceHandle fh6(6);
    FaceHandle fh7(7);
    FaceHandle fh8(8);
    FaceHandle fh9(9);
    FaceHandle fh10(10);

    fprop[fh0] = 11;
    fprop[fh1] = 10;
    fprop[fh2] = 9;
    fprop[fh3] = 8;
    fprop[fh4] = 7;
    fprop[fh5] = 6;
    fprop[fh6] = 5;
    fprop[fh7] = 4;
    fprop[fh8] = 3;
    fprop[fh9] = 2;
    fprop[fh10] = 1;

    status[VertexHandle(0)].set_deleted(true);

    status.garbage_collection(false);

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(11u, mesh_.n_vertices());
    EXPECT_EQ(17u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_faces());

    std::set<int> fprops_i;
    for(FaceIter f_it = mesh_.f_iter(); f_it.valid(); ++f_it) {
        fprops_i.insert(fprop[f_it->idx()]);
    }

    EXPECT_EQ(0, fprops_i.count(11));
    EXPECT_EQ(1, fprops_i.count(10));
    EXPECT_EQ(1, fprops_i.count(9));
    EXPECT_EQ(0, fprops_i.count(8));
    EXPECT_EQ(0, fprops_i.count(7));
    EXPECT_EQ(1, fprops_i.count(6));
    EXPECT_EQ(1, fprops_i.count(5));
    EXPECT_EQ(1, fprops_i.count(4));
    EXPECT_EQ(1, fprops_i.count(3));
    EXPECT_EQ(1, fprops_i.count(2));
    EXPECT_EQ(1, fprops_i.count(1));
}

TEST_F(HexahedralMeshBase, GarbageCollectionTestProps2) {

    generateHexahedralMesh(mesh_);

    StatusAttrib status(mesh_);

    FacePropertyT<int> fprop = mesh_.request_face_property<int>("FProp");

    FaceHandle fh0(0);
    FaceHandle fh1(1);
    FaceHandle fh2(2);
    FaceHandle fh3(3);
    FaceHandle fh4(4);
    FaceHandle fh5(5);
    FaceHandle fh6(6);
    FaceHandle fh7(7);
    FaceHandle fh8(8);
    FaceHandle fh9(9);
    FaceHandle fh10(10);

    fprop[fh0] = 11;
    fprop[fh1] = 10;
    fprop[fh2] = 9;
    fprop[fh3] = 8;
    fprop[fh4] = 7;
    fprop[fh5] = 6;
    fprop[fh6] = 5;
    fprop[fh7] = 4;
    fprop[fh8] = 3;
    fprop[fh9] = 2;
    fprop[fh10] = 1;

    status[FaceHandle(0)].set_deleted(true);

    status.garbage_collection(false);

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(12u, mesh_.n_vertices());
    EXPECT_EQ(20u, mesh_.n_edges());
    EXPECT_EQ(10u, mesh_.n_faces());

    std::set<int> fprops_i;
    for(FaceIter f_it = mesh_.f_iter(); f_it.valid(); ++f_it) {
        fprops_i.insert(fprop[f_it->idx()]);
    }

    EXPECT_EQ(0, fprops_i.count(11));
    EXPECT_EQ(1, fprops_i.count(10));
    EXPECT_EQ(1, fprops_i.count(9));
    EXPECT_EQ(1, fprops_i.count(8));
    EXPECT_EQ(1, fprops_i.count(7));
    EXPECT_EQ(1, fprops_i.count(6));
    EXPECT_EQ(1, fprops_i.count(5));
    EXPECT_EQ(1, fprops_i.count(4));
    EXPECT_EQ(1, fprops_i.count(3));
    EXPECT_EQ(1, fprops_i.count(2));
    EXPECT_EQ(1, fprops_i.count(1));
}

TEST_F(HexahedralMeshBase, HalfEdgeFetchFunction1) {

    generateHexahedralMesh(mesh_);

    VertexHandle v0(0);
    VertexHandle v1(1);

    VertexHandle v2(2);
    VertexHandle v3(3);

    VertexHandle v5(5);
    VertexHandle v6(6);
    VertexHandle v7(7);

    HalfEdgeHandle he0 = mesh_.halfedge(v0, v1);
    HalfEdgeHandle he5 = mesh_.halfedge(v3, v2);
    HalfEdgeHandle he10 = mesh_.halfedge(v5, v6);
    HalfEdgeHandle heInv = mesh_.halfedge(v5, v7);

    EXPECT_EQ(HalfEdgeHandle(0), he0);
    EXPECT_EQ(HalfEdgeHandle(5), he5);
    EXPECT_EQ(HalfEdgeHandle(10), he10);
    EXPECT_EQ(HexahedralMesh::InvalidHalfEdgeHandle, heInv);
}

TEST_F(HexahedralMeshBase, HalfFaceFetchFunction1) {

    generateHexahedralMesh(mesh_);

    HalfEdgeHandle he0(0);
    HalfEdgeHandle he2(2);
    HalfEdgeHandle he4(4);

    std::vector<HalfEdgeHandle> hes;
    hes.push_back(he0); hes.push_back(he2);

    HalfFaceHandle hf0_0 = mesh_.halfface(hes);
    hes.clear();
    hes.push_back(he0); hes.push_back(he4);
    HalfFaceHandle hf0_1 = mesh_.halfface(hes);

    HalfEdgeHandle he16(16);
    HalfEdgeHandle he18(18);

    hes.clear();
    hes.push_back(he16); hes.push_back(he18);
    HalfFaceHandle hf4_0 = mesh_.halfface(hes);

    hes.clear();
    hes.push_back(he0); hes.push_back(he18);
    HalfFaceHandle hfInv = mesh_.halfface(hes);

    HalfEdgeHandle he17(17);
    HalfEdgeHandle he19(19);

    hes.clear();
    hes.push_back(he17); hes.push_back(he19);
    HalfFaceHandle hf5_0 = mesh_.halfface(hes);

    EXPECT_EQ(HalfFaceHandle(0), hf0_0);
    EXPECT_EQ(HalfFaceHandle(0), hf0_1);
    EXPECT_EQ(HalfFaceHandle(4), hf4_0);
    EXPECT_EQ(HalfFaceHandle(5), hf5_0);
    EXPECT_EQ(HexahedralMesh::InvalidHalfFaceHandle, hfInv);
}

TEST_F(HexahedralMeshBase, HalfFaceFetchFunction2) {

    generateHexahedralMesh(mesh_);

    VertexHandle v0(0);
    VertexHandle v1(1);
    VertexHandle v2(2);
    VertexHandle v3(3);
    VertexHandle v4(4);
    VertexHandle v5(5);
    VertexHandle v6(6);
    VertexHandle v7(7);

    std::vector<VertexHandle> vs;
    vs.push_back(v0); vs.push_back(v1); vs.push_back(v2);
    HalfFaceHandle hf0 = mesh_.halfface(vs); vs.clear();

    vs.push_back(v2); vs.push_back(v1); vs.push_back(v0);
    HalfFaceHandle hf1 = mesh_.halfface(vs); vs.clear();

    vs.push_back(v2); vs.push_back(v1); vs.push_back(v5);
    HalfFaceHandle hf4 = mesh_.halfface(vs); vs.clear();

    vs.push_back(v6); vs.push_back(v5); vs.push_back(v4);
    HalfFaceHandle hf3 = mesh_.halfface(vs); vs.clear();

    vs.push_back(v4); vs.push_back(v5); vs.push_back(v6);
    HalfFaceHandle hf2 = mesh_.halfface(vs); vs.clear();

    vs.push_back(v0); vs.push_back(v1); vs.push_back(v4);
    HalfFaceHandle hfInv0 = mesh_.halfface(vs); vs.clear();

    vs.push_back(v0); vs.push_back(v1); vs.push_back(v6);
    HalfFaceHandle hfInv1 = mesh_.halfface(vs); vs.clear();

    EXPECT_EQ(HalfFaceHandle(0), hf0);
    EXPECT_EQ(HalfFaceHandle(1), hf1);
    EXPECT_EQ(HalfFaceHandle(4), hf4);
    EXPECT_EQ(HalfFaceHandle(3), hf3);
    EXPECT_EQ(HalfFaceHandle(2), hf2);
    EXPECT_EQ(HexahedralMesh::InvalidHalfFaceHandle, hfInv0);
    EXPECT_EQ(HexahedralMesh::InvalidHalfFaceHandle, hfInv1);
}

TEST_F(HexahedralMeshBase, AddCellViaVerticesFunction1) {

    generateHexahedralMesh(mesh_);

    StatusAttrib status(mesh_);

    status[FaceHandle(0)].set_deleted(true);

    status.garbage_collection(false);

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(10u, mesh_.n_faces());
    EXPECT_EQ(20u, mesh_.n_edges());
    EXPECT_EQ(12u, mesh_.n_vertices());

    std::vector<VertexHandle> vs;
    vs.push_back(VertexHandle(0));
    vs.push_back(VertexHandle(1));
    vs.push_back(VertexHandle(2));
    vs.push_back(VertexHandle(3));
    vs.push_back(VertexHandle(4));
    vs.push_back(VertexHandle(7));
    vs.push_back(VertexHandle(6));
    vs.push_back(VertexHandle(5));

    CellHandle ch = mesh_.add_cell(vs);

    EXPECT_NE(HexahedralMesh::InvalidCellHandle, ch);

    EXPECT_EQ(2u, mesh_.n_cells());
    EXPECT_EQ(11u, mesh_.n_faces());
    EXPECT_EQ(20u, mesh_.n_edges());
    EXPECT_EQ(12u, mesh_.n_vertices());
}

TEST_F(HexahedralMeshBase, AddCellViaVerticesFunction2) {

    generateHexahedralMesh(mesh_);

    StatusAttrib status(mesh_);

    status[FaceHandle(0)].set_deleted(true);

    status.garbage_collection(true);

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(6u, mesh_.n_faces());
    EXPECT_EQ(12u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_vertices());

    VertexHandle v0 = mesh_.add_vertex(Vec3d(0.0, 0.0, 0.0));
    VertexHandle v1 = mesh_.add_vertex(Vec3d(1.0, 0.0, 0.0));
    VertexHandle v2 = mesh_.add_vertex(Vec3d(1.0, 1.0, 0.0));
    VertexHandle v3 = mesh_.add_vertex(Vec3d(0.0, 1.0, 0.0));

    std::vector<VertexHandle> vs;
    vs.push_back(v0);
    vs.push_back(v1);
    vs.push_back(v2);
    vs.push_back(v3);
    vs.push_back(VertexHandle(0));
    vs.push_back(VertexHandle(3));
    vs.push_back(VertexHandle(2));
    vs.push_back(VertexHandle(1));

    CellHandle ch = mesh_.add_cell(vs);

    EXPECT_NE(HexahedralMesh::InvalidCellHandle, ch);

    EXPECT_EQ(2u, mesh_.n_cells());
    EXPECT_EQ(11u, mesh_.n_faces());
    EXPECT_EQ(20u, mesh_.n_edges());
    EXPECT_EQ(12u, mesh_.n_vertices());
}

#endif // INCLUDE GUARD
