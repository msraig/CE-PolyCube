#include <gtest/gtest.h>
#include <Unittests/unittests_common.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>

#include <iostream>

namespace {

class OpenMeshPropertyManager : public OpenMeshBase {

    protected:

        // This function is called before each test is run
        virtual void SetUp() {
        }

        // This function is called after all tests are through
        virtual void TearDown() {

            // Do some final stuff with the member data here...
        }

    // Member already defined in OpenMeshBase
    //Mesh mesh_;  
};

/*
 * ====================================================================
 * Define tests below
 * ====================================================================
 */

/*
 * Collapsing a tetrahedron
 */
TEST_F(OpenMeshPropertyManager, set_range_bool) {

  mesh_.clear();

  // Add some vertices
  Mesh::VertexHandle vhandle[4];

  vhandle[0] = mesh_.add_vertex(Mesh::Point(0, 0, 0));
  vhandle[1] = mesh_.add_vertex(Mesh::Point(0, 1, 0));
  vhandle[2] = mesh_.add_vertex(Mesh::Point(1, 1, 0));
  vhandle[3] = mesh_.add_vertex(Mesh::Point(0, 0, 1));

  // Add two faces
  std::vector<Mesh::VertexHandle> face_vhandles;

  face_vhandles.push_back(vhandle[0]);
  face_vhandles.push_back(vhandle[1]);
  face_vhandles.push_back(vhandle[2]);
  mesh_.add_face(face_vhandles);

  face_vhandles.clear();

  face_vhandles.push_back(vhandle[0]);
  face_vhandles.push_back(vhandle[2]);
  face_vhandles.push_back(vhandle[3]);
  mesh_.add_face(face_vhandles);

  OpenMesh::PropertyManager<
      OpenMesh::VPropHandleT<bool>, Mesh> pm_v_bool(mesh_, "pm_v_bool");
  pm_v_bool.set_range(mesh_.vertices_begin(), mesh_.vertices_end(), false);
  for (int i = 0; i < 4; ++i)
      ASSERT_FALSE(pm_v_bool[vhandle[i]]);
  pm_v_bool.set_range(mesh_.vertices_begin(), mesh_.vertices_end(), true);
  for (int i = 0; i < 4; ++i)
      ASSERT_TRUE(pm_v_bool[vhandle[i]]);

  OpenMesh::PropertyManager<
      OpenMesh::EPropHandleT<bool>, Mesh> pm_e_bool(mesh_, "pm_e_bool");
  pm_e_bool.set_range(mesh_.edges_begin(), mesh_.edges_end(), false);
  for (Mesh::EdgeIter e_it = mesh_.edges_begin(), f_end = mesh_.edges_end();
          e_it != f_end; ++e_it)
      ASSERT_FALSE(pm_e_bool[*e_it]);
  pm_e_bool.set_range(mesh_.edges_begin(), mesh_.edges_end(), true);
  for (Mesh::EdgeIter e_it = mesh_.edges_begin(), f_end = mesh_.edges_end();
          e_it != f_end; ++e_it)
      ASSERT_TRUE(pm_e_bool[*e_it]);

  OpenMesh::PropertyManager<
      OpenMesh::FPropHandleT<bool>, Mesh> pm_f_bool(mesh_, "pm_f_bool");
  pm_f_bool.set_range(mesh_.faces_begin(), mesh_.faces_end(), false);
  for (Mesh::FaceIter f_it = mesh_.faces_begin(), f_end = mesh_.faces_end();
          f_it != f_end; ++f_it)
      ASSERT_FALSE(pm_f_bool[*f_it]);
  pm_f_bool.set_range(mesh_.faces_begin(), mesh_.faces_end(), true);
  for (Mesh::FaceIter f_it = mesh_.faces_begin(), f_end = mesh_.faces_end();
          f_it != f_end; ++f_it)
      ASSERT_TRUE(pm_f_bool[*f_it]);
}








}
