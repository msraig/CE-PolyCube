import unittest
import openmesh

class DeleteFaceTriangleMesh(unittest.TestCase):

    def test_delete_half_triangle_mesh_cube_no_edge_status(self):
        self.mesh = openmesh.TriMesh()
        self.vhandle = []
        
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1, -1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1, -1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1,  1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1,  1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1, -1, -1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1, -1, -1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1,  1, -1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1,  1, -1)))
        
        # Add six faces to form a cube
        face_vhandles = []
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)
        
        #=======================
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[7])
        face_vhandles.append(self.vhandle[6])
        face_vhandles.append(self.vhandle[5])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[7])
        face_vhandles.append(self.vhandle[5])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        #=======================
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[4])
        face_vhandles.append(self.vhandle[5])
        self.mesh.add_face(face_vhandles)
        
        #=======================
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[5])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[5])
        face_vhandles.append(self.vhandle[6])
        self.mesh.add_face(face_vhandles)
        
        #=======================
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[6])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[6])
        face_vhandles.append(self.vhandle[7])
        self.mesh.add_face(face_vhandles)
        
        #=======================
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[7])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[7])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        # Test setup:
        #
        #    3 ======== 2
        #   /          /|
        #  /          / |      z
        # 0 ======== 1  |      |
        # |          |  |      |   y
        # |  7       |  6      |  /
        # |          | /       | /
        # |          |/        |/
        # 4 ======== 5         -------> x
        
        # Check setup
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_halfedges(), 36)
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        # =====================================================
        # Now we delete half of the mesh
        # =====================================================
        self.mesh.request_face_status()
        self.mesh.request_vertex_status()
        self.mesh.request_halfedge_status()
        
        n_face_to_delete = self.mesh.n_faces() / 2
        
        # Check the variable
        self.assertEqual(n_face_to_delete, 6)
        
        for i in range(int(n_face_to_delete)):
            self.mesh.delete_face(self.mesh.face_handle(i))
            
        # =====================================================
        # Check config afterwards
        # =====================================================
        
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_halfedges(), 36)
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        # =====================================================
        # Cleanup and recheck
        # =====================================================
        
        self.mesh.garbage_collection()
        
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_halfedges(), 36)
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 6)
        
    def test_delete_half_triangle_mesh_cube_with_edge_status(self):
        self.mesh = openmesh.TriMesh()
        self.vhandle = []
        
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1, -1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1, -1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1,  1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1,  1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1, -1, -1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1, -1, -1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1,  1, -1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1,  1, -1)))
        
        # Add six faces to form a cube
        face_vhandles = []
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)
        
        #=======================
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[7])
        face_vhandles.append(self.vhandle[6])
        face_vhandles.append(self.vhandle[5])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[7])
        face_vhandles.append(self.vhandle[5])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        #=======================
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[4])
        face_vhandles.append(self.vhandle[5])
        self.mesh.add_face(face_vhandles)
        
        #=======================
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[5])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[5])
        face_vhandles.append(self.vhandle[6])
        self.mesh.add_face(face_vhandles)
        
        #=======================
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[6])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[6])
        face_vhandles.append(self.vhandle[7])
        self.mesh.add_face(face_vhandles)
        
        #=======================
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[7])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[7])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        # Test setup:
        #
        #    3 ======== 2
        #   /          /|
        #  /          / |      z
        # 0 ======== 1  |      |
        # |          |  |      |   y
        # |  7       |  6      |  /
        # |          | /       | /
        # |          |/        |/
        # 4 ======== 5         -------> x
        
        # Check setup
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_halfedges(), 36)
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        # =====================================================
        # Now we delete half of the mesh
        # =====================================================
        self.mesh.request_face_status()
        self.mesh.request_vertex_status()
        self.mesh.request_edge_status()
        self.mesh.request_halfedge_status()
        
        n_face_to_delete = self.mesh.n_faces() / 2
        
        # Check the variable
        self.assertEqual(n_face_to_delete, 6)
        
        for i in range(int(n_face_to_delete)):
            self.mesh.delete_face(self.mesh.face_handle(i))
            
        # =====================================================
        # Check config afterwards
        # =====================================================
        
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_halfedges(), 36)
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        # =====================================================
        # Cleanup and recheck
        # =====================================================
        
        self.mesh.garbage_collection()
        
        self.assertEqual(self.mesh.n_edges(), 13)
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 6)
        
    def test_deletete_half_poly_mesh_cube_without_edge_status(self):
        self.mesh = openmesh.PolyMesh()
        self.vhandle = []
        
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1, -1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1, -1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1,  1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1,  1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1, -1, -1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1, -1, -1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1,  1, -1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1,  1, -1)))
        
        # Add six faces to form a cube
        face_vhandles = []
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[7])
        face_vhandles.append(self.vhandle[6])
        face_vhandles.append(self.vhandle[5])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[4])
        face_vhandles.append(self.vhandle[5])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[5])
        face_vhandles.append(self.vhandle[6])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[6])
        face_vhandles.append(self.vhandle[7])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[7])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        # Test setup:
        #
        #    3 ======== 2
        #   /          /|
        #  /          / |      z
        # 0 ======== 1  |      |
        # |          |  |      |   y
        # |  7       |  6      |  /
        # |          | /       | /
        # |          |/        |/
        # 4 ======== 5         -------> x
        
        # Check setup
        self.assertEqual(self.mesh.n_edges(), 12)
        self.assertEqual(self.mesh.n_halfedges(), 24)
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 6)
        
        # =====================================================
        # Now we delete half of the mesh
        # =====================================================
        self.mesh.request_face_status()
        self.mesh.request_vertex_status()
        self.mesh.request_halfedge_status()
        
        n_face_to_delete = self.mesh.n_faces() / 2
        
        # Check the variable
        self.assertEqual(n_face_to_delete, 3)
        
        for i in range(int(n_face_to_delete)):
            self.mesh.delete_face(self.mesh.face_handle(i))
            
        # =====================================================
        # Check config afterwards
        # =====================================================
        
        self.assertEqual(self.mesh.n_edges(), 12)
        self.assertEqual(self.mesh.n_halfedges(), 24)
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 6)
        
        # =====================================================
        # Cleanup and recheck
        # =====================================================
        
        self.mesh.garbage_collection()
        
        self.assertEqual(self.mesh.n_edges(), 12)
        self.assertEqual(self.mesh.n_halfedges(), 24)
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 3)
        
    def test_deletete_half_poly_mesh_cube_with_edge_status(self):
        self.mesh = openmesh.PolyMesh()
        self.vhandle = []
        
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1, -1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1, -1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1,  1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1,  1,  1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1, -1, -1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1, -1, -1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1,  1, -1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1,  1, -1)))
        
        # Add six faces to form a cube
        face_vhandles = []
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[7])
        face_vhandles.append(self.vhandle[6])
        face_vhandles.append(self.vhandle[5])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[4])
        face_vhandles.append(self.vhandle[5])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[5])
        face_vhandles.append(self.vhandle[6])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[6])
        face_vhandles.append(self.vhandle[7])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[7])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        # Test setup:
        #
        #    3 ======== 2
        #   /          /|
        #  /          / |      z
        # 0 ======== 1  |      |
        # |          |  |      |   y
        # |  7       |  6      |  /
        # |          | /       | /
        # |          |/        |/
        # 4 ======== 5         -------> x
        
        # Check setup
        self.assertEqual(self.mesh.n_edges(), 12)
        self.assertEqual(self.mesh.n_halfedges(), 24)
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 6)
        
        # =====================================================
        # Now we delete half of the mesh
        # =====================================================
        self.mesh.request_face_status()
        self.mesh.request_vertex_status()
        self.mesh.request_edge_status()
        self.mesh.request_halfedge_status()
        
        n_face_to_delete = self.mesh.n_faces() / 2
        
        # Check the variable
        self.assertEqual(n_face_to_delete, 3)
        
        for i in range(int(n_face_to_delete)):
            self.mesh.delete_face(self.mesh.face_handle(i))
            
        # =====================================================
        # Check config afterwards
        # =====================================================
        
        self.assertEqual(self.mesh.n_edges(), 12)
        self.assertEqual(self.mesh.n_halfedges(), 24)
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 6)
        
        # =====================================================
        # Cleanup and recheck
        # =====================================================
        
        self.mesh.garbage_collection()
        
        self.assertEqual(self.mesh.n_edges(), 10)
        self.assertEqual(self.mesh.n_halfedges(), 20)
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 3)
        

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(DeleteFaceTriangleMesh)
    unittest.TextTestRunner(verbosity=2).run(suite)
