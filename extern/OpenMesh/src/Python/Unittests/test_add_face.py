import unittest
import openmesh

class AddFace(unittest.TestCase):

    def test_add_triangles_to_trimesh(self):
        self.mesh = openmesh.TriMesh()
        self.vhandle = []
        
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 0, 0)))

        # Add two faces
        face_vhandles = []
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[0])
        
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[3])
        
        self.mesh.add_face(face_vhandles)

        # Test setup:
        #  1 === 2
        #  |   / |    
        #  |  /  |
        #  | /   |
        #  0 === 3
        
        # Check setup
        self.assertEqual(self.mesh.n_vertices(), 4)
        self.assertEqual(self.mesh.n_faces(), 2)

    def test_add_quad_to_trimesh(self):
        self.mesh = openmesh.TriMesh()
        self.vhandle = []
        
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 0, 0)))

        # Add two faces
        face_vhandles = []
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[3])
        
        self.mesh.add_face(face_vhandles)

        # Test setup:
        #  1 === 2
        #  |   / |    
        #  |  /  |
        #  | /   |
        #  0 === 3
        
        # Check setup
        self.assertEqual(self.mesh.n_vertices(), 4)
        self.assertEqual(self.mesh.n_faces(), 2)
        
    def test_create_triangle_mesh_cube(self):
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
        
    def test_create_strange_config(self):
        self.mesh = openmesh.TriMesh()
        self.vhandle = []
        
        #    2 x-----------x 1
        #       \         /
        #        \       /
        #         \     /
        #          \   /
        #           \ /
        #          0 x ---x 6
        #           /|\   |
        #          / | \  |
        #         /  |  \ |
        #        /   |   \|
        #       x----x    x
        #       3    4    5
        
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 1)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2, 2, 2)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(3, 3, 3)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(4, 4, 4)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(5, 5, 5)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(6, 6, 6)))
        
        self.mesh.add_face(self.vhandle[0], self.vhandle[1], self.vhandle[2])
        self.mesh.add_face(self.vhandle[0], self.vhandle[3], self.vhandle[4])
        self.mesh.add_face(self.vhandle[0], self.vhandle[5], self.vhandle[6])
        
        # non-manifold!
        invalid_fh = self.mesh.add_face(self.vhandle[3], self.vhandle[0], self.vhandle[4])

        # Check setup
        self.assertEqual(self.mesh.n_vertices(), 7)
        self.assertEqual(self.mesh.n_faces(), 3)
        self.assertEqual(invalid_fh, openmesh.TriMesh.InvalidFaceHandle)
        
    def test_add_quad_to_polymesh(self):
        self.mesh = openmesh.PolyMesh()
        self.vhandle = []
        
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 0, 0)))

        # Add one face
        face_vhandles = []
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[3])
        
        self.mesh.add_face(face_vhandles)

        # Test setup:
        #  1 === 2
        #  |     |    
        #  |     |
        #  |     |
        #  0 === 3
        
        # Check setup
        self.assertEqual(self.mesh.n_vertices(), 4)
        self.assertEqual(self.mesh.n_faces(), 1)
        
    def test_create_poly_mesh_cube(self):
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


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(AddFace)
    unittest.TextTestRunner(verbosity=2).run(suite)
