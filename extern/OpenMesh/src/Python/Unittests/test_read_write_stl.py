import unittest
import openmesh

class ReadWriteSTL(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()

    def test_load_simple_stl_file(self):
        ok = openmesh.read_mesh(self.mesh, "cube1.stl")
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 7526)
        self.assertEqual(self.mesh.n_edges(), 22572)
        self.assertEqual(self.mesh.n_faces(), 15048)
    
    def test_load_simple_stl_file_with_normals(self):
        self.mesh.request_face_normals()
        
        options = openmesh.Options()
        options += openmesh.Options.FaceNormal
        
        ok = openmesh.read_mesh(self.mesh, "cube1.stl", options)
        
        self.assertTrue(ok)
        
        self.assertAlmostEqual(self.mesh.normal(self.mesh.face_handle(0))[0], -0.038545)
        self.assertAlmostEqual(self.mesh.normal(self.mesh.face_handle(0))[1], -0.004330)
        self.assertAlmostEqual(self.mesh.normal(self.mesh.face_handle(0))[2],  0.999247)
        
        self.assertEqual(self.mesh.n_vertices(), 7526)
        self.assertEqual(self.mesh.n_edges(), 22572)
        self.assertEqual(self.mesh.n_faces(), 15048)
        
        self.mesh.release_face_normals()

    def test_load_simple_stl_binary_file(self):
        ok = openmesh.read_mesh(self.mesh, "cube1Binary.stl")
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 7526)
        self.assertEqual(self.mesh.n_edges(), 22572)
        self.assertEqual(self.mesh.n_faces(), 15048)

    def test_load_simple_stl_binary_file_with_normals(self):
        self.mesh.request_face_normals()
        
        options = openmesh.Options()
        options += openmesh.Options.FaceNormal
        options += openmesh.Options.Binary
        
        ok = openmesh.read_mesh(self.mesh, "cube1Binary.stl", options)
        
        self.assertTrue(ok)
        
        self.assertTrue(options.is_binary())
        self.assertTrue(options.face_has_normal())
        self.assertFalse(options.vertex_has_normal())
        
        self.assertAlmostEqual(self.mesh.normal(self.mesh.face_handle(0))[0], -0.038545, 5)
        self.assertAlmostEqual(self.mesh.normal(self.mesh.face_handle(0))[1], -0.004330, 5)
        self.assertAlmostEqual(self.mesh.normal(self.mesh.face_handle(0))[2],  0.999247, 5)
        
        self.assertEqual(self.mesh.n_vertices(), 7526)
        self.assertEqual(self.mesh.n_edges(), 22572)
        self.assertEqual(self.mesh.n_faces(), 15048)

        self.mesh.release_face_normals()


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ReadWriteSTL)
    unittest.TextTestRunner(verbosity=2).run(suite)
