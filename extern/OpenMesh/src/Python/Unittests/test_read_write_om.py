import unittest
import openmesh
import os

class ReadWriteOM(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()

    def test_load_simple_om_force_vertex_colors_although_not_available(self):
        self.mesh.request_vertex_colors()
        
        file_name = "cube-minimal.om"
        
        options = openmesh.Options()
        options += openmesh.Options.VertexColor
        
        ok = openmesh.read_mesh(self.mesh, file_name, options)
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        self.assertEqual(self.mesh.n_halfedges(), 36)

        self.assertFalse(options.vertex_has_normal())
        self.assertFalse(options.vertex_has_texcoord())
        self.assertFalse(options.vertex_has_color())

    def test_load_simple_om_with_texcoords(self):
        self.mesh.request_vertex_texcoords2D()

        options = openmesh.Options()
        options += openmesh.Options.VertexTexCoord
        
        ok = openmesh.read_mesh(self.mesh, "cube-minimal-texCoords.om", options)
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(0))[0], 10.0)
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(0))[1], 10.0)
        
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(2))[0], 6.0)
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(2))[1], 6.0)
        
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(4))[0], 9.0)
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(4))[1], 9.0)
        
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(7))[0], 12.0)
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(7))[1], 12.0)
        
        self.assertFalse(options.vertex_has_normal())
        self.assertTrue(options.vertex_has_texcoord())
        self.assertFalse(options.vertex_has_color())
        
        self.mesh.release_vertex_texcoords2D()

    def test_load_simple_om_with_vertex_colors(self):
        self.mesh.request_vertex_colors()
        
        options = openmesh.Options()
        options += openmesh.Options.VertexColor
        
        ok = openmesh.read_mesh(self.mesh, "cube-minimal-vertexColors.om", options)
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[0], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[2], 0.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[0], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[2], 0.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[2], 1.0)
        
        self.assertFalse(options.vertex_has_normal())
        self.assertFalse(options.vertex_has_texcoord())
        self.assertTrue(options.vertex_has_color())
        
        self.mesh.release_vertex_colors()

    def test_write_triangle(self):
        filename = "triangle-minimal.om";
        
        # Generate data
        v1 = self.mesh.add_vertex(openmesh.Vec3d(1.0, 0.0, 0.0))
        v2 = self.mesh.add_vertex(openmesh.Vec3d(0.0, 1.0, 0.0))
        v3 = self.mesh.add_vertex(openmesh.Vec3d(0.0, 0.0, 1.0))
        self.mesh.add_face(v1, v2, v3)
        
        # Save
        ok = openmesh.write_mesh(self.mesh, filename)
        self.assertTrue(ok)
        
        # Reset
        self.mesh.clear()
        
        # Load
        ok = openmesh.read_mesh(self.mesh, filename)
        self.assertTrue(ok)
        
        # Compare
        self.assertEqual(self.mesh.n_vertices(), 3)
        self.assertEqual(self.mesh.n_edges(), 3)
        self.assertEqual(self.mesh.n_faces(), 1)
        
        self.assertEqual(self.mesh.point(v1), openmesh.Vec3d(1.0, 0.0, 0.0))
        self.assertEqual(self.mesh.point(v2), openmesh.Vec3d(0.0, 1.0, 0.0))
        self.assertEqual(self.mesh.point(v3), openmesh.Vec3d(0.0, 0.0, 1.0))
        
        # Cleanup
        os.remove(filename)

    def test_write_triangle_vertex_integer_color(self):
        self.mesh.request_vertex_colors()
            
        options = openmesh.Options()
        options += openmesh.Options.VertexColor
        options += openmesh.Options.ColorFloat

        filename = "triangle-minimal-ColorsPerVertex.om"
            
        # Generate data
        v1 = self.mesh.add_vertex(openmesh.Vec3d(1.0, 0.0, 0.0))
        v2 = self.mesh.add_vertex(openmesh.Vec3d(0.0, 1.0, 0.0))
        v3 = self.mesh.add_vertex(openmesh.Vec3d(0.0, 0.0, 1.0))
        self.mesh.add_face(v1, v2, v3)
        
        c1 = openmesh.Vec4f(0.00, 0.00, 0.50, 1.00)
        c2 = openmesh.Vec4f(0.25, 0.00, 0.00, 1.00)
        c3 = openmesh.Vec4f(0.00, 0.75, 0.00, 1.00)
    
        self.mesh.set_color(v1, c1)
        self.mesh.set_color(v2, c2)
        self.mesh.set_color(v3, c3)
            
        # Save
        ok = openmesh.write_mesh(self.mesh, filename, options)
        self.assertTrue(ok)
            
        self.mesh.release_vertex_colors()
            
        # Load
        cmpMesh = openmesh.TriMesh()
        cmpMesh.request_vertex_colors()
        ok = openmesh.read_mesh(cmpMesh, filename, options)
        self.assertTrue(ok)
            
        self.assertTrue(cmpMesh.has_vertex_colors())
            
        # Compare
        self.assertEqual(self.mesh.n_vertices(), 3)
        self.assertEqual(self.mesh.n_edges(), 3)
        self.assertEqual(self.mesh.n_faces(), 1)
        
        self.assertEqual(cmpMesh.point(v1), openmesh.Vec3d(1.0, 0.0, 0.0))
        self.assertEqual(cmpMesh.point(v2), openmesh.Vec3d(0.0, 1.0, 0.0))
        self.assertEqual(cmpMesh.point(v3), openmesh.Vec3d(0.0, 0.0, 1.0))
        
        self.assertAlmostEqual(cmpMesh.color(v1)[0], c1[0], 2)
        self.assertAlmostEqual(cmpMesh.color(v1)[1], c1[1], 2)
        self.assertAlmostEqual(cmpMesh.color(v1)[2], c1[2], 2)
        self.assertAlmostEqual(cmpMesh.color(v1)[3], c1[3], 2)
        
        self.assertAlmostEqual(cmpMesh.color(v2)[0], c2[0], 2)
        self.assertAlmostEqual(cmpMesh.color(v2)[1], c2[1], 2)
        self.assertAlmostEqual(cmpMesh.color(v2)[2], c2[2], 2)
        self.assertAlmostEqual(cmpMesh.color(v2)[3], c2[3], 2)
        
        self.assertAlmostEqual(cmpMesh.color(v3)[0], c3[0], 2)
        self.assertAlmostEqual(cmpMesh.color(v3)[1], c3[1], 2)
        self.assertAlmostEqual(cmpMesh.color(v3)[2], c3[2], 2)
        self.assertAlmostEqual(cmpMesh.color(v3)[3], c3[3], 2)
        
        # Clean up
        cmpMesh.release_vertex_colors()
        os.remove(filename)

    # TODO property tests


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ReadWriteOM)
    unittest.TextTestRunner(verbosity=2).run(suite)
