import unittest
import openmesh

class ReadWritePLY(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()

    def test_load_simple_point_ply_file_with_bad_encoding(self):
        ok = openmesh.read_mesh(self.mesh, "pointCloudBadEncoding.ply")
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 10)
        self.assertEqual(self.mesh.n_edges(), 0)
        self.assertEqual(self.mesh.n_faces(), 0)

    def test_load_simple_point_ply_file_with_good_encoding(self):
        ok = openmesh.read_mesh(self.mesh, "pointCloudGoodEncoding.ply")
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 10)
        self.assertEqual(self.mesh.n_edges(), 0)
        self.assertEqual(self.mesh.n_faces(), 0)

    def test_load_simple_ply(self):
        ok = openmesh.read_mesh(self.mesh, "cube-minimal.ply")
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)

    def test_load_simple_ply_force_vertex_colors_although_not_available(self):
        self.mesh.request_vertex_colors()
        
        file_name = "cube-minimal.ply"
        
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

    def test_load_simple_ply_with_vertex_colors(self):
        self.mesh.request_vertex_colors()
        
        file_name = "cube-minimal.ply"
        
        options = openmesh.Options()
        options += openmesh.Options.VertexColor
        
        ok = openmesh.read_mesh(self.mesh, "cube-minimal-vertexColors.ply", options)
        
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

    def test_load_ply_from_mesh_lab_with_vertex_colors(self):
        self.mesh.request_vertex_colors()
    
        options = openmesh.Options()
        options += openmesh.Options.VertexColor
    
        ok = openmesh.read_mesh(self.mesh, "meshlab.ply", options)
    
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
    
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[2], 1.0)
        
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

    def test_write_and_read_binary_ply_with_vertex_colors(self):
        self.mesh.request_vertex_colors()
        
        options = openmesh.Options()
        options += openmesh.Options.VertexColor
        
        ok = openmesh.read_mesh(self.mesh, "meshlab.ply", options)
        
        self.assertTrue(ok)

        options += openmesh.Options.Binary

        ok = openmesh.write_mesh(self.mesh, "meshlab_binary.ply", options)
        self.assertTrue(ok)

        self.mesh.clear

        ok = openmesh.read_mesh(self.mesh, "meshlab_binary.ply", options)
        self.assertTrue(ok)

        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)

        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[2], 1.0)

        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[2], 1.0)

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

    def test_write_and_read_ply_with_float_vertex_colors(self):
        self.mesh.request_vertex_colors()
        
        options = openmesh.Options()
        options += openmesh.Options.VertexColor
        
        ok = openmesh.read_mesh(self.mesh, "meshlab.ply", options)
        
        self.assertTrue(ok)
        
        options += openmesh.Options.ColorFloat
    
        ok = openmesh.write_mesh(self.mesh, "meshlab_float.ply", options)
        self.assertTrue(ok)
        
        self.mesh.clear
        ok = openmesh.read_mesh(self.mesh, "meshlab_float.ply", options)
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[2], 1.0)
        
        self.assertFalse(options.vertex_has_normal())
        self.assertFalse(options.vertex_has_texcoord())
        self.assertTrue(options.vertex_has_color())
        self.assertTrue(options.color_is_float())
        
        self.mesh.release_vertex_colors()

    def test_write_and_read_binary_ply_with_float_vertex_colors(self):
        self.mesh.request_vertex_colors()
        
        options = openmesh.Options()
        options += openmesh.Options.VertexColor
        
        ok = openmesh.read_mesh(self.mesh, "meshlab.ply", options)
        
        self.assertTrue(ok)
        
        options += openmesh.Options.ColorFloat
        options += openmesh.Options.Binary
        
        ok = openmesh.write_mesh(self.mesh, "meshlab_binary_float.ply", options)
        self.assertTrue(ok)
        
        self.mesh.clear
        ok = openmesh.read_mesh(self.mesh, "meshlab_binary_float.ply", options)
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[2], 1.0)
        
        self.assertFalse(options.vertex_has_normal())
        self.assertFalse(options.vertex_has_texcoord())
        self.assertTrue(options.vertex_has_color())
        self.assertTrue(options.color_is_float())
        self.assertTrue(options.is_binary())
        
        self.mesh.release_vertex_colors()

    def test_load_simple_ply_with_texcoords(self):
        self.mesh.request_vertex_texcoords2D()
        
        options = openmesh.Options()
        options += openmesh.Options.VertexTexCoord
        
        ok = openmesh.read_mesh(self.mesh, "cube-minimal-texCoords.ply", options)
        
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

    def test_load_simple_ply_with_normals(self):
        self.mesh.request_vertex_normals()
        
        options = openmesh.Options()
        options += openmesh.Options.VertexNormal
        
        ok = openmesh.read_mesh(self.mesh, "cube-minimal-normals.ply", options)
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.assertTrue(options.vertex_has_normal())
        self.assertFalse(options.vertex_has_texcoord())
        self.assertFalse(options.vertex_has_color())
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(0))[0], 0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(0))[2], 1.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(3))[0], 1.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(3))[1], 0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(3))[2], 0.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(4))[0], 1.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(4))[1], 0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(4))[2], 1.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(7))[0], 1.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(7))[1], 1.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(7))[2], 2.0)
        
        self.mesh.release_vertex_normals()


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ReadWritePLY)
    unittest.TextTestRunner(verbosity=2).run(suite)
