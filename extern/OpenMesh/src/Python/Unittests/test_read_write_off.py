import unittest
import openmesh

class ReadWriteOFF(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()

    def test_load_simple_off_file(self):
        ok = openmesh.read_mesh(self.mesh, "cube1.off")
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 7526)
        self.assertEqual(self.mesh.n_edges(), 22572)
        self.assertEqual(self.mesh.n_faces(), 15048)

    def test_write_and_read_vertex_colors_to_and_from_off_file(self):
        self.mesh.request_vertex_colors()

        self.mesh.add_vertex(openmesh.Vec3d(0, 0, 1))
        self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0))
        self.mesh.add_vertex(openmesh.Vec3d(0, 1, 1))
        self.mesh.add_vertex(openmesh.Vec3d(1, 0, 1))
        
        # Using the default openmesh Python color type
        testColor = openmesh.Vec4f(1.0, 0.5, 0.25, 1.0)
        
        # Setting colors (different from black)
        for v in self.mesh.vertices():
            self.mesh.set_color(v, testColor)
        
        # Check if the colors are correctly set
        count = 0
        for v in self.mesh.vertices():
            color = self.mesh.color(v)
            if color[0] != testColor[0] or color[1] != testColor[1] or color[2] != testColor[2]:
                count += 1

        self.assertEqual(count, 0)

        options = openmesh.Options()
        options += openmesh.Options.VertexColor
        options += openmesh.Options.ColorFloat
        
        openmesh.write_mesh(self.mesh, "temp.off", options)
        openmesh.read_mesh(self.mesh, "temp.off", options)

        # Check if vertices still have the same color
        count = 0
        for v in self.mesh.vertices():
            color = self.mesh.color(v)
            if color[0] != testColor[0] or color[1] != testColor[1] or color[2] != testColor[2]:
                count += 1

        self.assertEqual(count, 0)

        self.mesh.release_vertex_colors()

    def test_write_and_read_float_vertex_colors_to_and_from_off_file(self):
        self.mesh.request_vertex_colors()
        
        options = openmesh.Options(openmesh.Options.VertexColor)
        
        ok = openmesh.read_mesh(self.mesh, "meshlab.ply", options)

        self.assertTrue(ok)

        options.clear()
        options += openmesh.Options.VertexColor
        options += openmesh.Options.ColorFloat

        # Write the mesh
        ok = openmesh.write_mesh(self.mesh, "cube_floating.off", options)
        self.assertTrue(ok)
        ok = openmesh.read_mesh(self.mesh, "cube_floating.off", options)
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

    def test_write_and_read_binary_float_vertex_colors_to_and_from_off_file(self):
        self.mesh.request_vertex_colors()
        
        options = openmesh.Options(openmesh.Options.VertexColor)
        
        ok = openmesh.read_mesh(self.mesh, "meshlab.ply", options)
        
        self.assertTrue(ok)
        
        options.clear()
        options += openmesh.Options.VertexColor
        options += openmesh.Options.Binary
        options += openmesh.Options.ColorFloat
        
        # Write the mesh
        ok = openmesh.write_mesh(self.mesh, "cube_floating_binary.off", options)
        self.assertTrue(ok)
        self.mesh.clear()
        options.clear()
        options += openmesh.Options.VertexColor
        options += openmesh.Options.Binary
        options += openmesh.Options.ColorFloat
        ok = openmesh.read_mesh(self.mesh, "cube_floating_binary.off", options)
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
        self.assertFalse(options.face_has_color())
        self.assertTrue(options.vertex_has_color())
        self.assertTrue(options.color_is_float())
        self.assertTrue(options.is_binary())
        
        self.mesh.release_vertex_colors()


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ReadWriteOFF)
    unittest.TextTestRunner(verbosity=2).run(suite)
