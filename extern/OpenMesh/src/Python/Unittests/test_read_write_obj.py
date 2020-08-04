import unittest
import openmesh

class ReadWriteOBJ(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()

    def test_load_simple_obj(self):
        ok = openmesh.read_mesh(self.mesh, "cube-minimal.obj")
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)

    def test_load_simple_obj_check_halfedge_and_vertex_normals(self):
        self.mesh.request_halfedge_normals()
        self.mesh.request_vertex_normals()
        
        options = openmesh.Options()
        options += openmesh.Options.VertexNormal
        
        file_name = "cube-minimal.obj"
        
        ok = openmesh.read_mesh(self.mesh, file_name, options)
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        self.assertEqual(self.mesh.n_halfedges(), 36)
        
        # =====================================================
        # Check vertex normals
        # =====================================================
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(0))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(0))[1], -1.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(0))[2],  0.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(3))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(3))[1],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(3))[2],  1.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(4))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(4))[1], -1.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(4))[2],  0.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(7))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(7))[1],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(7))[2],  1.0)
        
        # =====================================================
        # Check halfedge normals
        # =====================================================
        
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle( 0))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle( 0))[1],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle( 0))[2], -1.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(10))[0], -1.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(10))[1],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(10))[2],  0.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(19))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(19))[1],  1.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(19))[2],  0.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(24))[0],  1.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(24))[1],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(24))[2],  0.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(30))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(30))[1], -1.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(30))[2],  0.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(35))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(35))[1],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(35))[2],  1.0)
        
        self.mesh.release_vertex_normals()
        self.mesh.release_halfedge_normals()

    def test_load_simple_obj_force_vertex_colors_although_not_available(self):
        self.mesh.request_vertex_colors()
        
        file_name = "cube-minimal.obj"
        
        options = openmesh.Options()
        options += openmesh.Options.VertexColor
        
        ok = openmesh.read_mesh(self.mesh, file_name, options)
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        self.assertEqual(self.mesh.n_halfedges(), 36)

    def test_load_simple_obj_check_texcoords(self):
        self.mesh.request_halfedge_texcoords2D()
        
        options = openmesh.Options()
        options += openmesh.Options.FaceTexCoord
        
        file_name = "cube-minimal-texCoords.obj"
        
        ok = openmesh.read_mesh(self.mesh, file_name, options)
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.texcoord2D(self.mesh.halfedge_handle( 0))[0], 1.0)
        self.assertEqual(self.mesh.texcoord2D(self.mesh.halfedge_handle( 0))[1], 1.0)
        
        self.assertEqual(self.mesh.texcoord2D(self.mesh.halfedge_handle(10))[0], 3.0)
        self.assertEqual(self.mesh.texcoord2D(self.mesh.halfedge_handle(10))[1], 3.0)
        
        self.assertEqual(self.mesh.texcoord2D(self.mesh.halfedge_handle(19))[0], 6.0)
        self.assertEqual(self.mesh.texcoord2D(self.mesh.halfedge_handle(19))[1], 6.0)
        
        self.assertEqual(self.mesh.texcoord2D(self.mesh.halfedge_handle(24))[0], 7.0)
        self.assertEqual(self.mesh.texcoord2D(self.mesh.halfedge_handle(24))[1], 7.0)
        
        self.assertEqual(self.mesh.texcoord2D(self.mesh.halfedge_handle(30))[0], 9.0)
        self.assertEqual(self.mesh.texcoord2D(self.mesh.halfedge_handle(30))[1], 9.0)
        
        self.assertEqual(self.mesh.texcoord2D(self.mesh.halfedge_handle(35))[0], 12.0)
        self.assertEqual(self.mesh.texcoord2D(self.mesh.halfedge_handle(35))[1], 12.0)
        
        self.mesh.release_halfedge_texcoords2D()

    def test_load_obj_with_material(self):
        self.mesh.request_face_colors()

        options = openmesh.Options()
        options += openmesh.Options.FaceColor

        file_name = "square_material.obj"

        ok = openmesh.read_mesh(self.mesh, file_name, options)

        self.assertTrue(ok)

        fh = self.mesh.face_handle(self.mesh.halfedge_handle(0))

        self.assertTrue(fh.is_valid())

        self.assertAlmostEqual(self.mesh.color(fh)[0], 0.5, 2)
        self.assertAlmostEqual(self.mesh.color(fh)[1], 0.5, 2)
        self.assertAlmostEqual(self.mesh.color(fh)[2], 0.5, 2)

        self.mesh.release_face_colors()

    def test_load_simple_obj_with_vertex_colors_after_vertices(self):
        self.mesh.request_vertex_colors()
        
        options = openmesh.Options()
        options += openmesh.Options.VertexColor
        
        file_name = "cube-minimal-vertex-colors-after-vertex-definition.obj"
        
        ok = openmesh.read_mesh(self.mesh, file_name, options)
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[2], 0.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[1], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[0], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[2], 0.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[0], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[1], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[2], 1.0)
        
        self.mesh.release_vertex_colors()

    def test_load_simple_obj_with_vertex_colors_as_vc_lines(self):
        self.mesh.request_vertex_colors()
        
        options = openmesh.Options()
        options += openmesh.Options.VertexColor
        
        file_name = "cube-minimal-vertex-colors-as-vc-lines.obj"
        
        ok = openmesh.read_mesh(self.mesh, file_name, options)
        
        self.assertTrue(ok)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[2], 0.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[1], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[0], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[2], 0.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[0], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[1], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[2], 1.0)
        
        self.mesh.release_vertex_colors()


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ReadWriteOBJ)
    unittest.TextTestRunner(verbosity=2).run(suite)
