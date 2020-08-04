import unittest
import openmesh

class TriMeshIterators(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()
        self.vhandle = []

    def test_vertex_iter(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 0, 0)))
        
        # Add two faces
        self.mesh.add_face(self.vhandle[2], self.vhandle[1], self.vhandle[0])
        self.mesh.add_face(self.vhandle[2], self.vhandle[0], self.vhandle[3])
        
        # Test setup:
        #  1 === 2
        #  |   / |    
        #  |  /  |
        #  | /   |
        #  0 === 3
        
        v_it = self.mesh.vertices()
        
        self.assertEqual(v_it.__next__().idx(), 0)
        self.assertEqual(v_it.__next__().idx(), 1)
        self.assertEqual(v_it.__next__().idx(), 2)
        self.assertEqual(v_it.__next__().idx(), 3)
        
        self.assertRaises(StopIteration, v_it.__next__)
        
    def test_vertex_iter_start_position(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 0, 0)))
        
        # Add two faces
        self.mesh.add_face(self.vhandle[2], self.vhandle[1], self.vhandle[0])
        self.mesh.add_face(self.vhandle[2], self.vhandle[0], self.vhandle[3])
        
        # Test setup:
        #  1 === 2
        #  |   / |    
        #  |  /  |
        #  | /   |
        #  0 === 3
        
        v_it = openmesh.VertexIter(self.mesh, self.mesh.vertex_handle(2))
        
        self.assertEqual(v_it.__next__().idx(), 2)
        self.assertEqual(v_it.__next__().idx(), 3)
        
        self.assertRaises(StopIteration, v_it.__next__)
        
    def test_edge_iter(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 0, 0)))
        
        # Add two faces
        self.mesh.add_face(self.vhandle[2], self.vhandle[1], self.vhandle[0])
        self.mesh.add_face(self.vhandle[2], self.vhandle[0], self.vhandle[3])
        
        # Test setup:
        #  1 === 2
        #  |   / |    
        #  |  /  |
        #  | /   |
        #  0 === 3
        
        e_it = self.mesh.edges()
        
        e = e_it.__next__()
        self.assertEqual(e.idx(), 0)
        
        he = self.mesh.halfedge_handle(e, 0)
        self.assertEqual(self.mesh.to_vertex_handle(he).idx(), 1)
        self.assertEqual(self.mesh.from_vertex_handle(he).idx(), 2)
        he = self.mesh.halfedge_handle(e, 1)
        self.assertEqual(self.mesh.to_vertex_handle(he).idx(), 2)
        self.assertEqual(self.mesh.from_vertex_handle(he).idx(), 1)
        
        e = e_it.__next__()
        self.assertEqual(e.idx(), 1)
        
        he = self.mesh.halfedge_handle(e, 0)
        self.assertEqual(self.mesh.to_vertex_handle(he).idx(), 0)
        self.assertEqual(self.mesh.from_vertex_handle(he).idx(), 1)
        he = self.mesh.halfedge_handle(e, 1)
        self.assertEqual(self.mesh.to_vertex_handle(he).idx(), 1)
        self.assertEqual(self.mesh.from_vertex_handle(he).idx(), 0)
        
        e = e_it.__next__()
        self.assertEqual(e.idx(), 2)
        
        he = self.mesh.halfedge_handle(e, 0)
        self.assertEqual(self.mesh.to_vertex_handle(he).idx(), 2)
        self.assertEqual(self.mesh.from_vertex_handle(he).idx(), 0)
        he = self.mesh.halfedge_handle(e, 1)
        self.assertEqual(self.mesh.to_vertex_handle(he).idx(), 0)
        self.assertEqual(self.mesh.from_vertex_handle(he).idx(), 2)
        
        e = e_it.__next__()
        self.assertEqual(e.idx(), 3)
        
        he = self.mesh.halfedge_handle(e, 0)
        self.assertEqual(self.mesh.to_vertex_handle(he).idx(), 3)
        self.assertEqual(self.mesh.from_vertex_handle(he).idx(), 0)
        he = self.mesh.halfedge_handle(e, 1)
        self.assertEqual(self.mesh.to_vertex_handle(he).idx(), 0)
        self.assertEqual(self.mesh.from_vertex_handle(he).idx(), 3)
        
        e = e_it.__next__()
        self.assertEqual(e.idx(), 4)
        
        he = self.mesh.halfedge_handle(e, 0)
        self.assertEqual(self.mesh.to_vertex_handle(he).idx(), 2)
        self.assertEqual(self.mesh.from_vertex_handle(he).idx(), 3)
        he = self.mesh.halfedge_handle(e, 1)
        self.assertEqual(self.mesh.to_vertex_handle(he).idx(), 3)
        self.assertEqual(self.mesh.from_vertex_handle(he).idx(), 2)
    
    def test_halfedge_iter_skipping(self):
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
        self.mesh.add_face(self.vhandle[0], self.vhandle[1], self.vhandle[3])
        self.mesh.add_face(self.vhandle[1], self.vhandle[2], self.vhandle[3])
        self.mesh.add_face(self.vhandle[7], self.vhandle[6], self.vhandle[5])
        self.mesh.add_face(self.vhandle[7], self.vhandle[5], self.vhandle[4])
        self.mesh.add_face(self.vhandle[1], self.vhandle[0], self.vhandle[4])
        self.mesh.add_face(self.vhandle[1], self.vhandle[4], self.vhandle[5])
        self.mesh.add_face(self.vhandle[2], self.vhandle[1], self.vhandle[5])
        self.mesh.add_face(self.vhandle[2], self.vhandle[5], self.vhandle[6])
        self.mesh.add_face(self.vhandle[3], self.vhandle[2], self.vhandle[6])
        self.mesh.add_face(self.vhandle[3], self.vhandle[6], self.vhandle[7])
        self.mesh.add_face(self.vhandle[0], self.vhandle[3], self.vhandle[7])
        self.mesh.add_face(self.vhandle[0], self.vhandle[7], self.vhandle[4])
        
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
        
        # Run over all halfedges
        heCounter = 0

        self.mesh.request_face_status()
        self.mesh.request_vertex_status()
        self.mesh.request_halfedge_status()
        
        # Get second edge
        eh = self.mesh.edge_handle(2)
        
        # Delete one edge
        self.mesh.delete_edge(eh)
        
        # Check setup ( No garbage collection, so nothing should change!)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_halfedges(), 36)
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        # =====================================================
        # Check skipping iterator
        # =====================================================
        
        ok_4 = True
        ok_5 = True
        
        count = 0
        
        for he in self.mesh.shalfedges():
            if he.idx() == 4:
                ok_4 = False
            if he.idx() == 5:
                ok_5 = False
            count += 1
        
        self.assertEqual(count, 34)
        self.assertTrue(ok_4)
        self.assertTrue(ok_5)
        
        # =====================================================
        # Check non skipping iterator
        # =====================================================
        
        ok_4 = False
        ok_5 = False
        
        count = 0
        
        for he in self.mesh.halfedges():
            if he.idx() == 4:
                ok_4 = True
            if he.idx() == 5:
                ok_5 = True
            count += 1
        
        self.assertEqual(count, 36)
        self.assertTrue(ok_4)
        self.assertTrue(ok_5)
        
    def test_halfedge_iter_skipping_low_level(self):
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
        self.mesh.add_face(self.vhandle[0], self.vhandle[1], self.vhandle[3])
        self.mesh.add_face(self.vhandle[1], self.vhandle[2], self.vhandle[3])
        self.mesh.add_face(self.vhandle[7], self.vhandle[6], self.vhandle[5])
        self.mesh.add_face(self.vhandle[7], self.vhandle[5], self.vhandle[4])
        self.mesh.add_face(self.vhandle[1], self.vhandle[0], self.vhandle[4])
        self.mesh.add_face(self.vhandle[1], self.vhandle[4], self.vhandle[5])
        self.mesh.add_face(self.vhandle[2], self.vhandle[1], self.vhandle[5])
        self.mesh.add_face(self.vhandle[2], self.vhandle[5], self.vhandle[6])
        self.mesh.add_face(self.vhandle[3], self.vhandle[2], self.vhandle[6])
        self.mesh.add_face(self.vhandle[3], self.vhandle[6], self.vhandle[7])
        self.mesh.add_face(self.vhandle[0], self.vhandle[3], self.vhandle[7])
        self.mesh.add_face(self.vhandle[0], self.vhandle[7], self.vhandle[4])
        
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
        
        # Run over all halfedges
        heCounter = 0

        self.mesh.request_face_status()
        self.mesh.request_vertex_status()
        self.mesh.request_halfedge_status()
        
        # Get second edge
        eh = self.mesh.edge_handle(2)
        
        # Delete one edge
        self.mesh.delete_edge(eh)
        
        # Check setup ( No garbage collection, so nothing should change!)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_halfedges(), 36)
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        # =====================================================
        # Try to add low level edge with invalid incidents and
        # check skipping iterator
        # =====================================================
        
        # Add a low level edge without handles
        eh_test = self.mesh.edge_handle(self.mesh.new_edge(openmesh.VertexHandle(), openmesh.VertexHandle()))
        
        count = 0
        found_4  = False
        found_5  = False
        found_36 = False
        found_37 = False
        
        for he in self.mesh.shalfedges():
            if he.idx() == 4:
                found_4 = True
            if he.idx() == 5:
                found_5 = True
            if he.idx() == 36:
                found_36 = True
            if he.idx() == 37:
                found_37 = True
            count += 1
            
        self.assertEqual(count, 36)
        self.assertFalse(found_4)
        self.assertFalse(found_5)
        self.assertTrue(found_36)
        self.assertTrue(found_37)
        
        # =====================================================
        # Try to delete one edge with invalid incidents and
        # check skipping iterator
        # =====================================================
        
        # Delete one edge and recheck (Halfedges 4 and 5)
        self.mesh.delete_edge(eh_test)
        
        count = 0
        found_4  = False
        found_5  = False
        found_36 = False
        found_37 = False
        
        for he in self.mesh.shalfedges():
            if he.idx() == 4:
                found_4 = True
            if he.idx() == 5:
                found_5 = True
            if he.idx() == 36:
                found_36 = True
            if he.idx() == 37:
                found_37 = True
            count += 1
        
        self.assertEqual(count, 34)
        self.assertFalse(found_4)
        self.assertFalse(found_5)
        self.assertFalse(found_36)
        self.assertFalse(found_37)
        
    def test_face_iter_empty_mesh_one_deleted_face(self):
        # Request delete_face capability
        self.mesh.request_vertex_status()
        self.mesh.request_edge_status()
        self.mesh.request_face_status()
        
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 0)))
        
        # Add one face
        fh = self.mesh.add_face(self.vhandle[2], self.vhandle[1], self.vhandle[0])
        
        is_delete_isolated_vertex = False
        self.mesh.delete_face(fh, is_delete_isolated_vertex)
        
        # Test setup:
        #  1 === 2
        #  |   /
        #  |  /
        #  | /
        #  0
        
        # Normal iterators
        f_it = self.mesh.faces()
        
        self.assertEqual(f_it.__next__().idx(), 0)
        self.assertRaises(StopIteration, f_it.__next__)
        
        # Same with skipping iterators
        f_it = self.mesh.sfaces()
        
        self.assertRaises(StopIteration, f_it.__next__)
        

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TriMeshIterators)
    unittest.TextTestRunner(verbosity=2).run(suite)
