import unittest
import openmesh

class TriMeshCirculatorVertexIHalfEdge(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()

        # Add some vertices
        self.vhandle = []

        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0,-1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2,-1, 0)))
        
        # Add four faces
        self.mesh.add_face(self.vhandle[0], self.vhandle[1], self.vhandle[2])
        self.mesh.add_face(self.vhandle[1], self.vhandle[3], self.vhandle[4])
        self.mesh.add_face(self.vhandle[0], self.vhandle[3], self.vhandle[1])
        self.mesh.add_face(self.vhandle[2], self.vhandle[1], self.vhandle[4])

        '''
        Test setup:
            0 ==== 2
            |\  0 /|
            | \  / |
            |2  1 3|
            | /  \ |
            |/  1 \|
            3 ==== 4
        Starting halfedge is 1->4
        '''

    def test_vertex_incoming_halfedge_without_holes_increment(self):
        # Iterate around vertex 1 at the middle
        vih_it = openmesh.VertexIHalfedgeIter(self.mesh, self.vhandle[1])
        heh = vih_it.__next__()
        self.assertEqual(heh.idx(), 10)
        self.assertEqual(self.mesh.face_handle(heh).idx(), 1)
        heh = vih_it.__next__()
        self.assertEqual(heh.idx(), 7)
        self.assertEqual(self.mesh.face_handle(heh).idx(), 2)
        heh = vih_it.__next__()
        self.assertEqual(heh.idx(), 0)
        self.assertEqual(self.mesh.face_handle(heh).idx(), 0)
        heh = vih_it.__next__()
        self.assertEqual(heh.idx(), 3)
        self.assertEqual(self.mesh.face_handle(heh).idx(), 3)
        self.assertRaises(StopIteration, vih_it.__next__)
        
    def test_vertex_incoming_halfedge_boundary_increment(self):
        # Iterate around vertex 2 at the boundary
        vih_it = openmesh.VertexIHalfedgeIter(self.mesh, self.vhandle[2])
        heh = vih_it.__next__()
        self.assertEqual(heh.idx(), 14)
        self.assertEqual(self.mesh.face_handle(heh).idx(), 3)
        heh = vih_it.__next__()
        self.assertEqual(heh.idx(), 2)
        self.assertEqual(self.mesh.face_handle(heh).idx(), 0)
        heh = vih_it.__next__()
        self.assertEqual(heh.idx(), 5)
        self.assertEqual(self.mesh.face_handle(heh).idx(), -1)
        self.assertRaises(StopIteration, vih_it.__next__)
        
    def test_vertex_incoming_halfedge_dereference_increment(self):
        # Iterate around vertex 1 at the middle
        vih_it = openmesh.VertexIHalfedgeIter(self.mesh, self.vhandle[1])
        heh = vih_it.__next__()
        eh = self.mesh.edge_handle(heh)
        vh = self.mesh.to_vertex_handle(heh)
        self.assertEqual(heh.idx(), 10)
        self.assertEqual(eh.idx(), 5)
        self.assertEqual(vh.idx(), 1)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TriMeshCirculatorVertexIHalfEdge)
    unittest.TextTestRunner(verbosity=2).run(suite)
