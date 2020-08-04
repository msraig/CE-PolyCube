import unittest
import openmesh

class TriMeshCirculatorVertexVertex(unittest.TestCase):

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
        Starting vertex is 1->4
        '''

    def test_vertex_vertex_increment(self):
        # Iterate around vertex 1 at the middle
        vv_it = openmesh.VertexVertexIter(self.mesh, self.vhandle[1])
        self.assertEqual(vv_it.__next__().idx(), 4)
        self.assertEqual(vv_it.__next__().idx(), 3)
        self.assertEqual(vv_it.__next__().idx(), 0)
        self.assertEqual(vv_it.__next__().idx(), 2)
        self.assertRaises(StopIteration, vv_it.__next__)

    def test_vertex_vertex_boundary_increment(self):
        # Iterate around vertex 2 at the boundary
        vv_it = openmesh.VertexVertexIter(self.mesh, self.vhandle[2])
        self.assertEqual(vv_it.__next__().idx(), 4)
        self.assertEqual(vv_it.__next__().idx(), 1)
        self.assertEqual(vv_it.__next__().idx(), 0)
        self.assertRaises(StopIteration, vv_it.__next__)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TriMeshCirculatorVertexVertex)
    unittest.TextTestRunner(verbosity=2).run(suite)
