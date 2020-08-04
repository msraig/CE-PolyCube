import unittest
import openmesh

class TriMeshCirculatorFaceEdge(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()

        # Add some vertices
        self.vhandle = []

        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(3, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(4, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2,-1, 0)))
        
        # Add four faces
        self.mesh.add_face(self.vhandle[0], self.vhandle[1], self.vhandle[2])
        self.mesh.add_face(self.vhandle[2], self.vhandle[1], self.vhandle[3])
        self.mesh.add_face(self.vhandle[2], self.vhandle[3], self.vhandle[4])
        self.mesh.add_face(self.vhandle[1], self.vhandle[5], self.vhandle[3])

        '''
        Test setup:
            0 ------ 2 ------ 4
             \      / \      /
              \  0 /   \  2 /
               \  /  1  \  /
                1 ------- 3
                 \       /
                  \  3  /
                   \   /
                    \ /
                     5
        '''

    def test_face_edge_iter_without_holes_increment(self):
        # Iterate around face 1 at the middle
        fe_it = openmesh.FaceEdgeIter(self.mesh, self.mesh.face_handle(1))
        self.assertEqual(fe_it.__next__().idx(), 4)
        self.assertEqual(fe_it.__next__().idx(), 1)
        self.assertEqual(fe_it.__next__().idx(), 3)
        self.assertRaises(StopIteration, fe_it.__next__)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TriMeshCirculatorFaceEdge)
    unittest.TextTestRunner(verbosity=2).run(suite)