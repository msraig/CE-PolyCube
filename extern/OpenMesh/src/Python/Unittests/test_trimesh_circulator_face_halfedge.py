import unittest
import openmesh

class TriMeshCirculatorFaceHalfEdge(unittest.TestCase):

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

    def test_face_halfedge_iter_without_holes_increment(self):
        # Iterate around face 1 at the middle
        fh_it = openmesh.FaceHalfedgeIter(self.mesh, self.mesh.face_handle(1))
        self.assertEqual(fh_it.__next__().idx(), 8)
        self.assertEqual(fh_it.__next__().idx(), 3)
        self.assertEqual(fh_it.__next__().idx(), 6)
        self.assertRaises(StopIteration, fh_it.__next__)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TriMeshCirculatorFaceHalfEdge)
    unittest.TextTestRunner(verbosity=2).run(suite)