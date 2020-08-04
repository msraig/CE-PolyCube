import unittest
import openmesh

class TrimeshCirculatorHalfedgeLoop(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()
        self.vhandle = []
    
    def test_halfedge_loop_with_face(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0,  1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1,  0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2,  1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(3,  0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(4,  1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2, -1, 0)))

        # Add four faces
        face_vhandles = []

        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[5])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)

        # Test setup:
        #
        # edge x => halfedge x/x+1
        # i.e. edge 0 => halfedge 0/1
        #
        # 0 --4--- 2 ------ 4
        #  \      / \      /
        #   0  0 2   6  2 /
        #    \  /  1  \  /
        #     1 ---8--- 3
        #      \       /
        #       \  3  /
        #        \   /
        #         \ /
        #          5
        
        # Circle around face 1
        hl_it = self.mesh.hl(self.mesh.halfedge_handle(3))

        self.assertEqual(hl_it.__next__().idx(), 3)
        self.assertEqual(hl_it.__next__().idx(), 6)
        self.assertEqual(hl_it.__next__().idx(), 8)
        self.assertRaises(StopIteration, hl_it.__next__)

    def test_halfedge_loop_without_face(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0,  1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1,  0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2,  1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(3,  0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(4,  1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2, -1, 0)))
        
        # Add three faces
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[5])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)

        # Test setup:
        #
        # H => hole (no face)
        # fx => face #x
        # edge 0 => halfedge 0/1
        #
        # 0 --4--- 2 -10--- 4
        #  \      / \      /
        #   0 f0 2   6 f2 8
        #    \  /  H  \  /
        #     1 ---16---3
        #      \       /
        #      12 f3 14
        #        \   /
        #         \ /
        #          5

        # Circle around the hole
        hl_it = self.mesh.hl(self.mesh.halfedge_handle(3))

        self.assertEqual(hl_it.__next__().idx(), 3)
        self.assertEqual(hl_it.__next__().idx(), 17)
        self.assertEqual(hl_it.__next__().idx(), 7)
        self.assertRaises(StopIteration, hl_it.__next__)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TrimeshCirculatorHalfedgeLoop)
    unittest.TextTestRunner(verbosity=2).run(suite)
