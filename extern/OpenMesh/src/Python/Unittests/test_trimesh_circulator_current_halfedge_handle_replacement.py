import unittest
import openmesh

class TrimeshCirculatorCurrentHalfedgeHandleReplacement(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()
        self.vhandle = []
    
    def test_dereference(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0,  1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1,  0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2,  1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, -1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2, -1, 0)))

        # Add four faces
        face_vhandles = []

        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[1])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)

        # Test setup:
        #  0 ==== 2
        #  |\  0 /|
        #  | \  / |
        #  |2  1 3|
        #  | /  \ |
        #  |/  1 \|
        #  3 ==== 4
        # Starting vertex is 1->4
        
        # output from fh_it.current_halfedge_handle()
        current_halfedge_handles = [4, 0, 2, 10, 6, 8, 1, 12, 7, 14, 3, 11]
        
        i = 0
        for f in self.mesh.faces():
            for he in self.mesh.fh(f):
                self.assertEqual(he.idx(), current_halfedge_handles[i])
                i += 1

    def test_vv_iter(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0,  1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1,  0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2,  1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, -1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2, -1, 0)))
        
        # Add four faces
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[1])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)

        # Test setup:
        #  0 ==== 2
        #  |\  0 /|
        #  | \  / |
        #  |2  1 3|
        #  | /  \ |
        #  |/  1 \|
        #  3 ==== 4
        # Starting vertex is 1->4

        # output from vv_it.current_halfedge_handle()
        current_halfedge_handles = [5, 0, 12, 11, 6, 1, 2, 15, 3, 4, 13, 7, 8, 9, 10, 14]

        eh0 = []
        eh1 = []

        i = 0

        for v in self.mesh.vertices():
            for vv in self.mesh.vv(v):
                he = openmesh.HalfedgeHandle(current_halfedge_handles[i])
                eh0.append(self.mesh.edge_handle(he))
                i += 1
        for v in self.mesh.vertices():
            for he in self.mesh.voh(v):
                eh1.append(self.mesh.edge_handle(he))

        self.assertEqual(len(eh0), len(eh1))
        for i in range(len(eh0)):
            self.assertEqual(eh0[i], eh1[i])

    def test_fe_iter(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0,  1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1,  0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2,  1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, -1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2, -1, 0)))
        
        # Add four faces
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[1])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        # Test setup:
        #  0 ==== 2
        #  |\  0 /|
        #  | \  / |
        #  |2  1 3|
        #  | /  \ |
        #  |/  1 \|
        #  3 ==== 4
        # Starting vertex is 1->4
        
        # output from fe_it.current_halfedge_handle()
        current_halfedge_handles = [4, 0, 2, 10, 6, 8, 1, 12, 7, 14, 3, 11]
        
        heh0 = []
        heh1 = []
        
        i = 0
        
        for f in self.mesh.faces():
            for e in self.mesh.fe(f):
                heh0.append(openmesh.HalfedgeHandle(current_halfedge_handles[i]))
                i += 1
        for f in self.mesh.faces():
            for he in self.mesh.fh(f):
                heh1.append(he)
        
        self.assertEqual(len(heh0), len(heh1))
        for i in range(len(heh0)):
            self.assertEqual(heh0[i], heh1[i])

    def test_vf_iter_boundary(self):
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
        # 0 ------ 2 ------ 4
        #  \      / \      /
        #   \  0 /   \  1 /
        #    \  /     \  /
        #     1 ------- 3
        #      \       /
        #       \  2  /
        #        \   /
        #         \ /
        #          5
        
        # output from fe_it.current_halfedge_handle()
        current_halfedge_handles = [0, 2, 12, 4, 6, 8, 16, 10, 14]
        
        fh0 = []
        fh1 = []
        
        i = 0
        
        for v in self.mesh.vertices():
            for f in self.mesh.vf(v):
                he = openmesh.HalfedgeHandle(current_halfedge_handles[i])
                fh0.append(self.mesh.face_handle(he))
                i += 1
        for v in self.mesh.vertices():
            for f in self.mesh.vf(v):
                fh1.append(f)
        
        self.assertEqual(len(fh0), len(fh1))
        for i in range(len(fh0)):
            self.assertEqual(fh0[i], fh1[i])


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TrimeshCirculatorCurrentHalfedgeHandleReplacement)
    unittest.TextTestRunner(verbosity=2).run(suite)
