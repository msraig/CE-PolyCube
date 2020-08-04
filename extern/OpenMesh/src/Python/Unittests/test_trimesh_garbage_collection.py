import unittest
import openmesh

class TriMeshGarbageCollection(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()
        
        self.mesh.request_vertex_status()
        self.mesh.request_edge_status()
        self.mesh.request_halfedge_status()
        self.mesh.request_face_status()

        # Add some vertices
        self.vhandle = []
        
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
        
    def test_standard_garbage_collection(self):
        # Check setup
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.mesh.delete_vertex(self.vhandle[0])
        
        # Check setup
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.mesh.garbage_collection()
        
        # Check setup
        self.assertEqual(self.mesh.n_vertices(), 7)
        self.assertEqual(self.mesh.n_faces(), 8)
        
    def test_tracked_garbage_collection(self):
        # Check setup
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 12)

        #==================================================
        # Create lists containing the current handles
        #==================================================

        vertexHandles = []
        for v in self.mesh.vertices():
            vertexHandles.append(v)

        halfedgeHandles = []
        for he in self.mesh.halfedges():
            halfedgeHandles.append(he)

        faceHandles = []
        for f in self.mesh.faces():
            faceHandles.append(f)

        # Deleting vertex 0
        self.mesh.delete_vertex(self.vhandle[0])

        # Check setup
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_faces(), 12)

        self.mesh.garbage_collection(vertexHandles, halfedgeHandles, faceHandles, True, True, True)

        # Check setup
        self.assertEqual(self.mesh.n_vertices(), 7)
        self.assertEqual(self.mesh.n_faces(), 8)

        # Check setup of vertices
        self.assertEqual(vertexHandles[0].idx(), -1)
        self.assertEqual(vertexHandles[1].idx(),  1)
        self.assertEqual(vertexHandles[2].idx(),  2)
        self.assertEqual(vertexHandles[3].idx(),  3)
        self.assertEqual(vertexHandles[4].idx(),  4)
        self.assertEqual(vertexHandles[5].idx(),  5)
        self.assertEqual(vertexHandles[6].idx(),  6)
        self.assertEqual(vertexHandles[7].idx(),  0)

        # Check setup of halfedge handles
        self.assertEqual(halfedgeHandles[0 ].idx(), -1)
        self.assertEqual(halfedgeHandles[1 ].idx(), -1)
        self.assertEqual(halfedgeHandles[2 ].idx(),  2)
        self.assertEqual(halfedgeHandles[3 ].idx(),  3)
        self.assertEqual(halfedgeHandles[4 ].idx(), -1)
        self.assertEqual(halfedgeHandles[5 ].idx(), -1)
        self.assertEqual(halfedgeHandles[6 ].idx(),  6)
        self.assertEqual(halfedgeHandles[7 ].idx(),  7)
        self.assertEqual(halfedgeHandles[8 ].idx(),  8)
        self.assertEqual(halfedgeHandles[9 ].idx(),  9)
        self.assertEqual(halfedgeHandles[10].idx(), 10)
        self.assertEqual(halfedgeHandles[11].idx(), 11)
        self.assertEqual(halfedgeHandles[12].idx(), 12)
        self.assertEqual(halfedgeHandles[13].idx(), 13)
        self.assertEqual(halfedgeHandles[14].idx(), 14)
        self.assertEqual(halfedgeHandles[15].idx(), 15)
        self.assertEqual(halfedgeHandles[16].idx(), 16)
        self.assertEqual(halfedgeHandles[17].idx(), 17)
        self.assertEqual(halfedgeHandles[18].idx(), 18)
        self.assertEqual(halfedgeHandles[19].idx(), 19)
        self.assertEqual(halfedgeHandles[20].idx(), -1)
        self.assertEqual(halfedgeHandles[21].idx(), -1)
        self.assertEqual(halfedgeHandles[22].idx(), 22)
        self.assertEqual(halfedgeHandles[23].idx(), 23)
        self.assertEqual(halfedgeHandles[24].idx(), 24)
        self.assertEqual(halfedgeHandles[25].idx(), 25)
        self.assertEqual(halfedgeHandles[26].idx(), 26)
        self.assertEqual(halfedgeHandles[27].idx(), 27)
        self.assertEqual(halfedgeHandles[28].idx(), 20)
        self.assertEqual(halfedgeHandles[29].idx(), 21)
        self.assertEqual(halfedgeHandles[30].idx(),  4)
        self.assertEqual(halfedgeHandles[31].idx(),  5)
        self.assertEqual(halfedgeHandles[32].idx(),  0)
        self.assertEqual(halfedgeHandles[33].idx(),  1)
        self.assertEqual(halfedgeHandles[34].idx(), -1)
        self.assertEqual(halfedgeHandles[35].idx(), -1)

        # Check setup of faces
        self.assertEqual(faceHandles[0 ].idx(), -1)
        self.assertEqual(faceHandles[1 ].idx(),  1)
        self.assertEqual(faceHandles[2 ].idx(),  2)
        self.assertEqual(faceHandles[3 ].idx(),  3)
        self.assertEqual(faceHandles[4 ].idx(), -1)
        self.assertEqual(faceHandles[5 ].idx(),  5)
        self.assertEqual(faceHandles[6 ].idx(),  6)
        self.assertEqual(faceHandles[7 ].idx(),  7)
        self.assertEqual(faceHandles[8 ].idx(),  4)
        self.assertEqual(faceHandles[9 ].idx(),  0)
        self.assertEqual(faceHandles[10].idx(), -1)
        self.assertEqual(faceHandles[11].idx(), -1)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TriMeshGarbageCollection)
    unittest.TextTestRunner(verbosity=2).run(suite)
