import unittest
import openmesh

class BoundaryTriangleMesh(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()

        # Add some vertices
        self.vhandle = []

        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0,-1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2,-1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(3, 0, 0)))
        
        # Single point
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0,-2, 0)))
        
        # Add five faces
        self.fhandle = []
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        self.fhandle.append(self.mesh.add_face(face_vhandles))
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[4])
        self.fhandle.append(self.mesh.add_face(face_vhandles))
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[1])
        self.fhandle.append(self.mesh.add_face(face_vhandles))
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[4])
        self.fhandle.append(self.mesh.add_face(face_vhandles))
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[5])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[4])
        self.fhandle.append(self.mesh.add_face(face_vhandles))

        # Test setup:
        #  0 ==== 2
        #  |\  0 /|\
        #  | \  / | \
        #  |2  1 3|4 5
        #  | /  \ | /
        #  |/  1 \|/
        #  3 ==== 4
        #
        # Vertex 6 single

    def test_boundary_vertex(self):
        self.assertTrue (self.mesh.is_boundary(self.vhandle[0]))
        self.assertFalse(self.mesh.is_boundary(self.vhandle[1]))
        self.assertTrue (self.mesh.is_boundary(self.vhandle[2]))
        self.assertTrue (self.mesh.is_boundary(self.vhandle[3]))
        self.assertTrue (self.mesh.is_boundary(self.vhandle[4]))
        self.assertTrue (self.mesh.is_boundary(self.vhandle[5]))
        
        self.assertTrue (self.mesh.is_boundary(self.vhandle[6]))
        
    def test_boundary_face(self):
        self.assertTrue (self.mesh.is_boundary(self.fhandle[0]))
        self.assertTrue (self.mesh.is_boundary(self.fhandle[1]))
        self.assertTrue (self.mesh.is_boundary(self.fhandle[2]))
        self.assertFalse(self.mesh.is_boundary(self.fhandle[3]))
        self.assertTrue (self.mesh.is_boundary(self.fhandle[4]))
        

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(BoundaryTriangleMesh)
    unittest.TextTestRunner(verbosity=2).run(suite)
