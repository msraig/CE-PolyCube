import unittest
import openmesh

class SplitCopy(unittest.TestCase):

    def test_split_copy_triangle_mesh(self):
        self.mesh = openmesh.TriMesh()
        self.vhandle = []
        
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0.25, 0.25, 0)))

        # Add one face
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[0])
        
        fh = self.mesh.add_face(face_vhandles)

        # Test setup:
        #  1 === 2
        #  |   /
        #  |  /
        #  | /
        #  0
        
        # Set property
        fprop_int = openmesh.FPropHandle()
        self.mesh.add_property(fprop_int)
        self.mesh.set_property(fprop_int, fh, 999)

        # Split face with new vertex
        self.mesh.split_copy(fh, self.vhandle[3])

        # Check setup
        for f in self.mesh.faces():
            self.assertEqual(self.mesh.property(fprop_int, f), 999)

    def test_split_copy_polymesh(self):
        self.mesh = openmesh.PolyMesh()
        self.vhandle = []
        
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0.5, 0.5, 0)))
        
        # Add one face
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[3])
        
        fh = self.mesh.add_face(face_vhandles)
        
        # Test setup:
        #  1 === 2
        #  |     |
        #  |     |
        #  |     |
        #  0 === 3
        
        # Set property
        fprop_int = openmesh.FPropHandle()
        self.mesh.add_property(fprop_int)
        self.mesh.set_property(fprop_int, fh, 999)
        
        # Split face with new vertex
        self.mesh.split_copy(fh, self.vhandle[4])
        
        # Check setup
        for f in self.mesh.faces():
            self.assertEqual(self.mesh.property(fprop_int, f), 999)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(SplitCopy)
    unittest.TextTestRunner(verbosity=2).run(suite)
