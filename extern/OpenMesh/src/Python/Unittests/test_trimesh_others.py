import unittest
import openmesh

from math import pi, fabs

class Others(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()
        self.vhandle = []
    
    def test_is_estimated_feature_edge(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 1)))

        # Add four faces
        face_vhandles = []

        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[0])
        self.mesh.add_face(face_vhandles)

        # ===============================================
        # Setup complete
        # ===============================================


        # Check one Request only vertex normals
        # Face normals are required for vertex and halfedge normals, so
        # that prevent access to non existing properties are in place

        self.mesh.request_vertex_normals()
        self.mesh.request_halfedge_normals()
        self.mesh.request_face_normals()

        # Automatically compute all normals
        # As only vertex normals are requested and no face normals, this will compute nothing.
        self.mesh.update_normals()

        he = self.mesh.halfedges().__next__()

        self.assertTrue(self.mesh.is_estimated_feature_edge(he, 0.0))
        self.assertTrue(self.mesh.is_estimated_feature_edge(he, 0.125 * pi))
        self.assertTrue(self.mesh.is_estimated_feature_edge(he, 0.250 * pi))
        self.assertTrue(self.mesh.is_estimated_feature_edge(he, 0.375 * pi))
        self.assertTrue(self.mesh.is_estimated_feature_edge(he, 0.500 * pi))
        self.assertFalse(self.mesh.is_estimated_feature_edge(he, 0.625 * pi))
        self.assertFalse(self.mesh.is_estimated_feature_edge(he, 0.750 * pi))
        self.assertFalse(self.mesh.is_estimated_feature_edge(he, 0.875 * pi))
        self.assertFalse(self.mesh.is_estimated_feature_edge(he, 1.000 * pi))

    def test_is_estimated_feature_edge(self):
        # Test setup:
        #  1 -- 2
        #  |  / |
        #  | /  |
        #  0 -- 3

        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 0, 0)))

        # Add two faces
        face_vhandles = []

        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)

        # ===============================================
        # Setup complete
        # ===============================================

        he = self.mesh.halfedge_handle(4)

        self.assertEqual(self.mesh.to_vertex_handle(he).idx(), 0)
        self.assertEqual(self.mesh.from_vertex_handle(he).idx(), 2)
        self.assertEqual(self.mesh.edge_handle(he).idx(), 2)

        eh = self.mesh.edge_handle(he)
        self.assertEqual(self.mesh.calc_dihedral_angle(eh), 0.0)

        # Modify point
        tmp = (openmesh.Vec3d(0.0, 0.0, -1.0) + openmesh.Vec3d(1.0, 1.0, -1.0)) * 0.5
        self.mesh.set_point(self.vhandle[2], tmp)

        difference = fabs(1.36944 - self.mesh.calc_dihedral_angle(eh))

        self.assertTrue(difference < 0.00001)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Others)
    unittest.TextTestRunner(verbosity=2).run(suite)
