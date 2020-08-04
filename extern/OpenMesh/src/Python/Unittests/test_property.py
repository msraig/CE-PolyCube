import unittest
import openmesh

class Property(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()
        self.vhandle = []

    def test_vertex_property_copy_properties_int(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 0, 0)))
        
        # Add two faces
        face_vhandles = []
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[0])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)
        
        # Test setup:
        #  1 === 2
        #  |   / |
        #  |  /  |
        #  | /   |
        #  0 === 3
        
        # Check setup
        self.assertEqual(self.mesh.n_vertices(), 4)
        self.assertEqual(self.mesh.n_faces(), 2)
        
        # Add a vertex property
        intHandle = openmesh.VPropHandle()
        self.assertFalse(self.mesh.get_property_handle(intHandle, "intProp"))
        self.mesh.add_property(intHandle, "intProp")
        self.assertTrue(self.mesh.get_property_handle(intHandle, "intProp"))
        
        # Fill property
        for vh in self.mesh.vertices():
            self.mesh.set_property(intHandle, vh, vh.idx())
        
        # Check if property it is ok
        v_it = self.mesh.vertices()
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 0)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 1)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 2)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 3)
        
        # Check vertex positions
        v_it = self.mesh.vertices()
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 0)
        self.assertEqual(self.mesh.point(vh)[1], 0)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 0)
        self.assertEqual(self.mesh.point(vh)[1], 1)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 1)
        self.assertEqual(self.mesh.point(vh)[1], 1)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 1)
        self.assertEqual(self.mesh.point(vh)[1], 0)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        
        # Copy from vertex 1 to 0, with skipping build in properties
        self.mesh.copy_all_properties(self.vhandle[1], self.vhandle[0])
        v_it = self.mesh.vertices()
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 0)
        self.assertEqual(self.mesh.point(vh)[1], 0)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 0)
        self.assertEqual(self.mesh.point(vh)[1], 1)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 1)
        self.assertEqual(self.mesh.point(vh)[1], 1)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 1)
        self.assertEqual(self.mesh.point(vh)[1], 0)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        v_it = self.mesh.vertices()
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 1)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 1)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 2)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 3)

        # Copy from vertex 2 to 3, including build in properties
        self.mesh.copy_all_properties(self.vhandle[2], self.vhandle[3], True)
        v_it = self.mesh.vertices()
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 0)
        self.assertEqual(self.mesh.point(vh)[1], 0)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 0)
        self.assertEqual(self.mesh.point(vh)[1], 1)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 1)
        self.assertEqual(self.mesh.point(vh)[1], 1)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 1)
        self.assertEqual(self.mesh.point(vh)[1], 1)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        v_it = self.mesh.vertices()
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 1)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 1)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 2)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 2)
        
    def test_check_status_properties_halfedge_edge_all_deleted(self):
        self.mesh.request_vertex_status()
        self.mesh.request_face_status()
        self.mesh.request_halfedge_status()
        self.mesh.request_edge_status()

        # Define positions
        p1 = openmesh.Vec3d(0, 0, 0)
        p2 = openmesh.Vec3d(0, 1, 0)
        p3 = openmesh.Vec3d(1, 1, 0)
        p4 = openmesh.Vec3d(0, 0, 1)
        
        # Add some vertices
        vh1 = self.mesh.add_vertex(p1)
        vh2 = self.mesh.add_vertex(p2)
        vh3 = self.mesh.add_vertex(p3)
        vh4 = self.mesh.add_vertex(p4)

        # Add some faces
        f1 = self.mesh.add_face(vh1, vh3, vh2)
        f2 = self.mesh.add_face(vh1, vh2, vh4)
        f3 = self.mesh.add_face(vh2, vh3, vh4)
        f4 = self.mesh.add_face(vh3, vh1, vh4)

        # Delete all faces
        self.mesh.delete_face(f1)
        self.mesh.delete_face(f2)
        self.mesh.delete_face(f3)
        self.mesh.delete_face(f4)

        for heh in self.mesh.halfedges():
            self.assertTrue(self.mesh.status(self.mesh.edge_handle(heh)).deleted())
            self.assertTrue(self.mesh.status(heh).deleted())

    def test_copy_all_properties_vertex_after_remove_of_property(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 0, 0)))
        
        # Add two faces
        face_vhandles = []
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[0])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)
        
        # Test setup:
        #  1 === 2
        #  |   / |
        #  |  /  |
        #  | /   |
        #  0 === 3
        
        # Check setup
        self.assertEqual(self.mesh.n_vertices(), 4)
        self.assertEqual(self.mesh.n_faces(), 2)
        
        # Add a double vertex property
        doubleHandle = openmesh.VPropHandle()
        self.assertFalse(self.mesh.get_property_handle(doubleHandle, "doubleProp"))
        self.mesh.add_property(doubleHandle, "doubleProp")
        self.assertTrue(self.mesh.get_property_handle(doubleHandle, "doubleProp"))
        
        # Add a int vertex property
        intHandle = openmesh.VPropHandle()
        self.assertFalse(self.mesh.get_property_handle(intHandle, "intProp"))
        self.mesh.add_property(intHandle, "intProp")
        self.assertTrue(self.mesh.get_property_handle(intHandle, "intProp"))
        
        # Now remove the double property again
        self.mesh.remove_property(doubleHandle)
        
        # Fill int property
        for vh in self.mesh.vertices():
            self.mesh.set_property(intHandle, vh, vh.idx())
        
        # Check if property it is ok
        v_it = self.mesh.vertices()
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 0)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 1)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 2)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 3)

        # Check vertex positions
        v_it = self.mesh.vertices()
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 0)
        self.assertEqual(self.mesh.point(vh)[1], 0)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 0)
        self.assertEqual(self.mesh.point(vh)[1], 1)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 1)
        self.assertEqual(self.mesh.point(vh)[1], 1)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 1)
        self.assertEqual(self.mesh.point(vh)[1], 0)
        self.assertEqual(self.mesh.point(vh)[2], 0)

        # Copy from vertex 1 to 0, with skipping build in properties
        self.mesh.copy_all_properties(self.vhandle[1], self.vhandle[0])
        v_it = self.mesh.vertices()
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 0)
        self.assertEqual(self.mesh.point(vh)[1], 0)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 0)
        self.assertEqual(self.mesh.point(vh)[1], 1)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 1)
        self.assertEqual(self.mesh.point(vh)[1], 1)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 1)
        self.assertEqual(self.mesh.point(vh)[1], 0)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        v_it = self.mesh.vertices()
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 1)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 1)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 2)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 3)

        # Copy from vertex 2 to 3, including build in properties
        self.mesh.copy_all_properties(self.vhandle[2], self.vhandle[3], True)
        v_it = self.mesh.vertices()
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 0)
        self.assertEqual(self.mesh.point(vh)[1], 0)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 0)
        self.assertEqual(self.mesh.point(vh)[1], 1)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 1)
        self.assertEqual(self.mesh.point(vh)[1], 1)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        vh = v_it.next()
        self.assertEqual(self.mesh.point(vh)[0], 1)
        self.assertEqual(self.mesh.point(vh)[1], 1)
        self.assertEqual(self.mesh.point(vh)[2], 0)
        v_it = self.mesh.vertices()
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 1)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 1)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 2)
        self.assertEqual(self.mesh.property(intHandle, v_it.next()), 2)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Property)
    unittest.TextTestRunner(verbosity=2).run(suite)
