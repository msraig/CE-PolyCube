import unittest
import openmesh

class Collapse(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()
        self.vhandle = []
    
    def test_collapse_quad_with_center(self):
        
        # 0--------1
        # |\      /|
        # | \    / |
        # |  \  /  |
        # |    2   |
        # |  /  \  |
        # | /    \ |
        # 3--------4

        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 2, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(2, 2, 0)))

        # Add four faces
        face_vhandles = []

        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[2])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[4])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[4])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        self.mesh.add_face(face_vhandles)

        self.mesh.request_vertex_status()
        self.mesh.request_edge_status()
        self.mesh.request_face_status()

        # Get the halfedge
        v2v1 = self.mesh.find_halfedge(self.vhandle[2], self.vhandle[1])

        self.assertTrue(v2v1.is_valid())
        self.assertTrue(self.mesh.is_collapse_ok(v2v1))

        # Execute it as a crash test
        self.mesh.collapse(v2v1)

    def test_collapse_tetrahedron_complex(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(0, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(1, 0, 0)))
        
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
        
        self.mesh.request_vertex_status()
        self.mesh.request_edge_status()
        self.mesh.request_face_status()
        
        v0v1 = self.mesh.halfedge_handle(0)
        v1v0 = self.mesh.opposite_halfedge_handle(v0v1)

        v1vL = self.mesh.next_halfedge_handle(v0v1)
        vLv1 = self.mesh.opposite_halfedge_handle(v1vL)
        vLv0 = self.mesh.next_halfedge_handle(v1vL)
        v0vL = self.mesh.opposite_halfedge_handle(vLv0)

        vLvR = self.mesh.next_halfedge_handle(v0vL)
        vRvL = self.mesh.opposite_halfedge_handle(vLvR)

        v0vR = self.mesh.next_halfedge_handle(v1v0)
        vRv0 = self.mesh.opposite_halfedge_handle(v0vR)
        vRv1 = self.mesh.next_halfedge_handle(v0vR)
        v1vR = self.mesh.opposite_halfedge_handle(vRv1)

        v0 = self.mesh.from_vertex_handle(v0v1)
        v1 = self.mesh.to_vertex_handle(v0v1)
        vL = self.mesh.to_vertex_handle(self.mesh.next_halfedge_handle(v0v1))
        vR = self.mesh.to_vertex_handle(self.mesh.next_halfedge_handle(v1v0))

        # ===================================================================
        # Check preconditions
        # ===================================================================

        self.assertTrue(self.mesh.is_collapse_ok(v0v1))
        self.assertTrue(self.mesh.is_collapse_ok(v1v0))

        # Test the Vertex indices
        self.assertEqual(v0.idx(), 0)
        self.assertEqual(v1.idx(), 1)
        self.assertEqual(vL.idx(), 2)
        self.assertEqual(vR.idx(), 3)

        # Check the halfedges
        self.assertEqual(v0v1.idx(), 0)
        self.assertEqual(v1v0.idx(), 1)
        
        self.assertEqual(v1vL.idx(), 2)
        self.assertEqual(vLv1.idx(), 3)
        self.assertEqual(vLv0.idx(), 4)
        self.assertEqual(v0vL.idx(), 5)
        
        self.assertEqual(vLvR.idx(), 6)
        self.assertEqual(vRvL.idx(), 7)
        
        self.assertEqual(vRv0.idx(), 8)
        self.assertEqual(v0vR.idx(), 9)
        
        self.assertEqual(v1vR.idx(), 10)
        self.assertEqual(vRv1.idx(), 11)

        # ===================================================================
        # Execute collapse
        # ===================================================================

        self.mesh.collapse(v0v1)

        # ===================================================================
        # Check configuration afterwards
        # ===================================================================

        # Now the configuration should look like this:
        # The numbers at the side denote the halfedges
        #          1
        #         / \
        #        /   \
        #      //    \\
        #     3/2    11\10
        #     //        \\
        #    /    6-->   \
        #   2 ----------- 3
        #       <--7

        self.assertEqual(self.mesh.n_faces(), 4)

        # Check if the right vertices got deleted
        self.assertTrue(self.mesh.status(self.mesh.face_handle(0)).deleted())
        self.assertFalse(self.mesh.status(self.mesh.face_handle(1)).deleted())
        self.assertFalse(self.mesh.status(self.mesh.face_handle(2)).deleted())
        self.assertTrue(self.mesh.status(self.mesh.face_handle(3)).deleted())

        # Check the vertices of the two remaining faces
        fh_1 = self.mesh.face_handle(1)
        fh_2 = self.mesh.face_handle(2)

        fv_it = self.mesh.fv(fh_1)

        self.assertTrue(fv_it.__next__().idx(), 1)
        self.assertTrue(fv_it.__next__().idx(), 2)
        self.assertTrue(fv_it.__next__().idx(), 3)

        fv_it = self.mesh.fv(fh_2)

        self.assertTrue(fv_it.__next__().idx(), 2)
        self.assertTrue(fv_it.__next__().idx(), 1)
        self.assertTrue(fv_it.__next__().idx(), 3)

        # Get the first halfedge of face 1
        fh_1_he = self.mesh.halfedge_handle(fh_1)

        self.assertEqual(fh_1_he.idx(), 11)
        self.assertEqual(self.mesh.to_vertex_handle(fh_1_he).idx(), 1)

        next = self.mesh.next_halfedge_handle(fh_1_he)
        self.assertEqual(next.idx(), 2)
        self.assertEqual(self.mesh.to_vertex_handle(next).idx(), 2)

        next = self.mesh.next_halfedge_handle(next)
        self.assertEqual(next.idx(), 6)
        self.assertEqual(self.mesh.to_vertex_handle(next).idx(), 3)

        # Get the first halfedge of face 2
        fh_2_he = self.mesh.halfedge_handle(fh_2)

        self.assertEqual(fh_2_he.idx(), 7)
        self.assertEqual(self.mesh.to_vertex_handle(fh_2_he).idx(), 2)

        next = self.mesh.next_halfedge_handle(fh_2_he)
        self.assertEqual(next.idx(), 3)
        self.assertEqual(self.mesh.to_vertex_handle(next).idx(), 1)

        next = self.mesh.next_halfedge_handle(next)
        self.assertEqual(next.idx(), 10)
        self.assertEqual(self.mesh.to_vertex_handle(next).idx(), 3)

        # Vertex 1 outgoing
        voh_it = self.mesh.voh(self.mesh.vertex_handle(1))
        self.assertEqual(voh_it.__next__().idx(), 10)
        self.assertEqual(voh_it.__next__().idx(), 2)
        self.assertRaises(StopIteration, voh_it.__next__)

        # Vertex 2 outgoing
        voh_it = self.mesh.voh(self.mesh.vertex_handle(2))
        self.assertEqual(voh_it.__next__().idx(), 3)
        self.assertEqual(voh_it.__next__().idx(), 6)
        self.assertRaises(StopIteration, voh_it.__next__)

        # Vertex 2 outgoing
        voh_it = self.mesh.voh(self.mesh.vertex_handle(3))
        self.assertEqual(voh_it.__next__().idx(), 11)
        self.assertEqual(voh_it.__next__().idx(), 7)
        self.assertRaises(StopIteration, voh_it.__next__)

        # ===================================================================
        # Cleanup
        # ===================================================================
        self.mesh.garbage_collection()

        # ===================================================================
        # Check configuration afterwards
        # ===================================================================

        # Now the configuration should look like this:
        # The numbers at the side denote the halfedges
        #          0
        #         / \
        #        /   \
        #      //    \\
        #     4/5     0\1
        #     //        \\
        #    /    3-->   \
        #   2 ----------- 1
        #       <--2

        self.assertEqual(self.mesh.n_faces(), 2)

        # Check the vertices of the two remaining faces
        fh_0 = self.mesh.face_handle(0)
        fh_1 = self.mesh.face_handle(1)

        fv_it = self.mesh.fv(fh_0)

        self.assertEqual(fv_it.__next__().idx(), 2)
        self.assertEqual(fv_it.__next__().idx(), 1)
        self.assertEqual(fv_it.__next__().idx(), 0)

        fv_it = self.mesh.fv(fh_1)

        self.assertEqual(fv_it.__next__().idx(), 1)
        self.assertEqual(fv_it.__next__().idx(), 2)
        self.assertEqual(fv_it.__next__().idx(), 0)

        # Get the first halfedge of face 1
        fh_0_he = self.mesh.halfedge_handle(fh_0)

        self.assertEqual(fh_0_he.idx(), 5)
        self.assertEqual(self.mesh.to_vertex_handle(fh_0_he).idx(), 2)

        next = self.mesh.next_halfedge_handle(fh_0_he)
        self.assertEqual(next.idx(), 3)
        self.assertEqual(self.mesh.to_vertex_handle(next).idx(), 1)

        next = self.mesh.next_halfedge_handle(next)
        self.assertEqual(next.idx(), 0)
        self.assertEqual(self.mesh.to_vertex_handle(next).idx(), 0)

        # Get the first halfedge of face 1
        fh_1_he = self.mesh.halfedge_handle(fh_1)

        self.assertEqual(fh_1_he.idx(), 1)
        self.assertEqual(self.mesh.to_vertex_handle(fh_1_he).idx(), 1)

        next = self.mesh.next_halfedge_handle(fh_1_he)
        self.assertEqual(next.idx(), 2)
        self.assertEqual(self.mesh.to_vertex_handle(next).idx(), 2)

        next = self.mesh.next_halfedge_handle(next)
        self.assertEqual(next.idx(), 4)
        self.assertEqual(self.mesh.to_vertex_handle(next).idx(), 0)

        # Vertex 0 outgoing
        voh_it = self.mesh.voh(self.mesh.vertex_handle(0))
        self.assertEqual(voh_it.__next__().idx(), 1)
        self.assertEqual(voh_it.__next__().idx(), 5)
        self.assertRaises(StopIteration, voh_it.__next__)

        # Vertex 1 outgoing
        voh_it = self.mesh.voh(self.mesh.vertex_handle(1))
        self.assertEqual(voh_it.__next__().idx(), 0)
        self.assertEqual(voh_it.__next__().idx(), 2)
        self.assertRaises(StopIteration, voh_it.__next__)

        # Vertex 2 outgoing
        voh_it = self.mesh.voh(self.mesh.vertex_handle(2))
        self.assertEqual(voh_it.__next__().idx(), 3)
        self.assertEqual(voh_it.__next__().idx(), 4)
        self.assertRaises(StopIteration, voh_it.__next__)

        self.assertFalse(self.mesh.is_collapse_ok(self.mesh.halfedge_handle(0)))
        self.assertFalse(self.mesh.is_collapse_ok(self.mesh.halfedge_handle(1)))
        self.assertFalse(self.mesh.is_collapse_ok(self.mesh.halfedge_handle(2)))
        self.assertFalse(self.mesh.is_collapse_ok(self.mesh.halfedge_handle(3)))
        self.assertFalse(self.mesh.is_collapse_ok(self.mesh.halfedge_handle(4)))
        self.assertFalse(self.mesh.is_collapse_ok(self.mesh.halfedge_handle(5)))

    def test_collapse_tetrahedron(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 0,  0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1,  0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 0, -1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 0,  1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1,  0, 0)))
        
        # Add six faces
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[4])
        face_vhandles.append(self.vhandle[2])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[4])
        face_vhandles.append(self.vhandle[0])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[4])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[0])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)

        self.mesh.request_vertex_status()
        self.mesh.request_edge_status()
        self.mesh.request_face_status()

        # =============================================
        # Collapse halfedge from 0 to 4
        # =============================================

        heh_collapse1 = self.mesh.halfedge_handle(0)

        self.assertEqual(self.mesh.to_vertex_handle(heh_collapse1).idx(), 4)
        self.assertEqual(self.mesh.from_vertex_handle(heh_collapse1).idx(), 0)

        self.assertTrue(self.mesh.is_collapse_ok(heh_collapse1))
        self.mesh.collapse(heh_collapse1)

        heh_collapse2 = self.mesh.halfedge_handle(2)

        self.assertEqual(self.mesh.to_vertex_handle(heh_collapse2).idx(), 2)
        self.assertEqual(self.mesh.from_vertex_handle(heh_collapse2).idx(), 4)

        self.assertTrue(self.mesh.is_collapse_ok(heh_collapse2))
        self.mesh.collapse(heh_collapse2)

        heh_collapse3 = self.mesh.halfedge_handle(6)

        self.assertEqual(self.mesh.to_vertex_handle(heh_collapse3).idx(), 2)
        self.assertEqual(self.mesh.from_vertex_handle(heh_collapse3).idx(), 3)

        self.assertFalse(self.mesh.is_collapse_ok(heh_collapse3))

    def test_large_collapse_halfedge(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 0,  1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 1,  0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 2,  1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 0, -1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 2, -1, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d(-1,  0, 0)))
        self.vhandle.append(self.mesh.add_vertex(openmesh.Vec3d( 3,  0, 0)))
        
        # Add six faces
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[5])
        face_vhandles.append(self.vhandle[1])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[5])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)

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

        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[4])
        face_vhandles.append(self.vhandle[6])
        self.mesh.add_face(face_vhandles)

        # Test setup:
        #    0 ==== 2
        #   / \    /|\
        #  /   \  / | \
        # 5 --- 1   |  6
        #  \   /  \ | /
        #   \ /    \|/
        #    3 ==== 4

        # Request the status bits
        self.mesh.request_vertex_status()
        self.mesh.request_edge_status()
        self.mesh.request_face_status()

        # =============================================
        # Collapse halfedge from 1 to 4
        # =============================================

        heh_collapse = openmesh.HalfedgeHandle()

        for he in self.mesh.halfedges():
            if self.mesh.from_vertex_handle(he).idx() == 1 and self.mesh.to_vertex_handle(he).idx() == 4:
                heh_collapse = he

        # Check our halfedge
        self.assertEqual(self.mesh.to_vertex_handle(heh_collapse).idx(), 4)
        self.assertEqual(self.mesh.from_vertex_handle(heh_collapse).idx(), 1)
        self.assertTrue(self.mesh.is_collapse_ok(heh_collapse))

        # Remember the end vertices
        vh_from = self.mesh.from_vertex_handle(heh_collapse)
        vh_to = self.mesh.to_vertex_handle(heh_collapse)

        # Collapse it
        self.mesh.collapse(heh_collapse)

        self.assertTrue(self.mesh.status(vh_from).deleted())
        self.assertFalse(self.mesh.status(vh_to).deleted())


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Collapse)
    unittest.TextTestRunner(verbosity=2).run(suite)
