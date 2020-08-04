import unittest
import openmesh

class VectorTest(unittest.TestCase):

    def test_compute_triangle_surface_with_cross_product(self):
        # vec1
        #  y
        #  |
        #  |
        #  |
        #  x------>x vec2
        
        vec1 = openmesh.Vec3d(0.0, 1.0, 0.0)
        vec2 = openmesh.Vec3d(1.0, 0.0, 0.0)

        area = 0.5 * openmesh.cross(vec1, vec2).norm()
        self.assertEqual(area, 0.5)

        area = 0.5 * (vec1 % vec2).norm()
        self.assertEqual(area, 0.5)

    def test_equality_operator_vec3d(self):
        vec1 = openmesh.Vec3d(0.0, 1.0, 0.0)
        vec2 = openmesh.Vec3d(1.0, 0.0, 0.0)
        vec3 = openmesh.Vec3d(1.0, 0.0, 0.0)
        
        self.assertFalse(vec1==vec2)
        self.assertTrue(vec3==vec2)
        
    def test_equality_operator_vec3f(self):
        vec1 = openmesh.Vec3f(0.0, 1.0, 0.0)
        vec2 = openmesh.Vec3f(1.0, 0.0, 0.0)
        vec3 = openmesh.Vec3f(1.0, 0.0, 0.0)
        
        self.assertFalse(vec1==vec2)
        self.assertTrue(vec3==vec2)     

    def test_abs_test(self):
        vec1 = openmesh.Vec3d(0.5, 0.5, -0.5)
        self.assertEqual(vec1.l8_norm(), 0.5)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(VectorTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
