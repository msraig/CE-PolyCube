#include <gtest/gtest.h>
#include <Unittests/unittests_common.hh>
#include <iostream>

namespace {

class OpenMeshVectorTest : public testing::Test {

    protected:

        // This function is called before each test is run
        virtual void SetUp() {
            
            // Do some initial stuff with the member data here...
        }

        // This function is called after all tests are through
        virtual void TearDown() {

            // Do some final stuff with the member data here...
        }

};



/*
 * ====================================================================
 * Define tests below
 * ====================================================================
 */

/* Compute surface area via cross product
 */
TEST_F(OpenMeshVectorTest, ComputeTriangleSurfaceWithCrossProduct) {


  //
  // vec1
  //  x
  //  |
  //  |
  //  |
  //  x------>x vec2
  //

  OpenMesh::Vec3d vec1(0.0,1.0,0.0);
  OpenMesh::Vec3d vec2(1.0,0.0,0.0);

  double area = 0.5 * cross(vec1,vec2).norm();
  EXPECT_EQ(0.5f , area ) << "Wrong area in cross product function";

  area = 0.5 * ( vec1 % vec2 ).norm();
  EXPECT_EQ(0.5f , area ) << "Wrong area in cross product operator";

}

/* Check OpenMesh Vector type abs function
 */
TEST_F(OpenMeshVectorTest, AbsTest) {

  OpenMesh::Vec3d vec1(0.5,0.5,-0.5);

  EXPECT_EQ( vec1.l8_norm() , 0.5f ) << "Wrong l8norm computation";

}

/* Compute surface area via cross product
 */
TEST_F(OpenMeshVectorTest, VectorCasting) {

  OpenMesh::Vec3d vecd(1.0,2.0,3.0);
  OpenMesh::Vec3f vecf = OpenMesh::vector_cast<OpenMesh::Vec3f>(vecd);
  EXPECT_EQ(1.f, vecf[0]) << "vector type cast failed on component 0";
  EXPECT_EQ(2.f, vecf[1]) << "vector type cast failed on component 1";
  EXPECT_EQ(3.f, vecf[2]) << "vector type cast failed on component 2";

  OpenMesh::Vec4d vecd4(40.0,30.0,20.0,10.0);
  vecd = OpenMesh::vector_cast<OpenMesh::Vec3d>(vecd4);
  EXPECT_EQ(40.0, vecd[0]) << "vector dimension cast failed on component 0";
  EXPECT_EQ(30.0, vecd[1]) << "vector dimension cast failed on component 1";
  EXPECT_EQ(20.0, vecd[2]) << "vector dimension cast failed on component 2";

}

#if __cplusplus > 199711L || defined(__GXX_EXPERIMENTAL_CXX0X__)
TEST_F(OpenMeshVectorTest, cpp11_constructors) {
    OpenMesh::Vec3d vec1 { 1.2, 2.0, 3.0 };

    EXPECT_EQ(1.2, vec1[0]);
    EXPECT_EQ(2.0, vec1[1]);
    EXPECT_EQ(3.0, vec1[2]);

    OpenMesh::Vec4f vec2 { 1.2f, 3.5f, 1.0f, 0.0f };

    EXPECT_EQ(1.2f, vec2[0]);
    EXPECT_EQ(3.5f, vec2[1]);
    EXPECT_EQ(1.0f, vec2[2]);
    EXPECT_EQ(0.0f, vec2[3]);

    OpenMesh::Vec4f vec2b { vec2 };

    EXPECT_EQ(1.2f, vec2b[0]);
    EXPECT_EQ(3.5f, vec2b[1]);
    EXPECT_EQ(1.0f, vec2b[2]);
    EXPECT_EQ(0.0f, vec2b[3]);

    OpenMesh::Vec4d vec4d { 1.23 };
    EXPECT_EQ(1.23, vec4d[0]);
    EXPECT_EQ(1.23, vec4d[1]);
    EXPECT_EQ(1.23, vec4d[2]);
    EXPECT_EQ(1.23, vec4d[3]);
}

TEST_F(OpenMeshVectorTest, cpp11_htmlColorLiteral) {
    static constexpr OpenMesh::Vec4f rose = 0xFFC7F1FF_htmlColor;

    const OpenMesh::Vec4f light_blue = 0x1FCFFFFF_htmlColor;
    EXPECT_LE((OpenMesh::Vec4f(0.1215686274f, 0.8117647058f, 1.0f, 1.0f)
        - light_blue).sqrnorm(), 1e-10);

    const auto light_blue_2 = 0x1FCFFFFF_htmlColor;
    // Check whether auto type deduction works as expected.
    static_assert(std::is_same<decltype(light_blue_2), decltype(light_blue)>
        ::value, "Bad type deduced from _htmlColor literal.");
    EXPECT_EQ(light_blue, light_blue_2);
}
#endif

}
