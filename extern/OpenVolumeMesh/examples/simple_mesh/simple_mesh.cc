// C++ includes
#include <iostream>
#include <vector>

// Include vector classes
#include <OpenVolumeMesh/Geometry/VectorT.hh>

// Include polyhedral mesh kernel
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

// Make some typedefs to facilitate your life
typedef OpenVolumeMesh::Geometry::Vec3f         Vec3f;
typedef OpenVolumeMesh::GeometryKernel<Vec3f>   PolyhedralMeshV3f;

int main(int _argc, char** _argv) {

    // Create mesh object
    PolyhedralMeshV3f myMesh;

    // Add eight vertices
    OpenVolumeMesh::VertexHandle v0 = myMesh.add_vertex(Vec3f(-1.0, 0.0, 0.0));
    OpenVolumeMesh::VertexHandle v1 = myMesh.add_vertex(Vec3f( 0.0, 0.0, 1.0));
    OpenVolumeMesh::VertexHandle v2 = myMesh.add_vertex(Vec3f( 1.0, 0.0, 0.0));
    OpenVolumeMesh::VertexHandle v3 = myMesh.add_vertex(Vec3f( 0.0, 0.0,-1.0));
    OpenVolumeMesh::VertexHandle v4 = myMesh.add_vertex(Vec3f( 0.0, 1.0, 0.0));
    
    std::vector<OpenVolumeMesh::VertexHandle> vertices;
    
    // Add faces
    vertices.push_back(v0); vertices.push_back(v1);vertices.push_back(v4);
    OpenVolumeMesh::FaceHandle f0 = myMesh.add_face(vertices);
    
    vertices.clear();
    vertices.push_back(v1); vertices.push_back(v2);vertices.push_back(v4);
    OpenVolumeMesh::FaceHandle f1 = myMesh.add_face(vertices);
    
    vertices.clear();
    vertices.push_back(v0); vertices.push_back(v1);vertices.push_back(v2);
    OpenVolumeMesh::FaceHandle f2 = myMesh.add_face(vertices);
    
    vertices.clear();
    vertices.push_back(v0); vertices.push_back(v4);vertices.push_back(v2);
    OpenVolumeMesh::FaceHandle f3 = myMesh.add_face(vertices);
    
    vertices.clear();
    vertices.push_back(v0); vertices.push_back(v4);vertices.push_back(v3);
    OpenVolumeMesh::FaceHandle f4 = myMesh.add_face(vertices);
    
    vertices.clear();
    vertices.push_back(v2); vertices.push_back(v3);vertices.push_back(v4);
    OpenVolumeMesh::FaceHandle f5 = myMesh.add_face(vertices);
    
    vertices.clear();
    vertices.push_back(v0); vertices.push_back(v2);vertices.push_back(v3);
    OpenVolumeMesh::FaceHandle f6 = myMesh.add_face(vertices);
    
    std::vector<OpenVolumeMesh::HalfFaceHandle> halffaces;
    
    // Add first tetrahedron
    halffaces.push_back(myMesh.halfface_handle(f0, 1));
    halffaces.push_back(myMesh.halfface_handle(f1, 1));
    halffaces.push_back(myMesh.halfface_handle(f2, 0)); 
    halffaces.push_back(myMesh.halfface_handle(f3, 1)); 
    myMesh.add_cell(halffaces);
    
    // Add second tetrahedron
    halffaces.clear();
    halffaces.push_back(myMesh.halfface_handle(f4, 1));
    halffaces.push_back(myMesh.halfface_handle(f5, 1));
    halffaces.push_back(myMesh.halfface_handle(f3, 0)); 
    halffaces.push_back(myMesh.halfface_handle(f6, 0)); 
    myMesh.add_cell(halffaces);

    // Print positions of vertices to std out
    for(OpenVolumeMesh::VertexIter v_it = myMesh.vertices_begin();
            v_it != myMesh.vertices_end(); ++v_it) {

        std::cout << "Position of vertex " << v_it->idx() << ": " <<
            myMesh.vertex(*v_it) << std::endl;
    }

    return 0;
}
