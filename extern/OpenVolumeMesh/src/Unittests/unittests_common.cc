/*
 * unittests_common.cc
 *
 *  Created on: Jan 10, 2012
 *      Author: kremer
 */

#include "unittests_common.hh"

void PolyhedralMeshBase::generatePolyhedralMesh(PolyhedralMesh& _mesh) {

    Vec3d p1(0.0, 0.0, 0.0);
    Vec3d p2(1.0, 0.0, 0.0);
    Vec3d p3(1.0, 1.0, 0.0);
    Vec3d p4(0.0, 1.0, 0.0);

    Vec3d p5(0.0, 0.0, 1.0);
    Vec3d p6(1.0, 0.0, 1.0);
    Vec3d p7(1.0, 1.0, 1.0);
    Vec3d p8(0.0, 1.0, 1.0);

    Vec3d p9(0.0, 0.0, 2.0);
    Vec3d p10(1.0, 0.0, 2.0);
    Vec3d p11(1.0, 1.0, 2.0);
    Vec3d p12(0.0, 1.0, 2.0);

    VertexHandle v1 = _mesh.add_vertex(p1);
    VertexHandle v2 = _mesh.add_vertex(p2);
    VertexHandle v3 = _mesh.add_vertex(p3);
    VertexHandle v4 = _mesh.add_vertex(p4);

    VertexHandle v5 = _mesh.add_vertex(p5);
    VertexHandle v6 = _mesh.add_vertex(p6);
    VertexHandle v7 = _mesh.add_vertex(p7);
    VertexHandle v8 = _mesh.add_vertex(p8);

    VertexHandle v9  = _mesh.add_vertex(p9);
    VertexHandle v10 = _mesh.add_vertex(p10);
    VertexHandle v11 = _mesh.add_vertex(p11);
    VertexHandle v12 = _mesh.add_vertex(p12);

    std::vector<VertexHandle> vertices;
    vertices.push_back(v1); vertices.push_back(v2);
    vertices.push_back(v3); vertices.push_back(v4);
    FaceHandle f1 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v5); vertices.push_back(v6);
    vertices.push_back(v7); vertices.push_back(v8);
    FaceHandle f2 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v2); vertices.push_back(v6);
    vertices.push_back(v7); vertices.push_back(v3);
    FaceHandle f3 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v1); vertices.push_back(v4);
    vertices.push_back(v8); vertices.push_back(v5);
    FaceHandle f4 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v4); vertices.push_back(v3);
    vertices.push_back(v7); vertices.push_back(v8);
    FaceHandle f5 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v6); vertices.push_back(v5);
    vertices.push_back(v1); vertices.push_back(v2);
    FaceHandle f6 = _mesh.add_face(vertices);

    // Add first cell
    std::vector<HalfFaceHandle> halffaces;
    halffaces.push_back(_mesh.halfface_handle(f1, 1)); halffaces.push_back(_mesh.halfface_handle(f2, 0));
    halffaces.push_back(_mesh.halfface_handle(f3, 1)); halffaces.push_back(_mesh.halfface_handle(f4, 1));
    halffaces.push_back(_mesh.halfface_handle(f5, 1)); halffaces.push_back(_mesh.halfface_handle(f6, 0));
    _mesh.add_cell(halffaces);

    vertices.clear();
    vertices.push_back(v9); vertices.push_back(v10);
    vertices.push_back(v11); vertices.push_back(v12);
    FaceHandle f7 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v10); vertices.push_back(v11);
    vertices.push_back(v7); vertices.push_back(v6);
    FaceHandle f8 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v5); vertices.push_back(v8);
    vertices.push_back(v12); vertices.push_back(v9);
    FaceHandle f9 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v8); vertices.push_back(v7);
    vertices.push_back(v11); vertices.push_back(v12);
    FaceHandle f10 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v10); vertices.push_back(v9);
    vertices.push_back(v5); vertices.push_back(v6);
    FaceHandle f11 = _mesh.add_face(vertices);

    halffaces.clear();
    halffaces.push_back(_mesh.halfface_handle(f2,  1)); halffaces.push_back(_mesh.halfface_handle(f7,  0));
    halffaces.push_back(_mesh.halfface_handle(f8,  1)); halffaces.push_back(_mesh.halfface_handle(f9,  1));
    halffaces.push_back(_mesh.halfface_handle(f10, 1)); halffaces.push_back(_mesh.halfface_handle(f11, 0));
    _mesh.add_cell(halffaces);
}

void HexahedralMeshBase::generateHexahedralMesh(HexahedralMesh& _mesh) {

    Vec3d p1(0.0, 0.0, 0.0);
    Vec3d p2(1.0, 0.0, 0.0);
    Vec3d p3(1.0, 1.0, 0.0);
    Vec3d p4(0.0, 1.0, 0.0);

    Vec3d p5(0.0, 0.0, 1.0);
    Vec3d p6(1.0, 0.0, 1.0);
    Vec3d p7(1.0, 1.0, 1.0);
    Vec3d p8(0.0, 1.0, 1.0);

    Vec3d p9(0.0, 0.0, 2.0);
    Vec3d p10(1.0, 0.0, 2.0);
    Vec3d p11(1.0, 1.0, 2.0);
    Vec3d p12(0.0, 1.0, 2.0);

    VertexHandle v1 = _mesh.add_vertex(p1);
    VertexHandle v2 = _mesh.add_vertex(p2);
    VertexHandle v3 = _mesh.add_vertex(p3);
    VertexHandle v4 = _mesh.add_vertex(p4);

    VertexHandle v5 = _mesh.add_vertex(p5);
    VertexHandle v6 = _mesh.add_vertex(p6);
    VertexHandle v7 = _mesh.add_vertex(p7);
    VertexHandle v8 = _mesh.add_vertex(p8);

    VertexHandle v9  = _mesh.add_vertex(p9);
    VertexHandle v10 = _mesh.add_vertex(p10);
    VertexHandle v11 = _mesh.add_vertex(p11);
    VertexHandle v12 = _mesh.add_vertex(p12);

    std::vector<VertexHandle> vertices;
    vertices.push_back(v1); vertices.push_back(v2);
    vertices.push_back(v3); vertices.push_back(v4);
    FaceHandle f0 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v5); vertices.push_back(v6);
    vertices.push_back(v7); vertices.push_back(v8);
    FaceHandle f1 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v2); vertices.push_back(v6);
    vertices.push_back(v7); vertices.push_back(v3);
    FaceHandle f2 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v1); vertices.push_back(v5);
    vertices.push_back(v8); vertices.push_back(v4);
    FaceHandle f3 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v1); vertices.push_back(v2);
    vertices.push_back(v6); vertices.push_back(v5);
    FaceHandle f4 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v4); vertices.push_back(v3);
    vertices.push_back(v7); vertices.push_back(v8);
    FaceHandle f5 = _mesh.add_face(vertices);

    // Add first cell
    std::vector<HalfFaceHandle> halffaces;
    halffaces.push_back(_mesh.halfface_handle(f0, 1)); halffaces.push_back(_mesh.halfface_handle(f1, 0));
    halffaces.push_back(_mesh.halfface_handle(f2, 1)); halffaces.push_back(_mesh.halfface_handle(f3, 0));
    halffaces.push_back(_mesh.halfface_handle(f4, 0)); halffaces.push_back(_mesh.halfface_handle(f5, 1));
    _mesh.add_cell(halffaces);

    vertices.clear();
    vertices.push_back(v9); vertices.push_back(v10);
    vertices.push_back(v11); vertices.push_back(v12);
    FaceHandle f6 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v6); vertices.push_back(v10);
    vertices.push_back(v11); vertices.push_back(v7);
    FaceHandle f7 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v5); vertices.push_back(v9);
    vertices.push_back(v12); vertices.push_back(v8);
    FaceHandle f8 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v5); vertices.push_back(v6);
    vertices.push_back(v10); vertices.push_back(v9);
    FaceHandle f9 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v8); vertices.push_back(v7);
    vertices.push_back(v11); vertices.push_back(v12);
    FaceHandle f10 = _mesh.add_face(vertices);

    halffaces.clear();
    halffaces.push_back(_mesh.halfface_handle(f1, 1)); halffaces.push_back(_mesh.halfface_handle(f6, 0));
    halffaces.push_back(_mesh.halfface_handle(f7, 1)); halffaces.push_back(_mesh.halfface_handle(f8, 0));
    halffaces.push_back(_mesh.halfface_handle(f9, 0)); halffaces.push_back(_mesh.halfface_handle(f10, 1));
    _mesh.add_cell(halffaces);
}
