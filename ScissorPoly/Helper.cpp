#pragma warning( disable : 4477 4018 4267 4244 4838)
#include <iostream>
#include "omp.h"
#include "Helper.h"
#include "SmallMat.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <list>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/Status.hh>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/shared_array_property_map.hpp>
#include <boost/graph/smallest_last_ordering.hpp>
#include <boost/graph/sequential_vertex_coloring.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#define epsilon 0.0001
#define PI 3.1415926536
#ifndef rd2
#define rd2 (1.414213562373095633972752693807)
#endif
#ifndef rd3
#define rd3 (1.732050807568877637265813973499)
#endif
#ifndef rd6
#define rd6 (2.449489742783178098197284074706)
#endif
using std::vector;
namespace ig
{
    namespace IVF
    {
        typedef Eigen::SparseMatrix<double> SpMat;
        typedef Eigen::Triplet<double> T;
        using Eigen::VectorXd;
        using Eigen::Matrix3d;
        void IVF_svd(const std::vector<Vector3d> &pts_ref, const std::vector<unsigned int> &elems, std::vector<Vector3d> &pts_flip, const std::vector<int> &boundary_tag, bool using_standard_tet, double min_singular_value, int max_iter)
        {
            int nv = pts_ref.size();
            int nc = elems.size() / 4;
            std::vector<Matrix3d> cell_S;
            std::vector<double> cell_volume;
            cell_S.resize(nc);
            cell_volume.resize(nc);
            Matrix3d VS, IS, diag, Q, U, V, R;
            diag << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
            double avg_edge_length = 0.0;
            if (!using_standard_tet)
                for (size_t c_id = 0; c_id < nc; c_id++)
                {
                    int id0, id1, id2, id3;
                    id0 = elems[4 * c_id];
                    id1 = elems[4 * c_id + 1];
                    id2 = elems[4 * c_id + 2];
                    id3 = elems[4 * c_id + 3];
                    CVec<double, 3> cp;
                    CVec<double, 3> cr;
                    CVec<double, 3> cs;
                    CVec<double, 3> ct;
                    cp = pts_ref[id0];
                    cr = pts_ref[id1];
                    cs = pts_ref[id2];
                    ct = pts_ref[id3];
                    cr = cr - cp;
                    cs = cs - cp;
                    ct = ct - cp;
                    VS(0, 0) = cr[0]; VS(1, 0) = cr[1]; VS(2, 0) = cr[2];
                    VS(0, 1) = cs[0]; VS(1, 1) = cs[1]; VS(2, 1) = cs[2];
                    VS(0, 2) = ct[0]; VS(1, 2) = ct[1]; VS(2, 2) = ct[2];
                    IS = VS.inverse();
                    cell_S[c_id] = IS;
                    cell_volume[c_id] = -VS.determinant();
                }
            else
            {
                VS(0, 0) = 1.0; VS(1, 0) = 0.0; VS(2, 0) = 0.0;
                VS(0, 1) = -1.0 / rd3; VS(1, 1) = 2.0 / rd3; VS(2, 1) = 0.0;
                VS(0, 2) = -1.0 / rd6; VS(1, 2) = -1.0 / rd6; VS(2, 2) = rd3 / rd2;
                for (size_t c_id = 0; c_id < nc; c_id++)
                {
                    cell_S[c_id] = VS;
                    cell_volume[c_id] = 1.0;
                }
                for (size_t c_id = 0; c_id < nc; c_id++)
                {
                    int id0, id1, id2, id3;
                    id0 = elems[4 * c_id];
                    id1 = elems[4 * c_id + 1];
                    id2 = elems[4 * c_id + 2];
                    id3 = elems[4 * c_id + 3];
                    CVec<double, 3> cp;
                    CVec<double, 3> cr;
                    CVec<double, 3> cs;
                    CVec<double, 3> ct;
                    cp = pts_flip[id0];
                    cr = pts_flip[id1];
                    cs = pts_flip[id2];
                    ct = pts_flip[id3];
                    avg_edge_length += (cp - cr).L2Norm();
                    avg_edge_length += (cp - cs).L2Norm();
                    avg_edge_length += (cp - ct).L2Norm();
                    avg_edge_length += (cr - cs).L2Norm();
                    avg_edge_length += (cr - ct).L2Norm();
                    avg_edge_length += (cs - ct).L2Norm();
                }
                avg_edge_length = avg_edge_length / (6 * nc);
                std::cout << "avg edge length: " << avg_edge_length << std::endl;
                for (size_t i = 0; i < pts_flip.size(); i++)
                {
                    pts_flip[i] = pts_flip[i] / avg_edge_length;
                }
            }
            std::vector < std::vector < double >> Cell_R(nc);
            for (int c_id = 0; c_id < nc; ++c_id)
            {
                Cell_R[c_id].resize(9);
            }
            int var_count = 0; std::vector<int> map_v(nv, -1);
            for (int j = 0; j < nv; ++j)
            {
                if (boundary_tag[j] != 1)
                {
                    map_v[j] = var_count;
                    ++var_count;
                }
            }
            std::cout << "var number: " << var_count << std::endl;
            std::vector<T> tripletList;
            SpMat A_eigen(nc * 3, var_count);
            std::vector<double> cs(4);
            for (int c_id = 0; c_id < nc; ++c_id)
            {
                double cv = 1.0;
                Eigen::Matrix3d& CS = cell_S[c_id];
                for (int j = 0; j < 3; ++j)
                {
                    cs[1] = CS(0, j)*cv; cs[2] = CS(1, j)*cv; cs[3] = CS(2, j)*cv;
                    cs[0] = -(cs[1] + cs[2] + cs[3]);
                    for (int k = 0; k < 4; ++k)
                    {
                        int cv_id = elems[4 * c_id + k];
                        int var_id = map_v[cv_id];
                        if (var_id >= 0)
                        {
                            tripletList.push_back(T(3 * c_id + j, var_id, cs[k]));
                        }
                    }
                }
            }
            A_eigen.setFromTriplets(tripletList.begin(), tripletList.end());
            SpMat ATA = A_eigen.transpose() * A_eigen;
            Eigen::SimplicialLLT<SpMat> eigen_solver;
            eigen_solver.compute(ATA);
            for (size_t iter = 0; iter < max_iter; iter++)
            {
                for (size_t c_id = 0; c_id < nc; c_id++)
                {
                    int id0, id1, id2, id3;
                    id0 = elems[4 * c_id];
                    id1 = elems[4 * c_id + 1];
                    id2 = elems[4 * c_id + 2];
                    id3 = elems[4 * c_id + 3];
                    CVec<double, 3> cp;
                    CVec<double, 3> cr;
                    CVec<double, 3> cs;
                    CVec<double, 3> ct;
                    cp = pts_flip[id0];
                    cr = pts_flip[id1];
                    cs = pts_flip[id2];
                    ct = pts_flip[id3];
                    cr = cr - cp;
                    cs = cs - cp;
                    ct = ct - cp;
                    VS(0, 0) = cr[0]; VS(1, 0) = cr[1]; VS(2, 0) = cr[2];
                    VS(0, 1) = cs[0]; VS(1, 1) = cs[1]; VS(2, 1) = cs[2];
                    VS(0, 2) = ct[0]; VS(1, 2) = ct[1]; VS(2, 2) = ct[2];
                    Q = VS * cell_S[c_id];
                    Eigen::JacobiSVD<Eigen::Matrix3d> svd(Q, Eigen::ComputeFullU | Eigen::ComputeFullV);
                    U = svd.matrixU(); V = svd.matrixV();
                    double sig0 = svd.singularValues()[0];
                    double sig1 = svd.singularValues()[1];
                    double sig2 = svd.singularValues()[2];
                    diag(0, 0) = sig0;
                    diag(1, 1) = sig1;
                    diag(2, 2) = sig2;
                    Eigen::Matrix3d mult = U * diag * V.transpose();
                    if (Q.determinant() > 0)
                    {
                        R = V * diag * U.transpose();
                    }
                    else
                    {
                        U(0, 2) = -U(0, 2); U(1, 2) = -U(1, 2); U(2, 2) = -U(2, 2);
                        diag(2, 2) = min_singular_value;
                        R = V * diag * U.transpose();
                    }
                    Cell_R[c_id][0] = R(0, 0); Cell_R[c_id][1] = R(0, 1); Cell_R[c_id][2] = R(0, 2);
                    Cell_R[c_id][3] = R(1, 0); Cell_R[c_id][4] = R(1, 1); Cell_R[c_id][5] = R(1, 2);
                    Cell_R[c_id][6] = R(2, 0); Cell_R[c_id][7] = R(2, 1); Cell_R[c_id][8] = R(2, 2);
                }
                for (int i = 0; i < 3; ++i)
                {
                    VectorXd b_eigen(nc * 3);
                    for (size_t j = 0; j < nc * 3; j++)
                    {
                        b_eigen(j) = 0.0;
                    }
                    for (int c_id = 0; c_id < nc; ++c_id)
                    {
                        double cv = 1.0;
                        Eigen::Matrix3d& CS = cell_S[c_id];
                        std::vector<double>& CR = Cell_R[c_id];
                        for (int j = 0; j < 3; ++j)
                        {
                            cs[1] = CS(0, j)*cv; cs[2] = CS(1, j)*cv; cs[3] = CS(2, j)*cv;
                            cs[0] = -(cs[1] + cs[2] + cs[3]);
                            for (int k = 0; k < 4; ++k)
                            {
                                int cv_id = elems[4 * c_id + k];
                                int var_id = map_v[cv_id];
                                if (var_id < 0)
                                {
                                    b_eigen(3 * c_id + j) += -cs[k] * pts_flip[cv_id][i];
                                }
                                else
                                {
                                }
                            }
                            b_eigen(3 * c_id + j) += cv * CR[i + 3 * j];
                        }
                    }
                    VectorXd ATb = A_eigen.transpose() * b_eigen;
                    VectorXd solution = eigen_solver.solve(ATb);
                    for (int j = 0; j < nv; ++j)
                    {
                        if (map_v[j] >= 0)
                        {
                            int var_id = map_v[j];
                            CVec<double, 3> np = pts_flip[j];
                            np[i] = solution(var_id);
                            pts_flip[j] = np;
                        }
                    }
                }
            }
            if (using_standard_tet)
            {
                for (size_t i = 0; i < pts_flip.size(); i++)
                {
                    pts_flip[i] = pts_flip[i] * avg_edge_length;
                }
            }
        }
    }
    namespace SimpleTriangulation
    {
        struct FaceInfo2
        {
            FaceInfo2() {}
            int nesting_level;
            bool in_domain() {
                return nesting_level % 2 == 1;
            }
        };
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef CGAL::Triangulation_vertex_base_2<K> Vb;
        typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K> Fbb;
        typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb> Fb;
        typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
        typedef CGAL::Exact_predicates_tag Itag;
        typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
        typedef CDT::Point Point;
        typedef CGAL::Polygon_2<K> Polygon_2;
        typedef Polygon_2::Vertex_iterator VertexIterator;
        void
            mark_domains(CDT& ct,
                CDT::Face_handle start,
                int index,
                std::list<CDT::Edge>& border)
        {
            if (start->info().nesting_level != -1) {
                return;
            }
            std::list<CDT::Face_handle> queue;
            queue.push_back(start);
            while (!queue.empty()) {
                CDT::Face_handle fh = queue.front();
                queue.pop_front();
                if (fh->info().nesting_level == -1) {
                    fh->info().nesting_level = index;
                    for (int i = 0; i < 3; i++) {
                        CDT::Edge e(fh, i);
                        CDT::Face_handle n = fh->neighbor(i);
                        if (n->info().nesting_level == -1) {
                            if (ct.is_constrained(e)) border.push_back(e);
                            else queue.push_back(n);
                        }
                    }
                }
            }
        }
        void
            mark_domains(CDT& cdt)
        {
            for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it) {
                it->info().nesting_level = -1;
            }
            std::list<CDT::Edge> border;
            mark_domains(cdt, cdt.infinite_face(), 0, border);
            while (!border.empty()) {
                CDT::Edge e = border.front();
                border.pop_front();
                CDT::Face_handle n = e.first->neighbor(e.second);
                if (n->info().nesting_level == -1) {
                    mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
                }
            }
        }
        bool triangulation(const std::vector<double>& coord_x, const std::vector<double>& coord_y, const std::vector<double>& coord_z, const std::vector<std::vector<int>> &input_faces, std::vector<std::vector<int>> &output_faces)
        {
            std::vector<double> x_2d, y_2d;
            double normal_x, normal_y, normal_z;
            ig::CGALHelper::CrossVectorNormalize(coord_x[0] - coord_x[1], coord_y[0] - coord_y[1], coord_z[0] - coord_z[1], coord_x[2] - coord_x[1], coord_y[2] - coord_y[1], coord_z[2] - coord_z[1], normal_x, normal_y, normal_z);
            ig::CGALHelper::localize3dvector(coord_x, coord_y, coord_z, normal_x, normal_y, normal_z, x_2d, y_2d);
            int vert_count = 0;
            std::vector<int> idx_array;
            CDT cdt;
            for (size_t i = 0; i < input_faces.size(); i++)
            {
                Polygon_2 temp_polygon;
                for (size_t j = 0; j < input_faces[i].size(); j++)
                {
                    temp_polygon.push_back(Point(x_2d[vert_count], y_2d[vert_count]));
                    vert_count++;
                    idx_array.push_back(input_faces[i][j]);
                }
                cdt.insert_constraint(temp_polygon.vertices_begin(), temp_polygon.vertices_end(), true);
            }
            mark_domains(cdt);
            std::map<Point, int> pointtoindex;
            int ind = 0;
            Point temp_point;
            for (CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin();
                vit != cdt.finite_vertices_end(); ++vit)
            {
                temp_point = vit->point();
                pointtoindex.insert(std::pair<Point, int>(temp_point, idx_array[ind++]));
            }
            if (cdt.number_of_vertices() != vert_count) return false;
            assert(cdt.number_of_vertices() == vert_count);
            output_faces.clear();
            std::vector<int> face;
            int count = 0;
            for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
                fit != cdt.finite_faces_end(); ++fit)
            {
                face.clear();
                if (fit->info().in_domain())
                {
                    ++count;
                    face.push_back(pointtoindex.at(fit->vertex(0)->point()));
                    face.push_back(pointtoindex.at(fit->vertex(1)->point()));
                    face.push_back(pointtoindex.at(fit->vertex(2)->point()));
                    output_faces.push_back(face);
                }
            }
            return true;
        }
        bool sortface_area(std::vector<double>& coord_x, std::vector<double>& coord_y, std::vector<double>& coord_z, std::vector<std::vector<int>> &input_faces)
        {
            if (input_faces.size() == 1)
                return false;
            std::vector<double> x_2d, y_2d;
            double normal_x, normal_y, normal_z;
            ig::CGALHelper::CrossVectorNormalize(coord_x[0] - coord_x[1], coord_y[0] - coord_y[1], coord_z[0] - coord_z[1], coord_x[2] - coord_x[1], coord_y[2] - coord_y[1], coord_z[2] - coord_z[1], normal_x, normal_y, normal_z);
            ig::CGALHelper::localize3dvector(coord_x, coord_y, coord_z, normal_x, normal_y, normal_z, x_2d, y_2d);
            std::vector<Polygon_2> polygon_array;
            int idx = 0;
            int max_idx = -1;
            double max_area = -1.f;
            std::vector<double> area_array;
            for (size_t i = 0; i < input_faces.size(); i++)
            {
                Polygon_2 temp_polygon;
                for (size_t j = 0; j < input_faces[i].size(); j++)
                {
                    temp_polygon.push_back(Point(x_2d[idx], y_2d[idx]));
                    idx++;
                }
                polygon_array.push_back(temp_polygon);
                double Area = abs(temp_polygon.area());
                if (Area > max_area)
                {
                    max_area = Area;
                    max_idx = i;
                }
                area_array.push_back(Area);
            }
            if (max_idx != -1 && max_idx != 0)
            {
                std::vector<int> temp_face = input_faces[0];
                input_faces[0] = input_faces[max_idx];
                input_faces[max_idx] = temp_face;
                return true;
            }
            return false;
        }
        typedef OpenMesh::PolyMesh_ArrayKernelT< OpenMesh::DefaultTraits> Mesh;
        bool merge_triangles(const std::vector<double>& coord_x, const std::vector<double>& coord_y, const std::vector<double>& coord_z, const std::vector<std::vector<int>> &input_faces, std::vector<std::vector<int>> &triangles)
        {
            if (input_faces.size() == 1)
                return false;
            std::vector<int> local_to_global;
            std::map<int, int> global_to_local;
            std::map<int, int> vert_to_boundary_map;
            std::vector<std::vector<int>> local_faces;
            std::vector<int> face;
            for (size_t i = 0; i < input_faces.size(); i++)
            {
                face.clear();
                for (size_t j = 0; j < input_faces[i].size(); j++)
                {
                    local_to_global.push_back(input_faces[i][j]);
                    global_to_local[input_faces[i][j]] = local_to_global.size() - 1;
                    vert_to_boundary_map[local_to_global.size() - 1] = i;
                    face.push_back(local_to_global.size() - 1);
                }
                local_faces.push_back(face);
            }
            Mesh mesh;
            mesh.request_face_status();
            mesh.request_edge_status();
            mesh.request_vertex_status();
            OpenMesh::VPropHandleT<int> tag;
            mesh.add_property(tag);
            std::vector<Mesh::VertexHandle> vhandle;
            std::vector<Mesh::FaceHandle> fhandle;
            vhandle.resize(coord_x.size());
            fhandle.resize(triangles.size());
            for (size_t i = 0; i < coord_x.size(); i++)
            {
                vhandle[i] = mesh.add_vertex(Mesh::Point(coord_x[i], coord_y[i], coord_z[i]));
                mesh.property(tag, vhandle[i]) = vert_to_boundary_map[i];
            }
            std::vector<Mesh::VertexHandle> tmp_face_vhandles;
            for (size_t i = 0; i < triangles.size(); i++)
            {
                tmp_face_vhandles.clear();
                for (size_t j = 0; j < triangles[i].size(); j++)
                {
                    int tmp_idx = global_to_local[triangles[i][j]];
                    tmp_face_vhandles.push_back(vhandle[tmp_idx]);
                }
                fhandle[i] = mesh.add_face(tmp_face_vhandles);
            }
            std::vector<std::vector<Mesh::HalfedgeHandle>> cand_hf;
            cand_hf.resize(input_faces.size() - 1);
            Mesh::HalfedgeIter he_it = mesh.halfedges_begin();
            Mesh::HalfedgeIter he_end = mesh.halfedges_end();
            for (; he_it != he_end; ++he_it)
            {
                int idx1 = mesh.from_vertex_handle((*he_it)).idx();
                int idx2 = mesh.to_vertex_handle((*he_it)).idx();
                if (vert_to_boundary_map[idx2] == 0 && vert_to_boundary_map[idx1] != 0)
                    cand_hf[vert_to_boundary_map[idx1] - 1].push_back((*he_it));
            }
            std::vector<Mesh::EdgeHandle> selected_eh;
            for (size_t i = 0; i < cand_hf.size(); i++)
            {
                if (cand_hf[i].size() < 2)
                {
                    std::cout << "Merge Error1 " << std::endl;
                    return false;
                }
                selected_eh.push_back(mesh.edge_handle(cand_hf[i][0]));
                int to_idx = mesh.to_vertex_handle(cand_hf[i][0]).idx();
                int from_idx = mesh.from_vertex_handle(cand_hf[i][0]).idx();
                int j = 1;
                for (; j < cand_hf[i].size(); j++)
                {
                    int temp_to_idx = mesh.to_vertex_handle(cand_hf[i][j]).idx();
                    int temp_from_idx = mesh.from_vertex_handle(cand_hf[i][j]).idx();
                    if (temp_to_idx != to_idx && temp_from_idx != from_idx)
                        break;
                }
                if (j == cand_hf[i].size())
                {
                    std::cout << "Merge Error2 " << std::endl;
                    return false;
                }
                selected_eh.push_back(mesh.edge_handle(cand_hf[i][j]));
            }
            std::cout << "selected edge: " << selected_eh.size() << std::endl;
            std::vector<Mesh::EdgeHandle> deleted_edge;
            for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
            {
                if (mesh.is_boundary(e_it))
                {
                    continue;
                }
                std::vector<Mesh::EdgeHandle>::iterator tmp_it = std::find(selected_eh.begin(), selected_eh.end(), (*e_it));
                if (tmp_it == selected_eh.end())
                {
                    deleted_edge.push_back((*e_it));
                }
            }
            std::cout << "deleted edge: " << deleted_edge.size() << std::endl;
            for (size_t i = 0; i < deleted_edge.size(); i++)
            {
                Mesh::HalfedgeHandle heh0 = mesh.halfedge_handle((deleted_edge[i]), 0);
                Mesh::HalfedgeHandle heh1 = mesh.halfedge_handle((deleted_edge[i]), 1);
                tmp_face_vhandles.clear();
                Mesh::HalfedgeHandle next = heh0;
                while (mesh.to_vertex_handle(next) != mesh.from_vertex_handle(heh0))
                {
                    tmp_face_vhandles.push_back(mesh.to_vertex_handle(next));
                    next = mesh.next_halfedge_handle(next);
                }
                next = heh1;
                while (mesh.to_vertex_handle(next) != mesh.from_vertex_handle(heh1))
                {
                    tmp_face_vhandles.push_back(mesh.to_vertex_handle(next));
                    next = mesh.next_halfedge_handle(next);
                }
                mesh.delete_edge(deleted_edge[i], false);
                mesh.add_face(tmp_face_vhandles);
            }
            mesh.garbage_collection();
            mesh.update_normals();
            std::vector<std::vector<int>> new_faces;
            for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
            {
                face.clear();
                for (Mesh::FaceVertexIter fv_it = mesh.fv_begin((*f_it)); fv_it != mesh.fv_end((*f_it)); ++fv_it)
                {
                    int idx = fv_it.handle().idx();
                    face.push_back(local_to_global[idx]);
                }
                new_faces.push_back(face);
            }
            triangles = new_faces;
            try {
                if (!OpenMesh::IO::write_mesh(mesh, "openmesh_output.off")) {
                    std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
                }
            }
            catch (std::exception& x)
            {
                std::cerr << x.what() << std::endl;
            }
            return true;
        }
        void construct(Mesh *mesh, const std::vector<double> &pts, const std::vector<int> &faces)
        {
            mesh->clear();
            int pts_size = pts.size() / 3;
            int face_size = faces.size() / 3;
            std::vector<Mesh::VertexHandle> vhandle;
            vhandle.resize(pts_size);
            for (int i = 0; i < pts_size; i++)
            {
                vhandle[i] = mesh->add_vertex(Mesh::Point(pts[3 * i], pts[3 * i + 1], pts[3 * i + 2]));
            }
            std::vector<Mesh::VertexHandle> tmp_face_vhandles;
            for (size_t i = 0; i < face_size; i++)
            {
                tmp_face_vhandles.clear();
                tmp_face_vhandles.push_back(vhandle[faces[3 * i]]);
                tmp_face_vhandles.push_back(vhandle[faces[3 * i + 1]]);
                tmp_face_vhandles.push_back(vhandle[faces[3 * i + 2]]);
                mesh->add_face(tmp_face_vhandles);
            }
            mesh->update_normals();
        }
        void calc_weights(Mesh *mesh, OpenMesh::EPropHandleT<Mesh::Scalar>& eweight_)
        {
            Mesh::VertexIter v_it, v_end(mesh->vertices_end());
            Mesh::EdgeIter e_it, e_end(mesh->edges_end());
            Mesh::VertexFaceIter vf_it;
            Mesh::FaceVertexIter fv_it;
            Mesh::HalfedgeHandle h0, h1, h2;
            Mesh::VertexHandle v0, v1;
            Mesh::Point p0, p1, p2, d0, d1;
            Mesh::Scalar w, area, b(0.99);
            for (e_it = mesh->edges_begin(); e_it != e_end; ++e_it)
            {
                w = 0.0;
                h0 = mesh->halfedge_handle(e_it.handle(), 0);
                v0 = mesh->to_vertex_handle(h0);
                p0 = mesh->point(v0);
                h1 = mesh->halfedge_handle(e_it.handle(), 1);
                v1 = mesh->to_vertex_handle(h1);
                p1 = mesh->point(v1);
                {
                    h2 = mesh->next_halfedge_handle(h0);
                    p2 = mesh->point(mesh->to_vertex_handle(h2));
                    d0 = (p0 - p2).normalize();
                    d1 = (p1 - p2).normalize();
                    w += 1.0 / tan(acos(std::max(-b, std::min(b, (d0 | d1)))));
                }
                {
                    h2 = mesh->next_halfedge_handle(h1);
                    p2 = mesh->point(mesh->to_vertex_handle(h2));
                    d0 = (p0 - p2).normalize();
                    d1 = (p1 - p2).normalize();
                    w += 1.0 / tan(acos(std::max(-b, std::min(b, (d0 | d1)))));
                }
                w = std::max(w, 0.0f);
                mesh->property(eweight_, e_it) = w;
            }
        }
        void Laplacian_smooth(Mesh &mesh_, OpenMesh::EPropHandleT<Mesh::Scalar>& eweight_, int max_iter, double damp_ratio)
        {
            mesh_.request_face_normals();
            Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
            Mesh::HalfedgeHandle h;
            Mesh::EdgeHandle e;
            Mesh::VertexVertexIter vv_it;
            Mesh::Point laplace(0.0, 0.0, 0.0);
            Mesh::Scalar w, ww;
            OpenMesh::VPropHandleT<Mesh::Point> vpos_;
            mesh_.add_property(vpos_);
            for (unsigned int iter = 0; iter < max_iter; ++iter)
            {
                for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
                {
                    laplace = Mesh::Point(0, 0, 0);
                    ww = 0.0;
                    if (!mesh_.is_boundary(v_it))
                    {
                        for (vv_it = mesh_.vv_iter(v_it); vv_it; ++vv_it)
                        {
                            h = vv_it.current_halfedge_handle();
                            e = mesh_.edge_handle(h);
                            w = mesh_.property(eweight_, e);
                            ww += w;
                            laplace += w * (mesh_.point(vv_it) - mesh_.point(v_it));
                        }
                        laplace /= ww;
                        laplace *= damp_ratio;
                    }
                    mesh_.property(vpos_, v_it) = mesh_.point(v_it) + laplace;
                }
                for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
                    mesh_.set_point(v_it, mesh_.property(vpos_, v_it));
            }
            mesh_.update_normals();
        }
        void Laplacian_smooth(std::vector<double> &pts, const std::vector<int>& faces, int iter, double damp_ratio)
        {
            Mesh mesh_;
            mesh_.request_face_status();
            mesh_.request_edge_status();
            mesh_.request_vertex_status();
            construct(&mesh_, pts, faces);
            OpenMesh::EPropHandleT<Mesh::Scalar> eweight_;
            mesh_.add_property(eweight_);
            calc_weights(&mesh_, eweight_);
            Laplacian_smooth(mesh_, eweight_, iter, damp_ratio);
            Mesh::VertexIter v_iter = mesh_.vertices_begin();
            Mesh::VertexIter v_end = mesh_.vertices_end();
            for (; v_iter != v_end; ++v_iter)
            {
                int idx = v_iter.handle().idx();
                pts[3 * idx] = mesh_.point(v_iter)[0];
                pts[3 * idx + 1] = mesh_.point(v_iter)[1];
                pts[3 * idx + 2] = mesh_.point(v_iter)[2];
            }
        }
    }
    namespace CGALHelper
    {
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef CGAL::Triangulation_vertex_base_2<K> Vb;
        typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
        typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
        typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
        typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
        typedef CDT::Vertex_handle Vertex_handle;
        typedef CDT::Point Point;
        void TriQuality(const QualityType qtype, const Vector3d& v0,
            const Vector3d& v1, const Vector3d& v2, double& quality) {
            double L01 = (v0 - v1).SquaredLength();
            double L02 = (v0 - v2).SquaredLength();
            double L12 = (v1 - v2).SquaredLength();
            double area = (double) 0.5 * ((v0 - v1).Cross(v0 - v2)).Length();
            if (qtype == SIN_THETA)
            {
                double theta[3];
                theta[0] = 4 * area * area / (L01 * L02);
                theta[1] = 4 * area * area / (L01 * L12);
                theta[2] = 4 * area * area / (L02 * L12);
                quality = (double)4.0 / (double)3.0 * *std::min_element(&theta[0], &theta[3]);
            }
            else if (qtype == RHO)
            {
                double l01 = sqrt(L01);
                double l02 = sqrt(L02);
                double l12 = sqrt(L12);
                quality = 16 * area * area / (l01 * l02 * l12 * (l01 + l02 + l12));
            }
            else if (qtype == ETA)
            {
                quality = 4 * sqrt((double) 3.0) * area / (L01 + L02 + L12);
            }
        }
        bool QuadQuality(const QualityType qtype, const Vector3d& v0,
            const Vector3d& v1, const Vector3d& v2,
            const Vector3d& v3, double& quality) {
            Vector3d V[4];
            V[0] = (v1 - v0).UnitCross(v2 - v0);
            V[1] = (v1 - v0).UnitCross(v3 - v0);
            V[2] = (v2 - v0).UnitCross(v3 - v0);
            V[3] = (v2 - v1).UnitCross(v3 - v1);
            double n[6];
            int count = 0;
            for (int i = 0; i < 3; i++) {
                for (int j = i + 1; j < 4; j++) {
                    n[count] = V[i].Dot(V[j]);
                    if (n[count] <= 0)
                        return false;
                    count++;
                }
            }
            double min_n = *std::min_element(&n[0], &n[6]);
            min_n = min_n < -1.f ? -1.f : min_n;
            min_n = min_n > 1.f ? 1.f : min_n;
            quality = (double) 1.0 - (double) 2.0 / M_PI * acos(min_n);
            double mu[4];
            TriQuality(qtype, v0, v1, v2, mu[0]);
            TriQuality(qtype, v0, v1, v3, mu[1]);
            TriQuality(qtype, v0, v2, v3, mu[2]);
            TriQuality(qtype, v1, v2, v3, mu[3]);
            quality *= *std::min_element(&mu[0], &mu[4]);
            if (qtype == SIN_THETA)
            {
                quality *= (double)1.5;
            }
            else if (qtype == RHO)
            {
                quality *= (sqrt(2.0) + (double)1.0)*(double)0.5;
            }
            else if (qtype == ETA)
            {
                quality *= (double)2.0*sqrt(3.0) / (double)3.0;
            }
            return true;
        }
        bool IsValidOrientedQuad(const Vector3d& v0, const Vector3d& v1,
            const Vector3d& v2, const Vector3d& v3) {
            Vector3d V[4];
            V[0] = (v1 - v0).Cross(v2 - v0);
            V[1] = (v1 - v0).Cross(v3 - v0);
            V[2] = (v2 - v0).Cross(v3 - v0);
            V[3] = (v2 - v1).Cross(v3 - v1);
            for (int i = 0; i < 3; i++) {
                for (int j = i + 1; j < 4; j++) {
                    if (V[i].Dot(V[j]) <= 0)
                        return false;
                }
            }
            return true;
        }
        double compute_cell_volume(const CVec<double, 3> &v1, const CVec<double, 3> &v2, const CVec<double, 3> &v3, const CVec<double, 3> &v4)
        {
            double v;
            CVec<double, 3> vec1, vec2, vec3, temp_vec;
            vec1 = v2 - v1;
            vec2 = v3 - v1;
            vec3 = v4 - v1;
            temp_vec = vec1.Cross(vec2);
            v = temp_vec.Dot(vec3);
            return v / 6.0;
        }
        double compute_quad_distortion(std::vector<Vector3d> &quad_array, std::vector<double> &face_distortion)
        {
            face_distortion.clear();
            Vector3d v0(-2.81941, 8.58536, -3.46881), v1(-1.91519, 9.02237, -3.12251), v2(-1.40247, 8.82325, -3.95416), v3(-2.30669, 8.38624, -4.30047);
            double quality = 0;
            QualityType qtype = RHO;
            double ave_distortion = 0;
            int length = quad_array.size() / 4;
            int accum = 0;
            for (size_t i = 0; i < length; i++)
            {
                double quality = 0;
                if (IsValidOrientedQuad(quad_array[4 * i], quad_array[4 * i + 1], quad_array[4 * i + 2], quad_array[4 * i + 3]))
                {
                    QuadQuality(qtype, quad_array[4 * i], quad_array[4 * i + 1], quad_array[4 * i + 2], quad_array[4 * i + 3], quality);
                    face_distortion.push_back(quality);
                    if (quality <= 1.f && quality >= -1.f)
                        ave_distortion += quality;
                    else
                    {
                        std::cout << "unvalid quality: " << quality << std::endl;
                        std::cout << "four verts: " << std::endl;
                        std::cout << quad_array[4 * i] << " _ " << quad_array[4 * i + 1] << " _ " << quad_array[4 * i + 2] << " _ " << quad_array[4 * i + 3] << std::endl;
                    }
                    accum++;
                }
                else
                {
                    face_distortion.push_back(0);
                    std::cout << "invalid quad result! " << std::endl;
                    std::cout << "invalid four verts: " << std::endl;
                    std::cout << quad_array[4 * i] << " _ " << quad_array[4 * i + 1] << " _ " << quad_array[4 * i + 2] << " _ " << quad_array[4 * i + 3] << std::endl;
                }
            }
            std::cout << "accum number: " << accum << std::endl;
            if (accum == 0)
                return 0;
            ave_distortion = ave_distortion / accum;
            return ave_distortion;
        }
        void localize3dvector(Vector3d input1, Vector3d input2, Vector3d normal, Vector2d& output1, Vector2d& output2)
        {
            normal = normal / normal.L2Norm();
            double theta = acos(normal[2]);
            double phi = acos(normal[0] / sqrt(1 - normal[2] * normal[2]));
            double phi_ref = asin(normal[1] / sqrt(1 - normal[2] * normal[2]));
            if (phi_ref < 0)
                phi = 2 * PI - phi;
            double sin_theta = sqrt(1 - normal[2] * normal[2]);
            double cos_theta = normal[2];
            if (abs(sin_theta) < epsilon)
            {
                output1[0] = input1[0];
                output1[1] = input1[1];
                output2[0] = input2[0];
                output2[1] = input2[1];
                return;
            }
            CMat<double, 3, 3> R_z(cos(phi), sin(phi), 0, -sin(phi), cos(phi), 0, 0, 0, 1);
            CMat<double, 3, 3> R_y(cos_theta, 0, -sin_theta, 0, 1, 0, sin_theta, 0, cos_theta);
            CMat<double, 3, 3> R_mult = R_y * R_z;
            CVec<double, 3> temp1 = R_mult * input1;
            CVec<double, 3> temp2 = R_mult * input2;
            output1[0] = temp1[0];
            output1[1] = temp1[1];
            output2[0] = temp2[0];
            output2[1] = temp2[1];
        }
        void localize3dvector(const std::vector<double> &coord_x, const std::vector<double> &coord_y, const std::vector<double> &coord_z, double normal_x, double normal_y, double normal_z, std::vector<double> &output_x, std::vector<double> &output_y)
        {
            Vector3d normal(normal_x, normal_y, normal_z);
            normal = normal / normal.L2Norm();
            double theta = acos(normal[2]);
            double phi = acos(normal[0] / sqrt(1 - normal[2] * normal[2]));
            double phi_ref = asin(normal[1] / sqrt(1 - normal[2] * normal[2]));
            if (phi_ref < 0)
                phi = 2 * PI - phi;
            double sin_theta = sqrt(1 - normal[2] * normal[2]);
            double cos_theta = normal[2];
            if (abs(sin_theta) < epsilon)
            {
                output_x = coord_x;
                output_y = coord_y;
                return;
            }
            CMat<double, 3, 3> R_z(cos(phi), sin(phi), 0, -sin(phi), cos(phi), 0, 0, 0, 1);
            CMat<double, 3, 3> R_y(cos_theta, 0, -sin_theta, 0, 1, 0, sin_theta, 0, cos_theta);
            CMat<double, 3, 3> R_mult = R_y * R_z;
            assert(coord_x.size() == coord_y.size() && coord_x.size() == coord_z.size());
            output_x.clear();
            output_y.clear();
            for (size_t i = 0; i < coord_x.size(); i++)
            {
                Vector3d input_vec(coord_x[i], coord_y[i], coord_z[i]);
                Vector3d output_vec = R_mult * input_vec;
                output_x.push_back(output_vec[0]);
                output_y.push_back(output_vec[1]);
            }
        }
        bool PointInPolygon(std::vector<Vector3d> &key_pts_array, Vector3d point)
        {
            bool is_in_polygon = true;
            double angle = 0;
            Vector3d dir1 = key_pts_array[0] - key_pts_array[1];
            Vector3d dir2 = key_pts_array[2] - key_pts_array[1];
            Vector3d face_normal = dir1.Cross(dir2);
            for (size_t i = 0; i < key_pts_array.size(); i++)
            {
                size_t cur = i;
                size_t next = (i + 1) % key_pts_array.size();
                Vector3d v1_o = key_pts_array[cur] - point;
                Vector3d v2_o = key_pts_array[next] - point;
                Vector2d v1, v2;
                localize3dvector(v1_o, v2_o, face_normal, v1, v2);
                double dot = v1[0] * v2[0] + v1[1] * v2[1];
                double det = v1[0] * v2[1] - v1[1] * v2[0];
                angle += atan2(det, dot);
            }
            if (abs(angle) < 6.f)
            {
                is_in_polygon = false;
            }
            return is_in_polygon;
        }
        void CrossVectorNormalize(double d1x, double d1y, double d1z, double d2x, double d2y, double d2z, double &rx, double &ry, double &rz)
        {
            Vector3d d1(d1x, d1y, d1z);
            Vector3d d2(d2x, d2y, d2z);
            Vector3d r = d1.Cross(d2);
            r.Normalize();
            if (r.L2Norm() < 0.1f)
            {
                std::cout << "Norm too small" << std::endl;
                std::cout << d1 << std::endl;
                std::cout << d2 << std::endl;
                std::cout << "d1: " << d1x << " " << d1y << " " << d1z << std::endl;
                std::cout << "d2: " << d2x << " " << d2y << " " << d2z << std::endl;
            }
            rx = r[0];
            ry = r[1];
            rz = r[2];
        }
        bool PointInPolygon(std::vector<double> &key_pts_array_x, std::vector<double> &key_pts_array_y, std::vector<double> &key_pts_array_z, double p_x, double p_y, double p_z)
        {
            assert((key_pts_array_x.size() == key_pts_array_y.size()) && (key_pts_array_z.size() == key_pts_array_y.size()) && (key_pts_array_x.size() == key_pts_array_z.size()));
            std::vector<Vector3d> key_pts_array;
            Vector3d point;
            point[0] = p_x;
            point[1] = p_y;
            point[2] = p_z;
            Vector3d temp_vec;
            int length = key_pts_array_x.size();
            for (size_t i = 0; i < length; i++)
            {
                temp_vec[0] = key_pts_array_x[i];
                temp_vec[1] = key_pts_array_y[i];
                temp_vec[2] = key_pts_array_z[i];
                key_pts_array.push_back(temp_vec);
            }
            return PointInPolygon(key_pts_array, point);
        }
        bool ToLeftTest(CVec<double, 3> p0, CVec<double, 3> p1, CVec<double, 3> p2)
        {
            double x1(p0[0]), y1(p0[1]), x2(p1[0]), y2(p1[1]), x3(p2[0]), y3(p2[1]);
            double area = x1 * y2 - y1 * x2 + x2 * y3 - y2 * x3 + x3 * y1 - y3 * x1;
            return (area > 0);
        }
        bool LineSegmentIntersectTest(Vector3d u0, Vector3d u1, Vector3d v0, Vector3d v1)
        {
            bool test1 = ToLeftTest(u0, u1, v0);
            bool test2 = ToLeftTest(u0, u1, v1);
            bool test3 = ToLeftTest(v0, v1, u0);
            bool test4 = ToLeftTest(v0, v1, u1);
            if (test1 + test2 == 1 && test3 + test4 == 1)
            {
                return true;
            }
            if (PointOnLineSegment(u0, u1, v0))
                return true;
            if (PointOnLineSegment(u0, u1, v1))
                return true;
            if (PointOnLineSegment(v0, v1, u0))
                return true;
            if (PointOnLineSegment(v0, v1, u1))
                return true;
            return false;
        }
        bool CutLegalTest(std::vector<Vector3d> &key_pts_array, Vector3d begin_point, Vector3d end_point)
        {
            if (!PointInPolygon(key_pts_array, end_point))
                return false;
            int begin_idx = -1;
            for (int i = 0; i < key_pts_array.size(); i++)
            {
                if ((key_pts_array[i] - begin_point).L2Norm() < EPSILON)
                {
                    begin_idx = i;
                    break;
                }
            }
            if (begin_idx == -1)
            {
                std::cout << "begin point not found in polygon !" << std::endl;
                return false;
            }
            int after_idx = (begin_idx + 1) % key_pts_array.size();
            int before_idx = (begin_idx - 1 + key_pts_array.size()) % key_pts_array.size();
            double cur_length = (end_point - begin_point).L2Norm();
            double before_length = (key_pts_array[before_idx] - begin_point).L2Norm();
            double after_length = (key_pts_array[after_idx] - begin_point).L2Norm();
            if (cur_length > before_length * 0.5 || cur_length > after_length * 0.5)
                return false;
            for (int i = 0; i < key_pts_array.size(); i++)
            {
                int after = (i + 1) % key_pts_array.size();
                if (i == begin_idx || after == begin_idx)
                    continue;
                if (LineSegmentIntersectTest(begin_point, end_point, key_pts_array[i], key_pts_array[after]))
                {
                    return false;
                }
            }
            return true;
        }
        bool PointOnLineSegment(CVec<double, 3> p0, CVec<double, 3> p1, CVec<double, 3> p2)
        {
            double x1(p0[0]), y1(p0[1]), x2(p1[0]), y2(p1[1]), x3(p2[0]), y3(p2[1]);
            double area = x1 * y2 - y1 * x2 + x2 * y3 - y2 * x3 + x3 * y1 - y3 * x1;
            if (abs(area) > epsilon)
                return false;
            double ratio = (x3 - x1) / (x2 - x1);
            if (ratio <= 1 && ratio >= 0)
                return true;
            return false;
        }
        void triangulation(std::vector<std::vector<int>> &cornerrecord, std::vector<std::vector<CVec<double, 3>>> &key_pts_array, std::vector<std::vector<int>> &faces, std::vector<size_t> &begin_pts_index, std::vector<CVec<double, 3>> &end_pts, double tri_area)
        {
            if (begin_pts_index.size() != end_pts.size())
            {
                std::cout << "select error " << std::endl;
                return;
            }
            std::cout << "corner number: " << begin_pts_index.size() << std::endl;
            std::vector<CVec<double, 3>> all_pts;
            std::vector<Point> polygon1;
            std::map<Point, int> point_to_corner;
            std::map<Point, int>::iterator it;
            CDT cdt;
            cornerrecord.clear();
            for (size_t j = 0; j < key_pts_array.size(); j++)
            {
                for (size_t i = 0; i < key_pts_array[j].size(); i++)
                {
                    all_pts.push_back(key_pts_array[j][i]);
                    polygon1.push_back((Point(key_pts_array[j][i][0], key_pts_array[j][i][1])));
                }
                size_t N_poly = polygon1.size();
                for (size_t i = 0; i < N_poly; i++)
                {
                    cdt.insert_constraint(polygon1[i], polygon1[(i + 1) % N_poly]);
                }
                polygon1.clear();
            }
            int N_poly = all_pts.size();
            CVec<double, 3> first, second;
            std::vector<int> onecorner;
            for (size_t i = 0; i < begin_pts_index.size(); i++)
            {
                first = all_pts[begin_pts_index[i]];
                second = end_pts[i];
                cdt.insert(Point(first[0], first[1]), Point(second[0], second[1]));
            }
            std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
            std::cout << "Meshing of triangulation..." << std::endl;
            std::list<Point> list_of_seeds;
            CVec<double, 3> p1, p2, p3, ave;
            if (key_pts_array.size() == 1)
            {
                CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0, tri_area));
            }
            else
            {
                size_t inner_size = key_pts_array[1].size();
                for (size_t i = 0; i < key_pts_array[1].size(); i++)
                {
                    p1 = key_pts_array[1][(i - 1 + inner_size) % inner_size];
                    p2 = key_pts_array[1][i];
                    p3 = key_pts_array[1][(i + 1) % inner_size];
                    if (ToLeftTest(p1, p2, p3))
                        break;
                }
                ave = (p1 + p2 + p3) / 3;
                list_of_seeds.push_back(Point(ave[0], ave[1]));
                CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),
                    Criteria(0, tri_area));
            }
            int N_allvert = cdt.number_of_vertices();
            all_pts.clear();
            std::map<Point, int> pointtoindex;
            std::map<int, int> cornervertmap;
            int ind = 0;
            Point temp_point;
            first[2] = 0;
            for (CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin();
                vit != cdt.finite_vertices_end(); ++vit)
            {
                temp_point = vit->point();
                cornervertmap[ind] = ind;
                pointtoindex.insert(std::pair<Point, int>(temp_point, ind++));
                first[0] = temp_point[0];
                first[1] = temp_point[1];
                all_pts.push_back(first);
            }
            for (size_t i = 0; i < begin_pts_index.size(); i++)
            {
                first = all_pts[begin_pts_index[i]];
                second = end_pts[i];
                std::vector<std::pair<int, double>> candidate_points_index;
                for (size_t j = 0; j < all_pts.size(); j++)
                {
                    if (ig::CGALHelper::PointOnLineSegment(first, second, all_pts[j]))
                    {
                        double temp_dis = (all_pts[j] - first).L2Norm();
                        candidate_points_index.push_back(std::pair<int, double>(j, temp_dis));
                    }
                }
                for (size_t j = 0; j < candidate_points_index.size(); j++)
                {
                    std::pair<int, double> temp_pair;
                    for (size_t k = j + 1; k < candidate_points_index.size(); k++)
                    {
                        if (candidate_points_index[k].second < candidate_points_index[j].second)
                        {
                            temp_pair = candidate_points_index[k];
                            candidate_points_index[k] = candidate_points_index[j];
                            candidate_points_index[j] = temp_pair;
                        }
                    }
                }
                onecorner.clear();
                for (size_t j = 0; j < candidate_points_index.size(); j++)
                {
                    onecorner.push_back(candidate_points_index[j].first);
                }
                cornerrecord.push_back(onecorner);
            }
            int vert_count = N_allvert;
            for (size_t i = 0; i < begin_pts_index.size(); i++)
            {
                double temp_split_num = cornerrecord[i].size() - 1;
                for (size_t j = 0; j < temp_split_num; j++)
                {
                    cornerrecord[i].push_back(vert_count + j);
                    cornervertmap[cornerrecord[i][temp_split_num - 1 - j]] = vert_count + j;
                    all_pts.push_back(all_pts[cornerrecord[i][temp_split_num - 1 - j]]);
                }
                std::cout << "corner size: " << cornerrecord[i].size() << std::endl;
                vert_count += temp_split_num;
            }
            std::vector<std::vector<std::vector<int>>> corner_record_no_rep;
            std::vector<std::vector<int>> one_corner_no_rep;
            std::map<int, int> rep_corner_idx_map;
            std::map<int, int>::iterator rep_corner_it;
            for (int i = 0; i < cornerrecord.size(); i++)
            {
                one_corner_no_rep.clear();
                int corner_idx = -1;
                for (int j = 0; j < corner_record_no_rep.size(); j++)
                {
                    if (corner_record_no_rep[j][0][0] == cornerrecord[i][0])
                    {
                        corner_idx = j;
                        break;
                    }
                }
                if (corner_idx == -1)
                {
                    one_corner_no_rep.push_back(cornerrecord[i]);
                    corner_record_no_rep.push_back(one_corner_no_rep);
                }
                else
                {
                    corner_record_no_rep[corner_idx].push_back(cornerrecord[i]);
                }
            }
            for (int i = 0; i < corner_record_no_rep.size(); i++)
            {
                if (corner_record_no_rep[i].size() != 1)
                {
                    int begin_idx = corner_record_no_rep[i][0][0];
                    rep_corner_idx_map[begin_idx] = i;
                    for (int j = 0; j < corner_record_no_rep[i].size(); j++)
                    {
                        int min_idx = j;
                        Vector3d begin_vec = all_pts[corner_record_no_rep[i][j][1]] - all_pts[begin_idx];
                        double min_value = 0;
                        for (int k = j + 1; k < corner_record_no_rep[i].size(); k++)
                        {
                            Vector3d temp_vec = all_pts[corner_record_no_rep[i][k][1]] - all_pts[begin_idx];
                            double cross_value = begin_vec.Cross(temp_vec)[2];
                            if (cross_value < min_value)
                            {
                                min_value = cross_value;
                                min_idx = k;
                            }
                        }
                        if (min_idx != j)
                        {
                            std::vector<int> temp_corner = corner_record_no_rep[i][min_idx];
                            corner_record_no_rep[i][min_idx] = corner_record_no_rep[i][j];
                            corner_record_no_rep[i][j] = temp_corner;
                        }
                    }
                    for (int j = 0; j < corner_record_no_rep[i].size(); j++)
                    {
                        if (j != 0)
                        {
                            corner_record_no_rep[i][j][0] = corner_record_no_rep[i][j - 1].back();
                            cornervertmap[corner_record_no_rep[i][j][0]] = corner_record_no_rep[i][j].back();
                        }
                        for (int k = 0; k < corner_record_no_rep[i][j].size(); k++)
                        {
                            std::cout << " " << corner_record_no_rep[i][j][k];
                        }
                        std::cout << std::endl;
                    }
                }
            }
            cornerrecord.clear();
            for (int i = 0; i < corner_record_no_rep.size(); i++)
            {
                {
                    for (int j = 0; j < corner_record_no_rep[i].size(); j++)
                    {
                        for (int k = 0; k < corner_record_no_rep[i][j].size(); k++)
                        {
                            int temp_idx = corner_record_no_rep[i][j][k];
                            point_to_corner[Point(all_pts[temp_idx][0], all_pts[temp_idx][1])] = cornerrecord.size();
                        }
                        cornerrecord.push_back(corner_record_no_rep[i][j]);
                    }
                }
            }
            std::cout << "corner output: " << std::endl;
            for (int i = 0; i < cornerrecord.size(); i++)
            {
                std::cout << "i: " << i << std::endl;
                for (int j = 0; j < cornerrecord[i].size(); j++)
                {
                    std::cout << " " << cornerrecord[i][j];
                }
                std::cout << std::endl;
            }
            faces.clear();
            std::vector<int> face;
            int cornermark;
            int cornerind;
            CDT::Edge_iterator eit;
            int count = 0;
            for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
                fit != cdt.finite_faces_end(); ++fit)
            {
                if (fit->is_in_domain())
                {
                    face.clear();
                    ++count;
                    Point p0 = fit->vertex(0)->point();
                    Point p1 = fit->vertex(1)->point();
                    Point p2 = fit->vertex(2)->point();
                    CVec<double, 3> p_ave((p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0, 0);
                    cornermark = 0;
                    int repcornermark = 0;
                    int idx0, idx1, idx2;
                    std::map<int, int>::iterator rep_corner_select;
                    idx0 = pointtoindex.at(fit->vertex(0)->point());
                    idx1 = pointtoindex.at(fit->vertex(1)->point());
                    idx2 = pointtoindex.at(fit->vertex(2)->point());
                    rep_corner_it = rep_corner_idx_map.find(idx0);
                    if (rep_corner_it != rep_corner_idx_map.end())
                    {
                        repcornermark++;
                        rep_corner_select = rep_corner_it;
                    }
                    rep_corner_it = rep_corner_idx_map.find(idx1);
                    if (rep_corner_it != rep_corner_idx_map.end())
                    {
                        repcornermark++;
                        rep_corner_select = rep_corner_it;
                    }
                    rep_corner_it = rep_corner_idx_map.find(idx2);
                    if (rep_corner_it != rep_corner_idx_map.end())
                    {
                        repcornermark++;
                        rep_corner_select = rep_corner_it;
                    }
                    if (repcornermark > 0)
                    {
                        std::cout << "rep corner mark: " << repcornermark << std::endl;
                        int rep_array_idx = rep_corner_select->second;
                        int tri_corner_idx = corner_record_no_rep[rep_array_idx][0][0];
                        std::cout << "tri corner idx: " << tri_corner_idx << std::endl;
                        Vector3d triangle_vec = p_ave - all_pts[tri_corner_idx];
                        Vector3d first_vec = all_pts[corner_record_no_rep[rep_array_idx][0][1]] - all_pts[tri_corner_idx];
                        if (first_vec.Cross(triangle_vec)[2] < 0)
                        {
                            face.push_back(idx0);
                            face.push_back(idx1);
                            face.push_back(idx2);
                        }
                        else
                        {
                            int k = 1;
                            for (; k < corner_record_no_rep[rep_array_idx].size(); k++)
                            {
                                Vector3d temp_vec = all_pts[corner_record_no_rep[rep_array_idx][k][1]] - all_pts[tri_corner_idx];
                                if (temp_vec.Cross(triangle_vec)[2] < 0)
                                    break;
                            }
                            k--;
                            std::vector<int> temp_corner = corner_record_no_rep[rep_array_idx][k];
                            std::vector<int>::iterator temp_it;
                            if (idx0 == tri_corner_idx)
                            {
                                idx0 = corner_record_no_rep[rep_array_idx][k].back();
                            }
                            else
                            {
                                temp_it = std::find(temp_corner.begin(), temp_corner.end(), idx0);
                                if (temp_it != temp_corner.end())
                                {
                                    idx0 = cornervertmap[idx0];
                                }
                            }
                            if (idx1 == tri_corner_idx)
                            {
                                idx1 = corner_record_no_rep[rep_array_idx][k].back();
                            }
                            else
                            {
                                temp_it = std::find(temp_corner.begin(), temp_corner.end(), idx1);
                                if (temp_it != temp_corner.end())
                                {
                                    idx1 = cornervertmap[idx1];
                                }
                            }
                            if (idx2 == tri_corner_idx)
                            {
                                idx2 = corner_record_no_rep[rep_array_idx][k].back();
                            }
                            else
                            {
                                temp_it = std::find(temp_corner.begin(), temp_corner.end(), idx2);
                                if (temp_it != temp_corner.end())
                                {
                                    idx2 = cornervertmap[idx2];
                                }
                            }
                            face.push_back(idx0);
                            face.push_back(idx1);
                            face.push_back(idx2);
                        }
                        faces.push_back(face);
                        continue;
                    }
                    int cornerind0(-1), cornerind1(-1), cornerind2(-1);
                    it = point_to_corner.find(p0);
                    if (it != point_to_corner.end())
                    {
                        cornermark++;
                        cornerind0 = point_to_corner.find(p0)->second;
                    }
                    it = point_to_corner.find(p1);
                    if (it != point_to_corner.end())
                    {
                        cornermark++;
                        cornerind1 = point_to_corner.find(p1)->second;
                    }
                    it = point_to_corner.find(p2);
                    if (it != point_to_corner.end())
                    {
                        cornermark++;
                        cornerind2 = point_to_corner.find(p2)->second;
                    }
                    if (cornermark < 1)
                    {
                        face.push_back(idx0);
                        face.push_back(idx1);
                        face.push_back(idx2);
                        faces.push_back(face);
                    }
                    else
                    {
                        cornerind = all_pts.size();
                        if (cornerind > cornerind0 && cornerind0 > -1)
                            cornerind = cornerind0;
                        if (cornerind > cornerind1 && cornerind1 > -1)
                            cornerind = cornerind1;
                        if (cornerind > cornerind2 && cornerind2 > -1)
                            cornerind = cornerind2;
                        if (ToLeftTest(all_pts[cornerrecord[cornerind][1]], all_pts[cornerrecord[cornerind][0]], p_ave))
                        {
                            face.push_back(pointtoindex.at(fit->vertex(0)->point()));
                            face.push_back(pointtoindex.at(fit->vertex(1)->point()));
                            face.push_back(pointtoindex.at(fit->vertex(2)->point()));
                            faces.push_back(face);
                        }
                        else
                        {
                            std::vector<int> temp_corner = cornerrecord[cornerind];
                            std::vector<int>::iterator temp_it;
                            temp_it = std::find(temp_corner.begin(), temp_corner.end(), idx0);
                            if (temp_it != temp_corner.end())
                            {
                                idx0 = cornervertmap[idx0];
                            }
                            temp_it = std::find(temp_corner.begin(), temp_corner.end(), idx1);
                            if (temp_it != temp_corner.end())
                            {
                                idx1 = cornervertmap[idx1];
                            }
                            temp_it = std::find(temp_corner.begin(), temp_corner.end(), idx2);
                            if (temp_it != temp_corner.end())
                            {
                                idx2 = cornervertmap[idx2];
                            }
                            face.push_back(idx0);
                            face.push_back(idx1);
                            face.push_back(idx2);
                            faces.push_back(face);
                        }
                    }
                }
            }
            key_pts_array.clear();
            key_pts_array.push_back(all_pts);
            std::cout << "There are " << all_pts.size() << " verts in the domain." << std::endl;
            std::cout << "There are " << count << " facets in the domain." << std::endl;
        }
        void triangulation(std::vector<std::vector<int>> &cornerrecord, std::vector<std::vector<CVec<double, 3>>> &key_pts_array, std::vector<std::vector<int>> &faces, std::vector<size_t> &begin_pts_index, std::vector<CVec<double, 3>> &end_pts, int& split_num, double tri_area)
        {
            std::vector<CVec<double, 3>> all_pts;
            std::vector<Point> polygon1;
            std::map<Point, int> point_to_corner;
            std::map<Point, int>::iterator it;
            CDT cdt;
            cornerrecord.clear();
            for (size_t j = 0; j < key_pts_array.size(); j++)
            {
                for (size_t i = 0; i < key_pts_array[j].size(); i++)
                {
                    all_pts.push_back(key_pts_array[j][i]);
                    polygon1.push_back((Point(key_pts_array[j][i][0], key_pts_array[j][i][1])));
                }
                size_t N_poly = polygon1.size();
                for (size_t i = 0; i < N_poly; i++)
                {
                    cdt.insert_constraint(polygon1[i], polygon1[(i + 1) % N_poly]);
                }
                polygon1.clear();
            }
            int N_poly = all_pts.size();
            CVec<double, 3> first, second;
            std::vector<int> onecorner;
            for (size_t i = 0; i < begin_pts_index.size(); i++)
            {
                onecorner.clear();
                onecorner.push_back(begin_pts_index[i]);
                first = all_pts[begin_pts_index[i]];
                second = end_pts[i];
                CVec<double, 3> step = (second - first) / split_num;
                for (size_t j = 0; j < split_num; j++)
                {
                    onecorner.push_back(N_poly + i * split_num + j);
                    point_to_corner[Point(first[0] + j * step[0], first[1] + j * step[1])] = i;
                    cdt.insert(Point(first[0] + j * step[0], first[1] + j * step[1]), Point(first[0] + (j + 1) * step[0], first[1] + (j + 1) * step[1]));
                }
                point_to_corner[Point(first[0] + split_num * step[0], first[1] + split_num * step[1])] = i;
                cornerrecord.push_back(onecorner);
            }
            std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
            std::cout << "Meshing of triangulation..." << std::endl;
            std::list<Point> list_of_seeds;
            CVec<double, 3> p1, p2, p3, ave;
            if (key_pts_array.size() == 1)
            {
                CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0, tri_area));
            }
            else
            {
                size_t inner_size = key_pts_array[1].size();
                for (size_t i = 0; i < key_pts_array[1].size(); i++)
                {
                    p1 = key_pts_array[1][(i - 1 + inner_size) % inner_size];
                    p2 = key_pts_array[1][i];
                    p3 = key_pts_array[1][(i + 1) % inner_size];
                    if (ToLeftTest(p1, p2, p3))
                        break;
                }
                ave = (p1 + p2 + p3) / 3;
                list_of_seeds.push_back(Point(ave[0], ave[1]));
                CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),
                    Criteria(0, tri_area));
            }
            std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
            int N_allvert = cdt.number_of_vertices();
            all_pts.clear();
            std::map<Point, int> pointtoindex;
            std::map<int, int> cornervertmap;
            int ind = 0;
            std::cout << "construct point map" << std::endl;
            Point temp_point;
            first[2] = 0;
            for (CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin();
                vit != cdt.finite_vertices_end(); ++vit)
            {
                temp_point = vit->point();
                cornervertmap[ind] = ind;
                pointtoindex.insert(std::pair<Point, int>(temp_point, ind++));
                first[0] = temp_point[0];
                first[1] = temp_point[1];
                all_pts.push_back(first);
            }
            for (size_t i = 0; i < begin_pts_index.size(); i++)
            {
                for (size_t j = 0; j < split_num; j++)
                {
                    cornerrecord[i].push_back(N_allvert + i * split_num + j);
                    cornervertmap[cornerrecord[i][split_num - 1 - j]] = N_allvert + i * split_num + j;
                    all_pts.push_back(all_pts[cornerrecord[i][split_num - 1 - j]]);
                }
            }
            faces.clear();
            std::vector<int> face;
            int cornermark;
            int cornerind;
            CDT::Edge_iterator eit;
            int count = 0;
            for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
                fit != cdt.finite_faces_end(); ++fit)
            {
                if (fit->is_in_domain())
                {
                    face.clear();
                    ++count;
                    Point p0 = fit->vertex(0)->point();
                    Point p1 = fit->vertex(1)->point();
                    Point p2 = fit->vertex(2)->point();
                    CVec<double, 3> p_ave((p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0, 0);
                    cornermark = 0;
                    it = point_to_corner.find(p0);
                    if (it != point_to_corner.end())
                    {
                        cornermark++;
                        cornerind = point_to_corner.find(p0)->second;
                    }
                    it = point_to_corner.find(p1);
                    if (it != point_to_corner.end())
                    {
                        cornermark++;
                        cornerind = point_to_corner.find(p1)->second;
                    }
                    it = point_to_corner.find(p2);
                    if (it != point_to_corner.end())
                    {
                        cornermark++;
                        cornerind = point_to_corner.find(p2)->second;
                    }
                    if (cornermark < 1)
                    {
                        face.push_back(pointtoindex.at(fit->vertex(0)->point()));
                        face.push_back(pointtoindex.at(fit->vertex(1)->point()));
                        face.push_back(pointtoindex.at(fit->vertex(2)->point()));
                        faces.push_back(face);
                    }
                    else
                    {
                        if (ToLeftTest(end_pts[cornerind], all_pts[begin_pts_index[cornerind]], p_ave))
                        {
                            face.push_back(pointtoindex.at(fit->vertex(0)->point()));
                            face.push_back(pointtoindex.at(fit->vertex(1)->point()));
                            face.push_back(pointtoindex.at(fit->vertex(2)->point()));
                            faces.push_back(face);
                        }
                        else
                        {
                            face.push_back(cornervertmap[pointtoindex.at(fit->vertex(0)->point())]);
                            face.push_back(cornervertmap[pointtoindex.at(fit->vertex(1)->point())]);
                            face.push_back(cornervertmap[pointtoindex.at(fit->vertex(2)->point())]);
                            faces.push_back(face);
                        }
                    }
                }
            }
            key_pts_array.clear();
            key_pts_array.push_back(all_pts);
            std::cout << "There are " << all_pts.size() << " verts in the domain." << std::endl;
            std::cout << "There are " << count << " facets in the domain." << std::endl;
        }
        void triangulation(std::vector<CVec<double, 3>> &key_pts, std::vector<std::vector<int>> &faces, std::vector<size_t> &begin_pts_index, std::vector<CVec<double, 3>> &end_pts, double tri_area, int split_num)
        {
            std::vector<Point> polygon1;
            CDT cdt;
            for (size_t i = 0; i < key_pts.size(); i++)
            {
                polygon1.push_back((Point(key_pts[i][0], key_pts[i][1])));
            }
            size_t N_poly = polygon1.size();
            for (size_t i = 0; i < N_poly; i++)
            {
                cdt.insert_constraint(polygon1[i], polygon1[(i + 1) % N_poly]);
            }
            CVec<double, 3> first, second;
            for (size_t i = 0; i < begin_pts_index.size(); i++)
            {
                first = key_pts[begin_pts_index[i]];
                second = end_pts[i];
                CVec<double, 3> step = (second - first) / split_num;
                for (size_t j = 0; j < split_num; j++)
                {
                    cdt.insert(Point(first[0] + j * step[0], first[1] + j * step[1]), Point(first[0] + (j + 1) * step[0], first[1] + (j + 1) * step[1]));
                }
            }
            std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
            std::cout << "Meshing of triangulation..." << std::endl;
            CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0, tri_area));
            std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
            key_pts.clear();
            std::map<Point, int> pointtoindex;
            int ind = 0;
            std::cout << "construct point map" << std::endl;
            Point temp_point;
            first[2] = 0;
            for (CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin();
                vit != cdt.finite_vertices_end(); ++vit)
            {
                temp_point = vit->point();
                pointtoindex.insert(std::pair<Point, int>(temp_point, ind++));
                first[0] = temp_point[0];
                first[1] = temp_point[1];
                key_pts.push_back(first);
            }
            faces.clear();
            std::vector<int> face;
            CDT::Edge_iterator eit;
            int count = 0;
            for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
                fit != cdt.finite_faces_end(); ++fit)
            {
                if (fit->is_in_domain())
                {
                    face.clear();
                    ++count;
                    face.push_back(pointtoindex.at(fit->vertex(0)->point()));
                    face.push_back(pointtoindex.at(fit->vertex(1)->point()));
                    face.push_back(pointtoindex.at(fit->vertex(2)->point()));
                    faces.push_back(face);
                }
            }
            std::cout << "There are " << count << " facets in the domain." << std::endl;
        }
        void triangulation(std::vector<CVec<double, 3>> &key_pts, std::vector<std::vector<int>> &faces, double tri_area)
        {
            std::vector<size_t> begin_pts_index;
            std::vector<CVec<double, 3>> end_pts;
            triangulation(key_pts, faces, begin_pts_index, end_pts, tri_area);
        }
        void triangulation(std::vector<std::vector<CVec<double, 3>>> &key_pts_array, std::vector<std::vector<int>> &faces, double tri_area)
        {
            std::vector<Point> polygon1;
            CDT cdt;
            for (size_t j = 0; j < key_pts_array.size(); j++)
            {
                for (size_t i = 0; i < key_pts_array[j].size(); i++)
                {
                    polygon1.push_back((Point(key_pts_array[j][i][0], key_pts_array[j][i][1])));
                }
                size_t N_poly = polygon1.size();
                for (size_t i = 0; i < N_poly; i++)
                {
                    cdt.insert_constraint(polygon1[i], polygon1[(i + 1) % N_poly]);
                }
                polygon1.clear();
            }
            std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
            std::cout << "Meshing of triangulation..." << std::endl;
            std::list<Point> list_of_seeds;
            CVec<double, 3> p1, p2, p3, ave, first;
            if (key_pts_array.size() == 1)
            {
                CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0, tri_area));
            }
            else
            {
                p1 = key_pts_array[1][0];
                p2 = key_pts_array[1][1];
                p3 = key_pts_array[1][2];
                ave = (p1 + p2 + p3) / 3;
                list_of_seeds.push_back(Point(ave[0], ave[1]));
                CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),
                    Criteria(0, tri_area));
            }
            std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
            std::vector<CVec<double, 3>> key_pts;
            std::map<Point, int> pointtoindex;
            int ind = 0;
            std::cout << "construct point map" << std::endl;
            Point temp_point;
            first[2] = 0;
            for (CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin();
                vit != cdt.finite_vertices_end(); ++vit)
            {
                temp_point = vit->point();
                pointtoindex.insert(std::pair<Point, int>(temp_point, ind++));
                first[0] = temp_point[0];
                first[1] = temp_point[1];
                key_pts.push_back(first);
            }
            faces.clear();
            std::vector<int> face;
            CDT::Edge_iterator eit;
            int count = 0;
            for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
                fit != cdt.finite_faces_end(); ++fit)
            {
                if (fit->is_in_domain())
                {
                    face.clear();
                    ++count;
                    face.push_back(pointtoindex.at(fit->vertex(0)->point()));
                    face.push_back(pointtoindex.at(fit->vertex(1)->point()));
                    face.push_back(pointtoindex.at(fit->vertex(2)->point()));
                    faces.push_back(face);
                }
            }
            key_pts_array.clear();
            key_pts_array.push_back(key_pts);
            std::cout << "There are " << count << " facets in the domain." << std::endl;
        }
        template <class T>
        void printnumber(T i)
        {
            std::cout << "print number " << i << std::endl;
        }
    }
    namespace CGALHelper
    {
        typedef CGAL::Simple_cartesian<double> Ker;
        typedef Ker::FT FT;
        typedef Ker::Point_3 Point3;
        typedef Ker::Point_2 Point2;
        typedef std::vector<FT> Scalar_vector;
        typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Ker> Triangle_coordinates;
        typedef Ker::Segment_3 Segment;
        typedef CGAL::Polyhedron_3<Ker, CGAL::Polyhedron_items_with_id_3> Polyhedron;
        typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
        typedef CGAL::AABB_traits<Ker, Primitive> Traits;
        typedef CGAL::AABB_tree<Traits> Tree;
        typedef Tree::Point_and_primitive_id Point_and_primitive_id;
        typedef Polyhedron::HalfedgeDS HalfedgeDS;
        template <class HDS>
        class Build_triangle : public CGAL::Modifier_base<HDS> {
        public:
            Build_triangle() {}
            void ConstructList(std::vector<Point3> PointList, std::vector<std::vector<int>> FaceList)
            {
                point_list = PointList;
                face_list = FaceList;
            }
            void operator()(HDS& hds) {
                CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
                B.begin_surface(point_list.size(), face_list.size());
                typedef typename HDS::Vertex Vertex;
                typedef typename Vertex::Point Point;
                for (int i = 0; i < point_list.size(); i++)
                {
                    B.add_vertex(point_list[i]);
                }
                for (int i = 0; i < face_list.size(); i++)
                {
                    B.begin_facet();
                    for (int j = 0; j < face_list[i].size(); j++)
                    {
                        B.add_vertex_to_facet(face_list[i][j]);
                    }
                    B.end_facet();
                }
                B.end_surface();
            }
            std::vector<Point3> point_list;
            std::vector<std::vector<int>> face_list;
        };
        void construct_barycentric_coordinate(std::vector<CVec<double, 3>>& source_points, std::vector<std::vector<int>>& source_faces, std::vector<CVec<double, 3>>&target_points, std::vector<CVec<double, 3>>&barycentric_coord, std::vector<int> &face_idx_list)
        {
            Polyhedron polyhedron;
            std::vector<Point3> PointList;
            std::vector<std::vector<int>> FaceList;
            std::vector<int> face;
            for (size_t i = 0; i < source_points.size(); i++)
            {
                PointList.push_back(Point3(source_points[i][0], source_points[i][1], source_points[i][2]));
            }
            for (size_t i = 0; i < source_faces.size(); i++)
            {
                face.clear();
                face.push_back(source_faces[i][0]);
                face.push_back(source_faces[i][1]);
                face.push_back(source_faces[i][2]);
                FaceList.push_back(face);
            }
            Build_triangle<HalfedgeDS> triangles;
            triangles.ConstructList(PointList, FaceList);
            polyhedron.delegate(triangles);
            int count = 0;
            for (Polyhedron::Face_iterator i = polyhedron.facets_begin(); i != polyhedron.facets_end(); ++i)
                i->id() = count++;
            count = 0;
            for (Polyhedron::Vertex_iterator viter = polyhedron.vertices_begin(); viter != polyhedron.vertices_end(); ++viter)
                viter->id() = count++;
            Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
            tree.accelerate_distance_queries();
            face_idx_list.clear();
            barycentric_coord.clear();
            for (size_t i = 0; i < target_points.size(); i++)
            {
                Point3 query(target_points[i][0], target_points[i][1], target_points[i][2]);
                Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
                Polyhedron::Face_handle f = pp.second;
                face_idx_list.push_back(f->id());
                Point3 p1 = f->halfedge()->vertex()->point();
                Point3 p2 = f->halfedge()->next()->vertex()->point();
                Point3 p3 = f->halfedge()->next()->next()->vertex()->point();
                Scalar_vector coordinates;
                const Point2 first_vertex(p1[0], p1[1]);
                const Point2 second_vertex(p2[0], p2[1]);
                const Point2 third_vertex(p3[0], p3[1]);
                Triangle_coordinates triangle_coordinates(first_vertex, second_vertex, third_vertex);
                Point2 querypoint(target_points[i][0], target_points[i][1]);
                coordinates.reserve(3);
                triangle_coordinates(querypoint, std::inserter(coordinates, coordinates.end()));
                barycentric_coord.push_back(CVec<double, 3>(coordinates[0], coordinates[1], coordinates[2]));
            }
            return;
        }
        void project_points_on_faces(const std::vector<double> &input_pts, const std::vector<int> &input_faces, std::vector<double> &proj_pts)
        {
            std::cout << "project points on faces" << std::endl;
            Polyhedron polyhedron;
            std::vector<Point3> PointList;
            std::vector<std::vector<int>> FaceList;
            std::vector<int> face;
            int input_pts_size = input_pts.size() / 3;
            int input_face_size = input_faces.size() / 3;
            for (size_t i = 0; i < input_pts_size; i++)
            {
                PointList.push_back(Point3(input_pts[3 * i], input_pts[3 * i + 1], input_pts[3 * i + 2]));
            }
            for (size_t i = 0; i < input_face_size; i++)
            {
                face.clear();
                face.push_back(input_faces[3 * i + 0]);
                face.push_back(input_faces[3 * i + 1]);
                face.push_back(input_faces[3 * i + 2]);
                FaceList.push_back(face);
            }
            Build_triangle<HalfedgeDS> triangles;
            triangles.ConstructList(PointList, FaceList);
            polyhedron.delegate(triangles);
            int count = 0;
            for (Polyhedron::Face_iterator i = polyhedron.facets_begin(); i != polyhedron.facets_end(); ++i)
                i->id() = count++;
            count = 0;
            for (Polyhedron::Vertex_iterator viter = polyhedron.vertices_begin(); viter != polyhedron.vertices_end(); ++viter)
                viter->id() = count++;
            Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
            tree.accelerate_distance_queries();
            int proj_pts_size = proj_pts.size() / 3;
            for (int i = 0; i < proj_pts_size; i++)
            {
                Point3 query(proj_pts[3 * i], proj_pts[3 * i + 1], proj_pts[3 * i + 2]);
                Point3 closest = tree.closest_point(query);
                proj_pts[3 * i] = closest[0];
                proj_pts[3 * i + 1] = closest[1];
                proj_pts[3 * i + 2] = closest[2];
            }
        }
        bool gauss_elimination(std::vector<std::vector<double>> &A)
        {
            int n = A.size();
            int m = A[0].size();
            int begin_idx = 0;
            for (int i = 0; i < n; i++) {
                double maxEl = abs(A[i][i]);
                int maxRow = i;
                begin_idx = i;
                for (int k = i + 1; k < n; k++) {
                    if (abs(A[k][begin_idx]) > maxEl) {
                        maxEl = abs(A[k][begin_idx]);
                        maxRow = k;
                    }
                }
                while (abs(A[maxRow][begin_idx]) < epsilon)
                {
                    begin_idx++;
                    maxEl = abs(A[i][begin_idx]);
                    for (int k = i + 1; k < n; k++) {
                        if (abs(A[k][begin_idx]) > maxEl) {
                            maxEl = abs(A[k][begin_idx]);
                            maxRow = k;
                        }
                    }
                }
                if (abs(A[maxRow][begin_idx]) > epsilon)
                {
                    for (int k = begin_idx; k < m; k++) {
                        double tmp = A[maxRow][k];
                        A[maxRow][k] = A[i][k];
                        A[i][k] = tmp;
                    }
                    for (int k = begin_idx + 1; k < m; k++)
                    {
                        A[i][k] = A[i][k] / A[i][begin_idx];
                    }
                    A[i][begin_idx] = 1;
                    for (int k = i + 1; k < n; k++) {
                        double c = -A[k][begin_idx] / A[i][begin_idx];
                        for (int j = begin_idx; j < m; j++) {
                            if (i == j) {
                                A[k][j] = 0;
                            }
                            else {
                                A[k][j] += c * A[i][j];
                            }
                        }
                    }
                    for (int k = 0; k < i; k++) {
                        double c = -A[k][begin_idx] / A[i][begin_idx];
                        for (int j = begin_idx; j < m; j++) {
                            if (i == j) {
                                A[k][j] = 0;
                            }
                            else {
                                A[k][j] += c * A[i][j];
                            }
                        }
                    }
                }
            }
            return true;
        }
        void cut_array(const std::vector<double> &input, double sigma_ratio, std::vector<double> &output, double &output_max, double &output_mean)
        {
            double mean, var, std_var;
            mean = 0;
            var = 0;
            for (int i = 0; i < input.size(); i++)
            {
                mean += input[i];
            }
            mean = mean / input.size();
            for (int i = 0; i < input.size(); i++)
            {
                var += (input[i] - mean) * (input[i] - mean);
            }
            var = var / (input.size() - 1);
            std_var = sqrt(var);
            std::cout << "input size: " << input.size() << std::endl;
            std::cout << "---------mean value: " << mean << std::endl;
            std::cout << "---------variance: " << var << std::endl;
            double lb = mean - sigma_ratio * std_var;
            double hb = mean + sigma_ratio * std_var;
            output_max = -1;
            output_mean = 0;
            output.clear();
            for (int i = 0; i < input.size(); i++)
            {
                if (input[i] > lb && input[i] < hb)
                {
                    output.push_back(input[i]);
                    if (output_max < input[i])
                        output_max = input[i];
                    output_mean += input[i];
                }
            }
            output_mean = output_mean / output.size();
            std::cout << "output size: " << output.size() << std::endl;
        }
        using namespace boost;
        void graph_coloring(int N, std::vector<std::pair<int, int>>& graph_color_edge, std::vector< std::vector<int> >& v_same_color, std::vector<int>& v_color)
        {
            typedef adjacency_list<listS, vecS, undirectedS> Graph;
            typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
            typedef graph_traits<Graph>::vertices_size_type vertices_size_type;
            typedef property_map<Graph, vertex_index_t>::const_type vertex_index_map;
            Graph g;
            for (int i = 0; i < N; ++i)
            {
                add_vertex(g);
            }
            int edge_size = graph_color_edge.size();
            for (int i = 0; i < edge_size; i++)
            {
                add_edge(graph_color_edge[i].first, graph_color_edge[i].second, g);
            }
            boost::vector_property_map<vertex_descriptor> order;
            smallest_last_vertex_ordering(g, order);
            std::vector<vertices_size_type> color_vec(num_vertices(g));
            iterator_property_map<vertices_size_type*, vertex_index_map>
                color(&color_vec.front(), get(vertex_index, g));
            vertices_size_type num_colors = sequential_vertex_coloring(g, order, color);
            v_same_color.clear(); v_same_color.resize(num_colors);
            for (int i = 0; i < N; ++i)
            {
                v_same_color[color_vec[i]].push_back(i);
            }
            graph_color_edge.clear();
            v_color.clear();
            for (int i = 0; i < N; i++)
            {
                v_color.push_back(color_vec[i] + 1);
            }
        }
        bool all_pair_shortest_path(int N, const std::vector<std::pair<int, int>>& input_edges, const std::vector<double> &input_weight, std::vector<std::vector<double>> &pair_paths)
        {
            typedef double t_weight;
            typedef boost::property<boost::edge_weight_t, t_weight> EdgeWeightProperty;
            typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                boost::no_property, EdgeWeightProperty> Graph;
            typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightMap;
            typedef boost::exterior_vertex_property<Graph, t_weight> DistanceProperty;
            typedef DistanceProperty::matrix_type DistanceMatrix;
            typedef DistanceProperty::matrix_map_type DistanceMatrixMap;
            Graph g;
            assert(input_edges.size() == input_weight.size());
            const int num_edges = input_edges.size();
            for (std::size_t k = 0; k < num_edges; ++k)
                boost::add_edge(input_edges[k].first, input_edges[k].second, input_weight[k], g);
            WeightMap weight_pmap = boost::get(boost::edge_weight, g);
            DistanceMatrix distances(num_vertices(g));
            DistanceMatrixMap dm(distances, g);
            bool valid = floyd_warshall_all_pairs_shortest_paths(g, dm,
                boost::weight_map(weight_pmap));
            if (!valid) {
                std::cerr << "Error - Negative cycle in matrix" << std::endl;
                return -1;
            }
            pair_paths.clear();
            pair_paths.resize(N, std::vector<double>(N, -1.0));
            assert(num_vertices(g) == N);
            for (std::size_t i = 0; i < num_vertices(g); ++i) {
                for (std::size_t j = i; j < num_vertices(g); ++j) {
                    pair_paths[i][j] = distances[i][j];
                    pair_paths[j][i] = distances[i][j];
                }
            }
            return 0;
        }
    }
}
