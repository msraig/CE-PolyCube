#pragma once
#include "SmallVec.h"
#include <vector>
#include <set>
#define EPSILON 0.001
namespace ig
{
    namespace IVF
    {
        void IVF_svd(const std::vector<Vector3d> &pts_ref, const std::vector<unsigned int> &elems, std::vector<Vector3d> &pts_flip, const std::vector<int> &boundary_tag, bool using_standard_tet, double min_singular_value = 0.0001, int max_iter = 1);
    }
    namespace SimpleTriangulation
    {
        bool triangulation(const std::vector<double>& coord_x, const std::vector<double>& coord_y, const std::vector<double>& coord_z, const std::vector<std::vector<int>> &faces, std::vector<std::vector<int>> &new_faces);
        bool sortface_area(std::vector<double>& coord_x, std::vector<double>& coord_y, std::vector<double>& coord_z, std::vector<std::vector<int>> &faces);
        bool merge_triangles(const std::vector<double>& coord_x, const std::vector<double>& coord_y, const std::vector<double>& coord_z, const std::vector<std::vector<int>> &input_faces, std::vector<std::vector<int>> &faces);
        void Laplacian_smooth(std::vector<double> &pts, const std::vector<int>& faces, int iter, double damp_ratio);
    }
    namespace CGALHelper
    {
        enum QualityType {
            SIN_THETA, RHO, ETA
        };
        template <typename T>
        bool line_segm_intersection(std::pair<T, T> &s0, std::pair<T, T> &s1)
        {
            T large_start = s0.first;
            T small_end = s0.second;
            if (s1.first > large_start) large_start = s1.first;
            if (s1.second < small_end) small_end = s1.second;
            return large_start < small_end;
        }
        double compute_cell_volume(const CVec<double, 3> &v1, const CVec<double, 3> &v2, const CVec<double, 3> &v3, const CVec<double, 3> &v4);
        double compute_quad_distortion(std::vector<Vector3d> &quad_array, std::vector<double> &face_distortion);
        void project_points_on_faces(const std::vector<double> &input_pts, const std::vector<int> &faces, std::vector<double> &proj_pts);
        void rays_triangles_intersection(const std::vector<std::pair<CVec<double, 3>, CVec<double, 3>>> &rays, const std::vector<CVec<double, 3>> &tri_pts, const std::vector<unsigned> &tri_faces, const std::vector<std::set<unsigned>> &ignore_faces, std::vector<CVec<double, 3>> &first_intersection_points, std::vector<std::vector<int>> *all_intersection_face_id = NULL);
        void triangulation(std::vector<CVec<double, 3>> &key_pts, std::vector<std::vector<int>> &faces, std::vector<size_t> &begin_pts_index, std::vector<CVec<double, 3>> &end_pts, double tri_area = 0.5, int split_num = 5);
        void triangulation(std::vector<std::vector<int>> &cornerrecord, std::vector<std::vector<CVec<double, 3>>> &key_pts_array, std::vector<std::vector<int>> &faces, std::vector<size_t> &begin_pts_index, std::vector<CVec<double, 3>> &end_pts, double tri_area = 0.5);
        void triangulation(std::vector<std::vector<int>> &cornerrecord, std::vector<std::vector<CVec<double, 3>>> &key_pts_array, std::vector<std::vector<int>> &faces, std::vector<size_t> &begin_pts_index, std::vector<CVec<double, 3>> &end_pts, int& split_num, double tri_area = 0.5);
        void triangulation(std::vector<CVec<double, 3>> &key_pts, std::vector<std::vector<int>> &faces, double tri_area = 0.5);
        void triangulation(std::vector<std::vector<CVec<double, 3>>> &key_pts_array, std::vector<std::vector<int>> &faces, double tri_area = 0.5);
        void localize3dvector(Vector3d input1, Vector3d input2, Vector3d normal, Vector2d& output1, Vector2d& output2);
        void localize3dvector(const std::vector<double> &coord_x, const std::vector<double> &coord_y, const std::vector<double> &coord_z, double normal_x, double normal_y, double normal_z, std::vector<double> &output_x, std::vector<double> &output_y);
        void CrossVectorNormalize(double d1x, double d1y, double d1z, double d2x, double d2y, double d2z, double &rx, double &ry, double &rz);
        bool PointInPolygon(std::vector<Vector3d> &key_pts_array, Vector3d point);
        bool PointInPolygon(std::vector<double> &key_pts_array_x, std::vector<double> &key_pts_array_y, std::vector<double> &key_pts_array_z, double p_x, double p_y, double p_z);
        bool ToLeftTest(CVec<double, 3> p0, CVec<double, 3> p1, CVec<double, 3> p2);
        bool PointOnLineSegment(CVec<double, 3> p0, CVec<double, 3> p1, CVec<double, 3> p2);
        bool LineSegmentIntersectTest(Vector3d u0, Vector3d u1, Vector3d v0, Vector3d v1);
        bool CutLegalTest(std::vector<Vector3d> &key_pts_array, Vector3d begin_point, Vector3d end_point);
        template <class MESH>
        void triangulation(MESH* mesh, double tri_area);
        bool gauss_elimination(std::vector<std::vector<double>> &A);
        void cut_array(const std::vector<double> &input, double sigma_ratio, std::vector<double> &output, double &output_max, double &output_mean);
        template <class MESH>
        void get_outer_boundary_vertices(MESH* mesh, std::vector<Vector3d> &boundary_array)
        {
            std::vector<std::vector<MESH::VertexHandle> >& boundaryvertices = mesh->get_boundaryvertices();
            std::vector<Vector3d> pts;
            {
                pts.clear();
                for (size_t i = 0; i < boundaryvertices[0].size(); i++)
                {
                    Vector3d temp_vec = CVec<double, 3>(boundaryvertices[0][i]->pos);
                    if (i == 0 || (pts.back() - temp_vec).L2Norm() > 0.1 * EPSILON)
                        pts.push_back(temp_vec);
                }
            }
            boundary_array = pts;
        }
        template<typename T>
        int compare(const T& left, const T& right) {
            if (left < right) {
                return -1;
            }
            if (right < left) {
                return 1;
            }
            return 0;
        }
        template<typename MESH>
        void triangulation(MESH* mesh, double tri_area)
        {
            std::vector<CVec<double, 3>> pts;
            std::vector<CVec<double, 3>> key_pts;
            std::vector<std::vector<int>> faces;
            std::vector<int> face;
            double edge_length;
            std::vector<std::vector<MESH::VertexHandle> >& boundaryvertices = mesh->get_boundaryvertices();
            std::vector<std::vector<CVec<double, 3>>> pts_array, key_pts_array;
            for (size_t j = 0; j < boundaryvertices.size(); j++)
            {
                pts.clear();
                for (size_t i = 0; i < boundaryvertices[j].size(); i++)
                {
                    pts.push_back(CVec<double, 3>(boundaryvertices[j][i]->pos));
                }
                pts_array.push_back(pts);
            }
            CVec<double, 3> normal = boundaryvertices[0][0]->normal;
            key_pts_array = pts_array;
            ig::CGALHelper::triangulation(key_pts_array[0], faces, tri_area);
            construct(mesh, key_pts_array[0], faces);
        }
        void graph_coloring(int N, std::vector<std::pair<int, int>>& graph_color_edge, std::vector< std::vector<int> >& v_same_color, std::vector<int>& v_color);
        bool all_pair_shortest_path(int N, const std::vector<std::pair<int, int>>& edges, const std::vector<double> &weight, std::vector<std::vector<double>> &pair_paths);
        template <typename MESH>
        void construct(MESH* mesh, std::vector<CVec<double, 3>> &points, std::vector<std::vector<int>> &faces)
        {
            if (points.size() == 0 || faces.size() == 0)
                return;
            mesh->clear_data();
            for (size_t i = 0; i < points.size(); ++i)
            {
                mesh->insert_vertex(points[i]);
            }
            std::vector<MESH::VertexHandle> s_faceid;
            int face_size = faces[0].size();
            s_faceid.resize(face_size);
            for (size_t i = 0; i < faces.size(); ++i)
            {
                for (int j = 0; j < face_size; ++j)
                {
                    MESH::VertexHandle hv = mesh->get_vertex(faces[i][j]);
                    s_faceid[j] = hv;
                }
                mesh->insert_face(s_faceid);
            }
            mesh->update_mesh();
        }
    }
}
