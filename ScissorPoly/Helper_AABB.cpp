#pragma warning( disable : 4477 4018 4267 4244 4838)
#include <iostream>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include "Helper.h"
using ig::CVec;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Ray_3 Ray;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;
typedef Tree::Primitive_id Primitive_id;
struct Skip {
    std::set<face_descriptor> f_set;
    Skip()
    {
    }
    Skip(const face_descriptor face)
    {
        f_set.insert(face);
    }
    Skip(const face_descriptor face1, const face_descriptor face2)
    {
        f_set.insert(face1);
        f_set.insert(face2);
    }
    Skip(const std::set<face_descriptor> &face_set)
    {
        f_set = face_set;
    }
    bool operator()(const face_descriptor& t) const
    {
        auto search = f_set.find(t);
        if (search != f_set.end())
            return true;
        return false;
    }
};
struct Intersection_visitor {
    typedef void result_type;
    void operator()(const Point& p) const
    {
        std::cout << p << std::endl;
    }
};
void ig::CGALHelper::rays_triangles_intersection(const std::vector<std::pair<CVec<double, 3>, CVec<double, 3>>>& rays, const std::vector<CVec<double, 3>>& tri_pts, const std::vector<unsigned>& tri_faces, const std::vector<std::set<unsigned>>& ignore_faces, std::vector<CVec<double, 3>>& first_intersection_points, std::vector<std::vector<int>> *all_intersection_face_id)
{
    if (!ignore_faces.empty())
    {
        assert(rays.size() == ignore_faces.size());
    }
    if (all_intersection_face_id != NULL)
    {
        all_intersection_face_id->clear();
        all_intersection_face_id->resize(rays.size());
    }
    first_intersection_points.clear();
    Mesh mesh;
    std::vector<vertex_descriptor> v_handles;
    for (size_t i = 0; i < tri_pts.size(); i++)
    {
        v_handles.push_back(mesh.add_vertex(K::Point_3(tri_pts[i][0], tri_pts[i][1], tri_pts[i][2])));
    }
    std::vector<face_descriptor> f_array;
    int n_faces = tri_faces.size() / 3;
    for (size_t i = 0; i < n_faces; i++)
    {
        f_array.push_back(mesh.add_face(v_handles[tri_faces[3 * i]], v_handles[tri_faces[3 * i + 1]], v_handles[tri_faces[3 * i + 2]]));
    }
    Tree tree(faces(mesh).first, faces(mesh).second, mesh);
    for (size_t i = 0; i < rays.size(); i++)
    {
        Point p(rays[i].first[0], rays[i].first[1], rays[i].first[2]);
        Vector v(rays[i].second[0], rays[i].second[1], rays[i].second[2]);
        Ray ray(p, v);
        std::set<face_descriptor> face_set;
        if (!ignore_faces.empty())
        {
            for (auto it = ignore_faces[i].begin(); it != ignore_faces[i].end(); ++it)
            {
                face_set.insert(f_array[*it]);
            }
        }
        Skip skip(face_set);
        Ray_intersection intersection = tree.first_intersection(ray, skip);
        if (intersection) {
            if (boost::get<Point>(&(intersection->first))) {
                const Point* p = boost::get<Point>(&(intersection->first));
                first_intersection_points.push_back(CVec<double, 3>((*p)[0], (*p)[1], (*p)[2]));
            }
            else
            {
                std::cout << intersection->first << std::endl;
                first_intersection_points.push_back(first_intersection_points.back());
            }
        }
        if (all_intersection_face_id != NULL)
        {
            std::list<Primitive_id> primitives;
            tree.all_intersected_primitives(ray, std::back_inserter(primitives));
            for (auto p : primitives)
            {
                (*all_intersection_face_id)[i].push_back(p);
            }
        }
    }
}
