#ifndef _POLYCUBE_CUT_HELPER_H_
#define _POLYCUBE_CUT_HELPER_H_

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <tetgen/tetgen.h>

#include <fstream>

constexpr auto EPSILON = 1e-12;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3> Polyhedron;
typedef Polyhedron::Edge_iterator Edge_iterator;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef Polyhedron::Vertex_handle Vertex_handle;
typedef Polyhedron::Face_handle Face_handle;
typedef K::Vector_3 Vector;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef K::Triangle_3 Triangle;
typedef K::Plane_3 Plane;

// quadrilateral facet
class CutFacet
{
public:
    std::array<unsigned int, 4> a_vertex; // indices of four corners
    std::array<Point, 4> a_point; // positions of four corners
	std::array<Vector, 4> a_edge_vector; // vectors of four edges
	std::array<double, 4> a_edge_length_2; // squared edge lengths of four edges
    std::array<bool, 2> a_is_interior;
    bool is_convex; // true if this is a convex cut
    Vector normal;
	Plane plane;
};

// C++-style tetgenio structure, assisting in constructing PLC

class PolygonHelper
{
public:
    std::vector<unsigned int> vertexlist;
};

class FacetHelper
{
public:
    int chart;
    std::vector<PolygonHelper> polygonlist;
};

class TetgenHelper
{
public:
    std::vector<Point> pointlist;
    std::vector<FacetHelper> facetlist;
};

class PolyCubeCutHelper
{
public:
    enum LOCATION_RESULT
    {
        ON_INTERIOR_EDGE_0,
        ON_INTERIOR_EDGE_1,
        ON_SEAM_EDGE,
		BORDER_OR_INSIDE,
        OUTSIDE
    };

	static void off2cgal(std::string filename, Polyhedron *polymesh);
    static void assign_cgal_id(Polyhedron *polymesh);
	static void cgal2tetgen_helper(const Polyhedron &polymesh, const std::vector<int> &v_face2chart, TetgenHelper *tetgen_helper);
	static void tetgen_helper2tetgen(const TetgenHelper &tetgen_helper, tetgenio *tetin);
	static void tetgen2vtk(tetgenio *tetout, std::string filename);
    static LOCATION_RESULT locate_point_in_cut_facet(const Point &p, const CutFacet &facet);
    static bool point_inside_quadrilateral(const Point &p, const std::array<Point, 4> &v_quad);

private:
	PolyCubeCutHelper() {}
	~PolyCubeCutHelper() {}

    static bool point_on_segment(const Point &q, const Point &pa, const Point &pb);
};

#endif