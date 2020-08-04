#ifndef DEFINED_TENSOR_H
#define DEFINED_TENSOR_H

#define BOOST_PARAMETER_MAX_ARITY 24

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include "UtilityDefinitions.h"
#include "MeshDefinition.h"
//#include "..\VolumeMeshViewer\MeshDefinition.h"
#include "ANN\ANN.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_segment_primitive.h>

typedef CGAL::Simple_cartesian<double> CGAL_Simple_Cartesian_Kernal;
typedef CGAL_Simple_Cartesian_Kernal::Point_3 CGAL_double_3_Point;

typedef CGAL_Simple_Cartesian_Kernal::Triangle_3 CGAL_3_Triangle;
typedef std::vector<CGAL_3_Triangle>::iterator CGAL_Triangle_Iterator;
typedef CGAL::AABB_triangle_primitive<CGAL_Simple_Cartesian_Kernal, CGAL_Triangle_Iterator> CGAL_AABB_Triangle_Primitive;
typedef CGAL::AABB_traits<CGAL_Simple_Cartesian_Kernal, CGAL_AABB_Triangle_Primitive> CGAL_AABB_triangle_traits;
typedef CGAL::AABB_tree<CGAL_AABB_triangle_traits> CGAL_AABB_Tree;

typedef CGAL_Simple_Cartesian_Kernal::Segment_3 CGAL_3_Segment;
typedef std::vector<CGAL_3_Segment>::iterator CGAL_Segment_Iterator;
typedef CGAL::AABB_segment_primitive<CGAL_Simple_Cartesian_Kernal, CGAL_Segment_Iterator> CGAL_AABB_Segment_Primitive;
typedef CGAL::AABB_traits<CGAL_Simple_Cartesian_Kernal, CGAL_AABB_Segment_Primitive> CGAL_AABB_Segment_traits;
typedef CGAL::AABB_tree<CGAL_AABB_Segment_traits> CGAL_AABB_Segment_Tree;

class defined_tensor
{
public:
	defined_tensor();
	~defined_tensor();
	void build_AABB_Tree_using_ref_mesh(SurfaceMesh* mesh_);
	void build_AABB_Tree_using_ref_mesh(VolumeMesh* mesh_);
	void build_AABB_Segment_Tree_using_ref_mesh(VolumeMesh* mesh_, const std::vector<int>& edge_feature);
	void project_on_ref_mesh(OpenVolumeMesh::Geometry::Vec3d& p); //will change the p's value
	void project_on_ref_mesh_with_guess_face(OpenVolumeMesh::Geometry::Vec3d& p, int& guess_face_id, const double& radius);
	void project_on_ref_mesh_with_metric(OpenVolumeMesh::Geometry::Vec3d& p, OpenVolumeMesh::Geometry::Vec6d& M); //will change the p's value
	void project_on_ref_edge(OpenVolumeMesh::Geometry::Vec3d& p); //will change the p's value
	void project_for_interior(OpenVolumeMesh::Geometry::Vec3d& p, OpenVolumeMesh::Geometry::Vec6d& M);

	void set_tensor_name(int tn){ m_tensor_name = tn;};
	int get_tensor_name() {return m_tensor_name;};
	void get_tensor(const OpenVolumeMesh::Geometry::Vec3d& p, OpenVolumeMesh::Geometry::Vec6d& M);

	//propagate the surface metric field to the interior
	void smmoth_metric_to_interior(VolumeMesh* mesh_);
	void save_tet_vertex_metric(const char* filename);
	void load_tet_vertex_metric(const char* filename);
	void load_ref_tet_mesh(VolumeMesh* mesh_);

	void set_info_for_project(std::vector<std::vector<int>>& opp_face_id_,
		std::vector<std::vector<OpenVolumeMesh::Geometry::Vec3d>>& face_pos_,
		std::vector<int>& visited_ref_face_id_, std::vector<int>& visited_ref_face_)
	{
		opp_face_id_ = opp_face_id;
		face_pos_ = face_pos;
		visited_ref_face_id_ = visited_ref_face_id;
		visited_ref_face_ = visited_ref_face;
	}

private:
	int m_tensor_name;
	CGAL_AABB_Tree* AABB_Tree;
	std::vector<CGAL_3_Triangle> triangle_vectors;
	CGAL_AABB_Segment_Tree* AABB_Segment_Tree;
	std::vector<std::vector<int>> opp_face_id; 
	std::vector<std::vector<OpenVolumeMesh::Geometry::Vec3d>> face_pos;
	std::vector<int> visited_ref_face_id; std::vector<int> visited_ref_face;

	std::vector<OpenVolumeMesh::Geometry::Vec6d> tet_vertex_metric;
	std::vector<OpenVolumeMesh::Geometry::Vec6d> bm_metric;
	VolumeMesh* ref_mesh_; SurfaceMesh* bm;
	void construct_boundary_mesh();
	void build_ANN_KD_tree();
	ANNkd_tree* kdTree;
	int find_nearest_cell_id(const OpenVolumeMesh::Geometry::Vec3d& p);
	bool check_in_tet(const OpenVolumeMesh::Geometry::Vec3d& p, const int& cel_id);
	void find_cell_bay_cor(const OpenVolumeMesh::Geometry::Vec3d& p,const int& cell_id, OpenVolumeMesh::Geometry::Vec4i& v_id, OpenVolumeMesh::Geometry::Vec4d& b);
	void save_boundary_mesh_metric();
};

struct local_frame
{
	local_frame()
		:e_x(OpenVolumeMesh::Geometry::Vec3d(1,0,0)), e_y(OpenVolumeMesh::Geometry::Vec3d(0,1,0)), n(OpenVolumeMesh::Geometry::Vec3d(0,0,1))
	{}
	~local_frame(){}

	void find_e_x_y()
	{
		if(std::abs(n[2]) >= std::abs(n[1]) && std::abs(n[2]) >= std::abs(n[0]))
		{
			e_x[0] = 1.0; e_x[1] = 1.0; e_x[2] = ( -n[0] - n[1] ) / n[2];
		}
		else if (std::abs(n[1]) >= std::abs(n[2]) && std::abs(n[1]) >= std::abs(n[0]))
		{
			e_x[0] = 1.0; e_x[2] = 1.0; e_x[1] = ( -n[0] - n[2] ) / n[1];
		}
		else 
		{
			e_x[1] = 1.0; e_x[2] = 1.0; e_x[0] = ( -n[2] - n[1] ) / n[0];
		}
		e_x.normalize();
		e_y = OpenVolumeMesh::Geometry::cross(n, e_x);
	}

	OpenVolumeMesh::Geometry::Vec3d e_x;
	OpenVolumeMesh::Geometry::Vec3d e_y;
	OpenVolumeMesh::Geometry::Vec3d n;
};

#endif // !DEFINED_TENSOR_H
