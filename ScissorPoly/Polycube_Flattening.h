#ifndef TETRAHEDRAL_DEFORMATION_H
#define TETRAHEDRAL_DEFORMATION_H
#include "MeshDefinition.h"
#include "Defined_Tensor.h"
#include <Eigen\Dense>
//#include <hash_map>
#include <unordered_map>
#include "TetStructure.h"
#define ROUNDING_TH 1e-6
using ig::Tetrahedron;
using ig::TetStructure;
using ig::TetVertex;
using ig::CVec;
class polycube_flattening_interface
{
public:
	polycube_flattening_interface();
	~polycube_flattening_interface();
	bool get_prepare_ok() { return prepare_ok; };
	void set_prepare_ok(bool ok) { prepare_ok = ok; };
	void set_hex_meshing_flag(bool flag) { hex_meshing_flag = flag; };
	void set_bd_face_chart_label(const std::vector<CVec<double, 3>> &bd_pts, const std::vector<unsigned> &bd_face, const std::vector<int> &bd_chart, const std::vector<int> &bd_label);
	void get_bd_face_chart_label(TetStructure<double> *tet_mesh_);
	void prepare_for_deformation(TetStructure<double> *tet_mesh_, const std::vector<double> &x_ori, const std::vector<double> &y_ori, const std::vector<double> &z_ori);
	void load_deformation_result(const std::vector<double> &coord_ori, TetStructure<double> *tet_mesh_);
	void set_cube_length(double d) { cube_len = d; }
	double get_cube_length() { return cube_len; }
	ig::CVec<double, 3> compute_volumetric_distortion(TetStructure<double> *tet_mesh_);
	ig::CVec<double, 3> compute_distortion(TetStructure<double> *tet_mesh_);
	void get_chart_distortion();
	void get_polycube_edge_layer_distortion(int n_layer, TetStructure<double> *tet_mesh_);
	void get_edge_length(const std::vector<std::pair<int, int>> &edgearray, std::vector<double> &edge_length);
	void compute_triangle_area();
	void set_equal_faces(const std::vector<std::vector<std::vector<int>>> &faces_array, const std::vector<int>& common_verts_idx, const std::vector<int>& cut_types, int n_vert, const std::vector<std::pair<int, int>> &chart_pair, const std::vector<std::pair<int, int>> &chart_pair_neighbor, const std::vector<int> &three_cut_common_vert, const std::vector<std::vector<std::array<unsigned int, 3>>> &three_cut_vert, const std::vector<std::array<int, 3>> &three_cut_adjacent_one_cut_index);
	void set_sigma_r(double r) { sigma_r = r; }
	void set_sigma_s(double s) { sigma_s = s; }
	void load_boundary_face_label(const std::vector<int> &label, const std::vector<int> &chart, TetStructure<double> *tet_mesh_);
	//feature part
	void load_feature_edges_vtk(const std::vector<std::pair<int, int>> &feature_pairs);
	bool feature_extraction(double select_ratio = 0.7, int min_feature_length = 3);
	bool repair_chartlabel_feature(int max_line_length = 8, int max_area_size = 20);
	bool save_chartlabel(const char* filename);
	void save_feature_edges_vtk_feaformat(const char* filename);
	void save_feature_edges_vtk_tfeformat(const char* filename, TetStructure<double>* tet_mesh_ori);
	void save_feature_edges_vtk_vtkformat(const char* filename);
	void save_error_feature_ptid_ptsformat(const char* filename);
	void save_feature_long_edges_vtk_vtkformat(const char* filename);
	void save_feature_long_edges_vtk_psfeformat(const char* filename); //include segment info
	void find_all_corner(TetStructure<double>* tet_mesh_);
	std::vector<double> get_chart_mean_value() { return chart_mean_value; }
	bool find_all_chart_value_equal_face_rounding(TetStructure<double> *tet_mesh_, double offset = 0.001);
	bool find_all_chart_value_equal_face_mine(TetStructure<double> *tet_mesh_, double offset = 0.001);
	bool find_all_chart_value_equal_face_mine(TetStructure<double> *tet_mesh_, const std::vector<double>& chart_mean_value_ori, double lambda = 0.5, bool add_constraint = false, bool set_min_diff = false, double offset = 0.001);
	bool find_all_chart_value_preprocessing(TetStructure<double> *tet_mesh_);
	//void find_all_chart_value_equal_face_mine(TetStructure<double> *tet_mesh_);
	void update_updownchart(TetStructure<double> *tet_mesh_);
	void save_polycube_para(const char* filename, TetStructure<double>* tet_mesh_);
	bool deform_ARAP_polycube_equal_face(TetStructure<double>* tet_mesh_, bool int_solver_flag = true, bool add_constraint = false, bool set_min_diff = false, double offset = 0.001, bool double_cube_length_flag = false); //no ori chart version
	bool deform_ARAP_polycube_equal_face(TetStructure<double>* tet_mesh_, const std::vector<double>& chart_mean_value_ori, double lambda = 0.5, bool add_constraint = false, bool set_min_diff = false, double offset = 0.001, bool double_cube_length_flag = false);
	void boundary_mapping_polycube_equal_face(TetStructure<double>* tet_mesh_); //using IVF
	//polycube varibale
	std::vector<std::vector<int>> edge_with_same_label;
	std::vector<std::vector<int>> key_edges;
	std::vector<std::vector<int>> polycube_edges; //0: first endpoint, 1: second endpoint, 2: left chart, 3: right chart
	std::vector<std::vector<int>> polycube_short_edges; //same size as polycube_edges, store the short edge idx of each long polycube edges
	std::vector<double> polycube_edge_distortion_iso;
	std::vector<double> polycube_edge_distortion_conf;
	std::vector<double> polycube_edge_distortion_vol;
	std::vector<double> polycube_edge_face_distortion;
	std::vector<double> polycube_edge_layer_distortion;
	std::vector<double> polycube_edge_length;
	std::vector<int> polycube_edge_idx2same_label_idx;
	//for feature edges
	std::vector<bool> feature_edge_flag_ori; //size: id2edge
	std::vector<bool> feature_edge_flag_final; //size: id2edge
	std::vector<int> feature_polycube_edge; //size: featured polycube edges, include those with only a segment on the polycube edge
	std::vector<std::pair<double, double>> feature_polycube_edge_segm; //size: feature polycube edges, segmentation of polycube edges, if the whole edge is occupied, the value would be [0,1]
	std::vector<std::vector<int>> group_feature_edge; //same as feature_polycube_edge, work for vtk
	std::vector<int> group_feature_edge_start; //start point of group_feature_edge[i], same size
	std::vector<std::vector<int>> feature_v2v;
	std::vector<std::vector<int>> feature_v2e;
	std::vector<int> valence_three_error_feature_pt_id;
	std::vector<std::vector<int>> pqedge_v2v; //polycube edge vertex to vertex, used for feature
	std::vector<std::vector<int>> pqedge_v2e;
	void assign_pos_mesh(TetStructure<double>* tet_mesh_, bool r_order = false);
	void get_coord(const std::vector<int>& all_vert, std::vector<double> &coord_x, std::vector<double> &coord_y, std::vector<double> &coord_z);
	std::vector<OpenVolumeMesh::Geometry::Vec3i> get_vertex_type() { return vertex_type; }
	std::vector<int>& get_polycube_chart_label()
	{
		return polycube_chart_label;
	}
private:
	void reorder_edge(std::vector<std::pair<int, int>> &one_edge, int start);
public:
	bool hex_meshing_flag;
	std::vector<CVec<double, 3>> bd_pts;
	std::vector<unsigned int> bd_faces;
	std::vector<int> bd_chart;
	std::vector<int> bd_label;
	bool update_updownchart_flag;
private:
	bool prepare_ok;
	//basic data
	std::vector< std::vector<std::vector<double> > > S_v;
	std::vector< std::vector<Eigen::Matrix3d > > S;
	std::vector< double > cell_volume;
	std::vector< std::vector<int> > cell_vertex;
	std::vector< Eigen::Matrix3d > cell_S;
	std::vector< std::vector<std::vector<int> > > cell_vertex_vertex;
	std::vector< std::vector<int>> vertex_cell;
	std::vector< std::vector<std::vector<int>>> vertex_cell_vertex;
	std::vector< std::vector<std::vector<double> > > vcv_S;
	std::vector< std::vector<std::vector<double> > > init_vcv_A;
	double mu_A;
	std::vector<double> dpx; std::vector<double> dpy; std::vector<double> dpz;
	std::vector<double> dpx_ori; std::vector<double> dpy_ori; std::vector<double> dpz_ori;
	int number_of_color; int max_vc_size; double alpha;
	std::vector< std::vector<int> > vertex_diff_color;
	std::vector<double > vertex_move_d;
	//deformation related
	std::vector<int> is_deformation_constraint;
	std::vector<int> is_polycube_handles;
	std::vector<int> deformation_handle_id;
	std::vector<int> deformation_v_id;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> deformation_new_p;
	std::vector<int> deformation_ok;
	std::vector<int> handles_region_id;
	std::vector<std::vector<int> > region_handles;
	std::vector<int> rotate_v_id;
	std::vector<double> rotate_radius;
	std::vector<double> rotate_angle;
	std::vector<int> rest_v_id;
	int moving_handle; int moving_handle_id;
	bool is_auto_deformation;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> error_color;
	std::vector<double> ori_edge_len;
	std::vector<OpenVolumeMesh::CellHandle> flipped_cell;
	std::vector<int> con_dis_above_th_vertex;
	std::vector<int> iso_dis_above_th_vertex;
	std::vector<int> dis_above_th_cell;
	double dis_th; double bound_c; double bound_v;
	bool compute_eps(VolumeMesh* mesh_);
	double eps;
	double asimptotic_epsilon; // Asymptotic value of h(sigma) when sigma tends to -inf
	double EpsSafetyFactor;      // "Safety" factor for epsilon
	double meps;         // Machine epsilon
	double EpsilonEfectivo; // Effective epsilon
	double umbral_factor; // Effective epsilon
	std::vector<OpenVolumeMesh::Geometry::Vec3i> bfv_id;
	std::vector< std::vector<int> > bvf_id;
	std::vector< std::vector<int> > bvv_id;
	std::vector< OpenVolumeMesh::Geometry::Vec3d > target_bfn;
	std::vector<int> bf_chart;
	std::vector<std::vector<int>> polycube_chart; //each chart which stores the faces
	std::vector<int> polycube_chart_label;//0: x, 1: -x, 2: y, 3: -y, 4: z, 5: -z
	std::vector<double> chart_mean_value; //the mean value of each chart
	std::vector<int> up_chart_id; std::vector<int> down_chart_id; std::vector<double> diff_up_down;
	std::vector<std::vector<int> > bef_id;
	std::vector<OpenVolumeMesh::Geometry::Vec3i> vertex_type;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> src_pos;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> init_pos;
	//for equal faces
	std::vector<int> vert_update_type; // 0 : no cut, 1: update part, 2: no update part
	//std::vector<int> vert_cut_type; // 0: no cut, for the update part of the cut, set vert_cut_type to be vert type.
	std::vector<int> vert_belong_cut_array;
	std::vector<std::pair<int, int>> cut_to_chart_pair;
	std::vector<std::pair<int, int>> cut_to_chart_pair_neighbor;
	std::vector<int> cut_types; // totally 12 types
	std::vector<int> cut_common_verts_idx;
	std::map<int, int> vert_pairs_map;
	std::vector<std::vector<int>> equal_triangles;
	std::vector<int> three_cut_common_vert;
	std::vector<std::vector<std::array<unsigned int, 3>>> three_cut_vert;
	std::vector<int> three_cut_vert_flag;
	std::vector<std::array<int, 3>> three_cut_adjacent_one_cut_index;
	//boundary edge information for tet mesh
	std::map<std::pair<int, int>, int> edge2id; //for one edge, first idx less than the second, work for vtk
	std::vector<std::pair<int, int>> id2edge;  //work for vtk
	const double type_matrix[12][3][3] = {
		{ {0, 1, 0}, { 1,0,0 }, { 0,0,1 } },
		{ {0, -1, 0}, { -1,0,0 }, { 0,0,1 } },
		{ {0, -1, 0}, { -1,0,0 }, { 0,0,1 } },
		{ {0, 1, 0}, { 1,0,0 }, { 0,0,1 } },
		{ {0, 0, 1}, { 0,1,0 }, { 1,0,0 } },
		{ {0, 0, -1}, { 0,1,0 }, { -1,0,0 } },
		{ {0, 0, -1}, { 0,1,0 }, { -1,0,0 } },
		{ {0, 0, 1}, { 0,1,0 }, { 1,0,0 } },
		{ {1, 0, 0}, { 0,0,1 }, { 0,1,0 } },
		{ {1, 0, 0}, { 0,0,-1 }, { 0,-1,0 } },
		{ {1, 0, 0}, { 0,0,-1 }, { 0,-1,0 } },
		{ {1, 0, 0}, { 0,0,1 }, { 0,1,0 } } };
	std::vector<double> all_iso_d; std::vector<double> all_con_d; std::vector<double> all_vol_d;
	std::vector<double> chart_iso_d; //store the isometric distortion of each chart
	std::vector<double> chart_con_d;
	std::vector<double> chart_vol_d;
	void build_AABB_Tree();
	CGAL_AABB_Tree* AABB_Tree;
	std::vector<CGAL_3_Triangle> triangle_vectors;
	std::vector<int> is_boundary_v;
	void project_on_ref_mesh(OpenVolumeMesh::Geometry::Vec3d& p);
	void project_on_ref_mesh(double& px, double& py, double& pz);
	bool tan_smooth;
	double avg_boundary_edge_length; int boundary_face_number;
	double sigma_s; double sigma_r;
	double cube_len;
};
#endif // !TETRAHEDRAL_DEFORMATION_H