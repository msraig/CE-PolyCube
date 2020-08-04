#ifndef MESH_POLYCUBE_ALGO_H
#define MESH_POLYCUBE_ALGO_H

#include "MeshDefinition.h"
#include <Eigen/Dense>
#include "Polycube_Deformation.h"
#include "Polycube_Flattening.h"
#include "SmallVec.h"
#include "TetStructure.h"
using ig::Tetrahedron;
using ig::TetStructure;
using ig::TetVertex;
using ig::CVec;

class ThreeCut;
class PolyCubeCutInterface;
class Polycube_Algorithm
{
public:
	Polycube_Algorithm();
	Polycube_Algorithm(VolumeMesh* mesh);
	~Polycube_Algorithm();

	//basic tet operation
	TetStructure<double>* copy_mesh(TetStructure<double>* tet_mesh);
	void get_vert_cell(const TetStructure<double>* tet_mesh, std::vector<CVec<double, 3>> &points, std::vector<unsigned int> &indices);
	void get_boundary_vert_face(TetStructure<double>* tet_mesh, std::vector<double> &points, std::vector<int> &faces, std::vector<int>& s2v);
	
	//preparation function
	void SetMesh(VolumeMesh* mesh)
	{
		mesh_ = mesh;
		pq_deform = new polycube_deformation_interface();
		pq_flattening = new polycube_flattening_interface();
		reset_aLL_state();
	}
	void SetMesh(TetStructure<double>* tet_mesh)
	{
		tet_mesh_ = tet_mesh;
		tet_mesh_unchanged = copy_mesh(tet_mesh);
		pq_deform = new polycube_deformation_interface();
		pq_flattening = new polycube_flattening_interface();
		reset_aLL_state_vtk();
	}
	void SetMeshPolycube(TetStructure<double>* tet_mesh);
	void SetSigmaR(double sigma_r_);
	void SetPolycubeStructure(SurfaceMesh *sm, OpenMesh::VPropHandleT<double> *v_distortion_pt = NULL)
	{
		polycube_structure = sm;
		v_distortion = v_distortion_pt;
	}
	void data_preparation_ovm();
	void reset_aLL_state();
	void reset_aLL_state_vtk();
	void get_boundary_avg_edge_length();
	void get_min_polycube_edge_length();
	void get_bounding_box();
	void initlize_pc();
	void initlize_pc_vtk();

	//cut function
	int auto_cut_batch(double cut_ratio = 0.2, double max_volume = 1.0, int batch_size = 1, double cutting_edge_ratio = 0.8, int distortion_type = 0, int ivf_iter = 10, double min_depth_ratio = 0.5, bool flags_cut_convex_edge = false, bool flag_sample_hex = true, int default_cutnum = -1, double tet_threhold = 1.8, bool random_cut = false);
	int auto_cut_batch_min_max_depth(double cut_ratio = 0.2, double max_volume = 1.0, int batch_size = 1, double cutting_edge_ratio = 0.8, int distortion_type = 0, int ivf_iter = 10, bool min_flag = true, double min_depth_ratio = 0.5, bool flags_cut_convex_edge = false, bool flag_sample_hex = true, int default_cutnum = -1, double tet_threshold = 1.8, bool random_cut = false);
	void randomize_array(std::vector<std::pair<unsigned int, unsigned int>>& pq_edge);
	ig::CVec<double, 3> cut(PolyCubeCutInterface *pci, std::vector<unsigned int> &v_cut_edge, std::vector<double> &cut_depth, double max_volume_ratio = 1.0, bool split_tet = false);
	bool is_edges_cuttable(PolyCubeCutInterface *pci, std::vector<unsigned int> &v_cut_edge, std::vector<double> &cut_depth, double max_volume_ratio = 1.0);
	bool get_available_cut(PolyCubeCutInterface *pci, double min_cut_depth, double max_volume, const std::vector<std::pair<unsigned, unsigned>> &input_edges, const std::set<std::pair<unsigned, unsigned>> &ignore_edge_set, std::vector<unsigned> &avail_edges);
	ig::CVec<double, 3> cut_and_deform(PolyCubeCutInterface *pci, std::vector<unsigned int> &v_cut_edge, std::vector<double> &cut_depth, double sigma_s = 1.0, double max_volume_ratio = 2.0, int ivf_iter = 10);
	void get_cut_depth_min(const std::vector<unsigned int> &cut_edge, std::vector<double> &cut_depth, double cut_ratio);
	void get_cut_depth(const std::vector<unsigned int> &cut_edge, std::vector<double> &cut_depth, double cut_ratio, const std::map<std::pair<unsigned, unsigned>, std::vector<std::pair<unsigned, unsigned>>>& map_associated_edge, int n_vert);
	bool reset_max_cut_depth();
	void get_max_cut_depth(const std::vector<unsigned int> &cut_edge, std::vector<double> &max_cut_depth, double cut_ratio);
	bool dist_rays2mesh(const std::vector<std::pair<CVec<double, 3>, CVec<double, 3>>> &rays, const std::vector<std::set<unsigned>> &ignore_faces, std::vector<double> &distance_array);
	std::pair<int, int> compute_flips();

	//deformation function
	void auto_deformation(double sigma_s = 1, int iter_ivf = 10);
	void do_deformation_constrained_vtk(double sigma_s);
	void do_deformation_constrained_continued();

	//for sample hex distortion
	void do_sample_hex();
	bool construct_original_tet();
	void get_corresponding_point(TetStructure<double>* source, TetStructure<double>* target, const std::vector<ig::CVec<double, 3>>& source_point, std::vector<ig::CVec<double, 3>> &target_point);

	//coloring for acceleration
	void load_vertex_color_vtk(const char* filename);
	bool auto_color();
	void construct_old_new_vert_map();

	//polyline processing
	void deformation_polylines(int max_iter, double sigma_s, bool assign_current = true, bool ns_flag = false, double smooth_factor = 0.1, double threshold = 0.08);
	void feature_aware_deformation(int max_iter, double sigma_s, double threshold = 0.03);
	
	//iversion free functions
	void do_IVF_seamless_moving(int smooth_iter = 20, int max_iter = 40);
	void do_IVF_seamless_moving_part(int smooth_iter = 20, int max_iter = 40, int neighbor_size = 3);
	void do_IVF_seamless_svd(int max_iter = 20);
	void do_IVF_seamless_svd_part(int max_iter = 20, int neighbor = 3);
	void do_IVF_standard(TetStructure<double>* tet_, bool fix_boundary = true, int smooth_iter = 20, int max_iter = 40);
	void do_optimize_polycube_PC(bool fix_boundary = true, int max_smooth_iter = 0, int max_iter = 20);
	
	//extract initial polycube boundary
	bool do_Extract_polycube_vtk(int distortion_type = 0, bool sample_hex = false);
	void do_construct_polycube_structure();

	//flattening function
	void flattening_prepare_vtk();
	bool do_flattening_constrained_vtk(int distortion_type = 0, double default_cube_length = -1.0, double lambda = 0.5, bool add_constraint = false, int IVF_polycube_type = 0, int set_min_diff = false, bool int_solver_flag = true, double offset = 0.001, bool double_cube_length_flag = false);
	
	//feature operation
	void repair_features(bool polyline_guidance = false, bool ns_flag = false, double smooth_factor = 0.1, double ratio_threshold = 0.08);
	bool get_cut_edge_map();
	bool refine_polycube_feature_edge_array(const std::map<int, int> &vertmap);

	//load and save function
	bool load_chart_label(const char* filename);
	void load_equal_faces(const char* filename);
	void load_equal_faces(const std::vector<int> &v_cut_type, const std::vector<int> &v_cut_section_chart, const std::vector<int> &v_seam_edge_vertex_id, const std::vector<std::vector<int>>& vv_congruent_face, int nv, std::vector<ThreeCut> &three_cut_info);//from cut result
	bool load_feature_edges(const char* filename);
	void save_polyline_vtkformat(const char* filename);
	void save_polycube_vtk(const char* filename);
	void save_polycube_scaled_vtk(const char* filename);
	void save_tet_vtk(const char* filename);
	void save_hex_length(const char* filename);
	void save_feature_edges_ovm(const char* filename);
	void save_feature_edges_ovm_vtkformat(const char* filename);
	bool save_polycube_feature_edges_pfeformat(const char* filename);
	void save_fix_chart_label(const char* filename);
	void save_equal_faces(const char* filename);
	void save_tet_mesh(const char* filename);
	void save_hexex(const char* filename, const TetStructure<double>* oritet, const TetStructure<double>* pqtet);
	void save_polycube_tetformat(const char* filename);
	void save_polycube_para(const char* filename);
	void save_flip_info(const char* filename);

	//data
	OpenMesh::VPropHandleT<double> *v_distortion;
	VolumeMesh* mesh_;
	TetStructure<double>* tet_mesh_;
	TetStructure<double>* tet_mesh_polycube_;
	bool polycube_exist_flag;
	TetStructure<double>* tet_mesh_unchanged;
	TetStructure<double>* tet_mesh_polycube_unchanged;
	SurfaceMesh *polycube_structure;
	std::vector<double> polycube_structure_v_distortion;
	std::map<std::pair<int, int>, double> polycube_structure_e_distortion;
	std::vector<std::pair<unsigned int, unsigned int>> polycube_edge;
	std::vector<std::pair<int, int>> polycube_edge_chart;
	std::vector<std::pair<int, int>> polycube_edge_label;
	std::vector<std::pair<unsigned int, unsigned int>> sorted_polycube_edge;
	std::vector<double> sorted_polycube_edge_distortion;
	ig::CVec<double, 3> init_distortion;
	ig::CVec<double, 3> init_volume_distortion;
	ig::CVec<double, 3> bb_Min, bb_Max;
	double avg_boundary_edge_length;
	double min_polycube_edge_length;  //now set as min{bb / SPLIT, avg_boundary_edge_length}
	double sigma_r;
	double final_cube_length;

	std::vector<double> polycube_edge_distortion;
	std::vector<double> polycube_edge_hex_distortion;
	std::vector<double> polycube_edge_length;
	std::vector<double> polycube_edge_neighbor_min_length;
	std::map<std::pair<unsigned, unsigned>, unsigned> polycube_edge_to_idx_map;
	std::vector<double> polycube_edge_max_cut_depth_one_cut; //max cut depth of with no self-intersection
	std::vector<double> polycube_edge_max_cut_depth_three_cut;

	std::vector<double> coord_ori;
	std::vector<Eigen::Matrix3d> R;
	std::vector<Eigen::Matrix3d> PT;
	std::vector<OpenVolumeMesh::Geometry::Vec4i> cv_id;
	std::vector<double> cell_volume;
	double avg_cell_volume;
	std::vector<OpenVolumeMesh::Geometry::Vec3i> bfv_id;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> bf_n;
	std::vector<int> fix_cell_boundary;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> src_pos;
	std::vector<int> fix_charts;
	std::vector<int> fix_labels;

	//equal face info
	std::vector<std::vector<std::vector<int>>> equal_faces_array;
	std::vector<int> cut_types;
	std::vector<int> common_verts;
	std::vector<std::pair<int, int>> chart_pair, chart_pair_neighbor;
	std::vector<int> three_cut_type_array;
	std::vector<int> three_cut_common_vert;
	std::vector<std::array<int, 3>> temp_three_cut_adjacent_one_cut_type;
	std::vector<std::array<int, 3>> three_cut_adjacent_one_cut_index;
	std::vector<std::vector<std::array<unsigned int, 3>>> three_cut_vert;

	//feature part
	std::vector<std::pair<int, int>> init_feature_edge_array;
	std::vector<std::pair<unsigned, unsigned>> init_feature_polycube_edge_array_ori;
	std::vector<std::pair<unsigned, unsigned>> init_feature_polycube_edge_array;
	std::vector<std::pair<unsigned, unsigned>> init_feature_polycube_edge_array_sorted; //each pair is sorted
	std::map<std::pair<unsigned int, unsigned int>, std::vector<std::pair<unsigned int, unsigned int>>> cut_edge_o2n;
	std::vector<bool> *feature_edge_flag;

	//cut part
	std::vector<CVec<double, 3>> cut_input_point;
	std::vector<unsigned int> cut_input_face;
	std::vector<int> cut_input_chart;
	std::vector<int> cut_input_label;
	std::vector<double> chart_mean_value_ori; //original chart_mean value
	//old new vert map
	std::vector<int> o2n; //old to new
	std::vector<std::array<int, 3>> n2o; //new to old
	bool ovm_color_status;
	bool data_prepare_status;
	bool feature_status;

	polycube_deformation_interface *pq_deform;
	polycube_flattening_interface *pq_flattening;
};

#endif