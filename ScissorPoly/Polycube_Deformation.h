#ifndef POLYCUBE_DEFORMATION_H
#define POLYCUBE_DEFORMATION_H
#include "MeshDefinition.h"
#include <Eigen\Dense>
#include <unordered_map>
#include "TetStructure.h"
using ig::Tetrahedron;
using ig::TetStructure;
using ig::TetVertex;
class polycube_deformation_interface
{
public:
	polycube_deformation_interface();
	~polycube_deformation_interface();
	bool get_prepare_ok(){ return prepare_ok; };
	void set_prepare_ok(bool ok){ prepare_ok = ok; };
	void prepare_for_deformation(VolumeMesh* mesh_);
	void prepare_for_deformation(TetStructure<double>* tet_mesh_);
	void load_ori_tet(TetStructure<double>* tet_mesh_);
	//compute distortion
	void compute_distortion(VolumeMesh* mesh_);
	ig::CVec<double, 3> compute_distortion(TetStructure<double>* tet_mesh_);
	ig::CVec<double, 3> compute_volumetric_distortion(TetStructure<double>* tet_mesh_);
	std::pair<double, double> compute_chart_normal_diff();
	void update_face_center();
	//deformation
	void exp_mips_deformation_refine_polycube_omp_feature_nf(int max_iter, int iter2, double energy_power, double angle_area_ratio, VolumeMesh *mesh_, bool use_xyz = true, bool use_RGNF = true, bool use_half_sigma_s = false, bool feature_flag = false, bool boundary_global_opt = false);
	void deformation_polylines(VolumeMesh* mesh_, bool assign_current = true, bool ns_flag = false, double smooth_factor = 0.1);
	bool check_polycube_planarity(double threshold = 0.03);
	bool check_polyline_orthogonality(double threhold = 0.03);
	void exp_mips_deformation_refine_polycube_omp_constrained(int max_iter, int iter2, double energy_power, double angle_area_ratio, int nv, bool use_xyz = true, bool use_RGNF = true, bool use_half_sigma_s = false, bool change_mu = true);
	void optimize_all_verts_nf();
	void optimize_all_polylines(bool assign_current = true);
	void set_sigma_r(double r) { sigma_r = r; }
	void set_sigma_s(double s) { sigma_s = s; }
	void set_fixed_chart_label(std::vector<int> &chart, std::vector<int> &label, bool is_ovm = true);
	void set_equal_faces(const std::vector<std::vector<std::vector<int>>> &faces_array, const std::vector<int>& common_verts_idx, const std::vector<int>& cut_types, std::vector<std::pair<int, int>> &chart_pair, std::vector<int> &three_cut_type, std::vector<int> &three_cut_common_vert, std::vector<std::vector<std::array<unsigned int, 3>>> &three_cut_vert, const std::vector<std::array<int, 3>> &temp_three_cut_adjacent_one_cut_type);
	void adjust_orientation(VolumeMesh* mesh_, bool feature_flag = false);
	void adjust_orientation_feature(VolumeMesh* mesh_, bool polyline_flag = false);
	void adjust_orientation(TetStructure<double>* tet_mesh_);
	void assign_pos_mesh(VolumeMesh* mesh_, bool r_order = false);
	void assign_pos_feature(std::vector<double> &fx, std::vector<double> &fy, std::vector<double> &fz, bool r_order = false);
	void assign_pos_mesh(TetStructure<double>* tet_mesh_, bool r_order = false);
	void set_vertex_color(int n_color, const std::vector<int>& v_color, int nv);
	void set_feature_edges_ovm(const std::vector<std::pair<int, int>> &edge_pairs, VolumeMesh * mesh_);
	void build_feature_edges_connectivity(VolumeMesh * mesh_, bool refine_long_edge = false, int min_edge_size = 5);
	void repair_features(VolumeMesh * mesh_, bool polyline_guidance = false);
	void feature_long_edge_redefinition(int min_edge_size = 5);
	void corner_redefinition_polyline();
	bool deformation_prepare_OK;
	bool polyline_ns_flag;
private:
	//feature part
	bool check_and_repair_feature(bool repair_status, bool polyline_guidance = false);
	void refine_corner_and_long_edge(int min_edge_size);
	
	bool compute_exp_misp_energy_refine_polycube(const int& vc_size, const double& old_e, double& new_e,
		const double* posx, const double* posy, const double* posz,
		const std::vector<std::vector<double> >& vc_S, double* exp_vec,
		const double& npx, const double& npy, const double& npz,
		double alpha, double beta,
		const double& ga, const int& vf_size, 
		const double* vf_px, const double* vf_py, const double* vf_pz, const double* mu);
	bool compute_exp_misp_energy_refine_polycube(const int& vc_size, const double& old_e, double& new_e,
		const double* posx, const double* posy, const double* posz,
		const std::vector<std::vector<double> >& vc_S, double* exp_vec,
		const double& npx, const double& npy, const double& npz,
		double alpha, double beta,
		const double& ga, const int& vf_size,
		const double* vf_px, const double* vf_py, const double* vf_pz, const double* mu, const std::vector<int> &feature_face_pair, const std::vector<std::pair<int, int>> &feature_neighber_vert_localidx_pairidx);
	bool compute_exp_misp_energy_refine_polycube(const int& vc_size, const double& old_e, double& new_e,
		const double* posx, const double* posy, const double* posz,
		const std::vector<std::vector<double> >& vc_S, double* exp_vec,
		const double& npx, const double& npy, const double& npz,
		double alpha, double beta,
		const double& ga, const int& vf_size,
		const double* vf_px, const double* vf_py, const double* vf_pz, const double* mu,
		const int& vc_size_pair, const double& old_e_pair, double& new_e_pair,
		const double* posx_pair, const double* posy_pair, const double* posz_pair,
		const std::vector<std::vector<double> >& vc_S_pair, double* exp_vec_pair,
		const double& npx_pair, const double& npy_pair, const double& npz_pair,
		const double& ga_pair, const int& vf_size_pair,
		const double* vf_px_pair, const double* vf_py_pair, const double* vf_pz_pair, const double* mu_pair);
    bool compute_exp_misp_energy_refine_polycube_three_cut(
        const std::array<double, 3> &a_vc_size, const double &old_e, std::array<double, 3> &a_new_e,
        const std::array<double *, 3> &a_posx, const std::array<double *, 3> &a_posy, const std::array<double *, 3> &a_posz,
        const std::array<std::vector<std::vector<double>>, 3> &avv_vc_S, std::array<double *, 3> &a_exp_vec,
        const std::array<Eigen::Vector3d, 3> &a_new_coord,
        const double &alpha, const double &beta, const double &ga, const std::array<int, 3> &a_vf_size,
        const std::array<double *, 3> &a_vf_px, const std::array<double *, 3> a_vf_py, const std::array<double *, 3> &a_vf_pz, const std::array<double *, 3> a_mu
    );
		bool local_check_negative_volume4(const int& vc_size, const double* posx, const double* posy, const double* posz,
		const double* vc_n_cross_x, const double* vc_n_cross_y, const double* vc_n_cross_z,
		const double& npx, const double& npy, const double& npz);
    Eigen::Vector3d rotate_vector(const double &theta, const Eigen::Vector3d &axis, const Eigen::Vector3d &x);
    Eigen::Vector3d project2plane(const Eigen::Vector3d &normal, const Eigen::Vector3d &d);
    void compute_vertex_update_info(
        const bool update_all, const bool bv_flag,
        const int &v_id,
        const double *&p_dpx, const double *&p_dpy, const double *&p_dpz,
        const double &alpha, const double &beta,
        std::vector<std::vector<double>> &vc_S,
        double *&p_vf_pos_x, double *&p_vf_pos_y, double *&p_vf_pos_z,
        double *&p_vc_pos_x, double *&p_vc_pos_y, double *&p_vc_pos_z,
        double *&p_vc_n_cross_x, double *&p_vc_n_cross_y, double *&p_vc_n_cross_z,
        double *&exp_vec, double *&mu_vec,
        double *&gx_vec, double *&gy_vec, double *&gz_vec,
        double &local_energy_all, double &normal_e_all,
        double &gx_new, double &gy_new, double &gz_new,
        double &min_radius, int &bvf_size
    );
	void exp_mips_deformation_refine_one_polycube(int v_id, double angle_area_ratio, double energy_power, const int& omp_id, double ga = 1.0, bool update_all = true, bool feature_flag = false);
	void exp_mips_deformation_refine_one_polycube_normal_constrained(int v_id, double angle_area_ratio, double energy_power, const int& omp_id, double ga = 1.0, bool update_all = true);
	void exp_mips_deformation_refine_one_polycube_normal_constrained_single_line(int v_id, double angle_area_ratio, double energy_power, const int& omp_id, const double move_dir[3], double ga = 1.0, bool update_all = true);
	void exp_mips_deformation_refine_one_polycube_normal_constrained_equal_face(int v_id, double angle_area_ratio, double energy_power, const int& omp_id, const  double sym_M[3][3], int common_vert, double ga = 1.0, bool update_all = true);
	void exp_mips_deformation_refine_one_polycube_normal_constrained_three_cut(
        const int &v_id,
        const double &angle_area_ratio,
        const double &energy_power,
        const int &omp_id,
        const std::array<int, 3> &adjacent_one_cut_type,
        const int &common_vert,
        const double &ga,
        const bool update_all
    );
	
	void get_pair_coordinate(double x, double y, double z, const double M[3][3], double center_x, double center_y, double center_z, double &x_pair, double & y_pair, double &z_pair);
	double get_min_edge_length()
	{
		return min_boundary_edge_length;
	}
	double get_ave_edge_length()
	{
		return avg_boundary_edge_length;
	}
	void check_big_distortion_cell(
		const int& vc_size, std::vector<int>& large_distortion_flag,
		const double* posx, const double* posy, const double* posz,
		const std::vector<std::vector<double> >& vc_S,
		const double& npx, const double& npy, const double& npz);
public:
	bool prepare_ok;
	//basic data
	std::vector< double > cell_volume;
	std::vector< std::vector<int> > cell_vertex;
	std::vector< std::vector<int> > cell_vertex_right_order;
	std::vector< Eigen::Matrix3d > cell_S;
	std::vector< Eigen::Matrix3d > cell_S_ori;
	std::vector< std::vector<std::vector<int> > > cell_vertex_vertex;
	std::vector< std::vector<int>> vertex_cell;
	std::vector< std::vector<std::vector<int>>> vertex_cell_vertex;
	std::vector< std::vector<std::vector<double> > > vcv_S;
	std::vector<double> dpx; std::vector<double> dpy; std::vector<double> dpz;
	std::vector<double> dpx_feature; std::vector<double> dpy_feature; std::vector<double> dpz_feature; //same size as dpx
	std::vector<double> dpx_backup; std::vector<double> dpy_backup; std::vector<double> dpz_backup;
	int number_of_color; int max_vc_size;
	std::vector< std::vector<int> > vertex_diff_color;
	std::vector<std::vector<std::vector<double> >> vc_pos_x_omp;
	std::vector<std::vector<std::vector<double> >> vc_pos_y_omp;
	std::vector<std::vector<std::vector<double> >> vc_pos_z_omp;
	std::vector<std::vector<double >> vc_pos_x2_omp;
	std::vector<std::vector<double >> vc_pos_y2_omp;
	std::vector<std::vector<double >> vc_pos_z2_omp;
	std::vector<std::vector<double >> vc_pos_x3_omp;
	std::vector<std::vector<double >> vc_pos_y3_omp;
	std::vector<std::vector<double >> vc_pos_z3_omp;
	std::vector<std::vector<std::vector<double> >> vc_S_omp;
	std::vector<std::vector<double >> vc_n_cross_x_omp;
	std::vector<std::vector<double >> vc_n_cross_y_omp;
	std::vector<std::vector<double >> vc_n_cross_z_omp;
	std::vector<std::vector<double >> exp_vec_omp;
	std::vector<std::vector<double >> gx_vec_omp;
	std::vector<std::vector<double >> gy_vec_omp;
	std::vector<std::vector<double >> gz_vec_omp;
	std::vector<std::vector<double >> mu_vec_omp;
	std::vector<std::vector<int >> large_dis_flag_omp;
	
	std::vector<int> is_polycube_handles;
	std::vector<OpenVolumeMesh::CellHandle> flipped_cell;
	std::vector<OpenVolumeMesh::Geometry::Vec3i> bfv_id;
	std::vector< std::vector<int> > bvf_id; //only contains boundary cells
	std::vector<int> a2b; //all verts to boundary verts
	std::vector<int> boundary_verts;
	std::vector< std::vector<int> > bvv_id;
	std::vector< OpenVolumeMesh::Geometry::Vec3d > target_bfn;
	std::vector<int> error_edges;
	std::vector< OpenVolumeMesh::Geometry::Vec3d > face_center; //size equal to #cells
	std::vector<std::vector<int>> chart_tet_idx;
	std::vector<int> bf_chart;
	std::vector<std::vector<int>> polycube_chart; //each chart which stores the faces
	std::vector<int> polycube_chart_label;//0 for x, 1 for y, 2 for z, for each chart
	std::vector<double> chart_mean_value; //the mean value of each chart
	std::vector<int> up_chart_id; std::vector<int> down_chart_id; std::vector<double> diff_up_down;
	std::vector<std::vector<int> > bef_id;
	std::vector<std::vector<int> > bfe_id;
	std::vector<OpenVolumeMesh::Geometry::Vec3i> vertex_type;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> src_pos;
	double avg_boundary_edge_length;
	double min_boundary_edge_length;
	int boundary_face_number; 
	int boundary_vertface_number; //only initialize in prepare funciton for ovm format
	std::vector<int> is_bv;
	double sigma_s; double sigma_r; int RGNF_iter_count;
	void filter_boundary_normal_RGNF(double sigma_s_, double sigma_r_, int iter_count_, bool xyz_flag = false);
	void filter_boundary_normal_ring(double sigma_s_, int iter_count_, bool xyz_flag = false);
	void filter_boundary_normal_feature(double sigma_s_, double sigma_r_, int level);
	void filter_feature_edge_ring(int n_level, double ss, bool consider_face_normal = true);
	void get_polyline_normal();
	void set_feature_edge_array();
	void set_feature_edge_length();
	bool target_normal_check_repair(bool repair_flag = false);
	bool check_normal_difference(double th = 45);
	void compute_face_normal(std::vector< OpenVolumeMesh::Geometry::Vec3d >& bfn);
	double compute_energy_ratio_filter_deform();
	
	std::vector<double> local_energy_ratio;
	std::vector<int> change_big_flag;
	std::vector< OpenVolumeMesh::Geometry::Vec3d > last_bfn;
	std::vector<int> distortion_big_count;
	std::vector<int> distortion_big_cell;
	std::vector<int> distortion_big_count_update;
	std::vector<int> distortion_big_cell_update;
	std::vector<std::vector<int>> distortion_big_count_omp;
	std::vector<std::vector<int>> vc_count_omp;
	double bound_K;
	//for seamless deformation
	std::vector<int> fixed_chart;
	std::vector<int> fixed_label;
	std::vector<int> cut_face_tet_idx;
	std::vector<std::vector<int>> equal_triangles;
	std::vector<int> vert_update_type; // 0 : no cut, 1: update part, 2: no update part, 3: three_cut;
	//std::vector<int> vert_cut_type; // 0: no cut, for the update part of the cut, set vert_cut_type to be vert type.
	std::vector<int> vert_belong_cut_array;
	std::vector<int> cut_types; // totally 12 types
	std::vector<int> cut_common_verts_idx;
	std::map<int, int> vert_pairs_map;
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
	const double type_four_matrix[12][3] = {
		{0,0,1},
		{0,0,1},
		{0,0,1},
		{0,0,1},
		{0,1,0},
		{0,1,0},
		{0,1,0},
		{0,1,0},
		{1,0,0},
		{1,0,0},
		{1,0,0},
		{1,0,0},
	};
	//three cut_info
    static const std::array<Eigen::Vector3d, 8> a_three_cut_opposite_dir;
	std::vector<int> vert_belong_three_cut_array; //the idx of a vert in std::vector<ThreeCut>
	std::vector<int> three_cut_types;
    std::vector<std::array<int, 3>> three_cut_adjacent_one_cut_type;
	std::vector<int> three_cut_common_vertex_idx;
	std::map<int, std::pair<int, int>> three_cut_vert_to_pair;  //the two neightbors of a vertex can be find by this map
	//feature edge info
	bool feature_edge_change_flag;
	std::vector<bool> feature_edge_flag;
	std::vector<bool> corner_vert_flag;
	std::vector<std::vector<int>> feature_v2e; //global coordinats
	std::vector<std::vector<int>> feature_v2v; 
	std::vector<std::vector<std::pair<int, int>>> feature_neighbor_vert2cellpair; //first face containging vert, second face the other direction
	std::vector<OpenVolumeMesh::Geometry::Vec3d> target_fea_n; //size equal to ne
	std::vector<OpenVolumeMesh::Geometry::Vec3d> polyline_target_n;
	std::vector<int> feature_edge_array;
	std::vector<double> feature_edge_length;	//same size as feature_edge_array
	std::vector<std::pair<int, int>> feature_e2v; //consistent feature edges
	std::vector<std::vector<int>> feature_neighbor; //neighbor of edges at most 2 neighbors, cut by valence 3 vertex
	std::vector<std::vector<int>> long_feature_edge;
	std::vector<int> feature_edge2long;
	std::vector<int> feature_g2l; //local to global
	std::vector<int> feature_l2g; //global to local
	std::vector<double> lambda_array;
	std::vector<double> lambda_array_polyline;
	
};
#endif // !TETRAHEDRAL_DEFORMATION_H