#ifndef POLYCUBE_BOUNDARY_MAPPING_IVF_H
#define POLYCUBE_BOUNDARY_MAPPING_IVF_H
#include "MeshDefinition.h"
#include "TetStructure.h"
#include "SmallVec.h"
#include <Eigen\Dense>
#include <Eigen/Sparse>

using ig::Tetrahedron;
using ig::TetStructure;
using ig::TetVertex;
using ig::CVec;

class polycube_boundary_mapping_IVF
{
public:
	polycube_boundary_mapping_IVF();
	~polycube_boundary_mapping_IVF();
	void reset();

	void optimize_IVF(VolumeMesh* mesh_, int chart_label,
		const std::vector<OpenVolumeMesh::Geometry::Vec3d>& src_pos,
		const std::vector<OpenVolumeMesh::Geometry::Vec3i>& bfv_id, 
		const std::vector< std::vector<int> >& bvf_id,
		const std::vector<int>& map_v, const int& var_v_count, const std::vector<int>& map_c, const int& var_c_count);

	void optimize_IVF(TetStructure<double> *tet_mesh_, int chart_label,
		const std::vector<OpenVolumeMesh::Geometry::Vec3d>& src_pos,
		const std::vector<OpenVolumeMesh::Geometry::Vec3i>& bfv_id,
		const std::vector< std::vector<int> >& bvf_id,
		const std::vector<int>& map_v, const int& var_v_count, const std::vector<int>& map_c, const int& var_c_count);
	//void optimize_ivf_block(Mesh* mesh_);

private:
	std::vector<OpenVolumeMesh::Geometry::Vec3d> fe1;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> fe2;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> bp;

	Eigen::SparseMatrix<double> CTC;
	Eigen::SparseMatrix<double> CT;
	Eigen::Matrix<double, Eigen::Dynamic, 1> CB;
	double CBTCB;

	Eigen::SparseMatrix<double> hessian;
	Eigen::Matrix<double, Eigen::Dynamic, 1> gradient;
	Eigen::Matrix<double, Eigen::Dynamic, 1> prevSolution;
	Eigen::Matrix<double, Eigen::Dynamic, 1> solution;
	Eigen::Matrix<double, Eigen::Dynamic, 1> step;
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Upper> lltSolver;
	double lambda; double min_lambda; double max_lambda;
	double ivf_alpha; double min_ivf_alpha; double max_ivf_alpha;
	double step_size; double min_step_size; double max_step_size;
	double current_ivf_energy;
	double current_dis_energy; double pre_current_dis_energy;
	double amips_s; int var_num;

	void update_ivf_alpha();
	bool initilize_affine_trans(VolumeMesh* mesh_, int chart_label,
		const std::vector<OpenVolumeMesh::Geometry::Vec3d>& src_pos,
		const std::vector<OpenVolumeMesh::Geometry::Vec3i>& bfv_id,
		const std::vector< std::vector<int> >& bvf_id,
		const std::vector<int>& map_v, const int& var_v_count, const std::vector<int>& map_c, const int& var_c_count);
	
	bool initilize_affine_trans(TetStructure<double>* tet_mesh_, int chart_label,
		const std::vector<OpenVolumeMesh::Geometry::Vec3d>& src_pos,
		const std::vector<OpenVolumeMesh::Geometry::Vec3i>& bfv_id,
		const std::vector< std::vector<int> >& bvf_id,
		const std::vector<int>& map_v, const int& var_v_count, const std::vector<int>& map_c, const int& var_c_count);



	void optimize_AT(const std::vector<int>& map_c); //affine transformation
	double compute_energy_derivative_AT(const std::vector<int>& map_c, const Eigen::Matrix<double, Eigen::Dynamic, 1>& x);
	double compute_only_energy_AT(const std::vector<int>& map_c, const Eigen::Matrix<double, Eigen::Dynamic, 1>& x);
	void construct_face_vertex_pos(VolumeMesh* mesh_, int chart_label,
		const std::vector<OpenVolumeMesh::Geometry::Vec3d>& src_pos,
		const std::vector< std::vector<int> >& bvf_id,
		const std::vector<int>& map_c, const std::vector<int>& map_v);

	void construct_face_vertex_pos(TetStructure<double> *tet_mesh_, int chart_label,
		const std::vector<OpenVolumeMesh::Geometry::Vec3d>& src_pos,
		const std::vector< std::vector<int> >& bvf_id,
		const std::vector<int>& map_c, const std::vector<int>& map_v);
};

#endif // !POLYCUBE_BOUNDARY_MAPPING_IVF_H

