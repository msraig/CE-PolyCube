#pragma warning( disable : 4267)
#ifdef _DEBUG
#pragma comment(lib, "SUSd")
#pragma comment(lib, "SUS_standardd")
#pragma comment(lib, "SUS_polycubed")
#else
#pragma comment(lib, "SUS")
#pragma comment(lib, "SUS_standard")
#pragma comment(lib, "SUS_polycube")
#endif
#include <fstream>
#include <string>
#include <array>
#include <algorithm>  
#include <vector>     
#include <ctime>      
#include <cstdlib>    
#include "sus_lib.h"
#include "sus_lib_standard.h"
#include "sus_lib_polycube.h"
#include "Helper.h"
#include "Polycube_Algorithm.h"
#include "Sparse_Matrix.h"
#include "Sparse_Solver.h"
#include "PolyCubeCutInterface.h"
#include "Hex.h"
#include "TetOperation.h"
#define MIN_CUT_LENGTH 0.001
#define INIT_SIGMA_R 0.5
#define N_SPLIT 90.0 // used for determine min_polycube_edge_length
#define FEATURE_FLAG 1
#define CUT_EDGE_MAP_EPSILON 1e-6
#define DEFAULT_TET_QUALITY 2.0
using ig::TetOperation;
Polycube_Algorithm::Polycube_Algorithm()
{
	mesh_ = NULL;
	tet_mesh_ = NULL;
	tet_mesh_polycube_ = NULL;
	tet_mesh_unchanged = NULL;
	tet_mesh_polycube_unchanged = NULL;
	avg_boundary_edge_length = -1.0;
	min_polycube_edge_length = -1.0;
	sigma_r = INIT_SIGMA_R;
	final_cube_length = -1.0;
	//sigma_r = 1;
	init_distortion[0] = -1.0;
	init_distortion[1] = -1.0;
	init_distortion[2] = -1.0;
	init_volume_distortion[0] = -1.0;
	init_volume_distortion[1] = -1.0;
	init_volume_distortion[2] = -1.0;
	v_distortion = NULL;
	feature_edge_flag = NULL;
	reset_aLL_state_vtk();
	ovm_color_status = false;
	data_prepare_status = false;
	polycube_exist_flag = false;
	feature_status = FEATURE_FLAG;
}
Polycube_Algorithm::Polycube_Algorithm(VolumeMesh* mesh)
	: mesh_(mesh)
{
	ovm_color_status = false;
	data_prepare_status = false;
	polycube_exist_flag = false;
	feature_status = FEATURE_FLAG;
	reset_aLL_state_vtk();
	data_preparation_ovm();
}
Polycube_Algorithm::~Polycube_Algorithm()
{
	if (tet_mesh_unchanged) delete tet_mesh_unchanged;
	if (tet_mesh_polycube_unchanged) delete tet_mesh_polycube_unchanged;
	if (pq_deform) delete pq_deform;
}
void Polycube_Algorithm::reset_aLL_state()
{
	R.clear(); PT.clear();
	cell_volume.clear(); fix_cell_boundary.clear();
	cv_id.clear(); bfv_id.clear(); bf_n.clear();
	initlize_pc(); 
}
void Polycube_Algorithm::reset_aLL_state_vtk()
{
	R.clear(); PT.clear();
	cell_volume.clear(); fix_cell_boundary.clear();
	cv_id.clear(); bfv_id.clear(); bf_n.clear();
	fix_cell_boundary.clear();
	src_pos.clear();
	fix_charts.clear();
	fix_labels.clear();
	equal_faces_array.clear();
	initlize_pc_vtk();
	feature_edge_flag = NULL;
}
void Polycube_Algorithm::get_boundary_avg_edge_length()
{
	//if calculated in the previous step, then return directly
	if (avg_boundary_edge_length > 0.0 || bfv_id.size() == 0)
		return;
	avg_boundary_edge_length = 0.0; OpenVolumeMesh::Geometry::Vec3d p0, p1, p2; double count = 0.0;
	for (int i = 0; i < bfv_id.size(); ++i)
	{
		OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
		if (one_bfv_id[0] < 0) continue;
		p0 = src_pos[one_bfv_id[0]];
		p1 = src_pos[one_bfv_id[1]];
		p2 = src_pos[one_bfv_id[2]];
		avg_boundary_edge_length += (p0 - p1).norm();
		avg_boundary_edge_length += (p1 - p2).norm();
		avg_boundary_edge_length += (p2 - p0).norm();
		count += 3.0;
	}
	avg_boundary_edge_length /= count;
}
void Polycube_Algorithm::get_bounding_box()
{
	//get bb for tet_mesh_
	if (tet_mesh_ == NULL)
	{
		return;
	}
	std::vector<TetVertex<double>* > &tetra_vertices = tet_mesh_->tetra_vertices;
	bb_Min = tetra_vertices[0]->pos;
	bb_Max = tetra_vertices[0]->pos;
	int n_vert = (int)tetra_vertices.size();
	for (size_t i = 1; i < n_vert; i++)
	{
		ig::CVec<double, 3> tmp_pos = tetra_vertices[i]->pos;
		for (size_t j = 0; j < 3; j++)
		{
			if (bb_Min[j] > tmp_pos[j])
				bb_Min[j] = tmp_pos[j];
			if (bb_Max[j] < tmp_pos[j])
				bb_Max[j] = tmp_pos[j];
		}
	}
}
void Polycube_Algorithm::get_min_polycube_edge_length()
{
	//min polycube edge length: min(pq edge length, avg_boundary_edge_length, bb/SPLIT)
	if (min_polycube_edge_length > 0.0 || bfv_id.size() == 0)
		return;
	if (polycube_edge_length.size() == 0)
		return;
	min_polycube_edge_length = polycube_edge_length[0];
	for (size_t i = 1; i < polycube_edge_length.size(); i++)
	{
		if (min_polycube_edge_length > polycube_edge_length[i])
			min_polycube_edge_length = polycube_edge_length[i];
	}
	if (avg_boundary_edge_length > 0.0)
	{
		if (min_polycube_edge_length > avg_boundary_edge_length)
			min_polycube_edge_length = avg_boundary_edge_length;
	}
	//finally set to avg_boundary_edge_length
	min_polycube_edge_length = avg_boundary_edge_length;
	get_bounding_box();
	double bb[3];
	for (size_t i = 0; i < 3; i++)
	{
		bb[i] = (bb_Max[i] - bb_Min[i]) / N_SPLIT;
	}
	double max_bb = std::max(std::max(bb[0], bb[1]), bb[2]);
	if (min_polycube_edge_length > max_bb)
		min_polycube_edge_length = max_bb;
}
void Polycube_Algorithm::initlize_pc()
{
	//initialize the affine matrix
	if (!mesh_) return;
	int nc = mesh_->n_cells();
	R.resize(nc); PT.resize(nc); cell_volume.resize(nc); cv_id.resize(nc);
	OpenVolumeMesh::Geometry::Vec4i one_cv_id; Eigen::Matrix3d p;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> cv_p(4);
	for (OpenVolumeMesh::CellIter c_it = mesh_->cells_begin(); c_it != mesh_->cells_end(); ++c_it)
	{
		int c_id = c_it->idx(); int cv_count = 0;
		for (OpenVolumeMesh::CellVertexIter cv_it = mesh_->cv_iter(c_it.cur_handle()); cv_it; ++cv_it)
		{
			one_cv_id[cv_count] = cv_it->idx();
			cv_p[cv_count] = mesh_->vertex(*cv_it);
			++cv_count;
		}
		cv_id[c_id] = one_cv_id;
		p(0, 0) = cv_p[1][0] - cv_p[0][0]; p(1, 0) = cv_p[1][1] - cv_p[0][1]; p(2, 0) = cv_p[1][2] - cv_p[0][2];
		p(0, 1) = cv_p[2][0] - cv_p[0][0]; p(1, 1) = cv_p[2][1] - cv_p[0][1]; p(2, 1) = cv_p[2][2] - cv_p[0][2];
		p(0, 2) = cv_p[3][0] - cv_p[0][0]; p(1, 2) = cv_p[3][1] - cv_p[0][1]; p(2, 2) = cv_p[3][2] - cv_p[0][2];
		cell_volume[c_id] = std::abs(p.determinant());
		PT[c_id] = p.inverse();
	}
	bfv_id.clear(); bfv_id.resize(nc, OpenVolumeMesh::Geometry::Vec3i(-1, -1, -1));
	for (OpenVolumeMesh::FaceIter f_it = mesh_->faces_begin(); f_it != mesh_->faces_end(); ++f_it)
	{
		OpenVolumeMesh::HalfFaceHandle hfh0 = mesh_->halfface_handle(*f_it, 0);
		OpenVolumeMesh::HalfFaceHandle hfh1 = mesh_->halfface_handle(*f_it, 1);
		OpenVolumeMesh::CellHandle ch0 = mesh_->incident_cell(hfh0);
		OpenVolumeMesh::CellHandle ch1 = mesh_->incident_cell(hfh1);
		if (ch0 == VolumeMesh::InvalidCellHandle && ch1 != VolumeMesh::InvalidCellHandle)
		{
			int hfv_count = 0; int c_id = ch1.idx();
			for (OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh0); hfv_it; ++hfv_it)
			{
				int hfv_id = hfv_it->idx();
				bfv_id[c_id][hfv_count] = hfv_id;
				++hfv_count;
			}
		}
		else if (ch1 == VolumeMesh::InvalidCellHandle && ch0 != VolumeMesh::InvalidCellHandle)
		{
			int hfv_count = 0; int c_id = ch0.idx();
			for (OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh1); hfv_it; ++hfv_it)
			{
				int hfv_id = hfv_it->idx();
				bfv_id[c_id][hfv_count] = hfv_id;
				++hfv_count;
			}
		}
		else if (ch1 == VolumeMesh::InvalidCellHandle && ch0 == VolumeMesh::InvalidCellHandle)
		{
			printf("Error : Both two halffaces have no cells!\n");
		}
	}
	src_pos.resize(mesh_->n_vertices());
	for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		src_pos[v_it->idx()] = mesh_->vertex(*v_it);
	}
}
void Polycube_Algorithm::initlize_pc_vtk()
{
	if (!tet_mesh_) return;
	int nc = (int)tet_mesh_->tetras.size();
	std::vector<Tetrahedron<double>*> &tetras = tet_mesh_->tetras;
	R.resize(nc); PT.resize(nc); cell_volume.resize(nc); cv_id.resize(nc);
	OpenVolumeMesh::Geometry::Vec4i one_cv_id; Eigen::Matrix3d p;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> cv_p(4);
	avg_cell_volume = 0.0;
	int neg_cells = 0;
	for (size_t i = 0; i < nc; i++)
	{
		int c_id = (int)i; int cv_count = 0;
		for (; cv_count < 4; cv_count++)
		{
			one_cv_id[cv_count] = tet_mesh_->tetras[i]->vertex[cv_count]->id;
			cv_p[cv_count][0] = tetras[i]->vertex[cv_count]->pos[0];
			cv_p[cv_count][1] = tetras[i]->vertex[cv_count]->pos[1];
			cv_p[cv_count][2] = tetras[i]->vertex[cv_count]->pos[2];
		}
		cv_id[c_id] = one_cv_id;
		p(0, 0) = cv_p[1][0] - cv_p[0][0]; p(1, 0) = cv_p[1][1] - cv_p[0][1]; p(2, 0) = cv_p[1][2] - cv_p[0][2];
		p(0, 1) = cv_p[2][0] - cv_p[0][0]; p(1, 1) = cv_p[2][1] - cv_p[0][1]; p(2, 1) = cv_p[2][2] - cv_p[0][2];
		p(0, 2) = cv_p[3][0] - cv_p[0][0]; p(1, 2) = cv_p[3][1] - cv_p[0][1]; p(2, 2) = cv_p[3][2] - cv_p[0][2];
		cell_volume[c_id] = (p.determinant());
		if (cell_volume[c_id] < 0)
		{
			neg_cells++;
		}
		avg_cell_volume += cell_volume[c_id];
		PT[c_id] = p.inverse();
	}
	avg_cell_volume = avg_cell_volume / nc;
	std::cout << "Avg_cell_volume: " << avg_cell_volume << std::endl;
	std::cout << "Negative cell number: " << neg_cells << std::endl;
	bfv_id.clear(); bfv_id.resize(nc, OpenVolumeMesh::Geometry::Vec3i(-1, -1, -1));
	//bfv stores the information of boundary face vertices
	static int s_tet_id[4][3] = { { 1, 2, 3}, {2, 0, 3}, {0, 1, 3}, {1, 0, 2} };
	for (size_t i = 0; i < nc; i++)
	{
		int c_id = (int)i;
		if (!tetras[i]->is_onboundary())
		{
			continue;
		}
		int face_id = -1;
		for (size_t j = 0; j < 4; j++)
		{
			TetVertex<double>* tv_temp = tetras[i]->vertex[j];
			if (!(tv_temp->boundary))
			{
				face_id = (int)j;
				break;
			}
		}
		if (face_id == -1)
		{
			std::cout << "four verts on boundary" << std::endl;
			return;
		}
		for (size_t hfv_count = 0; hfv_count < 3; hfv_count++)
		{
			int hfv_id = (int)tetras[i]->vertex[s_tet_id[face_id][hfv_count]]->id;
			bfv_id[c_id][hfv_count] = hfv_id;
		}
	}
	int n_verts = (int)tet_mesh_->tetra_vertices.size();
	src_pos.clear();
	src_pos.resize(n_verts);
	std::vector<TetVertex<double>* > &tetra_vertices = tet_mesh_->tetra_vertices;
	for (size_t i = 0; i < n_verts; i++)
	{
		src_pos[i][0] = tetra_vertices[i]->pos[0];
		src_pos[i][1] = tetra_vertices[i]->pos[1];
		src_pos[i][2] = tetra_vertices[i]->pos[2];
	}
	get_boundary_avg_edge_length();
}
bool Polycube_Algorithm::construct_original_tet()
{
	//construct tet_mesh based on tet_mesh_polycube based on polycube coordinate
	std::vector<CVec<double, 3>> points_polycube, points_ori, points_changed;
	std::vector<unsigned int> indices_polycube, indices_ori;
	get_vert_cell(tet_mesh_polycube_, points_polycube, indices_polycube);
	get_vert_cell(tet_mesh_unchanged, points_ori, indices_ori);
	int nv = (int)tet_mesh_polycube_->tetra_vertices.size();
	const int s_tet_id[4][3] = { { 2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2} };
	for (size_t i = 0; i < nv; i++)
	{
		Tetrahedron<double> *mt = tet_mesh_polycube_unchanged->locate_point(points_polycube[i]);
		if (mt == NULL)
		{
			std::cout << "NULL pointer when locating polycube!" << std::endl;
			return false;
		}
		CVec<double, 3> p[4];
		p[0] = mt->vertex[0]->pos;
		p[1] = mt->vertex[1]->pos;
		p[2] = mt->vertex[2]->pos;
		p[3] = mt->vertex[3]->pos;
		int idx[4];
		idx[0] = (int)mt->vertex[0]->id;
		idx[1] = (int)mt->vertex[1]->id;
		idx[2] = (int)mt->vertex[2]->id;
		idx[3] = (int)mt->vertex[3]->id;
		double v_all;
		double v[4];
		double w[4];
		v_all = (ig::CGALHelper::compute_cell_volume(p[0], p[1], p[2], p[3]));
		for (size_t j = 0; j < 4; j++)
		{
			v[j] = (ig::CGALHelper::compute_cell_volume(p[s_tet_id[j][0]], p[s_tet_id[j][1]], p[s_tet_id[j][2]], points_polycube[i]));
		}
		if (std::abs(v_all) < 1.0e-14)
		{
			for (size_t j = 0; j < 4; j++)
			{
				w[j] = 0.25;
			}
		}
		else
		{
			for (size_t j = 0; j < 4; j++)
			{
				w[j] = v[j] / v_all;
			}
		}
		CVec<double, 3> temp_point(0, 0, 0);
		for (size_t j = 0; j < 4; j++)
		{
			temp_point = temp_point + points_ori[idx[j]] * w[j];
		}
		points_changed.push_back(temp_point);
	}
	tet_mesh_->load_tet(points_changed, indices_polycube);
	return true;
}
void Polycube_Algorithm::get_corresponding_point(TetStructure<double>* source, TetStructure<double>* target, const std::vector<ig::CVec<double, 3>>& source_point, std::vector<ig::CVec<double, 3>> &target_point)
{
	const static int s_tet_id[4][3] = { { 2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2} };
	std::vector<CVec<double, 3>> points_ori, points_changed;
	std::vector<unsigned int> indices_ori;
	get_vert_cell(target, points_ori, indices_ori);
	target_point.clear();
	int np = (int)source_point.size();
	for (size_t i = 0; i < np; i++)
	{
		Tetrahedron<double> *mt = source->locate_point(source_point[i]);
		if (mt == NULL)
		{
			std::cout << "NULL pointer when locating polycube!" << std::endl;
			return;
		}
		CVec<double, 3> p[4];
		p[0] = mt->vertex[0]->pos;
		p[1] = mt->vertex[1]->pos;
		p[2] = mt->vertex[2]->pos;
		p[3] = mt->vertex[3]->pos;
		int idx[4];
		idx[0] = (int)mt->vertex[0]->id;
		idx[1] = (int)mt->vertex[1]->id;
		idx[2] = (int)mt->vertex[2]->id;
		idx[3] = (int)mt->vertex[3]->id;
		double v_all;
		double v[4];
		double w[4];
		v_all = (ig::CGALHelper::compute_cell_volume(p[0], p[1], p[2], p[3]));
		for (size_t j = 0; j < 4; j++)
		{
			v[j] = (ig::CGALHelper::compute_cell_volume(p[s_tet_id[j][0]], p[s_tet_id[j][1]], p[s_tet_id[j][2]], source_point[i]));
			w[j] = v[j] / v_all;
		}
		if (v_all < 0.0000001 * EPSILON)
		{
			for (size_t j = 0; j < 4; j++)
			{
				w[j] = 0.25;
			}
		}
		CVec<double, 3> temp_point(0, 0, 0);
		for (size_t j = 0; j < 4; j++)
		{
			temp_point = temp_point + points_ori[idx[j]] * w[j];
		}
		target_point.push_back(temp_point);
	}
}
void Polycube_Algorithm::SetSigmaR(double sigma_r_)
{
	sigma_r = sigma_r_;
}
void Polycube_Algorithm::SetMeshPolycube(TetStructure<double>* tet_mesh)
{
	tet_mesh_polycube_ = tet_mesh;
	//tet_mesh and tet_mesh_polycube_should share the same structure
	assert(tet_mesh_->tetra_vertices.size() == tet_mesh_polycube_->tetra_vertices.size() && tet_mesh_->tetras.size() == tet_mesh_polycube_->tetras.size());
	//first index should be the same
	assert(tet_mesh_->tetras[0]->vertex[0]->id == tet_mesh_polycube_->tetras[0]->vertex[0]->id);
	tet_mesh_polycube_unchanged = copy_mesh(tet_mesh_polycube_);
	tet_mesh_polycube_unchanged->build_kdtree();
	(*pq_deform).prepare_for_deformation(tet_mesh_polycube_);
	compute_flips();
	polycube_exist_flag = true;
}
int Polycube_Algorithm::auto_cut_batch(double cut_depth_ratio, double max_volume, int batch_size, double cutting_edge_ratio, int distortion_type, int ivf_iter, double min_depth_ratio, bool flags_cut_convex_edge, bool flag_sample_hex, int default_cutnum, double tet_threshold, bool random_cut)
{
	//return the number of cut
	//number of cutting edges: # of available edges * cutting_edge_ratio
	if (!tet_mesh_polycube_)
		return -1;
	if (!do_Extract_polycube_vtk(distortion_type, flag_sample_hex)) return -2;
	double depth_decay_ratio = 0.5;
	int truncated_idx = 0;
	for (size_t i = 0; i < sorted_polycube_edge_distortion.size(); i++)
	{
		if (sorted_polycube_edge_distortion[i] < tet_threshold)
		{
			truncated_idx = (int)i;
			break;
		}
	}
	if (random_cut)
	{
		randomize_array(sorted_polycube_edge);
	}
	else
	{
		if (truncated_idx != 0)
		{
			sorted_polycube_edge.resize(truncated_idx);
		}
	}
	std::vector<ig::CVec<double, 3>> v_input_point = cut_input_point;
	std::vector<unsigned int> v_input_face = cut_input_face;
	std::vector<int> v_input_chart = cut_input_chart;
	std::vector<int> v_input_label = cut_input_label;
	PolyCubeCutInterface pci(v_input_point, v_input_face, v_input_chart, v_input_label);
	std::vector<unsigned int> v_cuttable_edge;
	std::vector<unsigned int> avail_edges;
	std::set<std::pair<unsigned, unsigned>> v_cuttable_edge_set;
	pci.get_cuttable_edge(v_cuttable_edge);
	std::map<std::pair<unsigned int, unsigned int>, bool> map_is_convex_edge;
	pci.get_convex_map(map_is_convex_edge);
	//update one_cut_depth
	double sqrt_two = sqrt(2.0);
	for (size_t i = 0; i < polycube_edge.size(); i++)
	{
		if (polycube_edge_max_cut_depth_one_cut[i] > polycube_edge_neighbor_min_length[i])
			polycube_edge_max_cut_depth_one_cut[i] = polycube_edge_neighbor_min_length[i];
	}
	for (size_t i = 0; i < v_cuttable_edge.size() / 2; i++)
	{
		unsigned id0 = v_cuttable_edge[2 * i];
		unsigned id1 = v_cuttable_edge[2 * i + 1];
		if (id0 < id1)
		{
			v_cuttable_edge_set.insert(std::pair<unsigned, unsigned>(id0, id1));
		}
		else
		{
			v_cuttable_edge_set.insert(std::pair<unsigned, unsigned>(id1, id0));
		}
	}
	std::vector<double> v_cut_depth;
	std::vector<double> v_max_cut_depth;
	std::map<std::pair<unsigned, unsigned>, std::vector<std::pair<unsigned, unsigned>>> map_associated_edge, map_affected_edge; //associated[a]: edges while is affected by cutting a, affected[a]: edges which wil affect a
	pci.get_associated_edge(map_associated_edge);
	//initialize affected
	for (size_t i = 0; i < polycube_edge.size(); i++)
	{
		std::vector<std::pair<unsigned, unsigned>> empty;
		unsigned min_id, max_id;
		min_id = polycube_edge[i].first;
		max_id = polycube_edge[i].second;
		std::pair<unsigned, unsigned> edge_pair(min_id, max_id);
		map_affected_edge[edge_pair] = empty;
	}
	assert(polycube_edge.size() == map_associated_edge.size());
	//construct affected map
	for (size_t i = 0; i < polycube_edge.size(); i++)
	{
		unsigned min_id, max_id;
		min_id = polycube_edge[i].first;
		max_id = polycube_edge[i].second;
		std::pair<unsigned, unsigned> edge_pair(min_id, max_id);
		std::vector<std::pair<unsigned, unsigned>> &assoc = map_associated_edge[edge_pair];
		for (size_t j = 0; j < assoc.size(); j++)
		{
			unsigned assoc0 = assoc[j].first;
			unsigned assoc1 = assoc[j].second;
			std::pair<unsigned, unsigned> assoc_pair(assoc0, assoc1);
			std::vector<std::pair<unsigned, unsigned>> &affected = map_affected_edge[assoc_pair];
			auto it = std::find(affected.begin(), affected.end(), edge_pair);
			if (it == affected.end())
			{
				//not found
				affected.push_back(edge_pair);
			}
		}
	}
	//min cut depth should be at least half the size of hex cube size, this bound is secure
	//this variable is set both in function auto_cut and select_cut
	double min_cut_depth = min_depth_ratio * sigma_r * min_polycube_edge_length;
	int max_edge_num = (int)sorted_polycube_edge.size();
	//get edges which can be cut deeply enough
	std::vector<unsigned int> v_cut_edge_all, v_cut_edge_all_ori, v_cut_edge;
	std::vector<double> v_cut_depth_all, v_cut_depth_all_ori;
	//v_cut_edge.clear();
	for (size_t i = 0; i < sorted_polycube_edge.size(); i++)
	{
		//only consider edge which is already long enough
		unsigned int id0 = sorted_polycube_edge[i].first;
		unsigned int id1 = sorted_polycube_edge[i].second;
		v_cut_edge_all_ori.push_back(sorted_polycube_edge[i].first);
		v_cut_edge_all_ori.push_back(sorted_polycube_edge[i].second);
	}
	get_cut_depth_min(v_cut_edge_all_ori, v_cut_depth_all_ori, cut_depth_ratio);
	std::cout << "min cut depth: " << min_cut_depth << std::endl;
	std::set<std::pair<unsigned, unsigned>> ignored_edges; //edges with min cut depth smalller than min_cut_detph and the ones who affected such edges
	//construct ignored_edges
	for (size_t i = 0; i < sorted_polycube_edge.size(); i++)
	{
		if (v_cut_depth_all_ori[i] < min_cut_depth)
		{
			//insert this edge and the affected to ignore_edges
			auto search = ignored_edges.find(sorted_polycube_edge[i]);
			if (search == ignored_edges.end())
			{
				//not found
				ignored_edges.insert(sorted_polycube_edge[i]);
				std::vector<std::pair<unsigned, unsigned>> &affected = map_affected_edge[sorted_polycube_edge[i]];
				for (size_t j = 0; j < affected.size(); j++)
				{
					ignored_edges.insert(affected[j]);
				}
			}
		}
	}
	int avail_edge_count = 0;
	//remove convex and concave edge from ignored_edges
	std::vector<std::pair<unsigned, unsigned>> convex_concave;
	for (auto it = ignored_edges.begin(); it != ignored_edges.end(); ++it)
	{
		auto cuttable_it = v_cuttable_edge_set.find(*it);
		if (cuttable_it != v_cuttable_edge_set.end())
		{
			//convex or concave
			//1 cut depth must be large enough
			unsigned tmp_idx = polycube_edge_to_idx_map[(*it)];
			double tmp_depth = polycube_edge_max_cut_depth_one_cut[tmp_idx] * cut_depth_ratio;
			if (tmp_depth > min_cut_depth)
			{
				convex_concave.push_back((*it));
			}
		}
	}
	for (size_t i = 0; i < convex_concave.size(); i++)
	{
		ignored_edges.erase(convex_concave[i]);
	}
	int n_ignored_edges_first_step = (int)ignored_edges.size();
	//ignore edges modified by feature edges here
	if (!init_feature_polycube_edge_array_sorted.empty())
	{
		for (size_t i = 0; i < init_feature_polycube_edge_array_sorted.size(); i++)
		{
			ignored_edges.insert(init_feature_polycube_edge_array_sorted[i]);
		}
	}
	get_available_cut(&pci, min_cut_depth, max_volume, sorted_polycube_edge, ignored_edges, avail_edges);
	if (n_ignored_edges_first_step == 0 && avail_edges.size() == 0)
		return -3;
	v_cut_depth.clear();
	v_cut_edge.clear();
	v_cut_depth_all.clear();
	v_cut_edge_all.clear();
	int final_cut_candidate_num = (int)(cutting_edge_ratio * (avail_edges.size() / 2));
	if (default_cutnum > 0 && default_cutnum <= avail_edges.size() / 2)
	{
		final_cut_candidate_num = default_cutnum;
	}
	std::vector<unsigned> v_cut_edge_success;
	std::vector<double> v_cut_depth_success;
	for (size_t i = 0; i < final_cut_candidate_num; i++)
	{
		v_cut_edge_all.push_back(avail_edges[2 * i]);
		v_cut_edge_all.push_back(avail_edges[2 * i + 1]);
	}
	get_cut_depth(v_cut_edge_all, v_cut_depth_all, cut_depth_ratio, map_associated_edge, (int)cut_input_point.size());
	if (batch_size <= 0 || batch_size > final_cut_candidate_num)
	{
		batch_size = final_cut_candidate_num;
	}
	//int n_batch = (v_cut_depth.size() / batch_size) + 1;
	int n_batch = (int)(std::ceil((double)final_cut_candidate_num / (double)batch_size));
	if (final_cut_candidate_num > 0)
	{
		for (size_t i = 0; i < n_batch; i++)
		{
			//insert elem
			int n_add = batch_size;
			if (batch_size * (i + 1) > final_cut_candidate_num)
			{
				n_add = (int)(final_cut_candidate_num - batch_size * i);
			}
			assert(n_add <= batch_size);
			for (size_t j = 0; j < n_add; j++)
			{
				int tmp_idx = (int)(i * batch_size + j);
				v_cut_edge.push_back(v_cut_edge_all[2 * tmp_idx]);
				v_cut_edge.push_back(v_cut_edge_all[2 * tmp_idx + 1]);
				v_cut_depth.push_back(v_cut_depth_all[tmp_idx]);
			}
			if (pci.is_equivalent_cut(v_cut_edge, v_cut_depth, v_cut_edge_success, v_cut_depth_success))
			{
				//equavalent cut
				continue;
			}
			bool cut_success = true;
			//cut
			ig::CVec<double, 3> cur_distortion(-1, -1, -1);
			while (cur_distortion[0] < 0.0 || cur_distortion[distortion_type] > init_distortion[distortion_type])
			{
				cur_distortion = cut_and_deform(&pci, v_cut_edge, v_cut_depth, max_volume, ivf_iter);
				int change_flag = -1;
				//not success
				if (cur_distortion[0] < 0.0 || cur_distortion[distortion_type] > init_distortion[distortion_type])
				{
					change_flag = 0;
					for (size_t k = 0; k < v_cut_depth.size(); k++)
					{
						v_cut_depth[k] = v_cut_depth[k] * depth_decay_ratio;
						if (v_cut_depth[k] <= min_cut_depth)
						{
							v_cut_depth[k] = min_cut_depth;
						}
						else
						{
							change_flag++;
						}
					}
				}
				if (change_flag == 0) // cut illegal but nothing changes
				{
					cut_success = false;
					break;
				}
			}
			if (cut_success)
			{
				v_cut_edge_success = v_cut_edge;
				v_cut_depth_success = v_cut_depth;
			}
			else
			{
				v_cut_edge = v_cut_edge_success;
				v_cut_depth = v_cut_depth_success;
			}
		}
	}
	if (flags_cut_convex_edge)
	{
		//convex edges in ignored_edges
		for (size_t i = 0; i < v_cut_depth_all_ori.size(); i++)
		{
			std::pair<unsigned, unsigned> tmp_edge_pair(v_cut_edge_all_ori[2 * i], v_cut_edge_all_ori[2 * i + 1]);
			auto it = ignored_edges.find(tmp_edge_pair);
			if (it != ignored_edges.end())
			{
				//not the end
				if (map_is_convex_edge[tmp_edge_pair])
				{
					//convex
					v_cut_edge.push_back(tmp_edge_pair.first);
					v_cut_edge.push_back(tmp_edge_pair.second);
					v_cut_depth.push_back(std::max(v_cut_depth_all_ori[i], min_cut_depth));
					if (pci.is_equivalent_cut(v_cut_edge, v_cut_depth, v_cut_edge_success, v_cut_depth_success))
					{
						continue;
					}
					ig::CVec<double, 3> cur_inner_distortion(-1, -1, -1);
					cur_inner_distortion = cut_and_deform(&pci, v_cut_edge, v_cut_depth, max_volume);
					if (cur_inner_distortion[0] < 0.0)
					{
						//not success
						v_cut_edge.pop_back();
						v_cut_edge.pop_back();
						v_cut_depth.pop_back();
					}
					else
					{
						if (cur_inner_distortion[0] < init_distortion[0])
						{
							//success
							v_cut_edge_success = v_cut_edge;
							v_cut_depth_success = v_cut_depth;
						}
						else
						{
							//current depth is not the best, edge select another depth
							if (v_cut_depth_all_ori[i] > min_cut_depth)
								v_cut_depth.back() = min_cut_depth;
							v_cut_edge_success = v_cut_edge;
							v_cut_depth_success = v_cut_depth;
						}
					}
				}
			}
		}
	}
	if (v_cut_edge.size() == 0)
		return 0;
	bool split_flag = true;
	//cut(&pci, v_cut_edge, v_cut_depth, max_volume, split_flag);
	cut(&pci, v_cut_edge, v_cut_depth, 0.1, split_flag);
	assert(tet_mesh_->tetras.size() == tet_mesh_polycube_->tetras.size());
	//refinement of polycube mesh might need here
	//pci.get_cut_edge_map(cut_edge_o2n);
	get_cut_edge_map();
	return (int)v_cut_edge.size() / 2;
}
bool Polycube_Algorithm::get_cut_edge_map()
{
	if (cut_input_point.empty()) return false;
	//get cut_edge_o2n, tet_mesh_polycube should not be deformed
	cut_edge_o2n.clear();
	std::vector<std::pair<CVec<double, 3>, CVec<double, 3>>> init_feature_polycube_edge_sorted_pts;
	for (size_t i = 0; i < init_feature_polycube_edge_array_sorted.size(); i++)
	{
		unsigned id0 = init_feature_polycube_edge_array_sorted[i].first;
		unsigned id1 = init_feature_polycube_edge_array_sorted[i].second;
		init_feature_polycube_edge_sorted_pts.push_back(std::pair<CVec<double, 3>, CVec<double, 3>>(cut_input_point[id0], cut_input_point[id1]));
	}
	polycube_flattening_interface pq_flatten_tmp;
	int nv = (int)tet_mesh_->tetra_vertices.size();
	std::vector<double> coord_ori_tmp;
	coord_ori_tmp.clear();
	std::vector<TetVertex<double>* > tetra_vertices = tet_mesh_->tetra_vertices;
	for (size_t i = 0; i < nv; i++)
	{
		double x, y, z;
		x = tetra_vertices[i]->pos[0];
		y = tetra_vertices[i]->pos[1];
		z = tetra_vertices[i]->pos[2];
		coord_ori_tmp.push_back(x);
		coord_ori_tmp.push_back(y);
		coord_ori_tmp.push_back(z);
	}
	pq_flatten_tmp.load_deformation_result(coord_ori_tmp, tet_mesh_polycube_);
	int m[6] = { 0, 3, 1, 4, 2, 5 };
	std::vector<int> fix_labels_changed;
	for (size_t i = 0; i < fix_labels.size(); i++)
	{
		if (fix_labels[i] >= 0)
			fix_labels_changed.push_back(m[fix_labels[i]]);
		else
			fix_labels_changed.push_back(-1);
	}
	pq_flatten_tmp.load_boundary_face_label(fix_labels_changed, fix_charts, tet_mesh_polycube_);
	assert(pq_flatten_tmp.polycube_edges.size() >= init_feature_polycube_edge_array_sorted.size());
	for (size_t i = 0; i < init_feature_polycube_edge_array_sorted.size(); i++)
	{
		std::vector<std::pair<unsigned int, unsigned int>> second_array;
		//construct second_array
		std::pair<CVec<double, 3>, CVec<double, 3>> init_two_point = init_feature_polycube_edge_sorted_pts[i];
		for (size_t j = 0; j < pq_flatten_tmp.polycube_edges.size(); j++)
		{
			std::pair<CVec<double, 3>, CVec<double, 3>> tmp_pt_pair(tet_mesh_polycube_->tetra_vertices[pq_flatten_tmp.polycube_edges[j][0]]->pos, tet_mesh_polycube_->tetra_vertices[pq_flatten_tmp.polycube_edges[j][1]]->pos);
			if ((init_two_point.first - tmp_pt_pair.first).L2Norm() < CUT_EDGE_MAP_EPSILON && (init_two_point.second - tmp_pt_pair.second).L2Norm() < CUT_EDGE_MAP_EPSILON)
			{
				second_array.push_back(std::pair<unsigned int, unsigned int>(pq_flatten_tmp.polycube_edges[j][0], pq_flatten_tmp.polycube_edges[j][1]));
			}
			else if ((init_two_point.first - tmp_pt_pair.second).L2Norm() < CUT_EDGE_MAP_EPSILON && (init_two_point.second - tmp_pt_pair.first).L2Norm() < CUT_EDGE_MAP_EPSILON)
			{
				second_array.push_back(std::pair<unsigned int, unsigned int>(pq_flatten_tmp.polycube_edges[j][1], pq_flatten_tmp.polycube_edges[j][0]));
			}
		}
		assert(second_array.size() == 1 || second_array.size() == 2);
		cut_edge_o2n[init_feature_polycube_edge_array_sorted[i]] = second_array;
	}
	return true;
}
bool Polycube_Algorithm::save_polycube_feature_edges_pfeformat(const char* filename)
{
	if (cut_edge_o2n.empty())
	{
		//no cut
		pq_flattening->feature_polycube_edge.clear();
		std::map<std::pair<unsigned, unsigned>, int> pe2id;
		for (int i = 0; i < pq_flattening->polycube_edges.size(); i++)
		{
			unsigned id0 = (unsigned)pq_flattening->polycube_edges[i][0];
			unsigned id1 = (unsigned)pq_flattening->polycube_edges[i][1];
			pe2id[std::pair<unsigned, unsigned>(id0, id1)] = i;
			pe2id[std::pair<unsigned, unsigned>(id1, id0)] = i;
		}
		std::ofstream ofs(filename);
		ofs << "Grouped PolyCube Features" << std::endl;
		ofs << init_feature_polycube_edge_array_ori.size() << std::endl;
		for (size_t i = 0; i < init_feature_polycube_edge_array_ori.size(); i++)
		{
			ofs << "1\n";
			ofs << init_feature_polycube_edge_array_ori[i].first << " " << init_feature_polycube_edge_array_ori[i].second << std::endl;
			auto pqeid0 = pe2id.find(init_feature_polycube_edge_array_ori[i]);
			assert(pqeid0 != pe2id.end());
			if (pqeid0 == pe2id.end()) return false;
			pq_flattening->feature_polycube_edge.push_back(pqeid0->second);
		}
		ofs.close();
	}
	else
	{
		pq_flattening->feature_polycube_edge.clear();
		std::map<std::pair<unsigned, unsigned>, int> pe2id;
		for (int i = 0; i < pq_flattening->polycube_edges.size(); i++)
		{
			unsigned id0 = (unsigned)pq_flattening->polycube_edges[i][0];
			unsigned id1 = (unsigned)pq_flattening->polycube_edges[i][1];
			pe2id[std::pair<unsigned, unsigned>(id0, id1)] = i;
			pe2id[std::pair<unsigned, unsigned>(id1, id0)] = i;
		}
		std::vector<std::vector<std::pair<unsigned, unsigned>>> grouped_polycube_features;
		for (size_t i = 0; i < init_feature_polycube_edge_array_sorted.size(); i++)
		{
			//cut or not
			std::vector<std::pair<unsigned, unsigned>> one_group_feature;
			auto cutit = cut_edge_o2n.find(init_feature_polycube_edge_array_sorted[i]);
			assert(cutit != cut_edge_o2n.end());
			if (cutit != cut_edge_o2n.end())
			{
				//cut
				//assert(cutit->second.size() == 2);
				for (size_t j = 0; j < cutit->second.size(); j++)
				{
					std::pair<unsigned, unsigned> pair0 = cutit->second[j];
					if (init_feature_polycube_edge_array[i].first != init_feature_polycube_edge_array_sorted[i].first)
					{
						//swap
						std::swap(pair0.first, pair0.second);
					}
					auto pqeid0 = pe2id.find(pair0);
					assert(pqeid0 != pe2id.end());
					if (pqeid0 == pe2id.end()) return false;
					pq_flattening->feature_polycube_edge.push_back(pqeid0->second);
					one_group_feature.push_back(pair0);
				}
			}
			else
			{
				return false;
			}
			grouped_polycube_features.push_back(one_group_feature);
		}
		const std::vector<std::vector<int>> &polycube_edges = pq_flattening->polycube_edges;
		std::vector<bool> polycube_vertex_flag(tet_mesh_->tetra_vertices.size(), false);
		for (size_t i = 0; i < polycube_edges.size(); i++)
		{
			polycube_vertex_flag[polycube_edges[i][0]] = true;
			polycube_vertex_flag[polycube_edges[i][1]] = true;
		}
		std::vector<int> vert_o2n(tet_mesh_->tetra_vertices.size(), -1);
		int count = 0;
		for (size_t i = 0; i < tet_mesh_->tetra_vertices.size(); i++)
		{
			if (polycube_vertex_flag[i])
				vert_o2n[i] = count++;
		}
		//save feature file
		std::ofstream ofs(filename);
		ofs << "Grouped PolyCube Features" << std::endl;
		ofs << grouped_polycube_features.size() << std::endl;
		for (size_t i = 0; i < grouped_polycube_features.size(); i++)
		{
			ofs << grouped_polycube_features[i].size() << std::endl;
			for (size_t j = 0; j < grouped_polycube_features[i].size(); j++)
			{
				ofs << vert_o2n[grouped_polycube_features[i][j].first] << " " << vert_o2n[grouped_polycube_features[i][j].second] << std::endl;
			}
		}
		ofs.close();
	}
	return true;
}
void Polycube_Algorithm::randomize_array(std::vector<std::pair<unsigned int, unsigned int>>& pq_edge)
{
	int edgesize = (int)pq_edge.size();
	std::srand(unsigned(std::time(0)));
	std::vector<int> myvector;
	// set some values:
	for (int i = 0; i < edgesize; ++i) myvector.push_back(i); // 1 2 3 4 5 6 7 8 9
	// using built-in random generator:
	std::random_shuffle(myvector.begin(), myvector.end());
	std::vector<std::pair<unsigned int, unsigned int>> pq_edge_new;
	for (size_t i = 0; i < edgesize; i++)
	{
		pq_edge_new.push_back(pq_edge[myvector[i]]);
	}
	pq_edge.clear();
	pq_edge = pq_edge_new;
}
int Polycube_Algorithm::auto_cut_batch_min_max_depth(double cut_ratio, double max_volume, int batch_size, double cutting_edge_ratio, int distortion_type, int ivf_iter, bool min_depth_flag, double min_depth_ratio, bool flags_cut_convex_edge, bool flag_sample_hex, int default_cutnum, double tet_threshold, bool random_cut)
{
	//return the number of cut
	//number of cutting edges: #of available edges * cutting_edge_ratio
	if (!tet_mesh_polycube_)
		return -1;
	if (!do_Extract_polycube_vtk(distortion_type, flag_sample_hex)) return -2;
	double depth_decay_ratio = 0.5;
	int truncated_idx = 0;
	for (size_t i = 0; i < sorted_polycube_edge_distortion.size(); i++)
	{
		if (sorted_polycube_edge_distortion[i] < tet_threshold)
		{
			truncated_idx = (int)i;
			break;
		}
	}
	if (random_cut)
	{
		randomize_array(sorted_polycube_edge);
	}
	if (truncated_idx != 0)
	{
		sorted_polycube_edge.resize(truncated_idx);
	}
	std::vector<ig::CVec<double, 3>> v_input_point = cut_input_point;
	std::vector<unsigned int> v_input_face = cut_input_face;
	std::vector<int> v_input_chart = cut_input_chart;
	std::vector<int> v_input_label = cut_input_label;
	PolyCubeCutInterface pci(v_input_point, v_input_face, v_input_chart, v_input_label);
	std::vector<unsigned int> v_cuttable_edge;
	std::vector<unsigned int> avail_edges;
	std::set<std::pair<unsigned, unsigned>> v_cuttable_edge_set;
	pci.get_cuttable_edge(v_cuttable_edge);
	std::map<std::pair<unsigned int, unsigned int>, bool> map_is_convex_edge;
	pci.get_convex_map(map_is_convex_edge);
	//update one_cut_depth
	double sqrt_two = sqrt(2.0);
	for (size_t i = 0; i < polycube_edge.size(); i++)
	{
		if (polycube_edge_max_cut_depth_one_cut[i] > polycube_edge_neighbor_min_length[i])
			polycube_edge_max_cut_depth_one_cut[i] = polycube_edge_neighbor_min_length[i];
	}
	for (size_t i = 0; i < v_cuttable_edge.size() / 2; i++)
	{
		unsigned id0 = v_cuttable_edge[2 * i];
		unsigned id1 = v_cuttable_edge[2 * i + 1];
		if (id0 < id1)
		{
			v_cuttable_edge_set.insert(std::pair<unsigned, unsigned>(id0, id1));
		}
		else
		{
			v_cuttable_edge_set.insert(std::pair<unsigned, unsigned>(id1, id0));
		}
	}
	std::vector<double> v_cut_depth;
	std::vector<double> v_max_cut_depth;
	std::map<std::pair<unsigned, unsigned>, std::vector<std::pair<unsigned, unsigned>>> map_associated_edge, map_affected_edge;
	pci.get_associated_edge(map_associated_edge);
	//initialize affected
	for (size_t i = 0; i < polycube_edge.size(); i++)
	{
		std::vector<std::pair<unsigned, unsigned>> empty;
		unsigned min_id, max_id;
		min_id = polycube_edge[i].first;
		max_id = polycube_edge[i].second;
		std::pair<unsigned, unsigned> edge_pair(min_id, max_id);
		map_affected_edge[edge_pair] = empty;
	}
	assert(polycube_edge.size() == map_associated_edge.size());
	//construct affected map
	for (size_t i = 0; i < polycube_edge.size(); i++)
	{
		unsigned min_id, max_id;
		min_id = polycube_edge[i].first;
		max_id = polycube_edge[i].second;
		std::pair<unsigned, unsigned> edge_pair(min_id, max_id);
		std::vector<std::pair<unsigned, unsigned>> &assoc = map_associated_edge[edge_pair];
		for (size_t j = 0; j < assoc.size(); j++)
		{
			unsigned assoc0 = assoc[j].first;
			unsigned assoc1 = assoc[j].second;
			std::pair<unsigned, unsigned> assoc_pair(assoc0, assoc1);
			std::vector<std::pair<unsigned, unsigned>> &affected = map_affected_edge[assoc_pair];
			auto it = std::find(affected.begin(), affected.end(), edge_pair);
			if (it == affected.end())
			{
				//not found
				affected.push_back(edge_pair);
			}
		}
	}
	//output
	//min cut depth should be at least half the size of hex cube size, this bound is secure
	//this variable is set both in function auto_cut and select_cut
	double min_cut_depth = min_depth_ratio * sigma_r * min_polycube_edge_length;
	int max_edge_num = (int)sorted_polycube_edge.size();
	//get edges which can be cut deeply enough
	std::vector<unsigned int> v_cut_edge_all, v_cut_edge_all_ori, v_cut_edge;
	std::vector<double> v_cut_depth_all, v_cut_depth_all_ori;
	for (size_t i = 0; i < sorted_polycube_edge.size(); i++)
	{
		//only consider edge which is already long enough
		unsigned int id0 = sorted_polycube_edge[i].first;
		unsigned int id1 = sorted_polycube_edge[i].second;
		v_cut_edge_all_ori.push_back(sorted_polycube_edge[i].first);
		v_cut_edge_all_ori.push_back(sorted_polycube_edge[i].second);
	}
	//std::vector<unsigned int> v_cut_edge{ 2, 1, 0, 3 };
	get_cut_depth_min(v_cut_edge_all_ori, v_cut_depth_all_ori, cut_ratio);
	std::cout << "min cut depth: " << min_cut_depth << std::endl;
	std::set<std::pair<unsigned, unsigned>> ignored_edges;
	//construct ignored_edges
	for (size_t i = 0; i < sorted_polycube_edge.size(); i++)
	{
		if (v_cut_depth_all_ori[i] < min_cut_depth)
		{
			//insert this edge and the affected to ignore_edges
			auto search = ignored_edges.find(sorted_polycube_edge[i]);
			if (search == ignored_edges.end())
			{
				//not found
				ignored_edges.insert(sorted_polycube_edge[i]);
				std::vector<std::pair<unsigned, unsigned>> &affected = map_affected_edge[sorted_polycube_edge[i]];
				for (size_t j = 0; j < affected.size(); j++)
				{
					ignored_edges.insert(affected[j]);
				}
			}
		}
	}
	//ignore edges with vertex of valence more than 3
	std::vector<std::set<unsigned>> vert_con_chart(v_input_point.size());
	for (size_t i = 0; i < v_input_face.size() / 3; i++)
	{
		unsigned tmpid[3] = { v_input_face[3 * i], v_input_face[3 * i + 1], v_input_face[3 * i + 2] };
		vert_con_chart[tmpid[0]].insert(v_input_chart[i]);
		vert_con_chart[tmpid[0]].insert(v_input_chart[i]);
		vert_con_chart[tmpid[1]].insert(v_input_chart[i]);
		vert_con_chart[tmpid[1]].insert(v_input_chart[i]);
		vert_con_chart[tmpid[2]].insert(v_input_chart[i]);
		vert_con_chart[tmpid[2]].insert(v_input_chart[i]);
	}
	for (size_t i = 0; i < vert_con_chart.size(); i++)
	{
		if (vert_con_chart[i].size() > 3)
		{
			for (auto nbv : vert_con_chart[i])
			{
				unsigned minid = std::min((unsigned)i, nbv);
				unsigned maxid = std::max((unsigned)i, nbv);
				std::pair<unsigned, unsigned> tmppair(minid, maxid);
				auto search = ignored_edges.find(tmppair);
				if (search == ignored_edges.end())
				{
					ignored_edges.insert(tmppair);
					std::vector<std::pair<unsigned, unsigned>> &affected = map_affected_edge[tmppair];
					for (size_t j = 0; j < affected.size(); j++)
					{
						ignored_edges.insert(affected[j]);
					}
				}
			}
		}
	}
	int avail_edge_count = 0;
	//remove convex and concave edge from ignored_edges
	std::vector<std::pair<unsigned, unsigned>> convex_concave;
	for (auto it = ignored_edges.begin(); it != ignored_edges.end(); ++it)
	{
		auto cuttable_it = v_cuttable_edge_set.find(*it);
		if (cuttable_it != v_cuttable_edge_set.end())
		{
			//convex or concave
			//1 cut depth must be large enough
			unsigned tmp_idx = polycube_edge_to_idx_map[(*it)];
			double tmp_depth = polycube_edge_max_cut_depth_one_cut[tmp_idx] * cut_ratio;
			if (tmp_depth > min_cut_depth)
			{
				convex_concave.push_back((*it));
			}
		}
	}
	for (size_t i = 0; i < convex_concave.size(); i++)
	{
		ignored_edges.erase(convex_concave[i]);
	}
	int n_ignored_edges_first_step = (int)ignored_edges.size();
	//ignore edges modified by feature edges here
	if (!init_feature_polycube_edge_array_sorted.empty())
	{
		for (size_t i = 0; i < init_feature_polycube_edge_array_sorted.size(); i++)
		{
			ignored_edges.insert(init_feature_polycube_edge_array_sorted[i]);
		}
	}
	get_available_cut(&pci, min_cut_depth, max_volume, sorted_polycube_edge, ignored_edges, avail_edges);
	//if first step ignored edge and available edges are empty, which means all polycube edges are feature edges
	if (n_ignored_edges_first_step == 0 && avail_edges.size() == 0)
		return -3;
	//depth need to be adjusted here
	v_cut_depth.clear();
	v_cut_edge.clear();
	v_cut_depth_all.clear();
	v_cut_edge_all.clear();
	int final_cut_candidate_num = (int)(cutting_edge_ratio * (avail_edges.size() / 2));
	if (default_cutnum > 0 && default_cutnum <= avail_edges.size() / 2)
	{
		final_cut_candidate_num = default_cutnum;
	}
	std::vector<unsigned> v_cut_edge_success;
	std::vector<double> v_cut_depth_success;
	v_cut_edge_success.clear();
	v_cut_depth_success.clear();
	for (size_t i = 0; i < final_cut_candidate_num; i++)
	{
		v_cut_edge_all.push_back(avail_edges[2 * i]);
		v_cut_edge_all.push_back(avail_edges[2 * i + 1]);
	}
	get_cut_depth(v_cut_edge_all, v_cut_depth_all, cut_ratio, map_associated_edge, (int)cut_input_point.size());
	if (batch_size <= 0 || batch_size > final_cut_candidate_num)
	{
		batch_size = final_cut_candidate_num;
	}
	//int n_batch = (v_cut_depth.size() / batch_size) + 1;
	int n_batch = (int)(std::ceil((double)final_cut_candidate_num / (double)batch_size));
	if (final_cut_candidate_num > 0)
	{
		for (size_t i = 0; i < n_batch; i++)
		{
			//insert elem
			int n_add = batch_size;
			if (batch_size * (i + 1) > final_cut_candidate_num)
			{
				n_add = (int)(final_cut_candidate_num - batch_size * i);
			}
			assert(n_add <= batch_size);
			for (size_t j = 0; j < n_add; j++)
			{
				int tmp_idx = (int)(i * batch_size + j);
				v_cut_edge.push_back(v_cut_edge_all[2 * tmp_idx]);
				v_cut_edge.push_back(v_cut_edge_all[2 * tmp_idx + 1]);
				v_cut_depth.push_back(v_cut_depth_all[tmp_idx]);
			}
			//equal
			if (pci.is_equivalent_cut(v_cut_edge, v_cut_depth, v_cut_edge_success, v_cut_depth_success))
			{
				continue;
			}
			bool cut_success = true;
			std::vector<double> max_cut_depth(n_add, -1.0);
			//cut
			ig::CVec<double, 3> cur_inner_distortion(-1, -1, -1);
			while (cur_inner_distortion[0] < 0.0 || cur_inner_distortion[distortion_type] > init_distortion[distortion_type])
			{
				cur_inner_distortion = cut_and_deform(&pci, v_cut_edge, v_cut_depth, max_volume, ivf_iter);
				int change_flag = -1;
				if (cur_inner_distortion[0] > 0.0)
				{
					if (max_cut_depth[0] < 0.0)
					{
						for (size_t j = 0; j < n_add; j++)
						{
							int tmp_idx = (int)(i * batch_size + j);
							//v_cut_depth.push_back(v_cut_depth_all[tmp_idx]);
							max_cut_depth[j] = v_cut_depth[tmp_idx]; //init cut depth
						}
					}
				}
				//not success
				if (cur_inner_distortion[0] < 0.0 || cur_inner_distortion[distortion_type] > init_distortion[distortion_type])
				{
					change_flag = 0;
					for (size_t k = 0; k < v_cut_depth.size(); k++)
					{
						v_cut_depth[k] = v_cut_depth[k] * depth_decay_ratio;
						if (v_cut_depth[k] <= min_cut_depth)
						{
							v_cut_depth[k] = min_cut_depth;
						}
						else
						{
							change_flag++;
						}
					}
					
				}
				if (change_flag == 0) // cut illegal but nothing changes
				{
					cut_success = false;
					break;
				}
			}
			if (cut_success)
			{
				v_cut_edge_success = v_cut_edge;
				v_cut_depth_success = v_cut_depth;
			}
			else
			{
				for (size_t j = 0; j < n_add; j++)
				{
					int tmp_idx = (int)(i * batch_size + j);
					if (!min_depth_flag && max_cut_depth[j] > 0.0)
					{
						//max cut depth
						v_cut_depth[tmp_idx] = max_cut_depth[j];
					}
					else
					{
						v_cut_depth[tmp_idx] = min_cut_depth;
					}
				}
				v_cut_edge_success = v_cut_edge;
				v_cut_depth_success = v_cut_depth;
			}
		}
	}
	if (flags_cut_convex_edge)
	{
		for (size_t i = 0; i < v_cut_depth_all_ori.size(); i++)
		{
			std::pair<unsigned, unsigned> tmp_edge_pair(v_cut_edge_all_ori[2 * i], v_cut_edge_all_ori[2 * i + 1]);
			auto it = ignored_edges.find(tmp_edge_pair);
			if (it != ignored_edges.end())
			{
				//not the end
				if (map_is_convex_edge[tmp_edge_pair])
				{
					//convex
					v_cut_edge.push_back(tmp_edge_pair.first);
					v_cut_edge.push_back(tmp_edge_pair.second);
					v_cut_depth.push_back(std::max(v_cut_depth_all_ori[i], min_cut_depth));
					if (pci.is_equivalent_cut(v_cut_edge, v_cut_depth, v_cut_edge_success, v_cut_depth_success))
					{
						continue;
					}
					ig::CVec<double, 3> cur_inner_distortion(-1, -1, -1);
					cur_inner_distortion = cut_and_deform(&pci, v_cut_edge, v_cut_depth, max_volume);
					if (cur_inner_distortion[0] < 0.0)
					{
						//not success
						v_cut_edge.pop_back();
						v_cut_edge.pop_back();
						v_cut_depth.pop_back();
					}
					else
					{
						if (cur_inner_distortion[0] < init_distortion[0])
						{
							//success
							v_cut_edge_success = v_cut_edge;
							v_cut_depth_success = v_cut_depth;
						}
						else
						{
							//current depth is not the best, edge select another depth
							if (v_cut_depth_all_ori[i] > min_cut_depth)
								v_cut_depth.back() = min_cut_depth;
							v_cut_edge_success = v_cut_edge;
							v_cut_depth_success = v_cut_depth;
						}
					}
				}
			}
		}
	}
	if (v_cut_edge.size() == 0)
		return 0;
	bool split_flag = true;
	cut(&pci, v_cut_edge, v_cut_depth, max_volume, split_flag);
	assert(tet_mesh_->tetras.size() == tet_mesh_polycube_->tetras.size());
	get_cut_edge_map();
	return (int)v_cut_edge.size() / 2;
}
ig::CVec<double, 3> Polycube_Algorithm::cut(PolyCubeCutInterface *pci, std::vector<unsigned int> &v_cut_edge, std::vector<double> &v_cut_depth, double max_volume_ratio, bool split_tet)
{
	ig::CVec<double, 3> cur_distortion(-1, -1, -1);
	int cut_num = (int)v_cut_edge.size() / 2;
	std::vector<unsigned int> v_cut_vertex;
	double max_volume = max_volume_ratio * avg_cell_volume;
	double tet_quality = DEFAULT_TET_QUALITY;
	std::vector<ig::CVec<double, 3>> v_output_point;
	std::vector<unsigned int> v_output_tet;
	std::vector<int> v_output_chart;
	std::vector<int> v_output_label;
	std::vector<int> v_cut_type;
	std::vector<int> v_cut_section_chart;
	std::vector<int> v_seam_edge_vertex_id;
	std::vector<std::vector<int>> vv_congruent_face;
	std::vector<ThreeCut> v_three_cut;
	double decay_ratio = 0.8; //cut depth decay ratio
	bool flag = true;
	try
	{
		pci->cut(
			v_cut_edge, v_cut_depth, max_volume, tet_quality,
			v_output_point, v_output_tet, v_output_chart, v_output_label, v_cut_type,
			v_cut_section_chart, v_seam_edge_vertex_id, vv_congruent_face, v_three_cut
		);
	}
	catch (const std::string &s)
	{
		flag = false;
		std::cerr << s;
		return cur_distortion;
	}
	//polycube needs to be updated
	tet_mesh_polycube_->load_tet(v_output_point, v_output_tet);
	if (split_tet)
	{
		//split tet, boundary kept unchanged update v_output_chart, v_output_label
		int max_chart_id = -1;
		for (size_t i = 0; i < v_output_chart.size(); i++)
		{
			if (max_chart_id < v_output_chart[i])
				max_chart_id = v_output_chart[i];
		}
		std::vector<int> chart2label(max_chart_id + 1, -1);
		for (size_t i = 0; i < v_output_chart.size(); i++)
		{
			if (v_output_chart[i] >= 0)
			{
				assert(v_output_label[i] >= 0);
				chart2label[v_output_chart[i]] = v_output_label[i];
			}
		}
		//the following two variable kept unchanged
		std::vector<std::vector<int>> all_boundary_faces_before; //sorted
		std::map<std::vector<int>, int> bf2chart_before;
		std::vector<double> bdpts;
		std::vector<int> bdfaces_before, bdfaces_after, s2v;
		std::map<std::vector<int>, int> bf2cell_before;
		TetOperation<double> tet_op;
		tet_op.get_boundary_vert_face(tet_mesh_polycube_, bdpts, bdfaces_before, s2v, true, &bf2cell_before);
		//set bf2chart_before
		for (auto it : bf2cell_before)
		{
			bf2chart_before[it.first] = v_output_chart[it.second];
		}
		//split tet
		tet_op.set_tet_structure(tet_mesh_polycube_);
		tet_op.split_inner_edge_two_boundary_point(*tet_mesh_polycube_);
		//update chart label
		std::map<std::vector<int>, int> bf2cell_after;
		tet_op.get_boundary_vert_face(tet_mesh_polycube_, bdpts, bdfaces_after, s2v, true, &bf2cell_after);
		assert(bdfaces_before.size() == bdfaces_after.size());
		v_output_chart.clear();
		v_output_label.clear();
		v_output_chart.resize(tet_mesh_polycube_->tetras.size(), -1);
		v_output_label.resize(tet_mesh_polycube_->tetras.size(), -1);
		for (auto it : bf2cell_after)
		{
			auto findit = bf2chart_before.find(it.first);
			assert(findit != bf2chart_before.end());
			int chart = findit->second;
			v_output_chart[it.second] = chart;
			v_output_label[it.second] = chart2label[chart];
		}
	}
	if (!construct_original_tet()) return cur_distortion;
	printf("Prepare Deformation.......\n");
	delete pq_deform;
	delete pq_flattening;
	pq_deform = new polycube_deformation_interface();
	pq_flattening = new polycube_flattening_interface();
	reset_aLL_state_vtk();
	(*pq_deform).prepare_for_deformation(tet_mesh_polycube_);
	coord_ori.clear();
	fix_charts.clear();
	fix_labels.clear();
	fix_charts = v_output_chart;
	fix_labels = v_output_label;
	(*pq_deform).set_fixed_chart_label(fix_charts, fix_labels, false);
	load_equal_faces(v_cut_type, v_cut_section_chart, v_seam_edge_vertex_id, vv_congruent_face, (int)tet_mesh_polycube_->tetra_vertices.size(), v_three_cut);
	compute_flips();
	cur_distortion = (*pq_deform).compute_distortion(tet_mesh_polycube_);
	return cur_distortion;
}
bool Polycube_Algorithm::is_edges_cuttable(PolyCubeCutInterface *pci, std::vector<unsigned int> &v_cut_edge, std::vector<double> &v_cut_depth, double max_volume_ratio)
{
	//only decide whether a cut is available
	bool cut_success = false;
	int cut_num = (int)v_cut_edge.size() / 2;
	std::vector<unsigned int> v_cut_vertex;
	double max_volume = max_volume_ratio * avg_cell_volume;
	double tet_quality = DEFAULT_TET_QUALITY;
	//double tet_quality = 2.0;
	std::vector<ig::CVec<double, 3>> v_output_point;
	std::vector<unsigned int> v_output_tet;
	std::vector<int> v_output_chart;
	std::vector<int> v_output_label;
	std::vector<int> v_cut_type;
	std::vector<int> v_cut_section_chart;
	std::vector<int> v_seam_edge_vertex_id;
	std::vector<std::vector<int>> vv_congruent_face;
	std::vector<ThreeCut> v_three_cut;
	double decay_ratio = 0.8; //cut depth decay ratio
	bool flag = true;
	try
	{
		pci->cut_test(
			v_cut_edge, v_cut_depth, max_volume, tet_quality
		);
	}
	catch (const std::string &s)
	{
		flag = false;
		std::cerr << s;
		//return;
		return cut_success;
	}
	
	return true;
}
bool Polycube_Algorithm::get_available_cut(PolyCubeCutInterface *pci, double min_cut_depth, double max_volume_ratio, const std::vector<std::pair<unsigned, unsigned>> &input_edges, const std::set<std::pair<unsigned, unsigned>> &ignore_edge_set, std::vector<unsigned> &avail_edges)
{
	//set all depth to min_cut_depth, get avail_edges from input_edges
	avail_edges.clear();
	std::vector<double> v_cut_depth;
	int cut_num = (int)input_edges.size();
	std::vector<unsigned int> v_cut_vertex;
	double max_volume = max_volume_ratio * avg_cell_volume;
	double tet_quality = DEFAULT_TET_QUALITY;
	std::vector<ig::CVec<double, 3>> v_output_point;
	std::vector<unsigned int> v_output_tet;
	std::vector<int> v_output_chart;
	std::vector<int> v_output_label;
	std::vector<int> v_cut_type;
	std::vector<int> v_cut_section_chart;
	std::vector<int> v_seam_edge_vertex_id;
	std::vector<std::vector<int>> vv_congruent_face;
	std::vector<ThreeCut> v_three_cut;
	std::vector<unsigned> success_cut_edge;
	std::vector<double> success_cut_depth;
	pci->is_equivalent_cut(avail_edges, v_cut_depth, success_cut_edge, success_cut_depth);
	for (size_t i = 0; i < cut_num; i++)
	{
		unsigned id1 = input_edges[i].first;
		unsigned id2 = input_edges[i].second;
		auto it = ignore_edge_set.find(input_edges[i]);
		if (it != ignore_edge_set.end())
		{
			//found in the set
			continue;
		}
		avail_edges.push_back(id1);
		avail_edges.push_back(id2);
		v_cut_depth.push_back(min_cut_depth);
		if (pci->is_equivalent_cut(avail_edges, v_cut_depth, success_cut_edge, success_cut_depth))
		{
			//equivalent cut
			continue;
		}
		bool flag = true;
		try
		{
			pci->cut_test(
				avail_edges, v_cut_depth, max_volume, tet_quality
			);
		}
		catch (const std::string &s)
		{
			flag = false;
			std::cerr << s;
			//return;
		}
		if (!flag)
		{
			//not success
			avail_edges.pop_back();
			avail_edges.pop_back();
			v_cut_depth.pop_back();
		}
		else
		{
			//cut success
			success_cut_edge = avail_edges;
			success_cut_depth = v_cut_depth;
		}
	}
	if (avail_edges.size() == 0)
		return false;
	else
		return true;
}
ig::CVec<double, 3> Polycube_Algorithm::cut_and_deform(PolyCubeCutInterface *pci, std::vector<unsigned int> &v_cut_edge, std::vector<double> &v_cut_depth, double sigma_s, double max_volume_ratio, int ivf_iter)
{
	ig::CVec<double, 3> cur_distortion(-1, -1, -1);
	int cut_num = (int)v_cut_edge.size() / 2;
	std::vector<unsigned int> v_cut_vertex;
	double max_volume = max_volume_ratio * avg_cell_volume;
	double tet_quality = DEFAULT_TET_QUALITY;
	std::vector<ig::CVec<double, 3>> v_output_point;
	std::vector<unsigned int> v_output_tet;
	std::vector<int> v_output_chart;
	std::vector<int> v_output_label;
	std::vector<int> v_cut_type;
	std::vector<int> v_cut_section_chart;
	std::vector<int> v_seam_edge_vertex_id;
	std::vector<std::vector<int>> vv_congruent_face;
	std::vector<ThreeCut> v_three_cut;
	double decay_ratio = 0.8; //cut depth decay ratio
	bool flag = true;
	try
	{
		pci->cut(
			v_cut_edge, v_cut_depth, max_volume, tet_quality,
			v_output_point, v_output_tet, v_output_chart, v_output_label, v_cut_type,
			v_cut_section_chart, v_seam_edge_vertex_id, vv_congruent_face, v_three_cut
		);
	}
	catch (const std::string &s)
	{
		flag = false;
		std::cerr << s;
		return cur_distortion;
	}
	tet_mesh_polycube_->load_tet(v_output_point, v_output_tet);
	if (!construct_original_tet())
	{
		std::cout << "!!!!construct ori tet failed" << std::endl;
		return cur_distortion;
	}
	printf("Prepare Deformation.......\n");
	delete pq_deform;
	delete pq_flattening;
	pq_deform = new polycube_deformation_interface();
	pq_flattening = new polycube_flattening_interface();
	reset_aLL_state_vtk();
	(*pq_deform).prepare_for_deformation(tet_mesh_polycube_);
	coord_ori.clear();
	fix_charts.clear();
	fix_labels.clear();
	fix_charts = v_output_chart;
	fix_labels = v_output_label;
	(*pq_deform).set_fixed_chart_label(fix_charts, fix_labels, false);
	load_equal_faces(v_cut_type, v_cut_section_chart, v_seam_edge_vertex_id, vv_congruent_face, (int)v_output_point.size(), v_three_cut);
	compute_flips();
	auto_deformation(sigma_s, ivf_iter);
	cur_distortion = (*pq_deform).compute_distortion(tet_mesh_polycube_);
	return cur_distortion;
}
void Polycube_Algorithm::auto_deformation(double sigma_s, int iter_ivf)
{
	//tet_mesh_polycube required
	if (tet_mesh_polycube_ == NULL) return;
	do_IVF_seamless_moving(0, iter_ivf);
	std::pair<int, int> flip_pair = compute_flips();
	if (flip_pair.first != 0)
	{
		//flips not removed totally
		do_IVF_seamless_svd(10);
		flip_pair = compute_flips();
		if (flip_pair.first != 0)
		{
			do_IVF_seamless_moving_part(0, iter_ivf, 0);
		}
	}
	pq_deform->load_ori_tet(tet_mesh_);
	init_distortion = pq_deform->compute_distortion(tet_mesh_polycube_);
	double avg_diff_th = 0.01;
	double max_diff_th = 0.2;
	int max_iter = 2, iter = 0;
	std::pair<double, double> avg_max_diff;
	do_deformation_constrained_vtk(sigma_s);
	avg_max_diff = pq_deform->compute_chart_normal_diff();
	while (avg_max_diff.first > avg_diff_th || avg_max_diff.second > max_diff_th)
	{
		do_deformation_constrained_continued();
		avg_max_diff = pq_deform->compute_chart_normal_diff();
		iter++;
		if (iter > max_iter)
		{
			break;
		}
	}
	pq_deform->compute_distortion(tet_mesh_polycube_);
	do_optimize_polycube_PC(true, 1, 2);
	compute_flips();
}
void Polycube_Algorithm::get_cut_depth_min(const std::vector<unsigned int> &cut_edge, std::vector<double> &cut_depth, double cut_ratio)
{
	cut_depth.clear();
	if (cut_edge.size() == 0)
		return;
	int cut_edge_num = (int)cut_edge.size() / 2;
	for (size_t i = 0; i < cut_edge_num; i++)
	{
		unsigned edge_idx = -1;
		edge_idx = polycube_edge_to_idx_map[std::pair<unsigned, unsigned>(cut_edge[2 * i], cut_edge[2 * i + 1])];
		assert(edge_idx >= 0);
		double init_depth = polycube_edge_max_cut_depth_three_cut[edge_idx];
		if (init_depth > polycube_edge_max_cut_depth_one_cut[edge_idx])
			init_depth = polycube_edge_max_cut_depth_one_cut[edge_idx];
		double depth = init_depth * cut_ratio;
		cut_depth.push_back(depth);
	}
}
void Polycube_Algorithm::get_cut_depth(const std::vector<unsigned int> &cut_edge, std::vector<double> &cut_depth, double cut_ratio, const std::map<std::pair<unsigned, unsigned>, std::vector<std::pair<unsigned, unsigned>>>& map_associated_edge, int n_vert)
{
	cut_depth.clear();
	if (cut_edge.size() == 0)
		return;
	//first input edges should be extended
	std::set<std::pair<unsigned, unsigned>> extended_edges;
	int cut_edge_num = (int)cut_edge.size() / 2;
	for (size_t i = 0; i < cut_edge_num; i++)
	{
		std::pair<unsigned, unsigned> tmp_edge(cut_edge[2 * i], cut_edge[2 * i + 1]);
		extended_edges.insert(tmp_edge);
		const std::vector<std::pair<unsigned, unsigned>> &associate_edges = map_associated_edge.at(tmp_edge);
		for (size_t j = 0; j < associate_edges.size(); j++)
		{
			extended_edges.insert(associate_edges[j]);
		}
	}
	//divide edges into one-cut edge and three-cut edge
	std::map<std::pair<unsigned, unsigned>, bool> map_is_one_cut;
	std::vector<int> v_vert_count(n_vert, 0);
	for (auto it : extended_edges)
	{
		unsigned id1 = it.first;
		unsigned id2 = it.second;
		v_vert_count[id1]++;
		v_vert_count[id2]++;
	}
	for (auto it : extended_edges)
	{
		unsigned id1 = it.first;
		unsigned id2 = it.second;
		if (v_vert_count[id1] == 1 && v_vert_count[id2] == 1)
		{
			//one cut
			map_is_one_cut[it] = true;
		}
		else
		{
			map_is_one_cut[it] = false;
		}
	}
	//construct map_edge_to_distance
	for (size_t i = 0; i < cut_edge_num; i++)
	{
		std::pair<unsigned, unsigned> tmp_edge(cut_edge[2 * i], cut_edge[2 * i + 1]);
		unsigned edge_idx = polycube_edge_to_idx_map[tmp_edge];
		if (map_is_one_cut[tmp_edge])
			cut_depth.push_back(polycube_edge_max_cut_depth_one_cut[edge_idx] * cut_ratio);
		else
			cut_depth.push_back(polycube_edge_max_cut_depth_three_cut[edge_idx] * cut_ratio);
	}
}
bool Polycube_Algorithm::reset_max_cut_depth()
{
	polycube_edge_max_cut_depth_one_cut.clear();
	polycube_edge_max_cut_depth_three_cut.clear();
	//polycube_edge_max_cut_depth_one_cut.resize(polycube_edge.size(), -1.0);
	int n_sample = 10; // larger than 2
	static CVec<double, 3> n[6];
	for (size_t i = 0; i < 6; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			n[i][j] = 0.0;
		}
	}
	n[0][0] = -1.0;
	n[1][0] = 1.0;
	n[2][1] = -1.0;
	n[3][1] = 1.0;
	n[4][2] = -1.0;
	n[5][2] = 1.0;
	std::vector<std::pair<CVec<double, 3>, CVec<double, 3>>> rays;
	std::vector<std::set<unsigned>> ignore_faces;
	int n_faces = (int)cut_input_face.size() / 3;
	//vertex face map
	std::vector<std::set<unsigned>> vf_array;
	vf_array.resize(cut_input_point.size());
	for (size_t i = 0; i < n_faces; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			vf_array[cut_input_face[3 * i + j]].insert((unsigned)i);
		}
	}
	for (size_t i = 0; i < polycube_edge.size(); i++)
	{
		std::vector<ig::CVec<double, 3>> begin_points;
		begin_points.resize(n_sample);
		CVec<double, 3> p0, p1, diff;
		p0 = cut_input_point[polycube_edge[i].first];
		p1 = cut_input_point[polycube_edge[i].second];
		diff = (p1 - p0) / (1.0 * (n_sample + 1));
		//begin_points[0] = p0;
		for (size_t j = 0; j < n_sample; j++)
		{
			begin_points[j] = p0 + (j + 1) * diff;
		}
		CVec<double, 3> dir;
		dir = n[polycube_edge_label[i].first] + n[polycube_edge_label[i].second];
		//dir = -dir;
		for (size_t j = 0; j < n_sample; j++)
		{
			rays.push_back(std::pair<CVec<double, 3>, CVec<double, 3>>(begin_points[j], dir));
		}
		unsigned int id0 = polycube_edge[i].first;
		unsigned int id1 = polycube_edge[i].second;
		std::vector<unsigned> ignore_face_one_dir;
		for (size_t j = 0; j < n_faces; j++)
		{
			int n_equal = 0;
			for (size_t k = 0; k < 3; k++)
			{
				if (id0 == cut_input_face[3 * j + k])
					n_equal++;
				if (id1 == cut_input_face[3 * j + k])
					n_equal++;
			}
			if (n_equal == 2)
				ignore_face_one_dir.push_back((unsigned)j);
		}
		assert(ignore_face_one_dir.size() == 2);
		std::set<unsigned> one_set;
		one_set.insert(ignore_face_one_dir[0]);
		one_set.insert(ignore_face_one_dir[1]);
		for (size_t j = 0; j < n_sample; j++)
		{
			ignore_faces.push_back(one_set);
		}
	}
	//rays and ignore faces finished.
	std::vector<double> distances;
	if (!dist_rays2mesh(rays, ignore_faces, distances))
		return false;
	for (size_t i = 0; i < polycube_edge.size(); i++)
	{
		double min_dist = distances[i * n_sample];
		for (size_t j = 1; j < n_sample; j++)
		{
			if (min_dist > distances[i * n_sample + j])
				min_dist = distances[i * n_sample + j];
		}
		polycube_edge_max_cut_depth_one_cut.push_back(min_dist);
	}
	//set polycube_edge_max_cut_depth_three_cut
	std::vector<std::vector<int>> vert_label, vert_chart;
	vert_label.resize(cut_input_point.size());
	vert_chart.resize(cut_input_point.size());
	for (size_t i = 0; i < polycube_edge.size(); i++)
	{
		unsigned int id0 = polycube_edge[i].first;
		unsigned int id1 = polycube_edge[i].second;
		vert_label[id0].push_back(polycube_edge_label[i].first);
		vert_label[id0].push_back(polycube_edge_label[i].second);
		vert_label[id1].push_back(polycube_edge_label[i].first);
		vert_label[id1].push_back(polycube_edge_label[i].second);
	}
	for (size_t i = 0; i < n_faces; i++)
	{
		vert_chart[cut_input_face[3 * i]].push_back((int)i);
		vert_chart[cut_input_face[3 * i + 1]].push_back((int)i);
		vert_chart[cut_input_face[3 * i + 2]].push_back((int)i);
	}
	rays.clear();
	ignore_faces.clear();
	distances.clear();
	for (size_t i = 0; i < cut_input_point.size(); i++)
	{
		assert(vert_label[i].size() == 6);
		CVec<double, 3> dir(0.0, 0.0, 0.0);
		for (size_t j = 0; j < 6; j++)
		{
			dir = dir + n[vert_label[i][j]];
		}
		rays.push_back(std::pair<CVec<double, 3>, CVec<double, 3>>(cut_input_point[i], dir));
		std::set<unsigned> one_set;
		for (size_t j = 0; j < vert_chart[i].size(); j++)
		{
			one_set.insert((unsigned)vert_chart[i][j]);
		}
		ignore_faces.push_back(one_set);
	}
	if (!dist_rays2mesh(rays, ignore_faces, distances))
		return false;
	polycube_edge_max_cut_depth_three_cut.resize(polycube_edge.size(), -1.0);
	double ratio = sqrt(2.0 / 3.0);
	for (size_t i = 0; i < polycube_edge.size(); i++)
	{
		unsigned int id0 = polycube_edge[i].first;
		unsigned int id1 = polycube_edge[i].second;
		double tmp_min = std::min(distances[id0], distances[id1]);
		polycube_edge_max_cut_depth_three_cut[i] = tmp_min * ratio;
	}
	return true;
}
bool Polycube_Algorithm::dist_rays2mesh(const std::vector<std::pair<CVec<double, 3>, CVec<double, 3>>> &rays, const std::vector<std::set<unsigned>> &ignore_faces, std::vector<double> &distance_array)
{
	distance_array.clear();
	std::vector<CVec<double, 3>> first_intersection_points;
	ig::CGALHelper::rays_triangles_intersection(rays, cut_input_point, cut_input_face, ignore_faces, first_intersection_points);
	if (rays.size() != first_intersection_points.size())
		return false;
	for (size_t i = 0; i < rays.size(); i++)
	{
		distance_array.push_back((rays[i].first - first_intersection_points[i]).L2Norm());
	}
	return true;
}
void Polycube_Algorithm::save_polycube_vtk(const char* filename)
{
	tet_mesh_polycube_->save_vtk(filename);
}
void Polycube_Algorithm::save_polycube_scaled_vtk(const char* filename)
{
	//tet_mesh_polycube_->save_vtk(filename);
	double cube_len = final_cube_length;
	if (final_cube_length < 0.0)
		cube_len = 1.0;
	const std::vector<Tetrahedron<double>*> &tetras = tet_mesh_polycube_->tetras;
	const std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_polycube_->tetra_vertices;
	std::ofstream vtkfile(filename);
	if (!vtkfile.is_open()) return;
	vtkfile.precision(16);
	vtkfile << "# vtk DataFile Version 2.0" << std::endl
		<< "Tetrahedral Mesh" << std::endl
		<< "ASCII" << std::endl
		<< "DATASET UNSTRUCTURED_GRID\n"
		<< "POINTS " << tetra_vertices.size() << " double\n";
	for (size_t i = 0; i < tetra_vertices.size(); i++)
	{
		vtkfile << std::scientific << tetra_vertices[i]->pos / cube_len << std::endl;
	}
	vtkfile << "CELLS " << tetras.size() << " " << tetras.size() * 5 << "\n";
	for (size_t i = 0; i < tetras.size(); i++)
	{
		vtkfile << "4 " << tetras[i]->vertex[0]->id << " "
			<< tetras[i]->vertex[1]->id << " "
			<< tetras[i]->vertex[2]->id << " "
			<< tetras[i]->vertex[3]->id << std::endl;
	}
	vtkfile << "CELL_TYPES " << tetras.size() << std::endl;
	for (size_t i = 0; i < tetras.size(); i++)
		vtkfile << 10 << std::endl;
	vtkfile.close();
}
void Polycube_Algorithm::save_tet_vtk(const char* filename)
{
	if (tet_mesh_)
	{
		tet_mesh_->save_vtk(filename);
	}
	else if (mesh_)
	{
		//ovm to vtk
		std::vector<double> x, y, z;
		std::vector<std::vector<int>> indices;
		for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
		{
			int v_id = v_it->idx();
			OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(*v_it);
			x.push_back(p[0]);
			y.push_back(p[1]);
			z.push_back(p[2]);
		}
		std::vector<OpenVolumeMesh::HalfFaceHandle> hff_vec;
		for (OpenVolumeMesh::CellIter c_it = mesh_->cells_begin(); c_it != mesh_->cells_end(); ++c_it)
		{
			double cv_count = 0.0; int c_id = c_it->idx();
			std::vector<int> tmp_indices;
			const std::vector<OpenVolumeMesh::VertexHandle>& vertices_ = mesh_->cell(*c_it).vertices();
			hff_vec = mesh_->cell(*c_it).halffaces();
			assert(hff_vec.size() == 4);
			OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hff_vec[0]);
			std::set<int> face_set;
			for (hfv_it; hfv_it; ++hfv_it)
			{
				tmp_indices.push_back(hfv_it->idx());
				face_set.insert(hfv_it->idx());
			}
			for (unsigned i = 0; i < vertices_.size(); ++i)
			{
				int v_id = vertices_[i];
				auto it = face_set.find(v_id);
				if (it == face_set.end())
				{
					tmp_indices.push_back(v_id);
				}
			}
			assert(tmp_indices.size() == 4);
			indices.push_back(tmp_indices);
		}
		std::ofstream outputfile(filename);
		outputfile << "# vtk DataFile Version 3.0\n"
			<< "Volume mesh\n"
			<< "ASCII\n"
			<< "DATASET UNSTRUCTURED_GRID\n";
		int n_point = (int)x.size();
		int n_cell = (int)indices.size();
		outputfile << "POINTS " << n_point << " double" << std::endl;
		for (size_t i = 0; i < n_point; i++)
		{
			outputfile << x[i] << " " << y[i] << " " << z[i] << std::endl;
		}
		outputfile << "CELLS " << n_cell << " " << 5 * n_cell << std::endl;
		for (size_t i = 0; i < n_cell; i++)
		{
			outputfile << "4 ";
			for (size_t j = 0; j < 4; j++)
			{
				outputfile << indices[i][j] << " ";
			}
			outputfile << std::endl;
		}
		outputfile << "CELL_TYPES " << n_cell << std::endl;
		for (size_t i = 0; i < n_cell; i++)
		{
			outputfile << "10" << std::endl;
		}
		outputfile.close();
	}
}
void Polycube_Algorithm::save_hex_length(const char* filename)
{
	if (final_cube_length < 0.0)
	{
		std::string tmp(filename);
		std::ofstream output_file(tmp + "_error");
		output_file.close();
		return;
	}
	std::ofstream output_file(filename);
	output_file << final_cube_length << std::endl;
	output_file.close();
}
void Polycube_Algorithm::get_max_cut_depth(const std::vector<unsigned int> &cut_edge, std::vector<double> &max_cut_depth, double cut_ratio)
{
	max_cut_depth.clear();
	if (cut_edge.size() == 0)
		return;
	int cut_edge_num = (int)cut_edge.size() / 2;
	for (size_t i = 0; i < cut_edge_num; i++)
	{
		int edge_idx = -1;
		for (size_t j = 0; j < polycube_edge.size(); j++)
		{
			if ((polycube_edge[j].first == cut_edge[2 * i] && polycube_edge[j].second == cut_edge[2 * i + 1]) || (polycube_edge[j].first == cut_edge[2 * i + 1] && polycube_edge[j].second == cut_edge[2 * i]))
			{
				edge_idx = (int)j;
				break;
			}
		}
		assert(edge_idx >= 0);
		double depth = polycube_edge_max_cut_depth_one_cut[edge_idx] * cut_ratio;
		max_cut_depth.push_back(depth);
	}
}
std::pair<int, int> Polycube_Algorithm::compute_flips()
{
	//compute flips of ori mesh and polycube mesh
	int flip_ori = 0;
	int flip_polycube = 0;
	int flip_ori_in = 0;
	int flip_polycube_in = 0;
	const std::vector<Tetrahedron<double>*> &tetras_ori = tet_mesh_->tetras;
	const std::vector<Tetrahedron<double>*> &tetras_polycube = tet_mesh_polycube_->tetras;
	int nc = (int)tetras_ori.size();
	for (size_t i = 0; i < nc; i++)
	{
		if (tetras_ori[i]->compute_tet_volume() < 0)
			flip_ori++;
		if (tetras_polycube[i]->compute_tet_volume() < 0)
			flip_polycube++;
	}
	for (size_t i = 0; i < nc; i++)
	{
		//compute interior flips
		TetVertex<double> **tet_vert_ = tetras_ori[i]->vertex;
		bool boundary_flag = false;
		for (unsigned i = 0; i < 4; ++i)
		{
			if (tet_vert_[i]->boundary)
				boundary_flag = true;
		}
		if (boundary_flag)
			continue;
		if (tetras_ori[i]->compute_tet_volume() < 0)
			flip_ori_in++;
		if (tetras_polycube[i]->compute_tet_volume() < 0)
			flip_polycube_in++;
	}
	std::cout << "origin mesh flip: " << flip_ori << std::endl;
	std::cout << "polycube mesh flip: " << flip_polycube << std::endl;
	std::cout << "origin mesh interior flip: " << flip_ori_in << std::endl;
	std::cout << "polycube mesh interior flip: " << flip_polycube_in << std::endl;
	return (std::pair<int, int>(flip_ori, flip_polycube));
}
void Polycube_Algorithm::get_vert_cell(const TetStructure<double>* tet_mesh, std::vector<CVec<double, 3>> &points, std::vector<unsigned int> &indices)
{
	points.clear();
	indices.clear();
	const std::vector<Tetrahedron<double>*> &tetras = tet_mesh->tetras;
	const std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh->tetra_vertices;
	int nv = (int)tetra_vertices.size();
	int nc = (int)tetras.size();
	for (size_t i = 0; i < nv; i++)
	{
		CVec<double, 3> temp = tetra_vertices[i]->pos;
		points.push_back(temp);
	}
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			unsigned id = (unsigned)tetras[i]->vertex[j]->id;
			indices.push_back(id);
		}
	}
}
void Polycube_Algorithm::get_boundary_vert_face(TetStructure<double>* tet_mesh, std::vector<double> &points, std::vector<int> &faces, std::vector<int> &s2v)
{
	points.clear();
	faces.clear();
	s2v.clear();
	//VolumeSurfaceVertex.clear();
	std::vector<int> tmp_face;
	std::map<int, int> VolumeSurfaceVertex;
	std::map<int, int>::iterator map_it;
	int indexOnSurface = 0;
	int num_tet = (int)tet_mesh->tetras.size();
	std::vector<Tetrahedron<double>*> tetras = tet_mesh->get_tetras();
	static int s_tet_id[4][3] = { { 1, 2, 3}, {2, 0, 3}, {0, 1, 3}, {1, 0, 2} }; //to outer part
	for (size_t i = 0; i < num_tet; i++)
	{
		if (!tetras[i]->is_onboundary())
		{
			continue;
		}
		int face_id = -1;
		for (size_t j = 0; j < 4; j++)
		{
			TetVertex<double>* tv_temp = tetras[i]->vertex[j];
			if (!(tv_temp->boundary))
			{
				face_id = (int)j;
				break;
			}
		}
		if (face_id == -1)
		{
			std::cout << "read vtk error" << std::endl;
			return;
		}
		tmp_face.clear();
		for (size_t j = 0; j < 3; j++)
		{
			map_it = VolumeSurfaceVertex.find((int)tetras[i]->vertex[s_tet_id[face_id][j]]->id);
			if (map_it == VolumeSurfaceVertex.end())
			{
				VolumeSurfaceVertex.insert(std::pair<int, int>(tetras[i]->vertex[s_tet_id[face_id][j]]->id, indexOnSurface));
				s2v.push_back((int)tetras[i]->vertex[s_tet_id[face_id][j]]->id);
				double px, py, pz;
				px = tetras[i]->vertex[s_tet_id[face_id][j]]->pos[0];
				py = tetras[i]->vertex[s_tet_id[face_id][j]]->pos[1];
				pz = tetras[i]->vertex[s_tet_id[face_id][j]]->pos[2];
				SurfaceMesh::Point v2(px, py, pz);
				points.push_back(px);
				points.push_back(py);
				points.push_back(pz);
				tmp_face.push_back(indexOnSurface);
				++indexOnSurface;
			}
			else
			{
				tmp_face.push_back(map_it->second);
			}
		}
		faces.push_back(tmp_face[0]);
		faces.push_back(tmp_face[1]);
		faces.push_back(tmp_face[2]);
	}
}
TetStructure<double>* Polycube_Algorithm::copy_mesh(TetStructure<double>* tet_mesh)
{
	TetStructure<double>* tet_new = new TetStructure<double>;
	std::vector<CVec<double, 3>> points;
	std::vector<unsigned int> indices;
	get_vert_cell(tet_mesh, points, indices);
	tet_new->load_tet(points, indices);
	return tet_new;
}
void Polycube_Algorithm::save_fix_chart_label(const char* filename)
{
	std::ofstream ofs;
	ofs.open(filename);
	int nc = (int)fix_charts.size();
	for (size_t i = 0; i < nc; i++)
	{
		ofs << i << " " << fix_charts[i] << " " << fix_labels[i] << std::endl;
	}
	ofs.close();
}
void Polycube_Algorithm::save_equal_faces(const char* filename)
{
	std::ofstream ofs;
	ofs.open(filename);
	int n_normal_cut, n_three_cut;
	n_normal_cut = (int)cut_types.size();
	ofs << n_normal_cut << std::endl;
	for (size_t i = 0; i < n_normal_cut; i++)
	{
		ofs << cut_types[i] << std::endl;
		ofs << chart_pair_neighbor[i].first << " " << chart_pair[i].first << " " << chart_pair[i].second << " " << chart_pair_neighbor[i].second << std::endl;
		ofs << common_verts[2 * i] << " " << common_verts[2 * i + 1] << std::endl;
		ofs << equal_faces_array[i].size() << std::endl;
		for (size_t j = 0; j < equal_faces_array[i].size(); j++)
		{
			for (size_t k = 0; k < 6; k++)
			{
				ofs << equal_faces_array[i][j][k] << " ";
			}
			ofs << std::endl;
		}
	}
	n_three_cut = (int)three_cut_type_array.size();
	ofs << n_three_cut << std::endl;
	for (size_t i = 0; i < n_three_cut; i++)
	{
		ofs << three_cut_type_array[i] << std::endl;
		ofs << three_cut_common_vert[i] << std::endl;
		ofs << three_cut_vert[i].size() << std::endl;
		for (size_t j = 0; j < three_cut_vert[i].size(); j++)
		{
			ofs << three_cut_vert[i][j][0] << " " << three_cut_vert[i][j][1] << " " << three_cut_vert[i][j][2] << std::endl;
		}
		for (size_t j = 0; j < 3; j++)
		{
			ofs << temp_three_cut_adjacent_one_cut_type[i][j] << " ";
		}
		ofs << std::endl;
		for (size_t j = 0; j < 3; j++)
		{
			ofs << three_cut_adjacent_one_cut_index[i][j] << " ";
		}
		ofs << std::endl;
	}
	ofs.close();
}
void Polycube_Algorithm::save_tet_mesh(const char* filename)
{
	//save tet_mesh_ in tet format
	if (tet_mesh_->tetra_vertices.size() != 0)
	{
		//save tet
		FILE* vtk_file = fopen(filename, "w");
		if (!vtk_file)
		{
			std::cout << "can't write file" << std::endl;
			return;
		}
		const std::vector<Tetrahedron<double>*> &tetras = tet_mesh_->tetras;
		const std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
		int nv = (int)tetra_vertices.size();
		int nc = (int)tetras.size();
		fprintf(vtk_file, "%d vertices\n", (int)tet_mesh_->tetra_vertices.size());
		fprintf(vtk_file, "%d cells\n", nc);
		for (size_t i = 0; i < nv; i++)
		{
			fprintf(vtk_file, "%20.19f %20.19f %20.19f\n", tetra_vertices[i]->pos[0], tetra_vertices[i]->pos[1], tetra_vertices[i]->pos[2]);
		}
		for (size_t i = 0; i < nc; i++)
		{
			fprintf(vtk_file, "4 %d %d %d %d\n", (int)tetras[i]->vertex[0]->id, (int)tetras[i]->vertex[1]->id, (int)tetras[i]->vertex[2]->id, (int)tetras[i]->vertex[3]->id);
		}
		fclose(vtk_file);
	}
}
void Polycube_Algorithm::save_hexex(const char* filename, const TetStructure<double>* oritet, const TetStructure<double>* pqtet)
{
	double cube_length = pq_flattening->get_cube_length();
	std::vector<CVec<double, 3>> pts_pq, pts_ori;
	std::vector<unsigned int> indices_pq, indices_ori;
	get_vert_cell(pqtet, pts_pq, indices_pq);
	get_vert_cell(oritet, pts_ori, indices_ori);
	assert(indices_ori[0] == indices_pq[0]);
	//enlarge pts
	for (size_t i = 0; i < pts_pq.size(); i++)
	{
		CVec<double, 3> q = pts_pq[i] / cube_length;
		for (size_t k = 0; k < 3; k++)
		{
			int tmp = (int)q[k];
			if (abs(tmp - q[k]) < ROUNDING_TH)
				q[k] = (double)tmp;
		}
		pts_pq[i] = q;
	}
	std::ofstream ofs(filename);
	ofs << pts_ori.size() << std::endl;
	for (size_t i = 0; i < pts_ori.size(); i++)
	{
		ofs << pts_ori[i] << std::endl;
	}
	ofs << indices_ori.size() / 4 << std::endl;
	for (size_t i = 0; i < indices_ori.size() / 4; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			ofs << indices_ori[4 * i + j] << ' ';
		}
		for (size_t j = 0; j < 4; j++)
		{
			unsigned int tmpid = indices_ori[4 * i + j];
			for (size_t k = 0; k < 3; k++)
			{
				ofs << std::scientific << pts_pq[tmpid][k] << ' ';
			}
		}
		ofs << std::endl;
	}
	ofs.close();
}
	
void Polycube_Algorithm::save_polycube_tetformat(const char* filename)
{
	//save tet_mesh_ in tet format to feed QuadHexGen
	if (tet_mesh_polycube_->tetra_vertices.size() != 0)
	{
		//save tet
		FILE* vtk_file = fopen(filename, "w");
		if (!vtk_file)
		{
			std::cout << "can't write file" << std::endl;
			return;
		}
		const std::vector<Tetrahedron<double>*> &tetras = tet_mesh_polycube_->tetras;
		const std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_polycube_->tetra_vertices;
		int nv = (int)tetra_vertices.size();
		int nc = (int)tetras.size();
		fprintf(vtk_file, "%d vertices\n", (int)tet_mesh_polycube_->tetra_vertices.size());
		fprintf(vtk_file, "%d cells\n", nc);
		for (size_t i = 0; i < nv; i++)
		{
			fprintf(vtk_file, "%20.19f %20.19f %20.19f\n", tetra_vertices[i]->pos[0], tetra_vertices[i]->pos[1], tetra_vertices[i]->pos[2]);
		}
		for (size_t i = 0; i < nc; i++)
		{
			fprintf(vtk_file, "4 %d %d %d %d\n", (int)tetras[i]->vertex[0]->id, (int)tetras[i]->vertex[1]->id, (int)tetras[i]->vertex[2]->id, (int)tetras[i]->vertex[3]->id);
		}
		fclose(vtk_file);
		std::cout << "output tet file success" << std::endl;
	}
}
void Polycube_Algorithm::save_polycube_para(const char* filename)
{
	if (pq_flattening)
		(*pq_flattening).save_polycube_para(filename, tet_mesh_polycube_);
}
void Polycube_Algorithm::save_flip_info(const char* filename)
{
	double cube_len = final_cube_length;
	if (final_cube_length < 0.0)
		cube_len = 1.0;
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
	int n_equal_vert_pair = 0;
	int n_chart_pair = (int)equal_faces_array.size();
	int nv = (int)tet_mesh_->tetra_vertices.size();
	std::vector<std::vector<std::pair<int, int>>> v_vert_pair;
	v_vert_pair.resize(n_chart_pair);
	std::vector<int> three_cut_vert_type(nv, 0);
	//update three cut vert type
	for (size_t i = 0; i < three_cut_vert.size(); i++)
	{
		for (size_t j = 0; j < three_cut_vert[i].size(); j++)
		{
			unsigned id0, id1, id2;
			id0 = three_cut_vert[i][j][0];
			id1 = three_cut_vert[i][j][1];
			id2 = three_cut_vert[i][j][2];
			three_cut_vert_type[id0] = 3;
			three_cut_vert_type[id1] = 2;
			three_cut_vert_type[id2] = 2;
		}
	}
	for (size_t i = 0; i < n_chart_pair; i++)
	{
		std::vector<int> v_vert_type(nv, 0);
		std::vector<int> v_vert_map(nv, 0);
		int tmp_n = 0;
		for (size_t j = 0; j < equal_faces_array[i].size(); j++)
		{
			assert(equal_faces_array[i][j].size() == 6);
			int id[6];
			id[0] = equal_faces_array[i][j][0];
			id[1] = equal_faces_array[i][j][1];
			id[2] = equal_faces_array[i][j][2];
			id[3] = equal_faces_array[i][j][3];
			id[4] = equal_faces_array[i][j][4];
			id[5] = equal_faces_array[i][j][5];
			for (size_t k = 0; k < 3; k++)
			{
				if (!(three_cut_vert_type[id[k]] == 2 && three_cut_vert_type[id[k + 3]] == 2))
				{
					if (id[k] != id[k + 3])
					{
						v_vert_type[id[k]] = 1;
						v_vert_type[id[k + 3]] = 2;
						v_vert_map[id[k]] = id[k + 3];
					}
					if (three_cut_vert_type[id[k]] == 2 && three_cut_vert_type[id[k + 3]] == 3)
					{
						v_vert_type[id[k]] = 2;
						v_vert_type[id[k + 3]] = 1;
						v_vert_map[id[k + 3]] = id[k];
					}
				}
			}
		}
		for (size_t j = 0; j < nv; j++)
		{
			if (v_vert_type[j] == 1)
			{
				tmp_n++;
				v_vert_pair[i].push_back(std::pair<int, int>(j, v_vert_map[j]));
			}
		}
		n_equal_vert_pair = n_equal_vert_pair + tmp_n;
	}
	std::vector<TetVertex<double>* > tetra_vertices = tet_mesh_polycube_->tetra_vertices;
	std::ofstream ofs(filename);
	ofs << n_equal_vert_pair << std::endl;
	for (size_t i = 0; i < n_chart_pair; i++)
	{
		int one_chart_size = (int)v_vert_pair[i].size();
		Eigen::Matrix3d p;
		Eigen::Vector3d v;
		int type = cut_types[i];
		for (size_t j = 0; j < 3; j++)
		{
			for (size_t k = 0; k < 3; k++)
			{
				p(j, k) = type_matrix[type][j][k];
			}
		}
		int common_vert_idx = common_verts[2 * i];
		v(0) = tetra_vertices[common_vert_idx]->pos[0] / cube_len;
		v(1) = tetra_vertices[common_vert_idx]->pos[1] / cube_len;
		v(2) = tetra_vertices[common_vert_idx]->pos[2] / cube_len;
		Eigen::Vector3d b = v - p * v;
		for (size_t a = 0; a < one_chart_size; a++)
		{
			ofs << v_vert_pair[i][a].first << " " << v_vert_pair[i][a].second << std::endl;
			for (size_t j = 0; j < 3; j++)
			{
				for (size_t k = 0; k < 3; k++)
				{
					//p(j, k) = type_matrix[i][j][k];
					ofs << type_matrix[type][j][k] << " ";
				}
			}
			for (size_t j = 0; j < 3; j++)
			{
				ofs << b(j) << " ";
			}
			ofs << std::endl;
		}
	}
	ofs << nv << std::endl;
	std::vector<OpenVolumeMesh::Geometry::Vec3i> v_type = (*pq_flattening).get_vertex_type();
	for (size_t i = 0; i < nv; i++)
	{
		ig::CVec<int, 3> tmp(0, 0, 0);
		if (v_type[i][0] >= 0)
		{
			tmp[0] = 1;
		}
		if (v_type[i][1] >= 0)
		{
			tmp[1] = 1;
		}
		if (v_type[i][2] >= 0)
		{
			tmp[2] = 1;
		}
		for (size_t j = 0; j < 3; j++)
		{
			ofs << tmp[j] << " ";
		}
		ofs << std::endl;
	}
	ofs.close();
}
void Polycube_Algorithm::data_preparation_ovm()
{
	if (data_prepare_status)
		return;
	if (mesh_ == NULL) return;
	data_prepare_status = true;
	(*pq_deform).prepare_for_deformation(mesh_);
	if (feature_status && !(*pq_deform).feature_edge_flag.empty())
	{
		(*pq_deform).adjust_orientation_feature(mesh_);
		//(*pd_flatten).adjust_orientation(mesh_, feature_flag);
	}
	else
	{
		(*pq_deform).adjust_orientation(mesh_);
	}
	coord_ori.clear();
	for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		int v_id = v_it->idx();
		OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(*v_it);
		coord_ori.push_back(p[0]);
		coord_ori.push_back(p[1]);
		coord_ori.push_back(p[2]);
	}
	//normal set here
	//sigma_s is settled as 1 here
	double sigma_s = 1;
	double sr = 0.3;
	(*pq_deform).build_feature_edges_connectivity(mesh_);
	(*pq_deform).filter_boundary_normal_feature(sigma_s, sr, 3);
	(*pq_deform).target_normal_check_repair();
#ifdef POLYLINE
	(*pd_normal).filter_feature_edge_ring(2, sigma_s, false);
#else
	(*pq_deform).filter_feature_edge_ring(2, sigma_s);
#endif
}
void Polycube_Algorithm::load_vertex_color_vtk(const char* filename)
{
	FILE* f_vc = fopen(filename, "r");
	//return;
	int nv = (int)tet_mesh_->tetra_vertices.size(); char buf[4096];
	int v_id_, c_; int number_of_color = 0;
	std::vector<int> vertex_color(nv, 0);
	char v_id[128]; char c[128]; int v_count = 0;
	while (!feof(f_vc))
	{
		fgets(buf, 4096, f_vc);
		if (v_count < nv)
		{
			sscanf(buf, "%s %s", v_id, c);
			v_id_ = atoi(v_id); c_ = atoi(c);
			vertex_color[v_id_] = c_;
			if (c_ > number_of_color) number_of_color = c_;
			++v_count;
		}
	}
	fclose(f_vc);
	printf("Number Of Colors: %d\n", number_of_color);
	(*pq_deform).set_vertex_color(number_of_color, vertex_color, (int)tet_mesh_->tetra_vertices.size());
	//prepare for deformation
	printf("Prepare Deformation.......\n");
	(*pq_deform).prepare_for_deformation(tet_mesh_);
	//set coord_ori here
	coord_ori.clear();
	std::vector<TetVertex<double>* > tetra_vertices = tet_mesh_->tetra_vertices;
	for (size_t i = 0; i < nv; i++)
	{
		double x, y, z;
		x = tetra_vertices[i]->pos[0];
		y = tetra_vertices[i]->pos[1];
		z = tetra_vertices[i]->pos[2];
		coord_ori.push_back(x);
		coord_ori.push_back(y);
		coord_ori.push_back(z);
	}
}
void Polycube_Algorithm::construct_old_new_vert_map()
{
	if (o2n.size() == tet_mesh_->tetra_vertices.size())
	{
		//already constructed
		return;
	}
	const std::vector<Tetrahedron<double>*> &tetras = tet_mesh_->tetras;
	const std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
	int n_points = (int)tetra_vertices.size();
	int n_cells = (int)tetras.size();
	o2n.clear();
	n2o.clear();
	o2n.resize(n_points, -1);
	std::array<int, 3> tmp_array;
	tmp_array[0] = -1;
	tmp_array[1] = -1;
	tmp_array[2] = -1;
	n2o.resize(n_points, std::array<int, 3>{ {-1, -1, -1}});
	int idx = 0;
	for (size_t i = 0; i < three_cut_vert.size(); i++)
	{
		for (size_t j = 0; j < three_cut_vert[i].size(); j++)
		{
			tmp_array[0] = three_cut_vert[i][j][0];
			tmp_array[1] = three_cut_vert[i][j][1];
			tmp_array[2] = three_cut_vert[i][j][2];
			o2n[tmp_array[0]] = idx;
			o2n[tmp_array[1]] = idx;
			o2n[tmp_array[2]] = idx;
			n2o[idx] = tmp_array;
			idx++;
		}
	}
	//equal_face_part
	tmp_array[2] = -1;
	for (size_t i = 0; i < equal_faces_array.size(); i++)
	{
		for (size_t j = 0; j < equal_faces_array[i].size(); j++)
		{
			for (size_t k = 0; k < 3; k++)
			{
				if (equal_faces_array[i][j][k] != equal_faces_array[i][j][k + 3])
				{
					//not equal
					assert((o2n[equal_faces_array[i][j][k]] != -1) == (o2n[equal_faces_array[i][j][k + 3]] != -1));
					if (o2n[equal_faces_array[i][j][k]] == -1)
					{
						//new vert here
						int id0 = equal_faces_array[i][j][k];
						int id1 = equal_faces_array[i][j][k + 3];
						tmp_array[0] = id0;
						tmp_array[1] = id1;
						o2n[id0] = idx;
						o2n[id1] = idx;
						n2o[idx] = tmp_array;
						idx++;
					}
				}
			}
		}
	}
	//other vert
	tmp_array[1] = -1;
	tmp_array[2] = -1;
	for (size_t i = 0; i < n_points; i++)
	{
		if (o2n[i] == -1)
		{
			tmp_array[0] = (int)i;
			o2n[i] = idx;
			n2o[idx] = tmp_array;
			idx++;
		}
	}
	n2o.resize(idx);
}
bool Polycube_Algorithm::auto_color()
{
	if (tet_mesh_)
	{
		const std::vector<Tetrahedron<double>*> &tetras = tet_mesh_->tetras;
		const std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
		int n_points = (int)tetra_vertices.size();
		int n_cells = (int)tetras.size();
		if (equal_faces_array.size() == 0)
		{
			static int s_tet_id[4][3] = { { 2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2} };
			std::vector<std::vector<int>> vert_edges;
			std::vector<int>::iterator iter;
			vert_edges.resize(n_points);
			for (size_t i = 0; i < n_cells; i++)
			{
				for (size_t j = 0; j < 4; j++)
				{
					//four verts
					//int center = cells[i][j];
					int center = (int)tetras[i]->vertex[j]->id;
					int idx[3];
					idx[0] = (int)tetras[i]->vertex[s_tet_id[j][0]]->id;
					idx[1] = (int)tetras[i]->vertex[s_tet_id[j][1]]->id;
					idx[2] = (int)tetras[i]->vertex[s_tet_id[j][2]]->id;
					//only add verts with bigger idx
					for (size_t k = 0; k < 3; k++)
					{
						if (idx[k] < center)
							continue;
						iter = std::find(vert_edges[center].begin(), vert_edges[center].end(), idx[k]);
						if (iter == vert_edges[center].end())
						{
							//no find
							vert_edges[center].push_back(idx[k]);
						}
					}
				}
			}
			std::vector<std::pair<int, int>> edges;
			std::pair<int, int> temp_pair;
			//collect edges
			for (size_t i = 0; i < n_points; i++)
			{
				temp_pair.first = (int)i;
				for (size_t j = 0; j < vert_edges[i].size(); j++)
				{
					temp_pair.second = vert_edges[i][j];
					edges.push_back(temp_pair);
				}
			}
			std::vector<std::vector<int>> v_same_color;
			std::vector<int> v_color;
			ig::CGALHelper::graph_coloring(n_points, edges, v_same_color, v_color);
			(*pq_deform).set_vertex_color((int)v_same_color.size(), v_color, (int)tet_mesh_->tetra_vertices.size());
		}
		else {
			construct_old_new_vert_map();
			int idx = (int)n2o.size();
			//map construction finished
			static int s_tet_id[4][3] = { { 2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2} };
			std::vector<std::vector<int>> vert_edges;
			std::vector<int>::iterator iter;
			vert_edges.resize(idx);
			for (size_t i = 0; i < n_cells; i++)
			{
				for (size_t j = 0; j < 4; j++)
				{
					//four verts
					int center = o2n[tetras[i]->vertex[j]->id];
					int tmp_idx[3];
					tmp_idx[0] = o2n[tetras[i]->vertex[s_tet_id[j][0]]->id];
					tmp_idx[1] = o2n[tetras[i]->vertex[s_tet_id[j][1]]->id];
					tmp_idx[2] = o2n[tetras[i]->vertex[s_tet_id[j][2]]->id];
					//only add verts with bigger idx
					for (size_t k = 0; k < 3; k++)
					{
						if (tmp_idx[k] < center)
							continue;
						iter = std::find(vert_edges[center].begin(), vert_edges[center].end(), tmp_idx[k]);
						if (iter == vert_edges[center].end())
						{
							//no find
							vert_edges[center].push_back(tmp_idx[k]);
						}
					}
				}
			}
			std::vector<std::pair<int, int>> edges;
			for (size_t i = 0; i < idx; i++)
			{
				std::pair<int, int> temp_pair;
				temp_pair.first = (int)i;
				for (size_t j = 0; j < vert_edges[i].size(); j++)
				{
					temp_pair.second = vert_edges[i][j];
					edges.push_back(temp_pair);
				}
			}
			std::vector<std::vector<int>> v_same_color;
			std::vector<int> v_color;
			std::vector<int> v_color_ori;
			v_color_ori.resize(n_points);
			ig::CGALHelper::graph_coloring(n_points, edges, v_same_color, v_color);
			//construct v_color_ori
			for (size_t i = 0; i < idx; i++)
			{
				int old_id0, old_id1, old_id2;
				old_id0 = n2o[i][0];
				old_id1 = n2o[i][1];
				old_id2 = n2o[i][2];
				assert(old_id0 != -1);
				v_color_ori[old_id0] = v_color[i];
				if (old_id1 != -1)
				{
					v_color_ori[old_id1] = v_color[i];
				}
				if (old_id2 != -1)
				{
					v_color_ori[old_id2] = v_color[i];
				}
			}
			(*pq_deform).set_vertex_color((int)v_same_color.size(), v_color_ori, (int)tet_mesh_->tetra_vertices.size());
		}
		return true;
	}
	//for ovm file
	if (mesh_)
	{
		//Only consider the case with no equal faces
		int n_points = mesh_->n_vertices();
		std::vector<std::pair<int, int>> edges;
		for (auto it = mesh_->edges_begin(); it != mesh_->edges_end(); ++it)
		{
			//two endpoints of edge
			OpenVolumeMesh::HalfEdgeHandle hfh = mesh_->halfedge_handle(*it, 0);
			OpenVolumeMesh::OpenVolumeMeshEdge ovme = mesh_->edge(*it);
			std::pair<int, int> one_edge;
			one_edge.first = ovme.from_vertex().idx();;
			one_edge.second = ovme.to_vertex().idx();
			edges.push_back(one_edge);
		}
		std::vector<std::vector<int>> v_same_color;
		std::vector<int> v_color;
		ig::CGALHelper::graph_coloring(n_points, edges, v_same_color, v_color);
		(*pq_deform).set_vertex_color((int)v_same_color.size(), v_color, (int)mesh_->n_vertices());
		return true;
	}
	return true;
}
bool Polycube_Algorithm::load_chart_label(const char* filename)
{
	fix_charts.clear();
	fix_labels.clear();
	std::ifstream ifs;
	ifs.open(filename);
	if (ifs.fail())
		return false;
	int nc;
	if (mesh_)
	{
		data_preparation_ovm();
		nc = mesh_->n_cells();
	}
	else
		nc = (int)tet_mesh_->tetras.size();
	int max_chart = -1;
	for (size_t i = 0; i < nc; i++)
	{
		//end of file reached
		if (ifs.eof())
			return false;
		int temp_idx, temp_chart, temp_label;
		ifs >> temp_idx >> temp_chart >> temp_label;
		fix_charts.push_back(temp_chart);
		fix_labels.push_back(temp_label);
		max_chart = max_chart < temp_chart ? temp_chart : max_chart;
	}
	if (mesh_)
		(*pq_deform).set_fixed_chart_label(fix_charts, fix_labels);
	else
		(*pq_deform).set_fixed_chart_label(fix_charts, fix_labels, false);
	ifs.close();
	return true;
}
//only for one cut version
void Polycube_Algorithm::load_equal_faces(const char* filename)
{
	(*pq_deform).prepare_for_deformation(tet_mesh_polycube_);
	coord_ori.clear();
	equal_faces_array.clear();
	cut_types.clear();
	common_verts.clear();
	chart_pair.clear();
	chart_pair_neighbor.clear();
	three_cut_type_array.clear();
	three_cut_common_vert.clear();
	temp_three_cut_adjacent_one_cut_type.clear();
	three_cut_adjacent_one_cut_index.clear();
	three_cut_vert.clear();
	//equal_triangle_size.clear();
	std::vector<std::vector<int>> equal_faces;
	std::ifstream ifs;
	std::string line;
	ifs.open(filename);
	int count = 0;
	int n_normal_cut, n_three_cut;
	n_normal_cut = 0;
	n_three_cut = 0;
	ifs >> n_normal_cut;
	for (size_t iter = 0; iter < n_normal_cut; iter++)
	{
		int type_id = -1;
		ifs >> type_id;
		equal_faces.clear();
		cut_types.push_back(type_id);
		std::pair<int, int> temp_pair, temp_pair_neighbor;
		ifs >> temp_pair_neighbor.first >> temp_pair.first >> temp_pair.second >> temp_pair_neighbor.second;
		chart_pair.push_back(temp_pair);
		chart_pair_neighbor.push_back(temp_pair_neighbor);
		int temp_common_vert;
		ifs >> temp_common_vert;
		common_verts.push_back(temp_common_vert);
		ifs >> temp_common_vert;
		common_verts.push_back(temp_common_vert);
		int n_lines;
		//equal_triangle_size.push_back(n_lines);
		ifs >> n_lines;
		for (size_t i = 0; i < n_lines; i++)
		{
			//temp_face.clear();
			int id1, id2, id3, id4, id5, id6;
			ifs >> id1 >> id2 >> id3 >> id4 >> id5 >> id6;
			std::vector<int> temp_face{ id1, id2, id3, id4, id5, id6 };
			equal_faces.push_back(temp_face);
		}
		equal_faces_array.push_back(equal_faces);
		type_id = -1;
		
	}
	int n_vert;
	if (mesh_)
		n_vert = mesh_->n_vertices();
	else
		n_vert = (int)tet_mesh_->tetra_vertices.size();
	ifs >> n_three_cut;
	for (size_t i = 0; i < n_three_cut; i++)
	{
		int cut_type, common_vert;
		ifs >> cut_type;
		ifs >> common_vert;
		three_cut_type_array.push_back(cut_type);
		three_cut_common_vert.push_back(common_vert);
		int n_line;
		ifs >> n_line;
		std::vector<std::array<unsigned int, 3>> one_three_cut;
		std::array<unsigned, 3> temp_array;
		for (size_t j = 0; j < n_line; j++)
		{
			ifs >> temp_array[0] >> temp_array[1] >> temp_array[2];
			one_three_cut.push_back(temp_array);
		}
		three_cut_vert.push_back(one_three_cut);
		std::array<int, 3> temp_array_int;
		ifs >> temp_array_int[0] >> temp_array_int[1] >> temp_array_int[2];
		temp_three_cut_adjacent_one_cut_type.push_back(temp_array_int);
		ifs >> temp_array_int[0] >> temp_array_int[1] >> temp_array_int[2];
		three_cut_adjacent_one_cut_index.push_back(temp_array_int);
	}
	(*pq_deform).set_equal_faces(equal_faces_array, common_verts, cut_types, chart_pair, three_cut_type_array, three_cut_common_vert, three_cut_vert, temp_three_cut_adjacent_one_cut_type);
	(*pq_flattening).set_equal_faces(equal_faces_array, common_verts, cut_types, n_vert, chart_pair, chart_pair_neighbor, three_cut_common_vert, three_cut_vert, three_cut_adjacent_one_cut_index);
	ifs.close();
}
void Polycube_Algorithm::load_equal_faces(const std::vector<int> &v_cut_type, const std::vector<int> &v_cut_section_chart, const std::vector<int> &v_seam_edge_vertex_id, const std::vector<std::vector<int>>& vv_congruent_face, int nv, std::vector<ThreeCut> &three_cut_info)
{
	equal_faces_array.clear();
	cut_types.clear();
	common_verts.clear();
	chart_pair.clear();
	chart_pair_neighbor.clear();
	three_cut_type_array.clear();
	three_cut_common_vert.clear();
	temp_three_cut_adjacent_one_cut_type.clear();
	three_cut_adjacent_one_cut_index.clear();
	three_cut_vert.clear();
	cut_types = v_cut_type;
	common_verts = v_seam_edge_vertex_id;
	std::vector<std::vector<int>> equal_faces;
	for (size_t i = 0; i < vv_congruent_face.size(); i++)
	{
		int n_faces = (int)vv_congruent_face[i].size() / 6;
		equal_faces.clear();
		for (size_t j = 0; j < n_faces; j++)
		{
			std::vector<int> temp{ vv_congruent_face[i][6 * j], vv_congruent_face[i][6 * j + 1], vv_congruent_face[i][6 * j + 2], vv_congruent_face[i][6 * j + 3], vv_congruent_face[i][6 * j + 4], vv_congruent_face[i][6 * j + 5] };
			equal_faces.push_back(temp);
		}
		equal_faces_array.push_back(equal_faces);
	}
	//deform part need to be loaded
	int chart_pair_size = (int)v_cut_section_chart.size() / 4;
	for (size_t i = 0; i < chart_pair_size; i++)
	{
		chart_pair.push_back(std::pair<int, int>(v_cut_section_chart[4 * i + 1], v_cut_section_chart[4 * i + 2]));
		chart_pair_neighbor.push_back(std::pair<int, int>(v_cut_section_chart[4 * i], v_cut_section_chart[4 * i + 3]));
	}
	//three cut information
	three_cut_vert.clear();
	temp_three_cut_adjacent_one_cut_type.resize(three_cut_info.size());
	for (size_t i = 0; i < three_cut_info.size(); i++)
	{
		three_cut_type_array.push_back(three_cut_info[i].cut_type);
		three_cut_common_vert.push_back(three_cut_info[i].shared_vertex_id);
		three_cut_vert.push_back(three_cut_info[i].va_corresponding_vertex);
		three_cut_adjacent_one_cut_index.push_back(three_cut_info[i].a_adjacent_one_cut_index);
		for (int j = 0; j < 3; ++j)
			temp_three_cut_adjacent_one_cut_type[i][j] = v_cut_type[three_cut_info[i].a_adjacent_one_cut_index[j]];
	}
	(*pq_deform).set_equal_faces(equal_faces_array, v_seam_edge_vertex_id, v_cut_type, chart_pair, three_cut_type_array, three_cut_common_vert, three_cut_vert, temp_three_cut_adjacent_one_cut_type);
	//to be update
	(*pq_flattening).set_equal_faces(equal_faces_array, v_seam_edge_vertex_id, v_cut_type, nv, chart_pair, chart_pair_neighbor, three_cut_common_vert, three_cut_vert, three_cut_adjacent_one_cut_index);
}
void Polycube_Algorithm::feature_aware_deformation(int max_iter, double sigma_s, double threshold)
{
#ifdef POLYLINE
	deformation_polylines(max_iter, sigma_s);
	return;
#endif
	data_preparation_ovm();
	//bool feature_flag = true;
	bool boundary_global_opt = false; //should always be false
	(*pq_deform).set_sigma_s(sigma_s);
	for (int j = 0; j < max_iter; ++j)
	{
		(*pq_deform).exp_mips_deformation_refine_polycube_omp_feature_nf(1, 1, 1.0, 0.5, mesh_, false, true, false, feature_status, boundary_global_opt);
		(*pq_deform).assign_pos_mesh(mesh_);
		if (pq_deform->check_polycube_planarity(threshold))
		{
			break;
		}
		if (feature_status && !pq_deform->feature_edge_flag.empty())
		{
			(*pq_deform).adjust_orientation(mesh_);
		}
		else
		{
			(*pq_deform).adjust_orientation(mesh_, feature_status);
		}
	}
	(*pq_deform).compute_distortion(mesh_);
	(*pq_deform).update_face_center();
	return;
}
void Polycube_Algorithm::deformation_polylines(int max_iter, double sigma_s, bool assign_current, bool ns_flag, double smooth_factor, double ratio_threshold)
{
	data_preparation_ovm();
	(*pq_deform).set_sigma_s(sigma_s);
	for (int j = 0; j < max_iter; ++j)
	{
		(*pq_deform).deformation_polylines(mesh_, assign_current, ns_flag, smooth_factor);
		(*pq_deform).assign_pos_mesh(mesh_);
		//angle check
		if ((*pq_deform).check_polyline_orthogonality(ratio_threshold))
			break;
		if (feature_status && !(*pq_deform).feature_edge_flag.empty())
		{
			(*pq_deform).adjust_orientation_feature(mesh_, !assign_current); //assign_current and polyline_flag is df
		}
	}
	//refine corners and long edge
	(*pq_deform).corner_redefinition_polyline();
	(*pq_deform).feature_long_edge_redefinition();
	(*pq_deform).compute_distortion(mesh_);
	(*pq_deform).update_face_center();
	return;
}
void Polycube_Algorithm::save_polyline_vtkformat(const char* filename)
{
	if (pq_deform->dpx_feature.empty())
	{
		std::cout << "No PolyLine" << std::endl;
		return;
	}
	std::vector<int> used_vert_map(pq_deform->dpx_feature.size(), -1);
	std::vector<int> reverse_vert_map(pq_deform->dpx_feature.size(), -1);
	size_t vert_count = 0;
	const std::vector<std::pair<int, int>> &feature_e2v = pq_deform->feature_e2v;
	size_t edge_count = 0;
	for (size_t i = 0; i < feature_e2v.size(); i++)
	{
		if (feature_e2v[i].first == -1) continue;
		edge_count++;
		int id0 = feature_e2v[i].first, id1 = feature_e2v[i].second;
		if (used_vert_map[id0] == -1)
		{
			used_vert_map[id0] = (int)vert_count;
			reverse_vert_map[vert_count] = id0;
			vert_count++;
		}
		if (used_vert_map[id1] == -1)
		{
			used_vert_map[id1] = (int)vert_count;
			reverse_vert_map[vert_count] = id1;
			vert_count++;
		}
	}
	std::ofstream outputfile(filename);
	outputfile << "# vtk DataFile Version 3.0\n"
		<< "mesh vtk data\n"
		<< "ASCII\n"
		<< "DATASET POLYDATA\n";
	//int n_point = dpx.size();
	//test code below
	vert_count = pq_deform->dpx_feature.size();
	for (size_t i = 0; i < vert_count; i++)
	{
		reverse_vert_map[i] = (int)i;
		used_vert_map[i] = (int)i;
	}
	outputfile << "POINTS " << vert_count << " double" << std::endl;
	for (size_t i = 0; i < vert_count; i++)
	{
		int oid = reverse_vert_map[i];
		outputfile << pq_deform->dpx_feature[oid] << " " << pq_deform->dpy_feature[oid] << " " << pq_deform->dpz_feature[oid] << std::endl;
	}
	outputfile << "LINES " << edge_count << " " << 3 * edge_count << std::endl;
	for (size_t i = 0; i < feature_e2v.size(); i++)
	{
		if (feature_e2v[i].first == -1) continue;
		//edge_count++;
		int id0 = feature_e2v[i].first, id1 = feature_e2v[i].second;
		outputfile << "2 " << used_vert_map[id0] << " " << used_vert_map[id1] << std::endl;
	}
	outputfile.close();
}
void Polycube_Algorithm::do_IVF_seamless_moving(int smooth_iter, int max_iter)
{
	//optimize original mesh based on polycube with seamless constraints
	if (!tet_mesh_polycube_)
		return;
	construct_old_new_vert_map();
	std::vector<double> nodes_polycube;
	std::vector<double> nodes_ori;
	std::vector<int> ref_nodes;
	std::vector<unsigned int> elems, elems_new;
	std::vector<int> node_normal_tag;
	std::vector<double> output_coord;
	std::vector<CVec<double, 3>> pts_polycube, pts_ori, pts_polycube_new, pts_ori_new;
	int nv = (int)tet_mesh_polycube_->tetra_vertices.size();
	int nc = (int)tet_mesh_polycube_->tetras.size();
	get_vert_cell(tet_mesh_polycube_, pts_polycube, elems);
	get_vert_cell(tet_mesh_, pts_ori, elems);
	assert(pts_polycube.size() == pts_ori.size());
	TetStructure<double> tet_mesh_new;
	//construct tet_mesh_new
	int nv_new = (int)n2o.size();
	pts_polycube_new.resize(nv_new);
	pts_ori_new.resize(nv_new);
	elems_new.resize(elems.size());
	for (int i = 0; i < nv; i++)
	{
		pts_polycube_new[o2n[i]] = pts_polycube[i];
		pts_ori_new[o2n[i]] = pts_ori[i];
	}
	for (int i = 0; i < 4 * nc; i++)
	{
		elems_new[i] = o2n[elems[i]];
	}
	tet_mesh_new.load_tet(pts_ori_new, elems_new);
	nodes_polycube.resize(3 * nv_new);
	nodes_ori.resize(3 * nv_new);
	ref_nodes.resize(nv_new, 0);
	node_normal_tag.resize(3 * nv_new, 0);
	for (size_t i = 0; i < nv_new; i++)
	{
		nodes_polycube[3 * i + 0] = pts_polycube_new[i][0];
		nodes_polycube[3 * i + 1] = pts_polycube_new[i][1];
		nodes_polycube[3 * i + 2] = pts_polycube_new[i][2];
		nodes_ori[3 * i + 0] = pts_ori_new[i][0];
		nodes_ori[3 * i + 1] = pts_ori_new[i][1];
		nodes_ori[3 * i + 2] = pts_ori_new[i][2];
		if (tet_mesh_new.tetra_vertices[i]->boundary)
		{
			//ref_nodes[i] = 1;
			node_normal_tag[3 * i] = 1;
		}
	}
	std::vector<int> ori_faces, boundary_faces;
	std::vector<double> ori_pts, boundary_pts;
	std::vector<int> s2v, s2v_ori;
	get_boundary_vert_face(&tet_mesh_new, boundary_pts, boundary_faces, s2v);
	get_boundary_vert_face(tet_mesh_unchanged, ori_pts, ori_faces, s2v_ori);
	SUS_proj_laplacian(nodes_ori, ref_nodes, elems_new, node_normal_tag, nodes_polycube, ori_pts, ori_faces, output_coord, s2v, boundary_faces, smooth_iter, max_iter);
	std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
	for (size_t i = 0; i < nv_new; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			int idx = n2o[i][j];
			if (idx >= 0)
			{
				//update pos
				tetra_vertices[idx]->pos[0] = output_coord[3 * i + 0];
				tetra_vertices[idx]->pos[1] = output_coord[3 * i + 1];
				tetra_vertices[idx]->pos[2] = output_coord[3 * i + 2];
			}
		}
	}
}
void Polycube_Algorithm::do_IVF_seamless_moving_part(int smooth_iter, int max_iter, int neighbor_size)
{
	//optimize original mesh based on polycube with seamless constraints
	if (!tet_mesh_polycube_)
		return;
	construct_old_new_vert_map();
	//here _new means that the seamless points are merged
	std::vector<double> nodes_polycube;
	std::vector<double> nodes_ori;
	std::vector<int> ref_nodes;
	std::vector<unsigned int> elems, elems_new;
	std::vector<int> node_normal_tag;
	std::vector<double> output_coord;
	std::vector<CVec<double, 3>> pts_polycube, pts_ori, pts_polycube_new, pts_ori_new;
	int nv = (int)tet_mesh_polycube_->tetra_vertices.size();
	int nc = (int)tet_mesh_polycube_->tetras.size();
	get_vert_cell(tet_mesh_polycube_, pts_polycube, elems);
	get_vert_cell(tet_mesh_, pts_ori, elems);
	assert(pts_polycube.size() == pts_ori.size());
	TetStructure<double> tet_mesh_new;
	//construct tet_mesh_new
	int nv_new = (int)n2o.size();
	pts_polycube_new.resize(nv_new);
	pts_ori_new.resize(nv_new);
	elems_new.resize(elems.size());
	for (int i = 0; i < nv; i++)
	{
		pts_polycube_new[o2n[i]] = pts_polycube[i];
		pts_ori_new[o2n[i]] = pts_ori[i];
	}
	for (int i = 0; i < 4 * nc; i++)
	{
		elems_new[i] = o2n[elems[i]];
	}
	tet_mesh_new.load_tet(pts_ori_new, elems_new);
	nodes_polycube.resize(3 * nv_new);
	nodes_ori.resize(3 * nv_new);
	ref_nodes.resize(nv_new, 1);
	node_normal_tag.resize(3 * nv_new, 0);
	for (size_t i = 0; i < nv_new; i++)
	{
		nodes_polycube[3 * i + 0] = pts_polycube_new[i][0];
		nodes_polycube[3 * i + 1] = pts_polycube_new[i][1];
		nodes_polycube[3 * i + 2] = pts_polycube_new[i][2];
		nodes_ori[3 * i + 0] = pts_ori_new[i][0];
		nodes_ori[3 * i + 1] = pts_ori_new[i][1];
		nodes_ori[3 * i + 2] = pts_ori_new[i][2];
		if (tet_mesh_new.tetra_vertices[i]->boundary)
		{
			//ref_nodes[i] = 1;
			node_normal_tag[3 * i] = 1;
		}
	}
	const std::vector<Tetrahedron<double>*> &tetras_new = tet_mesh_new.tetras;
	const std::vector<TetVertex<double>*> &tetra_vertices_new = tet_mesh_new.tetra_vertices;
	std::vector<std::vector<int>> vert_edges;
	vert_edges.resize(nv_new);
	std::vector<int> flip_vert_idx, new_flip_vert_idx;
	for (size_t i = 0; i < nc; i++)
	{
		if (tetras_new[i]->compute_tet_volume() < 0.0)
		{
			for (size_t j = 0; j < 4; j++)
			{
				flip_vert_idx.push_back((int)tetras_new[i]->vertex[j]->id);
			}
		}
	}
	static int s_tet_id[4][3] = { { 2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2} };
	std::vector<int>::iterator iter;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			//four verts
			int center = (int)tetras_new[i]->vertex[j]->id;
			int idx[3];
			idx[0] = (int)tetras_new[i]->vertex[s_tet_id[j][0]]->id;
			idx[1] = (int)tetras_new[i]->vertex[s_tet_id[j][1]]->id;
			idx[2] = (int)tetras_new[i]->vertex[s_tet_id[j][2]]->id;
			for (size_t k = 0; k < 3; k++)
			{
				iter = std::find(vert_edges[center].begin(), vert_edges[center].end(), idx[k]);
				if (iter == vert_edges[center].end())
				{
					//no find
					vert_edges[center].push_back(idx[k]);
				}
			}
		}
	}
	for (size_t i = 0; i < flip_vert_idx.size(); i++)
	{
		ref_nodes[flip_vert_idx[i]] = 0;
	}
	for (size_t iter = 0; iter < neighbor_size; iter++)
	{
		new_flip_vert_idx.clear();
		for (size_t i = 0; i < flip_vert_idx.size(); i++)
		{
			int center = flip_vert_idx[i];
			for (size_t j = 0; j < vert_edges[center].size(); j++)
			{
				new_flip_vert_idx.push_back(vert_edges[center][j]);
			}
		}
		flip_vert_idx = new_flip_vert_idx;
		for (size_t i = 0; i < flip_vert_idx.size(); i++)
		{
			ref_nodes[flip_vert_idx[i]] = 0;
		}
	}
	std::vector<int> ori_faces, boundary_faces;
	std::vector<double> ori_pts, boundary_pts;
	std::vector<int> s2v, s2v_ori;
	get_boundary_vert_face(&tet_mesh_new, boundary_pts, boundary_faces, s2v);
	get_boundary_vert_face(tet_mesh_unchanged, ori_pts, ori_faces, s2v_ori);
	std::vector<int> elems_new_int;
	for (size_t i = 0; i < elems_new.size(); i++)
	{
		elems_new_int.push_back(elems_new[i]);
	}
	SUSSpace::SUS_standard(nodes_ori, ref_nodes, elems_new_int, output_coord, smooth_iter, max_iter);
	std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
	for (size_t i = 0; i < nv_new; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			int idx = n2o[i][j];
			if (idx >= 0)
			{
				//update pos
				tetra_vertices[idx]->pos[0] = output_coord[3 * i + 0];
				tetra_vertices[idx]->pos[1] = output_coord[3 * i + 1];
				tetra_vertices[idx]->pos[2] = output_coord[3 * i + 2];
			}
		}
	}
}
void Polycube_Algorithm::do_IVF_seamless_svd(int max_iter)
{
	//optimize original mesh based on polycube with seamless constraints
	if (!tet_mesh_polycube_)
		return;
	bool fix_boundary = false;
	bool using_standard_tet = true;
	double min_singular_value = 0.001;
	construct_old_new_vert_map();
	//here _new means that the seamless points are merged
	std::vector<double> nodes_polycube;
	std::vector<double> nodes_ori;
	std::vector<int> ref_nodes;
	std::vector<unsigned int> elems, elems_new;
	std::vector<int> node_normal_tag;
	std::vector<double> output_coord;
	std::vector<CVec<double, 3>> pts_polycube, pts_ori, pts_polycube_new, pts_ori_new;
	int nv = (int)tet_mesh_polycube_->tetra_vertices.size();
	int nc = (int)tet_mesh_polycube_->tetras.size();
	get_vert_cell(tet_mesh_polycube_, pts_polycube, elems);
	get_vert_cell(tet_mesh_, pts_ori, elems);
	assert(pts_polycube.size() == pts_ori.size());
	TetStructure<double> tet_mesh_new;
	//construct tet_mesh_new
	int nv_new = (int)n2o.size();
	pts_polycube_new.resize(nv_new);
	pts_ori_new.resize(nv_new);
	elems_new.resize(elems.size());
	for (int i = 0; i < nv; i++)
	{
		pts_polycube_new[o2n[i]] = pts_polycube[i];
		pts_ori_new[o2n[i]] = pts_ori[i];
	}
	for (int i = 0; i < 4 * nc; i++)
	{
		elems_new[i] = o2n[elems[i]];
	}
	//always fix boundary 
	tet_mesh_new.load_tet(pts_ori_new, elems_new);
	nodes_polycube.resize(3 * nv_new);
	nodes_ori.resize(3 * nv_new);
	ref_nodes.resize(nv_new, 0);
	ref_nodes[0] = 1;
	node_normal_tag.resize(3 * nv_new, 0);
	if (fix_boundary)
	{
		ref_nodes[0] = 0;
		for (size_t i = 0; i < nv_new; i++)
		{
			if (tet_mesh_new.tetra_vertices[i]->boundary)
			{
				ref_nodes[i] = 1;
			}
		}
	}
	ig::IVF::IVF_svd(pts_polycube_new, elems_new, pts_ori_new, ref_nodes, using_standard_tet, min_singular_value, max_iter);
	std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
	for (size_t i = 0; i < nv_new; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			int idx = n2o[i][j];
			if (idx >= 0)
			{
				//update pos
				tetra_vertices[idx]->pos = pts_ori_new[i];
			}
		}
	}
	compute_flips();
}
void Polycube_Algorithm::do_IVF_seamless_svd_part(int max_iter, int neighbor_size)
{
	//optimize original mesh based on polycube with seamless constraints
	if (!tet_mesh_polycube_)
		return;
	bool fix_boundary = false;
	bool using_standard_tet = true;
	double min_singular_value = 0.001;
	construct_old_new_vert_map();
	//here _new means that the seamless points are merged
	std::vector<double> nodes_polycube;
	std::vector<double> nodes_ori;
	std::vector<int> ref_nodes;
	std::vector<unsigned int> elems, elems_new;
	std::vector<int> node_normal_tag;
	std::vector<double> output_coord;
	std::vector<CVec<double, 3>> pts_polycube, pts_ori, pts_polycube_new, pts_ori_new;
	int nv = (int)tet_mesh_polycube_->tetra_vertices.size();
	int nc = (int)tet_mesh_polycube_->tetras.size();
	get_vert_cell(tet_mesh_polycube_, pts_polycube, elems);
	get_vert_cell(tet_mesh_, pts_ori, elems);
	assert(pts_polycube.size() == pts_ori.size());
	TetStructure<double> tet_mesh_new;
	//construct tet_mesh_new
	int nv_new = (int)n2o.size();
	pts_polycube_new.resize(nv_new);
	pts_ori_new.resize(nv_new);
	elems_new.resize(elems.size());
	for (int i = 0; i < nv; i++)
	{
		pts_polycube_new[o2n[i]] = pts_polycube[i];
		pts_ori_new[o2n[i]] = pts_ori[i];
	}
	for (int i = 0; i < 4 * nc; i++)
	{
		elems_new[i] = o2n[elems[i]];
	}
	//always fix boundary 
	tet_mesh_new.load_tet(pts_ori_new, elems_new);
	nodes_polycube.resize(3 * nv_new);
	nodes_ori.resize(3 * nv_new);
	ref_nodes.resize(nv_new, 1);
	//ref_nodes[0] = 1;
	node_normal_tag.resize(3 * nv_new, 0);
	if (fix_boundary)
	{
		ref_nodes[0] = 0;
		for (size_t i = 0; i < nv_new; i++)
		{
			if (tet_mesh_new.tetra_vertices[i]->boundary)
			{
				ref_nodes[i] = 1;
			}
		}
	}
	//construct ref_nodes
	const std::vector<Tetrahedron<double>*> &tetras_new = tet_mesh_new.tetras;
	const std::vector<TetVertex<double>*> &tetra_vertices_new = tet_mesh_new.tetra_vertices;
	std::vector<std::vector<int>> vert_edges;
	vert_edges.resize(nv_new);
	std::vector<int> flip_vert_idx, new_flip_vert_idx;
	for (size_t i = 0; i < nc; i++)
	{
		if (tetras_new[i]->compute_tet_volume() < 0.0)
		{
			for (size_t j = 0; j < 4; j++)
			{
				flip_vert_idx.push_back((int)tetras_new[i]->vertex[j]->id);
			}
		}
	}
	static int s_tet_id[4][3] = { { 2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2} };
	std::vector<int>::iterator iter;
	for (size_t i = 0; i < nc; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			//four verts
			int center = (int)tetras_new[i]->vertex[j]->id;
			int idx[3];
			idx[0] = (int)tetras_new[i]->vertex[s_tet_id[j][0]]->id;
			idx[1] = (int)tetras_new[i]->vertex[s_tet_id[j][1]]->id;
			idx[2] = (int)tetras_new[i]->vertex[s_tet_id[j][2]]->id;
			for (size_t k = 0; k < 3; k++)
			{
				iter = std::find(vert_edges[center].begin(), vert_edges[center].end(), idx[k]);
				if (iter == vert_edges[center].end())
				{
					//no find
					vert_edges[center].push_back(idx[k]);
				}
			}
		}
	}
	for (size_t i = 0; i < flip_vert_idx.size(); i++)
	{
		ref_nodes[flip_vert_idx[i]] = 0;
	}
	for (size_t iter = 0; iter < neighbor_size; iter++)
	{
		new_flip_vert_idx.clear();
		for (size_t i = 0; i < flip_vert_idx.size(); i++)
		{
			int center = flip_vert_idx[i];
			for (size_t j = 0; j < vert_edges[center].size(); j++)
			{
				new_flip_vert_idx.push_back(vert_edges[center][j]);
			}
		}
		flip_vert_idx = new_flip_vert_idx;
		for (size_t i = 0; i < flip_vert_idx.size(); i++)
		{
			ref_nodes[flip_vert_idx[i]] = 0;
		}
	}
	ig::IVF::IVF_svd(pts_polycube_new, elems_new, pts_ori_new, ref_nodes, using_standard_tet, min_singular_value, max_iter);
	std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
	for (size_t i = 0; i < nv_new; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			int idx = n2o[i][j];
			if (idx >= 0)
			{
				//update pos
				tetra_vertices[idx]->pos = pts_ori_new[i];
			}
		}
	}
	compute_flips();
}
void Polycube_Algorithm::do_IVF_standard(TetStructure<double>* tet_, bool fix_boundary, int smooth_iter, int max_iter)
{
	if (!tet_)
		return;
	//optimize origin mesh based on polycube mesh
	std::vector<double> nodes_polycube;
	std::vector<int> ref_nodes;
	std::vector<unsigned int> elems;
	std::vector<double> output_coord;
	std::vector<CVec<double, 3>> pts_polycube, pts_ori;
	int nv = (int)tet_->tetra_vertices.size();
	int nc = (int)tet_->tetras.size();
	nodes_polycube.resize(3 * nv);
	ref_nodes.resize(nv, 0);
	get_vert_cell(tet_, pts_polycube, elems);
	const std::vector<Tetrahedron<double>*> &tetras = tet_->tetras;
	const std::vector<TetVertex<double>*> &tetra_vertices = tet_->tetra_vertices;
	for (size_t i = 0; i < nv; i++)
	{
		nodes_polycube[3 * i + 0] = pts_polycube[i][0];
		nodes_polycube[3 * i + 1] = pts_polycube[i][1];
		nodes_polycube[3 * i + 2] = pts_polycube[i][2];
		if (tetra_vertices[i]->boundary)
		{
			ref_nodes[i] = 1;
			//node_normal_tag[3 * i] = 1;
		}
	}
	std::vector<int> elems_int;
	for (size_t i = 0; i < elems.size(); i++)
	{
		elems_int.push_back(elems[i]);
	}
	std::vector<int> normal_tag;
	normal_tag.resize(3 * nv, 0);
	if (!fix_boundary)
	{
		//change normal tag and ref node
		std::vector<OpenVolumeMesh::Geometry::Vec3i> v_type = (*pq_flattening).get_vertex_type();
		if (v_type.size() != nv)
		{
			flattening_prepare_vtk();
			v_type = (*pq_flattening).get_vertex_type();
		}
		for (size_t i = 0; i < nv; i++)
		{
			if (v_type[i][0] >= 0)
			{
				normal_tag[3 * i] = 1;
			}
			if (v_type[i][1] >= 0)
			{
				normal_tag[3 * i + 1] = 1;
			}
			if (v_type[i][2] >= 0)
			{
				normal_tag[3 * i + 2] = 1;
			}
			ref_nodes[i] = 0;
		}
		//add seamless constraints here
		for (size_t i = 0; i < equal_faces_array.size(); i++)
		{
			for (size_t j = 0; j < equal_faces_array[i].size(); j++)
			{
				for (size_t k = 0; k < equal_faces_array[i][j].size(); k++)
				{
					ref_nodes[equal_faces_array[i][j][k]] = 1;
				}
			}
		}
	}
	SUSSpace::SUS_standard_polycube(nodes_polycube, ref_nodes, elems_int, normal_tag, output_coord, smooth_iter, max_iter);
	for (size_t i = 0; i < nv; i++)
	{
		pts_polycube[i][0] = output_coord[3 * i + 0];
		pts_polycube[i][1] = output_coord[3 * i + 1];
		pts_polycube[i][2] = output_coord[3 * i + 2];
	}
	tet_->load_tet(pts_polycube, elems);
	
}
void Polycube_Algorithm::do_optimize_polycube_PC(bool fix_boundary, int max_smooth_iter, int max_iter)
{
	//fix boundary, optimize polycube based on tet mesh
	if (!polycube_exist_flag)
		return;
	//optimize origin mesh based on polycube mesh
	std::vector<double> nodes_polycube;
	std::vector<double> nodes_ori;
	std::vector<int> ref_nodes;
	std::vector<unsigned int> elems;
	std::vector<int> node_normal_tag;
	std::vector<double> output_coord;
	std::vector<CVec<double, 3>> pts_polycube, pts_ori;
	int nv = (int)tet_mesh_polycube_->tetra_vertices.size();
	int nc = (int)tet_mesh_polycube_->tetras.size();
	nodes_polycube.resize(3 * nv);
	nodes_ori.resize(3 * nv);
	ref_nodes.resize(nv, 0);
	node_normal_tag.resize(3 * nv, 0);
	get_vert_cell(tet_mesh_polycube_, pts_polycube, elems);
	get_vert_cell(tet_mesh_, pts_ori, elems);
	assert(pts_polycube.size() == pts_ori.size());
	const std::vector<Tetrahedron<double>*> &tetras = tet_mesh_->tetras;
	const std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
	for (size_t i = 0; i < nv; i++)
	{
		nodes_polycube[3 * i + 0] = pts_polycube[i][0];
		nodes_polycube[3 * i + 1] = pts_polycube[i][1];
		nodes_polycube[3 * i + 2] = pts_polycube[i][2];
		nodes_ori[3 * i + 0] = pts_ori[i][0];
		nodes_ori[3 * i + 1] = pts_ori[i][1];
		nodes_ori[3 * i + 2] = pts_ori[i][2];
		if (tetra_vertices[i]->boundary)
		{
			ref_nodes[i] = 1;
		}
	}
	//equal_face_constrait;
	if (!fix_boundary)
	{
		//change normal_tag and ref_node
		std::vector<OpenVolumeMesh::Geometry::Vec3i> v_type = (*pq_flattening).get_vertex_type();
		for (size_t i = 0; i < nv; i++)
		{
			if (v_type[i][0] >= 0)
			{
				node_normal_tag[3 * i] = 1;
			}
			if (v_type[i][1] >= 0)
			{
				node_normal_tag[3 * i + 1] = 1;
			}
			if (v_type[i][2] >= 0)
			{
				node_normal_tag[3 * i + 2] = 1;
			}
			ref_nodes[i] = 0;
		}
		for (size_t i = 0; i < equal_faces_array.size(); i++)
		{
			for (size_t j = 0; j < equal_faces_array[i].size(); j++)
			{
				for (size_t k = 0; k < equal_faces_array[i][j].size(); k++)
				{
					ref_nodes[equal_faces_array[i][j][k]] = 1;
				}
			}
		}
	}
	SUSPolycube::SUS(nodes_polycube, ref_nodes, elems, node_normal_tag, nodes_ori, output_coord, max_smooth_iter, max_iter);
	//after this, ori_pts will be normalized
	for (size_t i = 0; i < nv; i++)
	{
		pts_polycube[i][0] = output_coord[3 * i + 0];
		pts_polycube[i][1] = output_coord[3 * i + 1];
		pts_polycube[i][2] = output_coord[3 * i + 2];
	}
	tet_mesh_polycube_->load_tet(pts_polycube, elems);
	//uppdate_vertex_position_signal();
	(*pq_deform).assign_pos_mesh(tet_mesh_polycube_, true);
	//(*pd_flatten).assign_pos_mesh(tet_mesh_polycube_, true);
	(*pq_deform).compute_distortion(tet_mesh_polycube_);
	compute_flips();
}
bool Polycube_Algorithm::do_Extract_polycube_vtk(int distortion_type, bool sample_hex)
{
	if (!tet_mesh_polycube_)
		tet_mesh_polycube_ = tet_mesh_;
	//if not prepared
	flattening_prepare_vtk();
	if (chart_mean_value_ori.size() == 0)
	{
		(*pq_flattening).find_all_corner(tet_mesh_polycube_);
		chart_mean_value_ori = (*pq_flattening).get_chart_mean_value();
	}
	(*pq_flattening).assign_pos_mesh(tet_mesh_polycube_, true);
	std::vector<std::vector<int>> edges = (*pq_flattening).polycube_edges;
	std::vector<double> polycube_edge_distortion_conf, polycube_edge_distortion_vol;
	(*pq_flattening).get_polycube_edge_layer_distortion(3, tet_mesh_);
	polycube_edge_distortion = (*pq_flattening).polycube_edge_layer_distortion;
	polycube_edge_distortion_conf = (*pq_flattening).polycube_edge_distortion_conf;
	polycube_edge_distortion_vol = (*pq_flattening).polycube_edge_distortion_vol;
	if (distortion_type == 1)
		polycube_edge_distortion = polycube_edge_distortion_conf;
	else if (distortion_type == 2)
		polycube_edge_distortion = polycube_edge_distortion_vol;
	polycube_edge_length = (*pq_flattening).polycube_edge_length;
	get_min_polycube_edge_length();
	if (sample_hex)
	{
		do_sample_hex();
		polycube_edge_distortion = polycube_edge_hex_distortion;
	}
	//get cut depth of each polycube edge here
	polycube_edge_neighbor_min_length.clear();
	polycube_edge_neighbor_min_length.resize(edges.size());
	polycube_edge_to_idx_map.clear();
	for (int i = 0; i < edges.size(); i++)
	{
		std::vector<int> cand_idx;
		int id0, id1;
		id0 = edges[i][0];
		id1 = edges[i][1];
		for (int j = 0; j < edges.size(); j++)
		{
			if (j == i)
				continue;
			if (edges[j][0] == id0 || edges[j][0] == id1 || edges[j][1] == id0 || edges[j][1] == id1)
				cand_idx.push_back(j);
		}
		//not suitable for valence 3
		double min_length = polycube_edge_length[cand_idx[0]];
		for (int j = 0; j < cand_idx.size(); j++)
		{
			if (min_length > polycube_edge_length[cand_idx[j]])
				min_length = polycube_edge_length[cand_idx[j]];
		}
		polycube_edge_neighbor_min_length[i] = min_length;
	}
	std::vector<int> chart_label_map = (*pq_flattening).get_polycube_chart_label();
	std::vector<std::vector<std::pair<int, int>>> edges_with_same_chart;
	int max_chart_idx = 0;
	for (size_t i = 0; i < edges.size(); i++)
	{
		if (max_chart_idx < edges[i][2])
			max_chart_idx = edges[i][2];
		if (max_chart_idx < edges[i][3])
			max_chart_idx = edges[i][3];
	}
	edges_with_same_chart.resize(max_chart_idx + 1);
	for (size_t i = 0; i < edges.size(); i++)
	{
		int chart_idx1 = edges[i][2];
		int chart_idx2 = edges[i][3];
		std::pair<int, int> point_pair(edges[i][0], edges[i][1]);
		edges_with_same_chart[chart_idx1].push_back(point_pair);
		edges_with_same_chart[chart_idx2].push_back(point_pair);
	}
	std::vector<std::vector<std::vector<int>>> faces_array;
	std::vector<int> faces_chart;
	//construct face and compute face_normal
	for (size_t i = 0; i < edges_with_same_chart.size(); i++)
	{
		std::vector<std::pair<int, int>> &one_chart = edges_with_same_chart[i];
		std::list<std::pair<int, int>> one_chart_list(one_chart.begin(), one_chart.end());
		std::vector<int> temp_face;
		std::vector<std::vector<int>> temp_faces;
		int next_idx;
		while (!one_chart_list.empty())
		{
			temp_face.clear();
			temp_face.push_back(one_chart_list.front().first);
			temp_face.push_back(one_chart_list.front().second);
			one_chart_list.pop_front();
			next_idx = -1;
			int max_try_time = (int)one_chart_list.size();
			int try_time = 0;
			while (!one_chart_list.empty() && try_time < max_try_time && next_idx != temp_face[0])
			{
				for (std::list<std::pair<int, int>>::iterator it = one_chart_list.begin(); it != one_chart_list.end(); ++it)
				{
					if ((*it).first == temp_face.back())
					{
						next_idx = (*it).second;
						one_chart_list.erase(it);
						break;
					}
					if ((*it).second == temp_face.back())
					{
						next_idx = (*it).first;
						one_chart_list.erase(it);
						break;
					}
				}
				try_time++;
				if (next_idx != temp_face.back() && next_idx != -1)
					temp_face.push_back(next_idx);
			}
			if (temp_face[0] == temp_face.back())
			{
				//only insert face which are close
				temp_faces.push_back(temp_face);
			}
		}
		faces_array.push_back(temp_faces);
	}
	std::map<int, int> all_vert_no_dup_map;
	std::vector<int> all_vert;
	int all_vert_idx = 0;
	std::map<int, int>::iterator all_vert_it;
	int total_face_number = 0;
	for (size_t i = 0; i < faces_array.size(); i++)
	{
		//remove the last element
		for (size_t j = 0; j < faces_array[i].size(); j++)
		{
			for (size_t k = 0; k < faces_array[i][j].size(); k++)
			{
				all_vert_it = all_vert_no_dup_map.find(faces_array[i][j][k]);
				if (all_vert_it == all_vert_no_dup_map.end())
				{
					//add new map element;
					all_vert.push_back(faces_array[i][j][k]);
					all_vert_no_dup_map[faces_array[i][j][k]] = all_vert_idx;
					all_vert_idx++;
				}
			}
			total_face_number++;
			if (faces_array[i][j][0] == faces_array[i][j].back())
				faces_array[i][j].pop_back();
		}
	}
	std::vector<double> coord_x, coord_y, coord_z;
	(*pq_flattening).get_coord(all_vert, coord_x, coord_y, coord_z);
	const double chart_normals[6][3] = { {1.f,0.f,0.f},{-1.f,0.f,0.f},{0.f,1.f,0.f},{0.f,-1.f,0.f},{0.f,0.f,1.f},{0.f,0.f,-1.f} };
	assert(faces_array.size() == edges_with_same_chart.size());
	std::vector<std::pair<int, int>>::iterator edge_it;
	//create face
	for (size_t i = 0; i < faces_array.size(); i++)
	{
		if (faces_array[i].size() == 0)
		{
			continue;
		}
		//faces with holes
		std::vector<double> face_coord_x, face_coord_y, face_coord_z;
		for (size_t k = 0; k < faces_array[i].size(); k++)
		{
			int face_size = (int)faces_array[i][k].size();
			for (size_t j = 0; j < face_size; j++)
			{
				face_coord_x.push_back(coord_x[all_vert_no_dup_map[faces_array[i][k][j]]]);
				face_coord_y.push_back(coord_y[all_vert_no_dup_map[faces_array[i][k][j]]]);
				face_coord_z.push_back(coord_z[all_vert_no_dup_map[faces_array[i][k][j]]]);
			}
		}
		//construct polygons
		std::vector<std::vector<int>> new_faces;
		if (ig::SimpleTriangulation::sortface_area(face_coord_x, face_coord_y, face_coord_z, faces_array[i]))
		{
			//recompute face coord
			face_coord_x.clear();
			face_coord_y.clear();
			face_coord_z.clear();
			for (size_t k = 0; k < faces_array[i].size(); k++)
			{
				int face_size = (int)faces_array[i][k].size();
				for (size_t j = 0; j < face_size; j++)
				{
					face_coord_x.push_back(coord_x[all_vert_no_dup_map[faces_array[i][k][j]]]);
					face_coord_y.push_back(coord_y[all_vert_no_dup_map[faces_array[i][k][j]]]);
					face_coord_z.push_back(coord_z[all_vert_no_dup_map[faces_array[i][k][j]]]);
				}
			}
		}
		bool triangulation_success = ig::SimpleTriangulation::triangulation(face_coord_x, face_coord_y, face_coord_z, faces_array[i], new_faces);
		if (!triangulation_success) return false;
		double normal_x, normal_y, normal_z;
		double dir1_x, dir1_y, dir1_z;
		double dir2_x, dir2_y, dir2_z;
		int normal_begin_idx = new_faces[0][0];
		int next_idx = new_faces[0][1];
		int nextnext_idx = new_faces[0][2];
		dir1_x = coord_x[all_vert_no_dup_map[normal_begin_idx]] - coord_x[all_vert_no_dup_map[next_idx]];
		dir1_y = coord_y[all_vert_no_dup_map[normal_begin_idx]] - coord_y[all_vert_no_dup_map[next_idx]];
		dir1_z = coord_z[all_vert_no_dup_map[normal_begin_idx]] - coord_z[all_vert_no_dup_map[next_idx]];
		dir2_x = coord_x[all_vert_no_dup_map[nextnext_idx]] - coord_x[all_vert_no_dup_map[next_idx]];
		dir2_y = coord_y[all_vert_no_dup_map[nextnext_idx]] - coord_y[all_vert_no_dup_map[next_idx]];
		dir2_z = coord_z[all_vert_no_dup_map[nextnext_idx]] - coord_z[all_vert_no_dup_map[next_idx]];
		ig::CGALHelper::CrossVectorNormalize(dir1_x, dir1_y, dir1_z, dir2_x, dir2_y, dir2_z, normal_x, normal_y, normal_z);
		int face_label = chart_label_map[i];
		//check if face normal are the same with normal 
		const double *chart_normal = chart_normals[face_label];
		double vector_dot = chart_normal[0] * normal_x + chart_normal[1] * normal_y + chart_normal[2] * normal_z;
		if (vector_dot < 0)
		{
			//same direction
			for (size_t j = 0; j < new_faces.size(); j++)
			{
				//std::reverse(new_faces[j].begin(), new_faces[j].end());
				faces_chart.push_back((int)i);
			}
			faces_array[i] = new_faces;
		}
		else
		{
			for (size_t j = 0; j < new_faces.size(); j++)
			{
				std::reverse(new_faces[j].begin(), new_faces[j].end());
				faces_chart.push_back((int)i);
			}
			faces_array[i] = new_faces;
		}
	}
	total_face_number = 0;
	for (size_t i = 0; i < faces_array.size(); i++)
	{
		total_face_number += (int)faces_array[i].size();
	}
	cut_input_point.clear();
	cut_input_face.clear();
	cut_input_chart.clear();
	cut_input_label.clear();
	for (size_t i = 0; i < all_vert.size(); i++)
	{
		CVec<double, 3> temp(coord_x[i], coord_y[i], coord_z[i]);
		cut_input_point.push_back(temp);
	}
	for (size_t i = 0; i < faces_array.size(); i++)
	{
		for (size_t j = 0; j < faces_array[i].size(); j++)
		{
			int face_size = (int)faces_array[i][j].size();
			for (size_t k = 0; k < faces_array[i][j].size(); k++)
			{
				cut_input_face.push_back(all_vert_no_dup_map[faces_array[i][j][k]]);
			}
		}
	}
	//feature edge renew here
	if (!init_feature_edge_array.empty() && init_feature_polycube_edge_array_sorted.empty())
	{
		if (!refine_polycube_feature_edge_array(all_vert_no_dup_map)) return false;
	}
	//construct edge distortion and vertex distortion
	polycube_structure_e_distortion.clear();
	polycube_structure_v_distortion.clear();
	polycube_structure_v_distortion.resize(all_vert_no_dup_map.size(), 0);
	polycube_edge.clear();
	polycube_edge_label.clear();
	polycube_edge_chart.clear();
	for (size_t i = 0; i < edges.size(); i++)
	{
		int id1 = all_vert_no_dup_map[edges[i][0]];
		int id2 = all_vert_no_dup_map[edges[i][1]];
		int label1 = chart_label_map[edges[i][2]];
		int label2 = chart_label_map[edges[i][3]];
		if (id1 < id2)
			polycube_edge.push_back(std::pair<unsigned int, unsigned int>((unsigned)id1, (unsigned)id2));
		else
			polycube_edge.push_back(std::pair<unsigned int, unsigned int>((unsigned)id2, (unsigned)id1));
		polycube_edge_to_idx_map[polycube_edge.back()] = (unsigned)i;
		polycube_edge_chart.push_back(std::pair<int, int>(edges[i][2], edges[i][3]));
		polycube_edge_label.push_back(std::pair<int, int>(label1, label2));
		if (id1 < id2)
			polycube_structure_e_distortion[std::pair<int, int>(id1, id2)] = polycube_edge_distortion[i];
		else
			polycube_structure_e_distortion[std::pair<int, int>(id2, id1)] = polycube_edge_distortion[i];
		polycube_structure_v_distortion[id1] += polycube_edge_distortion[i];
		polycube_structure_v_distortion[id2] += polycube_edge_distortion[i];
	}
	//refine polycube_edge_neighbor_min_length
	std::vector<std::set<unsigned>> vv_search;
	int n_cut_input_face = (int)cut_input_face.size() / 3;
	vv_search.resize(cut_input_point.size());
	for (size_t i = 0; i < n_cut_input_face; i++)
	{
		unsigned tmp_face[3];
		for (size_t j = 0; j < 3; j++)
		{
			tmp_face[j] = cut_input_face[3 * i + j];
		}
		for (size_t j = 0; j < 3; j++)
		{
			unsigned tmp_center = tmp_face[j];
			unsigned tmp_neighbor1 = tmp_face[(j + 1) % 3];
			unsigned tmp_neighbor2 = tmp_face[(j + 2) % 3];
			vv_search[tmp_center].insert(tmp_neighbor1);
			vv_search[tmp_center].insert(tmp_neighbor2);
		}
	}
	for (size_t i = 0; i < polycube_edge.size(); i++)
	{
		unsigned id1 = polycube_edge[i].first;
		unsigned id2 = polycube_edge[i].second;
		for (auto it = vv_search[id1].begin(); it != vv_search[id1].end(); ++it)
		{
			if (*it != id2)
			{
				//not polycube edge
				double tmp_length = (cut_input_point[id1] - cut_input_point[*it]).L2Norm();
				if (tmp_length < polycube_edge_neighbor_min_length[i])
					polycube_edge_neighbor_min_length[i] = tmp_length;
			}
		}
		for (auto it = vv_search[id2].begin(); it != vv_search[id2].end(); ++it)
		{
			if (*it != id1)
			{
				//not polycube edge
				double tmp_length = (cut_input_point[id2] - cut_input_point[*it]).L2Norm();
				if (tmp_length < polycube_edge_neighbor_min_length[i])
					polycube_edge_neighbor_min_length[i] = tmp_length;
			}
		}
	}
	double max_v_distortion = polycube_structure_v_distortion[0];
	double min_v_distortion = polycube_structure_v_distortion[0];
	for (size_t i = 0; i < polycube_structure_v_distortion.size(); i++)
	{
		if (max_v_distortion < polycube_structure_v_distortion[i])
			max_v_distortion = polycube_structure_v_distortion[i];
		if (min_v_distortion > polycube_structure_v_distortion[i])
			min_v_distortion = polycube_structure_v_distortion[i];
	}
	double distortion_diff = max_v_distortion - min_v_distortion;
	if (distortion_diff > 0.0001 * EPSILON)
		for (size_t i = 0; i < polycube_structure_v_distortion.size(); i++)
		{
			polycube_structure_v_distortion[i] = (polycube_structure_v_distortion[i] - min_v_distortion) / distortion_diff;
		}
	else
		for (size_t i = 0; i < polycube_structure_v_distortion.size(); i++)
		{
			polycube_structure_v_distortion[i] = 0.0;
		}
	//sort polycube edge according to distortion
	sorted_polycube_edge.clear();
	sorted_polycube_edge_distortion.clear();
	sorted_polycube_edge = polycube_edge;
	sorted_polycube_edge_distortion = polycube_edge_distortion;
	//bubble sort
	for (int i = 0; i < polycube_edge.size(); i++)
	{
		int change_id = i;
		for (int j = i + 1; j < polycube_edge.size(); j++)
		{
			if (sorted_polycube_edge_distortion[change_id] < sorted_polycube_edge_distortion[j])
				change_id = j;
		}
		if (change_id != i)
		{
			//change element
			double temp_dist = sorted_polycube_edge_distortion[change_id];
			sorted_polycube_edge_distortion[change_id] = sorted_polycube_edge_distortion[i];
			sorted_polycube_edge_distortion[i] = temp_dist;
			std::pair<unsigned int, unsigned int> temp_pair = sorted_polycube_edge[change_id];
			sorted_polycube_edge[change_id] = sorted_polycube_edge[i];
			sorted_polycube_edge[i] = temp_pair;
		}
	}
	//construct polycube structure here
	do_construct_polycube_structure();
	cut_input_chart = faces_chart;
	for (size_t i = 0; i < faces_chart.size(); i++)
	{
		int chart = faces_chart[i];
		int label = chart_label_map[chart];
		cut_input_label.push_back(label);
	}
	pq_flattening->set_bd_face_chart_label(cut_input_point, cut_input_face, cut_input_chart, cut_input_label);
	//reset polycube_edge_max_cut_depth_one_cut
	if (!reset_max_cut_depth())
		return false;
	return true;
}
bool Polycube_Algorithm::refine_polycube_feature_edge_array(const std::map<int, int> &vertmap)
{
	//refine init_feature_polycube_edge_array, init_feature_polycube_edge_array_sorted
	std::vector<std::pair<unsigned, unsigned>> new_feature_polycube;
	for (size_t i = 0; i < init_feature_polycube_edge_array.size(); i++)
	{
		auto it0 = vertmap.find(init_feature_polycube_edge_array[i].first);
		auto it1 = vertmap.find(init_feature_polycube_edge_array[i].second);
		assert(it0 != vertmap.end() && it1 != vertmap.end());
		if (it0 == vertmap.end() || it1 == vertmap.end()) return false;
		unsigned id0 = it0->second, id1 = it1->second;
		new_feature_polycube.push_back(std::pair<unsigned, unsigned>(id0, id1));
		if (id0 > id1) std::swap(id0, id1);
		init_feature_polycube_edge_array_sorted.push_back(std::pair<unsigned, unsigned>(id0, id1));
	}
	assert(new_feature_polycube.size() == init_feature_polycube_edge_array.size());
	init_feature_polycube_edge_array = new_feature_polycube;
	return true;
}
void Polycube_Algorithm::do_sample_hex()
{
	int quality_type = 1;
	//after hex-flattening, sample points along each polycube edge
	polycube_edge_hex_distortion.clear();
	//normal vector
	ig::CVec<double, 3> n[6];
	ig::CVec<int, 3> dir[3][3];
	for (size_t i = 0; i < 3; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			for (size_t k = 0; k < 3; k++)
			{
				dir[i][j][k] = 0;
			}
		}
	}
	dir[0][0][1] = -1;
	dir[0][1][2] = -1;
	dir[0][2][1] = -1;
	dir[0][2][2] = -1;
	dir[1][0][0] = -1;
	dir[1][1][2] = -1;
	dir[1][2][0] = -1;
	dir[1][2][2] = -1;
	dir[2][0][0] = -1;
	dir[2][1][1] = -1;
	dir[2][2][0] = -1;
	dir[2][2][1] = -1;
	for (size_t i = 0; i < 6; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			n[i][j] = 0;
		}
	}
	n[0][0] = -1;
	n[1][0] = 1;
	n[2][1] = -1;
	n[3][1] = 1;
	n[4][2] = -1;
	n[5][2] = 1;
	//construct polycube edge hex distortion
	std::vector<std::vector<int>> edges = (*pq_flattening).polycube_edges;
	std::vector<int> polycube_chart_label = (*pq_flattening).get_polycube_chart_label();
	double cube_len = min_polycube_edge_length * sigma_r;
	//output coord
	for (size_t i = 0; i < edges.size(); i++)
		//for (size_t i = 0; i < 1; i++)
	{
		int avail_hex_count = 0;
		double avg_hex_quality = 0.0;
		int id1 = edges[i][0];
		int id2 = edges[i][1];
		ig::CVec<double, 3> pos1 = tet_mesh_polycube_->tetra_vertices[id1]->pos;
		ig::CVec<double, 3> pos2 = tet_mesh_polycube_->tetra_vertices[id2]->pos;
		int chart1 = edges[i][2];
		int chart2 = edges[i][3];
		int label1 = polycube_chart_label[chart1];
		int label2 = polycube_chart_label[chart2];
		int fix_axis = 3 - label1 / 2 - label2 / 2;
		ig::CVec<int, 3> begin[4], end[4];
		//direct rounding
		if (pos1[fix_axis] < pos2[fix_axis])
		{
			begin[0][0] = (int)round(pos1[0] / cube_len);
			begin[0][1] = (int)round(pos1[1] / cube_len);
			begin[0][2] = (int)round(pos1[2] / cube_len);
			end[0][0] = (int)round(pos2[0] / cube_len);
			end[0][1] = (int)round(pos2[1] / cube_len);
			end[0][2] = (int)round(pos2[2] / cube_len);
		}
		else
		{
			begin[0][0] = (int)round(pos2[0] / cube_len);
			begin[0][1] = (int)round(pos2[1] / cube_len);
			begin[0][2] = (int)round(pos2[2] / cube_len);
			end[0][0] = (int)round(pos1[0] / cube_len);
			end[0][1] = (int)round(pos1[1] / cube_len);
			end[0][2] = (int)round(pos1[2] / cube_len);
		}
		//first select four begin points
		int total_length = end[0][fix_axis] - begin[0][fix_axis] - 1;
		for (size_t j = 0; j < 3; j++)
		{
			begin[j + 1] = begin[0] + dir[fix_axis][j];
		}
		//iterate over all hexes
		for (size_t j = 0; j < total_length; j++)
		{
			int tmp_count = 0;
			for (size_t k = 0; k < 4; k++)
			{
				//use tmp center to decide whether a point is inside the volume
				ig::CVec<double, 3> tmp_center;
				tmp_center[0] = cube_len * (begin[k][0] + 0.5);
				tmp_center[1] = cube_len * (begin[k][1] + 0.5);
				tmp_center[2] = cube_len * (begin[k][2] + 0.5);
				tmp_center[fix_axis] = cube_len * (begin[k][fix_axis] + j * 1.0 + 0.5);
				ig::CVec<double, 3> tmp_center_int = tmp_center / cube_len;
				if (tet_mesh_polycube_->is_inside(tmp_center))
				{
					tmp_count++;
					//construct hex and map it to tet_mesh
					std::vector<ig::CVec<double, 3>> hex_elems, hex_ori_elems;
					hex_elems.resize(8);
					hex_ori_elems.resize(8);
					hex_elems[0] = tmp_center - ig::CVec<double, 3>(cube_len * 0.5, cube_len * 0.5, cube_len * 0.5);
					//hex_elems[0] = tmp_center;
					hex_elems[1] = hex_elems[0] + cube_len * ig::CVec<double, 3>(1.0, 0, 0);
					hex_elems[2] = hex_elems[1] + cube_len * ig::CVec<double, 3>(0, 1.0, 0);
					hex_elems[3] = hex_elems[0] + cube_len * ig::CVec<double, 3>(0, 1.0, 0);
					hex_elems[4] = hex_elems[0] + cube_len * ig::CVec<double, 3>(0, 0, 1.0);
					hex_elems[5] = hex_elems[4] + cube_len * ig::CVec<double, 3>(1.0, 0, 0);
					hex_elems[6] = hex_elems[5] + cube_len * ig::CVec<double, 3>(0, 1.0, 0);
					hex_elems[7] = hex_elems[4] + cube_len * ig::CVec<double, 3>(0, 1.0, 0);
					Hex hex;
					get_corresponding_point(tet_mesh_polycube_, tet_mesh_, hex_elems, hex_ori_elems);
					hex.elems = hex_ori_elems;
					double norm_v = hex.compute_quality(quality_type) / (cube_len * cube_len * cube_len);
					avg_hex_quality += norm_v + 1.0 / norm_v;
				}
			}
			avail_hex_count += tmp_count;
		}
		avg_hex_quality /= (avail_hex_count * 1.0);
		polycube_edge_hex_distortion.push_back(avg_hex_quality);
	}
}
void Polycube_Algorithm::do_deformation_constrained_vtk(double sigma_s)
{
	long start_t = clock();
	if (!tet_mesh_polycube_)
	{
		std::cout << "polycube pointer NULL" << std::endl;
		tet_mesh_polycube_ = tet_mesh_;
		(*pq_deform).prepare_for_deformation(tet_mesh_);
	}
	auto_color();
	long start_t1 = clock();
	long diff1 = start_t1 - start_t;
	double t1 = (double)(diff1) / CLOCKS_PER_SEC;
	std::cout << "coloring time: " << t1 << std::endl;
	//fix normal and fix triangle area
	(*pq_deform).set_sigma_s(sigma_s);
	for (int j = 0; j < 5; ++j)
	{
		//change_mu is set true in the first round, then set to false in the following rounds
		(*pq_deform).exp_mips_deformation_refine_polycube_omp_constrained(1, 15, 1.0, 0.5, (int)tet_mesh_polycube_->tetra_vertices.size(), false, true, false, true);
		(*pq_deform).assign_pos_mesh(tet_mesh_polycube_);
	}
	//no labeling here, or a simplified version of labeling
	(*pq_deform).compute_distortion(tet_mesh_polycube_);
	pq_deform->compute_chart_normal_diff();
	long end_t = clock();
	long diff = end_t - start_t;
	double t = (double)(diff) / CLOCKS_PER_SEC;
	printf("Polycube Map : %f s\n", t);
	return;
}
void Polycube_Algorithm::do_deformation_constrained_continued()
{
	//continue deform for 10 iteration, not considering flattening
	if (tet_mesh_polycube_)
	{
		for (int j = 0; j < 10; ++j)
		{
			(*pq_deform).exp_mips_deformation_refine_polycube_omp_constrained(1, 15, 1.0, 0.5, (int)tet_mesh_polycube_->tetra_vertices.size(), false, true, false, false);
			(*pq_deform).assign_pos_mesh(tet_mesh_polycube_);
		}
		//no labeling here, or a simplified version of labeling
		(*pq_deform).compute_distortion(tet_mesh_polycube_);
		pq_deform->compute_chart_normal_diff();
	}
	else
	{
		for (int j = 0; j < 10; ++j)
		{
			(*pq_deform).exp_mips_deformation_refine_polycube_omp_constrained(1, 15, 1.0, 0.5, mesh_->n_vertices(), false, true, false);
			(*pq_deform).assign_pos_mesh(mesh_);
		}
		(*pq_deform).compute_distortion(mesh_);
	}
	return;
}
void Polycube_Algorithm::do_construct_polycube_structure()
{
	if (cut_input_point.size() == 0 || cut_input_face.size() == 0)
		return;
	if (v_distortion == NULL)
		return;
	if (polycube_structure->n_vertices() != 0)
		polycube_structure->clear();
	polycube_structure->request_face_normals();
	polycube_structure->add_property(*v_distortion);
	std::vector<SurfaceMesh::VertexHandle> vhandle;
	std::vector<SurfaceMesh::FaceHandle> fhandle;
	vhandle.resize(cut_input_point.size());
	fhandle.resize(cut_input_face.size() / 2);
	for (size_t i = 0; i < cut_input_point.size(); i++)
	{
		vhandle[i] = polycube_structure->add_vertex(SurfaceMesh::Point(cut_input_point[i][0], cut_input_point[i][1], cut_input_point[i][2]));
		polycube_structure->property(*v_distortion, vhandle[i]) = polycube_structure_v_distortion[i];
	}
	std::vector<SurfaceMesh::VertexHandle> tmp_face_vhandles;
	int n_face = (int)cut_input_face.size() / 3;
	for (size_t i = 0; i < n_face; i++)
	{
		tmp_face_vhandles.clear();
		for (size_t j = 0; j < 3; j++)
		{
			tmp_face_vhandles.push_back(vhandle[cut_input_face[3 * i + j]]);
		}
		fhandle[i] = polycube_structure->add_face(tmp_face_vhandles);
	}
	polycube_structure->update_normals();
}
void Polycube_Algorithm::flattening_prepare_vtk()
{
	if (tet_mesh_ == NULL) return;
	if (!tet_mesh_polycube_)
		tet_mesh_polycube_ = tet_mesh_;
	if (min_polycube_edge_length < 0.0)
	{
		get_boundary_avg_edge_length();
		min_polycube_edge_length = avg_boundary_edge_length;
		get_bounding_box();
		double bb[3];
		for (size_t i = 0; i < 3; i++)
		{
			bb[i] = (bb_Max[i] - bb_Min[i]) / N_SPLIT;
		}
		double max_bb = std::max(std::max(bb[0], bb[1]), bb[2]);
		if (min_polycube_edge_length > max_bb)
			min_polycube_edge_length = max_bb;
	}
	int nv = (int)tet_mesh_->tetra_vertices.size();
	if (coord_ori.size() != 3 * nv)
	{
		coord_ori.clear();
		std::vector<TetVertex<double>* > tetra_vertices = tet_mesh_->tetra_vertices;
		//flattening based on tet_mesh_polycube
		for (size_t i = 0; i < nv; i++)
		{
			double x, y, z;
			x = tetra_vertices[i]->pos[0];
			y = tetra_vertices[i]->pos[1];
			z = tetra_vertices[i]->pos[2];
			coord_ori.push_back(x);
			coord_ori.push_back(y);
			coord_ori.push_back(z);
		}
	}
	//here tet_mesh_ should be deformed mesh
	(*pq_flattening).load_deformation_result(coord_ori, tet_mesh_polycube_);
	(*pq_flattening).set_sigma_r(sigma_r);
	//in the following part, fix_labels has the following correspondence:
	//0 x, 1 y, 2 z, 3 -x, 4 -y, 5 -z
	int m[6] = { 0, 3, 1, 4, 2, 5 };
	std::vector<int> fix_labels_changed;
	for (size_t i = 0; i < fix_labels.size(); i++)
	{
		if (fix_labels[i] >= 0)
			fix_labels_changed.push_back(m[fix_labels[i]]);
		else
			fix_labels_changed.push_back(-1);
	}
	(*pq_flattening).compute_distortion(tet_mesh_polycube_);
	(*pq_flattening).load_boundary_face_label(fix_labels_changed, fix_charts, tet_mesh_polycube_);
	if (init_distortion[0] < 0)
		init_distortion = (*pq_flattening).compute_distortion(tet_mesh_polycube_);
	if (init_volume_distortion[0] < 0)
		init_volume_distortion = (*pq_flattening).compute_volumetric_distortion(tet_mesh_polycube_);
}
bool Polycube_Algorithm::do_flattening_constrained_vtk(int distortion_type, double default_cube_length, double lambda, bool add_constraint, int IVF_polycube_type, int set_min_diff, bool int_solver_flag, double offset, bool double_cube_length_flag)
{
	double sigma_r_backup = sigma_r;
	//if hex meshing defined, polycube must be provided
	//if deformation not prepared
	if (!pq_deform->deformation_prepare_OK)
	{
		(*pq_deform).prepare_for_deformation(tet_mesh_);
	}
	if (polycube_exist_flag)
	{
		if (min_polycube_edge_length < 0.0)
		{
			//set min_polycube_edge_length here
			if (!do_Extract_polycube_vtk(distortion_type)) return false;
		}
		else
		{
			flattening_prepare_vtk();
		}
	}
	else
	{
		flattening_prepare_vtk();
	}
	double initial_cube_length = min_polycube_edge_length * sigma_r;
	if (default_cube_length < 0.0)
		pq_flattening->set_cube_length(initial_cube_length);
	else
	{
		pq_flattening->set_cube_length(default_cube_length);
		initial_cube_length = default_cube_length;
	}
	std::vector<double> coord;
	TetStructure<double>* temp_tet_ptr = tet_mesh_;
	tet_mesh_ = tet_mesh_polycube_;
	if (!tet_mesh_polycube_)
	{
		tet_mesh_polycube_ = tet_mesh_;
	}
	tet_mesh_ = temp_tet_ptr;
	temp_tet_ptr = NULL;
	//flattening part
	double energy_power, area_angle_ratio;
	energy_power = 1.0;
	area_angle_ratio = 0.5;
	//temporarily set max iteration time as 10
	int max_iter = 10;
	long start_t = clock();
	(*pq_flattening).assign_pos_mesh(tet_mesh_polycube_);
	(*pq_flattening).find_all_corner(tet_mesh_polycube_);
	(*pq_flattening).compute_triangle_area();
	int shrink_count = 0;
	while (1)
	{
		bool flatten_flag;
		if (chart_mean_value_ori.size() == 0)
			flatten_flag = (*pq_flattening).deform_ARAP_polycube_equal_face(tet_mesh_polycube_, int_solver_flag, add_constraint, set_min_diff, offset, double_cube_length_flag);
		else
			flatten_flag = (*pq_flattening).deform_ARAP_polycube_equal_face(tet_mesh_polycube_, chart_mean_value_ori, lambda, add_constraint, set_min_diff, offset, double_cube_length_flag);
		if (flatten_flag)
			break;
		initial_cube_length = initial_cube_length * 0.5;
		std::cout << "Adjusting cube_length: " << initial_cube_length << std::endl;
		pq_flattening->set_cube_length(initial_cube_length);
		shrink_count++;
		if (shrink_count >= 20)
		{
			std::cout << "flatten failed: " << std::endl;
			return false;
		}
	}
	(*pq_flattening).compute_triangle_area();
	long end_t = clock();
	long diff = end_t - start_t;
	double t = (double)(diff) / CLOCKS_PER_SEC;
	printf("Flatten Boundary : %f s\n", t);
	(*pq_flattening).compute_distortion(tet_mesh_polycube_);
	//compute distortion
	if (IVF_polycube_type == 0)
	{
		do_IVF_standard(tet_mesh_polycube_, false, 0, 20);
	}
	else if (IVF_polycube_type == 1)
	{
		do_optimize_polycube_PC(true, 0, 20);
	}
	if (polycube_exist_flag)
	{
		do_optimize_polycube_PC(true, 1, 2);
	}
	(*pq_flattening).assign_pos_mesh(tet_mesh_polycube_, true);
	(*pq_flattening).compute_triangle_area();
	(*pq_flattening).compute_distortion(tet_mesh_polycube_);
	//update mesh information in (*pd_normal) and (*pd_flatten)
	(*pq_deform).assign_pos_mesh(tet_mesh_polycube_, true);
	compute_flips();
	final_cube_length = initial_cube_length;
	return true;
}
bool Polycube_Algorithm::load_feature_edges(const char* filename)
{
	init_feature_edge_array.clear();
	init_feature_polycube_edge_array.clear();
	std::ifstream ifs(filename);
	if (ifs.is_open())
	{
		//grouped or not
		std::string s;
		std::getline(ifs, s);
		if (s == "Grouped Tet Feature")
		{
			//group version
			int ng = 0;
			ifs >> ng;
			for (size_t i = 0; i < ng; i++)
			{
				int one_group_size = 0;
				ifs >> one_group_size;
				for (size_t j = 0; j < one_group_size; j++)
				{
					std::pair<int, int> one_edge;
					ifs >> one_edge.first >> one_edge.second;
					init_feature_edge_array.push_back(one_edge);
				}
			}
			//feature polycube edge part
			if (!ifs.eof())
			{
				std::getline(ifs, s);
				int np = 0;
				ifs >> np;
				for (size_t i = 0; i < np; i++)
				{
					std::pair<int, int> one_edge;
					ifs >> one_edge.first >> one_edge.second;
					init_feature_polycube_edge_array.push_back(one_edge);
				}
			}
		}
		else
		{
			int n_fea = std::stoi(s);
			for (size_t i = 0; i < n_fea; i++)
			{
				std::pair<int, int> one_edge;
				ifs >> one_edge.first >> one_edge.second;
				init_feature_edge_array.push_back(one_edge);
			}
			if (!ifs.eof())
			{
				//get line twice
				do
				{
					std::getline(ifs, s);
					if (ifs.eof()) break;
				} while (s != "Feature PolyCube Edge");
				int np = 0;
				ifs >> np;
				for (size_t i = 0; i < np; i++)
				{
					std::pair<unsigned int, unsigned> one_edge;
					ifs >> one_edge.first >> one_edge.second;
					init_feature_polycube_edge_array.push_back(one_edge);
				}
			}
			init_feature_polycube_edge_array_ori = init_feature_polycube_edge_array;
		}
	}
	else
	{
		return false;
	}
	ifs.close();
	//set feature edges
	(*pq_deform).set_feature_edges_ovm(init_feature_edge_array, mesh_);
	//feature edges also set for pd_flatten
	feature_edge_flag = &pq_deform->feature_edge_flag;
	data_prepare_status = false;
	data_preparation_ovm();
	return true;
}
void Polycube_Algorithm::save_feature_edges_ovm(const char* filename)
{
	//save file according to feature edge flag
	if ((*feature_edge_flag).empty()) return;
	std::vector<int> feature_edge_array;
	for (size_t i = 0; i < (*feature_edge_flag).size(); i++)
	{
		if ((*feature_edge_flag)[i] == true) feature_edge_array.push_back((int)i);
	}
	std::ofstream ofs(filename);
	ofs << feature_edge_array.size() << std::endl;
	for (size_t i = 0; i < feature_edge_array.size(); i++)
	{
		OpenVolumeMesh::EdgeHandle eh = OpenVolumeMesh::EdgeHandle(feature_edge_array[i]);
		OpenVolumeMesh::OpenVolumeMeshEdge ovme = mesh_->edge(eh);
		int id0 = ovme.from_vertex();
		int id1 = ovme.to_vertex();
		ofs << id0 << " " << id1 << std::endl;
	}
	ofs.close();
}
void Polycube_Algorithm::save_feature_edges_ovm_vtkformat(const char* filename)
{
	//save file according to feature edge flag
	if ((*feature_edge_flag).empty()) return;
	std::vector<int> feature_edge_array;
	for (size_t i = 0; i < (*feature_edge_flag).size(); i++)
	{
		if ((*feature_edge_flag)[i] == true) feature_edge_array.push_back((int)i);
	}
	std::vector<double> x, y, z;
	for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		int v_id = v_it->idx();
		OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(*v_it);
		x.push_back(p[0]);
		y.push_back(p[1]);
		z.push_back(p[2]);
	}
	std::ofstream outputfile(filename);
	outputfile << "# vtk DataFile Version 3.0\n"
		<< "mesh vtk data\n"
		<< "ASCII\n"
		<< "DATASET POLYDATA\n";
	int n_point = (int)x.size();
	outputfile << "POINTS " << n_point << " double" << std::endl;
	for (size_t i = 0; i < n_point; i++)
	{
		outputfile << x[i] << " " << y[i] << " " << z[i] << std::endl;
	}
	std::vector<int> feature_point_color(n_point, 0);
	outputfile << "LINES " << feature_edge_array.size() << " " << 3 * feature_edge_array.size() << std::endl;
	for (size_t i = 0; i < feature_edge_array.size(); i++)
	{
		OpenVolumeMesh::EdgeHandle eh = OpenVolumeMesh::EdgeHandle(feature_edge_array[i]);
		OpenVolumeMesh::OpenVolumeMeshEdge ovme = mesh_->edge(eh);
		int id0 = ovme.from_vertex();
		int id1 = ovme.to_vertex();
		feature_point_color[id0] = 1;
		feature_point_color[id1] = 1;
		outputfile << "2 " << id0 << " " << id1 << std::endl;
	}
	outputfile << "POINT_DATA " << n_point << "\n"
		<< "SCALARS V_Scalars int\nLOOKUP_TABLE V_Table" << std::endl;
	for (size_t i = 0; i < n_point; i++)
	{
		outputfile << feature_point_color[i] << std::endl;
	}
	outputfile.close();
}
void Polycube_Algorithm::repair_features(bool polyline_guidance, bool ns_flag, double smooth_factor, double threshold)
{
	//use polyline as guidance
	if (polyline_guidance && pq_deform->polyline_target_n.empty())
	{
		deformation_polylines(3, 1.0, false, ns_flag, smooth_factor, threshold);
	}
	pq_deform->repair_features(mesh_, polyline_guidance);
}
