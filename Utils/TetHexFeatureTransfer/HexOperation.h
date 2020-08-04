#pragma once

//Created by Haoxiang Guo 2019.9.20 version 1
//ghx17@mails.tsinghua.edu.cn
//The only exterior dependency are Eigen

#include <map>
#include "HexStructure.h"

namespace HexStructure
{
	bool read_hex_mesh_vtk(Mesh &m, char * path);
	void initial_hex_mesh(const vector<ig::CVec<Real, 3>> &pts, const vector<std::vector<Int>> &idx, Mesh &m);
	void get_hex_vert_cell(const Mesh &m, vector<ig::CVec<Real, 3>> &pts, vector<std::vector<Int>> &idx);
	void write_hex_mesh_vtk(const Mesh &m, char * path);
	void write_singularity_vtk(const Mesh &m, const Singularity &s, char * path);
	void write_metamesh_vtk(const Mesh &m, const MetaMesh & mm, char * path);
	
	void build_connectivity(Mesh &m);
	void extract_singularity(Mesh &m, Singularity &s);
	bool extract_metamesh(Mesh &m, Singularity &s, MetaMesh &mm);

	bool extract_sheet(const Mesh &m, const vector<Int> &cand, vector<vector<Int>> &output, vector<std::pair<Int, Int>> &edge2sheet);

	void get_consistent_rotation(const Mesh &m, const vector<Int> &sheet, vector<Int>& rotationflag, bool refine = false);
	Int face_cell_consistency(const Face &f, const Cell &c);
	Int get_boundary_face_dir(const Mesh &m, Int bf_id);
	bool get_neighbor_face_consistency(const Mesh &m, Int id1, Int id2);

	void compute_cell_quality(Mesh &m, int quality_type);
	
}


