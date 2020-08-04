#pragma once

#include <vector>
#include "SmallVec.h"

//definition for single hex

using ig::CVec;

class Hex
{
public:
	Hex();
	Hex(const std::vector<CVec<double, 3>> &elem_pts) { elems = elem_pts; }
	~Hex();
	double compute_quality(int type);

	bool check_cube_quality(const ig::CVec<double, 3>& A, const ig::CVec<double, 3>& B, const ig::CVec<double, 3>& C, const ig::CVec<double, 3>& D, const ig::CVec<double, 3>& E, const ig::CVec<double, 3>& F, const ig::CVec<double, 3>& G, const ig::CVec<double, 3>& H, double& quality);
	
	bool check_spatial_quad_quality(const ig::CVec<double, 3>& v0, const ig::CVec<double, 3>& v1, const ig::CVec<double, 3>& v2, const ig::CVec<double, 3>& v3, double& quality);

	void compute_tet_volume_quality(const ig::CVec<double, 3>& v0, const ig::CVec<double, 3>& v1, const ig::CVec<double, 3>& v2, const ig::CVec<double, 3>& v3, double& subvolume, double& quality);

	double compute_cell_volume(const CVec<double, 3> &v1, const CVec<double, 3> &v2, const CVec<double, 3> &v3, const CVec<double, 3> &v4);

	std::vector<ig::CVec<double, 3>> elems;


};
