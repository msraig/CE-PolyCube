#include <algorithm>
#include <math.h>
#include "Hex.h"
//#include "Helper.h"

Hex::Hex()
{
	;
}

Hex::~Hex()
{
	;
}

double Hex::compute_quality(int type)
{
	if (elems.size() != 8)
		return -1;

	//compute quality based on diff type
	//0: scaled jacobian
	//1: volume
	//2: quality
	if (type == 0)
	{
		//compute scaled jacobian
		static const unsigned int corner_vert[][4] = { {0,1,3,4 }, { 1,2,0,5 },
		{ 3,0,2,7}, { 2,3,1,6 }, { 4,5,0,7 }, { 5,1,4,6 }, {7,6,4,3}, {6,5,7,2} };

		double min_corner_jacobian = 2.0;
		for (int j = 0; j < 8; j++)
		{
			ig::CVec<double, 3> e1 = elems[corner_vert[j][1]] -elems[corner_vert[j][0]];
			ig::CVec<double, 3> e2 = elems[corner_vert[j][2]] -elems[corner_vert[j][0]];
			ig::CVec<double, 3> e3 = elems[corner_vert[j][3]] -elems[corner_vert[j][0]];
			e1.Normalize(); e2.Normalize(); e3.Normalize();
			double tmp_corner_jacobian = e3.Dot(e1.Cross(e2));
			if (tmp_corner_jacobian < min_corner_jacobian)
			{
				min_corner_jacobian = tmp_corner_jacobian;
			}
		}
		return min_corner_jacobian;


	}
	else if (type == 1)
	{
		//compute volumes here
		static const unsigned int face_id[][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 },
		{ 7, 6, 2, 3 }, { 5, 4, 0, 1 }, { 5, 1, 2, 6 }, { 4, 7, 3, 0 } };

		ig::CVec<double, 3> hex_center(0, 0, 0);
		for (int j = 0; j < 8; j++)
		{
			hex_center += elems[j];
		}
		hex_center *= 0.125;
		//m_hex_centers[i] = hex_center;

		ig::CVec<double, 3> face_center[6];
		for (int j = 0; j < 6; j++)
		{
			face_center[j] = ig::CVec<double, 3>(0, 0, 0);
			for (int k = 0; k < 4; k++)
			{
				//face_center[j] += vertices[elements[i]->indices[face_id[j][k]]];
				face_center[j] += elems[face_id[j][k]];
			}
			face_center[j] *= 0.25;
		}

		double hex_vol = 0.0;
		for (int j = 0; j < 6; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				ig::CVec<double, 3> pt1 = elems[face_id[j][k]];
				ig::CVec<double, 3> pt2 = elems[face_id[j][(k + 1) % 4]];
				double tet_vol = compute_cell_volume(pt2, pt1, face_center[j], hex_center);
				hex_vol += tet_vol;
			}
		}

		
		if (!isnan(hex_vol))
			return hex_vol;
		else
			return -1.0;
	}
	else if (type == 2)
	{
		//cube quality
		double tmp_quality = -1.0;
		bool suc = check_cube_quality(elems[0], elems[1], elems[2], elems[3], elems[4], elems[5], elems[6], elems[7], tmp_quality);
		if (suc)
			if (!isnan(tmp_quality))
				return tmp_quality;

	}

	
	return -1.0;
	
}

bool Hex::check_cube_quality(const ig::CVec<double, 3>& A, const ig::CVec<double, 3>& B, const ig::CVec<double, 3>& C, const ig::CVec<double, 3>& D, const ig::CVec<double, 3>& E, const ig::CVec<double, 3>& F, const ig::CVec<double, 3>& G, const ig::CVec<double, 3>& H, double& quality)
{
	int m_qualitytype = 0;
	//check the validity of the cube
	double planar_quality[6];
	if (!check_spatial_quad_quality(A, B, C, D,
		planar_quality[0]))
		return false;
	if (!check_spatial_quad_quality(E, F, G, H,
		planar_quality[1]))
		return false;
	if (!check_spatial_quad_quality(A, B, F, E,
		planar_quality[2]))
		return false;
	if (!check_spatial_quad_quality(D, C, G, H,
		planar_quality[3]))
		return false;
	if (!check_spatial_quad_quality(B, C, G, F,
		planar_quality[4]))
		return false;
	if (!check_spatial_quad_quality(A, D, H, E,
		planar_quality[5]))
		return false;
	double subvolume[8], tet_quality[8];
	compute_tet_volume_quality(A, B, D, E, subvolume[0],
		tet_quality[0]);
	compute_tet_volume_quality(B, C, A, F, subvolume[1],
		tet_quality[1]);
	compute_tet_volume_quality(C, D, B, G, subvolume[2],
		tet_quality[2]);
	compute_tet_volume_quality(D, A, C, H, subvolume[3],
		tet_quality[3]);
	compute_tet_volume_quality(E, H, F, A, subvolume[4],
		tet_quality[4]);
	compute_tet_volume_quality(F, E, G, B, subvolume[5],
		tet_quality[5]);
	compute_tet_volume_quality(G, F, H, C, subvolume[6],
		tet_quality[6]);
	compute_tet_volume_quality(H, G, E, D, subvolume[7],
		tet_quality[7]);
	//for (unsigned int i = 0; i < 7; i++) {
	//	if (subvolume[i] * subvolume[i + 1] <= 0) {
	//		return false;
	//	}
	//}
	//if (subvolume[7] * subvolume[0] <= 0) {
	//	return false;
	//}
	double quality_factor[3] = { 1.609, 1.365, 1.190 };
	quality = (*std::min_element(&planar_quality[0], &planar_quality[6]))
		* (*std::min_element(&tet_quality[0], &tet_quality[8]));
	//if (quality * quality_factor[m_qualitytype]< m_quality_bound)
	//{
	//	return false;
	//}
	return true;
}

bool Hex::check_spatial_quad_quality(const ig::CVec<double, 3>& v0, const ig::CVec<double, 3>& v1, const ig::CVec<double, 3>& v2, const ig::CVec<double, 3>& v3, double& quality)
{
	ig::CVec<double, 3> V[4];
	V[0] = (v1 - v0).UnitCross(v2 - v0);
	V[1] = (v1 - v0).UnitCross(v3 - v0);
	V[2] = (v2 - v0).UnitCross(v3 - v0);
	V[3] = (v2 - v1).UnitCross(v3 - v1);
	double n[6];
	unsigned int count = 0;
	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = i + 1; j < 4; j++)
		{
			n[count] = V[i].Dot(V[j]);
			//if (n[count] <= 0)
			//	return false;
			count++;
		}
	}
	/*quality =  1.0 -  2.0 / Math<Real>::PI
		* Math<Real>::ACos(*std::min_element(&n[0], &n[6]));*/
	quality = 1.0 - 2.0 / M_PI * acos(*std::min_element(&n[0], &n[6]));
	return true;
}

void Hex::compute_tet_volume_quality(const ig::CVec<double, 3>& v0, const ig::CVec<double, 3>& v1, const ig::CVec<double, 3>& v2, const ig::CVec<double, 3>& v3, double& subvolume, double& quality)
{
	int m_qualitytype = 1;
	const ig::CVec<double, 3> v10 = v1 - v0;
	const ig::CVec<double, 3> v20 = v2 - v0;
	const ig::CVec<double, 3> v30 = v3 - v0;
	const ig::CVec<double, 3> v12 = v1 - v2;
	const ig::CVec<double, 3> v13 = v1 - v3;
	const ig::CVec<double, 3> v23 = v2 - v3;
	subvolume = v30.Dot(v10.Cross(v20)) / 6.0;
	//Reference: "Shape Measures for Quadrilaterals, Pyramids, Wedges and Hexahedra", Barry Joe
	if (m_qualitytype == 0) {
		//(1) shape measure: solid angle
		double L10 = v10.Length(), L20 = v20.Length(), L30 = v30.Length(), L12 =
			v12.Length(), L13 = v13.Length(), L23 = v23.Length();
		double theta[4];
		double pos_volume = abs(subvolume);
		theta[0] = 12 * pos_volume / sqrt((L10 + L20 + L12) * (L10 + L20 - L12)
			* (L10 + L30 + L13) * (L10 + L30 - L13) * (L20 + L30 + L23)
			* (L20 + L30 - L23));
		theta[1] = 12 * pos_volume / sqrt((L10 + L12 + L20) * (L10 + L12 - L20)
			* (L10 + L13 + L30) * (L10 + L13 - L30) * (L12 + L13 + L23)
			* (L12 + L13 - L23));
		theta[2] = 12 * pos_volume / sqrt((L20 + L12 + L10) * (L20 + L12 - L10)
			* (L12 + L23 + L13) * (L12 + L23 - L13) * (L23 + L20 + L30)
			* (L23 + L20 - L30));
		theta[3] = 12 * pos_volume / sqrt((L30 + L13 + L10) * (L30 + L13 - L10)
			* (L13 + L23 + L12) * (L13 + L23 - L12) * (L23 + L30 + L20)
			* (L23 + L30 - L20));
		quality = 9.0 / sqrt(6) * *std::min_element(&theta[0], &theta[4]);
	}
	else if (m_qualitytype == 1) {
		//(2) shape measure: \rho = 3 r_{in} \ r_{circ}
		double L10 = v10.Length(), L20 = v20.Length(), L30 = v30.Length(), L12 =
			v12.Length(), L13 = v13.Length(), L23 = v23.Length();
		double a = L10 * L23, b = L20 * L13, c = L30 * L12;
		double s0 =  0.5 * (v12.Cross(v13)).Length();
		double s1 =  0.5 * (v20.Cross(v23)).Length();
		double s2 =  0.5 * (v10.Cross(v13)).Length();
		double s3 =  0.5 * (v10.Cross(v12)).Length();
		quality = 216 * subvolume * subvolume / ((s0 + s1 + s2 + s3) * sqrt((a
			+ b + c) * (a + b - c) * (a + c - b) * (b + c - a)));
	}
	else if (m_qualitytype == 2) {
		//(3) shape measure \eta
		quality = 12 * pow(abs(3 * subvolume),
			(2.0 / 3.0)) / (v10.SquaredLength()
				+ v20.SquaredLength() + v30.SquaredLength()
				+ v12.SquaredLength() + v13.SquaredLength()
				+ v23.SquaredLength());
	}
}


double Hex::compute_cell_volume(const CVec<double, 3> &v1, const CVec<double, 3> &v2, const CVec<double, 3> &v3, const CVec<double, 3> &v4)
{
	double v;
	CVec<double, 3> vec1, vec2, vec3, temp_vec;
	vec1 = v2 - v1;
	vec2 = v3 - v1;
	vec3 = v4 - v1;

	temp_vec = vec1.Cross(vec2);
	v = temp_vec.Dot(vec3);
	return v / 6.0;

}