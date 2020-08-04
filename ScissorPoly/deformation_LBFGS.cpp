#pragma warning( disable : 4477 4018 4267 4244 4838)
#include "deformation_LBFGS.h"
#include <omp.h>
#define EPSL 1e-12
#define BOUND 100
#define LEN_LB 1e-8
using std::vector;
void evalfunc_de_LBFGS(int N, double* x, double *prev_x, double* f, double* g, void* user_pointer)
{
	deform_pointer *dp = (deform_pointer *)(user_pointer);
	const vector<double> &dpx = dp->pdi->dpx;
	const vector<double> &dpy = dp->pdi->dpy;
	const vector<double> &dpz = dp->pdi->dpz;
	const vector<int> &boundary_verts = dp->pdi->boundary_verts;
	const vector<OpenVolumeMesh::Geometry::Vec3i> &bfv_id = dp->pdi->bfv_id;
	const std::vector<std::vector<int>> &vertex_cell = dp->pdi->vertex_cell;
	const std::vector<std::vector<std::vector<int>>> &vertex_cell_vertex = dp->pdi->vertex_cell_vertex;
	const std::vector< OpenVolumeMesh::Geometry::Vec3d > &target_bfn = dp->pdi->target_bfn;
	const std::vector< std::vector<int> > &bvv_id = dp->pdi->bvv_id;
	const std::vector< std::vector<int> > &bvf_id = dp->pdi->bvf_id;
	const std::vector<int> &a2b = dp->pdi->a2b;
	const std::vector<std::vector<int>> &cell_vertex = dp->pdi->cell_vertex; //not using right order version, cause inverse is not in right order
	const std::vector< std::vector<std::vector<double> > > &vcv_S = dp->pdi->vcv_S;
	const std::vector< Eigen::Matrix3d > &cell_S = dp->pdi->cell_S;
	const std::vector<double> &lambda_array = dp->pdi->lambda_array;
	int boundary_face_number = dp->pdi->boundary_face_number;
	int boundary_vertface_number = dp->pdi->boundary_vertface_number;
	const vector<vector<int>> &bef_id = dp->pdi->bef_id;
	//feature edge part
	const std::vector<bool> &feature_edge_flag = dp->pdi->feature_edge_flag;
	const std::vector<std::vector<int>> &feature_v2e = dp->pdi->feature_v2e;
	const std::vector<std::vector<std::pair<int, int>>> &feature_neighbor_vert2cellpair = dp->pdi->feature_neighbor_vert2cellpair; //first face containging vert, second face the other direction
	//std::vector<std::vector<int>> feature_faces;
	const std::vector<OpenVolumeMesh::Geometry::Vec3d> &target_fea_n = dp->pdi->target_fea_n;
	//int n_feature_edges;
	const std::vector<int> &feature_edges = dp->pdi->feature_edge_array;
	const std::vector<std::pair<int, int>> &feature_e2v = dp->pdi->feature_e2v; //consistent feature edges
	const std::vector<std::vector<int>> &feature_neighbor = dp->pdi->feature_neighbor; //neighbor of edges at most 2 neighbors, cut by valence 3 vertex
	//const std::vector<OpenVolumeMesh::Geometry::Vec3d> &target_feature_edge_normal = dp->pdi->target_feature_edge_normal;
	int v_id_1, v_id_2, v_id_3, k;
	double D00, D10, D20, D01, D11, D21, D02, D12, D22;
	double C00, C01, C02, C10, C11, C12, C20, C21, C22;
	double A00, A01, A02, A10, A11, A12, A20, A21, A22;
	double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
	double AF, AF_05, I_AF_05, AF2, AF_I, AF_I_05, I_AF_I_05, exp_k;
	double D_AF_x, D_AF_y, D_AF_z, D_AF2_x, D_AF2_y, D_AF2_z, D_AF_I_x, D_AF_I_y, D_AF_I_z, D_det_A_x, D_det_A_y, D_det_A_z, inv_det_A_2_05;
	double d_A00_x, d_A10_y, d_A20_z, d_A01_x, d_A11_y, d_A21_z, d_A02_x, d_A12_y, d_A22_z;
	double dvex, dvey, dvez, tmp_g, dgx, dgy, dgz, dex, dey, dez, e, mips_e, volume_e, d_mips_e_x, d_mips_e_y, d_mips_e_z;
	double min_radius = 1e30;
	double alpha = 0.5; double beta = 0.5;
	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	double len;
	int bv_size = boundary_verts.size();
	assert(bv_size > 0);
	*f = 0.0;
	for (size_t i = 0; i < 3 * bv_size; i++)
	{
		g[i] = 0.0;
	}
	//consider all cells with at least one element on the boundary
	std::vector<int> cell_flag(bfv_id.size(), -1);
	double coeff_iso = 1.0 / (1.0 * boundary_vertface_number);
	//isometric part
	bv_size = 0;
	for (size_t i = 0; i < bv_size; i++)
	{
		int bid = boundary_verts[i];
		assert(bvf_id[bid].size() > 0);
		assert(vertex_cell[bid].size() == vertex_cell_vertex[bid].size());
		int vc_size = vertex_cell[bid].size();
		vector<double> local_cell_energy(vc_size, 0.0);
		double x0 = x[3 * i], y0 = x[3 * i + 1], z0 = x[3 * i + 2];
		//isometric part
		vector<double> exp_vec(vc_size, 0.0), gx_vec(vc_size, 0.0), gy_vec(vc_size, 0.0), gz_vec(vc_size, 0.0);
		for (size_t j = 0; j < vc_size; j++)
		{
			//decide local_cell_energy and gradient here
			const int* vv_id = vertex_cell_vertex[bid][j].data();
			const double* s_data = vcv_S[bid][j].data();
			k = 3 * j;
			v_id_1 = vv_id[0];
			v_id_2 = vv_id[1];
			v_id_3 = vv_id[2];
			x1 = dpx[v_id_1]; x2 = dpx[v_id_2]; x3 = dpx[v_id_3];
			y1 = dpy[v_id_1]; y2 = dpy[v_id_2]; y3 = dpy[v_id_3];
			z1 = dpz[v_id_1]; z2 = dpz[v_id_2]; z3 = dpz[v_id_3];
			//redefine on boundary needed
			if (a2b[v_id_1] != -1)
			{
				//vid1 on the boundary
				int tmp_bid = a2b[v_id_1];
				x1 = x[3 * tmp_bid];
				y1 = x[3 * tmp_bid + 1];
				z1 = x[3 * tmp_bid + 2];
			}
			if (a2b[v_id_2] != -1)
			{
				//vid1 on the boundary
				int tmp_bid = a2b[v_id_2];
				x2 = x[3 * tmp_bid];
				y2 = x[3 * tmp_bid + 1];
				z2 = x[3 * tmp_bid + 2];
			}
			if (a2b[v_id_3] != -1)
			{
				//vid1 on the boundary
				int tmp_bid = a2b[v_id_3];
				x3 = x[3 * tmp_bid];
				y3 = x[3 * tmp_bid + 1];
				z3 = x[3 * tmp_bid + 2];
			}
			/*p_vc_pos_x[k + 0] = x1; p_vc_pos_x[k + 1] = x2; p_vc_pos_x[k + 2] = x3;
			p_vc_pos_y[k + 0] = y1; p_vc_pos_y[k + 1] = y2; p_vc_pos_y[k + 2] = y3;
			p_vc_pos_z[k + 0] = z1; p_vc_pos_z[k + 1] = z2; p_vc_pos_z[k + 2] = z3;
			p_vc_n_cross_x[i] = (y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1);
			p_vc_n_cross_y[i] = (x3 - x1)*(z2 - z1) - (x2 - x1)*(z3 - z1);
			p_vc_n_cross_z[i] = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);*/
			D00 = x1 - x0; D10 = y1 - y0; D20 = z1 - z0;
			D01 = x2 - x0; D11 = y2 - y0; D21 = z2 - z0;
			D02 = x3 - x0; D12 = y3 - y0; D22 = z3 - z0;
			len = D00 * D00 + D10 * D10 + D20 * D20; if (len < min_radius) min_radius = len;
			len = D01 * D01 + D11 * D11 + D21 * D21; if (len < min_radius) min_radius = len;
			len = D02 * D02 + D12 * D12 + D22 * D22; if (len < min_radius) min_radius = len;
			C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
			C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
			C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];
			A00 = C00 * D00 + C10 * D01 + C20 * D02; d_A00_x = -C00 - C10 - C20;
			A10 = C00 * D10 + C10 * D11 + C20 * D12; d_A10_y = -C00 - C10 - C20;
			A20 = C00 * D20 + C10 * D21 + C20 * D22; d_A20_z = -C00 - C10 - C20;
			A01 = C01 * D00 + C11 * D01 + C21 * D02; d_A01_x = -C01 - C11 - C21;
			A11 = C01 * D10 + C11 * D11 + C21 * D12; d_A11_y = -C01 - C11 - C21;
			A21 = C01 * D20 + C11 * D21 + C21 * D22; d_A21_z = -C01 - C11 - C21;
			A02 = C02 * D00 + C12 * D01 + C22 * D02; d_A02_x = -C02 - C12 - C22;
			A12 = C02 * D10 + C12 * D11 + C22 * D12; d_A12_y = -C02 - C12 - C22;
			A22 = C02 * D20 + C12 * D21 + C22 * D22; d_A22_z = -C02 - C12 - C22;
			A2_00 = A00 * A00 + A10 * A10 + A20 * A20;
			A2_01 = A00 * A01 + A10 * A11 + A20 * A21;
			A2_02 = A00 * A02 + A10 * A12 + A20 * A22;
			A2_11 = A01 * A01 + A11 * A11 + A21 * A21;
			A2_12 = A01 * A02 + A11 * A12 + A21 * A22;
			A2_22 = A02 * A02 + A12 * A12 + A22 * A22;
			AF = A2_00 + A2_11 + A2_22; AF_05 = std::sqrt(AF); I_AF_05 = 0.5 / AF_05;
			AF2 = A2_00 * A2_00 + A2_01 * A2_01 * 2 + A2_02 * A2_02 * 2 + A2_11 * A2_11 + A2_12 * A2_12 * 2 + A2_22 * A2_22;
			AF_I = (AF*AF - AF2)*0.5; AF_I_05 = std::sqrt(AF_I); I_AF_I_05 = 0.5 / AF_I_05;
			det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
			i_det_A = 1.0 / det_A;
			D_AF_x = (2.0*A00*d_A00_x + 2.0*A01*d_A01_x + 2.0*A02*d_A02_x);
			D_AF_y = (2.0*A10*d_A10_y + 2.0*A11*d_A11_y + 2.0*A12*d_A12_y);
			D_AF_z = (2.0*A20*d_A20_z + 2.0*A21*d_A21_z + 2.0*A22*d_A22_z);
			D_AF2_x = (2.0*A2_00*2.0*A00*d_A00_x + 4.0*A2_01*(A01*d_A00_x + A00 * d_A01_x) + 4.0*A2_02*(A02*d_A00_x + A00 * d_A02_x)
				+ 2.0*A2_11*2.0*A01*d_A01_x + 2.0*A2_22*2.0*A02*d_A02_x + 4.0*A2_12*(A01*d_A02_x + A02 * d_A01_x));
			D_AF2_y = (2.0*A2_00*2.0*A10*d_A10_y + 4.0*A2_01*(A11*d_A10_y + A10 * d_A11_y) + 4.0*A2_02*(A12*d_A10_y + A10 * d_A12_y)
				+ 2.0*A2_11*2.0*A11*d_A11_y + 2.0*A2_22*2.0*A12*d_A12_y + 4.0*A2_12*(A11*d_A12_y + A12 * d_A11_y));
			D_AF2_z = (2.0*A2_00*2.0*A20*d_A20_z + 4.0*A2_01*(A21*d_A20_z + A20 * d_A21_z) + 4.0*A2_02*(A22*d_A20_z + A20 * d_A22_z)
				+ 2.0*A2_11*2.0*A21*d_A21_z + 2.0*A2_22*2.0*A22*d_A22_z + 4.0*A2_12*(A21*d_A22_z + A22 * d_A21_z));
			D_AF_I_x = (AF*D_AF_x - 0.5*D_AF2_x)*I_AF_I_05; D_AF_I_y = (AF*D_AF_y - 0.5*D_AF2_y)*I_AF_I_05; D_AF_I_z = (AF*D_AF_z - 0.5*D_AF2_z)*I_AF_I_05;
			D_AF_x *= I_AF_05; D_AF_y *= I_AF_05; D_AF_z *= I_AF_05;
			D_det_A_x = d_A00_x * A11 * A22 + d_A01_x * A12 * A20 + d_A02_x * A10 * A21 - d_A02_x * A11 * A20 - d_A01_x * A10 * A22 - d_A00_x * A12 * A21;
			D_det_A_y = A00 * d_A11_y * A22 + A01 * d_A12_y * A20 + A02 * d_A10_y * A21 - A02 * d_A11_y * A20 - A01 * d_A10_y * A22 - A00 * d_A12_y * A21;
			D_det_A_z = A00 * A11 * d_A22_z + A01 * A12 * d_A20_z + A02 * A10 * d_A21_z - A02 * A11 * d_A20_z - A01 * A10 * d_A22_z - A00 * A12 * d_A21_z;
			inv_det_A_2_05 = 0.5 *(1.0 - i_det_A * i_det_A);
			dvex = D_det_A_x * inv_det_A_2_05;
			dvey = D_det_A_y * inv_det_A_2_05;
			dvez = D_det_A_z * inv_det_A_2_05;
			tmp_g = AF_05 * AF_I_05;
			dgx = D_AF_x * AF_I_05 + AF_05 * D_AF_I_x;
			dgy = D_AF_y * AF_I_05 + AF_05 * D_AF_I_y;
			dgz = D_AF_z * AF_I_05 + AF_05 * D_AF_I_z;
			dex = (dgx *i_det_A - tmp_g * D_det_A_x *i_det_A*i_det_A);
			dey = (dgy *i_det_A - tmp_g * D_det_A_y *i_det_A*i_det_A);
			dez = (dgz *i_det_A - tmp_g * D_det_A_z *i_det_A*i_det_A);
			e = tmp_g * i_det_A;
			mips_e = (e*e - 1.0)*0.125;
			//mips_e = e/3.0;
			volume_e = 0.5*(det_A + i_det_A);
			exp_k = (alpha*mips_e + beta * volume_e);
			if (exp_k > 60) exp_k = 60;
			exp_vec[j] = exp_k;
			/*exp_e = std::exp(exp_k);
			local_energy += exp_e;*/
			d_mips_e_x = e * 0.25 * dex;
			d_mips_e_y = e * 0.25 * dey;
			d_mips_e_z = e * 0.25 * dez;
			/*d_mips_e_x = dex/3.0;
			d_mips_e_y = dey/3.0;
			d_mips_e_z = dez/3.0;*/
			gx_vec[j] = alpha * d_mips_e_x + beta * dvex;
			gy_vec[j] = alpha * d_mips_e_y + beta * dvey;
			gz_vec[j] = alpha * d_mips_e_z + beta * dvez;
#if 0
//debug below
			double diff[3];
			double tmp_f = 0.0;
			//x
			//std::cout << "x: " << std::endl;
			//std::cout << gx_vec[j] << std::endl;
			
			x0 = x[3 * i] + EPSL;
			
			D00 = x1 - x0; D10 = y1 - y0; D20 = z1 - z0;
			D01 = x2 - x0; D11 = y2 - y0; D21 = z2 - z0;
			D02 = x3 - x0; D12 = y3 - y0; D22 = z3 - z0;
			len = D00 * D00 + D10 * D10 + D20 * D20; if (len < min_radius) min_radius = len;
			len = D01 * D01 + D11 * D11 + D21 * D21; if (len < min_radius) min_radius = len;
			len = D02 * D02 + D12 * D12 + D22 * D22; if (len < min_radius) min_radius = len;
			C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
			C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
			C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];
			A00 = C00 * D00 + C10 * D01 + C20 * D02; d_A00_x = -C00 - C10 - C20;
			A10 = C00 * D10 + C10 * D11 + C20 * D12; d_A10_y = -C00 - C10 - C20;
			A20 = C00 * D20 + C10 * D21 + C20 * D22; d_A20_z = -C00 - C10 - C20;
			A01 = C01 * D00 + C11 * D01 + C21 * D02; d_A01_x = -C01 - C11 - C21;
			A11 = C01 * D10 + C11 * D11 + C21 * D12; d_A11_y = -C01 - C11 - C21;
			A21 = C01 * D20 + C11 * D21 + C21 * D22; d_A21_z = -C01 - C11 - C21;
			A02 = C02 * D00 + C12 * D01 + C22 * D02; d_A02_x = -C02 - C12 - C22;
			A12 = C02 * D10 + C12 * D11 + C22 * D12; d_A12_y = -C02 - C12 - C22;
			A22 = C02 * D20 + C12 * D21 + C22 * D22; d_A22_z = -C02 - C12 - C22;
			A2_00 = A00 * A00 + A10 * A10 + A20 * A20;
			A2_01 = A00 * A01 + A10 * A11 + A20 * A21;
			A2_02 = A00 * A02 + A10 * A12 + A20 * A22;
			A2_11 = A01 * A01 + A11 * A11 + A21 * A21;
			A2_12 = A01 * A02 + A11 * A12 + A21 * A22;
			A2_22 = A02 * A02 + A12 * A12 + A22 * A22;
			AF = A2_00 + A2_11 + A2_22; AF_05 = std::sqrt(AF); I_AF_05 = 0.5 / AF_05;
			AF2 = A2_00 * A2_00 + A2_01 * A2_01 * 2 + A2_02 * A2_02 * 2 + A2_11 * A2_11 + A2_12 * A2_12 * 2 + A2_22 * A2_22;
			AF_I = (AF*AF - AF2)*0.5; AF_I_05 = std::sqrt(AF_I); I_AF_I_05 = 0.5 / AF_I_05;
			det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
			i_det_A = 1.0 / det_A;
			D_AF_x = (2.0*A00*d_A00_x + 2.0*A01*d_A01_x + 2.0*A02*d_A02_x);
			D_AF_y = (2.0*A10*d_A10_y + 2.0*A11*d_A11_y + 2.0*A12*d_A12_y);
			D_AF_z = (2.0*A20*d_A20_z + 2.0*A21*d_A21_z + 2.0*A22*d_A22_z);
			D_AF2_x = (2.0*A2_00*2.0*A00*d_A00_x + 4.0*A2_01*(A01*d_A00_x + A00 * d_A01_x) + 4.0*A2_02*(A02*d_A00_x + A00 * d_A02_x)
				+ 2.0*A2_11*2.0*A01*d_A01_x + 2.0*A2_22*2.0*A02*d_A02_x + 4.0*A2_12*(A01*d_A02_x + A02 * d_A01_x));
			D_AF2_y = (2.0*A2_00*2.0*A10*d_A10_y + 4.0*A2_01*(A11*d_A10_y + A10 * d_A11_y) + 4.0*A2_02*(A12*d_A10_y + A10 * d_A12_y)
				+ 2.0*A2_11*2.0*A11*d_A11_y + 2.0*A2_22*2.0*A12*d_A12_y + 4.0*A2_12*(A11*d_A12_y + A12 * d_A11_y));
			D_AF2_z = (2.0*A2_00*2.0*A20*d_A20_z + 4.0*A2_01*(A21*d_A20_z + A20 * d_A21_z) + 4.0*A2_02*(A22*d_A20_z + A20 * d_A22_z)
				+ 2.0*A2_11*2.0*A21*d_A21_z + 2.0*A2_22*2.0*A22*d_A22_z + 4.0*A2_12*(A21*d_A22_z + A22 * d_A21_z));
			D_AF_I_x = (AF*D_AF_x - 0.5*D_AF2_x)*I_AF_I_05; D_AF_I_y = (AF*D_AF_y - 0.5*D_AF2_y)*I_AF_I_05; D_AF_I_z = (AF*D_AF_z - 0.5*D_AF2_z)*I_AF_I_05;
			D_AF_x *= I_AF_05; D_AF_y *= I_AF_05; D_AF_z *= I_AF_05;
			D_det_A_x = d_A00_x * A11 * A22 + d_A01_x * A12 * A20 + d_A02_x * A10 * A21 - d_A02_x * A11 * A20 - d_A01_x * A10 * A22 - d_A00_x * A12 * A21;
			D_det_A_y = A00 * d_A11_y * A22 + A01 * d_A12_y * A20 + A02 * d_A10_y * A21 - A02 * d_A11_y * A20 - A01 * d_A10_y * A22 - A00 * d_A12_y * A21;
			D_det_A_z = A00 * A11 * d_A22_z + A01 * A12 * d_A20_z + A02 * A10 * d_A21_z - A02 * A11 * d_A20_z - A01 * A10 * d_A22_z - A00 * A12 * d_A21_z;
			inv_det_A_2_05 = 0.5 *(1.0 - i_det_A * i_det_A);
			dvex = D_det_A_x * inv_det_A_2_05;
			dvey = D_det_A_y * inv_det_A_2_05;
			dvez = D_det_A_z * inv_det_A_2_05;
			tmp_g = AF_05 * AF_I_05;
			dgx = D_AF_x * AF_I_05 + AF_05 * D_AF_I_x;
			dgy = D_AF_y * AF_I_05 + AF_05 * D_AF_I_y;
			dgz = D_AF_z * AF_I_05 + AF_05 * D_AF_I_z;
			dex = (dgx *i_det_A - tmp_g * D_det_A_x *i_det_A*i_det_A);
			dey = (dgy *i_det_A - tmp_g * D_det_A_y *i_det_A*i_det_A);
			dez = (dgz *i_det_A - tmp_g * D_det_A_z *i_det_A*i_det_A);
			e = tmp_g * i_det_A;
			mips_e = (e*e - 1.0)*0.125;
			//mips_e = e/3.0;
			volume_e = 0.5*(det_A + i_det_A);
			exp_k = (alpha*mips_e + beta * volume_e);
			if (exp_k > 60) exp_k = 60;
			tmp_f = exp_k;
			//std::cout << (tmp_f - exp_vec[j]) / EPSL << std::endl;
			diff[0] = (tmp_f - exp_vec[j]) / EPSL;
			//y
			//std::cout << "y: " << std::endl;
			//std::cout << gy_vec[j] << std::endl;
			tmp_f = 0.0;
			x0 = x[3 * i];
			y0 = x[3 * i + 1] + EPSL;
			D00 = x1 - x0; D10 = y1 - y0; D20 = z1 - z0;
			D01 = x2 - x0; D11 = y2 - y0; D21 = z2 - z0;
			D02 = x3 - x0; D12 = y3 - y0; D22 = z3 - z0;
			len = D00 * D00 + D10 * D10 + D20 * D20; if (len < min_radius) min_radius = len;
			len = D01 * D01 + D11 * D11 + D21 * D21; if (len < min_radius) min_radius = len;
			len = D02 * D02 + D12 * D12 + D22 * D22; if (len < min_radius) min_radius = len;
			C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
			C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
			C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];
			A00 = C00 * D00 + C10 * D01 + C20 * D02; d_A00_x = -C00 - C10 - C20;
			A10 = C00 * D10 + C10 * D11 + C20 * D12; d_A10_y = -C00 - C10 - C20;
			A20 = C00 * D20 + C10 * D21 + C20 * D22; d_A20_z = -C00 - C10 - C20;
			A01 = C01 * D00 + C11 * D01 + C21 * D02; d_A01_x = -C01 - C11 - C21;
			A11 = C01 * D10 + C11 * D11 + C21 * D12; d_A11_y = -C01 - C11 - C21;
			A21 = C01 * D20 + C11 * D21 + C21 * D22; d_A21_z = -C01 - C11 - C21;
			A02 = C02 * D00 + C12 * D01 + C22 * D02; d_A02_x = -C02 - C12 - C22;
			A12 = C02 * D10 + C12 * D11 + C22 * D12; d_A12_y = -C02 - C12 - C22;
			A22 = C02 * D20 + C12 * D21 + C22 * D22; d_A22_z = -C02 - C12 - C22;
			A2_00 = A00 * A00 + A10 * A10 + A20 * A20;
			A2_01 = A00 * A01 + A10 * A11 + A20 * A21;
			A2_02 = A00 * A02 + A10 * A12 + A20 * A22;
			A2_11 = A01 * A01 + A11 * A11 + A21 * A21;
			A2_12 = A01 * A02 + A11 * A12 + A21 * A22;
			A2_22 = A02 * A02 + A12 * A12 + A22 * A22;
			AF = A2_00 + A2_11 + A2_22; AF_05 = std::sqrt(AF); I_AF_05 = 0.5 / AF_05;
			AF2 = A2_00 * A2_00 + A2_01 * A2_01 * 2 + A2_02 * A2_02 * 2 + A2_11 * A2_11 + A2_12 * A2_12 * 2 + A2_22 * A2_22;
			AF_I = (AF*AF - AF2)*0.5; AF_I_05 = std::sqrt(AF_I); I_AF_I_05 = 0.5 / AF_I_05;
			det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
			i_det_A = 1.0 / det_A;
			D_AF_x = (2.0*A00*d_A00_x + 2.0*A01*d_A01_x + 2.0*A02*d_A02_x);
			D_AF_y = (2.0*A10*d_A10_y + 2.0*A11*d_A11_y + 2.0*A12*d_A12_y);
			D_AF_z = (2.0*A20*d_A20_z + 2.0*A21*d_A21_z + 2.0*A22*d_A22_z);
			D_AF2_x = (2.0*A2_00*2.0*A00*d_A00_x + 4.0*A2_01*(A01*d_A00_x + A00 * d_A01_x) + 4.0*A2_02*(A02*d_A00_x + A00 * d_A02_x)
				+ 2.0*A2_11*2.0*A01*d_A01_x + 2.0*A2_22*2.0*A02*d_A02_x + 4.0*A2_12*(A01*d_A02_x + A02 * d_A01_x));
			D_AF2_y = (2.0*A2_00*2.0*A10*d_A10_y + 4.0*A2_01*(A11*d_A10_y + A10 * d_A11_y) + 4.0*A2_02*(A12*d_A10_y + A10 * d_A12_y)
				+ 2.0*A2_11*2.0*A11*d_A11_y + 2.0*A2_22*2.0*A12*d_A12_y + 4.0*A2_12*(A11*d_A12_y + A12 * d_A11_y));
			D_AF2_z = (2.0*A2_00*2.0*A20*d_A20_z + 4.0*A2_01*(A21*d_A20_z + A20 * d_A21_z) + 4.0*A2_02*(A22*d_A20_z + A20 * d_A22_z)
				+ 2.0*A2_11*2.0*A21*d_A21_z + 2.0*A2_22*2.0*A22*d_A22_z + 4.0*A2_12*(A21*d_A22_z + A22 * d_A21_z));
			D_AF_I_x = (AF*D_AF_x - 0.5*D_AF2_x)*I_AF_I_05; D_AF_I_y = (AF*D_AF_y - 0.5*D_AF2_y)*I_AF_I_05; D_AF_I_z = (AF*D_AF_z - 0.5*D_AF2_z)*I_AF_I_05;
			D_AF_x *= I_AF_05; D_AF_y *= I_AF_05; D_AF_z *= I_AF_05;
			D_det_A_x = d_A00_x * A11 * A22 + d_A01_x * A12 * A20 + d_A02_x * A10 * A21 - d_A02_x * A11 * A20 - d_A01_x * A10 * A22 - d_A00_x * A12 * A21;
			D_det_A_y = A00 * d_A11_y * A22 + A01 * d_A12_y * A20 + A02 * d_A10_y * A21 - A02 * d_A11_y * A20 - A01 * d_A10_y * A22 - A00 * d_A12_y * A21;
			D_det_A_z = A00 * A11 * d_A22_z + A01 * A12 * d_A20_z + A02 * A10 * d_A21_z - A02 * A11 * d_A20_z - A01 * A10 * d_A22_z - A00 * A12 * d_A21_z;
			inv_det_A_2_05 = 0.5 *(1.0 - i_det_A * i_det_A);
			dvex = D_det_A_x * inv_det_A_2_05;
			dvey = D_det_A_y * inv_det_A_2_05;
			dvez = D_det_A_z * inv_det_A_2_05;
			tmp_g = AF_05 * AF_I_05;
			dgx = D_AF_x * AF_I_05 + AF_05 * D_AF_I_x;
			dgy = D_AF_y * AF_I_05 + AF_05 * D_AF_I_y;
			dgz = D_AF_z * AF_I_05 + AF_05 * D_AF_I_z;
			dex = (dgx *i_det_A - tmp_g * D_det_A_x *i_det_A*i_det_A);
			dey = (dgy *i_det_A - tmp_g * D_det_A_y *i_det_A*i_det_A);
			dez = (dgz *i_det_A - tmp_g * D_det_A_z *i_det_A*i_det_A);
			e = tmp_g * i_det_A;
			mips_e = (e*e - 1.0)*0.125;
			//mips_e = e/3.0;
			volume_e = 0.5*(det_A + i_det_A);
			exp_k = (alpha*mips_e + beta * volume_e);
			if (exp_k > 60) exp_k = 60;
			tmp_f = exp_k;
			//std::cout << (tmp_f - exp_vec[j]) / EPSL << std::endl;
			diff[1] = (tmp_f - exp_vec[j]) / EPSL;
			
			//z
			//std::cout << "z: " << std::endl;
			//std::cout << gz_vec[j] << std::endl;
			tmp_f = 0.0;
			x0 = x[3 * i];
			y0 = x[3 * i + 1];
			z0 = x[3 * i + 2] + EPSL;
			D00 = x1 - x0; D10 = y1 - y0; D20 = z1 - z0;
			D01 = x2 - x0; D11 = y2 - y0; D21 = z2 - z0;
			D02 = x3 - x0; D12 = y3 - y0; D22 = z3 - z0;
			len = D00 * D00 + D10 * D10 + D20 * D20; if (len < min_radius) min_radius = len;
			len = D01 * D01 + D11 * D11 + D21 * D21; if (len < min_radius) min_radius = len;
			len = D02 * D02 + D12 * D12 + D22 * D22; if (len < min_radius) min_radius = len;
			C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
			C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
			C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];
			A00 = C00 * D00 + C10 * D01 + C20 * D02; d_A00_x = -C00 - C10 - C20;
			A10 = C00 * D10 + C10 * D11 + C20 * D12; d_A10_y = -C00 - C10 - C20;
			A20 = C00 * D20 + C10 * D21 + C20 * D22; d_A20_z = -C00 - C10 - C20;
			A01 = C01 * D00 + C11 * D01 + C21 * D02; d_A01_x = -C01 - C11 - C21;
			A11 = C01 * D10 + C11 * D11 + C21 * D12; d_A11_y = -C01 - C11 - C21;
			A21 = C01 * D20 + C11 * D21 + C21 * D22; d_A21_z = -C01 - C11 - C21;
			A02 = C02 * D00 + C12 * D01 + C22 * D02; d_A02_x = -C02 - C12 - C22;
			A12 = C02 * D10 + C12 * D11 + C22 * D12; d_A12_y = -C02 - C12 - C22;
			A22 = C02 * D20 + C12 * D21 + C22 * D22; d_A22_z = -C02 - C12 - C22;
			A2_00 = A00 * A00 + A10 * A10 + A20 * A20;
			A2_01 = A00 * A01 + A10 * A11 + A20 * A21;
			A2_02 = A00 * A02 + A10 * A12 + A20 * A22;
			A2_11 = A01 * A01 + A11 * A11 + A21 * A21;
			A2_12 = A01 * A02 + A11 * A12 + A21 * A22;
			A2_22 = A02 * A02 + A12 * A12 + A22 * A22;
			AF = A2_00 + A2_11 + A2_22; AF_05 = std::sqrt(AF); I_AF_05 = 0.5 / AF_05;
			AF2 = A2_00 * A2_00 + A2_01 * A2_01 * 2 + A2_02 * A2_02 * 2 + A2_11 * A2_11 + A2_12 * A2_12 * 2 + A2_22 * A2_22;
			AF_I = (AF*AF - AF2)*0.5; AF_I_05 = std::sqrt(AF_I); I_AF_I_05 = 0.5 / AF_I_05;
			det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
			i_det_A = 1.0 / det_A;
			D_AF_x = (2.0*A00*d_A00_x + 2.0*A01*d_A01_x + 2.0*A02*d_A02_x);
			D_AF_y = (2.0*A10*d_A10_y + 2.0*A11*d_A11_y + 2.0*A12*d_A12_y);
			D_AF_z = (2.0*A20*d_A20_z + 2.0*A21*d_A21_z + 2.0*A22*d_A22_z);
			D_AF2_x = (2.0*A2_00*2.0*A00*d_A00_x + 4.0*A2_01*(A01*d_A00_x + A00 * d_A01_x) + 4.0*A2_02*(A02*d_A00_x + A00 * d_A02_x)
				+ 2.0*A2_11*2.0*A01*d_A01_x + 2.0*A2_22*2.0*A02*d_A02_x + 4.0*A2_12*(A01*d_A02_x + A02 * d_A01_x));
			D_AF2_y = (2.0*A2_00*2.0*A10*d_A10_y + 4.0*A2_01*(A11*d_A10_y + A10 * d_A11_y) + 4.0*A2_02*(A12*d_A10_y + A10 * d_A12_y)
				+ 2.0*A2_11*2.0*A11*d_A11_y + 2.0*A2_22*2.0*A12*d_A12_y + 4.0*A2_12*(A11*d_A12_y + A12 * d_A11_y));
			D_AF2_z = (2.0*A2_00*2.0*A20*d_A20_z + 4.0*A2_01*(A21*d_A20_z + A20 * d_A21_z) + 4.0*A2_02*(A22*d_A20_z + A20 * d_A22_z)
				+ 2.0*A2_11*2.0*A21*d_A21_z + 2.0*A2_22*2.0*A22*d_A22_z + 4.0*A2_12*(A21*d_A22_z + A22 * d_A21_z));
			D_AF_I_x = (AF*D_AF_x - 0.5*D_AF2_x)*I_AF_I_05; D_AF_I_y = (AF*D_AF_y - 0.5*D_AF2_y)*I_AF_I_05; D_AF_I_z = (AF*D_AF_z - 0.5*D_AF2_z)*I_AF_I_05;
			D_AF_x *= I_AF_05; D_AF_y *= I_AF_05; D_AF_z *= I_AF_05;
			D_det_A_x = d_A00_x * A11 * A22 + d_A01_x * A12 * A20 + d_A02_x * A10 * A21 - d_A02_x * A11 * A20 - d_A01_x * A10 * A22 - d_A00_x * A12 * A21;
			D_det_A_y = A00 * d_A11_y * A22 + A01 * d_A12_y * A20 + A02 * d_A10_y * A21 - A02 * d_A11_y * A20 - A01 * d_A10_y * A22 - A00 * d_A12_y * A21;
			D_det_A_z = A00 * A11 * d_A22_z + A01 * A12 * d_A20_z + A02 * A10 * d_A21_z - A02 * A11 * d_A20_z - A01 * A10 * d_A22_z - A00 * A12 * d_A21_z;
			inv_det_A_2_05 = 0.5 *(1.0 - i_det_A * i_det_A);
			dvex = D_det_A_x * inv_det_A_2_05;
			dvey = D_det_A_y * inv_det_A_2_05;
			dvez = D_det_A_z * inv_det_A_2_05;
			tmp_g = AF_05 * AF_I_05;
			dgx = D_AF_x * AF_I_05 + AF_05 * D_AF_I_x;
			dgy = D_AF_y * AF_I_05 + AF_05 * D_AF_I_y;
			dgz = D_AF_z * AF_I_05 + AF_05 * D_AF_I_z;
			dex = (dgx *i_det_A - tmp_g * D_det_A_x *i_det_A*i_det_A);
			dey = (dgy *i_det_A - tmp_g * D_det_A_y *i_det_A*i_det_A);
			dez = (dgz *i_det_A - tmp_g * D_det_A_z *i_det_A*i_det_A);
			e = tmp_g * i_det_A;
			mips_e = (e*e - 1.0)*0.125;
			//mips_e = e/3.0;
			volume_e = 0.5*(det_A + i_det_A);
			exp_k = (alpha*mips_e + beta * volume_e);
			if (exp_k > 60) exp_k = 60;
			tmp_f = exp_k;
			//std::cout << (tmp_f - exp_vec[j]) / EPSL << std::endl;
			diff[2] = (tmp_f - exp_vec[j]) / EPSL;
			tmp_f = 0.0;
			z0 = x[3 * i + 2];
#endif
		}
		for (size_t j = 0; j < vc_size; j++)
		{
			double exp_e = std::exp(exp_vec[j]);
			if (cell_flag[vertex_cell[bid][j]] == -1)
			{
				//not marked yet
				cell_flag[vertex_cell[bid][j]] = 1;
				*f += coeff_iso * exp_e;
			}
			//derivative
			g[3 * i]	 += coeff_iso * gx_vec[j] * exp_e;
			g[3 * i + 1] += coeff_iso * gy_vec[j] * exp_e;
			g[3 * i + 2] += coeff_iso * gz_vec[j] * exp_e;
		}
	}
	cell_flag.clear();
	cell_flag.resize(bfv_id.size(), -1);
	//normal part
	double lambda = lambda_array[0] / (1.0 * boundary_face_number); //uniform weight used here
	bv_size = boundary_verts.size();
	for (size_t p = 0; p < bv_size; p++)
	{	
		int bid = boundary_verts[p];
		double mu = 1.0; int bvf_size = 0; int bvv_size = 0;
		//double* p_vf_pos_x = NULL; double* p_vf_pos_y = NULL; double* p_vf_pos_z = NULL;
		double normal_e = 0.0; double smooth_e = 0;
		double ne_gx = 0.0; double ne_gy = 0.0; double ne_gz = 0.0;
		double se_gx = 0.0; double se_gy = 0.0; double se_gz = 0.0;
		const std::vector<int>& one_vv_id = bvv_id[bid];
		const std::vector<int>& one_vf_id = bvf_id[bid];
		bvv_size = one_vv_id.size(); bvf_size = bvv_size;
		vector<double> normal_e_array(bvv_size, 0.0);
		double x0 = x[3 * p], y0 = x[3 * p + 1], z0 = x[3 * p + 2];
		for (int i = 0; i < bvv_size; ++i)
		{
			int j = (i + 1) % bvv_size; int k = (i + 2) % bvv_size;
			/*x1 = dpx[one_vv_id[i]]; y1 = dpy[one_vv_id[i]]; z1 = dpz[one_vv_id[i]];
			x2 = dpx[one_vv_id[j]]; y2 = dpy[one_vv_id[j]]; z2 = dpz[one_vv_id[j]];
			x3 = dpx[one_vv_id[k]]; y3 = dpy[one_vv_id[k]]; z3 = dpz[one_vv_id[k]];*/
			v_id_1 = one_vv_id[i];
			v_id_2 = one_vv_id[j];
			v_id_3 = one_vv_id[k];
			assert(a2b[v_id_1] != -1);
			assert(a2b[v_id_2] != -1);
			assert(a2b[v_id_3] != -1);
			if (a2b[v_id_1] != -1)
			{
				//vid1 on the boundary
				int tmp_bid = a2b[v_id_1];
				x1 = x[3 * tmp_bid];
				y1 = x[3 * tmp_bid + 1];
				z1 = x[3 * tmp_bid + 2];
			}
			if (a2b[v_id_2] != -1)
			{
				//vid1 on the boundary
				int tmp_bid = a2b[v_id_2];
				x2 = x[3 * tmp_bid];
				y2 = x[3 * tmp_bid + 1];
				z2 = x[3 * tmp_bid + 2];
			}
			if (a2b[v_id_3] != -1)
			{
				//vid1 on the boundary
				int tmp_bid = a2b[v_id_3];
				x3 = x[3 * tmp_bid];
				y3 = x[3 * tmp_bid + 1];
				z3 = x[3 * tmp_bid + 2];
			}
			/*p_vf_pos_x[4 * i + 0] = x1; p_vf_pos_x[4 * i + 1] = x2; p_vf_pos_x[4 * i + 2] = x3;
			p_vf_pos_y[4 * i + 0] = y1; p_vf_pos_y[4 * i + 1] = y2; p_vf_pos_y[4 * i + 2] = y3;
			p_vf_pos_z[4 * i + 0] = z1; p_vf_pos_z[4 * i + 1] = z2; p_vf_pos_z[4 * i + 2] = z3;*/
			const OpenVolumeMesh::Geometry::Vec3d& tn = target_bfn[one_vf_id[i]];
			//p_vf_pos_x[4 * i + 3] = tn[0]; p_vf_pos_y[4 * i + 3] = tn[1]; p_vf_pos_z[4 * i + 3] = tn[2];
			double nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
			double ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
			double nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);
			double len2 = nx * nx + ny * ny + nz * nz;  double len = std::sqrt(len2);
			if (len < LEN_LB)
			{
				len = LEN_LB;
				len2 = len * len;
			}
			
			//double E_ne = 1.0 + 1.0e-8 - (nx*tn[0] + ny*tn[1] + nz*tn[2])/len;//polycube energy
			double E_ne = 0.5*((nx / len - tn[0])*(nx / len - tn[0]) + (ny / len - tn[1])*(ny / len - tn[1]) + (nz / len - tn[2])*(nz / len - tn[2]));//polycube energy
			
			////local_energy_ratio[one_vf_id[i]];
			//normal_e += lambda * E_ne;
			normal_e_array[i] = lambda * E_ne;
			double d_nx_x0 = 0.0;     double d_nx_y0 = z1 - z2; double d_nx_z0 = y2 - y1;
			double d_ny_x0 = z2 - z1; double d_ny_y0 = 0.0;     double d_ny_z0 = x1 - x2;
			double d_nz_x0 = y1 - y2; double d_nz_y0 = x2 - x1; double d_nz_z0 = 0.0;
			double d_len2_x0 = (nx*d_nx_x0 + ny * d_ny_x0 + nz * d_nz_x0); //miss 2
			double d_len2_y0 = (nx*d_nx_y0 + ny * d_ny_y0 + nz * d_nz_y0);
			double d_len2_z0 = (nx*d_nx_z0 + ny * d_ny_z0 + nz * d_nz_z0);
			double d_len_x0 = d_len2_x0 / len; double d_len_y0 = d_len2_y0 / len; double d_len_z0 = d_len2_z0 / len;
			double d_nx_len_x0 = d_nx_x0 / len - nx * d_len_x0 / (len*len);
			double d_nx_len_y0 = d_nx_y0 / len - nx * d_len_y0 / (len*len);
			double d_nx_len_z0 = d_nx_z0 / len - nx * d_len_z0 / (len*len);
			double d_ny_len_x0 = d_ny_x0 / len - ny * d_len_x0 / (len*len);
			double d_ny_len_y0 = d_ny_y0 / len - ny * d_len_y0 / (len*len);
			double d_ny_len_z0 = d_ny_z0 / len - ny * d_len_z0 / (len*len);
			double d_nz_len_x0 = d_nz_x0 / len - nz * d_len_x0 / (len*len);
			double d_nz_len_y0 = d_nz_y0 / len - nz * d_len_y0 / (len*len);
			double d_nz_len_z0 = d_nz_z0 / len - nz * d_len_z0 / (len*len);
			ne_gx += lambda * (d_nx_len_x0*(nx / len - tn[0]) + d_ny_len_x0 * (ny / len - tn[1]) + d_nz_len_x0 * (nz / len - tn[2]));
			ne_gy += lambda * (d_nx_len_y0*(nx / len - tn[0]) + d_ny_len_y0 * (ny / len - tn[1]) + d_nz_len_y0 * (nz / len - tn[2]));
			ne_gz += lambda * (d_nx_len_z0*(nx / len - tn[0]) + d_ny_len_z0 * (ny / len - tn[1]) + d_nz_len_z0 * (nz / len - tn[2]));
			
#if 0
//debug part
			double tmpx, tmpy, tmpz;
			tmpx = lambda * (d_nx_len_x0*(nx / len - tn[0]) + d_ny_len_x0 * (ny / len - tn[1]) + d_nz_len_x0 * (nz / len - tn[2]));
			tmpy = lambda * (d_nx_len_y0*(nx / len - tn[0]) + d_ny_len_y0 * (ny / len - tn[1]) + d_nz_len_y0 * (nz / len - tn[2]));
			tmpz = lambda * (d_nx_len_z0*(nx / len - tn[0]) + d_ny_len_z0 * (ny / len - tn[1]) + d_nz_len_z0 * (nz / len - tn[2]));
			//x
			double tmp_f = 0.0;
			
			tmp_f = 0.0;
			x0 = x[3 * p] + EPSL, y0 = x[3 * p + 1], z0 = x[3 * p + 2];
			
			nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
			ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
			nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);
			len2 = nx * nx + ny * ny + nz * nz; len = std::sqrt(len2);
			if (len < LEN_LB)
			{
				len = LEN_LB;
				len2 = len * len;
			}
			//double E_ne = 1.0 + 1.0e-8 - (nx*tn[0] + ny*tn[1] + nz*tn[2])/len;//polycube energy
			E_ne = 0.5*((nx / len - tn[0])*(nx / len - tn[0]) + (ny / len - tn[1])*(ny / len - tn[1]) + (nz / len - tn[2])*(nz / len - tn[2]));//polycube energy
			////local_energy_ratio[one_vf_id[i]];
			//normal_e += lambda * E_ne;
			//normal_e_array[i] = lambda * E_ne;
			tmp_f = lambda * E_ne;
			
			std::cout << "x: " << std::endl;
			std::cout << tmpx << std::endl;
			std::cout << (tmp_f - normal_e_array[i]) / EPSL << std::endl;
			//y
			tmp_f = 0.0;
			x0 = x[3 * p], y0 = x[3 * p + 1] + EPSL, z0 = x[3 * p + 2];
			nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
			ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
			nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);
			len2 = nx * nx + ny * ny + nz * nz; len = std::sqrt(len2);
			if (len < LEN_LB)
			{
				len = LEN_LB;
				len2 = len * len;
			}
			//double E_ne = 1.0 + 1.0e-8 - (nx*tn[0] + ny*tn[1] + nz*tn[2])/len;//polycube energy
			E_ne = 0.5*((nx / len - tn[0])*(nx / len - tn[0]) + (ny / len - tn[1])*(ny / len - tn[1]) + (nz / len - tn[2])*(nz / len - tn[2]));//polycube energy
			////local_energy_ratio[one_vf_id[i]];
			//normal_e += lambda * E_ne;
			//normal_e_array[i] = lambda * E_ne;
			tmp_f = lambda * E_ne;
			std::cout << "y: " << std::endl;
			std::cout << tmpy << std::endl;
			std::cout << (tmp_f - normal_e_array[i]) / EPSL << std::endl;
			//z
			tmp_f = 0.0;
			std::cout << "z0: " << std::setprecision(16) << z0 << std::endl;
			x0 = x[3 * p]; 
			y0 = x[3 * p + 1];
			z0 = x[3 * p + 2] + EPSL;
			//z0 = x[3 * p + 2] + 0.0001;
			std::cout << "z0: " << std::setprecision(16) << z0 << std::endl;
			std::cout << "nx: " << std::setprecision(16) << nx << std::endl;
			std::cout << "ny: " << std::setprecision(16) << ny << std::endl;
			std::cout << "nz: " << std::setprecision(16) << nz << std::endl;
			nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
			ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
			nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);
			std::cout << "nx: " << std::setprecision(16) << nx << std::endl;
			std::cout << "ny: " << std::setprecision(16) << ny << std::endl;
			std::cout << "nz: " << std::setprecision(16) << nz << std::endl;
			std::cout << "len: " << std::setprecision(16)  << len << std::endl;
			len2 = nx * nx + ny * ny + nz * nz; len = std::sqrt(len2);
			if (len < LEN_LB)
			{
				len = LEN_LB;
				len2 = len * len;
			}
			std::cout << "len: " << std::setprecision(16) << len << std::endl;
			
			//double E_ne = 1.0 + 1.0e-8 - (nx*tn[0] + ny*tn[1] + nz*tn[2])/len;//polycube energy
			E_ne = 0.5*((nx / len - tn[0])*(nx / len - tn[0]) + (ny / len - tn[1])*(ny / len - tn[1]) + (nz / len - tn[2])*(nz / len - tn[2]));//polycube energy
			////local_energy_ratio[one_vf_id[i]];
			//normal_e += lambda * E_ne;
			//normal_e_array[i] = lambda * E_ne;
			tmp_f = lambda * E_ne;
			std::cout << "z: " << std::endl;
			std::cout << std::setprecision(16) <<  tmpz << std::endl;
			std::cout << std::setprecision(16) << (tmp_f - normal_e_array[i]) / EPSL << std::endl;
			tmp_f = 0.0;
			z0 = x[3 * p + 2];
			
#endif
		}
		for (size_t i = 0; i < bvf_size; i++)
		{
			if (cell_flag[bvf_id[p][i]] == -1)
			{
				cell_flag[bvf_id[p][i]] = 1;
				*f += normal_e_array[i];
			}
		}
		g[3 * p] += ne_gx;
		g[3 * p + 1] += ne_gy;
		g[3 * p + 2] += ne_gz;
		
	}
	
	//feature preserving
	double coeff_fea = lambda_array[0] / (1.0 * feature_edges.size());
	int feature_edge_size = feature_edges.size();
#if 0
	//original vertion of energy
	feature_edge_size = 0; //not using feature edge
	for (size_t p = 0; p < feature_edge_size; p++)
	{
		int eid = feature_edges[p];
		int bf[2];
		assert(bef_id[eid].size() == 2);
		bf[0] = bef_id[eid][0];
		bf[1] = bef_id[eid][1];
		int v0, v1, v2(-1), v3(-1);
		v0 = feature_e2v[eid].first;
		v1 = feature_e2v[eid].second;
		for (size_t i = 0; i < 3; i++)
		{
			//set v2
			int after = (i + 1) % 3;
			for (size_t j = 0; j < 2; j++)
			{
				if (bfv_id[bf[j]][i] == v0 && bfv_id[bf[j]][after] == v1)
				{
					int afterafter = (after + 1) % 3;
					v2 = bfv_id[bf[j]][afterafter];
				}
				if (bfv_id[bf[j]][i] == v1 && bfv_id[bf[j]][after] == v0)
				{
					int afterafter = (after + 1) % 3;
					v3 = bfv_id[bf[j]][afterafter];
				}
			}
			
		}
		assert(v2 != -1 && v3 != -1);
		int v0_id[2] = { v0, v0 };
		int v1_id[2] = { v3, v1 };
		int v2_id[2] = { v1, v2 };
		double nx[2], ny[2], nz[2];
		double nxn[2], nyn[2], nzn[2];
		double x1[2], x2[2], y1[2], y2[2], z1[2], z2[2];
		double x0[2], y0[2], z0[2];
		double len2[2], len[2];
		double dx_n1n2, dy_n1n2, dz_n1n2;
		//derivative to v0
		for (size_t j = 0; j < 2; j++)
		{
			x0[j] = x[3 * a2b[v0_id[j]] + 0];
			y0[j] = x[3 * a2b[v0_id[j]] + 1];
			z0[j] = x[3 * a2b[v0_id[j]] + 2];
			x1[j] = x[3 * a2b[v1_id[j]] + 0];
			y1[j] = x[3 * a2b[v1_id[j]] + 1];
			z1[j] = x[3 * a2b[v1_id[j]] + 2];
			x2[j] = x[3 * a2b[v2_id[j]] + 0];
			y2[j] = x[3 * a2b[v2_id[j]] + 1];
			z2[j] = x[3 * a2b[v2_id[j]] + 2];
			nx[j] = (y0[j] - y1[j])*(z0[j] - z2[j]) - (y0[j] - y2[j])*(z0[j] - z1[j]);
			ny[j] = (x0[j] - x2[j])*(z0[j] - z1[j]) - (x0[j] - x1[j])*(z0[j] - z2[j]);
			nz[j] = (x0[j] - x1[j])*(y0[j] - y2[j]) - (x0[j] - x2[j])*(y0[j] - y1[j]);
			len2[j] = nx[j] * nx[j] + ny[j] * ny[j] + nz[j] * nz[j];
			len[j] = std::sqrt(len2[j]);
			nxn[j] = nx[j] / len[j];
			nyn[j] = ny[j] / len[j];
			nzn[j] = nz[j] / len[j];
		}
		double n1n2 = nxn[0] * nxn[1] + nyn[0] * nyn[1] + nzn[0] * nzn[1];
		//energy term below
		double n1n22 = n1n2 * n1n2;
		//normal_e += tmp_mu * n1n22;
		*f += coeff_fea * n1n22;
		double dx_nxn[2], dy_nyn[2], dz_nzn[2], dy_nxn[2], dz_nxn[2], dx_nyn[2], dz_nyn[2], dx_nzn[2], dy_nzn[2];
		//calculate these item
		for (size_t j = 0; j < 2; j++)
		{
			//coordinates of 3 points, first dimension index, second dimension index
			double tmp_coord[3][3] = {
				{x0[j], y0[j], z0[j]},
				{x1[j], y1[j], z1[j]},
				{x2[j], y2[j], z2[j]},
			};
			double normal_coord[3] = { nx[j], ny[j], nz[j] };
			//der v0
			double der_normal[3][3] = {  //derivative corresponds to coordinates, first dimension: normal axis, second dimension: derivative axis
				{0.0, z1[j] - z2[j], y2[j] - y1[j] },
				{z2[j] - z1[j], 0.0, x1[j] - x2[j] },
				{y1[j] - y2[j], x2[j] - x1[j], 0.0 },
			};
			double der_len[3];
			for (size_t a = 0; a < 3; a++)
			{
				der_len[a] = (der_normal[0][a] * normal_coord[0] + der_normal[1][a] * normal_coord[1] + der_normal[2][a] * normal_coord[2]) / len[j];
			}
			double der_normaln[3][3];
			for (size_t a = 0; a < 3; a++)
			{
				for (size_t b = 0; b < 3; b++)
				{
					der_normaln[a][b] = (der_normal[a][b] * len[j] - normal_coord[a] * der_len[b]) / len2[j];
					//der_normaln[a][b] = 0.f;
				}
			}
			dx_nxn[j] = der_normaln[0][0]; dy_nxn[j] = der_normaln[0][1]; dz_nxn[j] = der_normaln[0][2];
			dx_nyn[j] = der_normaln[1][0]; dy_nyn[j] = der_normaln[1][1]; dz_nyn[j] = der_normaln[1][2];
			dx_nzn[j] = der_normaln[2][0]; dy_nzn[j] = der_normaln[2][1]; dz_nzn[j] = der_normaln[2][2];
		}
		//derivative to be calculated
		
		dx_n1n2 = nyn[0] * dx_nyn[1] + dx_nyn[0] * nyn[1] + nzn[0] * dx_nzn[1] + dx_nzn[0] * nzn[1] + nxn[0] * dx_nxn[1] + dx_nxn[0] * nxn[1];
		dy_n1n2 = nxn[0] * dy_nxn[1] + dy_nxn[0] * nxn[1] + nzn[0] * dy_nzn[1] + dy_nzn[0] * nzn[1] + nyn[0] * dy_nyn[1] + dy_nyn[0] * nyn[1];
		dz_n1n2 = nxn[0] * dz_nxn[1] + dz_nxn[0] * nxn[1] + nyn[0] * dz_nyn[1] + dz_nyn[0] * nyn[1] + nzn[0] * dz_nzn[1] + dz_nzn[0] * nzn[1];
		g[ 3 * a2b[v0] + 0] += 2 * coeff_fea * n1n2 * dx_n1n2;
		g[ 3 * a2b[v0] + 1] += 2 * coeff_fea * n1n2 * dy_n1n2;
		g[ 3 * a2b[v0] + 2] += 2 * coeff_fea * n1n2 * dz_n1n2;
		//derivative of v1
		v0_id[0] = v1;
		v0_id[1] = v1;
		v1_id[0] = v0;
		v1_id[1] = v2;
		v2_id[0] = v3;
		v2_id[1] = v0;
		double len_v1[2];
		for (size_t j = 0; j < 2; j++)
		{
			x0[j] = x[3 * a2b[v0_id[j]] + 0];
			y0[j] = x[3 * a2b[v0_id[j]] + 1];
			z0[j] = x[3 * a2b[v0_id[j]] + 2];
			x1[j] = x[3 * a2b[v1_id[j]] + 0];
			y1[j] = x[3 * a2b[v1_id[j]] + 1];
			z1[j] = x[3 * a2b[v1_id[j]] + 2];
			x2[j] = x[3 * a2b[v2_id[j]] + 0];
			y2[j] = x[3 * a2b[v2_id[j]] + 1];
			z2[j] = x[3 * a2b[v2_id[j]] + 2];
			/*nx[j] = (y0[j] - y1[j])*(z0[j] - z2[j]) - (y0[j] - y2[j])*(z0[j] - z1[j]);
			ny[j] = (x0[j] - x2[j])*(z0[j] - z1[j]) - (x0[j] - x1[j])*(z0[j] - z2[j]);
			nz[j] = (x0[j] - x1[j])*(y0[j] - y2[j]) - (x0[j] - x2[j])*(y0[j] - y1[j]);
			len2[j] = nx[j] * nx[j] + ny[j] * ny[j] + nz[j] * nz[j];
			len_v1[j] = std::sqrt(len2[j]);
			nxn[j] = nx[j] / len_v1[j];
			nyn[j] = ny[j] / len_v1[j];
			nzn[j] = nz[j] / len_v1[j];*/
		}
		
		//assert(abs(len[0] - len_v1[0]) < 0.00001 && abs(len[1] - len_v1[1]) < 0.00001);
		for (size_t j = 0; j < 2; j++)
		{
			//coordinates of 3 points, first dimension index, second dimension index
			double tmp_coord[3][3] = {
				{x0[j], y0[j], z0[j]},
				{x1[j], y1[j], z1[j]},
				{x2[j], y2[j], z2[j]},
			};
			double normal_coord[3] = { nx[j], ny[j], nz[j] };
			//der v0
			double der_normal[3][3] = {  //derivative corresponds to coordinates, first dimension: normal axis, second dimension: derivative axis
				{0.0, z1[j] - z2[j], y2[j] - y1[j] },
				{z2[j] - z1[j], 0.0, x1[j] - x2[j] },
				{y1[j] - y2[j], x2[j] - x1[j], 0.0 },
			};
			double der_len[3];
			for (size_t a = 0; a < 3; a++)
			{
				der_len[a] = (der_normal[0][a] + der_normal[1][a] + der_normal[2][a]) / len[j];
			}
			double der_normaln[3][3];
			for (size_t a = 0; a < 3; a++)
			{
				for (size_t b = 0; b < 3; b++)
				{
					der_normaln[a][b] = (der_normal[a][b] * len[j] - normal_coord[a] * der_len[b]) / len2[j];
					//der_normaln[a][b] = 0.f;
				}
			}
			dx_nxn[j] = der_normaln[0][0]; dy_nxn[j] = der_normaln[0][1]; dz_nxn[j] = der_normaln[0][2];
			dx_nyn[j] = der_normaln[1][0]; dy_nyn[j] = der_normaln[1][1]; dz_nyn[j] = der_normaln[1][2];
			dx_nzn[j] = der_normaln[2][0]; dy_nzn[j] = der_normaln[2][1]; dz_nzn[j] = der_normaln[2][2];
		}
		//derivative to be calculated
		dx_n1n2 = nyn[0] * dx_nyn[1] + dx_nyn[0] * nyn[1] + nzn[0] * dx_nzn[1] + dx_nzn[0] * nzn[1] + nxn[0] * dx_nxn[1] + dx_nxn[0] * nxn[1];
		dy_n1n2 = nxn[0] * dy_nxn[1] + dy_nxn[0] * nxn[1] + nzn[0] * dy_nzn[1] + dy_nzn[0] * nzn[1] + nyn[0] * dy_nyn[1] + dy_nyn[0] * nyn[1];
		dz_n1n2 = nxn[0] * dz_nxn[1] + dz_nxn[0] * nxn[1] + nyn[0] * dz_nyn[1] + dz_nyn[0] * nyn[1] + nzn[0] * dz_nzn[1] + dz_nzn[0] * nzn[1];
		g[3 * a2b[v1] + 0] += 2 * coeff_fea * n1n2 * dx_n1n2;
		g[3 * a2b[v1] + 1] += 2 * coeff_fea * n1n2 * dy_n1n2;
		g[3 * a2b[v1] + 2] += 2 * coeff_fea * n1n2 * dz_n1n2;
		
		//derivative of v2
		v0_id[1] = v2;
		v1_id[1] = v0;
		v2_id[1] = v1;
		for (size_t j = 1; j < 2; j++)
		{
			x0[j] = x[3 * a2b[v0_id[j]] + 0];
			y0[j] = x[3 * a2b[v0_id[j]] + 1];
			z0[j] = x[3 * a2b[v0_id[j]] + 2];
			x1[j] = x[3 * a2b[v1_id[j]] + 0];
			y1[j] = x[3 * a2b[v1_id[j]] + 1];
			z1[j] = x[3 * a2b[v1_id[j]] + 2];
			x2[j] = x[3 * a2b[v2_id[j]] + 0];
			y2[j] = x[3 * a2b[v2_id[j]] + 1];
			z2[j] = x[3 * a2b[v2_id[j]] + 2];
		}
		for (size_t j = 1; j < 2; j++)
		{
			//coordinates of 3 points, first dimension index, second dimension index
			double tmp_coord[3][3] = {
				{x0[j], y0[j], z0[j]},
				{x1[j], y1[j], z1[j]},
				{x2[j], y2[j], z2[j]},
			};
			double normal_coord[3] = { nx[j], ny[j], nz[j] };
			//der v0
			double der_normal[3][3] = {  //derivative corresponds to coordinates, first dimension: normal axis, second dimension: derivative axis
				{0.0, z1[j] - z2[j], y2[j] - y1[j] },
				{z2[j] - z1[j], 0.0, x1[j] - x2[j] },
				{y1[j] - y2[j], x2[j] - x1[j], 0.0 },
			};
			double der_len[3];
			for (size_t a = 0; a < 3; a++)
			{
				der_len[a] = (der_normal[0][a] + der_normal[1][a] + der_normal[2][a]) / len[j];
			}
			double der_normaln[3][3];
			for (size_t a = 0; a < 3; a++)
			{
				for (size_t b = 0; b < 3; b++)
				{
					der_normaln[a][b] = (der_normal[a][b] * len[j] - normal_coord[a] * der_len[b]) / len2[j];
					//der_normaln[a][b] = 0.f;
				}
			}
			dx_nxn[j] = der_normaln[0][0]; dy_nxn[j] = der_normaln[0][1]; dz_nxn[j] = der_normaln[0][2];
			dx_nyn[j] = der_normaln[1][0]; dy_nyn[j] = der_normaln[1][1]; dz_nyn[j] = der_normaln[1][2];
			dx_nzn[j] = der_normaln[2][0]; dy_nzn[j] = der_normaln[2][1]; dz_nzn[j] = der_normaln[2][2];
		}
		dx_n1n2 = nxn[0] * dx_nxn[1] + nyn[0] * dx_nyn[1] + nzn[0] * dx_nzn[1];
		dy_n1n2 = nxn[0] * dy_nxn[1] + nyn[0] * dy_nyn[1] + nzn[0] * dy_nzn[1];
		dz_n1n2 = nxn[0] * dz_nxn[1] + nyn[0] * dz_nyn[1] + nzn[0] * dz_nzn[1];
		g[3 * a2b[v2] + 0] += 2 * coeff_fea * n1n2 * dx_n1n2;
		g[3 * a2b[v2] + 1] += 2 * coeff_fea * n1n2 * dy_n1n2;
		g[3 * a2b[v2] + 2] += 2 * coeff_fea * n1n2 * dz_n1n2;
		//derivative of v3
		v0_id[0] = v3;
		v1_id[0] = v1;
		v2_id[0] = v0;
		for (size_t j = 0; j < 1; j++)
		{
			x0[j] = x[3 * a2b[v0_id[j]] + 0];
			y0[j] = x[3 * a2b[v0_id[j]] + 1];
			z0[j] = x[3 * a2b[v0_id[j]] + 2];
			x1[j] = x[3 * a2b[v1_id[j]] + 0];
			y1[j] = x[3 * a2b[v1_id[j]] + 1];
			z1[j] = x[3 * a2b[v1_id[j]] + 2];
			x2[j] = x[3 * a2b[v2_id[j]] + 0];
			y2[j] = x[3 * a2b[v2_id[j]] + 1];
			z2[j] = x[3 * a2b[v2_id[j]] + 2];
		}
		for (size_t j = 0; j < 1; j++)
		{
			//coordinates of 3 points, first dimension index, second dimension index
			double tmp_coord[3][3] = {
				{x0[j], y0[j], z0[j]},
				{x1[j], y1[j], z1[j]},
				{x2[j], y2[j], z2[j]},
			};
			double normal_coord[3] = { nx[j], ny[j], nz[j] };
			//der v0
			double der_normal[3][3] = {  //derivative corresponds to coordinates, first dimension: normal axis, second dimension: derivative axis
				{0.0, z1[j] - z2[j], y2[j] - y1[j] },
				{z2[j] - z1[j], 0.0, x1[j] - x2[j] },
				{y1[j] - y2[j], x2[j] - x1[j], 0.0 },
			};
			double der_len[3];
			for (size_t a = 0; a < 3; a++)
			{
				der_len[a] = (der_normal[0][a] + der_normal[1][a] + der_normal[2][a]) / len[j];
			}
			double der_normaln[3][3];
			for (size_t a = 0; a < 3; a++)
			{
				for (size_t b = 0; b < 3; b++)
				{
					der_normaln[a][b] = (der_normal[a][b] * len[j] - normal_coord[a] * der_len[b]) / len2[j];
					//der_normaln[a][b] = 0.f;
				}
			}
			dx_nxn[j] = der_normaln[0][0]; dy_nxn[j] = der_normaln[0][1]; dz_nxn[j] = der_normaln[0][2];
			dx_nyn[j] = der_normaln[1][0]; dy_nyn[j] = der_normaln[1][1]; dz_nyn[j] = der_normaln[1][2];
			dx_nzn[j] = der_normaln[2][0]; dy_nzn[j] = der_normaln[2][1]; dz_nzn[j] = der_normaln[2][2];
		}
		dx_n1n2 = nxn[1] * dx_nxn[0] + nyn[1] * dx_nyn[0] + nzn[1] * dx_nzn[0];
		dy_n1n2 = nxn[1] * dy_nxn[0] + nyn[1] * dy_nyn[0] + nzn[1] * dy_nzn[0];
		dz_n1n2 = nxn[1] * dz_nxn[0] + nyn[1] * dz_nyn[0] + nzn[1] * dz_nzn[0];
		g[3 * a2b[v3] + 0] += 2 * coeff_fea * n1n2 * dx_n1n2;
		g[3 * a2b[v3] + 1] += 2 * coeff_fea * n1n2 * dy_n1n2;
		g[3 * a2b[v3] + 2] += 2 * coeff_fea * n1n2 * dz_n1n2;
		
	}
	
#endif
	feature_edge_size = 0;
	for (size_t p = 0; p < feature_edge_size; p++)
	{
		int eid = feature_edges[p];
		int bf[2];
		assert(bef_id[eid].size() == 2);
		bf[0] = bef_id[eid][0];
		bf[1] = bef_id[eid][1];
		int v0, v1, v2(-1), v3(-1);
		v0 = feature_e2v[eid].first;
		v1 = feature_e2v[eid].second;
		for (size_t i = 0; i < 3; i++)
		{
			//set v2
			int after = (i + 1) % 3;
			for (size_t j = 0; j < 2; j++)
			{
				if (bfv_id[bf[j]][i] == v0 && bfv_id[bf[j]][after] == v1)
				{
					int afterafter = (after + 1) % 3;
					v2 = bfv_id[bf[j]][afterafter];
				}
				if (bfv_id[bf[j]][i] == v1 && bfv_id[bf[j]][after] == v0)
				{
					int afterafter = (after + 1) % 3;
					v3 = bfv_id[bf[j]][afterafter];
				}
			}
		}
		assert(v2 != -1 && v3 != -1);
		//two face test
		std::set<int> allverts;
		for (size_t i = 0; i < 2; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				allverts.insert(bfv_id[bf[i]][j]);
			}
		}
		assert(allverts.size() == 4);
		int v0_id[2] = { v0, v0 };
		int v1_id[2] = { v3, v1 };
		int v2_id[2] = { v1, v2 };
		double nx[2], ny[2], nz[2];
		double nxn[2], nyn[2], nzn[2];
		double x1[2], x2[2], y1[2], y2[2], z1[2], z2[2];
		double x0[2], y0[2], z0[2];
		double len2[2], len[2];
		//double dx_n1n2, dy_n1n2, dz_n1n2;
		double dnx[2][3][3], dny[2][3][3], dnz[2][3][3];
		//derivative to v0
		for (size_t j = 0; j < 2; j++)
		{
			x0[j] = x[3 * a2b[v0_id[j]] + 0];
			y0[j] = x[3 * a2b[v0_id[j]] + 1];
			z0[j] = x[3 * a2b[v0_id[j]] + 2];
			x1[j] = x[3 * a2b[v1_id[j]] + 0];
			y1[j] = x[3 * a2b[v1_id[j]] + 1];
			z1[j] = x[3 * a2b[v1_id[j]] + 2];
			x2[j] = x[3 * a2b[v2_id[j]] + 0];
			y2[j] = x[3 * a2b[v2_id[j]] + 1];
			z2[j] = x[3 * a2b[v2_id[j]] + 2];
			nx[j] = (y0[j] - y1[j])*(z0[j] - z2[j]) - (y0[j] - y2[j])*(z0[j] - z1[j]);
			ny[j] = (x0[j] - x2[j])*(z0[j] - z1[j]) - (x0[j] - x1[j])*(z0[j] - z2[j]);
			nz[j] = (x0[j] - x1[j])*(y0[j] - y2[j]) - (x0[j] - x2[j])*(y0[j] - y1[j]);
			dnx[j][0][0] = 0.0; //x0
			dnx[j][1][0] = 0.0;	//x1
			dnx[j][2][0] = 0.0; //x2
			dnx[j][0][1] = z1[j] - z2[j]; //y0
			dnx[j][1][1] = z2[j] - z0[j];
			dnx[j][2][1] = z0[j] - z1[j];
			dnx[j][0][2] = y2[j] - y1[j]; //z0
			dnx[j][1][2] = y0[j] - y2[j];
			dnx[j][2][2] = y1[j] - y0[j];
			dny[j][0][0] = z2[j] - z1[j]; //x0
			dny[j][1][0] = z0[j] - z2[j]; //x1
			dny[j][2][0] = z1[j] - z0[j]; //x2
			dny[j][0][1] = 0.0; //y0
			dny[j][1][1] = 0.0;
			dny[j][2][1] = 0.0;
			dny[j][0][2] = x1[j] - x2[j]; //z0
			dny[j][1][2] = x2[j] - x0[j];
			dny[j][2][2] = x0[j] - x1[j];
			dnz[j][0][0] = y1[j] - y2[j]; //x0
			dnz[j][1][0] = y2[j] - y0[j]; //x1
			dnz[j][2][0] = y0[j] - y1[j]; //x2
			dnz[j][0][1] = x2[j] - x1[j]; //y0
			dnz[j][1][1] = x0[j] - x2[j];
			dnz[j][2][1] = x1[j] - x0[j];
			dnz[j][0][2] = 0.0; //z0
			dnz[j][1][2] = 0.0;
			dnz[j][2][2] = 0.0;
			len2[j] = nx[j] * nx[j] + ny[j] * ny[j] + nz[j] * nz[j];
			len[j] = std::sqrt(len2[j]);
			if (len[j] < LEN_LB)
			{
				len[j] = LEN_LB;
				len2[j] = len[j] * len[j];
			}
			nxn[j] = nx[j] / len[j];
			nyn[j] = ny[j] / len[j];
			nzn[j] = nz[j] / len[j];
		}
		double n1n2 = nxn[0] * nxn[1] + nyn[0] * nyn[1] + nzn[0] * nzn[1];
		//energy term below
		double n1n22 = n1n2 * n1n2;
		//normal_e += tmp_mu * n1n22;
		*f += coeff_fea * n1n22;
		double dnxn[2][3][3], dnyn[2][3][3], dnzn[2][3][3]; //first dim: face, second dim: vert, third dim: axis
		double dlen[2][3][3];
		
		//calculate above
		for (size_t i = 0; i < 2; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				for (size_t k = 0; k < 3; k++)
				{
					dlen[i][j][k] = (nx[i] * dnx[i][j][k] + ny[i] * dny[i][j][k] + nz[i] * dnz[i][j][k]) / len[i];
				}
			}
		}
		for (size_t i = 0; i < 2; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				for (size_t k = 0; k < 3; k++)
				{
					dnxn[i][j][k] = (dnx[i][j][k] * len[i] - nx[i] * dlen[i][j][k]) / len2[i];
					dnyn[i][j][k] = (dny[i][j][k] * len[i] - ny[i] * dlen[i][j][k]) / len2[i];
					dnzn[i][j][k] = (dnz[i][j][k] * len[i] - nz[i] * dlen[i][j][k]) / len2[i];
				}
			}
		}
		double dn1n2[4][3]; //4 verts, 3 axis
		//3 axis
		for (size_t i = 0; i < 4; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				dn1n2[i][j] = 0.0;
			}
		}
		//v0: 0 0
		for (size_t i = 0; i < 3; i++)
		{
			dn1n2[0][i] = dnxn[0][0][i] * nxn[1] + nxn[0] * dnxn[1][0][i] + dnyn[0][0][i] * nyn[1] + nyn[0] * dnyn[1][0][i] + dnzn[0][0][i] * nzn[1] + nzn[0] * dnzn[1][0][i];
		}
		//v1: 2, 1
		for (size_t i = 0; i < 3; i++)
		{
			dn1n2[1][i] = dnxn[0][2][i] * nxn[1] + nxn[0] * dnxn[1][1][i] + dnyn[0][2][i] * nyn[1] + nyn[0] * dnyn[1][1][i] + dnzn[0][2][i] * nzn[1] + nzn[0] * dnzn[1][1][i];
		}
		//v2: -1 2
		for (size_t i = 0; i < 3; i++)
		{
			dn1n2[2][i] = nxn[0] * dnxn[1][2][i] + nyn[0] * dnyn[1][2][i] + nzn[0] * dnzn[1][2][i];
		}
		//v3: 1, -1
		for (size_t i = 0; i < 3; i++)
		{
			dn1n2[3][i] = dnxn[0][1][i] * nxn[1] + dnyn[0][1][i] * nyn[1] + dnzn[0][1][i] * nzn[1];
		}
		g[3 * a2b[v0] + 0] += 2 * coeff_fea * n1n2 * dn1n2[0][0];
		g[3 * a2b[v0] + 1] += 2 * coeff_fea * n1n2 * dn1n2[0][1];
		g[3 * a2b[v0] + 2] += 2 * coeff_fea * n1n2 * dn1n2[0][2];
						   
		g[3 * a2b[v1] + 0] += 2 * coeff_fea * n1n2 * dn1n2[1][0];
		g[3 * a2b[v1] + 1] += 2 * coeff_fea * n1n2 * dn1n2[1][1];
		g[3 * a2b[v1] + 2] += 2 * coeff_fea * n1n2 * dn1n2[1][2];
						   
		g[3 * a2b[v2] + 0] += 2 * coeff_fea * n1n2 * dn1n2[2][0];
		g[3 * a2b[v2] + 1] += 2 * coeff_fea * n1n2 * dn1n2[2][1];
		g[3 * a2b[v2] + 2] += 2 * coeff_fea * n1n2 * dn1n2[2][2];
						   
		g[3 * a2b[v3] + 0] += 2 * coeff_fea * n1n2 * dn1n2[3][0];
		g[3 * a2b[v3] + 1] += 2 * coeff_fea * n1n2 * dn1n2[3][1];
		g[3 * a2b[v3] + 2] += 2 * coeff_fea * n1n2 * dn1n2[3][2];
	}
	//feature normal constraint
	double coeff_fea_normal = lambda_array[0] / (1.0 * feature_edges.size());
	//assert(feature_edges.size() == target_feature_edge_normal.size());
	int n_fea_normal = feature_edges.size();
	n_fea_normal = 0;
	for (size_t p = 0; p < n_fea_normal; p++)
	{
		int eid = feature_edges[p];
		const OpenVolumeMesh::Geometry::Vec3d target_n = target_fea_n[eid];
		int aid0 = feature_e2v[eid].first, aid1 = feature_e2v[eid].second;
		int bid0 = a2b[aid0], bid1 = a2b[aid1];
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		double n[3];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		double len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		double len = sqrt(len2);
		assert(_finite(len) == 1 && _finite(len2) == 1);
		if (len < LEN_LB)
		{
			len = LEN_LB;
			len2 = len * len;
		}
		double nn[3];
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		assert(_finite(nn[0]) && _finite(nn[1]) && _finite(nn[2]));
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		double f_ori = 0;
		for (size_t i = 0; i < 3; i++)
		{
			*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			f_ori += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		double dlen[2][3]; //first dimension: vert id, second dim: axis
		for (size_t i = 0; i < 3; i++)
		{
			dlen[0][i] = -n[i] / len;
			dlen[1][i] = n[i] / len;
		}
		double diff[3];
		for (size_t i = 0; i < 3; i++)
		{
			diff[i] = nn[i] - target_n[i];
		}
		
		double ddiff0[3][3], ddiff1[3][3]; //second dimension: axis
		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				ddiff0[i][j] = -n[i] * dlen[0][j] / len2;
				ddiff1[i][j] = -n[i] * dlen[1][j] / len2;
			}
			ddiff0[i][i] += -1.0 / len;
			ddiff1[i][i] += 1.0 / len;
		}
		
		double df[2][3];
		for (size_t i = 0; i < 3; i++)
		{
			df[0][i] = 0.0; //i the axis
			df[1][i] = 0.0;
			for (size_t j = 0; j < 3; j++)
			{
				df[0][i] += 2 * coeff_fea_normal * diff[j] * ddiff0[j][i];
				df[1][i] += 2 * coeff_fea_normal * diff[j] * ddiff1[j][i];
			}
		}
		g[3 * bid0 + 0] += df[0][0];
		g[3 * bid0 + 1] += df[0][1];
		g[3 * bid0 + 2] += df[0][2];
		g[3 * bid1 + 0] += df[1][0];
		g[3 * bid1 + 1] += df[1][1];
		g[3 * bid1 + 2] += df[1][2];
#if 0
		//for debugging
		std::cout << "v0 grad: " << std::endl;
		for (size_t i = 0; i < 3; i++)
		{
			std::cout << df[0][i] << std::endl;
		}
		std::cout << "v1 grad: " << std::endl;
		for (size_t i = 0; i < 3; i++)
		{
			std::cout << df[1][i] << std::endl;
		}
		
		//difference method
		//v0x
		double tmpf = 0;
		x1 = x[3 * bid0] + EPSL;
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 x: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][0] - (tmpf - f_ori) / EPSL) < BOUND);
		
		//v0y
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1] + EPSL;
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 y: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][1] - (tmpf - f_ori) / EPSL) < BOUND);
		//v0z
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2] + EPSL;
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 z: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][2] - (tmpf - f_ori) / EPSL) < BOUND);
		//v1x
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1] + EPSL;
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 x: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][0] - (tmpf - f_ori) / EPSL) < BOUND);
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1] + EPSL;
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 y: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][1] - (tmpf - f_ori) / EPSL) < BOUND);
		//v1z
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2] + EPSL;
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 z: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][2] - (tmpf - f_ori) / EPSL) < BOUND);
		tmpf = 0;
#endif
		
	}
}
void evalfunc_de_LBFGS_all(int N, double* x, double *prev_x, double* f, double* g, void* user_pointer)
{
	deform_pointer *dp = (deform_pointer *)(user_pointer);
	const vector<double> &dpx = dp->pdi->dpx;
	const vector<double> &dpy = dp->pdi->dpy;
	const vector<double> &dpz = dp->pdi->dpz;
	const vector<int> &boundary_verts = dp->pdi->boundary_verts;
	const vector<OpenVolumeMesh::Geometry::Vec3i> &bfv_id = dp->pdi->bfv_id;
	const std::vector<std::vector<int>> &vertex_cell = dp->pdi->vertex_cell;
	const std::vector<std::vector<std::vector<int>>> &vertex_cell_vertex = dp->pdi->vertex_cell_vertex;
	const std::vector< OpenVolumeMesh::Geometry::Vec3d > &target_bfn = dp->pdi->target_bfn;
	const std::vector< std::vector<int> > &bvv_id = dp->pdi->bvv_id;
	const std::vector< std::vector<int> > &bvf_id = dp->pdi->bvf_id;
	const std::vector<int> &a2b = dp->pdi->a2b;
	//std::vector<std::vector<int>> &cell_vertex = dp->pdi->cell_vertex; //not using right order version, cause inverse is not in right order
	const std::vector< std::vector<std::vector<double> > > &vcv_S = dp->pdi->vcv_S;
	const std::vector< Eigen::Matrix3d > &cell_S = dp->pdi->cell_S;
	const std::vector<std::vector<int>> &cell_vertex = dp->pdi->cell_vertex;
	const std::vector<double> &lambda_array = dp->pdi->lambda_array;
	int boundary_face_number = dp->pdi->boundary_face_number;
	//int boundary_vertface_number = dp->pdi->boundary_vertface_number;
	int nc = dp->pdi->bfv_id.size();
	const vector<vector<int>> &bef_id = dp->pdi->bef_id;
	//feature edge part
	const std::vector<bool> &feature_edge_flag = dp->pdi->feature_edge_flag;
	const std::vector<std::vector<int>> &feature_v2e = dp->pdi->feature_v2e;
	const std::vector<std::vector<std::pair<int, int>>> &feature_neighbor_vert2cellpair = dp->pdi->feature_neighbor_vert2cellpair; //first face containging vert, second face the other direction
	//std::vector<std::vector<int>> feature_faces;
	const std::vector<OpenVolumeMesh::Geometry::Vec3d> &target_fea_n = dp->pdi->target_fea_n;
	//int n_feature_edges;
	const std::vector<int> &feature_edges = dp->pdi->feature_edge_array;
	const std::vector<std::pair<int, int>> &feature_e2v = dp->pdi->feature_e2v; //consistent feature edges
	const std::vector<std::vector<int>> &feature_neighbor = dp->pdi->feature_neighbor; //neighbor of edges at most 2 neighbors, cut by valence 3 vertex
	//const std::vector<OpenVolumeMesh::Geometry::Vec3d> &target_feature_edge_normal = dp->pdi->target_feature_edge_normal;
	/*int v_id_0, v_id_1, v_id_2, v_id_3, k, c_id;
	double D00, D10, D20, D01, D11, D21, D02, D12, D22;
	double C00, C01, C02, C10, C11, C12, C20, C21, C22;
	double A00, A01, A02, A10, A11, A12, A20, A21, A22;
	double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
	double AF, AF_05, I_AF_05, AF2, AF_I, AF_I_05, I_AF_I_05, exp_k, exp_e;
	double D_AF_x, D_AF_y, D_AF_z, D_AF2_x, D_AF2_y, D_AF2_z, D_AF_I_x, D_AF_I_y, D_AF_I_z, D_det_A_x, D_det_A_y, D_det_A_z, inv_det_A_2_05;
	double d_A00_x, d_A10_y, d_A20_z, d_A01_x, d_A11_y, d_A21_z, d_A02_x, d_A12_y, d_A22_z;
	double dvex, dvey, dvez, tmp_g, dgx, dgy, dgz, dex, dey, dez, e, mips_e, volume_e, d_mips_e_x, d_mips_e_y, d_mips_e_z;
	
	
	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	double len;*/
	//double min_radius = 1e30;
	double alpha = 0.5; //double beta = 0.5;
	int bv_size = boundary_verts.size();
	int nv = dp->pdi->bvf_id.size();
	assert(bv_size > 0);
	//std::vector<std::vector<double>> g_omp;
	unsigned int num_core = omp_get_max_threads();
	std::vector<std::vector<double>> g_omp(num_core, std::vector<double>(3 * nv, 0.0));
	double func = 0.0;
	*f = 0.0;
	for (size_t i = 0; i < 3 * nv; i++)
	{
		g[i] = 0.0;
	}
	
	//consider all cells with at least one element on the boundary
	std::vector<int> cell_flag(bfv_id.size(), -1);
	double coeff_iso = 1.0 / (1.0 * nc);
	//isometric part
	//nv set here
	//nv = 0;
//#pragma omp parallel for reduction(+:*f)
	//not parallelable
	for (int i = 0; i < 0; i++)
	{
		int tid = omp_get_thread_num();
		int id = i;
		//assert(bvf_id[bid].size() > 0);
		//assert(vertex_cell[bid].size() == vertex_cell_vertex[bid].size());
		double D00, D10, D20, D01, D11, D21, D02, D12, D22;
		double C00, C01, C02, C10, C11, C12, C20, C21, C22;
		double A00, A01, A02, A10, A11, A12, A20, A21, A22;
		double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
		double AF, AF_05, I_AF_05, AF2, AF_I, AF_I_05, I_AF_I_05, exp_k;
		double D_AF_x, D_AF_y, D_AF_z, D_AF2_x, D_AF2_y, D_AF2_z, D_AF_I_x, D_AF_I_y, D_AF_I_z, D_det_A_x, D_det_A_y, D_det_A_z, inv_det_A_2_05;
		double d_A00_x, d_A10_y, d_A20_z, d_A01_x, d_A11_y, d_A21_z, d_A02_x, d_A12_y, d_A22_z;
		double dvex, dvey, dvez, tmp_g, dgx, dgy, dgz, dex, dey, dez, e, mips_e, volume_e;
		int vc_size = vertex_cell[id].size();
		vector<double> local_cell_energy(vc_size, 0.0);
		double x0 = x[3 * id], y0 = x[3 * id + 1], z0 = x[3 * id + 2];
		//isometric part
		vector<double> exp_vec(vc_size, 0.0), gx_vec(vc_size, 0.0), gy_vec(vc_size, 0.0), gz_vec(vc_size, 0.0);
		for (size_t j = 0; j < vc_size; j++)
		{
			//decide local_cell_energy and gradient here
			const int* vv_id = vertex_cell_vertex[id][j].data();
			const double* s_data = vcv_S[id][j].data();
			int k = 3 * j;
			int v_id_1 = vv_id[0];
			int v_id_2 = vv_id[1];
			int v_id_3 = vv_id[2];
			double x1 = dpx[v_id_1],  x2 = dpx[v_id_2], x3 = dpx[v_id_3];
			double y1 = dpy[v_id_1],  y2 = dpy[v_id_2], y3 = dpy[v_id_3];
			double z1 = dpz[v_id_1],  z2 = dpz[v_id_2], z3 = dpz[v_id_3];
			//redefine the variables above
			x1 = x[3 * v_id_1 + 0];
			y1 = x[3 * v_id_1 + 1];
			z1 = x[3 * v_id_1 + 2];
			x2 = x[3 * v_id_2 + 0];
			y2 = x[3 * v_id_2 + 1];
			z2 = x[3 * v_id_2 + 2];
			x3 = x[3 * v_id_3 + 0];
			y3 = x[3 * v_id_3 + 1];
			z3 = x[3 * v_id_3 + 2];
			/*p_vc_pos_x[k + 0] = x1; p_vc_pos_x[k + 1] = x2; p_vc_pos_x[k + 2] = x3;
			p_vc_pos_y[k + 0] = y1; p_vc_pos_y[k + 1] = y2; p_vc_pos_y[k + 2] = y3;
			p_vc_pos_z[k + 0] = z1; p_vc_pos_z[k + 1] = z2; p_vc_pos_z[k + 2] = z3;
			p_vc_n_cross_x[i] = (y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1);
			p_vc_n_cross_y[i] = (x3 - x1)*(z2 - z1) - (x2 - x1)*(z3 - z1);
			p_vc_n_cross_z[i] = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);*/
			D00 = x1 - x0; D10 = y1 - y0; D20 = z1 - z0;
			D01 = x2 - x0; D11 = y2 - y0; D21 = z2 - z0;
			D02 = x3 - x0; D12 = y3 - y0; D22 = z3 - z0;
			C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
			C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
			C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];
			A00 = C00 * D00 + C10 * D01 + C20 * D02; d_A00_x = -C00 - C10 - C20;
			A10 = C00 * D10 + C10 * D11 + C20 * D12; d_A10_y = -C00 - C10 - C20;
			A20 = C00 * D20 + C10 * D21 + C20 * D22; d_A20_z = -C00 - C10 - C20;
			A01 = C01 * D00 + C11 * D01 + C21 * D02; d_A01_x = -C01 - C11 - C21;
			A11 = C01 * D10 + C11 * D11 + C21 * D12; d_A11_y = -C01 - C11 - C21;
			A21 = C01 * D20 + C11 * D21 + C21 * D22; d_A21_z = -C01 - C11 - C21;
			A02 = C02 * D00 + C12 * D01 + C22 * D02; d_A02_x = -C02 - C12 - C22;
			A12 = C02 * D10 + C12 * D11 + C22 * D12; d_A12_y = -C02 - C12 - C22;
			A22 = C02 * D20 + C12 * D21 + C22 * D22; d_A22_z = -C02 - C12 - C22;
			A2_00 = A00 * A00 + A10 * A10 + A20 * A20;
			A2_01 = A00 * A01 + A10 * A11 + A20 * A21;
			A2_02 = A00 * A02 + A10 * A12 + A20 * A22;
			A2_11 = A01 * A01 + A11 * A11 + A21 * A21;
			A2_12 = A01 * A02 + A11 * A12 + A21 * A22;
			A2_22 = A02 * A02 + A12 * A12 + A22 * A22;
			AF = A2_00 + A2_11 + A2_22; AF_05 = std::sqrt(AF); I_AF_05 = 0.5 / AF_05;
			AF2 = A2_00 * A2_00 + A2_01 * A2_01 * 2 + A2_02 * A2_02 * 2 + A2_11 * A2_11 + A2_12 * A2_12 * 2 + A2_22 * A2_22;
			AF_I = (AF*AF - AF2)*0.5; AF_I_05 = std::sqrt(AF_I); I_AF_I_05 = 0.5 / AF_I_05;
			det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
			i_det_A = 1.0 / det_A;
			D_AF_x = (2.0*A00*d_A00_x + 2.0*A01*d_A01_x + 2.0*A02*d_A02_x);
			D_AF_y = (2.0*A10*d_A10_y + 2.0*A11*d_A11_y + 2.0*A12*d_A12_y);
			D_AF_z = (2.0*A20*d_A20_z + 2.0*A21*d_A21_z + 2.0*A22*d_A22_z);
			D_AF2_x = (2.0*A2_00*2.0*A00*d_A00_x + 4.0*A2_01*(A01*d_A00_x + A00 * d_A01_x) + 4.0*A2_02*(A02*d_A00_x + A00 * d_A02_x)
				+ 2.0*A2_11*2.0*A01*d_A01_x + 2.0*A2_22*2.0*A02*d_A02_x + 4.0*A2_12*(A01*d_A02_x + A02 * d_A01_x));
			D_AF2_y = (2.0*A2_00*2.0*A10*d_A10_y + 4.0*A2_01*(A11*d_A10_y + A10 * d_A11_y) + 4.0*A2_02*(A12*d_A10_y + A10 * d_A12_y)
				+ 2.0*A2_11*2.0*A11*d_A11_y + 2.0*A2_22*2.0*A12*d_A12_y + 4.0*A2_12*(A11*d_A12_y + A12 * d_A11_y));
			D_AF2_z = (2.0*A2_00*2.0*A20*d_A20_z + 4.0*A2_01*(A21*d_A20_z + A20 * d_A21_z) + 4.0*A2_02*(A22*d_A20_z + A20 * d_A22_z)
				+ 2.0*A2_11*2.0*A21*d_A21_z + 2.0*A2_22*2.0*A22*d_A22_z + 4.0*A2_12*(A21*d_A22_z + A22 * d_A21_z));
			D_AF_I_x = (AF*D_AF_x - 0.5*D_AF2_x)*I_AF_I_05; D_AF_I_y = (AF*D_AF_y - 0.5*D_AF2_y)*I_AF_I_05; D_AF_I_z = (AF*D_AF_z - 0.5*D_AF2_z)*I_AF_I_05;
			D_AF_x *= I_AF_05; D_AF_y *= I_AF_05; D_AF_z *= I_AF_05;
			D_det_A_x = d_A00_x * A11 * A22 + d_A01_x * A12 * A20 + d_A02_x * A10 * A21 - d_A02_x * A11 * A20 - d_A01_x * A10 * A22 - d_A00_x * A12 * A21;
			D_det_A_y = A00 * d_A11_y * A22 + A01 * d_A12_y * A20 + A02 * d_A10_y * A21 - A02 * d_A11_y * A20 - A01 * d_A10_y * A22 - A00 * d_A12_y * A21;
			D_det_A_z = A00 * A11 * d_A22_z + A01 * A12 * d_A20_z + A02 * A10 * d_A21_z - A02 * A11 * d_A20_z - A01 * A10 * d_A22_z - A00 * A12 * d_A21_z;
			inv_det_A_2_05 = 0.5 *(1.0 - i_det_A * i_det_A);
			dvex = D_det_A_x * inv_det_A_2_05;
			dvey = D_det_A_y * inv_det_A_2_05;
			dvez = D_det_A_z * inv_det_A_2_05;
			tmp_g = AF_05 * AF_I_05;
			dgx = D_AF_x * AF_I_05 + AF_05 * D_AF_I_x;
			dgy = D_AF_y * AF_I_05 + AF_05 * D_AF_I_y;
			dgz = D_AF_z * AF_I_05 + AF_05 * D_AF_I_z;
			dex = (dgx *i_det_A - tmp_g * D_det_A_x *i_det_A*i_det_A);
			dey = (dgy *i_det_A - tmp_g * D_det_A_y *i_det_A*i_det_A);
			dez = (dgz *i_det_A - tmp_g * D_det_A_z *i_det_A*i_det_A);
			e = tmp_g * i_det_A;
			mips_e = (e*e - 1.0)*0.125;
			//mips_e = e/3.0;
			volume_e = 0.5*(det_A + i_det_A);
			exp_k = (alpha*mips_e + (1.0 - alpha) * volume_e);
			if (exp_k > 60) exp_k = 60;
			exp_vec[j] = exp_k;
			/*exp_e = std::exp(exp_k);
			local_energy += exp_e;*/
			double d_mips_e_x = e * 0.25 * dex;
			double d_mips_e_y = e * 0.25 * dey;
			double d_mips_e_z = e * 0.25 * dez;
			/*d_mips_e_x = dex/3.0;
			d_mips_e_y = dey/3.0;
			d_mips_e_z = dez/3.0;*/
			gx_vec[j] = alpha * d_mips_e_x + (1.0 - alpha) * dvex;
			gy_vec[j] = alpha * d_mips_e_y + (1.0 - alpha) * dvey;
			gz_vec[j] = alpha * d_mips_e_z + (1.0 - alpha) * dvez;
#if 0
			//debug below
			double diff[3];
			double tmp_f = 0.0;
			//x
			//std::cout << "x: " << std::endl;
			//std::cout << gx_vec[j] << std::endl;
			x0 = x[3 * i] + EPSL;
			D00 = x1 - x0; D10 = y1 - y0; D20 = z1 - z0;
			D01 = x2 - x0; D11 = y2 - y0; D21 = z2 - z0;
			D02 = x3 - x0; D12 = y3 - y0; D22 = z3 - z0;
			len = D00 * D00 + D10 * D10 + D20 * D20; if (len < min_radius) min_radius = len;
			len = D01 * D01 + D11 * D11 + D21 * D21; if (len < min_radius) min_radius = len;
			len = D02 * D02 + D12 * D12 + D22 * D22; if (len < min_radius) min_radius = len;
			C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
			C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
			C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];
			A00 = C00 * D00 + C10 * D01 + C20 * D02; d_A00_x = -C00 - C10 - C20;
			A10 = C00 * D10 + C10 * D11 + C20 * D12; d_A10_y = -C00 - C10 - C20;
			A20 = C00 * D20 + C10 * D21 + C20 * D22; d_A20_z = -C00 - C10 - C20;
			A01 = C01 * D00 + C11 * D01 + C21 * D02; d_A01_x = -C01 - C11 - C21;
			A11 = C01 * D10 + C11 * D11 + C21 * D12; d_A11_y = -C01 - C11 - C21;
			A21 = C01 * D20 + C11 * D21 + C21 * D22; d_A21_z = -C01 - C11 - C21;
			A02 = C02 * D00 + C12 * D01 + C22 * D02; d_A02_x = -C02 - C12 - C22;
			A12 = C02 * D10 + C12 * D11 + C22 * D12; d_A12_y = -C02 - C12 - C22;
			A22 = C02 * D20 + C12 * D21 + C22 * D22; d_A22_z = -C02 - C12 - C22;
			A2_00 = A00 * A00 + A10 * A10 + A20 * A20;
			A2_01 = A00 * A01 + A10 * A11 + A20 * A21;
			A2_02 = A00 * A02 + A10 * A12 + A20 * A22;
			A2_11 = A01 * A01 + A11 * A11 + A21 * A21;
			A2_12 = A01 * A02 + A11 * A12 + A21 * A22;
			A2_22 = A02 * A02 + A12 * A12 + A22 * A22;
			AF = A2_00 + A2_11 + A2_22; AF_05 = std::sqrt(AF); I_AF_05 = 0.5 / AF_05;
			AF2 = A2_00 * A2_00 + A2_01 * A2_01 * 2 + A2_02 * A2_02 * 2 + A2_11 * A2_11 + A2_12 * A2_12 * 2 + A2_22 * A2_22;
			AF_I = (AF*AF - AF2)*0.5; AF_I_05 = std::sqrt(AF_I); I_AF_I_05 = 0.5 / AF_I_05;
			det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
			i_det_A = 1.0 / det_A;
			D_AF_x = (2.0*A00*d_A00_x + 2.0*A01*d_A01_x + 2.0*A02*d_A02_x);
			D_AF_y = (2.0*A10*d_A10_y + 2.0*A11*d_A11_y + 2.0*A12*d_A12_y);
			D_AF_z = (2.0*A20*d_A20_z + 2.0*A21*d_A21_z + 2.0*A22*d_A22_z);
			D_AF2_x = (2.0*A2_00*2.0*A00*d_A00_x + 4.0*A2_01*(A01*d_A00_x + A00 * d_A01_x) + 4.0*A2_02*(A02*d_A00_x + A00 * d_A02_x)
				+ 2.0*A2_11*2.0*A01*d_A01_x + 2.0*A2_22*2.0*A02*d_A02_x + 4.0*A2_12*(A01*d_A02_x + A02 * d_A01_x));
			D_AF2_y = (2.0*A2_00*2.0*A10*d_A10_y + 4.0*A2_01*(A11*d_A10_y + A10 * d_A11_y) + 4.0*A2_02*(A12*d_A10_y + A10 * d_A12_y)
				+ 2.0*A2_11*2.0*A11*d_A11_y + 2.0*A2_22*2.0*A12*d_A12_y + 4.0*A2_12*(A11*d_A12_y + A12 * d_A11_y));
			D_AF2_z = (2.0*A2_00*2.0*A20*d_A20_z + 4.0*A2_01*(A21*d_A20_z + A20 * d_A21_z) + 4.0*A2_02*(A22*d_A20_z + A20 * d_A22_z)
				+ 2.0*A2_11*2.0*A21*d_A21_z + 2.0*A2_22*2.0*A22*d_A22_z + 4.0*A2_12*(A21*d_A22_z + A22 * d_A21_z));
			D_AF_I_x = (AF*D_AF_x - 0.5*D_AF2_x)*I_AF_I_05; D_AF_I_y = (AF*D_AF_y - 0.5*D_AF2_y)*I_AF_I_05; D_AF_I_z = (AF*D_AF_z - 0.5*D_AF2_z)*I_AF_I_05;
			D_AF_x *= I_AF_05; D_AF_y *= I_AF_05; D_AF_z *= I_AF_05;
			D_det_A_x = d_A00_x * A11 * A22 + d_A01_x * A12 * A20 + d_A02_x * A10 * A21 - d_A02_x * A11 * A20 - d_A01_x * A10 * A22 - d_A00_x * A12 * A21;
			D_det_A_y = A00 * d_A11_y * A22 + A01 * d_A12_y * A20 + A02 * d_A10_y * A21 - A02 * d_A11_y * A20 - A01 * d_A10_y * A22 - A00 * d_A12_y * A21;
			D_det_A_z = A00 * A11 * d_A22_z + A01 * A12 * d_A20_z + A02 * A10 * d_A21_z - A02 * A11 * d_A20_z - A01 * A10 * d_A22_z - A00 * A12 * d_A21_z;
			inv_det_A_2_05 = 0.5 *(1.0 - i_det_A * i_det_A);
			dvex = D_det_A_x * inv_det_A_2_05;
			dvey = D_det_A_y * inv_det_A_2_05;
			dvez = D_det_A_z * inv_det_A_2_05;
			tmp_g = AF_05 * AF_I_05;
			dgx = D_AF_x * AF_I_05 + AF_05 * D_AF_I_x;
			dgy = D_AF_y * AF_I_05 + AF_05 * D_AF_I_y;
			dgz = D_AF_z * AF_I_05 + AF_05 * D_AF_I_z;
			dex = (dgx *i_det_A - tmp_g * D_det_A_x *i_det_A*i_det_A);
			dey = (dgy *i_det_A - tmp_g * D_det_A_y *i_det_A*i_det_A);
			dez = (dgz *i_det_A - tmp_g * D_det_A_z *i_det_A*i_det_A);
			e = tmp_g * i_det_A;
			mips_e = (e*e - 1.0)*0.125;
			//mips_e = e/3.0;
			volume_e = 0.5*(det_A + i_det_A);
			exp_k = (alpha*mips_e + beta * volume_e);
			if (exp_k > 60) exp_k = 60;
			tmp_f = exp_k;
			//std::cout << (tmp_f - exp_vec[j]) / EPSL << std::endl;
			diff[0] = (tmp_f - exp_vec[j]) / EPSL;
			//y
			//std::cout << "y: " << std::endl;
			//std::cout << gy_vec[j] << std::endl;
			tmp_f = 0.0;
			x0 = x[3 * i];
			y0 = x[3 * i + 1] + EPSL;
			D00 = x1 - x0; D10 = y1 - y0; D20 = z1 - z0;
			D01 = x2 - x0; D11 = y2 - y0; D21 = z2 - z0;
			D02 = x3 - x0; D12 = y3 - y0; D22 = z3 - z0;
			len = D00 * D00 + D10 * D10 + D20 * D20; if (len < min_radius) min_radius = len;
			len = D01 * D01 + D11 * D11 + D21 * D21; if (len < min_radius) min_radius = len;
			len = D02 * D02 + D12 * D12 + D22 * D22; if (len < min_radius) min_radius = len;
			C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
			C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
			C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];
			A00 = C00 * D00 + C10 * D01 + C20 * D02; d_A00_x = -C00 - C10 - C20;
			A10 = C00 * D10 + C10 * D11 + C20 * D12; d_A10_y = -C00 - C10 - C20;
			A20 = C00 * D20 + C10 * D21 + C20 * D22; d_A20_z = -C00 - C10 - C20;
			A01 = C01 * D00 + C11 * D01 + C21 * D02; d_A01_x = -C01 - C11 - C21;
			A11 = C01 * D10 + C11 * D11 + C21 * D12; d_A11_y = -C01 - C11 - C21;
			A21 = C01 * D20 + C11 * D21 + C21 * D22; d_A21_z = -C01 - C11 - C21;
			A02 = C02 * D00 + C12 * D01 + C22 * D02; d_A02_x = -C02 - C12 - C22;
			A12 = C02 * D10 + C12 * D11 + C22 * D12; d_A12_y = -C02 - C12 - C22;
			A22 = C02 * D20 + C12 * D21 + C22 * D22; d_A22_z = -C02 - C12 - C22;
			A2_00 = A00 * A00 + A10 * A10 + A20 * A20;
			A2_01 = A00 * A01 + A10 * A11 + A20 * A21;
			A2_02 = A00 * A02 + A10 * A12 + A20 * A22;
			A2_11 = A01 * A01 + A11 * A11 + A21 * A21;
			A2_12 = A01 * A02 + A11 * A12 + A21 * A22;
			A2_22 = A02 * A02 + A12 * A12 + A22 * A22;
			AF = A2_00 + A2_11 + A2_22; AF_05 = std::sqrt(AF); I_AF_05 = 0.5 / AF_05;
			AF2 = A2_00 * A2_00 + A2_01 * A2_01 * 2 + A2_02 * A2_02 * 2 + A2_11 * A2_11 + A2_12 * A2_12 * 2 + A2_22 * A2_22;
			AF_I = (AF*AF - AF2)*0.5; AF_I_05 = std::sqrt(AF_I); I_AF_I_05 = 0.5 / AF_I_05;
			det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
			i_det_A = 1.0 / det_A;
			D_AF_x = (2.0*A00*d_A00_x + 2.0*A01*d_A01_x + 2.0*A02*d_A02_x);
			D_AF_y = (2.0*A10*d_A10_y + 2.0*A11*d_A11_y + 2.0*A12*d_A12_y);
			D_AF_z = (2.0*A20*d_A20_z + 2.0*A21*d_A21_z + 2.0*A22*d_A22_z);
			D_AF2_x = (2.0*A2_00*2.0*A00*d_A00_x + 4.0*A2_01*(A01*d_A00_x + A00 * d_A01_x) + 4.0*A2_02*(A02*d_A00_x + A00 * d_A02_x)
				+ 2.0*A2_11*2.0*A01*d_A01_x + 2.0*A2_22*2.0*A02*d_A02_x + 4.0*A2_12*(A01*d_A02_x + A02 * d_A01_x));
			D_AF2_y = (2.0*A2_00*2.0*A10*d_A10_y + 4.0*A2_01*(A11*d_A10_y + A10 * d_A11_y) + 4.0*A2_02*(A12*d_A10_y + A10 * d_A12_y)
				+ 2.0*A2_11*2.0*A11*d_A11_y + 2.0*A2_22*2.0*A12*d_A12_y + 4.0*A2_12*(A11*d_A12_y + A12 * d_A11_y));
			D_AF2_z = (2.0*A2_00*2.0*A20*d_A20_z + 4.0*A2_01*(A21*d_A20_z + A20 * d_A21_z) + 4.0*A2_02*(A22*d_A20_z + A20 * d_A22_z)
				+ 2.0*A2_11*2.0*A21*d_A21_z + 2.0*A2_22*2.0*A22*d_A22_z + 4.0*A2_12*(A21*d_A22_z + A22 * d_A21_z));
			D_AF_I_x = (AF*D_AF_x - 0.5*D_AF2_x)*I_AF_I_05; D_AF_I_y = (AF*D_AF_y - 0.5*D_AF2_y)*I_AF_I_05; D_AF_I_z = (AF*D_AF_z - 0.5*D_AF2_z)*I_AF_I_05;
			D_AF_x *= I_AF_05; D_AF_y *= I_AF_05; D_AF_z *= I_AF_05;
			D_det_A_x = d_A00_x * A11 * A22 + d_A01_x * A12 * A20 + d_A02_x * A10 * A21 - d_A02_x * A11 * A20 - d_A01_x * A10 * A22 - d_A00_x * A12 * A21;
			D_det_A_y = A00 * d_A11_y * A22 + A01 * d_A12_y * A20 + A02 * d_A10_y * A21 - A02 * d_A11_y * A20 - A01 * d_A10_y * A22 - A00 * d_A12_y * A21;
			D_det_A_z = A00 * A11 * d_A22_z + A01 * A12 * d_A20_z + A02 * A10 * d_A21_z - A02 * A11 * d_A20_z - A01 * A10 * d_A22_z - A00 * A12 * d_A21_z;
			inv_det_A_2_05 = 0.5 *(1.0 - i_det_A * i_det_A);
			dvex = D_det_A_x * inv_det_A_2_05;
			dvey = D_det_A_y * inv_det_A_2_05;
			dvez = D_det_A_z * inv_det_A_2_05;
			tmp_g = AF_05 * AF_I_05;
			dgx = D_AF_x * AF_I_05 + AF_05 * D_AF_I_x;
			dgy = D_AF_y * AF_I_05 + AF_05 * D_AF_I_y;
			dgz = D_AF_z * AF_I_05 + AF_05 * D_AF_I_z;
			dex = (dgx *i_det_A - tmp_g * D_det_A_x *i_det_A*i_det_A);
			dey = (dgy *i_det_A - tmp_g * D_det_A_y *i_det_A*i_det_A);
			dez = (dgz *i_det_A - tmp_g * D_det_A_z *i_det_A*i_det_A);
			e = tmp_g * i_det_A;
			mips_e = (e*e - 1.0)*0.125;
			//mips_e = e/3.0;
			volume_e = 0.5*(det_A + i_det_A);
			exp_k = (alpha*mips_e + beta * volume_e);
			if (exp_k > 60) exp_k = 60;
			tmp_f = exp_k;
			//std::cout << (tmp_f - exp_vec[j]) / EPSL << std::endl;
			diff[1] = (tmp_f - exp_vec[j]) / EPSL;
			//z
			//std::cout << "z: " << std::endl;
			//std::cout << gz_vec[j] << std::endl;
			tmp_f = 0.0;
			x0 = x[3 * i];
			y0 = x[3 * i + 1];
			z0 = x[3 * i + 2] + EPSL;
			D00 = x1 - x0; D10 = y1 - y0; D20 = z1 - z0;
			D01 = x2 - x0; D11 = y2 - y0; D21 = z2 - z0;
			D02 = x3 - x0; D12 = y3 - y0; D22 = z3 - z0;
			len = D00 * D00 + D10 * D10 + D20 * D20; if (len < min_radius) min_radius = len;
			len = D01 * D01 + D11 * D11 + D21 * D21; if (len < min_radius) min_radius = len;
			len = D02 * D02 + D12 * D12 + D22 * D22; if (len < min_radius) min_radius = len;
			C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
			C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
			C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];
			A00 = C00 * D00 + C10 * D01 + C20 * D02; d_A00_x = -C00 - C10 - C20;
			A10 = C00 * D10 + C10 * D11 + C20 * D12; d_A10_y = -C00 - C10 - C20;
			A20 = C00 * D20 + C10 * D21 + C20 * D22; d_A20_z = -C00 - C10 - C20;
			A01 = C01 * D00 + C11 * D01 + C21 * D02; d_A01_x = -C01 - C11 - C21;
			A11 = C01 * D10 + C11 * D11 + C21 * D12; d_A11_y = -C01 - C11 - C21;
			A21 = C01 * D20 + C11 * D21 + C21 * D22; d_A21_z = -C01 - C11 - C21;
			A02 = C02 * D00 + C12 * D01 + C22 * D02; d_A02_x = -C02 - C12 - C22;
			A12 = C02 * D10 + C12 * D11 + C22 * D12; d_A12_y = -C02 - C12 - C22;
			A22 = C02 * D20 + C12 * D21 + C22 * D22; d_A22_z = -C02 - C12 - C22;
			A2_00 = A00 * A00 + A10 * A10 + A20 * A20;
			A2_01 = A00 * A01 + A10 * A11 + A20 * A21;
			A2_02 = A00 * A02 + A10 * A12 + A20 * A22;
			A2_11 = A01 * A01 + A11 * A11 + A21 * A21;
			A2_12 = A01 * A02 + A11 * A12 + A21 * A22;
			A2_22 = A02 * A02 + A12 * A12 + A22 * A22;
			AF = A2_00 + A2_11 + A2_22; AF_05 = std::sqrt(AF); I_AF_05 = 0.5 / AF_05;
			AF2 = A2_00 * A2_00 + A2_01 * A2_01 * 2 + A2_02 * A2_02 * 2 + A2_11 * A2_11 + A2_12 * A2_12 * 2 + A2_22 * A2_22;
			AF_I = (AF*AF - AF2)*0.5; AF_I_05 = std::sqrt(AF_I); I_AF_I_05 = 0.5 / AF_I_05;
			det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
			i_det_A = 1.0 / det_A;
			D_AF_x = (2.0*A00*d_A00_x + 2.0*A01*d_A01_x + 2.0*A02*d_A02_x);
			D_AF_y = (2.0*A10*d_A10_y + 2.0*A11*d_A11_y + 2.0*A12*d_A12_y);
			D_AF_z = (2.0*A20*d_A20_z + 2.0*A21*d_A21_z + 2.0*A22*d_A22_z);
			D_AF2_x = (2.0*A2_00*2.0*A00*d_A00_x + 4.0*A2_01*(A01*d_A00_x + A00 * d_A01_x) + 4.0*A2_02*(A02*d_A00_x + A00 * d_A02_x)
				+ 2.0*A2_11*2.0*A01*d_A01_x + 2.0*A2_22*2.0*A02*d_A02_x + 4.0*A2_12*(A01*d_A02_x + A02 * d_A01_x));
			D_AF2_y = (2.0*A2_00*2.0*A10*d_A10_y + 4.0*A2_01*(A11*d_A10_y + A10 * d_A11_y) + 4.0*A2_02*(A12*d_A10_y + A10 * d_A12_y)
				+ 2.0*A2_11*2.0*A11*d_A11_y + 2.0*A2_22*2.0*A12*d_A12_y + 4.0*A2_12*(A11*d_A12_y + A12 * d_A11_y));
			D_AF2_z = (2.0*A2_00*2.0*A20*d_A20_z + 4.0*A2_01*(A21*d_A20_z + A20 * d_A21_z) + 4.0*A2_02*(A22*d_A20_z + A20 * d_A22_z)
				+ 2.0*A2_11*2.0*A21*d_A21_z + 2.0*A2_22*2.0*A22*d_A22_z + 4.0*A2_12*(A21*d_A22_z + A22 * d_A21_z));
			D_AF_I_x = (AF*D_AF_x - 0.5*D_AF2_x)*I_AF_I_05; D_AF_I_y = (AF*D_AF_y - 0.5*D_AF2_y)*I_AF_I_05; D_AF_I_z = (AF*D_AF_z - 0.5*D_AF2_z)*I_AF_I_05;
			D_AF_x *= I_AF_05; D_AF_y *= I_AF_05; D_AF_z *= I_AF_05;
			D_det_A_x = d_A00_x * A11 * A22 + d_A01_x * A12 * A20 + d_A02_x * A10 * A21 - d_A02_x * A11 * A20 - d_A01_x * A10 * A22 - d_A00_x * A12 * A21;
			D_det_A_y = A00 * d_A11_y * A22 + A01 * d_A12_y * A20 + A02 * d_A10_y * A21 - A02 * d_A11_y * A20 - A01 * d_A10_y * A22 - A00 * d_A12_y * A21;
			D_det_A_z = A00 * A11 * d_A22_z + A01 * A12 * d_A20_z + A02 * A10 * d_A21_z - A02 * A11 * d_A20_z - A01 * A10 * d_A22_z - A00 * A12 * d_A21_z;
			inv_det_A_2_05 = 0.5 *(1.0 - i_det_A * i_det_A);
			dvex = D_det_A_x * inv_det_A_2_05;
			dvey = D_det_A_y * inv_det_A_2_05;
			dvez = D_det_A_z * inv_det_A_2_05;
			tmp_g = AF_05 * AF_I_05;
			dgx = D_AF_x * AF_I_05 + AF_05 * D_AF_I_x;
			dgy = D_AF_y * AF_I_05 + AF_05 * D_AF_I_y;
			dgz = D_AF_z * AF_I_05 + AF_05 * D_AF_I_z;
			dex = (dgx *i_det_A - tmp_g * D_det_A_x *i_det_A*i_det_A);
			dey = (dgy *i_det_A - tmp_g * D_det_A_y *i_det_A*i_det_A);
			dez = (dgz *i_det_A - tmp_g * D_det_A_z *i_det_A*i_det_A);
			e = tmp_g * i_det_A;
			mips_e = (e*e - 1.0)*0.125;
			//mips_e = e/3.0;
			volume_e = 0.5*(det_A + i_det_A);
			exp_k = (alpha*mips_e + beta * volume_e);
			if (exp_k > 60) exp_k = 60;
			tmp_f = exp_k;
			//std::cout << (tmp_f - exp_vec[j]) / EPSL << std::endl;
			diff[2] = (tmp_f - exp_vec[j]) / EPSL;
			tmp_f = 0.0;
			z0 = x[3 * i + 2];
#endif
		}
		for (size_t j = 0; j < vc_size; j++)
		{
			double exp_e = std::exp(exp_vec[j]);
			if (cell_flag[vertex_cell[id][j]] == -1)
			{
				//not marked yet
				cell_flag[vertex_cell[id][j]] = 1;
				*f += coeff_iso * exp_e;
			}
			//derivative
			g[3 * i] += coeff_iso * gx_vec[j] * exp_e;
			g[3 * i + 1] += coeff_iso * gy_vec[j] * exp_e;
			g[3 * i + 2] += coeff_iso * gz_vec[j] * exp_e;
			
		}
	}
	//assert all cells have been colored
	/*for (size_t i = 0; i < nc; i++)
	{
		assert(cell_flag[i] == 1);
	}*/
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < nc; p++)
	{
		int tid = omp_get_thread_num();
		
		int vid[4];
		vid[0] = cell_vertex[p][0];
		vid[1] = cell_vertex[p][1];
		vid[2] = cell_vertex[p][2];
		vid[3] = cell_vertex[p][3];
		OpenVolumeMesh::Geometry::Vec3d cp_ori(x[3 * vid[0]], x[3 * vid[0] + 1], x[3 * vid[0] + 2]);
		OpenVolumeMesh::Geometry::Vec3d cr_ori(x[3 * vid[1]], x[3 * vid[1] + 1], x[3 * vid[1] + 2]);
		OpenVolumeMesh::Geometry::Vec3d cs_ori(x[3 * vid[2]], x[3 * vid[2] + 1], x[3 * vid[2] + 2]);
		OpenVolumeMesh::Geometry::Vec3d ct_ori(x[3 * vid[3]], x[3 * vid[3] + 1], x[3 * vid[3] + 2]);
		OpenVolumeMesh::Geometry::Vec3d cr = cr_ori - cp_ori;
		OpenVolumeMesh::Geometry::Vec3d cs = cs_ori - cp_ori;
		OpenVolumeMesh::Geometry::Vec3d ct = ct_ori - cp_ori;
		Eigen::Matrix3d T;
		T(0, 0) = cr[0]; T(1, 0) = cr[1]; T(2, 0) = cr[2];
		T(0, 1) = cs[0]; T(1, 1) = cs[1]; T(2, 1) = cs[2];
		T(0, 2) = ct[0]; T(1, 2) = ct[1]; T(2, 2) = ct[2];
		const Eigen::Matrix3d &invT0 = cell_S[p];
		Eigen::Matrix3d J = T * invT0;
		
		double det = J.determinant();
		if ((det) < 1.0e-10)
		{
			func = INFINITY;
			continue;
		}
		//double mips_e = (det * det - 1.0) * 0.125;
		double volume_e = 0.5 * (det + 1.0 / det);
		
		//local_f += 0.5 * (det + 1.0 / det) * (1 - alpha);
		double factor = 0.5 * (det - 1.0 / det) * (1 - alpha);
		Eigen::Matrix3d diffA;
		Eigen::Matrix3d invJ = J.inverse();
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				diffA(j, k) = factor * invJ(k, j);
		//conformal
		double c0 = 0, c1 = 0;
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
			{
				c0 += J(j, k) * J(j, k);
				c1 += invJ(j, k) * invJ(j, k);
			}
		double mips_e = 1.0 / 8.0 * (c0 * c1 - 1);
		double exp_k = alpha * mips_e + (1.0 - alpha) * volume_e;
		double exp_e = std::exp(exp_k);
		
		func += coeff_iso * exp_e;
		factor = 1.0 / 8.0 * alpha;
		Eigen::Matrix3d invJTinvJ = invJ.transpose() * invJ;
		Eigen::Matrix3d invJTinvJinvJT = invJTinvJ * invJ.transpose();
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
			{
				//Matrix3<Real> Q;
				//Q.MakeZero(); Q[j][k] = 1;
				//diffA[j][k] += factor * (c1 * Trace(Q.TransposeTimes(J) + J.TransposeTimes(Q)) -
				//	c0 * Trace(invJ.TransposeTimes(Q.TransposeTimesTranspose(invJ) + invJ*Q) * invJ)
				//	);
				double tmp0 = c1 * 2 * J(j, k);
				/*double tmp1 = invJ(0, j) * invJ(k, 0) + invJ(1, j) * invJ(k, 1) + invJ(2, j) * invJ(k, 2);
				double tmp2 = invJ(k, 0) * invJTinvJ(0, j) + invJ(k, 1) * invJTinvJ(1, j) + invJ(k, 2) * invJTinvJ(2, j);
				diffA(j, k) += factor * (tmp0 - c0 * (tmp1 + tmp2));*/
				double tmp1 = -2 * c0 * invJTinvJinvJT(j, k);
				diffA(j, k) += factor * (tmp0 + tmp1);
			}
		g_omp[tid][3 * vid[0]]	   -= coeff_iso * exp_e * (diffA(0, 0) * (invT0(0,0) + invT0(1,0) + invT0(2,0)) + diffA(0, 1) * (invT0(0, 1) + invT0(1, 1) + invT0(2, 1)) + diffA(0, 2) * (invT0(0, 2) + invT0(1, 2) + invT0(2, 2)));
		g_omp[tid][3 * vid[0] + 1] -= coeff_iso * exp_e * (diffA(1, 0) * (invT0(0,0) + invT0(1,0) + invT0(2,0)) + diffA(1, 1) * (invT0(0, 1) + invT0(1, 1) + invT0(2, 1)) + diffA(1, 2) * (invT0(0, 2) + invT0(1, 2) + invT0(2, 2)));
		g_omp[tid][3 * vid[0] + 2] -= coeff_iso * exp_e * (diffA(2, 0) * (invT0(0,0) + invT0(1,0) + invT0(2,0)) + diffA(2, 1) * (invT0(0, 1) + invT0(1, 1) + invT0(2, 1)) + diffA(2, 2) * (invT0(0, 2) + invT0(1, 2) + invT0(2, 2)));
		for (int j = 0; j < 3; j++)
		{
			g_omp[tid][3 * vid[j + 1]]     += coeff_iso * exp_e * (diffA(0, 0) * invT0(j, 0) + diffA(0, 1) * invT0(j, 1) + diffA(0, 2) * invT0(j, 2));
			g_omp[tid][3 * vid[j + 1] + 1] += coeff_iso * exp_e * (diffA(1, 0) * invT0(j, 0) + diffA(1, 1) * invT0(j, 1) + diffA(1, 2) * invT0(j, 2));
			g_omp[tid][3 * vid[j + 1] + 2] += coeff_iso * exp_e * (diffA(2, 0) * invT0(j, 0) + diffA(2, 1) * invT0(j, 1) + diffA(2, 2) * invT0(j, 2));
		}
		
#if 0
		//gradient check
		double gv0[3];
		gv0[0] = -coeff_iso * exp_e * (diffA(0, 0) * (invT0(0, 0) + invT0(1, 0) + invT0(2, 0)) + diffA(0, 1) * (invT0(0, 1) + invT0(1, 1) + invT0(2, 1)) + diffA(0, 2) * (invT0(0, 2) + invT0(1, 2) + invT0(2, 2)));
		gv0[1] = -coeff_iso * exp_e * (diffA(1, 0) * (invT0(0, 0) + invT0(1, 0) + invT0(2, 0)) + diffA(1, 1) * (invT0(0, 1) + invT0(1, 1) + invT0(2, 1)) + diffA(1, 2) * (invT0(0, 2) + invT0(1, 2) + invT0(2, 2)));
		gv0[2] = -coeff_iso * exp_e * (diffA(2, 0) * (invT0(0,0) + invT0(1,0) + invT0(2,0)) + diffA(2, 1) * (invT0(0, 1) + invT0(1, 1) + invT0(2, 1)) + diffA(2, 2) * (invT0(0, 2) + invT0(1, 2) + invT0(2, 2)));
		//the other 3
		double gv[3][3];
		for (size_t j = 0; j < 3; j++)
		{
			gv[j][0] = coeff_iso * exp_e * (diffA(0, 0) * invT0(j, 0) + diffA(0, 1) * invT0(j, 1) + diffA(0, 2) * invT0(j, 2));
			gv[j][1] = coeff_iso * exp_e * (diffA(1, 0) * invT0(j, 0) + diffA(1, 1) * invT0(j, 1) + diffA(1, 2) * invT0(j, 2));
			gv[j][2] = coeff_iso * exp_e * (diffA(2, 0) * invT0(j, 0) + diffA(2, 1) * invT0(j, 1) + diffA(2, 2) * invT0(j, 2));
		}
		
		
		
		double diffv0[3];
		Eigen::Matrix3d T_ori = T;
		if (1)
		{
			OpenVolumeMesh::Geometry::Vec3d cp = cp_ori;
			//x
			cp = cp + OpenVolumeMesh::Geometry::Vec3d(EPSL, 0.0, 0.0);
			cr = cr_ori - cp;
			cs = cs_ori - cp;
			ct = ct_ori - cp;
			T(0, 0) = cr[0]; T(1, 0) = cr[1]; T(2, 0) = cr[2];
			T(0, 1) = cs[0]; T(1, 1) = cs[1]; T(2, 1) = cs[2];
			T(0, 2) = ct[0]; T(1, 2) = ct[1]; T(2, 2) = ct[2];
			J = T * invT0;
			det = J.determinant();
			//double mips_e = (det * det - 1.0) * 0.125;
			volume_e = 0.5 * (det + 1.0 / det);
			//local_f += 0.5 * (det + 1.0 / det) * (1 - alpha);
			//Eigen::Matrix3d diffA;
			invJ = J.inverse();
			/*for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					diffA(j, k) = factor * invJ(k, j);*/
					//conformal
			c0 = 0, c1 = 0;
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
				{
					c0 += J(j, k) * J(j, k);
					c1 += invJ(j, k) * invJ(j, k);
				}
			mips_e = 1.0 / 8.0 * (c0 * c1 - 1);
			exp_k = alpha * mips_e + (1.0 - alpha) * volume_e;
			diffv0[0] = coeff_iso * (std::exp(exp_k) - exp_e) / EPSL;
			//y
			cp = cp_ori;
			cp = cp + OpenVolumeMesh::Geometry::Vec3d(0.0, EPSL, 0.0);
			cr = cr_ori - cp;
			cs = cs_ori - cp;
			ct = ct_ori - cp;
			T(0, 0) = cr[0]; T(1, 0) = cr[1]; T(2, 0) = cr[2];
			T(0, 1) = cs[0]; T(1, 1) = cs[1]; T(2, 1) = cs[2];
			T(0, 2) = ct[0]; T(1, 2) = ct[1]; T(2, 2) = ct[2];
			J = T * invT0;
			det = J.determinant();
			//double mips_e = (det * det - 1.0) * 0.125;
			volume_e = 0.5 * (det + 1.0 / det);
			//local_f += 0.5 * (det + 1.0 / det) * (1 - alpha);
			//Eigen::Matrix3d diffA;
			invJ = J.inverse();
			/*for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					diffA(j, k) = factor * invJ(k, j);*/
					//conformal
			c0 = 0, c1 = 0;
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
				{
					c0 += J(j, k) * J(j, k);
					c1 += invJ(j, k) * invJ(j, k);
				}
			mips_e = 1.0 / 8.0 * (c0 * c1 - 1);
			exp_k = alpha * mips_e + (1.0 - alpha) * volume_e;
			diffv0[1] = coeff_iso * (std::exp(exp_k) - exp_e) / EPSL;
			//z
			cp = cp_ori;
			cp = cp + OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, EPSL);
			cr = cr_ori - cp;
			cs = cs_ori - cp;
			ct = ct_ori - cp;
			T(0, 0) = cr[0]; T(1, 0) = cr[1]; T(2, 0) = cr[2];
			T(0, 1) = cs[0]; T(1, 1) = cs[1]; T(2, 1) = cs[2];
			T(0, 2) = ct[0]; T(1, 2) = ct[1]; T(2, 2) = ct[2];
			J = T * invT0;
			det = J.determinant();
			//double mips_e = (det * det - 1.0) * 0.125;
			volume_e = 0.5 * (det + 1.0 / det);
			//local_f += 0.5 * (det + 1.0 / det) * (1 - alpha);
			//Eigen::Matrix3d diffA;
			invJ = J.inverse();
			/*for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					diffA(j, k) = factor * invJ(k, j);*/
					//conformal
			c0 = 0, c1 = 0;
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
				{
					c0 += J(j, k) * J(j, k);
					c1 += invJ(j, k) * invJ(j, k);
				}
			mips_e = 1.0 / 8.0 * (c0 * c1 - 1);
			exp_k = alpha * mips_e + (1.0 - alpha) * volume_e;
			diffv0[2] = coeff_iso * (std::exp(exp_k) - exp_e) / EPSL;
			cp = cp_ori;
		}
		
		//the other 3 verts
		double diffv[3][3];
		for (size_t iter = 0; iter < 3; iter++)
		{
			for (size_t iter1 = 0; iter1 < 3; iter1++)
			{
				T = T_ori;
				T(iter1, iter) += EPSL;
				J = T * invT0;
				det = J.determinant();
				//double mips_e = (det * det - 1.0) * 0.125;
				volume_e = 0.5 * (det + 1.0 / det);
				//local_f += 0.5 * (det + 1.0 / det) * (1 - alpha);
				//Eigen::Matrix3d diffA;
				invJ = J.inverse();
				/*for (int j = 0; j < 3; j++)
					for (int k = 0; k < 3; k++)
						diffA(j, k) = factor * invJ(k, j);*/
						//conformal
				c0 = 0, c1 = 0;
				for (int j = 0; j < 3; j++)
					for (int k = 0; k < 3; k++)
					{
						c0 += J(j, k) * J(j, k);
						c1 += invJ(j, k) * invJ(j, k);
					}
				mips_e = 1.0 / 8.0 * (c0 * c1 - 1);
				exp_k = alpha * mips_e + (1.0 - alpha) * volume_e;
				diffv[iter][iter1] = coeff_iso * (std::exp(exp_k) - exp_e) / EPSL;
			}
		}
		T = T_ori;
		
#endif		
	}
	
	cell_flag.clear();
	cell_flag.resize(bfv_id.size(), -1);
	//normal part
	double lambda = lambda_array[0] / (1.0 * boundary_face_number); //uniform weight used here
	//normal part: cellwise
	//int nc = bfv_id.size();
	//nc = 0;
	int ncn = nc;
	ncn = 0;
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < ncn; p++)
	{
		if (bfv_id[p][0] < 0) continue;
		
		int tid = omp_get_thread_num();
		int v_id_0 = bfv_id[p][0];
		int v_id_1 = bfv_id[p][1];
		int v_id_2 = bfv_id[p][2];
		
		/*v_id_1 = one_vv_id[i];
		v_id_2 = one_vv_id[j];*/
		//v_id_3 = one_vv_id[k];
		assert(a2b[v_id_0] != -1);
		assert(a2b[v_id_1] != -1);
		assert(a2b[v_id_2] != -1);
		double x0, y0, z0;
		x0 = x[3 * v_id_0];
		y0 = x[3 * v_id_0 + 1];
		z0 = x[3 * v_id_0 + 2];
		double x1, x2, y1, y2, z1, z2;
		if (a2b[v_id_1] != -1)
		{
			//vid1 on the boundary
			int tmp_bid = a2b[v_id_1];
			x1 = x[3 * tmp_bid];
			y1 = x[3 * tmp_bid + 1];
			z1 = x[3 * tmp_bid + 2];
		}
		if (a2b[v_id_2] != -1)
		{
			//vid1 on the boundary
			int tmp_bid = a2b[v_id_2];
			x2 = x[3 * tmp_bid];
			y2 = x[3 * tmp_bid + 1];
			z2 = x[3 * tmp_bid + 2];
		}
		/*p_vf_pos_x[4 * i + 0] = x1; p_vf_pos_x[4 * i + 1] = x2; p_vf_pos_x[4 * i + 2] = x3;
		p_vf_pos_y[4 * i + 0] = y1; p_vf_pos_y[4 * i + 1] = y2; p_vf_pos_y[4 * i + 2] = y3;
		p_vf_pos_z[4 * i + 0] = z1; p_vf_pos_z[4 * i + 1] = z2; p_vf_pos_z[4 * i + 2] = z3;*/
		const OpenVolumeMesh::Geometry::Vec3d& tn = target_bfn[p];
		//p_vf_pos_x[4 * i + 3] = tn[0]; p_vf_pos_y[4 * i + 3] = tn[1]; p_vf_pos_z[4 * i + 3] = tn[2];
		double nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
		double ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
		double nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);
		double len2 = nx * nx + ny * ny + nz * nz;  double len = std::sqrt(len2);
		if (len < LEN_LB)
		{
			len = LEN_LB;
			len2 = len * len;
		}
		//double E_ne = 1.0 + 1.0e-8 - (nx*tn[0] + ny*tn[1] + nz*tn[2])/len;//polycube energy
		double E_ne = 0.5*((nx / len - tn[0])*(nx / len - tn[0]) + (ny / len - tn[1])*(ny / len - tn[1]) + (nz / len - tn[2])*(nz / len - tn[2]));//polycube energy
		//*f += lambda * E_ne;
		func += lambda * E_ne;
		//derivative with respect to 3 point
		double d_nx_x0 = 0.0;     double d_nx_y0 = z1 - z2; double d_nx_z0 = y2 - y1;
		double d_ny_x0 = z2 - z1; double d_ny_y0 = 0.0;     double d_ny_z0 = x1 - x2;
		double d_nz_x0 = y1 - y2; double d_nz_y0 = x2 - x1; double d_nz_z0 = 0.0;
		double d_len2_x0 = (nx*d_nx_x0 + ny * d_ny_x0 + nz * d_nz_x0); 
		double d_len2_y0 = (nx*d_nx_y0 + ny * d_ny_y0 + nz * d_nz_y0);
		double d_len2_z0 = (nx*d_nx_z0 + ny * d_ny_z0 + nz * d_nz_z0);
		double d_len_x0 = d_len2_x0 / len; double d_len_y0 = d_len2_y0 / len; double d_len_z0 = d_len2_z0 / len;
		double d_nx_len_x0 = d_nx_x0 / len - nx * d_len_x0 / (len*len);
		double d_nx_len_y0 = d_nx_y0 / len - nx * d_len_y0 / (len*len);
		double d_nx_len_z0 = d_nx_z0 / len - nx * d_len_z0 / (len*len);
		double d_ny_len_x0 = d_ny_x0 / len - ny * d_len_x0 / (len*len);
		double d_ny_len_y0 = d_ny_y0 / len - ny * d_len_y0 / (len*len);
		double d_ny_len_z0 = d_ny_z0 / len - ny * d_len_z0 / (len*len);
		double d_nz_len_x0 = d_nz_x0 / len - nz * d_len_x0 / (len*len);
		double d_nz_len_y0 = d_nz_y0 / len - nz * d_len_y0 / (len*len);
		double d_nz_len_z0 = d_nz_z0 / len - nz * d_len_z0 / (len*len);
		/*g[3 * v_id_0 + 0] += lambda * (d_nx_len_x0*(nx / len - tn[0]) + d_ny_len_x0 * (ny / len - tn[1]) + d_nz_len_x0 * (nz / len - tn[2]));
		g[3 * v_id_0 + 1] += lambda * (d_nx_len_y0*(nx / len - tn[0]) + d_ny_len_y0 * (ny / len - tn[1]) + d_nz_len_y0 * (nz / len - tn[2]));
		g[3 * v_id_0 + 2] += lambda * (d_nx_len_z0*(nx / len - tn[0]) + d_ny_len_z0 * (ny / len - tn[1]) + d_nz_len_z0 * (nz / len - tn[2]));
*/
		g_omp[tid][3 * v_id_0 + 0] += lambda * (d_nx_len_x0*(nx / len - tn[0]) + d_ny_len_x0 * (ny / len - tn[1]) + d_nz_len_x0 * (nz / len - tn[2]));
		g_omp[tid][3 * v_id_0 + 1] += lambda * (d_nx_len_y0*(nx / len - tn[0]) + d_ny_len_y0 * (ny / len - tn[1]) + d_nz_len_y0 * (nz / len - tn[2]));
		g_omp[tid][3 * v_id_0 + 2] += lambda * (d_nx_len_z0*(nx / len - tn[0]) + d_ny_len_z0 * (ny / len - tn[1]) + d_nz_len_z0 * (nz / len - tn[2]));
		double d_nx_x[2], d_nx_y[2], d_nx_z[2], d_ny_x[2], d_ny_y[2], d_ny_z[2], d_nz_x[2], d_nz_y[2], d_nz_z[2];
		
		//defination of above varaibles
		d_nx_x[0] = 0.0;		d_nx_y[0] = z2 - z0;	d_nx_z[0] = y0 - y2;
		d_ny_x[0] = z0 - z2;	d_ny_y[0] = 0.0;		d_ny_z[0] = x2 - x0;
		d_nz_x[0] = y2 - y0;	d_nz_y[0] = x0 - x2;	d_nz_z[0] = 0.0;
		d_nx_x[1] = 0.0;		d_nx_y[1] = z0 - z1;	d_nx_z[1] = y1 - y0;
		d_ny_x[1] = z1 - z0;	d_ny_y[1] = 0.0;		d_ny_z[1] = x0 - x1;
		d_nz_x[1] = y0 - y1;	d_nz_y[1] = x1 - x0;	d_nz_z[1] = 0.0;
		double d_len2_x[2];
		double d_len2_y[2];
		double d_len2_z[2];
		
		double d_len_x[2];
		double d_len_y[2];
		double d_len_z[2];
		
		double d_nx_len_x[2];
		double d_nx_len_y[2];
		double d_nx_len_z[2];
		
		double d_ny_len_x[2];
		double d_ny_len_y[2];
		double d_ny_len_z[2];
		
		double d_nz_len_x[2];
		double d_nz_len_y[2];
		double d_nz_len_z[2];
		for (size_t i = 0; i < 2; i++)
		{
			d_len2_x[i] = (nx*d_nx_x[i] + ny * d_ny_x[i] + nz * d_nz_x[i]);
			d_len2_y[i] = (nx*d_nx_y[i] + ny * d_ny_y[i] + nz * d_nz_y[i]);
			d_len2_z[i] = (nx*d_nx_z[i] + ny * d_ny_z[i] + nz * d_nz_z[i]);
			d_len_x[i] = d_len2_x[i] / len; d_len_y[i] = d_len2_y[i] / len; d_len_z[i] = d_len2_z[i] / len;
			d_nx_len_x[i] = d_nx_x[i] / len - nx * d_len_x[i] / (len*len);
			d_nx_len_y[i] = d_nx_y[i] / len - nx * d_len_y[i] / (len*len);
			d_nx_len_z[i] = d_nx_z[i] / len - nx * d_len_z[i] / (len*len);
			d_ny_len_x[i] = d_ny_x[i] / len - ny * d_len_x[i] / (len*len);
			d_ny_len_y[i] = d_ny_y[i] / len - ny * d_len_y[i] / (len*len);
			d_ny_len_z[i] = d_ny_z[i] / len - ny * d_len_z[i] / (len*len);
			d_nz_len_x[i] = d_nz_x[i] / len - nz * d_len_x[i] / (len*len);
			d_nz_len_y[i] = d_nz_y[i] / len - nz * d_len_y[i] / (len*len);
			d_nz_len_z[i] = d_nz_z[i] / len - nz * d_len_z[i] / (len*len);
		}
		/*g[3 * v_id_1 + 0] += lambda * (d_nx_len_x[0]*(nx / len - tn[0]) + d_ny_len_x[0] * (ny / len - tn[1]) + d_nz_len_x[0] * (nz / len - tn[2]));
		g[3 * v_id_1 + 1] += lambda * (d_nx_len_y[0]*(nx / len - tn[0]) + d_ny_len_y[0] * (ny / len - tn[1]) + d_nz_len_y[0] * (nz / len - tn[2]));
		g[3 * v_id_1 + 2] += lambda * (d_nx_len_z[0]*(nx / len - tn[0]) + d_ny_len_z[0] * (ny / len - tn[1]) + d_nz_len_z[0] * (nz / len - tn[2]));
		g[3 * v_id_2 + 0] += lambda * (d_nx_len_x[1] * (nx / len - tn[0]) + d_ny_len_x[1] * (ny / len - tn[1]) + d_nz_len_x[1] * (nz / len - tn[2]));
		g[3 * v_id_2 + 1] += lambda * (d_nx_len_y[1] * (nx / len - tn[0]) + d_ny_len_y[1] * (ny / len - tn[1]) + d_nz_len_y[1] * (nz / len - tn[2]));
		g[3 * v_id_2 + 2] += lambda * (d_nx_len_z[1] * (nx / len - tn[0]) + d_ny_len_z[1] * (ny / len - tn[1]) + d_nz_len_z[1] * (nz / len - tn[2]));
*/
		
		g_omp[tid][3 * v_id_1 + 0] += lambda * (d_nx_len_x[0]*(nx / len - tn[0]) + d_ny_len_x[0] * (ny / len - tn[1]) + d_nz_len_x[0] * (nz / len - tn[2]));
		g_omp[tid][3 * v_id_1 + 1] += lambda * (d_nx_len_y[0]*(nx / len - tn[0]) + d_ny_len_y[0] * (ny / len - tn[1]) + d_nz_len_y[0] * (nz / len - tn[2]));
		g_omp[tid][3 * v_id_1 + 2] += lambda * (d_nx_len_z[0]*(nx / len - tn[0]) + d_ny_len_z[0] * (ny / len - tn[1]) + d_nz_len_z[0] * (nz / len - tn[2]));
		g_omp[tid][3 * v_id_2 + 0] += lambda * (d_nx_len_x[1] * (nx / len - tn[0]) + d_ny_len_x[1] * (ny / len - tn[1]) + d_nz_len_x[1] * (nz / len - tn[2]));
		g_omp[tid][3 * v_id_2 + 1] += lambda * (d_nx_len_y[1] * (nx / len - tn[0]) + d_ny_len_y[1] * (ny / len - tn[1]) + d_nz_len_y[1] * (nz / len - tn[2]));
		g_omp[tid][3 * v_id_2 + 2] += lambda * (d_nx_len_z[1] * (nx / len - tn[0]) + d_ny_len_z[1] * (ny / len - tn[1]) + d_nz_len_z[1] * (nz / len - tn[2]));
	}
	//normal energy, normal vertical to edge
	int ncv = nc;
	//ncv = 0;
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < ncv; p++)
	{
		int tid = omp_get_thread_num();
		if (bfv_id[p][0] < 0) continue;
		int v_id_0 = bfv_id[p][0];
		int v_id_1 = bfv_id[p][1];
		int v_id_2 = bfv_id[p][2];
		/*v_id_1 = one_vv_id[i];
		v_id_2 = one_vv_id[j];*/
		//v_id_3 = one_vv_id[k];
		assert(a2b[v_id_0] != -1);
		assert(a2b[v_id_1] != -1);
		assert(a2b[v_id_2] != -1);
		double x0, y0, z0;
		x0 = x[3 * v_id_0];
		y0 = x[3 * v_id_0 + 1];
		z0 = x[3 * v_id_0 + 2];
		double x1, x2, y1, y2, z1, z2;
		if (a2b[v_id_1] != -1)
		{
			//vid1 on the boundary
			int tmp_bid = a2b[v_id_1];
			x1 = x[3 * tmp_bid];
			y1 = x[3 * tmp_bid + 1];
			z1 = x[3 * tmp_bid + 2];
		}
		if (a2b[v_id_2] != -1)
		{
			//vid1 on the boundary
			int tmp_bid = a2b[v_id_2];
			x2 = x[3 * tmp_bid];
			y2 = x[3 * tmp_bid + 1];
			z2 = x[3 * tmp_bid + 2];
		}
		/*p_vf_pos_x[4 * i + 0] = x1; p_vf_pos_x[4 * i + 1] = x2; p_vf_pos_x[4 * i + 2] = x3;
		p_vf_pos_y[4 * i + 0] = y1; p_vf_pos_y[4 * i + 1] = y2; p_vf_pos_y[4 * i + 2] = y3;
		p_vf_pos_z[4 * i + 0] = z1; p_vf_pos_z[4 * i + 1] = z2; p_vf_pos_z[4 * i + 2] = z3;*/
		const OpenVolumeMesh::Geometry::Vec3d& tn = target_bfn[p];
		//p_vf_pos_x[4 * i + 3] = tn[0]; p_vf_pos_y[4 * i + 3] = tn[1]; p_vf_pos_z[4 * i + 3] = tn[2];
		/*double nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
		double ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
		double nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);*/
		OpenVolumeMesh::Geometry::Vec3d T01(x0 - x1, y0 - y1, z0 - z1);
		OpenVolumeMesh::Geometry::Vec3d T02(x0 - x2, y0 - y2, z0 - z2);
		OpenVolumeMesh::Geometry::Vec3d T12(x1 - x2, y1 - y2, z1 - z2);
		
		double c01 = OpenVolumeMesh::Geometry::dot(T01, tn);
		double c02 = OpenVolumeMesh::Geometry::dot(T02, tn);
		double c12 = OpenVolumeMesh::Geometry::dot(T12, tn);
		double cc = c01 + c02;
		double cc1 = -c01 + c12;
		double cc2 = -c02 - c12; 
		//*f += lambda * 0.5 * (c01 * c01 + c02 * c02 + c12 * c12);
		func += lambda * 0.5 * (c01 * c01 + c02 * c02 + c12 * c12);
		/*g[3 * v_id_0 + 0] += lambda* cc  * tn[0];
		g[3 * v_id_0 + 1] += lambda* cc  * tn[1];
		g[3 * v_id_0 + 2] += lambda* cc  * tn[2];
		g[3 * v_id_1 + 0] += lambda* cc1 * tn[0];
		g[3 * v_id_1 + 1] += lambda* cc1 * tn[1];
		g[3 * v_id_1 + 2] += lambda* cc1 * tn[2];
		g[3 * v_id_2 + 0] += lambda* cc2 * tn[0];
		g[3 * v_id_2 + 1] += lambda* cc2 * tn[1];
		g[3 * v_id_2 + 2] += lambda* cc2 * tn[2];*/
		g_omp[tid][3 * v_id_0 + 0] += lambda * cc  * tn[0];
		g_omp[tid][3 * v_id_0 + 1] += lambda * cc  * tn[1];
		g_omp[tid][3 * v_id_0 + 2] += lambda * cc  * tn[2];
		g_omp[tid][3 * v_id_1 + 0] += lambda * cc1 * tn[0];
		g_omp[tid][3 * v_id_1 + 1] += lambda * cc1 * tn[1];
		g_omp[tid][3 * v_id_1 + 2] += lambda * cc1 * tn[2];
		g_omp[tid][3 * v_id_2 + 0] += lambda * cc2 * tn[0];
		g_omp[tid][3 * v_id_2 + 1] += lambda * cc2 * tn[1];
		g_omp[tid][3 * v_id_2 + 2] += lambda * cc2 * tn[2];
#if 0
		//debug
		double f_ori = lambda * 0.5 * (c01 * c01 + c02 * c02 + c12 * c12);
		double gv[3][3];
		gv[0][0] = lambda * cc  * tn[0];	gv[0][1] = lambda * cc  * tn[1];	gv[0][2] = lambda * cc  * tn[2];
		gv[1][0] = lambda * cc1  * tn[0];	gv[1][1] = lambda * cc1  * tn[1];	gv[1][2] = lambda * cc1  * tn[2];
		gv[2][0] = lambda * cc2  * tn[0];	gv[2][1] = lambda * cc2  * tn[1];	gv[2][2] = lambda * cc2  * tn[2];
		OpenVolumeMesh::Geometry::Vec3d v0_ori(x0, y0, z0);
		OpenVolumeMesh::Geometry::Vec3d v1_ori(x1, y1, z1);
		OpenVolumeMesh::Geometry::Vec3d v2_ori(x2, y2, z2);
		OpenVolumeMesh::Geometry::Vec3d v[3];
		double diff[3][3];
		
		//all vert at once
		for (size_t iter = 0; iter < 3; iter++)
		{
			for (size_t iter1 = 0; iter1 < 3; iter1++)
			{
				v[0] = v0_ori;
				v[1] = v1_ori;
				v[2] = v2_ori;
				v[iter][iter1] += EPSL;
				T01 = v[0] - v[1];
				T02 = v[0] - v[2];
				T12 = v[1] - v[2];
				c01 = OpenVolumeMesh::Geometry::dot(T01, tn);
				c02 = OpenVolumeMesh::Geometry::dot(T02, tn);
				c12 = OpenVolumeMesh::Geometry::dot(T12, tn);
				cc = c01 + c02;
				cc1 = -c01 + c12;
				cc2 = -c02 - c12;
				double f_cur = lambda * 0.5 * (c01 * c01 + c02 * c02 + c12 * c12);
				diff[iter][iter1] = (f_cur - f_ori) / EPSL;
			}
		}
		v[0] = v0_ori;
		
#endif
		
		//double E_ne = 1.0 + 1.0e-8 - (nx*tn[0] + ny*tn[1] + nz*tn[2])/len;//polycube energy
		// = 0.5*((nx / len - tn[0])*(nx / len - tn[0]) + (ny / len - tn[1])*(ny / len - tn[1]) + (nz / len - tn[2])*(nz / len - tn[2]));//polycube energy
		//define derivative
	}
	
	//feature preserving
	//double coeff_fea = lambda_array[0] / (1.0 * feature_edges.size() + LEN_LB);
	double coeff_fea = 0.5 * lambda_array[0] / (1.0 * feature_edges.size());
	int feature_edge_size = feature_edges.size();
#if 0
	//original vertion of energy
	feature_edge_size = 0; //not using feature edge
	for (size_t p = 0; p < feature_edge_size; p++)
	{
		int eid = feature_edges[p];
		int bf[2];
		assert(bef_id[eid].size() == 2);
		bf[0] = bef_id[eid][0];
		bf[1] = bef_id[eid][1];
		int v0, v1, v2(-1), v3(-1);
		v0 = feature_e2v[eid].first;
		v1 = feature_e2v[eid].second;
		for (size_t i = 0; i < 3; i++)
		{
			//set v2
			int after = (i + 1) % 3;
			for (size_t j = 0; j < 2; j++)
			{
				if (bfv_id[bf[j]][i] == v0 && bfv_id[bf[j]][after] == v1)
				{
					int afterafter = (after + 1) % 3;
					v2 = bfv_id[bf[j]][afterafter];
				}
				if (bfv_id[bf[j]][i] == v1 && bfv_id[bf[j]][after] == v0)
				{
					int afterafter = (after + 1) % 3;
					v3 = bfv_id[bf[j]][afterafter];
				}
			}
		}
		assert(v2 != -1 && v3 != -1);
		int v0_id[2] = { v0, v0 };
		int v1_id[2] = { v3, v1 };
		int v2_id[2] = { v1, v2 };
		double nx[2], ny[2], nz[2];
		double nxn[2], nyn[2], nzn[2];
		double x1[2], x2[2], y1[2], y2[2], z1[2], z2[2];
		double x0[2], y0[2], z0[2];
		double len2[2], len[2];
		double dx_n1n2, dy_n1n2, dz_n1n2;
		//derivative to v0
		for (size_t j = 0; j < 2; j++)
		{
			x0[j] = x[3 * a2b[v0_id[j]] + 0];
			y0[j] = x[3 * a2b[v0_id[j]] + 1];
			z0[j] = x[3 * a2b[v0_id[j]] + 2];
			x1[j] = x[3 * a2b[v1_id[j]] + 0];
			y1[j] = x[3 * a2b[v1_id[j]] + 1];
			z1[j] = x[3 * a2b[v1_id[j]] + 2];
			x2[j] = x[3 * a2b[v2_id[j]] + 0];
			y2[j] = x[3 * a2b[v2_id[j]] + 1];
			z2[j] = x[3 * a2b[v2_id[j]] + 2];
			nx[j] = (y0[j] - y1[j])*(z0[j] - z2[j]) - (y0[j] - y2[j])*(z0[j] - z1[j]);
			ny[j] = (x0[j] - x2[j])*(z0[j] - z1[j]) - (x0[j] - x1[j])*(z0[j] - z2[j]);
			nz[j] = (x0[j] - x1[j])*(y0[j] - y2[j]) - (x0[j] - x2[j])*(y0[j] - y1[j]);
			len2[j] = nx[j] * nx[j] + ny[j] * ny[j] + nz[j] * nz[j];
			len[j] = std::sqrt(len2[j]);
			nxn[j] = nx[j] / len[j];
			nyn[j] = ny[j] / len[j];
			nzn[j] = nz[j] / len[j];
		}
		double n1n2 = nxn[0] * nxn[1] + nyn[0] * nyn[1] + nzn[0] * nzn[1];
		//energy term below
		double n1n22 = n1n2 * n1n2;
		//normal_e += tmp_mu * n1n22;
		*f += coeff_fea * n1n22;
		double dx_nxn[2], dy_nyn[2], dz_nzn[2], dy_nxn[2], dz_nxn[2], dx_nyn[2], dz_nyn[2], dx_nzn[2], dy_nzn[2];
		//calculate these item
		for (size_t j = 0; j < 2; j++)
		{
			//coordinates of 3 points, first dimension index, second dimension index
			double tmp_coord[3][3] = {
				{x0[j], y0[j], z0[j]},
				{x1[j], y1[j], z1[j]},
				{x2[j], y2[j], z2[j]},
			};
			double normal_coord[3] = { nx[j], ny[j], nz[j] };
			//der v0
			double der_normal[3][3] = {  //derivative corresponds to coordinates, first dimension: normal axis, second dimension: derivative axis
				{0.0, z1[j] - z2[j], y2[j] - y1[j] },
				{z2[j] - z1[j], 0.0, x1[j] - x2[j] },
				{y1[j] - y2[j], x2[j] - x1[j], 0.0 },
			};
			double der_len[3];
			for (size_t a = 0; a < 3; a++)
			{
				der_len[a] = (der_normal[0][a] * normal_coord[0] + der_normal[1][a] * normal_coord[1] + der_normal[2][a] * normal_coord[2]) / len[j];
			}
			double der_normaln[3][3];
			for (size_t a = 0; a < 3; a++)
			{
				for (size_t b = 0; b < 3; b++)
				{
					der_normaln[a][b] = (der_normal[a][b] * len[j] - normal_coord[a] * der_len[b]) / len2[j];
					//der_normaln[a][b] = 0.f;
				}
			}
			dx_nxn[j] = der_normaln[0][0]; dy_nxn[j] = der_normaln[0][1]; dz_nxn[j] = der_normaln[0][2];
			dx_nyn[j] = der_normaln[1][0]; dy_nyn[j] = der_normaln[1][1]; dz_nyn[j] = der_normaln[1][2];
			dx_nzn[j] = der_normaln[2][0]; dy_nzn[j] = der_normaln[2][1]; dz_nzn[j] = der_normaln[2][2];
		}
		//derivative to be calculated
		dx_n1n2 = nyn[0] * dx_nyn[1] + dx_nyn[0] * nyn[1] + nzn[0] * dx_nzn[1] + dx_nzn[0] * nzn[1] + nxn[0] * dx_nxn[1] + dx_nxn[0] * nxn[1];
		dy_n1n2 = nxn[0] * dy_nxn[1] + dy_nxn[0] * nxn[1] + nzn[0] * dy_nzn[1] + dy_nzn[0] * nzn[1] + nyn[0] * dy_nyn[1] + dy_nyn[0] * nyn[1];
		dz_n1n2 = nxn[0] * dz_nxn[1] + dz_nxn[0] * nxn[1] + nyn[0] * dz_nyn[1] + dz_nyn[0] * nyn[1] + nzn[0] * dz_nzn[1] + dz_nzn[0] * nzn[1];
		g[3 * a2b[v0] + 0] += 2 * coeff_fea * n1n2 * dx_n1n2;
		g[3 * a2b[v0] + 1] += 2 * coeff_fea * n1n2 * dy_n1n2;
		g[3 * a2b[v0] + 2] += 2 * coeff_fea * n1n2 * dz_n1n2;
		//derivative of v1
		v0_id[0] = v1;
		v0_id[1] = v1;
		v1_id[0] = v0;
		v1_id[1] = v2;
		v2_id[0] = v3;
		v2_id[1] = v0;
		double len_v1[2];
		for (size_t j = 0; j < 2; j++)
		{
			x0[j] = x[3 * a2b[v0_id[j]] + 0];
			y0[j] = x[3 * a2b[v0_id[j]] + 1];
			z0[j] = x[3 * a2b[v0_id[j]] + 2];
			x1[j] = x[3 * a2b[v1_id[j]] + 0];
			y1[j] = x[3 * a2b[v1_id[j]] + 1];
			z1[j] = x[3 * a2b[v1_id[j]] + 2];
			x2[j] = x[3 * a2b[v2_id[j]] + 0];
			y2[j] = x[3 * a2b[v2_id[j]] + 1];
			z2[j] = x[3 * a2b[v2_id[j]] + 2];
			/*nx[j] = (y0[j] - y1[j])*(z0[j] - z2[j]) - (y0[j] - y2[j])*(z0[j] - z1[j]);
			ny[j] = (x0[j] - x2[j])*(z0[j] - z1[j]) - (x0[j] - x1[j])*(z0[j] - z2[j]);
			nz[j] = (x0[j] - x1[j])*(y0[j] - y2[j]) - (x0[j] - x2[j])*(y0[j] - y1[j]);
			len2[j] = nx[j] * nx[j] + ny[j] * ny[j] + nz[j] * nz[j];
			len_v1[j] = std::sqrt(len2[j]);
			nxn[j] = nx[j] / len_v1[j];
			nyn[j] = ny[j] / len_v1[j];
			nzn[j] = nz[j] / len_v1[j];*/
		}
		//assert(abs(len[0] - len_v1[0]) < 0.00001 && abs(len[1] - len_v1[1]) < 0.00001);
		for (size_t j = 0; j < 2; j++)
		{
			//coordinates of 3 points, first dimension index, second dimension index
			double tmp_coord[3][3] = {
				{x0[j], y0[j], z0[j]},
				{x1[j], y1[j], z1[j]},
				{x2[j], y2[j], z2[j]},
			};
			double normal_coord[3] = { nx[j], ny[j], nz[j] };
			//der v0
			double der_normal[3][3] = {  //derivative corresponds to coordinates, first dimension: normal axis, second dimension: derivative axis
				{0.0, z1[j] - z2[j], y2[j] - y1[j] },
				{z2[j] - z1[j], 0.0, x1[j] - x2[j] },
				{y1[j] - y2[j], x2[j] - x1[j], 0.0 },
			};
			double der_len[3];
			for (size_t a = 0; a < 3; a++)
			{
				der_len[a] = (der_normal[0][a] + der_normal[1][a] + der_normal[2][a]) / len[j];
			}
			double der_normaln[3][3];
			for (size_t a = 0; a < 3; a++)
			{
				for (size_t b = 0; b < 3; b++)
				{
					der_normaln[a][b] = (der_normal[a][b] * len[j] - normal_coord[a] * der_len[b]) / len2[j];
					//der_normaln[a][b] = 0.f;
				}
			}
			dx_nxn[j] = der_normaln[0][0]; dy_nxn[j] = der_normaln[0][1]; dz_nxn[j] = der_normaln[0][2];
			dx_nyn[j] = der_normaln[1][0]; dy_nyn[j] = der_normaln[1][1]; dz_nyn[j] = der_normaln[1][2];
			dx_nzn[j] = der_normaln[2][0]; dy_nzn[j] = der_normaln[2][1]; dz_nzn[j] = der_normaln[2][2];
		}
		//derivative to be calculated
		dx_n1n2 = nyn[0] * dx_nyn[1] + dx_nyn[0] * nyn[1] + nzn[0] * dx_nzn[1] + dx_nzn[0] * nzn[1] + nxn[0] * dx_nxn[1] + dx_nxn[0] * nxn[1];
		dy_n1n2 = nxn[0] * dy_nxn[1] + dy_nxn[0] * nxn[1] + nzn[0] * dy_nzn[1] + dy_nzn[0] * nzn[1] + nyn[0] * dy_nyn[1] + dy_nyn[0] * nyn[1];
		dz_n1n2 = nxn[0] * dz_nxn[1] + dz_nxn[0] * nxn[1] + nyn[0] * dz_nyn[1] + dz_nyn[0] * nyn[1] + nzn[0] * dz_nzn[1] + dz_nzn[0] * nzn[1];
		g[3 * a2b[v1] + 0] += 2 * coeff_fea * n1n2 * dx_n1n2;
		g[3 * a2b[v1] + 1] += 2 * coeff_fea * n1n2 * dy_n1n2;
		g[3 * a2b[v1] + 2] += 2 * coeff_fea * n1n2 * dz_n1n2;
		//derivative of v2
		v0_id[1] = v2;
		v1_id[1] = v0;
		v2_id[1] = v1;
		for (size_t j = 1; j < 2; j++)
		{
			x0[j] = x[3 * a2b[v0_id[j]] + 0];
			y0[j] = x[3 * a2b[v0_id[j]] + 1];
			z0[j] = x[3 * a2b[v0_id[j]] + 2];
			x1[j] = x[3 * a2b[v1_id[j]] + 0];
			y1[j] = x[3 * a2b[v1_id[j]] + 1];
			z1[j] = x[3 * a2b[v1_id[j]] + 2];
			x2[j] = x[3 * a2b[v2_id[j]] + 0];
			y2[j] = x[3 * a2b[v2_id[j]] + 1];
			z2[j] = x[3 * a2b[v2_id[j]] + 2];
		}
		for (size_t j = 1; j < 2; j++)
		{
			//coordinates of 3 points, first dimension index, second dimension index
			double tmp_coord[3][3] = {
				{x0[j], y0[j], z0[j]},
				{x1[j], y1[j], z1[j]},
				{x2[j], y2[j], z2[j]},
			};
			double normal_coord[3] = { nx[j], ny[j], nz[j] };
			//der v0
			double der_normal[3][3] = {  //derivative corresponds to coordinates, first dimension: normal axis, second dimension: derivative axis
				{0.0, z1[j] - z2[j], y2[j] - y1[j] },
				{z2[j] - z1[j], 0.0, x1[j] - x2[j] },
				{y1[j] - y2[j], x2[j] - x1[j], 0.0 },
			};
			double der_len[3];
			for (size_t a = 0; a < 3; a++)
			{
				der_len[a] = (der_normal[0][a] + der_normal[1][a] + der_normal[2][a]) / len[j];
			}
			double der_normaln[3][3];
			for (size_t a = 0; a < 3; a++)
			{
				for (size_t b = 0; b < 3; b++)
				{
					der_normaln[a][b] = (der_normal[a][b] * len[j] - normal_coord[a] * der_len[b]) / len2[j];
					//der_normaln[a][b] = 0.f;
				}
			}
			dx_nxn[j] = der_normaln[0][0]; dy_nxn[j] = der_normaln[0][1]; dz_nxn[j] = der_normaln[0][2];
			dx_nyn[j] = der_normaln[1][0]; dy_nyn[j] = der_normaln[1][1]; dz_nyn[j] = der_normaln[1][2];
			dx_nzn[j] = der_normaln[2][0]; dy_nzn[j] = der_normaln[2][1]; dz_nzn[j] = der_normaln[2][2];
		}
		dx_n1n2 = nxn[0] * dx_nxn[1] + nyn[0] * dx_nyn[1] + nzn[0] * dx_nzn[1];
		dy_n1n2 = nxn[0] * dy_nxn[1] + nyn[0] * dy_nyn[1] + nzn[0] * dy_nzn[1];
		dz_n1n2 = nxn[0] * dz_nxn[1] + nyn[0] * dz_nyn[1] + nzn[0] * dz_nzn[1];
		g[3 * a2b[v2] + 0] += 2 * coeff_fea * n1n2 * dx_n1n2;
		g[3 * a2b[v2] + 1] += 2 * coeff_fea * n1n2 * dy_n1n2;
		g[3 * a2b[v2] + 2] += 2 * coeff_fea * n1n2 * dz_n1n2;
		//derivative of v3
		v0_id[0] = v3;
		v1_id[0] = v1;
		v2_id[0] = v0;
		for (size_t j = 0; j < 1; j++)
		{
			x0[j] = x[3 * a2b[v0_id[j]] + 0];
			y0[j] = x[3 * a2b[v0_id[j]] + 1];
			z0[j] = x[3 * a2b[v0_id[j]] + 2];
			x1[j] = x[3 * a2b[v1_id[j]] + 0];
			y1[j] = x[3 * a2b[v1_id[j]] + 1];
			z1[j] = x[3 * a2b[v1_id[j]] + 2];
			x2[j] = x[3 * a2b[v2_id[j]] + 0];
			y2[j] = x[3 * a2b[v2_id[j]] + 1];
			z2[j] = x[3 * a2b[v2_id[j]] + 2];
		}
		for (size_t j = 0; j < 1; j++)
		{
			//coordinates of 3 points, first dimension index, second dimension index
			double tmp_coord[3][3] = {
				{x0[j], y0[j], z0[j]},
				{x1[j], y1[j], z1[j]},
				{x2[j], y2[j], z2[j]},
			};
			double normal_coord[3] = { nx[j], ny[j], nz[j] };
			//der v0
			double der_normal[3][3] = {  //derivative corresponds to coordinates, first dimension: normal axis, second dimension: derivative axis
				{0.0, z1[j] - z2[j], y2[j] - y1[j] },
				{z2[j] - z1[j], 0.0, x1[j] - x2[j] },
				{y1[j] - y2[j], x2[j] - x1[j], 0.0 },
			};
			double der_len[3];
			for (size_t a = 0; a < 3; a++)
			{
				der_len[a] = (der_normal[0][a] + der_normal[1][a] + der_normal[2][a]) / len[j];
			}
			double der_normaln[3][3];
			for (size_t a = 0; a < 3; a++)
			{
				for (size_t b = 0; b < 3; b++)
				{
					der_normaln[a][b] = (der_normal[a][b] * len[j] - normal_coord[a] * der_len[b]) / len2[j];
					//der_normaln[a][b] = 0.f;
				}
			}
			dx_nxn[j] = der_normaln[0][0]; dy_nxn[j] = der_normaln[0][1]; dz_nxn[j] = der_normaln[0][2];
			dx_nyn[j] = der_normaln[1][0]; dy_nyn[j] = der_normaln[1][1]; dz_nyn[j] = der_normaln[1][2];
			dx_nzn[j] = der_normaln[2][0]; dy_nzn[j] = der_normaln[2][1]; dz_nzn[j] = der_normaln[2][2];
		}
		dx_n1n2 = nxn[1] * dx_nxn[0] + nyn[1] * dx_nyn[0] + nzn[1] * dx_nzn[0];
		dy_n1n2 = nxn[1] * dy_nxn[0] + nyn[1] * dy_nyn[0] + nzn[1] * dy_nzn[0];
		dz_n1n2 = nxn[1] * dz_nxn[0] + nyn[1] * dz_nyn[0] + nzn[1] * dz_nzn[0];
		g[3 * a2b[v3] + 0] += 2 * coeff_fea * n1n2 * dx_n1n2;
		g[3 * a2b[v3] + 1] += 2 * coeff_fea * n1n2 * dy_n1n2;
		g[3 * a2b[v3] + 2] += 2 * coeff_fea * n1n2 * dz_n1n2;
	}
#endif
	//feature_edge_size = 0;
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < feature_edge_size; p++)
	{
		int tid = omp_get_thread_num();
		int eid = feature_edges[p];
		int bf[2];
		assert(bef_id[eid].size() == 2);
		bf[0] = bef_id[eid][0];
		bf[1] = bef_id[eid][1];
		int v0, v1, v2(-1), v3(-1);
		v0 = feature_e2v[eid].first;
		v1 = feature_e2v[eid].second;
		for (size_t i = 0; i < 3; i++)
		{
			//set v2
			int after = (i + 1) % 3;
			for (size_t j = 0; j < 2; j++)
			{
				if (bfv_id[bf[j]][i] == v0 && bfv_id[bf[j]][after] == v1)
				{
					int afterafter = (after + 1) % 3;
					v2 = bfv_id[bf[j]][afterafter];
				}
				if (bfv_id[bf[j]][i] == v1 && bfv_id[bf[j]][after] == v0)
				{
					int afterafter = (after + 1) % 3;
					v3 = bfv_id[bf[j]][afterafter];
				}
			}
		}
		assert(v2 != -1 && v3 != -1);
		//two face test
		std::set<int> allverts;
		for (size_t i = 0; i < 2; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				allverts.insert(bfv_id[bf[i]][j]);
			}
		}
		assert(allverts.size() == 4);
		int v0_id[2] = { v0, v0 };
		int v1_id[2] = { v3, v1 };
		int v2_id[2] = { v1, v2 };
		double nx[2], ny[2], nz[2];
		double nxn[2], nyn[2], nzn[2];
		double x1[2], x2[2], y1[2], y2[2], z1[2], z2[2];
		double x0[2], y0[2], z0[2];
		double len2[2], len[2];
		//double dx_n1n2, dy_n1n2, dz_n1n2;
		double dnx[2][3][3], dny[2][3][3], dnz[2][3][3];
		//derivative to v0
		for (size_t j = 0; j < 2; j++)
		{
			x0[j] = x[3 * a2b[v0_id[j]] + 0];
			y0[j] = x[3 * a2b[v0_id[j]] + 1];
			z0[j] = x[3 * a2b[v0_id[j]] + 2];
			x1[j] = x[3 * a2b[v1_id[j]] + 0];
			y1[j] = x[3 * a2b[v1_id[j]] + 1];
			z1[j] = x[3 * a2b[v1_id[j]] + 2];
			x2[j] = x[3 * a2b[v2_id[j]] + 0];
			y2[j] = x[3 * a2b[v2_id[j]] + 1];
			z2[j] = x[3 * a2b[v2_id[j]] + 2];
			nx[j] = (y0[j] - y1[j])*(z0[j] - z2[j]) - (y0[j] - y2[j])*(z0[j] - z1[j]);
			ny[j] = (x0[j] - x2[j])*(z0[j] - z1[j]) - (x0[j] - x1[j])*(z0[j] - z2[j]);
			nz[j] = (x0[j] - x1[j])*(y0[j] - y2[j]) - (x0[j] - x2[j])*(y0[j] - y1[j]);
			dnx[j][0][0] = 0.0; //x0
			dnx[j][1][0] = 0.0;	//x1
			dnx[j][2][0] = 0.0; //x2
			dnx[j][0][1] = z1[j] - z2[j]; //y0
			dnx[j][1][1] = z2[j] - z0[j];
			dnx[j][2][1] = z0[j] - z1[j];
			dnx[j][0][2] = y2[j] - y1[j]; //z0
			dnx[j][1][2] = y0[j] - y2[j];
			dnx[j][2][2] = y1[j] - y0[j];
			dny[j][0][0] = z2[j] - z1[j]; //x0
			dny[j][1][0] = z0[j] - z2[j]; //x1
			dny[j][2][0] = z1[j] - z0[j]; //x2
			dny[j][0][1] = 0.0; //y0
			dny[j][1][1] = 0.0;
			dny[j][2][1] = 0.0;
			dny[j][0][2] = x1[j] - x2[j]; //z0
			dny[j][1][2] = x2[j] - x0[j];
			dny[j][2][2] = x0[j] - x1[j];
			dnz[j][0][0] = y1[j] - y2[j]; //x0
			dnz[j][1][0] = y2[j] - y0[j]; //x1
			dnz[j][2][0] = y0[j] - y1[j]; //x2
			dnz[j][0][1] = x2[j] - x1[j]; //y0
			dnz[j][1][1] = x0[j] - x2[j];
			dnz[j][2][1] = x1[j] - x0[j];
			dnz[j][0][2] = 0.0; //z0
			dnz[j][1][2] = 0.0;
			dnz[j][2][2] = 0.0;
			len2[j] = nx[j] * nx[j] + ny[j] * ny[j] + nz[j] * nz[j];
			len[j] = std::sqrt(len2[j]);
			if (len[j] < LEN_LB)
			{
				len[j] = LEN_LB;
				len2[j] = len[j] * len[j];
			}
			nxn[j] = nx[j] / len[j];
			nyn[j] = ny[j] / len[j];
			nzn[j] = nz[j] / len[j];
		}
		double n1n2 = nxn[0] * nxn[1] + nyn[0] * nyn[1] + nzn[0] * nzn[1];
		//energy term below
		double n1n22 = n1n2 * n1n2;
		//normal_e += tmp_mu * n1n22;
		//*f += coeff_fea * n1n22;
		func += coeff_fea * n1n22;
		double dnxn[2][3][3], dnyn[2][3][3], dnzn[2][3][3]; //first dim: face, second dim: vert, third dim: axis
		double dlen[2][3][3];
		//calculate above
		for (size_t i = 0; i < 2; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				for (size_t k = 0; k < 3; k++)
				{
					dlen[i][j][k] = (nx[i] * dnx[i][j][k] + ny[i] * dny[i][j][k] + nz[i] * dnz[i][j][k]) / len[i];
				}
			}
		}
		for (size_t i = 0; i < 2; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				for (size_t k = 0; k < 3; k++)
				{
					dnxn[i][j][k] = (dnx[i][j][k] * len[i] - nx[i] * dlen[i][j][k]) / len2[i];
					dnyn[i][j][k] = (dny[i][j][k] * len[i] - ny[i] * dlen[i][j][k]) / len2[i];
					dnzn[i][j][k] = (dnz[i][j][k] * len[i] - nz[i] * dlen[i][j][k]) / len2[i];
				}
			}
		}
		double dn1n2[4][3]; //4 verts, 3 axis
		//3 axis
		for (size_t i = 0; i < 4; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				dn1n2[i][j] = 0.0;
			}
		}
		//v0: 0 0
		for (size_t i = 0; i < 3; i++)
		{
			dn1n2[0][i] = dnxn[0][0][i] * nxn[1] + nxn[0] * dnxn[1][0][i] + dnyn[0][0][i] * nyn[1] + nyn[0] * dnyn[1][0][i] + dnzn[0][0][i] * nzn[1] + nzn[0] * dnzn[1][0][i];
		}
		//v1: 2, 1
		for (size_t i = 0; i < 3; i++)
		{
			dn1n2[1][i] = dnxn[0][2][i] * nxn[1] + nxn[0] * dnxn[1][1][i] + dnyn[0][2][i] * nyn[1] + nyn[0] * dnyn[1][1][i] + dnzn[0][2][i] * nzn[1] + nzn[0] * dnzn[1][1][i];
		}
		//v2: -1 2
		for (size_t i = 0; i < 3; i++)
		{
			dn1n2[2][i] = nxn[0] * dnxn[1][2][i] + nyn[0] * dnyn[1][2][i] + nzn[0] * dnzn[1][2][i];
		}
		//v3: 1, -1
		for (size_t i = 0; i < 3; i++)
		{
			dn1n2[3][i] = dnxn[0][1][i] * nxn[1] + dnyn[0][1][i] * nyn[1] + dnzn[0][1][i] * nzn[1];
		}
		/*g[3 * a2b[v0] + 0] += 2 * coeff_fea * n1n2 * dn1n2[0][0];
		g[3 * a2b[v0] + 1] += 2 * coeff_fea * n1n2 * dn1n2[0][1];
		g[3 * a2b[v0] + 2] += 2 * coeff_fea * n1n2 * dn1n2[0][2];
		g[3 * a2b[v1] + 0] += 2 * coeff_fea * n1n2 * dn1n2[1][0];
		g[3 * a2b[v1] + 1] += 2 * coeff_fea * n1n2 * dn1n2[1][1];
		g[3 * a2b[v1] + 2] += 2 * coeff_fea * n1n2 * dn1n2[1][2];
		g[3 * a2b[v2] + 0] += 2 * coeff_fea * n1n2 * dn1n2[2][0];
		g[3 * a2b[v2] + 1] += 2 * coeff_fea * n1n2 * dn1n2[2][1];
		g[3 * a2b[v2] + 2] += 2 * coeff_fea * n1n2 * dn1n2[2][2];
		g[3 * a2b[v3] + 0] += 2 * coeff_fea * n1n2 * dn1n2[3][0];
		g[3 * a2b[v3] + 1] += 2 * coeff_fea * n1n2 * dn1n2[3][1];
		g[3 * a2b[v3] + 2] += 2 * coeff_fea * n1n2 * dn1n2[3][2];*/
		g_omp[tid][3 * a2b[v0] + 0] += 2 * coeff_fea * n1n2 * dn1n2[0][0];
		g_omp[tid][3 * a2b[v0] + 1] += 2 * coeff_fea * n1n2 * dn1n2[0][1];
		g_omp[tid][3 * a2b[v0] + 2] += 2 * coeff_fea * n1n2 * dn1n2[0][2];
		g_omp[tid][3 * a2b[v1] + 0] += 2 * coeff_fea * n1n2 * dn1n2[1][0];
		g_omp[tid][3 * a2b[v1] + 1] += 2 * coeff_fea * n1n2 * dn1n2[1][1];
		g_omp[tid][3 * a2b[v1] + 2] += 2 * coeff_fea * n1n2 * dn1n2[1][2];
		g_omp[tid][3 * a2b[v2] + 0] += 2 * coeff_fea * n1n2 * dn1n2[2][0];
		g_omp[tid][3 * a2b[v2] + 1] += 2 * coeff_fea * n1n2 * dn1n2[2][1];
		g_omp[tid][3 * a2b[v2] + 2] += 2 * coeff_fea * n1n2 * dn1n2[2][2];
		g_omp[tid][3 * a2b[v3] + 0] += 2 * coeff_fea * n1n2 * dn1n2[3][0];
		g_omp[tid][3 * a2b[v3] + 1] += 2 * coeff_fea * n1n2 * dn1n2[3][1];
		g_omp[tid][3 * a2b[v3] + 2] += 2 * coeff_fea * n1n2 * dn1n2[3][2];
	}
	//feature normal constraint
	//double coeff_fea_normal = lambda_array[0] / (1.0 * feature_edges.size() + LEN_LB);
	double coeff_fea_normal = 0.5 * lambda_array[0] / (1.0 * feature_edges.size());
	//assert(feature_edges.size() == target_feature_edge_normal.size());
	int n_fea_normal = feature_edges.size();
	//n_fea_normal = 0;
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < n_fea_normal; p++)
	{
		int tid = omp_get_thread_num();
		int eid = feature_edges[p];
		const OpenVolumeMesh::Geometry::Vec3d target_n = target_fea_n[eid];
		int aid0 = feature_e2v[eid].first, aid1 = feature_e2v[eid].second;
		int bid0 = a2b[aid0], bid1 = a2b[aid1];
		double x1, x2, y1, y2, z1, z2;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		double n[3];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		double len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		double len = sqrt(len2);
		assert(_finite(len) == 1 && _finite(len2) == 1);
		if (len < LEN_LB)
		{
			len = LEN_LB;
			len2 = len * len;
		}
		double nn[3];
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		assert(_finite(nn[0]) && _finite(nn[1]) && _finite(nn[2]));
		double f_ori = 0;
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			func += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			f_ori += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		double dlen[2][3]; //first dimension: vert id, second dim: axis
		for (size_t i = 0; i < 3; i++)
		{
			dlen[0][i] = -n[i] / len;
			dlen[1][i] = n[i] / len;
		}
		double diff[3];
		for (size_t i = 0; i < 3; i++)
		{
			diff[i] = nn[i] - target_n[i];
		}
		double ddiff0[3][3], ddiff1[3][3]; //second dimension: axis
		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				ddiff0[i][j] = -n[i] * dlen[0][j] / len2;
				ddiff1[i][j] = -n[i] * dlen[1][j] / len2;
			}
			ddiff0[i][i] += -1.0 / len;
			ddiff1[i][i] += 1.0 / len;
		}
		double df[2][3];
		for (size_t i = 0; i < 3; i++)
		{
			df[0][i] = 0.0; //i the axis
			df[1][i] = 0.0;
			for (size_t j = 0; j < 3; j++)
			{
				df[0][i] += 2 * coeff_fea_normal * diff[j] * ddiff0[j][i];
				df[1][i] += 2 * coeff_fea_normal * diff[j] * ddiff1[j][i];
			}
		}
		/*g[3 * bid0 + 0] += df[0][0];
		g[3 * bid0 + 1] += df[0][1];
		g[3 * bid0 + 2] += df[0][2];
		g[3 * bid1 + 0] += df[1][0];
		g[3 * bid1 + 1] += df[1][1];
		g[3 * bid1 + 2] += df[1][2];*/
		g_omp[tid][3 * bid0 + 0] += df[0][0];
		g_omp[tid][3 * bid0 + 1] += df[0][1];
		g_omp[tid][3 * bid0 + 2] += df[0][2];
		g_omp[tid][3 * bid1 + 0] += df[1][0];
		g_omp[tid][3 * bid1 + 1] += df[1][1];
		g_omp[tid][3 * bid1 + 2] += df[1][2];
#if 0
		//for debugging
		std::cout << "v0 grad: " << std::endl;
		for (size_t i = 0; i < 3; i++)
		{
			std::cout << df[0][i] << std::endl;
		}
		std::cout << "v1 grad: " << std::endl;
		for (size_t i = 0; i < 3; i++)
		{
			std::cout << df[1][i] << std::endl;
		}
		//difference method
		//v0x
		double tmpf = 0;
		x1 = x[3 * bid0] + EPSL;
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 x: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][0] - (tmpf - f_ori) / EPSL) < BOUND);
		//v0y
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1] + EPSL;
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 y: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][1] - (tmpf - f_ori) / EPSL) < BOUND);
		//v0z
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2] + EPSL;
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 z: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][2] - (tmpf - f_ori) / EPSL) < BOUND);
		//v1x
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1] + EPSL;
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 x: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][0] - (tmpf - f_ori) / EPSL) < BOUND);
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1] + EPSL;
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 y: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][1] - (tmpf - f_ori) / EPSL) < BOUND);
		//v1z
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2] + EPSL;
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 z: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][2] - (tmpf - f_ori) / EPSL) < BOUND);
		tmpf = 0;
#endif
	}
	//update gradient and function
	*f += func;
#pragma omp parallel for
	for (int i = 0; i < 3 * nv; i++)
	{
		for (size_t j = 0; j < num_core; j++)
		{
			g[i] += g_omp[j][i];
		}
	}
}
void evalfunc_de_LBFGS_all_nf(int N, double* x, double *prev_x, double* f, double* g, void* user_pointer)
{
	deform_pointer *dp = (deform_pointer *)(user_pointer);
	const vector<double> &dpx = dp->pdi->dpx;
	const vector<double> &dpy = dp->pdi->dpy;
	const vector<double> &dpz = dp->pdi->dpz;
	const vector<int> &boundary_verts = dp->pdi->boundary_verts;
	const vector<OpenVolumeMesh::Geometry::Vec3i> &bfv_id = dp->pdi->bfv_id;
	const std::vector<std::vector<int>> &vertex_cell = dp->pdi->vertex_cell;
	const std::vector<std::vector<std::vector<int>>> &vertex_cell_vertex = dp->pdi->vertex_cell_vertex;
	const std::vector< OpenVolumeMesh::Geometry::Vec3d > &target_bfn = dp->pdi->target_bfn;
	const std::vector< std::vector<int> > &bvv_id = dp->pdi->bvv_id;
	const std::vector< std::vector<int> > &bvf_id = dp->pdi->bvf_id;
	const std::vector<int> &a2b = dp->pdi->a2b;
	//std::vector<std::vector<int>> &cell_vertex = dp->pdi->cell_vertex; //not using right order version, cause inverse is not in right order
	const std::vector< std::vector<std::vector<double> > > &vcv_S = dp->pdi->vcv_S;
	const std::vector< Eigen::Matrix3d > &cell_S = dp->pdi->cell_S;
	const std::vector<std::vector<int>> &cell_vertex = dp->pdi->cell_vertex;
	const std::vector<double> &lambda_array = dp->pdi->lambda_array;
	int boundary_face_number = dp->pdi->boundary_face_number;
	//int boundary_vertface_number = dp->pdi->boundary_vertface_number;
	int nc = dp->pdi->bfv_id.size();
	const vector<vector<int>> &bef_id = dp->pdi->bef_id;
	//feature edge part
	const std::vector<bool> &feature_edge_flag = dp->pdi->feature_edge_flag;
	const std::vector<std::vector<int>> &feature_v2e = dp->pdi->feature_v2e;
	const std::vector<std::vector<std::pair<int, int>>> &feature_neighbor_vert2cellpair = dp->pdi->feature_neighbor_vert2cellpair; //first face containging vert, second face the other direction
	//std::vector<std::vector<int>> feature_faces;
	const std::vector<OpenVolumeMesh::Geometry::Vec3d> &target_fea_n = dp->pdi->target_fea_n;
	//int n_feature_edges;
	const std::vector<int> &feature_edges = dp->pdi->feature_edge_array;
	const std::vector<std::pair<int, int>> &feature_e2v = dp->pdi->feature_e2v; //consistent feature edges
	const std::vector<std::vector<int>> &feature_neighbor = dp->pdi->feature_neighbor; //neighbor of edges at most 2 neighbors, cut by valence 3 vertex
	const std::vector<int> &feature_g2l = dp->pdi->feature_g2l;
	const std::vector<int> &feature_l2g = dp->pdi->feature_l2g;
	const std::vector<std::vector<int>> &feature_v2v = dp->pdi->feature_v2v;
																					   //const std::vector<OpenVolumeMesh::Geometry::Vec3d> &target_feature_edge_normal = dp->pdi->target_feature_edge_normal;
	/*int v_id_0, v_id_1, v_id_2, v_id_3, k, c_id;
	double D00, D10, D20, D01, D11, D21, D02, D12, D22;
	double C00, C01, C02, C10, C11, C12, C20, C21, C22;
	double A00, A01, A02, A10, A11, A12, A20, A21, A22;
	double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
	double AF, AF_05, I_AF_05, AF2, AF_I, AF_I_05, I_AF_I_05, exp_k, exp_e;
	double D_AF_x, D_AF_y, D_AF_z, D_AF2_x, D_AF2_y, D_AF2_z, D_AF_I_x, D_AF_I_y, D_AF_I_z, D_det_A_x, D_det_A_y, D_det_A_z, inv_det_A_2_05;
	double d_A00_x, d_A10_y, d_A20_z, d_A01_x, d_A11_y, d_A21_z, d_A02_x, d_A12_y, d_A22_z;
	double dvex, dvey, dvez, tmp_g, dgx, dgy, dgz, dex, dey, dez, e, mips_e, volume_e, d_mips_e_x, d_mips_e_y, d_mips_e_z;
	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	double len;*/
	//double min_radius = 1e30;
	double alpha = 0.5; //double beta = 0.5;
	int bv_size = boundary_verts.size();
	int nv = dp->pdi->bvf_id.size();
	assert(bv_size > 0);
	//std::vector<std::vector<double>> g_omp;
	unsigned int num_core = omp_get_max_threads();
	std::vector<std::vector<double>> g_omp(num_core, std::vector<double>(3 * nv, 0.0));
	double func = 0.0;
	*f = 0.0;
	for (size_t i = 0; i < 3 * nv; i++)
	{
		g[i] = 0.0;
	}
	//consider all cells with at least one element on the boundary
	std::vector<int> cell_flag(bfv_id.size(), -1);
	double coeff_iso = 1.0 / (1.0 * nc);
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < nc; p++)
	{
		int tid = omp_get_thread_num();
		int vid[4];
		vid[0] = cell_vertex[p][0];
		vid[1] = cell_vertex[p][1];
		vid[2] = cell_vertex[p][2];
		vid[3] = cell_vertex[p][3];
		OpenVolumeMesh::Geometry::Vec3d cp_ori(x[3 * vid[0]], x[3 * vid[0] + 1], x[3 * vid[0] + 2]);
		OpenVolumeMesh::Geometry::Vec3d cr_ori(x[3 * vid[1]], x[3 * vid[1] + 1], x[3 * vid[1] + 2]);
		OpenVolumeMesh::Geometry::Vec3d cs_ori(x[3 * vid[2]], x[3 * vid[2] + 1], x[3 * vid[2] + 2]);
		OpenVolumeMesh::Geometry::Vec3d ct_ori(x[3 * vid[3]], x[3 * vid[3] + 1], x[3 * vid[3] + 2]);
		OpenVolumeMesh::Geometry::Vec3d cr = cr_ori - cp_ori;
		OpenVolumeMesh::Geometry::Vec3d cs = cs_ori - cp_ori;
		OpenVolumeMesh::Geometry::Vec3d ct = ct_ori - cp_ori;
		Eigen::Matrix3d T;
		T(0, 0) = cr[0]; T(1, 0) = cr[1]; T(2, 0) = cr[2];
		T(0, 1) = cs[0]; T(1, 1) = cs[1]; T(2, 1) = cs[2];
		T(0, 2) = ct[0]; T(1, 2) = ct[1]; T(2, 2) = ct[2];
		const Eigen::Matrix3d &invT0 = cell_S[p];
		Eigen::Matrix3d J = T * invT0;
		double det = J.determinant();
		if ((det) < 1.0e-10)
		{
			func = INFINITY;
			continue;
		}
		//double mips_e = (det * det - 1.0) * 0.125;
		double volume_e = 0.5 * (det + 1.0 / det);
		//local_f += 0.5 * (det + 1.0 / det) * (1 - alpha);
		double factor = 0.5 * (det - 1.0 / det) * (1 - alpha);
		Eigen::Matrix3d diffA;
		Eigen::Matrix3d invJ = J.inverse();
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				diffA(j, k) = factor * invJ(k, j);
		//conformal
		double c0 = 0, c1 = 0;
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
			{
				c0 += J(j, k) * J(j, k);
				c1 += invJ(j, k) * invJ(j, k);
			}
		double mips_e = 1.0 / 8.0 * (c0 * c1 - 1);
		double exp_k = alpha * mips_e + (1.0 - alpha) * volume_e;
		double exp_e = std::exp(exp_k);
		func += coeff_iso * exp_e;
		factor = 1.0 / 8.0 * alpha;
		Eigen::Matrix3d invJTinvJ = invJ.transpose() * invJ;
		Eigen::Matrix3d invJTinvJinvJT = invJTinvJ * invJ.transpose();
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
			{
				//Matrix3<Real> Q;
				//Q.MakeZero(); Q[j][k] = 1;
				//diffA[j][k] += factor * (c1 * Trace(Q.TransposeTimes(J) + J.TransposeTimes(Q)) -
				//	c0 * Trace(invJ.TransposeTimes(Q.TransposeTimesTranspose(invJ) + invJ*Q) * invJ)
				//	);
				double tmp0 = c1 * 2 * J(j, k);
				/*double tmp1 = invJ(0, j) * invJ(k, 0) + invJ(1, j) * invJ(k, 1) + invJ(2, j) * invJ(k, 2);
				double tmp2 = invJ(k, 0) * invJTinvJ(0, j) + invJ(k, 1) * invJTinvJ(1, j) + invJ(k, 2) * invJTinvJ(2, j);
				diffA(j, k) += factor * (tmp0 - c0 * (tmp1 + tmp2));*/
				double tmp1 = -2 * c0 * invJTinvJinvJT(j, k);
				diffA(j, k) += factor * (tmp0 + tmp1);
			}
		g_omp[tid][3 * vid[0]] -= coeff_iso * exp_e * (diffA(0, 0) * (invT0(0, 0) + invT0(1, 0) + invT0(2, 0)) + diffA(0, 1) * (invT0(0, 1) + invT0(1, 1) + invT0(2, 1)) + diffA(0, 2) * (invT0(0, 2) + invT0(1, 2) + invT0(2, 2)));
		g_omp[tid][3 * vid[0] + 1] -= coeff_iso * exp_e * (diffA(1, 0) * (invT0(0, 0) + invT0(1, 0) + invT0(2, 0)) + diffA(1, 1) * (invT0(0, 1) + invT0(1, 1) + invT0(2, 1)) + diffA(1, 2) * (invT0(0, 2) + invT0(1, 2) + invT0(2, 2)));
		g_omp[tid][3 * vid[0] + 2] -= coeff_iso * exp_e * (diffA(2, 0) * (invT0(0, 0) + invT0(1, 0) + invT0(2, 0)) + diffA(2, 1) * (invT0(0, 1) + invT0(1, 1) + invT0(2, 1)) + diffA(2, 2) * (invT0(0, 2) + invT0(1, 2) + invT0(2, 2)));
		for (int j = 0; j < 3; j++)
		{
			g_omp[tid][3 * vid[j + 1]] += coeff_iso * exp_e * (diffA(0, 0) * invT0(j, 0) + diffA(0, 1) * invT0(j, 1) + diffA(0, 2) * invT0(j, 2));
			g_omp[tid][3 * vid[j + 1] + 1] += coeff_iso * exp_e * (diffA(1, 0) * invT0(j, 0) + diffA(1, 1) * invT0(j, 1) + diffA(1, 2) * invT0(j, 2));
			g_omp[tid][3 * vid[j + 1] + 2] += coeff_iso * exp_e * (diffA(2, 0) * invT0(j, 0) + diffA(2, 1) * invT0(j, 1) + diffA(2, 2) * invT0(j, 2));
		}
#if 0
		//gradient check
		double gv0[3];
		gv0[0] = -coeff_iso * exp_e * (diffA(0, 0) * (invT0(0, 0) + invT0(1, 0) + invT0(2, 0)) + diffA(0, 1) * (invT0(0, 1) + invT0(1, 1) + invT0(2, 1)) + diffA(0, 2) * (invT0(0, 2) + invT0(1, 2) + invT0(2, 2)));
		gv0[1] = -coeff_iso * exp_e * (diffA(1, 0) * (invT0(0, 0) + invT0(1, 0) + invT0(2, 0)) + diffA(1, 1) * (invT0(0, 1) + invT0(1, 1) + invT0(2, 1)) + diffA(1, 2) * (invT0(0, 2) + invT0(1, 2) + invT0(2, 2)));
		gv0[2] = -coeff_iso * exp_e * (diffA(2, 0) * (invT0(0, 0) + invT0(1, 0) + invT0(2, 0)) + diffA(2, 1) * (invT0(0, 1) + invT0(1, 1) + invT0(2, 1)) + diffA(2, 2) * (invT0(0, 2) + invT0(1, 2) + invT0(2, 2)));
		//the other 3
		double gv[3][3];
		for (size_t j = 0; j < 3; j++)
		{
			gv[j][0] = coeff_iso * exp_e * (diffA(0, 0) * invT0(j, 0) + diffA(0, 1) * invT0(j, 1) + diffA(0, 2) * invT0(j, 2));
			gv[j][1] = coeff_iso * exp_e * (diffA(1, 0) * invT0(j, 0) + diffA(1, 1) * invT0(j, 1) + diffA(1, 2) * invT0(j, 2));
			gv[j][2] = coeff_iso * exp_e * (diffA(2, 0) * invT0(j, 0) + diffA(2, 1) * invT0(j, 1) + diffA(2, 2) * invT0(j, 2));
		}
		double diffv0[3];
		Eigen::Matrix3d T_ori = T;
		if (1)
		{
			OpenVolumeMesh::Geometry::Vec3d cp = cp_ori;
			//x
			cp = cp + OpenVolumeMesh::Geometry::Vec3d(EPSL, 0.0, 0.0);
			cr = cr_ori - cp;
			cs = cs_ori - cp;
			ct = ct_ori - cp;
			T(0, 0) = cr[0]; T(1, 0) = cr[1]; T(2, 0) = cr[2];
			T(0, 1) = cs[0]; T(1, 1) = cs[1]; T(2, 1) = cs[2];
			T(0, 2) = ct[0]; T(1, 2) = ct[1]; T(2, 2) = ct[2];
			J = T * invT0;
			det = J.determinant();
			//double mips_e = (det * det - 1.0) * 0.125;
			volume_e = 0.5 * (det + 1.0 / det);
			//local_f += 0.5 * (det + 1.0 / det) * (1 - alpha);
			//Eigen::Matrix3d diffA;
			invJ = J.inverse();
			/*for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					diffA(j, k) = factor * invJ(k, j);*/
					//conformal
			c0 = 0, c1 = 0;
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
				{
					c0 += J(j, k) * J(j, k);
					c1 += invJ(j, k) * invJ(j, k);
				}
			mips_e = 1.0 / 8.0 * (c0 * c1 - 1);
			exp_k = alpha * mips_e + (1.0 - alpha) * volume_e;
			diffv0[0] = coeff_iso * (std::exp(exp_k) - exp_e) / EPSL;
			//y
			cp = cp_ori;
			cp = cp + OpenVolumeMesh::Geometry::Vec3d(0.0, EPSL, 0.0);
			cr = cr_ori - cp;
			cs = cs_ori - cp;
			ct = ct_ori - cp;
			T(0, 0) = cr[0]; T(1, 0) = cr[1]; T(2, 0) = cr[2];
			T(0, 1) = cs[0]; T(1, 1) = cs[1]; T(2, 1) = cs[2];
			T(0, 2) = ct[0]; T(1, 2) = ct[1]; T(2, 2) = ct[2];
			J = T * invT0;
			det = J.determinant();
			//double mips_e = (det * det - 1.0) * 0.125;
			volume_e = 0.5 * (det + 1.0 / det);
			//local_f += 0.5 * (det + 1.0 / det) * (1 - alpha);
			//Eigen::Matrix3d diffA;
			invJ = J.inverse();
			/*for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					diffA(j, k) = factor * invJ(k, j);*/
					//conformal
			c0 = 0, c1 = 0;
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
				{
					c0 += J(j, k) * J(j, k);
					c1 += invJ(j, k) * invJ(j, k);
				}
			mips_e = 1.0 / 8.0 * (c0 * c1 - 1);
			exp_k = alpha * mips_e + (1.0 - alpha) * volume_e;
			diffv0[1] = coeff_iso * (std::exp(exp_k) - exp_e) / EPSL;
			//z
			cp = cp_ori;
			cp = cp + OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, EPSL);
			cr = cr_ori - cp;
			cs = cs_ori - cp;
			ct = ct_ori - cp;
			T(0, 0) = cr[0]; T(1, 0) = cr[1]; T(2, 0) = cr[2];
			T(0, 1) = cs[0]; T(1, 1) = cs[1]; T(2, 1) = cs[2];
			T(0, 2) = ct[0]; T(1, 2) = ct[1]; T(2, 2) = ct[2];
			J = T * invT0;
			det = J.determinant();
			//double mips_e = (det * det - 1.0) * 0.125;
			volume_e = 0.5 * (det + 1.0 / det);
			//local_f += 0.5 * (det + 1.0 / det) * (1 - alpha);
			//Eigen::Matrix3d diffA;
			invJ = J.inverse();
			/*for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					diffA(j, k) = factor * invJ(k, j);*/
					//conformal
			c0 = 0, c1 = 0;
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
				{
					c0 += J(j, k) * J(j, k);
					c1 += invJ(j, k) * invJ(j, k);
				}
			mips_e = 1.0 / 8.0 * (c0 * c1 - 1);
			exp_k = alpha * mips_e + (1.0 - alpha) * volume_e;
			diffv0[2] = coeff_iso * (std::exp(exp_k) - exp_e) / EPSL;
			cp = cp_ori;
		}
		//the other 3 verts
		double diffv[3][3];
		for (size_t iter = 0; iter < 3; iter++)
		{
			for (size_t iter1 = 0; iter1 < 3; iter1++)
			{
				T = T_ori;
				T(iter1, iter) += EPSL;
				J = T * invT0;
				det = J.determinant();
				//double mips_e = (det * det - 1.0) * 0.125;
				volume_e = 0.5 * (det + 1.0 / det);
				//local_f += 0.5 * (det + 1.0 / det) * (1 - alpha);
				//Eigen::Matrix3d diffA;
				invJ = J.inverse();
				/*for (int j = 0; j < 3; j++)
					for (int k = 0; k < 3; k++)
						diffA(j, k) = factor * invJ(k, j);*/
						//conformal
				c0 = 0, c1 = 0;
				for (int j = 0; j < 3; j++)
					for (int k = 0; k < 3; k++)
					{
						c0 += J(j, k) * J(j, k);
						c1 += invJ(j, k) * invJ(j, k);
					}
				mips_e = 1.0 / 8.0 * (c0 * c1 - 1);
				exp_k = alpha * mips_e + (1.0 - alpha) * volume_e;
				diffv[iter][iter1] = coeff_iso * (std::exp(exp_k) - exp_e) / EPSL;
			}
		}
		T = T_ori;
#endif		
	}
	cell_flag.clear();
	cell_flag.resize(bfv_id.size(), -1);
	//normal part
	double lambda = lambda_array[0] / (1.0 * boundary_face_number); //uniform weight used here
	//normal part: cellwise
	//int nc = bfv_id.size();
	//nc = 0;
	int ncn = nc;
	//ncn = 0;
	//normal flattening
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < ncn; p++)
	{
		if (bfv_id[p][0] < 0) continue;
		int tid = omp_get_thread_num();
		int v_id_0 = bfv_id[p][0];
		int v_id_1 = bfv_id[p][1];
		int v_id_2 = bfv_id[p][2];
		/*v_id_1 = one_vv_id[i];
		v_id_2 = one_vv_id[j];*/
		//v_id_3 = one_vv_id[k];
		assert(a2b[v_id_0] != -1);
		assert(a2b[v_id_1] != -1);
		assert(a2b[v_id_2] != -1);
		double x0, y0, z0;
		x0 = x[3 * v_id_0];
		y0 = x[3 * v_id_0 + 1];
		z0 = x[3 * v_id_0 + 2];
		double x1, x2, y1, y2, z1, z2;
		if (a2b[v_id_1] != -1)
		{
			//vid1 on the boundary
			int tmp_bid = a2b[v_id_1];
			x1 = x[3 * tmp_bid];
			y1 = x[3 * tmp_bid + 1];
			z1 = x[3 * tmp_bid + 2];
		}
		if (a2b[v_id_2] != -1)
		{
			//vid1 on the boundary
			int tmp_bid = a2b[v_id_2];
			x2 = x[3 * tmp_bid];
			y2 = x[3 * tmp_bid + 1];
			z2 = x[3 * tmp_bid + 2];
		}
		/*p_vf_pos_x[4 * i + 0] = x1; p_vf_pos_x[4 * i + 1] = x2; p_vf_pos_x[4 * i + 2] = x3;
		p_vf_pos_y[4 * i + 0] = y1; p_vf_pos_y[4 * i + 1] = y2; p_vf_pos_y[4 * i + 2] = y3;
		p_vf_pos_z[4 * i + 0] = z1; p_vf_pos_z[4 * i + 1] = z2; p_vf_pos_z[4 * i + 2] = z3;*/
		const OpenVolumeMesh::Geometry::Vec3d& tn = target_bfn[p];
		//p_vf_pos_x[4 * i + 3] = tn[0]; p_vf_pos_y[4 * i + 3] = tn[1]; p_vf_pos_z[4 * i + 3] = tn[2];
		double nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
		double ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
		double nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);
		double len2 = nx * nx + ny * ny + nz * nz;  double len = std::sqrt(len2);
		if (len < LEN_LB)
		{
			len = LEN_LB;
			len2 = len * len;
		}
		//double E_ne = 1.0 + 1.0e-8 - (nx*tn[0] + ny*tn[1] + nz*tn[2])/len;//polycube energy
		//double E_ne = 0.5*((nx / len - tn[0])*(nx / len - tn[0]) + (ny / len - tn[1])*(ny / len - tn[1]) + (nz / len - tn[2])*(nz / len - tn[2]));//polycube energy
		//*f += lambda * E_ne;
		double E_ne = 0.5 * (nx * nx * ny * ny / (len2 * len2) + ny * ny * nz * nz / (len2 * len2) + nz * nz * nx * nx / (len2 * len2));
		func += lambda * E_ne;
		//derivative with respect to 3 point
		double d_nx_x0 = 0.0;     double d_nx_y0 = z1 - z2; double d_nx_z0 = y2 - y1;
		double d_ny_x0 = z2 - z1; double d_ny_y0 = 0.0;     double d_ny_z0 = x1 - x2;
		double d_nz_x0 = y1 - y2; double d_nz_y0 = x2 - x1; double d_nz_z0 = 0.0;
		double d_len2_x0 = (nx*d_nx_x0 + ny * d_ny_x0 + nz * d_nz_x0);
		double d_len2_y0 = (nx*d_nx_y0 + ny * d_ny_y0 + nz * d_nz_y0);
		double d_len2_z0 = (nx*d_nx_z0 + ny * d_ny_z0 + nz * d_nz_z0);
		double d_len_x0 = d_len2_x0 / len; double d_len_y0 = d_len2_y0 / len; double d_len_z0 = d_len2_z0 / len;
		double d_nx_len_x0 = d_nx_x0 / len - nx * d_len_x0 / (len*len);
		double d_nx_len_y0 = d_nx_y0 / len - nx * d_len_y0 / (len*len);
		double d_nx_len_z0 = d_nx_z0 / len - nx * d_len_z0 / (len*len);
		double d_ny_len_x0 = d_ny_x0 / len - ny * d_len_x0 / (len*len);
		double d_ny_len_y0 = d_ny_y0 / len - ny * d_len_y0 / (len*len);
		double d_ny_len_z0 = d_ny_z0 / len - ny * d_len_z0 / (len*len);
		double d_nz_len_x0 = d_nz_x0 / len - nz * d_len_x0 / (len*len);
		double d_nz_len_y0 = d_nz_y0 / len - nz * d_len_y0 / (len*len);
		double d_nz_len_z0 = d_nz_z0 / len - nz * d_len_z0 / (len*len);
		/*g[3 * v_id_0 + 0] += lambda * (d_nx_len_x0*(nx / len - tn[0]) + d_ny_len_x0 * (ny / len - tn[1]) + d_nz_len_x0 * (nz / len - tn[2]));
		g[3 * v_id_0 + 1] += lambda * (d_nx_len_y0*(nx / len - tn[0]) + d_ny_len_y0 * (ny / len - tn[1]) + d_nz_len_y0 * (nz / len - tn[2]));
		g[3 * v_id_0 + 2] += lambda * (d_nx_len_z0*(nx / len - tn[0]) + d_ny_len_z0 * (ny / len - tn[1]) + d_nz_len_z0 * (nz / len - tn[2]));
*/
		/*g_omp[tid][3 * v_id_0 + 0] += lambda * (d_nx_len_x0*(nx / len - tn[0]) + d_ny_len_x0 * (ny / len - tn[1]) + d_nz_len_x0 * (nz / len - tn[2]));
		g_omp[tid][3 * v_id_0 + 1] += lambda * (d_nx_len_y0*(nx / len - tn[0]) + d_ny_len_y0 * (ny / len - tn[1]) + d_nz_len_y0 * (nz / len - tn[2]));
		g_omp[tid][3 * v_id_0 + 2] += lambda * (d_nx_len_z0*(nx / len - tn[0]) + d_ny_len_z0 * (ny / len - tn[1]) + d_nz_len_z0 * (nz / len - tn[2]));
*/
		double nxn0, nyn0, nzn0;
		nxn0 = nx / len;
		nyn0 = ny / len;
		nzn0 = nz / len;
		g_omp[tid][3 * v_id_0 + 0] += lambda * (d_nx_len_x0*(nxn0 * nyn0 * nyn0 + nxn0 * nzn0 * nzn0) + d_ny_len_x0 * (nyn0 * nzn0 * nzn0 + nyn0 * nxn0 * nxn0) + d_nz_len_x0 * (nzn0 * nxn0 * nxn0 + nzn0 * nyn0 * nyn0));
		g_omp[tid][3 * v_id_0 + 1] += lambda * (d_nx_len_y0*(nxn0 * nyn0 * nyn0 + nxn0 * nzn0 * nzn0) + d_ny_len_y0 * (nyn0 * nzn0 * nzn0 + nyn0 * nxn0 * nxn0) + d_nz_len_y0 * (nzn0 * nxn0 * nxn0 + nzn0 * nyn0 * nyn0));
		g_omp[tid][3 * v_id_0 + 2] += lambda * (d_nx_len_z0*(nxn0 * nyn0 * nyn0 + nxn0 * nzn0 * nzn0) + d_ny_len_z0 * (nyn0 * nzn0 * nzn0 + nyn0 * nxn0 * nxn0) + d_nz_len_z0 * (nzn0 * nxn0 * nxn0 + nzn0 * nyn0 * nyn0));
		double d_nx_x[2], d_nx_y[2], d_nx_z[2], d_ny_x[2], d_ny_y[2], d_ny_z[2], d_nz_x[2], d_nz_y[2], d_nz_z[2];
		//defination of above varaibles
		d_nx_x[0] = 0.0;		d_nx_y[0] = z2 - z0;	d_nx_z[0] = y0 - y2;
		d_ny_x[0] = z0 - z2;	d_ny_y[0] = 0.0;		d_ny_z[0] = x2 - x0;
		d_nz_x[0] = y2 - y0;	d_nz_y[0] = x0 - x2;	d_nz_z[0] = 0.0;
		d_nx_x[1] = 0.0;		d_nx_y[1] = z0 - z1;	d_nx_z[1] = y1 - y0;
		d_ny_x[1] = z1 - z0;	d_ny_y[1] = 0.0;		d_ny_z[1] = x0 - x1;
		d_nz_x[1] = y0 - y1;	d_nz_y[1] = x1 - x0;	d_nz_z[1] = 0.0;
		double d_len2_x[2];
		double d_len2_y[2];
		double d_len2_z[2];
		double d_len_x[2];
		double d_len_y[2];
		double d_len_z[2];
		double d_nx_len_x[2];
		double d_nx_len_y[2];
		double d_nx_len_z[2];
		double d_ny_len_x[2];
		double d_ny_len_y[2];
		double d_ny_len_z[2];
		double d_nz_len_x[2];
		double d_nz_len_y[2];
		double d_nz_len_z[2];
		for (size_t i = 0; i < 2; i++)
		{
			d_len2_x[i] = (nx*d_nx_x[i] + ny * d_ny_x[i] + nz * d_nz_x[i]);
			d_len2_y[i] = (nx*d_nx_y[i] + ny * d_ny_y[i] + nz * d_nz_y[i]);
			d_len2_z[i] = (nx*d_nx_z[i] + ny * d_ny_z[i] + nz * d_nz_z[i]);
			d_len_x[i] = d_len2_x[i] / len; d_len_y[i] = d_len2_y[i] / len; d_len_z[i] = d_len2_z[i] / len;
			d_nx_len_x[i] = d_nx_x[i] / len - nx * d_len_x[i] / (len*len);
			d_nx_len_y[i] = d_nx_y[i] / len - nx * d_len_y[i] / (len*len);
			d_nx_len_z[i] = d_nx_z[i] / len - nx * d_len_z[i] / (len*len);
			d_ny_len_x[i] = d_ny_x[i] / len - ny * d_len_x[i] / (len*len);
			d_ny_len_y[i] = d_ny_y[i] / len - ny * d_len_y[i] / (len*len);
			d_ny_len_z[i] = d_ny_z[i] / len - ny * d_len_z[i] / (len*len);
			d_nz_len_x[i] = d_nz_x[i] / len - nz * d_len_x[i] / (len*len);
			d_nz_len_y[i] = d_nz_y[i] / len - nz * d_len_y[i] / (len*len);
			d_nz_len_z[i] = d_nz_z[i] / len - nz * d_len_z[i] / (len*len);
		}
		/*g[3 * v_id_1 + 0] += lambda * (d_nx_len_x[0]*(nx / len - tn[0]) + d_ny_len_x[0] * (ny / len - tn[1]) + d_nz_len_x[0] * (nz / len - tn[2]));
		g[3 * v_id_1 + 1] += lambda * (d_nx_len_y[0]*(nx / len - tn[0]) + d_ny_len_y[0] * (ny / len - tn[1]) + d_nz_len_y[0] * (nz / len - tn[2]));
		g[3 * v_id_1 + 2] += lambda * (d_nx_len_z[0]*(nx / len - tn[0]) + d_ny_len_z[0] * (ny / len - tn[1]) + d_nz_len_z[0] * (nz / len - tn[2]));
		g[3 * v_id_2 + 0] += lambda * (d_nx_len_x[1] * (nx / len - tn[0]) + d_ny_len_x[1] * (ny / len - tn[1]) + d_nz_len_x[1] * (nz / len - tn[2]));
		g[3 * v_id_2 + 1] += lambda * (d_nx_len_y[1] * (nx / len - tn[0]) + d_ny_len_y[1] * (ny / len - tn[1]) + d_nz_len_y[1] * (nz / len - tn[2]));
		g[3 * v_id_2 + 2] += lambda * (d_nx_len_z[1] * (nx / len - tn[0]) + d_ny_len_z[1] * (ny / len - tn[1]) + d_nz_len_z[1] * (nz / len - tn[2]));
*/
		g_omp[tid][3 * v_id_1 + 0] += lambda * (d_nx_len_x[0] * (nxn0 * nyn0 * nyn0 + nxn0 * nzn0 * nzn0) + d_ny_len_x[0] * (nyn0 * nzn0 * nzn0 + nyn0 * nxn0 * nxn0) + d_nz_len_x[0] * (nzn0 * nxn0 * nxn0 + nzn0 * nyn0 * nyn0));
		g_omp[tid][3 * v_id_1 + 1] += lambda * (d_nx_len_y[0] * (nxn0 * nyn0 * nyn0 + nxn0 * nzn0 * nzn0) + d_ny_len_y[0] * (nyn0 * nzn0 * nzn0 + nyn0 * nxn0 * nxn0) + d_nz_len_y[0] * (nzn0 * nxn0 * nxn0 + nzn0 * nyn0 * nyn0));
		g_omp[tid][3 * v_id_1 + 2] += lambda * (d_nx_len_z[0] * (nxn0 * nyn0 * nyn0 + nxn0 * nzn0 * nzn0) + d_ny_len_z[0] * (nyn0 * nzn0 * nzn0 + nyn0 * nxn0 * nxn0) + d_nz_len_z[0] * (nzn0 * nxn0 * nxn0 + nzn0 * nyn0 * nyn0));
		g_omp[tid][3 * v_id_2 + 0] += lambda * (d_nx_len_x[1] * (nxn0 * nyn0 * nyn0 + nxn0 * nzn0 * nzn0) + d_ny_len_x[1] * (nyn0 * nzn0 * nzn0 + nyn0 * nxn0 * nxn0) + d_nz_len_x[1] * (nzn0 * nxn0 * nxn0 + nzn0 * nyn0 * nyn0));
		g_omp[tid][3 * v_id_2 + 1] += lambda * (d_nx_len_y[1] * (nxn0 * nyn0 * nyn0 + nxn0 * nzn0 * nzn0) + d_ny_len_y[1] * (nyn0 * nzn0 * nzn0 + nyn0 * nxn0 * nxn0) + d_nz_len_y[1] * (nzn0 * nxn0 * nxn0 + nzn0 * nyn0 * nyn0));
		g_omp[tid][3 * v_id_2 + 2] += lambda * (d_nx_len_z[1] * (nxn0 * nyn0 * nyn0 + nxn0 * nzn0 * nzn0) + d_ny_len_z[1] * (nyn0 * nzn0 * nzn0 + nyn0 * nxn0 * nxn0) + d_nz_len_z[1] * (nzn0 * nxn0 * nxn0 + nzn0 * nyn0 * nyn0));
	}
	//normal energy, normal vertical to edge
	int ncv = nc;
	ncv = 0;
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < ncv; p++)
	{
		int tid = omp_get_thread_num();
		if (bfv_id[p][0] < 0) continue;
		int v_id_0 = bfv_id[p][0];
		int v_id_1 = bfv_id[p][1];
		int v_id_2 = bfv_id[p][2];
		/*v_id_1 = one_vv_id[i];
		v_id_2 = one_vv_id[j];*/
		//v_id_3 = one_vv_id[k];
		assert(a2b[v_id_0] != -1);
		assert(a2b[v_id_1] != -1);
		assert(a2b[v_id_2] != -1);
		double x0, y0, z0;
		x0 = x[3 * v_id_0];
		y0 = x[3 * v_id_0 + 1];
		z0 = x[3 * v_id_0 + 2];
		double x1, x2, y1, y2, z1, z2;
		if (a2b[v_id_1] != -1)
		{
			//vid1 on the boundary
			int tmp_bid = a2b[v_id_1];
			x1 = x[3 * tmp_bid];
			y1 = x[3 * tmp_bid + 1];
			z1 = x[3 * tmp_bid + 2];
		}
		if (a2b[v_id_2] != -1)
		{
			//vid1 on the boundary
			int tmp_bid = a2b[v_id_2];
			x2 = x[3 * tmp_bid];
			y2 = x[3 * tmp_bid + 1];
			z2 = x[3 * tmp_bid + 2];
		}
		/*p_vf_pos_x[4 * i + 0] = x1; p_vf_pos_x[4 * i + 1] = x2; p_vf_pos_x[4 * i + 2] = x3;
		p_vf_pos_y[4 * i + 0] = y1; p_vf_pos_y[4 * i + 1] = y2; p_vf_pos_y[4 * i + 2] = y3;
		p_vf_pos_z[4 * i + 0] = z1; p_vf_pos_z[4 * i + 1] = z2; p_vf_pos_z[4 * i + 2] = z3;*/
		const OpenVolumeMesh::Geometry::Vec3d& tn = target_bfn[p];
		//p_vf_pos_x[4 * i + 3] = tn[0]; p_vf_pos_y[4 * i + 3] = tn[1]; p_vf_pos_z[4 * i + 3] = tn[2];
		/*double nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
		double ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
		double nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);*/
		OpenVolumeMesh::Geometry::Vec3d T01(x0 - x1, y0 - y1, z0 - z1);
		OpenVolumeMesh::Geometry::Vec3d T02(x0 - x2, y0 - y2, z0 - z2);
		OpenVolumeMesh::Geometry::Vec3d T12(x1 - x2, y1 - y2, z1 - z2);
		double c01 = OpenVolumeMesh::Geometry::dot(T01, tn);
		double c02 = OpenVolumeMesh::Geometry::dot(T02, tn);
		double c12 = OpenVolumeMesh::Geometry::dot(T12, tn);
		double cc = c01 + c02;
		double cc1 = -c01 + c12;
		double cc2 = -c02 - c12;
		//*f += lambda * 0.5 * (c01 * c01 + c02 * c02 + c12 * c12);
		func += lambda * 0.5 * (c01 * c01 + c02 * c02 + c12 * c12);
		/*g[3 * v_id_0 + 0] += lambda* cc  * tn[0];
		g[3 * v_id_0 + 1] += lambda* cc  * tn[1];
		g[3 * v_id_0 + 2] += lambda* cc  * tn[2];
		g[3 * v_id_1 + 0] += lambda* cc1 * tn[0];
		g[3 * v_id_1 + 1] += lambda* cc1 * tn[1];
		g[3 * v_id_1 + 2] += lambda* cc1 * tn[2];
		g[3 * v_id_2 + 0] += lambda* cc2 * tn[0];
		g[3 * v_id_2 + 1] += lambda* cc2 * tn[1];
		g[3 * v_id_2 + 2] += lambda* cc2 * tn[2];*/
		g_omp[tid][3 * v_id_0 + 0] += lambda * cc  * tn[0];
		g_omp[tid][3 * v_id_0 + 1] += lambda * cc  * tn[1];
		g_omp[tid][3 * v_id_0 + 2] += lambda * cc  * tn[2];
		g_omp[tid][3 * v_id_1 + 0] += lambda * cc1 * tn[0];
		g_omp[tid][3 * v_id_1 + 1] += lambda * cc1 * tn[1];
		g_omp[tid][3 * v_id_1 + 2] += lambda * cc1 * tn[2];
		g_omp[tid][3 * v_id_2 + 0] += lambda * cc2 * tn[0];
		g_omp[tid][3 * v_id_2 + 1] += lambda * cc2 * tn[1];
		g_omp[tid][3 * v_id_2 + 2] += lambda * cc2 * tn[2];
#if 0
		//debug
		double f_ori = lambda * 0.5 * (c01 * c01 + c02 * c02 + c12 * c12);
		double gv[3][3];
		gv[0][0] = lambda * cc  * tn[0];	gv[0][1] = lambda * cc  * tn[1];	gv[0][2] = lambda * cc  * tn[2];
		gv[1][0] = lambda * cc1  * tn[0];	gv[1][1] = lambda * cc1  * tn[1];	gv[1][2] = lambda * cc1  * tn[2];
		gv[2][0] = lambda * cc2  * tn[0];	gv[2][1] = lambda * cc2  * tn[1];	gv[2][2] = lambda * cc2  * tn[2];
		OpenVolumeMesh::Geometry::Vec3d v0_ori(x0, y0, z0);
		OpenVolumeMesh::Geometry::Vec3d v1_ori(x1, y1, z1);
		OpenVolumeMesh::Geometry::Vec3d v2_ori(x2, y2, z2);
		OpenVolumeMesh::Geometry::Vec3d v[3];
		double diff[3][3];
		//all vert at once
		for (size_t iter = 0; iter < 3; iter++)
		{
			for (size_t iter1 = 0; iter1 < 3; iter1++)
			{
				v[0] = v0_ori;
				v[1] = v1_ori;
				v[2] = v2_ori;
				v[iter][iter1] += EPSL;
				T01 = v[0] - v[1];
				T02 = v[0] - v[2];
				T12 = v[1] - v[2];
				c01 = OpenVolumeMesh::Geometry::dot(T01, tn);
				c02 = OpenVolumeMesh::Geometry::dot(T02, tn);
				c12 = OpenVolumeMesh::Geometry::dot(T12, tn);
				cc = c01 + c02;
				cc1 = -c01 + c12;
				cc2 = -c02 - c12;
				double f_cur = lambda * 0.5 * (c01 * c01 + c02 * c02 + c12 * c12);
				diff[iter][iter1] = (f_cur - f_ori) / EPSL;
			}
		}
		v[0] = v0_ori;
#endif
		//double E_ne = 1.0 + 1.0e-8 - (nx*tn[0] + ny*tn[1] + nz*tn[2])/len;//polycube energy
		// = 0.5*((nx / len - tn[0])*(nx / len - tn[0]) + (ny / len - tn[1])*(ny / len - tn[1]) + (nz / len - tn[2])*(nz / len - tn[2]));//polycube energy
		//define derivative
	}
	//feature preserving
	//double coeff_fea = lambda_array[0] / (1.0 * feature_edges.size() + LEN_LB);
	double coeff_fea = 0.5 * lambda_array[0] / (1.0 * feature_edges.size());
	int feature_edge_size = feature_edges.size();
#if 0
	//original vertion of energy
	feature_edge_size = 0; //not using feature edge
	for (size_t p = 0; p < feature_edge_size; p++)
	{
		int eid = feature_edges[p];
		int bf[2];
		assert(bef_id[eid].size() == 2);
		bf[0] = bef_id[eid][0];
		bf[1] = bef_id[eid][1];
		int v0, v1, v2(-1), v3(-1);
		v0 = feature_e2v[eid].first;
		v1 = feature_e2v[eid].second;
		for (size_t i = 0; i < 3; i++)
		{
			//set v2
			int after = (i + 1) % 3;
			for (size_t j = 0; j < 2; j++)
			{
				if (bfv_id[bf[j]][i] == v0 && bfv_id[bf[j]][after] == v1)
				{
					int afterafter = (after + 1) % 3;
					v2 = bfv_id[bf[j]][afterafter];
				}
				if (bfv_id[bf[j]][i] == v1 && bfv_id[bf[j]][after] == v0)
				{
					int afterafter = (after + 1) % 3;
					v3 = bfv_id[bf[j]][afterafter];
				}
			}
		}
		assert(v2 != -1 && v3 != -1);
		int v0_id[2] = { v0, v0 };
		int v1_id[2] = { v3, v1 };
		int v2_id[2] = { v1, v2 };
		double nx[2], ny[2], nz[2];
		double nxn[2], nyn[2], nzn[2];
		double x1[2], x2[2], y1[2], y2[2], z1[2], z2[2];
		double x0[2], y0[2], z0[2];
		double len2[2], len[2];
		double dx_n1n2, dy_n1n2, dz_n1n2;
		//derivative to v0
		for (size_t j = 0; j < 2; j++)
		{
			x0[j] = x[3 * a2b[v0_id[j]] + 0];
			y0[j] = x[3 * a2b[v0_id[j]] + 1];
			z0[j] = x[3 * a2b[v0_id[j]] + 2];
			x1[j] = x[3 * a2b[v1_id[j]] + 0];
			y1[j] = x[3 * a2b[v1_id[j]] + 1];
			z1[j] = x[3 * a2b[v1_id[j]] + 2];
			x2[j] = x[3 * a2b[v2_id[j]] + 0];
			y2[j] = x[3 * a2b[v2_id[j]] + 1];
			z2[j] = x[3 * a2b[v2_id[j]] + 2];
			nx[j] = (y0[j] - y1[j])*(z0[j] - z2[j]) - (y0[j] - y2[j])*(z0[j] - z1[j]);
			ny[j] = (x0[j] - x2[j])*(z0[j] - z1[j]) - (x0[j] - x1[j])*(z0[j] - z2[j]);
			nz[j] = (x0[j] - x1[j])*(y0[j] - y2[j]) - (x0[j] - x2[j])*(y0[j] - y1[j]);
			len2[j] = nx[j] * nx[j] + ny[j] * ny[j] + nz[j] * nz[j];
			len[j] = std::sqrt(len2[j]);
			nxn[j] = nx[j] / len[j];
			nyn[j] = ny[j] / len[j];
			nzn[j] = nz[j] / len[j];
		}
		double n1n2 = nxn[0] * nxn[1] + nyn[0] * nyn[1] + nzn[0] * nzn[1];
		//energy term below
		double n1n22 = n1n2 * n1n2;
		//normal_e += tmp_mu * n1n22;
		*f += coeff_fea * n1n22;
		double dx_nxn[2], dy_nyn[2], dz_nzn[2], dy_nxn[2], dz_nxn[2], dx_nyn[2], dz_nyn[2], dx_nzn[2], dy_nzn[2];
		//calculate these item
		for (size_t j = 0; j < 2; j++)
		{
			//coordinates of 3 points, first dimension index, second dimension index
			double tmp_coord[3][3] = {
				{x0[j], y0[j], z0[j]},
				{x1[j], y1[j], z1[j]},
				{x2[j], y2[j], z2[j]},
			};
			double normal_coord[3] = { nx[j], ny[j], nz[j] };
			//der v0
			double der_normal[3][3] = {  //derivative corresponds to coordinates, first dimension: normal axis, second dimension: derivative axis
				{0.0, z1[j] - z2[j], y2[j] - y1[j] },
				{z2[j] - z1[j], 0.0, x1[j] - x2[j] },
				{y1[j] - y2[j], x2[j] - x1[j], 0.0 },
			};
			double der_len[3];
			for (size_t a = 0; a < 3; a++)
			{
				der_len[a] = (der_normal[0][a] * normal_coord[0] + der_normal[1][a] * normal_coord[1] + der_normal[2][a] * normal_coord[2]) / len[j];
			}
			double der_normaln[3][3];
			for (size_t a = 0; a < 3; a++)
			{
				for (size_t b = 0; b < 3; b++)
				{
					der_normaln[a][b] = (der_normal[a][b] * len[j] - normal_coord[a] * der_len[b]) / len2[j];
					//der_normaln[a][b] = 0.f;
				}
			}
			dx_nxn[j] = der_normaln[0][0]; dy_nxn[j] = der_normaln[0][1]; dz_nxn[j] = der_normaln[0][2];
			dx_nyn[j] = der_normaln[1][0]; dy_nyn[j] = der_normaln[1][1]; dz_nyn[j] = der_normaln[1][2];
			dx_nzn[j] = der_normaln[2][0]; dy_nzn[j] = der_normaln[2][1]; dz_nzn[j] = der_normaln[2][2];
		}
		//derivative to be calculated
		dx_n1n2 = nyn[0] * dx_nyn[1] + dx_nyn[0] * nyn[1] + nzn[0] * dx_nzn[1] + dx_nzn[0] * nzn[1] + nxn[0] * dx_nxn[1] + dx_nxn[0] * nxn[1];
		dy_n1n2 = nxn[0] * dy_nxn[1] + dy_nxn[0] * nxn[1] + nzn[0] * dy_nzn[1] + dy_nzn[0] * nzn[1] + nyn[0] * dy_nyn[1] + dy_nyn[0] * nyn[1];
		dz_n1n2 = nxn[0] * dz_nxn[1] + dz_nxn[0] * nxn[1] + nyn[0] * dz_nyn[1] + dz_nyn[0] * nyn[1] + nzn[0] * dz_nzn[1] + dz_nzn[0] * nzn[1];
		g[3 * a2b[v0] + 0] += 2 * coeff_fea * n1n2 * dx_n1n2;
		g[3 * a2b[v0] + 1] += 2 * coeff_fea * n1n2 * dy_n1n2;
		g[3 * a2b[v0] + 2] += 2 * coeff_fea * n1n2 * dz_n1n2;
		//derivative of v1
		v0_id[0] = v1;
		v0_id[1] = v1;
		v1_id[0] = v0;
		v1_id[1] = v2;
		v2_id[0] = v3;
		v2_id[1] = v0;
		double len_v1[2];
		for (size_t j = 0; j < 2; j++)
		{
			x0[j] = x[3 * a2b[v0_id[j]] + 0];
			y0[j] = x[3 * a2b[v0_id[j]] + 1];
			z0[j] = x[3 * a2b[v0_id[j]] + 2];
			x1[j] = x[3 * a2b[v1_id[j]] + 0];
			y1[j] = x[3 * a2b[v1_id[j]] + 1];
			z1[j] = x[3 * a2b[v1_id[j]] + 2];
			x2[j] = x[3 * a2b[v2_id[j]] + 0];
			y2[j] = x[3 * a2b[v2_id[j]] + 1];
			z2[j] = x[3 * a2b[v2_id[j]] + 2];
			/*nx[j] = (y0[j] - y1[j])*(z0[j] - z2[j]) - (y0[j] - y2[j])*(z0[j] - z1[j]);
			ny[j] = (x0[j] - x2[j])*(z0[j] - z1[j]) - (x0[j] - x1[j])*(z0[j] - z2[j]);
			nz[j] = (x0[j] - x1[j])*(y0[j] - y2[j]) - (x0[j] - x2[j])*(y0[j] - y1[j]);
			len2[j] = nx[j] * nx[j] + ny[j] * ny[j] + nz[j] * nz[j];
			len_v1[j] = std::sqrt(len2[j]);
			nxn[j] = nx[j] / len_v1[j];
			nyn[j] = ny[j] / len_v1[j];
			nzn[j] = nz[j] / len_v1[j];*/
		}
		//assert(abs(len[0] - len_v1[0]) < 0.00001 && abs(len[1] - len_v1[1]) < 0.00001);
		for (size_t j = 0; j < 2; j++)
		{
			//coordinates of 3 points, first dimension index, second dimension index
			double tmp_coord[3][3] = {
				{x0[j], y0[j], z0[j]},
				{x1[j], y1[j], z1[j]},
				{x2[j], y2[j], z2[j]},
			};
			double normal_coord[3] = { nx[j], ny[j], nz[j] };
			//der v0
			double der_normal[3][3] = {  //derivative corresponds to coordinates, first dimension: normal axis, second dimension: derivative axis
				{0.0, z1[j] - z2[j], y2[j] - y1[j] },
				{z2[j] - z1[j], 0.0, x1[j] - x2[j] },
				{y1[j] - y2[j], x2[j] - x1[j], 0.0 },
			};
			double der_len[3];
			for (size_t a = 0; a < 3; a++)
			{
				der_len[a] = (der_normal[0][a] + der_normal[1][a] + der_normal[2][a]) / len[j];
			}
			double der_normaln[3][3];
			for (size_t a = 0; a < 3; a++)
			{
				for (size_t b = 0; b < 3; b++)
				{
					der_normaln[a][b] = (der_normal[a][b] * len[j] - normal_coord[a] * der_len[b]) / len2[j];
					//der_normaln[a][b] = 0.f;
				}
			}
			dx_nxn[j] = der_normaln[0][0]; dy_nxn[j] = der_normaln[0][1]; dz_nxn[j] = der_normaln[0][2];
			dx_nyn[j] = der_normaln[1][0]; dy_nyn[j] = der_normaln[1][1]; dz_nyn[j] = der_normaln[1][2];
			dx_nzn[j] = der_normaln[2][0]; dy_nzn[j] = der_normaln[2][1]; dz_nzn[j] = der_normaln[2][2];
		}
		//derivative to be calculated
		dx_n1n2 = nyn[0] * dx_nyn[1] + dx_nyn[0] * nyn[1] + nzn[0] * dx_nzn[1] + dx_nzn[0] * nzn[1] + nxn[0] * dx_nxn[1] + dx_nxn[0] * nxn[1];
		dy_n1n2 = nxn[0] * dy_nxn[1] + dy_nxn[0] * nxn[1] + nzn[0] * dy_nzn[1] + dy_nzn[0] * nzn[1] + nyn[0] * dy_nyn[1] + dy_nyn[0] * nyn[1];
		dz_n1n2 = nxn[0] * dz_nxn[1] + dz_nxn[0] * nxn[1] + nyn[0] * dz_nyn[1] + dz_nyn[0] * nyn[1] + nzn[0] * dz_nzn[1] + dz_nzn[0] * nzn[1];
		g[3 * a2b[v1] + 0] += 2 * coeff_fea * n1n2 * dx_n1n2;
		g[3 * a2b[v1] + 1] += 2 * coeff_fea * n1n2 * dy_n1n2;
		g[3 * a2b[v1] + 2] += 2 * coeff_fea * n1n2 * dz_n1n2;
		//derivative of v2
		v0_id[1] = v2;
		v1_id[1] = v0;
		v2_id[1] = v1;
		for (size_t j = 1; j < 2; j++)
		{
			x0[j] = x[3 * a2b[v0_id[j]] + 0];
			y0[j] = x[3 * a2b[v0_id[j]] + 1];
			z0[j] = x[3 * a2b[v0_id[j]] + 2];
			x1[j] = x[3 * a2b[v1_id[j]] + 0];
			y1[j] = x[3 * a2b[v1_id[j]] + 1];
			z1[j] = x[3 * a2b[v1_id[j]] + 2];
			x2[j] = x[3 * a2b[v2_id[j]] + 0];
			y2[j] = x[3 * a2b[v2_id[j]] + 1];
			z2[j] = x[3 * a2b[v2_id[j]] + 2];
		}
		for (size_t j = 1; j < 2; j++)
		{
			//coordinates of 3 points, first dimension index, second dimension index
			double tmp_coord[3][3] = {
				{x0[j], y0[j], z0[j]},
				{x1[j], y1[j], z1[j]},
				{x2[j], y2[j], z2[j]},
			};
			double normal_coord[3] = { nx[j], ny[j], nz[j] };
			//der v0
			double der_normal[3][3] = {  //derivative corresponds to coordinates, first dimension: normal axis, second dimension: derivative axis
				{0.0, z1[j] - z2[j], y2[j] - y1[j] },
				{z2[j] - z1[j], 0.0, x1[j] - x2[j] },
				{y1[j] - y2[j], x2[j] - x1[j], 0.0 },
			};
			double der_len[3];
			for (size_t a = 0; a < 3; a++)
			{
				der_len[a] = (der_normal[0][a] + der_normal[1][a] + der_normal[2][a]) / len[j];
			}
			double der_normaln[3][3];
			for (size_t a = 0; a < 3; a++)
			{
				for (size_t b = 0; b < 3; b++)
				{
					der_normaln[a][b] = (der_normal[a][b] * len[j] - normal_coord[a] * der_len[b]) / len2[j];
					//der_normaln[a][b] = 0.f;
				}
			}
			dx_nxn[j] = der_normaln[0][0]; dy_nxn[j] = der_normaln[0][1]; dz_nxn[j] = der_normaln[0][2];
			dx_nyn[j] = der_normaln[1][0]; dy_nyn[j] = der_normaln[1][1]; dz_nyn[j] = der_normaln[1][2];
			dx_nzn[j] = der_normaln[2][0]; dy_nzn[j] = der_normaln[2][1]; dz_nzn[j] = der_normaln[2][2];
		}
		dx_n1n2 = nxn[0] * dx_nxn[1] + nyn[0] * dx_nyn[1] + nzn[0] * dx_nzn[1];
		dy_n1n2 = nxn[0] * dy_nxn[1] + nyn[0] * dy_nyn[1] + nzn[0] * dy_nzn[1];
		dz_n1n2 = nxn[0] * dz_nxn[1] + nyn[0] * dz_nyn[1] + nzn[0] * dz_nzn[1];
		g[3 * a2b[v2] + 0] += 2 * coeff_fea * n1n2 * dx_n1n2;
		g[3 * a2b[v2] + 1] += 2 * coeff_fea * n1n2 * dy_n1n2;
		g[3 * a2b[v2] + 2] += 2 * coeff_fea * n1n2 * dz_n1n2;
		//derivative of v3
		v0_id[0] = v3;
		v1_id[0] = v1;
		v2_id[0] = v0;
		for (size_t j = 0; j < 1; j++)
		{
			x0[j] = x[3 * a2b[v0_id[j]] + 0];
			y0[j] = x[3 * a2b[v0_id[j]] + 1];
			z0[j] = x[3 * a2b[v0_id[j]] + 2];
			x1[j] = x[3 * a2b[v1_id[j]] + 0];
			y1[j] = x[3 * a2b[v1_id[j]] + 1];
			z1[j] = x[3 * a2b[v1_id[j]] + 2];
			x2[j] = x[3 * a2b[v2_id[j]] + 0];
			y2[j] = x[3 * a2b[v2_id[j]] + 1];
			z2[j] = x[3 * a2b[v2_id[j]] + 2];
		}
		for (size_t j = 0; j < 1; j++)
		{
			//coordinates of 3 points, first dimension index, second dimension index
			double tmp_coord[3][3] = {
				{x0[j], y0[j], z0[j]},
				{x1[j], y1[j], z1[j]},
				{x2[j], y2[j], z2[j]},
			};
			double normal_coord[3] = { nx[j], ny[j], nz[j] };
			//der v0
			double der_normal[3][3] = {  //derivative corresponds to coordinates, first dimension: normal axis, second dimension: derivative axis
				{0.0, z1[j] - z2[j], y2[j] - y1[j] },
				{z2[j] - z1[j], 0.0, x1[j] - x2[j] },
				{y1[j] - y2[j], x2[j] - x1[j], 0.0 },
			};
			double der_len[3];
			for (size_t a = 0; a < 3; a++)
			{
				der_len[a] = (der_normal[0][a] + der_normal[1][a] + der_normal[2][a]) / len[j];
			}
			double der_normaln[3][3];
			for (size_t a = 0; a < 3; a++)
			{
				for (size_t b = 0; b < 3; b++)
				{
					der_normaln[a][b] = (der_normal[a][b] * len[j] - normal_coord[a] * der_len[b]) / len2[j];
					//der_normaln[a][b] = 0.f;
				}
			}
			dx_nxn[j] = der_normaln[0][0]; dy_nxn[j] = der_normaln[0][1]; dz_nxn[j] = der_normaln[0][2];
			dx_nyn[j] = der_normaln[1][0]; dy_nyn[j] = der_normaln[1][1]; dz_nyn[j] = der_normaln[1][2];
			dx_nzn[j] = der_normaln[2][0]; dy_nzn[j] = der_normaln[2][1]; dz_nzn[j] = der_normaln[2][2];
		}
		dx_n1n2 = nxn[1] * dx_nxn[0] + nyn[1] * dx_nyn[0] + nzn[1] * dx_nzn[0];
		dy_n1n2 = nxn[1] * dy_nxn[0] + nyn[1] * dy_nyn[0] + nzn[1] * dy_nzn[0];
		dz_n1n2 = nxn[1] * dz_nxn[0] + nyn[1] * dz_nyn[0] + nzn[1] * dz_nzn[0];
		g[3 * a2b[v3] + 0] += 2 * coeff_fea * n1n2 * dx_n1n2;
		g[3 * a2b[v3] + 1] += 2 * coeff_fea * n1n2 * dy_n1n2;
		g[3 * a2b[v3] + 2] += 2 * coeff_fea * n1n2 * dz_n1n2;
	}
#endif
	//feature vertical
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < feature_edge_size; p++)
	{
		int tid = omp_get_thread_num();
		int eid = feature_edges[p];
		int bf[2];
		assert(bef_id[eid].size() == 2);
		bf[0] = bef_id[eid][0];
		bf[1] = bef_id[eid][1];
		int v0, v1, v2(-1), v3(-1);
		v0 = feature_e2v[eid].first;
		v1 = feature_e2v[eid].second;
		for (size_t i = 0; i < 3; i++)
		{
			//set v2
			int after = (i + 1) % 3;
			for (size_t j = 0; j < 2; j++)
			{
				if (bfv_id[bf[j]][i] == v0 && bfv_id[bf[j]][after] == v1)
				{
					int afterafter = (after + 1) % 3;
					v2 = bfv_id[bf[j]][afterafter];
				}
				if (bfv_id[bf[j]][i] == v1 && bfv_id[bf[j]][after] == v0)
				{
					int afterafter = (after + 1) % 3;
					v3 = bfv_id[bf[j]][afterafter];
				}
			}
		}
		assert(v2 != -1 && v3 != -1);
		//two face test
		std::set<int> allverts;
		for (size_t i = 0; i < 2; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				allverts.insert(bfv_id[bf[i]][j]);
			}
		}
		assert(allverts.size() == 4);
		int v0_id[2] = { v0, v0 };
		int v1_id[2] = { v3, v1 };
		int v2_id[2] = { v1, v2 };
		double nx[2], ny[2], nz[2];
		double nxn[2], nyn[2], nzn[2];
		double x1[2], x2[2], y1[2], y2[2], z1[2], z2[2];
		double x0[2], y0[2], z0[2];
		double len2[2], len[2];
		//double dx_n1n2, dy_n1n2, dz_n1n2;
		double dnx[2][3][3], dny[2][3][3], dnz[2][3][3];
		//derivative to v0
		for (size_t j = 0; j < 2; j++)
		{
			x0[j] = x[3 * a2b[v0_id[j]] + 0];
			y0[j] = x[3 * a2b[v0_id[j]] + 1];
			z0[j] = x[3 * a2b[v0_id[j]] + 2];
			x1[j] = x[3 * a2b[v1_id[j]] + 0];
			y1[j] = x[3 * a2b[v1_id[j]] + 1];
			z1[j] = x[3 * a2b[v1_id[j]] + 2];
			x2[j] = x[3 * a2b[v2_id[j]] + 0];
			y2[j] = x[3 * a2b[v2_id[j]] + 1];
			z2[j] = x[3 * a2b[v2_id[j]] + 2];
			nx[j] = (y0[j] - y1[j])*(z0[j] - z2[j]) - (y0[j] - y2[j])*(z0[j] - z1[j]);
			ny[j] = (x0[j] - x2[j])*(z0[j] - z1[j]) - (x0[j] - x1[j])*(z0[j] - z2[j]);
			nz[j] = (x0[j] - x1[j])*(y0[j] - y2[j]) - (x0[j] - x2[j])*(y0[j] - y1[j]);
			dnx[j][0][0] = 0.0; //x0
			dnx[j][1][0] = 0.0;	//x1
			dnx[j][2][0] = 0.0; //x2
			dnx[j][0][1] = z1[j] - z2[j]; //y0
			dnx[j][1][1] = z2[j] - z0[j];
			dnx[j][2][1] = z0[j] - z1[j];
			dnx[j][0][2] = y2[j] - y1[j]; //z0
			dnx[j][1][2] = y0[j] - y2[j];
			dnx[j][2][2] = y1[j] - y0[j];
			dny[j][0][0] = z2[j] - z1[j]; //x0
			dny[j][1][0] = z0[j] - z2[j]; //x1
			dny[j][2][0] = z1[j] - z0[j]; //x2
			dny[j][0][1] = 0.0; //y0
			dny[j][1][1] = 0.0;
			dny[j][2][1] = 0.0;
			dny[j][0][2] = x1[j] - x2[j]; //z0
			dny[j][1][2] = x2[j] - x0[j];
			dny[j][2][2] = x0[j] - x1[j];
			dnz[j][0][0] = y1[j] - y2[j]; //x0
			dnz[j][1][0] = y2[j] - y0[j]; //x1
			dnz[j][2][0] = y0[j] - y1[j]; //x2
			dnz[j][0][1] = x2[j] - x1[j]; //y0
			dnz[j][1][1] = x0[j] - x2[j];
			dnz[j][2][1] = x1[j] - x0[j];
			dnz[j][0][2] = 0.0; //z0
			dnz[j][1][2] = 0.0;
			dnz[j][2][2] = 0.0;
			len2[j] = nx[j] * nx[j] + ny[j] * ny[j] + nz[j] * nz[j];
			len[j] = std::sqrt(len2[j]);
			if (len[j] < LEN_LB)
			{
				len[j] = LEN_LB;
				len2[j] = len[j] * len[j];
			}
			nxn[j] = nx[j] / len[j];
			nyn[j] = ny[j] / len[j];
			nzn[j] = nz[j] / len[j];
		}
		double n1n2 = nxn[0] * nxn[1] + nyn[0] * nyn[1] + nzn[0] * nzn[1];
		//energy term below
		double n1n22 = n1n2 * n1n2;
		//normal_e += tmp_mu * n1n22;
		//*f += coeff_fea * n1n22;
		func += coeff_fea * n1n22;
		double dnxn[2][3][3], dnyn[2][3][3], dnzn[2][3][3]; //first dim: face, second dim: vert, third dim: axis
		double dlen[2][3][3];
		//calculate above
		for (size_t i = 0; i < 2; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				for (size_t k = 0; k < 3; k++)
				{
					dlen[i][j][k] = (nx[i] * dnx[i][j][k] + ny[i] * dny[i][j][k] + nz[i] * dnz[i][j][k]) / len[i];
				}
			}
		}
		for (size_t i = 0; i < 2; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				for (size_t k = 0; k < 3; k++)
				{
					dnxn[i][j][k] = (dnx[i][j][k] * len[i] - nx[i] * dlen[i][j][k]) / len2[i];
					dnyn[i][j][k] = (dny[i][j][k] * len[i] - ny[i] * dlen[i][j][k]) / len2[i];
					dnzn[i][j][k] = (dnz[i][j][k] * len[i] - nz[i] * dlen[i][j][k]) / len2[i];
				}
			}
		}
		double dn1n2[4][3]; //4 verts, 3 axis
		//3 axis
		for (size_t i = 0; i < 4; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				dn1n2[i][j] = 0.0;
			}
		}
		//v0: 0 0
		for (size_t i = 0; i < 3; i++)
		{
			dn1n2[0][i] = dnxn[0][0][i] * nxn[1] + nxn[0] * dnxn[1][0][i] + dnyn[0][0][i] * nyn[1] + nyn[0] * dnyn[1][0][i] + dnzn[0][0][i] * nzn[1] + nzn[0] * dnzn[1][0][i];
		}
		//v1: 2, 1
		for (size_t i = 0; i < 3; i++)
		{
			dn1n2[1][i] = dnxn[0][2][i] * nxn[1] + nxn[0] * dnxn[1][1][i] + dnyn[0][2][i] * nyn[1] + nyn[0] * dnyn[1][1][i] + dnzn[0][2][i] * nzn[1] + nzn[0] * dnzn[1][1][i];
		}
		//v2: -1 2
		for (size_t i = 0; i < 3; i++)
		{
			dn1n2[2][i] = nxn[0] * dnxn[1][2][i] + nyn[0] * dnyn[1][2][i] + nzn[0] * dnzn[1][2][i];
		}
		//v3: 1, -1
		for (size_t i = 0; i < 3; i++)
		{
			dn1n2[3][i] = dnxn[0][1][i] * nxn[1] + dnyn[0][1][i] * nyn[1] + dnzn[0][1][i] * nzn[1];
		}
		/*g[3 * a2b[v0] + 0] += 2 * coeff_fea * n1n2 * dn1n2[0][0];
		g[3 * a2b[v0] + 1] += 2 * coeff_fea * n1n2 * dn1n2[0][1];
		g[3 * a2b[v0] + 2] += 2 * coeff_fea * n1n2 * dn1n2[0][2];
		g[3 * a2b[v1] + 0] += 2 * coeff_fea * n1n2 * dn1n2[1][0];
		g[3 * a2b[v1] + 1] += 2 * coeff_fea * n1n2 * dn1n2[1][1];
		g[3 * a2b[v1] + 2] += 2 * coeff_fea * n1n2 * dn1n2[1][2];
		g[3 * a2b[v2] + 0] += 2 * coeff_fea * n1n2 * dn1n2[2][0];
		g[3 * a2b[v2] + 1] += 2 * coeff_fea * n1n2 * dn1n2[2][1];
		g[3 * a2b[v2] + 2] += 2 * coeff_fea * n1n2 * dn1n2[2][2];
		g[3 * a2b[v3] + 0] += 2 * coeff_fea * n1n2 * dn1n2[3][0];
		g[3 * a2b[v3] + 1] += 2 * coeff_fea * n1n2 * dn1n2[3][1];
		g[3 * a2b[v3] + 2] += 2 * coeff_fea * n1n2 * dn1n2[3][2];*/
		g_omp[tid][3 * a2b[v0] + 0] += 2 * coeff_fea * n1n2 * dn1n2[0][0];
		g_omp[tid][3 * a2b[v0] + 1] += 2 * coeff_fea * n1n2 * dn1n2[0][1];
		g_omp[tid][3 * a2b[v0] + 2] += 2 * coeff_fea * n1n2 * dn1n2[0][2];
		g_omp[tid][3 * a2b[v1] + 0] += 2 * coeff_fea * n1n2 * dn1n2[1][0];
		g_omp[tid][3 * a2b[v1] + 1] += 2 * coeff_fea * n1n2 * dn1n2[1][1];
		g_omp[tid][3 * a2b[v1] + 2] += 2 * coeff_fea * n1n2 * dn1n2[1][2];
		g_omp[tid][3 * a2b[v2] + 0] += 2 * coeff_fea * n1n2 * dn1n2[2][0];
		g_omp[tid][3 * a2b[v2] + 1] += 2 * coeff_fea * n1n2 * dn1n2[2][1];
		g_omp[tid][3 * a2b[v2] + 2] += 2 * coeff_fea * n1n2 * dn1n2[2][2];
		g_omp[tid][3 * a2b[v3] + 0] += 2 * coeff_fea * n1n2 * dn1n2[3][0];
		g_omp[tid][3 * a2b[v3] + 1] += 2 * coeff_fea * n1n2 * dn1n2[3][1];
		g_omp[tid][3 * a2b[v3] + 2] += 2 * coeff_fea * n1n2 * dn1n2[3][2];
	}
	//feature normal constraint, normal is decided by filtering or polyline 
	//double coeff_fea_normal = lambda_array[0] / (1.0 * feature_edges.size() + LEN_LB);
	double coeff_fea_normal = 0.5 * lambda_array[2] / (1.0 * feature_edges.size());
	//assert(feature_edges.size() == target_feature_edge_normal.size());
	int n_fea_normal = feature_edges.size();
	
	if (lambda_array[2] < LEN_LB)
		n_fea_normal = 0;
	//n_fea_normal = 0;
	//feature normal guidance
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < n_fea_normal; p++)
	{
		int tid = omp_get_thread_num();
		int eid = feature_edges[p];
		const OpenVolumeMesh::Geometry::Vec3d target_n = target_fea_n[eid];
		int aid0 = feature_e2v[eid].first, aid1 = feature_e2v[eid].second;
		int bid0 = a2b[aid0], bid1 = a2b[aid1];
		double x1, x2, y1, y2, z1, z2;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		double n[3];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		double len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		double len = sqrt(len2);
		assert(_finite(len) == 1 && _finite(len2) == 1);
		if (len < LEN_LB)
		{
			len = LEN_LB;
			len2 = len * len;
		}
		double nn[3];
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		assert(_finite(nn[0]) && _finite(nn[1]) && _finite(nn[2]));
		double f_ori = 0;
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			func += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			f_ori += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		double dlen[2][3]; //first dimension: vert id, second dim: axis
		for (size_t i = 0; i < 3; i++)
		{
			dlen[0][i] = -n[i] / len;
			dlen[1][i] = n[i] / len;
		}
		double diff[3];
		for (size_t i = 0; i < 3; i++)
		{
			diff[i] = nn[i] - target_n[i];
		}
		double ddiff0[3][3], ddiff1[3][3]; //second dimension: axis
		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				ddiff0[i][j] = -n[i] * dlen[0][j] / len2;
				ddiff1[i][j] = -n[i] * dlen[1][j] / len2;
			}
			ddiff0[i][i] += -1.0 / len;
			ddiff1[i][i] += 1.0 / len;
		}
		double df[2][3];
		for (size_t i = 0; i < 3; i++)
		{
			df[0][i] = 0.0; //i the axis
			df[1][i] = 0.0;
			for (size_t j = 0; j < 3; j++)
			{
				df[0][i] += 2 * coeff_fea_normal * diff[j] * ddiff0[j][i];
				df[1][i] += 2 * coeff_fea_normal * diff[j] * ddiff1[j][i];
			}
		}
		/*g[3 * bid0 + 0] += df[0][0];
		g[3 * bid0 + 1] += df[0][1];
		g[3 * bid0 + 2] += df[0][2];
		g[3 * bid1 + 0] += df[1][0];
		g[3 * bid1 + 1] += df[1][1];
		g[3 * bid1 + 2] += df[1][2];*/
		g_omp[tid][3 * bid0 + 0] += df[0][0];
		g_omp[tid][3 * bid0 + 1] += df[0][1];
		g_omp[tid][3 * bid0 + 2] += df[0][2];
		g_omp[tid][3 * bid1 + 0] += df[1][0];
		g_omp[tid][3 * bid1 + 1] += df[1][1];
		g_omp[tid][3 * bid1 + 2] += df[1][2];
#if 0
		//for debugging
		std::cout << "v0 grad: " << std::endl;
		for (size_t i = 0; i < 3; i++)
		{
			std::cout << df[0][i] << std::endl;
		}
		std::cout << "v1 grad: " << std::endl;
		for (size_t i = 0; i < 3; i++)
		{
			std::cout << df[1][i] << std::endl;
		}
		//difference method
		//v0x
		double tmpf = 0;
		x1 = x[3 * bid0] + EPSL;
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 x: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][0] - (tmpf - f_ori) / EPSL) < BOUND);
		//v0y
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1] + EPSL;
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 y: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][1] - (tmpf - f_ori) / EPSL) < BOUND);
		//v0z
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2] + EPSL;
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 z: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][2] - (tmpf - f_ori) / EPSL) < BOUND);
		//v1x
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1] + EPSL;
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 x: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][0] - (tmpf - f_ori) / EPSL) < BOUND);
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1] + EPSL;
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 y: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][1] - (tmpf - f_ori) / EPSL) < BOUND);
		//v1z
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2] + EPSL;
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 z: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][2] - (tmpf - f_ori) / EPSL) < BOUND);
		tmpf = 0;
#endif
	}
	int n_fea_normal_nf = feature_edges.size();
	if (lambda_array[2] < LEN_LB)
		n_fea_normal_nf = 0;
	//n_fea_normal_nf = 0;
	//feature normal flattening
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < n_fea_normal_nf; p++)
	{
		int tid = omp_get_thread_num();
		int eid = feature_edges[p];
		const OpenVolumeMesh::Geometry::Vec3d target_n = target_fea_n[eid];
		int aid0 = feature_e2v[eid].first, aid1 = feature_e2v[eid].second;
		int bid0 = a2b[aid0], bid1 = a2b[aid1];
		double x1, x2, y1, y2, z1, z2;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		double n[3];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		double len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		double len = sqrt(len2);
		assert(_finite(len) == 1 && _finite(len2) == 1);
		if (len < LEN_LB)
		{
			len = LEN_LB;
			len2 = len * len;
		}
		double nn[3];
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		assert(_finite(nn[0]) && _finite(nn[1]) && _finite(nn[2]));
		//double f_ori = 0;
		//for (size_t i = 0; i < 3; i++)
		//{
		//	//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		//	//func += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		//	func += 0.5 * coeff_fea_normal * ()
		//	//f_ori += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		//}
		func += coeff_fea_normal * (nn[0] * nn[0] * nn[1] * nn[1] + nn[2] * nn[2] * nn[1] * nn[1] + nn[2] * nn[2] * nn[0] * nn[0]);
		double dlen[2][3]; //first dimension: vert id, second dim: axis
		for (size_t i = 0; i < 3; i++)
		{
			dlen[0][i] = -n[i] / len;
			dlen[1][i] = n[i] / len;
		}
		//double diff[3];
		double ddiff0[3][3], ddiff1[3][3]; //second dimension: axis
		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				ddiff0[i][j] = -n[i] * dlen[0][j] / len2;
				ddiff1[i][j] = -n[i] * dlen[1][j] / len2;
			}
			ddiff0[i][i] += -1.0 / len;
			ddiff1[i][i] += 1.0 / len;
		}
		double df[2][3];
		for (int i = 0; i < 3; i++)
		{
			df[0][i] = 0.0; //i the axis
			df[1][i] = 0.0;
			
			for (int j = 0; j < 3; j++)
			{
				int id1 = (j + 1) % 3;
				int id2 = (j + 2) % 3;
				df[0][i] += 2 * coeff_fea_normal * (nn[j] * nn[id1] * nn[id1] + nn[j] * nn[id2] * nn[id2]) * ddiff0[j][i];
				df[1][i] += 2 * coeff_fea_normal * (nn[j] * nn[id1] * nn[id1] + nn[j] * nn[id2] * nn[id2]) * ddiff1[j][i];
			}
		}
		/*g[3 * bid0 + 0] += df[0][0];
		g[3 * bid0 + 1] += df[0][1];
		g[3 * bid0 + 2] += df[0][2];
		g[3 * bid1 + 0] += df[1][0];
		g[3 * bid1 + 1] += df[1][1];
		g[3 * bid1 + 2] += df[1][2];*/
		g_omp[tid][3 * bid0 + 0] += df[0][0];
		g_omp[tid][3 * bid0 + 1] += df[0][1];
		g_omp[tid][3 * bid0 + 2] += df[0][2];
		g_omp[tid][3 * bid1 + 0] += df[1][0];
		g_omp[tid][3 * bid1 + 1] += df[1][1];
		g_omp[tid][3 * bid1 + 2] += df[1][2];
#if 0
		//for debugging
		std::cout << "v0 grad: " << std::endl;
		for (size_t i = 0; i < 3; i++)
		{
			std::cout << df[0][i] << std::endl;
		}
		std::cout << "v1 grad: " << std::endl;
		for (size_t i = 0; i < 3; i++)
		{
			std::cout << df[1][i] << std::endl;
		}
		//difference method
		//v0x
		double tmpf = 0;
		x1 = x[3 * bid0] + EPSL;
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 x: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][0] - (tmpf - f_ori) / EPSL) < BOUND);
		//v0y
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1] + EPSL;
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 y: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][1] - (tmpf - f_ori) / EPSL) < BOUND);
		//v0z
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2] + EPSL;
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 z: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][2] - (tmpf - f_ori) / EPSL) < BOUND);
		//v1x
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1] + EPSL;
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 x: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][0] - (tmpf - f_ori) / EPSL) < BOUND);
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1] + EPSL;
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 y: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][1] - (tmpf - f_ori) / EPSL) < BOUND);
		//v1z
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2] + EPSL;
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 z: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][2] - (tmpf - f_ori) / EPSL) < BOUND);
		tmpf = 0;
#endif
	}
	//smooth term
	double coeff_smooth = 0.5 * lambda_array[1] / boundary_verts.size();
	int bnv = boundary_verts.size();
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < bnv; p++)
	{
		//if number is 1, then skip
		int tid = omp_get_thread_num();
		int center = boundary_verts[p];
		
		
		std::vector<int> lidx;
		int nn = 0;
		if (feature_g2l[center] != -1)
		{
			nn = feature_v2v[center].size();
			for (size_t i = 0; i < nn; i++)
			{
				lidx.push_back(feature_v2v[center][i]);
			}
		}
		else
		{
			//non-feature, all boundary
			for (size_t i = 0; i < bvv_id[center].size(); i++)
			{
				lidx.push_back(bvv_id[center][i]);
			}
			nn = lidx.size();
		}
		
		
		if (nn <= 1) continue;
		double xc, yc, zc, xn, yn, zn;
		xc = x[3 * center];
		yc = x[3 * center + 1];
		zc = x[3 * center + 2];
		xn = 0.0; yn = 0.0; zn = 0.0;
		for (size_t i = 0; i < nn; i++)
		{
			xn += 1.0 / nn * x[3 * lidx[i]];
			yn += 1.0 / nn * x[3 * lidx[i] + 1];
			zn += 1.0 / nn * x[3 * lidx[i] + 2];
		}
		func += coeff_smooth * ((xc - xn) * (xc - xn) + (yc - yn) * (yc - yn) + (zc - zn) * (zc - zn));
		g_omp[tid][3 * center + 0] += 2 * coeff_smooth * (xc - xn);
		g_omp[tid][3 * center + 1] += 2 * coeff_smooth * (yc - yn);
		g_omp[tid][3 * center + 2] += 2 * coeff_smooth * (zc - zn);
		/*g_omp[tid][3 * bid1 + 0] += df[1][0];
		g_omp[tid][3 * bid1 + 1] += df[1][1];
		g_omp[tid][3 * bid1 + 2] += df[1][2];*/
		for (size_t i = 0; i < nn; i++)
		{
			g_omp[tid][3 * lidx[i] + 0] += -2 * coeff_smooth / nn * (xc - xn);
			g_omp[tid][3 * lidx[i] + 1] += -2 * coeff_smooth / nn * (yc - yn);
			g_omp[tid][3 * lidx[i] + 2] += -2 * coeff_smooth / nn * (zc - zn);
		}
	}
	
	//update gradient and function
	*f += func;
#pragma omp parallel for
	for (int i = 0; i < 3 * nv; i++)
	{
		for (size_t j = 0; j < num_core; j++)
		{
			g[i] += g_omp[j][i];
		}
	}
}
void evalfunc_de_LBFGS_polylines(int N, double* x, double *prev_x, double* f, double* g, void* user_pointer)
{
	deform_pointer *dp = (deform_pointer *)(user_pointer);
	const vector<double> &dpx = dp->pdi->dpx;
	const vector<double> &dpy = dp->pdi->dpy;
	const vector<double> &dpz = dp->pdi->dpz;
	const vector<int> &boundary_verts = dp->pdi->boundary_verts;
	const vector<OpenVolumeMesh::Geometry::Vec3i> &bfv_id = dp->pdi->bfv_id;
	const std::vector<std::vector<int>> &vertex_cell = dp->pdi->vertex_cell;
	const std::vector<std::vector<std::vector<int>>> &vertex_cell_vertex = dp->pdi->vertex_cell_vertex;
	const std::vector< OpenVolumeMesh::Geometry::Vec3d > &target_bfn = dp->pdi->target_bfn;
	const std::vector< std::vector<int> > &bvv_id = dp->pdi->bvv_id;
	const std::vector< std::vector<int> > &bvf_id = dp->pdi->bvf_id;
	const std::vector<int> &a2b = dp->pdi->a2b;
	//std::vector<std::vector<int>> &cell_vertex = dp->pdi->cell_vertex; //not using right order version, cause inverse is not in right order
	const std::vector< std::vector<std::vector<double> > > &vcv_S = dp->pdi->vcv_S;
	const std::vector< Eigen::Matrix3d > &cell_S = dp->pdi->cell_S;
	const std::vector<std::vector<int>> &cell_vertex = dp->pdi->cell_vertex;
	const std::vector<double> &lambda_array = dp->pdi->lambda_array_polyline;
	int boundary_face_number = dp->pdi->boundary_face_number;
	//int boundary_vertface_number = dp->pdi->boundary_vertface_number;
	int nc = dp->pdi->bfv_id.size();
	const vector<vector<int>> &bef_id = dp->pdi->bef_id;
	//feature edge part
	const std::vector<bool> &feature_edge_flag = dp->pdi->feature_edge_flag;
	const std::vector<std::vector<int>> &feature_v2e = dp->pdi->feature_v2e;
	const std::vector<std::vector<std::pair<int, int>>> &feature_neighbor_vert2cellpair = dp->pdi->feature_neighbor_vert2cellpair; //first face containging vert, second face the other direction
	//std::vector<std::vector<int>> feature_faces;
	const std::vector<OpenVolumeMesh::Geometry::Vec3d> &target_fea_n = dp->pdi->target_fea_n;
	//int n_feature_edges;
	const std::vector<int> &feature_edges = dp->pdi->feature_edge_array;
	const std::vector<double> &feature_edge_length = dp->pdi->feature_edge_length;
	const std::vector<std::pair<int, int>> &feature_e2v = dp->pdi->feature_e2v; //consistent feature edges
	const std::vector<std::vector<int>> &feature_neighbor = dp->pdi->feature_neighbor; //neighbor of edges at most 2 neighbors, cut by valence 3 vertex
	
	const std::vector<int> &feature_g2l = dp->pdi->feature_g2l;
	const std::vector<int> &feature_l2g = dp->pdi->feature_l2g;
	const std::vector<std::vector<int>> &feature_v2v = dp->pdi->feature_v2v;
	bool polyline_ns_flag = dp->pdi->polyline_ns_flag;
																					   //const std::vector<OpenVolumeMesh::Geometry::Vec3d> &target_feature_edge_normal = dp->pdi->target_feature_edge_normal;
	/*int v_id_0, v_id_1, v_id_2, v_id_3, k, c_id;
	double D00, D10, D20, D01, D11, D21, D02, D12, D22;
	double C00, C01, C02, C10, C11, C12, C20, C21, C22;
	double A00, A01, A02, A10, A11, A12, A20, A21, A22;
	double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
	double AF, AF_05, I_AF_05, AF2, AF_I, AF_I_05, I_AF_I_05, exp_k, exp_e;
	double D_AF_x, D_AF_y, D_AF_z, D_AF2_x, D_AF2_y, D_AF2_z, D_AF_I_x, D_AF_I_y, D_AF_I_z, D_det_A_x, D_det_A_y, D_det_A_z, inv_det_A_2_05;
	double d_A00_x, d_A10_y, d_A20_z, d_A01_x, d_A11_y, d_A21_z, d_A02_x, d_A12_y, d_A22_z;
	double dvex, dvey, dvez, tmp_g, dgx, dgy, dgz, dex, dey, dez, e, mips_e, volume_e, d_mips_e_x, d_mips_e_y, d_mips_e_z;
	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	double len;*/
	//double min_radius = 1e30;
	double alpha = 0.5; //double beta = 0.5;
	int bv_size = boundary_verts.size();
	int nv = dp->pdi->bvf_id.size();
	assert(bv_size > 0);
	//std::vector<std::vector<double>> g_omp;
	unsigned int num_core = omp_get_max_threads();
	std::vector<std::vector<double>> g_omp(num_core, std::vector<double>(3 * nv, 0.0));
	double func = 0.0;
	*f = 0.0;
	for (size_t i = 0; i < 3 * feature_l2g.size(); i++)
	{
		g[i] = 0.0;
	}
	
	//feature edge length preserving preserving
	//Sigma ||len - ori_len||^2
	
	double coeff_fea = 0.5 / (1.0 * feature_edges.size());
	int n_fea = feature_edges.size();
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < n_fea; p++)
	{
		int tid = omp_get_thread_num();
		int eid = feature_edges[p];
		double ori_len = feature_edge_length[p];
		const OpenVolumeMesh::Geometry::Vec3d target_n = target_fea_n[eid];
		int aid0 = feature_e2v[eid].first, aid1 = feature_e2v[eid].second;
		//int bid0 = a2b[aid0], bid1 = a2b[aid1];
		int bid0 = feature_g2l[aid0], bid1 = feature_g2l[aid1];
		double x1, x2, y1, y2, z1, z2;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		double n[3];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		double len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		double len = sqrt(len2);
		assert(_finite(len) == 1 && _finite(len2) == 1);
		if (len < LEN_LB)
		{
			len = LEN_LB;
			len2 = len * len;
		}
		double nn[3];
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		assert(_finite(nn[0]) && _finite(nn[1]) && _finite(nn[2]));
		double f_ori = 0;
		func += coeff_fea * (len - ori_len) * (len - ori_len);
		double dlen[2][3]; //first dimension: vert id, second dim: axis
		for (size_t i = 0; i < 3; i++)
		{
			dlen[0][i] = -n[i] / len;
			dlen[1][i] = n[i] / len;
		}
		double df[2][3];
		for (size_t i = 0; i < 3; i++)
		{
			df[0][i] = 2 * coeff_fea * (len - ori_len) * dlen[0][i];
			df[1][i] = 2 * coeff_fea * (len - ori_len) * dlen[1][i];
		}
		/*g[3 * bid0 + 0] += df[0][0];
		g[3 * bid0 + 1] += df[0][1];
		g[3 * bid0 + 2] += df[0][2];
		g[3 * bid1 + 0] += df[1][0];
		g[3 * bid1 + 1] += df[1][1];
		g[3 * bid1 + 2] += df[1][2];*/
		g_omp[tid][3 * bid0 + 0] += df[0][0];
		g_omp[tid][3 * bid0 + 1] += df[0][1];
		g_omp[tid][3 * bid0 + 2] += df[0][2];
		g_omp[tid][3 * bid1 + 0] += df[1][0];
		g_omp[tid][3 * bid1 + 1] += df[1][1];
		g_omp[tid][3 * bid1 + 2] += df[1][2];
	}
	
	//feature normal constraint
	//double coeff_fea_normal = lambda_array[0] / (1.0 * feature_edges.size() + LEN_LB);
	double coeff_fea_normal = 0.5 * lambda_array[0] / (1.0 * feature_edges.size());
	//assert(feature_edges.size() == target_feature_edge_normal.size());
	int n_fea_normal = feature_edges.size();
	if (!polyline_ns_flag)
		n_fea_normal = 0;
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < n_fea_normal; p++)
	{
		int tid = omp_get_thread_num();
		int eid = feature_edges[p];
		const OpenVolumeMesh::Geometry::Vec3d target_n = target_fea_n[eid];
		int aid0 = feature_e2v[eid].first, aid1 = feature_e2v[eid].second;
		int bid0 = feature_g2l[aid0], bid1 = feature_g2l[aid1];
		double x1, x2, y1, y2, z1, z2;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		double n[3];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		double len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		double len = sqrt(len2);
		assert(_finite(len) == 1 && _finite(len2) == 1);
		if (len < LEN_LB)
		{
			len = LEN_LB;
			len2 = len * len;
		}
		double nn[3];
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		assert(_finite(nn[0]) && _finite(nn[1]) && _finite(nn[2]));
		double f_ori = 0;
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			func += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			f_ori += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		double dlen[2][3]; //first dimension: vert id, second dim: axis
		for (size_t i = 0; i < 3; i++)
		{
			dlen[0][i] = -n[i] / len;
			dlen[1][i] = n[i] / len;
		}
		double diff[3];
		for (size_t i = 0; i < 3; i++)
		{
			diff[i] = nn[i] - target_n[i];
		}
		double ddiff0[3][3], ddiff1[3][3]; //second dimension: axis
		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				ddiff0[i][j] = -n[i] * dlen[0][j] / len2;
				ddiff1[i][j] = -n[i] * dlen[1][j] / len2;
			}
			ddiff0[i][i] += -1.0 / len;
			ddiff1[i][i] += 1.0 / len;
		}
		double df[2][3];
		for (size_t i = 0; i < 3; i++)
		{
			df[0][i] = 0.0; //i the axis
			df[1][i] = 0.0;
			for (size_t j = 0; j < 3; j++)
			{
				df[0][i] += 2 * coeff_fea_normal * diff[j] * ddiff0[j][i];
				df[1][i] += 2 * coeff_fea_normal * diff[j] * ddiff1[j][i];
			}
		}
		/*g[3 * bid0 + 0] += df[0][0];
		g[3 * bid0 + 1] += df[0][1];
		g[3 * bid0 + 2] += df[0][2];
		g[3 * bid1 + 0] += df[1][0];
		g[3 * bid1 + 1] += df[1][1];
		g[3 * bid1 + 2] += df[1][2];*/
		g_omp[tid][3 * bid0 + 0] += df[0][0];
		g_omp[tid][3 * bid0 + 1] += df[0][1];
		g_omp[tid][3 * bid0 + 2] += df[0][2];
		g_omp[tid][3 * bid1 + 0] += df[1][0];
		g_omp[tid][3 * bid1 + 1] += df[1][1];
		g_omp[tid][3 * bid1 + 2] += df[1][2];
#if 0
		//for debugging
		std::cout << "v0 grad: " << std::endl;
		for (size_t i = 0; i < 3; i++)
		{
			std::cout << df[0][i] << std::endl;
		}
		std::cout << "v1 grad: " << std::endl;
		for (size_t i = 0; i < 3; i++)
		{
			std::cout << df[1][i] << std::endl;
		}
		//difference method
		//v0x
		double tmpf = 0;
		x1 = x[3 * bid0] + EPSL;
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 x: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][0] - (tmpf - f_ori) / EPSL) < BOUND);
		//v0y
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1] + EPSL;
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 y: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][1] - (tmpf - f_ori) / EPSL) < BOUND);
		//v0z
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2] + EPSL;
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 z: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][2] - (tmpf - f_ori) / EPSL) < BOUND);
		//v1x
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1] + EPSL;
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 x: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][0] - (tmpf - f_ori) / EPSL) < BOUND);
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1] + EPSL;
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 y: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][1] - (tmpf - f_ori) / EPSL) < BOUND);
		//v1z
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2] + EPSL;
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 z: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][2] - (tmpf - f_ori) / EPSL) < BOUND);
		tmpf = 0;
#endif
	}
	int n_fea_normal_nf = feature_edges.size();
	if (polyline_ns_flag)
		n_fea_normal_nf = 0;
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < n_fea_normal_nf; p++)
	{
		int tid = omp_get_thread_num();
		int eid = feature_edges[p];
		const OpenVolumeMesh::Geometry::Vec3d target_n = target_fea_n[eid];
		int aid0 = feature_e2v[eid].first, aid1 = feature_e2v[eid].second;
		int bid0 = feature_g2l[aid0], bid1 = feature_g2l[aid1];
		double x1, x2, y1, y2, z1, z2;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		double n[3];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		double len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		double len = sqrt(len2);
		assert(_finite(len) == 1 && _finite(len2) == 1);
		if (len < LEN_LB)
		{
			len = LEN_LB;
			len2 = len * len;
		}
		double nn[3];
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		assert(_finite(nn[0]) && _finite(nn[1]) && _finite(nn[2]));
		//double f_ori = 0;
		//for (size_t i = 0; i < 3; i++)
		//{
		//	//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		//	//func += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		//	func += 0.5 * coeff_fea_normal * ()
		//	//f_ori += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		//}
		func += coeff_fea_normal * (nn[0] * nn[0] * nn[1] * nn[1] + nn[2] * nn[2] * nn[1] * nn[1] + nn[2] * nn[2] * nn[0] * nn[0]);
		double dlen[2][3]; //first dimension: vert id, second dim: axis
		for (size_t i = 0; i < 3; i++)
		{
			dlen[0][i] = -n[i] / len;
			dlen[1][i] = n[i] / len;
		}
		//double diff[3];
		double ddiff0[3][3], ddiff1[3][3]; //second dimension: axis
		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				ddiff0[i][j] = -n[i] * dlen[0][j] / len2;
				ddiff1[i][j] = -n[i] * dlen[1][j] / len2;
			}
			ddiff0[i][i] += -1.0 / len;
			ddiff1[i][i] += 1.0 / len;
		}
		double df[2][3];
		for (int i = 0; i < 3; i++)
		{
			df[0][i] = 0.0; //i the axis
			df[1][i] = 0.0;
			for (int j = 0; j < 3; j++)
			{
				int id1 = (j + 1) % 3;
				int id2 = (j + 2) % 3;
				df[0][i] += 2 * coeff_fea_normal * (nn[j] * nn[id1] * nn[id1] + nn[j] * nn[id2] * nn[id2]) * ddiff0[j][i];
				df[1][i] += 2 * coeff_fea_normal * (nn[j] * nn[id1] * nn[id1] + nn[j] * nn[id2] * nn[id2]) * ddiff1[j][i];
			}
		}
		/*g[3 * bid0 + 0] += df[0][0];
		g[3 * bid0 + 1] += df[0][1];
		g[3 * bid0 + 2] += df[0][2];
		g[3 * bid1 + 0] += df[1][0];
		g[3 * bid1 + 1] += df[1][1];
		g[3 * bid1 + 2] += df[1][2];*/
		g_omp[tid][3 * bid0 + 0] += df[0][0];
		g_omp[tid][3 * bid0 + 1] += df[0][1];
		g_omp[tid][3 * bid0 + 2] += df[0][2];
		g_omp[tid][3 * bid1 + 0] += df[1][0];
		g_omp[tid][3 * bid1 + 1] += df[1][1];
		g_omp[tid][3 * bid1 + 2] += df[1][2];
#if 0
		//for debugging
		std::cout << "v0 grad: " << std::endl;
		for (size_t i = 0; i < 3; i++)
		{
			std::cout << df[0][i] << std::endl;
		}
		std::cout << "v1 grad: " << std::endl;
		for (size_t i = 0; i < 3; i++)
		{
			std::cout << df[1][i] << std::endl;
		}
		//difference method
		//v0x
		double tmpf = 0;
		x1 = x[3 * bid0] + EPSL;
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 x: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][0] - (tmpf - f_ori) / EPSL) < BOUND);
		//v0y
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1] + EPSL;
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 y: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][1] - (tmpf - f_ori) / EPSL) < BOUND);
		//v0z
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2] + EPSL;
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v0 z: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[0][2] - (tmpf - f_ori) / EPSL) < BOUND);
		//v1x
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1] + EPSL;
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 x: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][0] - (tmpf - f_ori) / EPSL) < BOUND);
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1] + EPSL;
		z2 = x[3 * bid1 + 2];
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 y: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][1] - (tmpf - f_ori) / EPSL) < BOUND);
		//v1z
		tmpf = 0.0;
		x1 = x[3 * bid0];
		y1 = x[3 * bid0 + 1];
		z1 = x[3 * bid0 + 2];
		x2 = x[3 * bid1];
		y2 = x[3 * bid1 + 1];
		z2 = x[3 * bid1 + 2] + EPSL;
		n[0] = x2 - x1;
		n[1] = y2 - y1;
		n[2] = z2 - z1;
		len2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		len = sqrt(len2);
		nn[0] = n[0] / len;
		nn[1] = n[1] / len;
		nn[2] = n[2] / len;
		//*f += (nn[0] - target_n[0]) * (nn - target_n[0]) + (nyn - target_n[1]) * (nyn - target_n[1]) + (nzn - target_n[2]) * (nzn - target_n[2]);
		for (size_t i = 0; i < 3; i++)
		{
			//*f += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
			tmpf += coeff_fea_normal * (nn[i] - target_n[i]) * (nn[i] - target_n[i]);
		}
		std::cout << "diff f v1 z: " << (tmpf - f_ori) / EPSL << std::endl;
		assert(abs(df[1][2] - (tmpf - f_ori) / EPSL) < BOUND);
		tmpf = 0;
#endif
	}
	double coeff_smooth = 0.5 * lambda_array[1] / (1.0 * feature_l2g.size());
	
	int n_fea_v = feature_l2g.size();
	if (polyline_ns_flag)
		n_fea_v = 0;
#pragma omp parallel for reduction(+:func)
	for (int p = 0; p < n_fea_v; p++)
	{
		//if number is 1, then skip
		int tid = omp_get_thread_num();
		int center = p;
		int gc = feature_l2g[center];
		int nn = feature_v2v[gc].size();
		if (nn == 1) continue;
		std::vector<int> lidx;
		for (size_t i = 0; i < nn; i++)
		{
			lidx.push_back(feature_g2l[feature_v2v[gc][i]]);
		}
		double xc, yc, zc, xn, yn, zn;
		xc = x[3 * center];
		yc = x[3 * center + 1];
		zc = x[3 * center + 2];
		xn = 0.0; yn = 0.0; zn = 0.0;
		for (size_t i = 0; i < nn; i++)
		{
			xn += 1.0 / nn * x[3 * lidx[i]];
			yn += 1.0 / nn * x[3 * lidx[i] + 1];
			zn += 1.0 / nn * x[3 * lidx[i] + 2];
		}
		
		func += coeff_smooth * ((xc - xn) * (xc - xn) + (yc - yn) * (yc - yn) + (zc - zn) * (zc - zn));
		
		g_omp[tid][3 * center + 0] += 2 * coeff_smooth * (xc - xn);
		g_omp[tid][3 * center + 1] += 2 * coeff_smooth * (yc - yn);
		g_omp[tid][3 * center + 2] += 2 * coeff_smooth * (zc - zn);
		/*g_omp[tid][3 * bid1 + 0] += df[1][0];
		g_omp[tid][3 * bid1 + 1] += df[1][1];
		g_omp[tid][3 * bid1 + 2] += df[1][2];*/
		for (size_t i = 0; i < nn; i++)
		{
			g_omp[tid][3 * lidx[i] + 0] += -2 * coeff_smooth / nn * (xc - xn);
			g_omp[tid][3 * lidx[i] + 1] += -2 * coeff_smooth / nn * (yc - yn);
			g_omp[tid][3 * lidx[i] + 2] += -2 * coeff_smooth / nn * (zc - zn);
		}
		
		
	}
	//update gradient and function
	*f += func;
#pragma omp parallel for
	for (int i = 0; i < 3 * feature_l2g.size(); i++)
	{
		for (size_t j = 0; j < num_core; j++)
		{
			g[i] += g_omp[j][i];
		}
	}
}
void newiteration_de_LBFGS(int iter, int call_iter, double *x, double* f, double *g, double* gnorm, void* user_pointer)
{
	printf("%d %d; %4.3e, %4.3e\n", iter, call_iter, f[0], gnorm[0]);
}
bool de_LBFGS(deform_pointer *dp, std::vector<double>& X)
{
	//user_pointer_ad up;
	//up.bfn = bfn;
	//LBFGS
	double parameter[20];
	int info[20];
	//initialize
	INIT_HLBFGS(parameter, info);
	parameter[2] = 0.9; //gtol, default value
	info[4] = 500; //number of iteration
	info[6] = 0; //update interval of Hessian
	info[7] = 0; //if with hessian 1, without 0
	info[10] = 0;
	int N = dp->pdi->boundary_verts.size() * 3;
	int M = 7;
	//m is the number of history value
	//n is the number of variables
	printf("-------------------------------\n");
	printf("start LBFGS\n");
	HLBFGS(N, M, &X[0], evalfunc_de_LBFGS, 0, HLBFGS_UPDATE_Hessian, newiteration_de_LBFGS, parameter, info, dp);
	return true;
}
bool de_LBFGS_all(deform_pointer *dp, std::vector<double>& X)
{
	//user_pointer_ad up;
	//up.bfn = bfn;
	//LBFGS
	double parameter[20];
	int info[20];
	//initialize
	INIT_HLBFGS(parameter, info);
	parameter[2] = 0.9; //gtol, default value
	info[4] = 500; //number of iteration
	info[6] = 0; //update interval of Hessian
	info[7] = 0; //if with hessian 1, without 0
	info[10] = 0;
	int N = dp->pdi->bvf_id.size() * 3;
	int M = 7;
	//m is the number of history value
	//n is the number of variables
	printf("-------------------------------\n");
	printf("start LBFGS\n");
	HLBFGS(N, M, &X[0], evalfunc_de_LBFGS_all, 0, HLBFGS_UPDATE_Hessian, newiteration_de_LBFGS, parameter, info, dp);
	return true;
}
bool de_LBFGS_all_nf(deform_pointer *dp, std::vector<double>& X)
{
	//user_pointer_ad up;
	//up.bfn = bfn;
	//LBFGS
	double parameter[20];
	int info[20];
	//initialize
	INIT_HLBFGS(parameter, info);
	parameter[2] = 0.9; //gtol, default value
	info[4] = 500; //number of iteration
	info[6] = 0; //update interval of Hessian
	info[7] = 0; //if with hessian 1, without 0
	info[10] = 0;
	int N = dp->pdi->bvf_id.size() * 3;
	int M = 7;
	//m is the number of history value
	//n is the number of variables
	printf("-------------------------------\n");
	printf("start LBFGS for normal flattening\n");
	HLBFGS(N, M, &X[0], evalfunc_de_LBFGS_all_nf, 0, HLBFGS_UPDATE_Hessian, newiteration_de_LBFGS, parameter, info, dp);
	return true;
}
bool de_LBFGS_polylines(deform_pointer *dp, std::vector<double>& X)
{
	//user_pointer_ad up;
	//up.bfn = bfn;
	//LBFGS
	double parameter[20];
	int info[20];
	//initialize
	INIT_HLBFGS(parameter, info);
	parameter[2] = 0.9; //gtol, default value
	info[4] = 500; //number of iteration
	info[6] = 0; //update interval of Hessian
	info[7] = 0; //if with hessian 1, without 0
	info[10] = 0;
	int N = dp->pdi->feature_l2g.size() * 3;
	int M = 7;
	//m is the number of history value
	//n is the number of variables
	printf("-------------------------------\n");
	printf("start LBFGS for normal flattening\n");
	HLBFGS(N, M, &X[0], evalfunc_de_LBFGS_polylines, 0, HLBFGS_UPDATE_Hessian, newiteration_de_LBFGS, parameter, info, dp);
	return true;
}