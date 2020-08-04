#include "Polycube_Boundary_Map_IVF.h"

polycube_boundary_mapping_IVF::polycube_boundary_mapping_IVF()
{
	reset();
}

polycube_boundary_mapping_IVF::~polycube_boundary_mapping_IVF()
{

}

void polycube_boundary_mapping_IVF::reset()
{
	lambda = 1e-1; min_lambda = 1e-5; max_lambda = 1e16;
	ivf_alpha = 1e5; min_ivf_alpha = 1e3; max_ivf_alpha = 1e16;
	step_size = 1.0; min_step_size = 1.0e-16; max_step_size = 1.0;
	amips_s = 2.0; var_num = 0;
}

void polycube_boundary_mapping_IVF::update_ivf_alpha()
{
	double s_ivf_alpha = ivf_alpha;
	ivf_alpha = max_ivf_alpha;

	ivf_alpha = current_dis_energy*min_ivf_alpha / current_ivf_energy;
	if (ivf_alpha > max_ivf_alpha) ivf_alpha = max_ivf_alpha;
	if (ivf_alpha < min_ivf_alpha) ivf_alpha = min_ivf_alpha;
	//if (ivf_alpha < s_ivf_alpha) ivf_alpha = s_ivf_alpha;
}

bool polycube_boundary_mapping_IVF::initilize_affine_trans(VolumeMesh* mesh_, int chart_label,
	const std::vector<OpenVolumeMesh::Geometry::Vec3d>& src_pos,
	const std::vector<OpenVolumeMesh::Geometry::Vec3i>& bfv_id,
	const std::vector< std::vector<int> >& bvf_id,
	const std::vector<int>& map_v, const int& var_v_count, const std::vector<int>& map_c, const int& var_c_count)
{
	var_num = 6 * var_c_count; int flipped_count = 0;
	std::vector<Eigen::Triplet<double> > f_triplets; solution.resize(var_num);
	std::vector<Eigen::Triplet<double> > c_triplets; std::vector<Eigen::Triplet<double> > triplets;
	int constraint_count = 0; std::vector<double> b;
	int nc = bfv_id.size(); fe1.resize(nc); fe2.resize(nc); bp.resize(nc);
	double q0_x, q0_y, q1_x, q1_y, q2_x, q2_y;
	for (int bf_id = 0; bf_id < nc; ++bf_id)
	{
		int var_c_id = map_c[bf_id];
		if (var_c_id < 0) continue;

		const OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[bf_id];
		const OpenVolumeMesh::Geometry::Vec3d& p0 = src_pos[ one_bfv_id[0] ];
		const OpenVolumeMesh::Geometry::Vec3d& p1 = src_pos[ one_bfv_id[1] ];
		const OpenVolumeMesh::Geometry::Vec3d& p2 = src_pos[ one_bfv_id[2] ];

		OpenVolumeMesh::Geometry::Vec3d e1 = (p1 - p0); double l01 = e1.norm(); e1.normalize();
		OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0); n.normalize();
		OpenVolumeMesh::Geometry::Vec3d e2 = OpenVolumeMesh::Geometry::cross(n, e1);
		bp[bf_id] = p0; fe1[bf_id] = e1; fe2[bf_id] = e2;

		double x10 = l01; double x20 = OpenVolumeMesh::Geometry::dot(p2 - p0, e1);
		double y10 = 0.0; double y20 = OpenVolumeMesh::Geometry::dot(p2 - p0, e2);

		const OpenVolumeMesh::Geometry::Vec3d& q0 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[0]));
		const OpenVolumeMesh::Geometry::Vec3d& q1 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[1]));
		const OpenVolumeMesh::Geometry::Vec3d& q2 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[2]));

		// chart_label: 0 for x, 1 for -x, 2 for y, 3 for -y, 4 for z, 5 for -z
		if (chart_label == 0)//y, z
		{
			q0_x = q0[1]; q0_y = q0[2]; 
			q1_x = q1[1]; q1_y = q1[2]; 
			q2_x = q2[1]; q2_y = q2[2];
		}
		else if (chart_label == 1) // y, z
		{
			q0_x = -q0[1]; q0_y = q0[2];
			q1_x = -q1[1]; q1_y = q1[2];
			q2_x = -q2[1]; q2_y = q2[2];
		}
		else if (chart_label == 2)// z, x
		{
			q0_x = q0[2]; q0_y = q0[0];
			q1_x = q1[2]; q1_y = q1[0];
			q2_x = q2[2]; q2_y = q2[0];
		}
		else if (chart_label == 3) // z, x
		{
			q0_x = -q0[2]; q0_y = q0[0];
			q1_x = -q1[2]; q1_y = q1[0];
			q2_x = -q2[2]; q2_y = q2[0];
		}
		else if (chart_label == 4)// x,y
		{
			q0_x = q0[0]; q0_y = q0[1];
			q1_x = q1[0]; q1_y = q1[1];
			q2_x = q2[0]; q2_y = q2[1];
		}
		else if (chart_label == 5) // x,y
		{
			q0_x = -q0[0]; q0_y = q0[1];
			q1_x = -q1[0]; q1_y = q1[1];
			q2_x = -q2[0]; q2_y = q2[1];
		}
		//affine matrix
		double det_ = 1.0 / (x10*y20 - x20*y10);
		double A_00 = ((q1_x - q0_x)*y20 - (q2_x - q0_x)*y10)*det_; double A_01 = (-(q1_x - q0_x)*x20 + (q2_x - q0_x)*x10)*det_;
		double A_10 = ((q1_y - q0_y)*y20 - (q2_y - q0_y)*y10)*det_; double A_11 = (-(q1_y - q0_y)*x20 + (q2_y - q0_y)*x10)*det_;
		Eigen::Matrix2d A; A << A_00, A_01, A_10, A_11;
		Eigen::JacobiSVD<Eigen::Matrix2d> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::Matrix2d U = svd.matrixU(); Eigen::Matrix2d V = svd.matrixV();
		//Eigen::Matrix2d rot = U*V.transpose();
		Eigen::Matrix2d rot = V*U.transpose();
		Eigen::Vector2d sv = svd.singularValues();
		double det = rot.determinant();
#if 1
		if (det < 0)
		{
			//A_00 = 1.0; A_01 = 0.0; A_10 = 0.0; A_11 = 1.0;
			if (sv(0) < sv(1))
			{
				U(0, 0) = -U(0, 0); U(1, 0) = -U(1, 0);
			}
			else
			{
				U(0, 1) = -U(0, 1); U(1, 1) = -U(1, 1);
			}
			//rot = U*V.transpose();
			rot = V*U.transpose();
			A_00 = rot(0, 0); A_01 = rot(0, 1); A_10 = rot(1, 0); A_11 = rot(1, 1);

			++flipped_count;
		}
		else
		{
			double err = (sv(0) / sv(1) + sv(1) / sv(0) + 1 / (sv(0)*sv(1)) + sv(0)*sv(1))*0.25;
			if (err > 10)
			{
				A_00 = rot(0, 0); A_01 = rot(0, 1); A_10 = rot(1, 0); A_11 = rot(1, 1);
			}
		}
#else
		double K = 10;
		if (det < 0)
		{
			if (sv(0) < sv(1))
			{
				V(0, 0) = -V(0, 0); V(1, 0) = -V(1, 0);
				sv(0) = -sv(0);
			}
			else
			{
				V(0, 1) = -V(0, 1); V(1, 1) = -V(1, 1);
				sv(1) = -sv(1);
			}
			
			Eigen::Matrix2d sv_m; sv_m.setZero(); sv_m(1, 1) = 1 / K;
			if (sv(0) > K) { sv_m(0, 0) = K; }
			else if (sv(0) < 1 / K) { sv_m(0, 0) = 1 / K; }
			else sv_m(0, 0) = sv(0);
			A = U*sv_m*V.transpose();
				//A_00 = rot(0, 0); A_01 = rot(0, 1); A_10 = rot(1, 0); A_11 = rot(1, 1);
			++flipped_count;
		}
		else
		{
			double err = sv(0) > 1.0 / sv(1) ? sv(0) : 1.0 / sv(1);
			
			if (err > K)
			{
				Eigen::Matrix2d sv_m; sv_m.setZero();
				if (sv(1) > K) { sv_m(0, 0) = K; sv_m(1, 1) = K; }
				else if (sv(0) < 1 / K) { sv_m(0, 0) = 1 / K; sv_m(1, 1) = 1 / K; }
				else if (sv(0) > K)
				{
					sv_m(0, 0) = K;
					if (sv(1) >= 1 / K && sv(1) <= K)
					{
						sv_m(1, 1) = sv(1);
					}
					else if (sv(1) < 1 / K)
					{
						sv_m(1, 1) = 1 / K;
					}
				}
				else if (sv(0) < K && sv(0) > 1 / K)
				{
					sv_m(0, 0) = sv(0);
					if (sv(1) >= 1 / K && sv(1) <= K)
					{
						sv_m(1, 1) = sv(1);
					}
					else if (sv(1) < 1 / K)
					{
						sv_m(1, 1) = 1 / K;
					}
				}
				A = U*sv_m*V.transpose();
			}
		}
		A_00 = A(0, 0); A_01 = A(0, 1); A_10 = A(1, 0); A_11 = A(1, 1);
#endif

		solution(6 * var_c_id + 0) = A_00; solution(6 * var_c_id + 1) = A_01;
		solution(6 * var_c_id + 2) = A_10; solution(6 * var_c_id + 3) = A_11;

		OpenVolumeMesh::Geometry::Vec2d np0(0.0, 0.0);
		OpenVolumeMesh::Geometry::Vec2d np1(A_00*x10, A_10*x10);
		OpenVolumeMesh::Geometry::Vec2d np2(A_00*x20 + A_01*y20, A_10*x20 + A_11*y20);

		//translation
		solution(6 * var_c_id + 4) = ((q0_x + q1_x + q2_x) - (np0[0] + np1[0] + np2[0])) / 3.0;
		solution(6 * var_c_id + 5) = ((q0_y + q1_y + q2_y) - (np0[1] + np1[1] + np2[1])) / 3.0;

		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				f_triplets.push_back(Eigen::Triplet<double>(6 * var_c_id + i, 6 * var_c_id + j, 1.0));
			}
		}

		int var_v_id = map_v[ one_bfv_id[0] ]; //0
		if (var_v_id < 0) //fixed point
		{
			//x
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 0, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 1, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 4, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 0, 0.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 1, 0.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 4, 1.0));
			b.push_back(q0_x);
			//y
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 2, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 3, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 5, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 2, 0.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 3, 0.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 5, 1.0));
			b.push_back(q0_y);

			constraint_count += 2;
		}
		var_v_id = map_v[one_bfv_id[1]]; //1
		if (var_v_id < 0) //fixed point
		{
			//x
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 0, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 1, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 4, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 0, x10));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 1, y10));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 4, 1.0));
			b.push_back(q1_x);
			//y
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 2, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 3, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 5, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 2, x10));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 3, y10));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 5, 1.0));
			b.push_back(q1_y);

			constraint_count += 2;
		}
		var_v_id = map_v[one_bfv_id[2]]; //2
		if (var_v_id < 0) //fixed point
		{
			//x
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 0, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 1, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 4, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 0, x20));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 1, y20));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 4, 1.0));
			b.push_back(q2_x);
			//y
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 2, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 3, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 5, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 2, x20));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 3, y20));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 5, 1.0));
			b.push_back(q2_y);

			constraint_count += 2;
		}
	}
	
	int nv = bvf_id.size();
	for (int bv_id = 0; bv_id < nv; ++bv_id)
	{
		int var_v_id = map_v[bv_id];
		if (var_v_id < 0) continue;
		
		const std::vector<int>& one_bvf_id = bvf_id[bv_id];
		const OpenVolumeMesh::Geometry::Vec3d& p = src_pos[bv_id];
		for (int i = 0; i < one_bvf_id.size() - 1; ++i)
		{
			int f0 = one_bvf_id[i]; int f1 = one_bvf_id[i + 1];
			int var_c_id0 = map_c[f0]; int var_c_id1 = map_c[f1];
			double x0 = OpenVolumeMesh::Geometry::dot(p - bp[f0], fe1[f0]);
			double y0 = OpenVolumeMesh::Geometry::dot(p - bp[f0], fe2[f0]);
			double x1 = OpenVolumeMesh::Geometry::dot(p - bp[f1], fe1[f1]);
			double y1 = OpenVolumeMesh::Geometry::dot(p - bp[f1], fe2[f1]);
			//printf("%f %f %f %f\n",x0,y0,x1,y1);

			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id0 + 0, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id0 + 1, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id0 + 4, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id1 + 0, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id1 + 1, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id1 + 4, 1.0));

			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id0 + 2, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id0 + 3, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id0 + 5, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id1 + 2, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id1 + 3, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id1 + 5, 1.0));

			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id0 + 0, x0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id0 + 1, y0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id0 + 4, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id1 + 0, -x1));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id1 + 1, -y1));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id1 + 4, -1.0));

			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id0 + 2, x0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id0 + 3, y0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id0 + 5, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id1 + 2, -x1));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id1 + 3, -y1));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id1 + 5, -1.0));

			b.push_back(0.0); b.push_back(0.0);
			constraint_count += 2;
		}
	}

	Eigen::SparseMatrix<double> C; C.resize(constraint_count, var_num);
	C.setFromTriplets(c_triplets.begin(), c_triplets.end());
	CT = C.transpose();
	CTC.resize(var_num, var_num);
	CTC = CT*C;
	Eigen::SparseMatrix<double> H;
	H.resize(var_num, var_num);
	H.setFromTriplets(f_triplets.begin(), f_triplets.end());

	hessian.resize(var_num, var_num);
	hessian = H + CTC;
	lltSolver.analyzePattern(hessian);

	C *= 0.0;
	C.setFromTriplets(triplets.begin(), triplets.end());

	CT = C.transpose();
	CTC = CT*C;
	CB.resize(constraint_count);
	for (int i = 0; i < constraint_count; ++i)
	{
		CB(i) = b[i];
	}
	CBTCB = CB.transpose() * CB;

	gradient.resize(var_num); prevSolution.resize(var_num); step.resize(var_num);

	return (flipped_count == 0);
}

bool polycube_boundary_mapping_IVF::initilize_affine_trans(TetStructure<double> *tet_mesh_, int chart_label,
	const std::vector<OpenVolumeMesh::Geometry::Vec3d>& src_pos,
	const std::vector<OpenVolumeMesh::Geometry::Vec3i>& bfv_id,
	const std::vector< std::vector<int> >& bvf_id,
	const std::vector<int>& map_v, const int& var_v_count, const std::vector<int>& map_c, const int& var_c_count)
{

	const std::vector<Tetrahedron<double>*> &tetras = tet_mesh_->tetras;
	std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;

	var_num = 6 * var_c_count; int flipped_count = 0;
	std::vector<Eigen::Triplet<double> > f_triplets; solution.resize(var_num);
	std::vector<Eigen::Triplet<double> > c_triplets; std::vector<Eigen::Triplet<double> > triplets;
	int constraint_count = 0; std::vector<double> b;
	int nc = bfv_id.size(); fe1.resize(nc); fe2.resize(nc); bp.resize(nc);
	double q0_x, q0_y, q1_x, q1_y, q2_x, q2_y;
	for (int bf_id = 0; bf_id < nc; ++bf_id)
	{
		int var_c_id = map_c[bf_id];
		if (var_c_id < 0) continue;

		const OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[bf_id];
		const OpenVolumeMesh::Geometry::Vec3d& p0 = src_pos[one_bfv_id[0]];
		const OpenVolumeMesh::Geometry::Vec3d& p1 = src_pos[one_bfv_id[1]];
		const OpenVolumeMesh::Geometry::Vec3d& p2 = src_pos[one_bfv_id[2]];

		OpenVolumeMesh::Geometry::Vec3d e1 = (p1 - p0); double l01 = e1.norm(); e1.normalize();
		OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0); n.normalize();
		OpenVolumeMesh::Geometry::Vec3d e2 = OpenVolumeMesh::Geometry::cross(n, e1);
		bp[bf_id] = p0; fe1[bf_id] = e1; fe2[bf_id] = e2;

		double x10 = l01; double x20 = OpenVolumeMesh::Geometry::dot(p2 - p0, e1);
		double y10 = 0.0; double y20 = OpenVolumeMesh::Geometry::dot(p2 - p0, e2);

		/*const OpenVolumeMesh::Geometry::Vec3d& q0 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[0]));
		const OpenVolumeMesh::Geometry::Vec3d& q1 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[1]));
		const OpenVolumeMesh::Geometry::Vec3d& q2 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[2]));*/
		
		CVec<double, 3> q0 = tetra_vertices[one_bfv_id[0]]->pos;
		CVec<double, 3> q1 = tetra_vertices[one_bfv_id[1]]->pos;
		CVec<double, 3> q2 = tetra_vertices[one_bfv_id[2]]->pos;
		

		// chart_label: 0 for x, 1 for -x, 2 for y, 3 for -y, 4 for z, 5 for -z
		if (chart_label == 0)//y, z
		{
			q0_x = q0[1]; q0_y = q0[2];
			q1_x = q1[1]; q1_y = q1[2];
			q2_x = q2[1]; q2_y = q2[2];
		}
		else if (chart_label == 1) // y, z
		{
			q0_x = -q0[1]; q0_y = q0[2];
			q1_x = -q1[1]; q1_y = q1[2];
			q2_x = -q2[1]; q2_y = q2[2];
		}
		else if (chart_label == 2)// z, x
		{
			q0_x = q0[2]; q0_y = q0[0];
			q1_x = q1[2]; q1_y = q1[0];
			q2_x = q2[2]; q2_y = q2[0];
		}
		else if (chart_label == 3) // z, x
		{
			q0_x = -q0[2]; q0_y = q0[0];
			q1_x = -q1[2]; q1_y = q1[0];
			q2_x = -q2[2]; q2_y = q2[0];
		}
		else if (chart_label == 4)// x,y
		{
			q0_x = q0[0]; q0_y = q0[1];
			q1_x = q1[0]; q1_y = q1[1];
			q2_x = q2[0]; q2_y = q2[1];
		}
		else if (chart_label == 5) // x,y
		{
			q0_x = -q0[0]; q0_y = q0[1];
			q1_x = -q1[0]; q1_y = q1[1];
			q2_x = -q2[0]; q2_y = q2[1];
		}
		//affine matrix
		double det_ = 1.0 / (x10*y20 - x20 * y10);
		double A_00 = ((q1_x - q0_x)*y20 - (q2_x - q0_x)*y10)*det_; double A_01 = (-(q1_x - q0_x)*x20 + (q2_x - q0_x)*x10)*det_;
		double A_10 = ((q1_y - q0_y)*y20 - (q2_y - q0_y)*y10)*det_; double A_11 = (-(q1_y - q0_y)*x20 + (q2_y - q0_y)*x10)*det_;
		Eigen::Matrix2d A; A << A_00, A_01, A_10, A_11;
		Eigen::JacobiSVD<Eigen::Matrix2d> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::Matrix2d U = svd.matrixU(); Eigen::Matrix2d V = svd.matrixV();
		//Eigen::Matrix2d rot = U*V.transpose();
		Eigen::Matrix2d rot = V * U.transpose();
		Eigen::Vector2d sv = svd.singularValues();
		double det = rot.determinant();
#if 1
		if (det < 0)
		{
			//A_00 = 1.0; A_01 = 0.0; A_10 = 0.0; A_11 = 1.0;
			if (sv(0) < sv(1))
			{
				U(0, 0) = -U(0, 0); U(1, 0) = -U(1, 0);
			}
			else
			{
				U(0, 1) = -U(0, 1); U(1, 1) = -U(1, 1);
			}
			//rot = U*V.transpose();
			rot = V * U.transpose();
			A_00 = rot(0, 0); A_01 = rot(0, 1); A_10 = rot(1, 0); A_11 = rot(1, 1);

			++flipped_count;
		}
		else
		{
			double err = (sv(0) / sv(1) + sv(1) / sv(0) + 1 / (sv(0)*sv(1)) + sv(0)*sv(1))*0.25;
			if (err > 10)
			{
				A_00 = rot(0, 0); A_01 = rot(0, 1); A_10 = rot(1, 0); A_11 = rot(1, 1);
			}
		}
#else
		double K = 10;
		if (det < 0)
		{
			if (sv(0) < sv(1))
			{
				V(0, 0) = -V(0, 0); V(1, 0) = -V(1, 0);
				sv(0) = -sv(0);
			}
			else
			{
				V(0, 1) = -V(0, 1); V(1, 1) = -V(1, 1);
				sv(1) = -sv(1);
			}

			Eigen::Matrix2d sv_m; sv_m.setZero(); sv_m(1, 1) = 1 / K;
			if (sv(0) > K) { sv_m(0, 0) = K; }
			else if (sv(0) < 1 / K) { sv_m(0, 0) = 1 / K; }
			else sv_m(0, 0) = sv(0);
			A = U * sv_m*V.transpose();
			//A_00 = rot(0, 0); A_01 = rot(0, 1); A_10 = rot(1, 0); A_11 = rot(1, 1);
			++flipped_count;
		}
		else
		{
			double err = sv(0) > 1.0 / sv(1) ? sv(0) : 1.0 / sv(1);

			if (err > K)
			{
				Eigen::Matrix2d sv_m; sv_m.setZero();
				if (sv(1) > K) { sv_m(0, 0) = K; sv_m(1, 1) = K; }
				else if (sv(0) < 1 / K) { sv_m(0, 0) = 1 / K; sv_m(1, 1) = 1 / K; }
				else if (sv(0) > K)
				{
					sv_m(0, 0) = K;
					if (sv(1) >= 1 / K && sv(1) <= K)
					{
						sv_m(1, 1) = sv(1);
					}
					else if (sv(1) < 1 / K)
					{
						sv_m(1, 1) = 1 / K;
					}
				}
				else if (sv(0) < K && sv(0) > 1 / K)
				{
					sv_m(0, 0) = sv(0);
					if (sv(1) >= 1 / K && sv(1) <= K)
					{
						sv_m(1, 1) = sv(1);
					}
					else if (sv(1) < 1 / K)
					{
						sv_m(1, 1) = 1 / K;
					}
				}
				A = U * sv_m*V.transpose();
			}
		}
		A_00 = A(0, 0); A_01 = A(0, 1); A_10 = A(1, 0); A_11 = A(1, 1);
#endif

		solution(6 * var_c_id + 0) = A_00; solution(6 * var_c_id + 1) = A_01;
		solution(6 * var_c_id + 2) = A_10; solution(6 * var_c_id + 3) = A_11;

		OpenVolumeMesh::Geometry::Vec2d np0(0.0, 0.0);
		OpenVolumeMesh::Geometry::Vec2d np1(A_00*x10, A_10*x10);
		OpenVolumeMesh::Geometry::Vec2d np2(A_00*x20 + A_01 * y20, A_10*x20 + A_11 * y20);

		//translation
		solution(6 * var_c_id + 4) = ((q0_x + q1_x + q2_x) - (np0[0] + np1[0] + np2[0])) / 3.0;
		solution(6 * var_c_id + 5) = ((q0_y + q1_y + q2_y) - (np0[1] + np1[1] + np2[1])) / 3.0;

		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				f_triplets.push_back(Eigen::Triplet<double>(6 * var_c_id + i, 6 * var_c_id + j, 1.0));
			}
		}

		int var_v_id = map_v[one_bfv_id[0]]; //0
		if (var_v_id < 0) //fixed point
		{
			//x
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 0, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 1, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 4, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 0, 0.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 1, 0.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 4, 1.0));
			b.push_back(q0_x);
			//y
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 2, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 3, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 5, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 2, 0.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 3, 0.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 5, 1.0));
			b.push_back(q0_y);

			constraint_count += 2;
		}
		var_v_id = map_v[one_bfv_id[1]]; //1
		if (var_v_id < 0) //fixed point
		{
			//x
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 0, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 1, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 4, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 0, x10));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 1, y10));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 4, 1.0));
			b.push_back(q1_x);
			//y
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 2, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 3, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 5, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 2, x10));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 3, y10));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 5, 1.0));
			b.push_back(q1_y);

			constraint_count += 2;
		}
		var_v_id = map_v[one_bfv_id[2]]; //2
		if (var_v_id < 0) //fixed point
		{
			//x
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 0, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 1, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 4, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 0, x20));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 1, y20));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id + 4, 1.0));
			b.push_back(q2_x);
			//y
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 2, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 3, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 5, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 2, x20));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 3, y20));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id + 5, 1.0));
			b.push_back(q2_y);

			constraint_count += 2;
		}
			}

	int nv = bvf_id.size();
	for (int bv_id = 0; bv_id < nv; ++bv_id)
	{
		int var_v_id = map_v[bv_id];
		if (var_v_id < 0) continue;

		const std::vector<int>& one_bvf_id = bvf_id[bv_id];
		const OpenVolumeMesh::Geometry::Vec3d& p = src_pos[bv_id];
		for (int i = 0; i < one_bvf_id.size() - 1; ++i)
		{
			int f0 = one_bvf_id[i]; int f1 = one_bvf_id[i + 1];
			int var_c_id0 = map_c[f0]; int var_c_id1 = map_c[f1];
			double x0 = OpenVolumeMesh::Geometry::dot(p - bp[f0], fe1[f0]);
			double y0 = OpenVolumeMesh::Geometry::dot(p - bp[f0], fe2[f0]);
			double x1 = OpenVolumeMesh::Geometry::dot(p - bp[f1], fe1[f1]);
			double y1 = OpenVolumeMesh::Geometry::dot(p - bp[f1], fe2[f1]);
			//printf("%f %f %f %f\n",x0,y0,x1,y1);

			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id0 + 0, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id0 + 1, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id0 + 4, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id1 + 0, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id1 + 1, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id1 + 4, 1.0));

			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id0 + 2, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id0 + 3, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id0 + 5, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id1 + 2, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id1 + 3, 1.0));
			c_triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id1 + 5, 1.0));

			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id0 + 0, x0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id0 + 1, y0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id0 + 4, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id1 + 0, -x1));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id1 + 1, -y1));
			triplets.push_back(Eigen::Triplet<double>(constraint_count, 6 * var_c_id1 + 4, -1.0));

			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id0 + 2, x0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id0 + 3, y0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id0 + 5, 1.0));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id1 + 2, -x1));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id1 + 3, -y1));
			triplets.push_back(Eigen::Triplet<double>(constraint_count + 1, 6 * var_c_id1 + 5, -1.0));

			b.push_back(0.0); b.push_back(0.0);
			constraint_count += 2;
		}
	}

	Eigen::SparseMatrix<double> C; C.resize(constraint_count, var_num);
	C.setFromTriplets(c_triplets.begin(), c_triplets.end());
	CT = C.transpose();
	CTC.resize(var_num, var_num);
	CTC = CT * C;
	Eigen::SparseMatrix<double> H;
	H.resize(var_num, var_num);
	H.setFromTriplets(f_triplets.begin(), f_triplets.end());

	hessian.resize(var_num, var_num);
	hessian = H + CTC;
	lltSolver.analyzePattern(hessian);

	C *= 0.0;
	C.setFromTriplets(triplets.begin(), triplets.end());

	CT = C.transpose();
	CTC = CT * C;
	CB.resize(constraint_count);
	for (int i = 0; i < constraint_count; ++i)
	{
		CB(i) = b[i];
	}
	CBTCB = CB.transpose() * CB;

	gradient.resize(var_num); prevSolution.resize(var_num); step.resize(var_num);

	return (flipped_count == 0);
		
}


//function value, gradient, hessian
double polycube_boundary_mapping_IVF::compute_energy_derivative_AT(const std::vector<int>& map_c, const Eigen::Matrix<double, Eigen::Dynamic, 1>& x)
{
	current_ivf_energy = x.transpose()*CTC*x;
	current_ivf_energy -= 2.0*x.transpose()*CT*CB;
	current_ivf_energy = std::abs(current_ivf_energy + CBTCB);
	
	gradient = 2.0*ivf_alpha*(CTC*x - CT*CB);
	
	hessian *= 0; current_dis_energy = 0.0; double s = amips_s;
	double distortion_e;
	double g_a1, g_a2, g_b1, g_b2;
	double h_a1_a1, h_a1_a2, h_a1_b1, h_a1_b2, h_a2_a2, h_a2_b1, h_a2_b2, h_b1_b1, h_b1_b2, h_b2_b2;
	int nc = map_c.size();
	for (int face_id = 0; face_id < nc; ++face_id)
	{
		int var_c_id = map_c[face_id];
		if (var_c_id < 0) continue;
		double a1 = x(6 * var_c_id + 0); double a2 = x(6 * var_c_id + 1);
		double b1 = x(6 * var_c_id + 2); double b2 = x(6 * var_c_id + 3);
		double det_J = a1*b2 - a2*b1;
		double mips_e = 0.5*(a1*a1 + a2*a2 + b1*b1 + b2*b2) / det_J;
		double area_e = 0.5*(det_J + 1.0 / det_J);
		distortion_e = 0.5*(mips_e + area_e);
		//AMIPS

		double k = s*distortion_e;
		if (k > 60) k = 60;
		double exp_e = std::exp(k);
		current_dis_energy += exp_e;



		g_a1 = b2 / 4.0 + a1 / (2.0 * det_J) - b2 / (4.0 * det_J*det_J) - (b2*mips_e) / (2.0*det_J);
		g_a2 = a2 / (2.0 * det_J) - b1 / 4.0 + b1 / (4.0 * det_J*det_J) + (b1*mips_e) / (2.0*det_J);
		g_b1 = a2 / (4.0 * det_J*det_J) - a2 / 4.0 + b1 / (2.0 * det_J) + (a2*mips_e) / (2.0*det_J);
		g_b2 = a1 / 4.0 - a1 / (4.0 * det_J*det_J) + b2 / (2.0 * det_J) - (a1*mips_e) / (2.0*det_J);

		h_a1_a1 = 1.0 / (2.0 * det_J) + b2*b2 / (2.0 * det_J*det_J*det_J) - (a1*b2) / (det_J*det_J) + (b2*b2 * mips_e) / (det_J*det_J);
		h_a1_a2 = (a1*b1) / (2.0 * det_J*det_J) - (a2*b2) / (2.0 * det_J*det_J) - (b1*b2) / (2.0 * det_J*det_J*det_J) - (b1*b2*mips_e) / (det_J*det_J);
		h_a1_b1 = (a1*a2) / (2.0 * det_J*det_J) - (a2*b2) / (2.0 * det_J*det_J*det_J) - (b1*b2) / (2.0 * det_J*det_J) - (a2*b2*mips_e) / (det_J*det_J);
		h_a1_b2 = (a1*b2) / (2.0 * det_J *det_J*det_J) - a1*a1 / (2.0 * det_J*det_J) - b2*b2 / (2.0 * det_J*det_J) - mips_e / (2.0 * det_J) - 1.0 / (4.0 * det_J*det_J) + (a1*b2*mips_e) / (det_J *det_J) + 1.0 / 4.0;

		h_a2_a2 = 1.0 / (2.0 * det_J) + b1 *b1 / (2.0 * det_J*det_J*det_J) + (a2*b1) / (det_J*det_J) + (b1*b1*mips_e) / (det_J*det_J);
		h_a2_b1 = 1.0 / (4.0 * det_J*det_J) + a2*a2 / (2.0 * det_J*det_J) + b1*b1 / (2.0 * det_J*det_J) + mips_e / (2.0 * det_J) + (a2*b1) / (2.0 * det_J *det_J*det_J) + (a2*b1*mips_e) / (det_J *det_J) - 1.0 / 4.0;
		h_a2_b2 = (b1*b2) / (2.0 * det_J*det_J) - (a1*b1) / (2.0 * det_J*det_J*det_J) - (a1*a2) / (2.0 * det_J*det_J) - (a1*b1*mips_e) / (det_J*det_J);

		h_b1_b1 = 1.0 / (2.0 * det_J) + a2*a2 / (2.0 * det_J*det_J*det_J) + (a2*b1) / (det_J*det_J) + (a2*a2 *mips_e) / (det_J*det_J);;
		h_b1_b2 = (a2*b2) / (2.0 * det_J*det_J) - (a1*b1) / (2.0 * det_J*det_J) - (a1*a2) / (2.0 * det_J*det_J*det_J) - (a1*a2*mips_e) / (det_J*det_J);
		h_b2_b2 = 1.0 / (2.0 * det_J) + a1*a1 / (2.0 * det_J*det_J*det_J) - (a1*b2) / (det_J*det_J) + (a1*a1 *mips_e) / (det_J*det_J);



		gradient(6 * var_c_id + 0) += exp_e*s*g_a1;
		gradient(6 * var_c_id + 1) += exp_e*s*g_a2;
		gradient(6 * var_c_id + 2) += exp_e*s*g_b1;
		gradient(6 * var_c_id + 3) += exp_e*s*g_b2;

		hessian.coeffRef(6 * var_c_id + 0, 6 * var_c_id + 0) = exp_e*s*h_a1_a1 + exp_e*s*s*g_a1*g_a1;
		hessian.coeffRef(6 * var_c_id + 0, 6 * var_c_id + 1) = exp_e*s*h_a1_a2 + exp_e*s*s*g_a1*g_a2;
		hessian.coeffRef(6 * var_c_id + 0, 6 * var_c_id + 2) = exp_e*s*h_a1_b1 + exp_e*s*s*g_a1*g_b1;
		hessian.coeffRef(6 * var_c_id + 0, 6 * var_c_id + 3) = exp_e*s*h_a1_b2 + exp_e*s*s*g_a1*g_b2;

		hessian.coeffRef(6 * var_c_id + 1, 6 * var_c_id + 1) = exp_e*s*h_a2_a2 + exp_e*s*s*g_a2*g_a2;
		hessian.coeffRef(6 * var_c_id + 1, 6 * var_c_id + 2) = exp_e*s*h_a2_b1 + exp_e*s*s*g_a2*g_b1;
		hessian.coeffRef(6 * var_c_id + 1, 6 * var_c_id + 3) = exp_e*s*h_a2_b2 + exp_e*s*s*g_a2*g_b2;

		hessian.coeffRef(6 * var_c_id + 2, 6 * var_c_id + 2) = exp_e*s*h_b1_b1 + exp_e*s*s*g_b1*g_b1;
		hessian.coeffRef(6 * var_c_id + 2, 6 * var_c_id + 3) = exp_e*s*h_b1_b2 + exp_e*s*s*g_b1*g_b2;

		hessian.coeffRef(6 * var_c_id + 3, 6 * var_c_id + 3) = exp_e*s*h_b2_b2 + exp_e*s*s*g_b2*g_b2;

		hessian.coeffRef(6 * var_c_id + 1, 6 * var_c_id + 0) = hessian.coeffRef(6 * var_c_id + 0, 6 * var_c_id + 1);
		hessian.coeffRef(6 * var_c_id + 2, 6 * var_c_id + 0) = hessian.coeffRef(6 * var_c_id + 0, 6 * var_c_id + 2);
		hessian.coeffRef(6 * var_c_id + 3, 6 * var_c_id + 0) = hessian.coeffRef(6 * var_c_id + 0, 6 * var_c_id + 3);
		hessian.coeffRef(6 * var_c_id + 2, 6 * var_c_id + 1) = hessian.coeffRef(6 * var_c_id + 1, 6 * var_c_id + 2);
		hessian.coeffRef(6 * var_c_id + 3, 6 * var_c_id + 1) = hessian.coeffRef(6 * var_c_id + 1, 6 * var_c_id + 3);
		hessian.coeffRef(6 * var_c_id + 3, 6 * var_c_id + 2) = hessian.coeffRef(6 * var_c_id + 2, 6 * var_c_id + 3);
	} 
	//if (!is_ivf_below_th)
	{
		hessian += 2.0*ivf_alpha*CTC;
	}

	return ivf_alpha*current_ivf_energy + current_dis_energy;
}

double polycube_boundary_mapping_IVF::compute_only_energy_AT(const std::vector<int>& map_c, const Eigen::Matrix<double, Eigen::Dynamic, 1>& x)
{
	current_dis_energy = 0.0; double s = amips_s; double distortion_e;
	int nc = map_c.size();
	for (int face_id = 0; face_id < nc; ++face_id)
	{
		int var_c_id = map_c[face_id];
		if (var_c_id < 0) continue;
		double a1 = x(6 * var_c_id + 0); double a2 = x(6 * var_c_id + 1);
		double b1 = x(6 * var_c_id + 2); double b2 = x(6 * var_c_id + 3);

		double det_J = a1*b2 - a2*b1;
		if (det_J < 1e-8) // min_det_j = 1e-8
		{
			current_dis_energy = std::numeric_limits<double>::infinity();
			return std::numeric_limits<double>::infinity();
		}

		double mips_e = 0.5*(a1*a1 + a2*a2 + b1*b1 + b2*b2) / det_J;
		
		{
			double area_e = 0.5*(det_J + 1.0 / det_J);
			distortion_e = 0.5*(mips_e + area_e);
		}

		
		{
			double k = s*distortion_e;
			if (k > 60) k = 60;
			double exp_e = std::exp(k);
			current_dis_energy += exp_e;
		}
	}

	current_ivf_energy = x.transpose()*CTC*x;
	current_ivf_energy -= 2.0*x.transpose()*CT*CB;
	current_ivf_energy = std::abs(current_ivf_energy + CBTCB);

	return ivf_alpha*current_ivf_energy + current_dis_energy;
}


void polycube_boundary_mapping_IVF::optimize_AT(const std::vector<int>& map_c )
{
	int iter_count = 0; step_size = 1.0; int max_iter = 120; min_ivf_alpha = 1e3;
	bool is_ivf_below_th = false; pre_current_dis_energy = 0;
	int nf = var_num / 6;
	while (iter_count < max_iter && step_size > min_step_size)
	{
		double old_e = compute_energy_derivative_AT(map_c, solution);

		lltSolver.factorize(hessian);
		int status = lltSolver.info();
		if (status == Eigen::Success)
		{
			lambda = min_lambda;
		}
		else
		{
			int run = 0;
			while (status != Eigen::Success && lambda < max_lambda)
			{
				for (int i = 0; i < nf; i++)
				{
					for (int j = 0; j < 4; ++j)
					{
						hessian.coeffRef(6 * i + j, 6 * i + j) += lambda;
					}
					hessian.coeffRef(6 * i + 4, 6 * i + 4) += min_lambda;
					hessian.coeffRef(6 * i + 5, 6 * i + 5) += min_lambda;
				}

				lltSolver.factorize(hessian);
				status = lltSolver.info();

				if (status != Eigen::Success)
				{
					if (lambda < max_lambda)
					{
						lambda *= 10;
					}
				}
				else
				{
					if (run == 0 && lambda > min_lambda)
					{
						lambda *= 0.1;
					}
				}
				run++;
			}
		}
		if (lambda < max_lambda)
		{
			// compute newton step direction
			step = lltSolver.solve(-gradient);
		}
		else
		{
			// gradient descent
			lambda *= 0.1;
			step = -gradient;
		}

		prevSolution = solution + step_size*step;
		double new_e = compute_only_energy_AT(map_c, prevSolution);
		if (new_e > old_e)
		{
			while (new_e > old_e)
			{
				if (step_size < min_step_size)
				{
					break;
				}
				step_size *= 0.5;
				prevSolution = solution + step_size*step;
				new_e = compute_only_energy_AT(map_c, prevSolution);
			}
		}
		else
		{
			step_size *= 2.0;
		}
		++iter_count;
		//printf("%d A:%e,L:%3.2e,S:%4.3e,I:%e,D:%e\n", iter_count, ivf_alpha, lambda, step_size, current_ivf_energy, current_dis_energy);
		
		solution = prevSolution;
		if (current_ivf_energy < 1e-12 && !is_ivf_below_th)
		{
			//max_iter = iter_count + 15;
			//min_ivf_alpha = 1e-3;
			//amips_s = 1.0;
			is_ivf_below_th = true;
			//compute_only_energy_AT(mesh_, prevSolution);
			//break;
		}
		else if (current_ivf_energy < 1e-12 && is_ivf_below_th)
		{
			if (std::abs(current_dis_energy - pre_current_dis_energy) / pre_current_dis_energy < 1e-8)
			{
				break;
			}
		}

		pre_current_dis_energy = current_dis_energy;
		update_ivf_alpha();
		if (is_ivf_below_th && ivf_alpha < 1e16)
		{
			ivf_alpha = 1e16;
		}

		if (step_size < 1e-12 || lambda > 1e10)
		{
			step_size = 1.0;
			lambda = min_lambda;
			//printf("--------------------------------------\n");
			//printf("reinitilize...............\n");
			min_ivf_alpha = 1e0;
			update_ivf_alpha();
		}
		else
		{
			min_ivf_alpha = 1e3;
		}
	}

	/*for (Mesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		if (!mesh_->data(v_it).get_new_pos_fixed())
		{
			const OpenMesh::Vec3d& p0 = mesh_->point(v_it); OpenMesh::Vec3d np(0, 0, 0); double count = 0;
			for (Mesh::VertexOHalfedgeIter voh_it = mesh_->voh_iter(v_it); voh_it; ++voh_it)
			{
				Mesh::FaceHandle fh = mesh_->face_handle(voh_it);
				if (fh != Mesh::InvalidFaceHandle)
				{
					int face_id = fh.idx();
					double a0 = solution(6 * face_id + 0); double a1 = solution(6 * face_id + 1);
					double b0 = solution(6 * face_id + 2); double b1 = solution(6 * face_id + 3);
					double t0 = solution(6 * face_id + 4); double t1 = solution(6 * face_id + 5);
					double qx = a0*p0[0] + a1*p0[1] + t0;
					double qy = b0*p0[0] + b1*p0[1] + t1;
					np[0] += qx; np[1] += qy; count += 1.0;
				}
			}
			np /= count; mesh_->data(v_it).set_New_Pos(np);
		}
	}*/
}

void polycube_boundary_mapping_IVF::construct_face_vertex_pos(VolumeMesh* mesh_, int chart_label,
	const std::vector<OpenVolumeMesh::Geometry::Vec3d>& src_pos,
	const std::vector< std::vector<int> >& bvf_id,
	const std::vector<int>& map_c, const std::vector<int>& map_v)
{
	int nv = map_v.size();
	for (int i = 0; i < nv; ++i)
	{
		int var_v_id = map_v[i];
		if (var_v_id < 0) continue;

		const std::vector<int>& one_bvf_id = bvf_id[i];
		const OpenVolumeMesh::Geometry::Vec3d& p = src_pos[i];
		double nx = 0.0; double ny = 0.0; double count = 0;
		for (int j = 0; j < one_bvf_id.size(); ++j)
		{
			double x0 = OpenVolumeMesh::Geometry::dot(p - bp[one_bvf_id[j]], fe1[one_bvf_id[j]]);
			double y0 = OpenVolumeMesh::Geometry::dot(p - bp[one_bvf_id[j]], fe2[one_bvf_id[j]]);
			int var_c_id = map_c[one_bvf_id[j]];
			double a0 = solution(6 * var_c_id + 0); double a1 = solution(6 * var_c_id + 1);
			double b0 = solution(6 * var_c_id + 2); double b1 = solution(6 * var_c_id + 3);
			double t0 = solution(6 * var_c_id + 4); double t1 = solution(6 * var_c_id + 5);
			double qx = a0*x0 + a1*y0 + t0;
			double qy = b0*x0 + b1*y0 + t1;

			nx += qx; ny += qy; count += 1.0;
		}
		nx /= count; ny /= count;

		OpenVolumeMesh::VertexHandle vh = OpenVolumeMesh::VertexHandle(i);
		OpenVolumeMesh::Geometry::Vec3d np = mesh_->vertex(vh);
		if (chart_label == 0) //y, z
		{
			np[1] = nx; np[2] = ny;
		}
		else if (chart_label == 1)
		{
			np[1] = -nx; np[2] = ny;
		}
		else if (chart_label == 2) // z, x
		{
			np[2] = nx; np[0] = ny;
		}
		else if (chart_label == 3)
		{
			np[2] = -nx; np[0] = ny;
		}
		else if (chart_label == 4) // x, y
		{
			np[0] = nx; np[1] = ny;
		}
		else if (chart_label == 5)
		{
			np[0] = -nx; np[1] = ny;
		}
		mesh_->set_vertex(vh, np);
	}
}

void polycube_boundary_mapping_IVF::construct_face_vertex_pos(TetStructure<double> *tet_mesh_, int chart_label,
	const std::vector<OpenVolumeMesh::Geometry::Vec3d>& src_pos,
	const std::vector< std::vector<int> >& bvf_id,
	const std::vector<int>& map_c, const std::vector<int>& map_v)
{
	const std::vector<Tetrahedron<double>*> &tetras = tet_mesh_->tetras;
    std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
	int nv = map_v.size();
	for (int i = 0; i < nv; ++i)
	{
		int var_v_id = map_v[i];
		if (var_v_id < 0) continue;

		const std::vector<int>& one_bvf_id = bvf_id[i];
		const OpenVolumeMesh::Geometry::Vec3d& p = src_pos[i];
		double nx = 0.0; double ny = 0.0; double count = 0;
		for (int j = 0; j < one_bvf_id.size(); ++j)
		{
			double x0 = OpenVolumeMesh::Geometry::dot(p - bp[one_bvf_id[j]], fe1[one_bvf_id[j]]);
			double y0 = OpenVolumeMesh::Geometry::dot(p - bp[one_bvf_id[j]], fe2[one_bvf_id[j]]);
			int var_c_id = map_c[one_bvf_id[j]];
			double a0 = solution(6 * var_c_id + 0); double a1 = solution(6 * var_c_id + 1);
			double b0 = solution(6 * var_c_id + 2); double b1 = solution(6 * var_c_id + 3);
			double t0 = solution(6 * var_c_id + 4); double t1 = solution(6 * var_c_id + 5);
			double qx = a0 * x0 + a1 * y0 + t0;
			double qy = b0 * x0 + b1 * y0 + t1;

			nx += qx; ny += qy; count += 1.0;
		}
		nx /= count; ny /= count;

		/*OpenVolumeMesh::VertexHandle vh = OpenVolumeMesh::VertexHandle(i);
		OpenVolumeMesh::Geometry::Vec3d np = mesh_->vertex(vh);*/
		CVec<double, 3> np = tetra_vertices[i]->pos;

		if (chart_label == 0) //y, z
		{
			np[1] = nx; np[2] = ny;
		}
		else if (chart_label == 1)
		{
			np[1] = -nx; np[2] = ny;
		}
		else if (chart_label == 2) // z, x
		{
			np[2] = nx; np[0] = ny;
		}
		else if (chart_label == 3)
		{
			np[2] = -nx; np[0] = ny;
		}
		else if (chart_label == 4) // x, y
		{
			np[0] = nx; np[1] = ny;
		}
		else if (chart_label == 5)
		{
			np[0] = -nx; np[1] = ny;
		}
		//mesh_->set_vertex(vh, np);
		tetra_vertices[i]->pos = np;
	}
}


void polycube_boundary_mapping_IVF::optimize_IVF(VolumeMesh* mesh_, int chart_label,
	const std::vector<OpenVolumeMesh::Geometry::Vec3d>& src_pos, 
	const std::vector<OpenVolumeMesh::Geometry::Vec3i>& bfv_id, const std::vector< std::vector<int> >& bvf_id,
	const std::vector<int>& map_v, const int& var_v_count, const std::vector<int>& map_c, const int& var_c_count)
{
	bool no_flipped = initilize_affine_trans(mesh_,chart_label,src_pos,bfv_id,bvf_id,map_v,var_v_count,map_c,var_c_count);
	if (no_flipped) return;

	compute_only_energy_AT(map_c, solution);
	if (current_dis_energy < 1e-12) return;
	update_ivf_alpha();
	//long start_t = clock();

	optimize_AT(map_c);

	/*long end_t = clock();
	long diff = end_t - start_t;
	double t = (double)(diff) / CLOCKS_PER_SEC;*/
	//printf("Optimize AT Time: %f s\n", t);

	construct_face_vertex_pos(mesh_, chart_label, src_pos, bvf_id, map_c, map_v);
}

void polycube_boundary_mapping_IVF::optimize_IVF(TetStructure<double> *tet_mesh_, int chart_label,
	const std::vector<OpenVolumeMesh::Geometry::Vec3d>& src_pos,
	const std::vector<OpenVolumeMesh::Geometry::Vec3i>& bfv_id, const std::vector< std::vector<int> >& bvf_id,
	const std::vector<int>& map_v, const int& var_v_count, const std::vector<int>& map_c, const int& var_c_count)
{
	bool no_flipped = initilize_affine_trans(tet_mesh_, chart_label, src_pos, bfv_id, bvf_id, map_v, var_v_count, map_c, var_c_count);
	if (no_flipped) return;

	compute_only_energy_AT(map_c, solution);
	if (current_dis_energy < 1e-12) return;
	update_ivf_alpha();
	//long start_t = clock();

	optimize_AT(map_c);

	/*long end_t = clock();
	long diff = end_t - start_t;
	double t = (double)(diff) / CLOCKS_PER_SEC;*/
	//printf("Optimize AT Time: %f s\n", t);

	construct_face_vertex_pos(tet_mesh_, chart_label, src_pos, bvf_id, map_c, map_v);
}
