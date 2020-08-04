#include "Defined_Tensor.h"
#include <Eigen\Dense>
#include <queue>

defined_tensor::defined_tensor()
{
	m_tensor_name = Utility_Definition::T_IDENTITY;
	AABB_Tree = NULL;
	AABB_Segment_Tree = NULL;
	bm = NULL;
	ref_mesh_ = NULL;
	kdTree = NULL;
}

defined_tensor::~defined_tensor()
{
	if(AABB_Tree) delete AABB_Tree;
	if(AABB_Segment_Tree) delete AABB_Segment_Tree;
	if(bm) delete bm;
	if(kdTree) {delete kdTree; annClose();}
}

void defined_tensor::get_tensor( const OpenVolumeMesh::Geometry::Vec3d& p, OpenVolumeMesh::Geometry::Vec6d& M )
{
	double det = 1.0; double t = 1.0;
	switch (m_tensor_name)
	{
	case Utility_Definition::T_IDENTITY:
		M[0] = 1.0; M[1] = 0.0; M[2] = 0.0; M[3] = 1.0; M[4] = 0.0; M[5] = 1.0;
		break;
	case Utility_Definition::T_10_1_1:
		M[0] = 1.0; M[1] = 0.0; M[2] = 0.0; M[3] = 1.0; M[4] = 0.0; M[5] = 100.0;
		break;
	case Utility_Definition::T_X2_Y2_Z4: //0.5*(x^2 + y^2 + z^4)
		M[0] = 1.0; M[1] = 0.0; M[2] = 0.0; M[3] = 1.0; M[4] = 0.0; M[5] = 6*p[2]*p[2];
		det = M[0]*M[3]*M[5] + M[1]*M[4]*M[2] + M[2]*M[1]*M[4] - M[2]*M[3]*M[2] - M[4]*M[4]*M[1]- M[5]*M[1]*M[1];
		t = std::pow( det, -0.2);
		M *= t;
		break;
	case Utility_Definition::T_EXP:
		M[0] = std::exp(p[0]*p[0]/10 + p[1]*p[1]/10 + p[2]*p[2]/10)/5 + (p[0]*p[0]*exp(p[0]*p[0]/10 + p[1]*p[1]/10 + p[2]*p[2]/10))/25;
		M[3] = std::exp(p[0]*p[0]/10 + p[1]*p[1]/10 + p[2]*p[2]/10)/5 + (p[1]*p[1]*exp(p[0]*p[0]/10 + p[1]*p[1]/10 + p[2]*p[2]/10))/25;
		M[5] = std::exp(p[0]*p[0]/10 + p[1]*p[1]/10 + p[2]*p[2]/10)/5 + (p[2]*p[2]*exp(p[0]*p[0]/10 + p[1]*p[1]/10 + p[2]*p[2]/10))/25;
		M[1] = (p[0]*p[1]*exp(p[0]*p[0]/10 + p[1]*p[1]/10 + p[2]*p[2]/10))/25;
		M[4] = (p[2]*p[1]*exp(p[0]*p[0]/10 + p[1]*p[1]/10 + p[2]*p[2]/10))/25;
		M[2] = (p[2]*p[0]*exp(p[0]*p[0]/10 + p[1]*p[1]/10 + p[2]*p[2]/10))/25;
		det = M[0]*M[3]*M[5] + M[1]*M[4]*M[2] + M[2]*M[1]*M[4] - M[2]*M[3]*M[2] - M[4]*M[4]*M[1]- M[5]*M[1]*M[1];
		t = std::pow( det, -0.2);
		M *= t;
		break;
	default:
		M[0] = 1.0; M[1] = 0.0; M[2] = 0.0; M[3] = 1.0; M[4] = 0.0; M[5] = 1.0;
		break;
	}

	if(Utility_Definition::T_SPHERICAL_SHOCK == m_tensor_name)
	{
		double t = 7;
		Eigen::Matrix3d R;
		double r = p.norm();
		OpenVolumeMesh::Geometry::Vec3d e_r = p / r;
		local_frame lf; lf.n = e_r; lf.find_e_x_y();
		R(0,0) =    e_r[0]; R(0,1) =    e_r[1]; R(0,2) =    e_r[2];
		R(1,0) = lf.e_y[0]; R(1,1) = lf.e_y[1]; R(1,2) = lf.e_y[2];
		R(2,0) = lf.e_x[0]; R(2,1) = lf.e_x[1]; R(2,2) = lf.e_x[2];

		Eigen::Matrix3d D; D.setZero();
		D(0,0) = 1.0/( 0.025 + 1.0 * (1.0 - std::exp( -0.01*std::abs(r*r - t*t) ) ) );
		D(1,1) = 2.0 / 2.0; D(2,2) = 2.0 / 2.0;
		Eigen::Matrix3d TM  = R.transpose() *D * D * R;

		M[0] = TM(0,0); M[1] = TM(0,1); M[2] = TM(0,2);
		M[3] = TM(1,1); M[4] = TM(1,2); M[5] = TM(2,2);
	}
	else if (Utility_Definition::T_CYLINDER_SHOCK == m_tensor_name)
	{
		double t = 7; Eigen::Matrix3d R;
		double r = std::sqrt(p[0]*p[0] + p[1]*p[1]);
		R(0,0) =  p[0]/r; R(0,1) = p[1]/r; R(0,2) = 0;
		R(1,0) = -p[1]/r; R(1,1) = p[0]/r; R(1,2) = 0;
		R(2,0) = 0; R(2,1) = 0; R(2,2) = 1;

		Eigen::Matrix3d D; D.setZero();
		D(0,0) = 2.0/( 0.1 + 2.0 * (1.0 - std::exp( -0.01*std::abs(r*r - t*t) ) ) );
		D(1,1) = 2.0 / 2.0; D(2,2) = 2.0 / 2.0;
		Eigen::Matrix3d TM  = R.transpose() *D * D * R;

		M[0] = TM(0,0); M[1] = TM(0,1); M[2] = TM(0,2);
		M[3] = TM(1,1); M[4] = TM(1,2); M[5] = TM(2,2);
	}

	if(Utility_Definition::T_PLANAR_SHOCK == m_tensor_name)
	{
		double t = 0.6;
		Eigen::Matrix3d D; D.setZero();
		D(0,0) = 1.0/( 0.0025 + 0.2 * (1.0 - std::exp( -1.0*std::abs(p[0] - t) ) ) );
		D(1,1) = 1.0 / 0.2; D(2,2) = 1.0 / 0.2;
		Eigen::Matrix3d TM = D * D;

		M[0] = TM(0,0); M[1] = TM(0,1); M[2] = TM(0,2);
		M[3] = TM(1,1); M[4] = TM(1,2); M[5] = TM(2,2);
	}

	if(Utility_Definition::T_SINE == m_tensor_name)
	{
		//printf("Sine tensor!!!!!!!\n");
		Eigen::Matrix3d R;
		OpenVolumeMesh::Geometry::Vec3d e_1( 2.0*std::cos(p[0]*6.0), 1.0, 0.0 ); e_1.normalize();
		OpenVolumeMesh::Geometry::Vec3d e_2(0.0, 1.0, 0.0);
		e_2 = OpenVolumeMesh::Geometry::cross(e_1, e_2); e_2.normalize();
		OpenVolumeMesh::Geometry::Vec3d e_3 = OpenVolumeMesh::Geometry::cross(e_1, e_2);
		
		R(0, 0) = e_1[0]; R(0, 1) = e_1[1]; R(0, 2) = e_1[2];
		R(1, 0) = e_2[0]; R(1, 1) = e_2[1]; R(1, 2) = e_2[2];
		R(2, 0) = e_3[0]; R(2, 1) = e_3[1]; R(2, 2) = e_3[2];

		//double R_det = R.determinant();

		Eigen::Matrix3d D; D.setZero();
		D(0, 0) = 100;
		D(1, 1) = 1.0; D(2, 2) = 1.0;
		D *= 10;
		Eigen::Matrix3d TM = R.transpose()* D * R;

		M[0] = TM(0, 0); M[1] = TM(0, 1); M[2] = TM(0, 2);
		M[3] = TM(1, 1); M[4] = TM(1, 2); M[5] = TM(2, 2);

		

		double D00 = 1000; double D11 = 10.0; double D22 = 10.0;

		/*M[0] = D00*e_1[0] * e_1[0] + D11*e_2[0] * e_2[0] + D22*e_3[0] * e_3[0];
		M[1] = e_1[0] * e_1[1] * D00 + e_2[0] * e_2[1] * D11 + e_3[0] * e_3[1] * D22;
		M[2] = e_1[0] * e_1[2] * D00 + e_2[0] * e_2[2] * D11 + e_3[0] * e_3[2] * D22;
		M[3] = D00*e_1[1] * e_1[1] + D11*e_2[1] * e_2[1] + D22*e_3[1] * e_3[1];
		M[4] = e_1[1] * e_1[2] * D00 + e_2[1] * e_2[2] * D11 + e_3[1] * e_3[2] * D22;
		M[5] = D00*e_1[2] * e_1[2] + D11*e_2[2] * e_2[2] + D22*e_3[2] * e_3[2];*/

		double M0 = D00*e_1[0] * e_1[0] + D11*e_2[0] * e_2[0] + D22*e_3[0] * e_3[0];
		double M1 = e_1[0] * e_1[1] * D00 + e_2[0] * e_2[1] * D11 + e_3[0] * e_3[1] * D22;
		double M2 = e_1[0] * e_1[2] * D00 + e_2[0] * e_2[2] * D11 + e_3[0] * e_3[2] * D22;
		double M3 = D00*e_1[1] * e_1[1] + D11*e_2[1] * e_2[1] + D22*e_3[1] * e_3[1];
		double M4 = e_1[1] * e_1[2] * D00 + e_2[1] * e_2[2] * D11 + e_3[1] * e_3[2] * D22;
		double M5 = D00*e_1[2] * e_1[2] + D11*e_2[2] * e_2[2] + D22*e_3[2] * e_3[2];

		/*M[0] = D00*e_1[0] * e_1[0] + D11*e_1[1] * e_1[1] + D22*e_1[2] * e_1[2];
		M[1] = e_1[0] * e_2[0] * D00 + e_1[1] * e_2[1] * D11 + e_1[2] * e_2[2] * D22;
		M[2] = e_1[0] * e_3[0] * D00 + e_1[1] * e_3[1] * D11 + e_1[2] * e_3[2] * D22;
		M[3] = D00*e_2[0] * e_2[0] + D11*e_2[1] * e_2[1] + D22*e_2[2] * e_2[2];
		M[4] = e_2[0] * e_3[0] * D00 + e_2[1] * e_3[1] * D11 + e_2[2] * e_3[2] * D22;
		M[5] = D00*e_3[0] * e_3[0] + D11*e_3[1] * e_3[1] + D22*e_3[2] * e_3[2];*/
	}

	if(Utility_Definition::T_SINK == m_tensor_name)
	{
		Eigen::Matrix3d R;
		double r = std::sqrt(p[0]*p[0] + p[1]*p[1]);
		R(0,0) =  p[0]/r; R(0,1) = p[1]/r; R(0,2) = 0;
		R(1,0) = -p[1]/r; R(1,1) = p[0]/r; R(1,2) = 0;
		R(2,0) = 0; R(2,1) = 0; R(2,2) = 1;

		Eigen::Matrix3d D; D.setZero();
		D(0,0) = 1.0;
		D(1,1) = 1.0 + r*5;
		D(2,2) = 1.0 + r*5;
		D *= 1.2/(0.5 + 1 - std::exp(-0.05*(r*r - 2.56)) );
		Eigen::Matrix3d TM = R.transpose() *D * D * R ;

		M[0] = TM(0,0); M[1] = TM(0,1); M[2] = TM(0,2);
		M[3] = TM(1,1); M[4] = TM(1,2); M[5] = TM(2,2);
	}

	if(Utility_Definition::T_TEST == m_tensor_name)
	{
		double t = 4;
		Eigen::Matrix3d R;
		double r = p.norm();
		OpenVolumeMesh::Geometry::Vec3d e_r = p / r;
		local_frame lf; lf.n = e_r; lf.find_e_x_y();
		R(0,0) =    e_r[0]; R(0,1) =    e_r[1]; R(0,2) =    e_r[2];
		R(1,0) = lf.e_y[0]; R(1,1) = lf.e_y[1]; R(1,2) = lf.e_y[2];
		R(2,0) = lf.e_x[0]; R(2,1) = lf.e_x[1]; R(2,2) = lf.e_x[2];

		Eigen::Matrix3d D; D.setZero();
		D(0,0) = 1.0/( 0.05 + 1.0 * (1.0 - std::exp( -0.01*std::abs(r*r - t*t) ) ) );
		D(1,1) = 1; D(2,2) = 1;
		D *= 1.25;
		Eigen::Matrix3d TM = R.transpose() *D * D * R ;

		M[0] = TM(0,0); M[1] = TM(0,1); M[2] = TM(0,2);
		M[3] = TM(1,1); M[4] = TM(1,2); M[5] = TM(2,2);
	}
}

void defined_tensor::build_AABB_Tree_using_ref_mesh(SurfaceMesh* mesh_)
{
	if(AABB_Tree) {delete AABB_Tree; AABB_Tree = NULL;}

	unsigned nf = mesh_->n_faces();
	if (nf == 0) return;
	opp_face_id.resize(nf); face_pos.resize(nf);
	std::vector<int> one_face_opp_face_id(3); std::vector<OpenVolumeMesh::Geometry::Vec3d> one_face_pos(3);
	visited_ref_face_id.resize(nf, -1); visited_ref_face.resize(nf, -1);

	std::vector<CGAL_double_3_Point> v_pos( mesh_->n_vertices() ); OpenMesh::Vec3d p;
	for(SurfaceMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		int vertex_id = v_it.handle().idx();
		p = mesh_->point( v_it );
		v_pos[vertex_id] = CGAL_double_3_Point( p[0], p[1], p[2] );
	}

	triangle_vectors.clear();
	triangle_vectors.resize(nf); int fv_id[3];
	for(SurfaceMesh::FaceIter f_it = mesh_->faces_begin(); f_it != mesh_->faces_end(); ++f_it)
	{
		SurfaceMesh::FaceVertexIter fv_it = mesh_->fv_iter(f_it);
		fv_id[0] = fv_it.handle().idx(); ++fv_it;
		fv_id[1] = fv_it.handle().idx(); ++fv_it;
		fv_id[2] = fv_it.handle().idx(); ++fv_it;
		one_face_pos[0] = OpenVolumeMesh::Geometry::Vec3d(v_pos[fv_id[0]].x(), v_pos[fv_id[0]].y(), v_pos[fv_id[0]].z());
		one_face_pos[1] = OpenVolumeMesh::Geometry::Vec3d(v_pos[fv_id[1]].x(), v_pos[fv_id[1]].y(), v_pos[fv_id[1]].z());
		one_face_pos[2] = OpenVolumeMesh::Geometry::Vec3d(v_pos[fv_id[2]].x(), v_pos[fv_id[2]].y(), v_pos[fv_id[2]].z());
		face_pos[f_it.handle().idx()] = one_face_pos;
		SurfaceMesh::FaceHalfedgeIter fhe_it = mesh_->fh_iter(f_it);
		one_face_opp_face_id[0] = mesh_->face_handle(mesh_->opposite_halfedge_handle(fhe_it)).idx();
		++fhe_it; one_face_opp_face_id[1] = mesh_->face_handle(mesh_->opposite_halfedge_handle(fhe_it)).idx();
		++fhe_it; one_face_opp_face_id[2] = mesh_->face_handle(mesh_->opposite_halfedge_handle(fhe_it)).idx();
		opp_face_id[f_it.handle().idx()] = one_face_opp_face_id;

		triangle_vectors[f_it.handle().idx()] = CGAL_3_Triangle(v_pos[fv_id[0]], v_pos[fv_id[1]], v_pos[fv_id[2]]);

		if (fv_it)
		{
			fv_id[1] = fv_id[2];
			fv_id[2] = fv_it.handle().idx();
			triangle_vectors.push_back(CGAL_3_Triangle(v_pos[fv_id[0]], v_pos[fv_id[1]], v_pos[fv_id[2]]));
		}
	}
	AABB_Tree = new CGAL_AABB_Tree(triangle_vectors.begin(), triangle_vectors.end());
	AABB_Tree->accelerate_distance_queries();

	printf("Finish Constructing AABB Tree.\n");
}

void defined_tensor::build_AABB_Tree_using_ref_mesh(VolumeMesh* mesh_)
{
	if (AABB_Tree) { delete AABB_Tree; AABB_Tree = NULL; }
	triangle_vectors.clear();
	for (OpenVolumeMesh::FaceIter f_it = mesh_->faces_begin(); f_it != mesh_->faces_end(); ++f_it)
	{
		if (!mesh_->is_boundary(*f_it)) continue;

		OpenVolumeMesh::HalfFaceHandle hfHandle = mesh_->halfface_handle(*f_it, 0);
		OpenVolumeMesh::HalfFaceHandle hfHandle2 = mesh_->halfface_handle(*f_it, 1);
		int cell_id = mesh_->incident_cell(hfHandle2).idx();
		if (cell_id < 0)
		{
			hfHandle2 = hfHandle;
			cell_id = mesh_->incident_cell(hfHandle2).idx();
			hfHandle = mesh_->halfface_handle(*f_it, 1);
		}

		OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfHandle);
		OpenVolumeMesh::Geometry::Vec3d p0 = mesh_->vertex(*hfv_it);
		++hfv_it; OpenVolumeMesh::Geometry::Vec3d p1 = mesh_->vertex(*hfv_it);
		++hfv_it; OpenVolumeMesh::Geometry::Vec3d p2 = mesh_->vertex(*hfv_it);
		
		CGAL_double_3_Point q0(p0[0], p0[1], p0[2]);
		CGAL_double_3_Point q1(p1[0], p1[1], p1[2]);
		CGAL_double_3_Point q2(p2[0], p2[1], p2[2]);

		triangle_vectors.push_back( CGAL_3_Triangle(q0, q1, q2) );
	}

	AABB_Tree = new CGAL_AABB_Tree(triangle_vectors.begin(), triangle_vectors.end());
	AABB_Tree->accelerate_distance_queries();

	printf("Finish Constructing AABB Tree.\n");
}

void defined_tensor::build_AABB_Segment_Tree_using_ref_mesh(VolumeMesh* mesh_, const std::vector<int>& edge_feature)
{
	if(AABB_Segment_Tree) {delete AABB_Segment_Tree; AABB_Segment_Tree = NULL;}
	if(mesh_->n_edges() == 0) return;

	std::vector<CGAL_3_Segment> segment_vectors;
	for(unsigned i=0;i<edge_feature.size();++i)
	{
		if(edge_feature[i] == 1)
		{
			OpenVolumeMesh::EdgeHandle eh(i);
			OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(eh);
			OpenVolumeMesh::VertexHandle vh_from = edge.from_vertex();
			OpenVolumeMesh::VertexHandle vh_to = edge.to_vertex();
			OpenVolumeMesh::Geometry::Vec3d p_from = mesh_->vertex(vh_from);
			OpenVolumeMesh::Geometry::Vec3d p_to = mesh_->vertex(vh_to);
			CGAL_double_3_Point pf(p_from[0], p_from[1], p_from[2]);
			CGAL_double_3_Point pt(p_to[0], p_to[1], p_to[2]);
			segment_vectors.push_back( CGAL_3_Segment(pf, pt) );
		}
	}
	AABB_Segment_Tree = new CGAL_AABB_Segment_Tree(segment_vectors.begin(), segment_vectors.end());
	AABB_Segment_Tree->accelerate_distance_queries();
	printf("Finish Constructing AABB Segment Tree.\n");
}

void defined_tensor::project_on_ref_mesh(OpenVolumeMesh::Geometry::Vec3d& p)
{
	CGAL_double_3_Point pos = AABB_Tree->closest_point( CGAL_double_3_Point(p[0], p[1], p[2]) );
	p = OpenVolumeMesh::Geometry::Vec3d( pos.x(), pos.y(), pos.z() );
}

//if input face id = -1, will change face_id , p
void defined_tensor::project_on_ref_mesh_with_guess_face(OpenVolumeMesh::Geometry::Vec3d& p, int& guess_face_id, const double& radius)
{
	if (guess_face_id == -1)
	{
		CGAL_AABB_Tree::Point_and_primitive_id point_primitive = AABB_Tree->closest_point_and_primitive(CGAL_double_3_Point(p[0], p[1], p[2]));
		CGAL_double_3_Point pos = point_primitive.first;
		CGAL_Triangle_Iterator it = point_primitive.second;
		guess_face_id = std::distance(triangle_vectors.begin(), it);
		p = OpenVolumeMesh::Geometry::Vec3d(pos.x(), pos.y(), pos.z());
	}
	else
	{
		int visited_face_count = 0; double* p_data = p.data();
		visited_ref_face_id[0] = guess_face_id; visited_ref_face[guess_face_id] = 1; ++visited_face_count;
		OpenVolumeMesh::Geometry::Vec3d np; OpenVolumeMesh::Geometry::Vec3d np_; double min_d = 1.0e30; int min_face_id = guess_face_id;
		double distance_th = 25.0*radius; int current_face_count = 0;
		while (current_face_count < visited_face_count)
		{
			//int face_id = Q.front(); Q.pop();
			int face_id = visited_ref_face_id[current_face_count]; ++current_face_count;
			std::vector<OpenVolumeMesh::Geometry::Vec3d>& one_face_pos = face_pos[face_id];
			//double d = distPointTriangleSquared(p, one_face_pos[0], one_face_pos[1], one_face_pos[2], np);
			double d = distPointTriangleSquared2(p_data, one_face_pos[0].data(), one_face_pos[1].data(), one_face_pos[2].data(), np.data());
			if (d < min_d) { min_d = d; min_face_id = face_id; np_ = np; }

			std::vector<int>& one_face_opp_face_id = opp_face_id[face_id];
			for (unsigned i = 0; i < 3; ++i)
			{
				if (one_face_opp_face_id[i] >= 0 && visited_ref_face[one_face_opp_face_id[i]] == -1)
				{
					std::vector<OpenVolumeMesh::Geometry::Vec3d>& opp_face_pos = face_pos[one_face_opp_face_id[i]];
					OpenVolumeMesh::Geometry::Vec3d& opp_face_pos1 = opp_face_pos[0];
					double d1 = (p[0] - opp_face_pos1[0])*(p[0] - opp_face_pos1[0]) + (p[1] - opp_face_pos1[1])*(p[1] - opp_face_pos1[1]) + (p[2] - opp_face_pos1[2])*(p[2] - opp_face_pos1[2]);
					if (d1 < distance_th)
					{
						visited_ref_face_id[visited_face_count] = one_face_opp_face_id[i];
						visited_ref_face[one_face_opp_face_id[i]] = 1;
						++visited_face_count;
					}
					else
					{
						OpenVolumeMesh::Geometry::Vec3d& opp_face_pos2 = opp_face_pos[1];
						double d2 = (p[0] - opp_face_pos2[0])*(p[0] - opp_face_pos2[0]) + (p[1] - opp_face_pos2[1])*(p[1] - opp_face_pos2[1]) + (p[2] - opp_face_pos2[2])*(p[2] - opp_face_pos2[2]);
						if (d2 < distance_th)
						{
							visited_ref_face_id[visited_face_count] = one_face_opp_face_id[i];
							visited_ref_face[one_face_opp_face_id[i]] = 1;
							++visited_face_count;
						}
						else
						{
							OpenVolumeMesh::Geometry::Vec3d& opp_face_pos3 = opp_face_pos[2];
							double d3 = (p[0] - opp_face_pos3[0])*(p[0] - opp_face_pos3[0]) + (p[1] - opp_face_pos3[1])*(p[1] - opp_face_pos3[1]) + (p[2] - opp_face_pos3[2])*(p[2] - opp_face_pos3[2]);
							if (d3 < distance_th)
							{
								visited_ref_face_id[visited_face_count] = one_face_opp_face_id[i];
								visited_ref_face[one_face_opp_face_id[i]] = 1;
								++visited_face_count;
							}
						}
					}
				}
			}
		}

		for (unsigned i = 0; i < visited_face_count; ++i)
		{
			visited_ref_face[visited_ref_face_id[i]] = -1;
		}

		guess_face_id = min_face_id;
		p = np_;
	}
	
}

void defined_tensor::project_on_ref_mesh_with_metric(OpenVolumeMesh::Geometry::Vec3d& p, OpenVolumeMesh::Geometry::Vec6d& M)
{
	CGAL_AABB_Tree::Point_and_primitive_id point_primitive = AABB_Tree->closest_point_and_primitive( CGAL_double_3_Point(p[0], p[1], p[2]) );
	CGAL_double_3_Point pos = point_primitive.first;
	CGAL_Triangle_Iterator it = point_primitive.second;
	unsigned face_id = std::distance( triangle_vectors.begin(), it);
	p = OpenVolumeMesh::Geometry::Vec3d( pos.x(), pos.y(), pos.z() );

	//printf("F : %d\n", face_id);

	OpenMesh::Vec3d p_(p[0], p[1], p[2]);
	SurfaceMesh::ConstFaceVertexIter fv_it = bm->cfv_iter( bm->face_handle(face_id) );
	OpenMesh::Vec3d p0 = bm->point(fv_it.handle()); int v0_id = fv_it.handle().idx();
	++fv_it; OpenMesh::Vec3d p1 = bm->point(fv_it.handle()); int v1_id = fv_it.handle().idx();
	++fv_it; OpenMesh::Vec3d p2 = bm->point(fv_it.handle()); int v2_id = fv_it.handle().idx();

	OpenMesh::Vec3d bc;
	if (!baryCoord(p_, p0, p1, p2, bc))
		bc[0] = bc[1] = bc[2] = 1.0/3.0;

	M = bc[0] *bm_metric[v0_id] + bc[1] *bm_metric[v1_id] + bc[2] *bm_metric[v2_id];
}

void defined_tensor::project_on_ref_edge(OpenVolumeMesh::Geometry::Vec3d& p)
{
	CGAL_double_3_Point pos = AABB_Segment_Tree->closest_point( CGAL_double_3_Point(p[0], p[1], p[2]) );
	p = OpenVolumeMesh::Geometry::Vec3d( pos.x(), pos.y(), pos.z() );
}

void defined_tensor::construct_boundary_mesh()
{
	OpenVolumeMesh::BoundaryFaceIter bf_it(ref_mesh_->bf_iter());
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = ref_mesh_->hfv_iter( ref_mesh_->halfface_handle(*bf_it,0) );

	std::vector<SurfaceMesh::VertexHandle> vertexHandleVec;
	std::vector<SurfaceMesh::VertexHandle> faceVertexHandle;
	std::map<int,int> VolumeSurfaceVertex;
	std::map<int,int>::iterator map_it;
	int indexOnSurface = 0;
	int InvalidHF =0;
	bm->clear();
	for (; bf_it.valid(); ++bf_it)
	{
		OpenVolumeMesh::HalfFaceHandle hfHandle = ref_mesh_->halfface_handle(*bf_it,0);
		OpenVolumeMesh::HalfFaceHandle hfHandle2 = ref_mesh_->halfface_handle(*bf_it,1);
		int cell_id = ref_mesh_->incident_cell(hfHandle2).idx();
		if( cell_id >= 0 )
		{
			hfHandle2 = hfHandle;
			cell_id = ref_mesh_->incident_cell(hfHandle2).idx();
			hfHandle = ref_mesh_->halfface_handle(*bf_it,1);
		}
		++InvalidHF;
		hfv_it = ref_mesh_->hfv_iter( hfHandle );
		faceVertexHandle.clear();
		for( hfv_it; hfv_it.valid(); ++hfv_it)
		{
			map_it = VolumeSurfaceVertex.find(hfv_it->idx()) ;
			if( map_it == VolumeSurfaceVertex.end() )
			{
				VolumeSurfaceVertex.insert(std::pair<int,int>( hfv_it->idx(), indexOnSurface ));
				OpenVolumeMesh::Geometry::Vec3d v1 = ref_mesh_->vertex(*hfv_it);
				SurfaceMesh::Point v2( v1[0], v1[1], v1[2] );
				vertexHandleVec.push_back( bm->add_vertex( v2 ) );
				bm->data(vertexHandleVec.back()).set_tet_vertex_id(hfv_it->idx());
				faceVertexHandle.push_back( vertexHandleVec.back() );
				++indexOnSurface;
			}
			else
			{
				faceVertexHandle.push_back(vertexHandleVec[map_it->second]);
			}
		}

		std::reverse(faceVertexHandle.begin(), faceVertexHandle.end());
		OpenMesh::FaceHandle fh = bm->add_face(faceVertexHandle);
		bm->data(fh).set_cell_id( cell_id );
	}

	bm->update_normals();
}

void defined_tensor::save_boundary_mesh_metric()
{
	FILE* f_b_m = fopen("A:\\Code\\VolumeMeshProcessing\\Models\\Tetrahedrization\\elphant\\boundary.off", "w");
	unsigned nv = bm->n_vertices(); unsigned nf = bm->n_faces();
	fprintf(f_b_m, "OFF");
	fprintf(f_b_m, "\n%d %d 0", nv, nf);
	for(SurfaceMesh::VertexIter v_it = bm->vertices_begin(); v_it != bm->vertices_end(); ++v_it)
	{
		OpenMesh::Vec3d p = bm->point(v_it);
		fprintf(f_b_m, "\n%20.19f %20.19f %20.19f", p[0], p[1], p[2]);
	}
	for(SurfaceMesh::FaceIter f_it = bm->faces_begin(); f_it != bm->faces_end(); ++f_it)
	{
		fprintf(f_b_m, "\n3");
		for( SurfaceMesh::FaceVertexIter fv_it = bm->fv_iter(f_it); fv_it; ++fv_it )
		{
			fprintf(f_b_m, " %d", fv_it.handle().idx());
		}
	}
	for(unsigned i=0;i<bm_metric.size();++i)
	{
		fprintf(f_b_m, "\n%20.19f %20.19f %20.19f %20.19f %20.19f %20.19f", bm_metric[i][0],bm_metric[i][1],bm_metric[i][2],bm_metric[i][3],bm_metric[i][4],bm_metric[i][5]);
	}

	fclose(f_b_m);
}

void defined_tensor::smmoth_metric_to_interior(VolumeMesh* mesh_)
{
	if(ref_mesh_) delete ref_mesh_;
	ref_mesh_ = mesh_;

	if(bm) delete bm;
	bm = new SurfaceMesh();
	bm->request_face_status();
	bm->request_vertex_status();
	bm->request_edge_status();
	bm->request_face_normals();
	bm->request_vertex_normals();

	construct_boundary_mesh();
	//printf("1\n");

	std::vector<double> k1, k2; std::vector<OpenMesh::Vec3d> d1, d2;
	compute_principal_curvature(bm,k1,k2,d1,d2); bm_metric.resize(k1.size());
	tet_vertex_metric.clear(); unsigned tet_nv = mesh_->n_vertices(); tet_vertex_metric.resize(tet_nv);
	std::vector<int> tet_vertex_flag(tet_nv, -1); std::vector<int> new_one_ring_vertex(tet_nv, -1); int new_one_ring_vertex_count = 0;
	Eigen::Matrix3d R, D, M; D.setZero();
	//printf("2 %d %d %d\n", k1.size(), bm->n_vertices(), tet_nv);
	for(unsigned i=0;i<k1.size();++i)
	{
		SurfaceMesh::VertexHandle vh = bm->vertex_handle(i);
		int tet_v_id = bm->data(vh).get_tet_vertex_id();
		//printf("%d\n", tet_v_id);
		//use curvature to construct metric
		OpenMesh::Vec3d n = OpenMesh::cross(d1[i], d2[i]); n.normalize();
		R(0,0) = d1[i][0]; R(1,0) = d1[i][1]; R(2,0) = d1[i][2];
		R(0,1) = d2[i][0]; R(1,1) = d2[i][1]; R(2,1) = d2[i][2];
		R(0,2) =     n[0]; R(1,2) =     n[1]; R(2,2) =     n[2];
		D(0,0) = std::abs(k1[i]) < 1.0e-4 ? 1.0e-4 : std::abs(k1[i]);
		D(1,1) = std::abs(k2[i]) < 1.0e-4 ? 1.0e-4 : std::abs(k2[i]);
		double max_k = std::abs(k1[i]) > std::abs(k2[i]) ? std::abs(k1[i]) : std::abs(k2[i]);
		D(2,2) = max_k < 1e-4 ? 1e-4 : max_k;
		//D(2,2) = 0.0;
		M = R*D*R.transpose();
		tet_vertex_metric[tet_v_id] = OpenVolumeMesh::Geometry::Vec6d(M(0,0), M(0,1), M(0,2), M(1,1), M(1,2), M(2,2));
		bm_metric[i] = tet_vertex_metric[tet_v_id];
		tet_vertex_flag[tet_v_id] = 1; OpenVolumeMesh::VertexHandle tet_vh(tet_v_id);
		for(OpenVolumeMesh::VertexOHalfEdgeIter voh_it = mesh_->voh_iter(tet_vh); voh_it; ++voh_it)
		{
			OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge( mesh_->edge_handle(*voh_it) );
			OpenVolumeMesh::VertexHandle f_vh = edge.from_vertex();
			OpenVolumeMesh::VertexHandle t_vh = edge.to_vertex();
			if(f_vh != tet_vh && !mesh_->is_boundary(f_vh))
			{
				new_one_ring_vertex[f_vh.idx()] = 1; ++new_one_ring_vertex_count;
			}
			else if(t_vh != tet_vh && !mesh_->is_boundary(t_vh))
			{
				new_one_ring_vertex[t_vh.idx()] = 1; ++new_one_ring_vertex_count;
			}
		}
	}

	//printf("3,1\n");

	save_boundary_mesh_metric();

	//printf("3,2\n");

	//init value
	while(new_one_ring_vertex_count > 0)
	{
		std::vector<int> old_v; old_v.reserve(new_one_ring_vertex_count);
		for(unsigned i=0;i<new_one_ring_vertex.size();++i)
		{
			if(new_one_ring_vertex[i] == 1)
			{
				old_v.push_back(i); OpenVolumeMesh::VertexHandle tet_vh(i);
				OpenVolumeMesh::Geometry::Vec6d ave_M(0,0,0,0,0,0); double sum_count = 0.0;
				for(OpenVolumeMesh::VertexOHalfEdgeIter voh_it = mesh_->voh_iter(tet_vh); voh_it; ++voh_it)
				{
					OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge( mesh_->edge_handle(*voh_it) );
					OpenVolumeMesh::VertexHandle f_vh = edge.from_vertex();
					OpenVolumeMesh::VertexHandle t_vh = edge.to_vertex();
					if(f_vh != tet_vh)
					{
						if(tet_vertex_flag[f_vh.idx()] == 1)
						{
							ave_M += tet_vertex_metric[f_vh.idx()]; sum_count += 1.0;
						}
					}
					else if(t_vh != tet_vh)
					{
						if(tet_vertex_flag[t_vh.idx()] == 1)
						{
							ave_M += tet_vertex_metric[t_vh.idx()]; sum_count += 1.0;
						}
					}
				}

				ave_M /= sum_count;
				tet_vertex_metric[i] = ave_M;
			}
		}

		for(unsigned i=0;i<old_v.size();++i)
		{
			tet_vertex_flag[old_v[i]] = 1; new_one_ring_vertex[old_v[i]] = -1;
		}

		new_one_ring_vertex_count = 0;
		for(unsigned i=0;i<old_v.size();++i)
		{
			OpenVolumeMesh::VertexHandle tet_vh(old_v[i]);
			for(OpenVolumeMesh::VertexOHalfEdgeIter voh_it = mesh_->voh_iter(tet_vh); voh_it; ++voh_it)
			{
				OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge( mesh_->edge_handle(*voh_it) );
				OpenVolumeMesh::VertexHandle f_vh = edge.from_vertex();
				OpenVolumeMesh::VertexHandle t_vh = edge.to_vertex();
				if(f_vh != tet_vh && tet_vertex_flag[f_vh.idx()] == -1)
				{
					new_one_ring_vertex[f_vh.idx()] = 1; ++new_one_ring_vertex_count;
				}
				else if(t_vh != tet_vh && tet_vertex_flag[t_vh.idx()] == -1)
				{
					new_one_ring_vertex[t_vh.idx()] = 1; ++new_one_ring_vertex_count;
				}
			}
		}
	}
	
	//printf("4\n");
	//smooth
	for(unsigned i=0;i<50;++i)
	{
		printf("%d ", i);
		for(OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
		{
			if(mesh_->is_boundary(*v_it)) continue;
			OpenVolumeMesh::Geometry::Vec6d ave_M(0,0,0,0,0,0); double sum_count = 0.0;
			for(OpenVolumeMesh::VertexOHalfEdgeIter voh_it = mesh_->voh_iter(*v_it); voh_it; ++voh_it)
			{
				OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge( mesh_->edge_handle(*voh_it) );
				OpenVolumeMesh::VertexHandle f_vh = edge.from_vertex();
				OpenVolumeMesh::VertexHandle t_vh = edge.to_vertex();
				if( f_vh != *v_it )
				{
					ave_M += tet_vertex_metric[f_vh.idx()]; sum_count += 1.0;
				}
				else if( t_vh != *v_it )
				{
					ave_M += tet_vertex_metric[t_vh.idx()]; sum_count += 1.0;
				}
			}
			ave_M /= sum_count;
			tet_vertex_metric[v_it->idx()] = ave_M;
		}
	}
	printf("\n");
}

void defined_tensor::save_tet_vertex_metric(const char* filename)
{
	FILE* f_v_m = fopen(filename, "w");

	fprintf(f_v_m, "%d", tet_vertex_metric.size());
	for(unsigned i=0;i<tet_vertex_metric.size();++i)
	{
		fprintf(f_v_m, "\n%20.19f %20.19f %20.19f %20.19f %20.19f %20.19f", tet_vertex_metric[i][0],tet_vertex_metric[i][1],tet_vertex_metric[i][2],tet_vertex_metric[i][3],tet_vertex_metric[i][4],tet_vertex_metric[i][5] );
	}

	fclose(f_v_m);
}

void defined_tensor::load_tet_vertex_metric(const char* filename)
{
	FILE* f_v_m = fopen(filename, "r");
	int nv = 0;
	char buf[4096];  fgets(buf, 4096, f_v_m);
	sscanf(buf, "%d", &nv); tet_vertex_metric.clear(); tet_vertex_metric.resize(nv); int v_count = 0;
	char h0[128]; char h1[128]; char h2[128]; char h4[128]; char h5[128]; char h8[128];
	while ( !feof(f_v_m) )
	{
		fgets(buf, 4096, f_v_m);
		sscanf(buf, "%s %s %s %s %s %s", h0, h1, h2, h4, h5, h8);
		tet_vertex_metric[v_count] = OpenVolumeMesh::Geometry::Vec6d(atof(h0), atof(h1),atof(h2), atof(h4), atof(h5), atof(h8));
		++v_count;
	}

	fclose(f_v_m);
}

void defined_tensor::load_ref_tet_mesh(VolumeMesh* mesh_)
{
	if(tet_vertex_metric.size() == 0)
	{
		printf("Please Load Metric First!!!\n");
		return;
	}
	if(ref_mesh_) delete ref_mesh_;
	ref_mesh_ = mesh_;
	if(bm) delete bm;
	bm = new SurfaceMesh();
	construct_boundary_mesh();
	bm_metric.resize(bm->n_vertices());
	for(unsigned i=0;i<bm->n_vertices(); ++i)
	{
		SurfaceMesh::VertexHandle vh = bm->vertex_handle(i);
		int tet_v_id = bm->data(vh).get_tet_vertex_id();
		bm_metric[i] = tet_vertex_metric[tet_v_id];
	}
	build_AABB_Tree_using_ref_mesh(bm);
	build_ANN_KD_tree();
}

void defined_tensor::build_ANN_KD_tree()
{
	int nc = ref_mesh_->n_cells();
	ANNpointArray dataPts = annAllocPts(nc, 3);
	for(OpenVolumeMesh::CellIter c_it = ref_mesh_->cells_begin(); c_it != ref_mesh_->cells_end();++c_it)
	{
		OpenVolumeMesh::Geometry::Vec3d c_c(0,0,0); double count = 0.0;
		for(OpenVolumeMesh::CellVertexIter cv_it = ref_mesh_->cv_iter(*c_it); cv_it; ++cv_it )
		{
			c_c += ref_mesh_->vertex(*cv_it); count += 1.0;
		}
		c_c /= count;
		int c_id = c_it->idx();
		dataPts[c_id][0] = c_c[0]; dataPts[c_id][1] = c_c[1]; dataPts[c_id][2] = c_c[2];
	}
	if(kdTree) delete kdTree;
	kdTree = new ANNkd_tree(dataPts, nc, 3);
}

bool defined_tensor::check_in_tet(const OpenVolumeMesh::Geometry::Vec3d& p, const int& cel_id)
{
	if(cel_id < 0) return false;

	OpenVolumeMesh::CellHandle ch(cel_id);
	OpenVolumeMesh::OpenVolumeMeshCell cell = ref_mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hfh_vec = cell.halffaces();
	for(unsigned i=0;i<hfh_vec.size();++i)
	{
		OpenVolumeMesh::HalfFaceVertexIter hfv_it = ref_mesh_->hfv_iter(hfh_vec[i]);
		OpenVolumeMesh::Geometry::Vec3d p0 = ref_mesh_->vertex(*hfv_it);
		++hfv_it; OpenVolumeMesh::Geometry::Vec3d p1 = ref_mesh_->vertex(*hfv_it);
		++hfv_it; OpenVolumeMesh::Geometry::Vec3d p2 = ref_mesh_->vertex(*hfv_it);

		double v = OpenVolumeMesh::Geometry::dot(p-p0, OpenVolumeMesh::Geometry::cross(p1-p0,p2-p0));
		if(v < 0) return false;
	}
	return true;
}

int defined_tensor::find_nearest_cell_id(const OpenVolumeMesh::Geometry::Vec3d& p)
{
	ANNpoint tp = annAllocPt(3); tp[0] = p[0]; tp[1] = p[1]; tp[2] = p[2];
	ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];
	kdTree->annkSearch(tp, 1, nnIdx, dists);
	int cell_id = nnIdx[0];
	delete [] nnIdx;

	std::vector<int> visited_cell_id(ref_mesh_->n_cells(), -1);
	visited_cell_id[cell_id] = 1;
	std::queue<int> Q; Q.push(cell_id);
	while (Q.size() > 0)
	{
		cell_id = Q.front(); Q.pop();
		if( check_in_tet(p, cell_id) ) return cell_id;
		OpenVolumeMesh::CellHandle ch(cell_id);
		OpenVolumeMesh::OpenVolumeMeshCell cell = ref_mesh_->cell(ch);
		std::vector<OpenVolumeMesh::HalfFaceHandle> hfh_vec = cell.halffaces();
		for(unsigned i=0;i<hfh_vec.size();++i)
		{
			OpenVolumeMesh::HalfFaceHandle hfh = ref_mesh_->opposite_halfface_handle(hfh_vec[i]);
			OpenVolumeMesh::CellHandle t_ch = ref_mesh_->incident_cell(hfh);
			if(t_ch != VolumeMesh::InvalidCellHandle && visited_cell_id[t_ch.idx()] == -1)
			{
				Q.push( t_ch.idx() ); visited_cell_id[t_ch.idx()] = 1;
			}
		}
	}

	return -1;
}

void defined_tensor::find_cell_bay_cor(const OpenVolumeMesh::Geometry::Vec3d& p,const int& cell_id, 
									   OpenVolumeMesh::Geometry::Vec4i& v_id, OpenVolumeMesh::Geometry::Vec4d& b)
{
	OpenVolumeMesh::CellHandle ch(cell_id);
	OpenVolumeMesh::OpenVolumeMeshCell cell = ref_mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hfh_vec = cell.halffaces(); double sum_b= 0.0;
	for(unsigned i=0;i<hfh_vec.size();++i)
	{
		OpenVolumeMesh::HalfFaceVertexIter hfv_it = ref_mesh_->hfv_iter(hfh_vec[i]);
		OpenVolumeMesh::Geometry::Vec3d p0 = ref_mesh_->vertex(*hfv_it); int v0 = hfv_it->idx();
		++hfv_it; OpenVolumeMesh::Geometry::Vec3d p1 = ref_mesh_->vertex(*hfv_it); int v1 = hfv_it->idx();
		++hfv_it; OpenVolumeMesh::Geometry::Vec3d p2 = ref_mesh_->vertex(*hfv_it); int v2 = hfv_it->idx();

		b[i] = OpenVolumeMesh::Geometry::dot(p - p0, OpenVolumeMesh::Geometry::cross(p1-p0,p2-p0));
		sum_b += b[i];
		unsigned j = (i+1) % hfh_vec.size();
		for( hfv_it = ref_mesh_->hfv_iter(hfh_vec[j]); hfv_it; ++hfv_it)
		{
			if(hfv_it->idx() != v0 &&hfv_it->idx() != v1 &&hfv_it->idx() != v2 )
			{
				v_id[i] = hfv_it->idx();
				break;
			}
		}
	}
	b /= sum_b;
}

void defined_tensor::project_for_interior(OpenVolumeMesh::Geometry::Vec3d& p, OpenVolumeMesh::Geometry::Vec6d& M)
{
	//printf("...........\n");
	int cell_id = find_nearest_cell_id(p);
	if(cell_id < 0) printf("%d\n", cell_id);

	OpenVolumeMesh::Geometry::Vec4i v_id; OpenVolumeMesh::Geometry::Vec4d b;
	find_cell_bay_cor(p, cell_id, v_id, b);
	M = tet_vertex_metric[v_id[0]] * b[0] + tet_vertex_metric[v_id[1]] * b[1] + tet_vertex_metric[v_id[2]] * b[2] + tet_vertex_metric[v_id[3]] * b[3];
}