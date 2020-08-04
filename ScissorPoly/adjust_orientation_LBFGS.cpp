#include "adjust_orientation_LBFGS.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
using Eigen::Matrix3d;
using Eigen::Vector3d;
void evalfunc_ad_LBFGS(int N, double* x, double *prev_x, double* f, double* g, void* user_pointer)
{
	double alpha = x[0]; double beta = x[1]; double gamma = x[2];
	double cos_alpha = std::cos(alpha); double sin_alpha = std::sin(alpha);
	double cos_beta = std::cos(beta); double sin_beta = std::sin(beta);
	double cos_gamma = std::cos(gamma); double sin_gamma = std::sin(gamma);
	double r00 = cos_alpha*cos_gamma - cos_beta*sin_alpha*sin_gamma;
	double r01 = sin_alpha*cos_gamma + cos_beta*cos_alpha*sin_gamma;
	double r02 = sin_beta*sin_gamma;
	double r10 = -cos_alpha*sin_gamma - cos_beta*sin_alpha*cos_gamma;
	double r11 = -sin_alpha*sin_gamma + cos_beta*cos_alpha*cos_gamma;
	double r12 = sin_beta*cos_gamma;
	double r20 = sin_beta*sin_alpha;
	double r21 = -sin_beta*cos_alpha;
	double r22 = cos_beta;
	
	user_pointer_ad* up = (user_pointer_ad*)(user_pointer);
	int nbf = (int)up->bfn.size(); *f = 0.0;
	g[0] = 0; g[1] = 0; g[2] = 0;
	for (int i = 0; i < nbf; ++i)
	{
		OpenVolumeMesh::Geometry::Vec3d& n = up->bfn[i];
		//OpenVolumeMesh::Geometry::Vec3d n(4.0, 5.0, 6.0);
		double nx = r00 * n[0] + r01 * n[1] + r02 * n[2];
		double ny = r10 * n[0] + r11 * n[1] + r12 * n[2];
		double nz = r20 * n[0] + r21 * n[1] + r22 * n[2];
		*f += (nx*ny)*(nx*ny) + (ny*nz)*(ny*nz) + (nz*nx)*(nz*nx);
		double d_nx_alpha = (-sin_alpha*cos_gamma - cos_beta*cos_alpha*sin_gamma) * n[0] + (cos_alpha*cos_gamma - cos_beta*sin_alpha*sin_gamma) * n[1] + 0 * n[2];
		double d_ny_alpha = (sin_alpha*sin_gamma - cos_beta*cos_alpha*cos_gamma) * n[0] + (-cos_alpha*sin_gamma - cos_beta*sin_alpha*cos_gamma) * n[1] + 0 * n[2];
		double d_nz_alpha = sin_beta*cos_alpha * n[0] + sin_beta*sin_alpha * n[1] + 0 * n[2];
		double d_nx_beta = (sin_beta*sin_alpha*sin_gamma) * n[0] + (-sin_beta*cos_alpha*sin_gamma) * n[1] + cos_beta*sin_gamma * n[2];
		double d_ny_beta = (sin_beta*sin_alpha*cos_gamma) * n[0] + (-sin_beta*cos_alpha*cos_gamma) * n[1] + cos_beta*cos_gamma * n[2];
		double d_nz_beta = cos_beta*sin_alpha *n[0] - cos_beta*cos_alpha * n[1] - sin_beta * n[2];
		double d_nx_gamma = (-cos_alpha*sin_gamma - cos_beta*sin_alpha*cos_gamma) * n[0] + (-sin_alpha*sin_gamma + cos_beta*cos_alpha*cos_gamma) * n[1] + sin_beta*cos_gamma * n[2];
		double d_ny_gamma = (-cos_alpha*cos_gamma + cos_beta*sin_alpha*sin_gamma) * n[0] + (-sin_alpha*cos_gamma - cos_beta*cos_alpha*sin_gamma) * n[1] + (-sin_beta*sin_gamma) * n[2];
		double d_nz_gamma = 0.0;
		g[0] += 2 * (nx*ny)*(d_nx_alpha*ny + nx*d_ny_alpha) + 2 * (ny*nz)*(d_ny_alpha*nz + ny*d_nz_alpha) + 2 * (nz*nx)*(d_nz_alpha*nx + nz*d_nx_alpha);
		g[1] += 2 * (nx*ny)*(d_nx_beta*ny + nx*d_ny_beta) + 2 * (ny*nz)*(d_ny_beta*nz + ny*d_nz_beta) + 2 * (nz*nx)*(d_nz_beta*nx + nz*d_nx_beta);
		g[2] += 2 * (nx*ny)*(d_nx_gamma*ny + nx*d_ny_gamma) + 2 * (ny*nz)*(d_ny_gamma*nz + ny*d_nz_gamma) + 2 * (nz*nx)*(d_nz_gamma*nx + nz*d_nx_gamma);
		std::cout << "g: " << g[0] << " " << g[1] << " " << g[2] << std::endl;
	}
}
void evalfunc_ad_LBFGS_xyz(int N, double* x, double *prev_x, double* f, double* g, void* user_pointer)
{
	//x y z axis rotation
	double cos_phi = std::cos(x[0]), sin_phi = std::sin(x[0]);
	double cos_theta = std::cos(x[1]), sin_theta = std::sin(x[1]);
	double cos_xi = std::cos(x[2]), sin_xi = std::sin(x[2]);
	Matrix3d Mtheta, Mphi, Msi, DMtheta, DMphi, DMsi;
	Mtheta.setZero(); Mphi.setZero(); Msi.setZero();
	Mphi(0,0) = 1; Mphi(1,1) = Mphi(2,2) = cos_phi, Mphi(1, 2) = sin_phi, Mphi(2, 1) = -sin_phi;
	Mtheta(1, 1) = 1; Mtheta(0,0) = Mtheta(2,2) = cos_theta, Mtheta(2,0) = sin_theta, Mtheta(0, 2) = -sin_theta;
	Msi(2, 2) = 1; Msi(0, 0) = Msi(1, 1) = cos_xi, Msi(0, 1) = sin_xi, Msi(1, 0) = -sin_xi;
	DMtheta.setZero(); DMphi.setZero(); DMsi.setZero();
	DMphi(1, 1) = DMphi(2, 2) = -sin_phi, DMphi(1, 2) = cos_phi, DMphi(2, 1) = -cos_phi;
	DMtheta(0, 0) = DMtheta(2, 2) = -sin_theta, DMtheta(2, 0) = cos_theta, DMtheta(0, 2) = -cos_theta;
	DMsi(0, 0) = DMsi(1, 1) = -sin_xi, DMsi(0, 1) = cos_xi, DMsi(1,0) = -cos_xi;
	Matrix3d M = Msi * Mphi * Mtheta, dMsi = DMsi * Mphi * Mtheta, dMphi = Msi * DMphi * Mtheta, dMtheta = Msi * Mphi * DMtheta;
	Vector3d DN, diffN_theta, diffN_phi, diffN_si;
	double xy, yz, zx;
	user_pointer_ad* up = (user_pointer_ad*)(user_pointer);
	int nbf = (int)up->bfn.size(); *f = 0.0;
	g[0] = 0; g[1] = 0; g[2] = 0;
	for (int i = 0; i < nbf; ++i)
	{
		OpenVolumeMesh::Geometry::Vec3d& n = up->bfn[i];
		Vector3d N(n[0], n[1], n[2]);
		//Vector3d  N(1, 2, 3);
		DN = M * N;
		diffN_theta = dMtheta * N;
		diffN_phi = dMphi * N;
		diffN_si = dMsi * N;
		xy = DN[0] * DN[1];
		yz = DN[1] * DN[2];
		zx = DN[2] * DN[0];
		*f += (xy * xy + yz * yz + zx * zx) / 2;
		g[0] += xy * (diffN_phi[0] * DN[1] + DN[0] * diffN_phi[1]) +
			yz * (diffN_phi[1] * DN[2] + DN[1] * diffN_phi[2]) +
			zx * (diffN_phi[2] * DN[0] + DN[2] * diffN_phi[0]);
		g[1] += xy * (diffN_theta[0] * DN[1] + DN[0] * diffN_theta[1]) +
			yz * (diffN_theta[1] * DN[2] + DN[1] * diffN_theta[2]) +
			zx * (diffN_theta[2] * DN[0] + DN[2] * diffN_theta[0]);
		g[2] += xy * (diffN_si[0] * DN[1] + DN[0] * diffN_si[1]) +
			yz * (diffN_si[1] * DN[2] + DN[1] * diffN_si[2]) +
			zx * (diffN_si[2] * DN[0] + DN[2] * diffN_si[0]);
		//std::cout << "f: " << *f << " g: " << g[0] << " " << g[1] << " " << g[2] << std::endl;
	}
	
	//double alpha = x[0]; double beta = x[1]; double gamma = x[2];
	//double cos_alpha = std::cos(alpha); double sin_alpha = std::sin(alpha);
	//double cos_beta = std::cos(beta); double sin_beta = std::sin(beta);
	//double cos_gamma = std::cos(gamma); double sin_gamma = std::sin(gamma);
	//double r00 = cos_alpha * cos_gamma - cos_beta * sin_alpha*sin_gamma;
	//double r01 = sin_alpha * cos_gamma + cos_beta * cos_alpha*sin_gamma;
	//double r02 = sin_beta * sin_gamma;
	//double r10 = -cos_alpha * sin_gamma - cos_beta * sin_alpha*cos_gamma;
	//double r11 = -sin_alpha * sin_gamma + cos_beta * cos_alpha*cos_gamma;
	//double r12 = sin_beta * cos_gamma;
	//double r20 = sin_beta * sin_alpha;
	//double r21 = -sin_beta * cos_alpha;
	//double r22 = cos_beta;
	//user_pointer_ad* up = (user_pointer_ad*)(user_pointer);
	//int nbf = up->bfn.size(); *f = 0.0;
	//g[0] = 0; g[1] = 0; g[2] = 0;
	//for (int i = 0; i < nbf; ++i)
	//{
	//	OpenVolumeMesh::Geometry::Vec3d& n = up->bfn[i];
	//	double nx = r00 * n[0] + r01 * n[1] + r02 * n[2];
	//	double ny = r10 * n[0] + r11 * n[1] + r12 * n[2];
	//	double nz = r20 * n[0] + r21 * n[1] + r22 * n[2];
	//	//*f += (nx*ny)*(nx*ny) + (ny*nz)*(ny*nz) + (nz*nx)*(nz*nx);
	//	*f += nx * nx * nx * nx + ny * ny * ny * ny + nz * nz * nz * nz;
	//	double d_nx_alpha = (-sin_alpha * cos_gamma - cos_beta * cos_alpha*sin_gamma) * n[0] + (cos_alpha*cos_gamma - cos_beta * sin_alpha*sin_gamma) * n[1] + 0 * n[2];
	//	double d_ny_alpha = (sin_alpha*sin_gamma - cos_beta * cos_alpha*cos_gamma) * n[0] + (-cos_alpha * sin_gamma - cos_beta * sin_alpha*cos_gamma) * n[1] + 0 * n[2];
	//	double d_nz_alpha = sin_beta * cos_alpha * n[0] + sin_beta * sin_alpha * n[1] + 0 * n[2];
	//	double d_nx_beta = (sin_beta*sin_alpha*sin_gamma) * n[0] + (-sin_beta * cos_alpha*sin_gamma) * n[1] + cos_beta * sin_gamma * n[2];
	//	double d_ny_beta = (sin_beta*sin_alpha*cos_gamma) * n[0] + (-sin_beta * cos_alpha*cos_gamma) * n[1] + cos_beta * cos_gamma * n[2];
	//	double d_nz_beta = cos_beta * sin_alpha *n[0] - cos_beta * cos_alpha * n[1] - sin_beta * n[2];
	//	double d_nx_gamma = (-cos_alpha * sin_gamma - cos_beta * sin_alpha*cos_gamma) * n[0] + (-sin_alpha * sin_gamma + cos_beta * cos_alpha*cos_gamma) * n[1] + sin_beta * cos_gamma * n[2];
	//	double d_ny_gamma = (-cos_alpha * cos_gamma + cos_beta * sin_alpha*sin_gamma) * n[0] + (-sin_alpha * cos_gamma - cos_beta * cos_alpha*sin_gamma) * n[1] + (-sin_beta * sin_gamma) * n[2];
	//	double d_nz_gamma = 0.0;
	//	/*g[0] += 2 * (nx*ny)*(d_nx_alpha*ny + nx * d_ny_alpha) + 2 * (ny*nz)*(d_ny_alpha*nz + ny * d_nz_alpha) + 2 * (nz*nx)*(d_nz_alpha*nx + nz * d_nx_alpha);
	//	g[1] += 2 * (nx*ny)*(d_nx_beta*ny + nx * d_ny_beta) + 2 * (ny*nz)*(d_ny_beta*nz + ny * d_nz_beta) + 2 * (nz*nx)*(d_nz_beta*nx + nz * d_nx_beta);
	//	g[2] += 2 * (nx*ny)*(d_nx_gamma*ny + nx * d_ny_gamma) + 2 * (ny*nz)*(d_ny_gamma*nz + ny * d_nz_gamma) + 2 * (nz*nx)*(d_nz_gamma*nx + nz * d_nx_gamma);*/
	//	g[0] += 4 * nx * nx * nx * (d_nx_alpha) + 4 * ny * ny * ny * (d_ny_alpha) + 4 * nz * nz * nz * (d_nz_alpha);
	//	g[1] += 4 * nx * nx * nx * (d_nx_beta) + 4 * ny * ny * ny * (d_ny_beta) + 4 * nz * nz * nz * (d_nz_beta);
	//	g[2] += 4 * nx * nx * nx * (d_nx_gamma) + 4 * ny * ny * ny * (d_ny_gamma) + 4 * nz * nz * nz * (d_nz_gamma);
	//}
}
void newiteration_ad_LBFGS(int iter, int call_iter, double *x, double* f, double *g, double* gnorm, void* user_pointer)
{
	printf("%d %d; %4.3e, %4.3e\n", iter, call_iter, f[0], gnorm[0]);
}
bool ad_LBFGS(std::vector<OpenVolumeMesh::Geometry::Vec3d>& bfn, std::vector<double>& X)
{
	user_pointer_ad up;
	up.bfn = bfn;
	//LBFGS
	double parameter[20];
	int info[20];
	//initialize
	INIT_HLBFGS(parameter, info);
	parameter[2] = 0.9;
	info[4] = 100; //number of iteration
	info[6] = 0;
	info[7] = 0; //if with hessian 1, without 0
	info[10] = 0;
	int N = 3; int M = 7;
	//m is the number of history value
	//n is the number of variables
	printf("-------------------------------\n");
	printf("start LBFGS\n");
	//function change to nx^4 + ny^4 + nz^ 4?
	//HLBFGS(N, M, &X[0], evalfunc_ad_LBFGS, 0, HLBFGS_UPDATE_Hessian, newiteration_ad_LBFGS, parameter, info, &up);
	HLBFGS(N, M, &X[0], evalfunc_ad_LBFGS_xyz, 0, HLBFGS_UPDATE_Hessian, newiteration_ad_LBFGS, parameter, info, &up);
	return true;
}
