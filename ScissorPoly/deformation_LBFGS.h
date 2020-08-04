#ifndef DEFORMATION_LBFGS_H
#define DEFORMATION_LBFGS_H
#include "MeshDefinition.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "Polycube_Deformation.h"
#include "HLBFGS.h"
struct deform_pointer
{
	polycube_deformation_interface* pdi;
};
void evalfunc_de_LBFGS(int N, double* x, double *prev_x, double* f, double* g, void* user_pointer);
void evalfunc_de_LBFGS_all(int N, double* x, double *prev_x, double* f, double* g, void* user_pointer);
void evalfunc_de_LBFGS_all_nf(int N, double* x, double *prev_x, double* f, double* g, void* user_pointer);
void evalfunc_de_LBFGS_polylines(int N, double* x, double *prev_x, double* f, double* g, void* user_pointer);
void newiteration_de_LBFGS(int iter, int call_iter, double *x, double* f, double *g, double* gnorm, void* user_pointer);
bool de_LBFGS(deform_pointer *dp, std::vector<double>& X);
bool de_LBFGS_all(deform_pointer *dp, std::vector<double>& X);
bool de_LBFGS_all_nf(deform_pointer *dp, std::vector<double>& X);
bool de_LBFGS_polylines(deform_pointer *dp, std::vector<double>& X);
#endif