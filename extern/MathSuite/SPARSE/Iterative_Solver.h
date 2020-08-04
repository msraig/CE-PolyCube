#pragma once

#include "Sparse_Matrix.h"

bool solve_by_CG(Sparse_Matrix *m_sparse_matrix, int max_iter_num, int *num_iter = 0, double tolerance = 1.0e-8);
bool solve_by_CG_with_jacobian_pred(Sparse_Matrix *m_sparse_matrix, int max_iter_num, int *num_iter = 0, double tolerance = 1.0e-8);

class ordering_tuple
{
public:
	size_t degree, index;
	double value;
	ordering_tuple(size_t d = 0, size_t ind = 0, double v = 0)
		:degree(d), index(ind), value(v) {}
};
bool ordering_tuple_compare(const ordering_tuple& p1, const ordering_tuple& p2);
Sparse_Matrix* cmkOrdering(Sparse_Matrix *m_sparse_matrix, std::vector<size_t>& iperm);
Sparse_Matrix* amdOrdering(Sparse_Matrix *m_sparse_matrix, std::vector<size_t>& iperm);
Sparse_Matrix* metisOrdering(Sparse_Matrix *m_sparse_matrix, std::vector<size_t>& iperm);

class ordering_pair
{
public:
	double value;
	size_t index;
	ordering_pair(double v = 0, size_t ind = 0) :value(v), index(ind) {}
};
bool first_ordering_sparse_vector(const Sparse_Entry& p1, const Sparse_Entry& p2);
bool second_ordering_sparse_vector(const Sparse_Entry& p1, const Sparse_Entry& p2);

Sparse_Matrix* incomplete_Cholesky_factorization(Sparse_Matrix *m_sparse_matrix); //A = GG^T
Sparse_Matrix* incomplete_Cholesky_factorization_by_filling(Sparse_Matrix *m_sparse_matrix, int fill_in_p = 0); //A = GG^T

void back_substitution_with_incomplete_Cholesky(Sparse_Matrix* chol, double* R, double *X);   //GG^T*X = R

bool solve_by_CG_with_incomplete_Cholesky_pred(Sparse_Matrix *m_sparse_matrix, Sparse_Matrix* chol, int max_iter_num, int *num_iter = 0, double tolerance = 1.0e-8);
bool solve_by_CG_with_incomplete_Cholesky_pred(Sparse_Matrix *m_sparse_matrix, int max_iter_num, int *num_iter = 0, double tolerance = 1.0e-8);
bool solve_by_CG_with_incomplete_Cholesky_pred_ordering(Sparse_Matrix *m_sparse_matrix, int max_iter_num, int *num_iter = 0, double tolerance = 1.0e-8);

//BiCGSTAB
int solve_by_BiCGSTAB(Sparse_Matrix *m_sparse_matrix, int max_iter_num, int *num_iter = 0, double tolerance = 1.0e-8);

//nonnegative optimization, lower_bound can be provided
double solve_by_BBNNLS(Sparse_Matrix *m_sparse_matrix, std::vector<double> &lower_bound, int max_iter_num, int *num_iter = 0);
