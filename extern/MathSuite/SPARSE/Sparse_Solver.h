#pragma once

#include "Sparse_Config.h"
#include "Sparse_Matrix.h"

/*! \addtogroup MathSuite
//@{
//! Solve a linear system by TAUCS library
/*!
*	\param m_sparse_matrix A CCS sparse matrix.
*	\return the solver status. If the matrix is not symmetric but not SPD, the solver will fail.
*/
#ifdef USE_TAUCS
bool solve_by_TAUCS(Sparse_Matrix *m_sparse_matrix, bool spd = true, bool cg_solver = false);
#endif

//! Solve a linear system by UMFPACK library
/*!
*	\param m_sparse_matrix a CCS non-symmetric sparse matrix
*	\return the solver status.
*/
#ifdef USE_SUITESPARSE
bool solve_by_UMFPACK(Sparse_Matrix *m_sparse_matrix);
bool inverse_power_method_by_UMFPACK(Sparse_Matrix* m_sparse_matrix, double* target_eigen_vec,
	int max_iter = 0, double tiny_value = 1.0e-9, double *min_eigenvalue = 0);
#endif

#ifdef USE_SUITESPARSE
bool solve_by_CHOLMOD(Sparse_Matrix *m_sparse_matrix);
#endif

#ifdef USE_ARPACK
//! Eigenvalue type for ARPACK
enum EIGENVALUE_TYPE
{
	ARPACK_LA, //!< largest algebraic value
	ARPACK_SA, //!< smallest algebraic value
	ARPACK_LM, //!< largest magnitude
	ARPACK_SM, //!< smallest magnitude
	ARPACK_LR, //!< largest real part
	ARPACK_SR, //!< smallest real part
	ARPACK_LI, //!< largest imaginary part
	ARPACK_SI, //!< smallest imaginary part
	ARPACK_BE //!< from both ends of spectrum
};

//! character array for ARPACK
static char arpack_type[][3] =
{
	"LA", "SA", "LM", "SM", "LR", "SR", "LI", "SI", "BE"
};

//! Compute SVD by ARPACK library
/*!
*	\param A A sparse matrix.
*	\param num_of_eigenvalue The number of required singular values.
*	\param eigentype The order type for eigenvalues, \sa EIGENVALUE_TYPE .
*	\param eigenvalues The array for singular values
*	\param eigenvectors The array for the right side of SVD
*/
bool solve_SVD_by_ARPACK(Sparse_Matrix *A,
	int num_of_eigenvalue,
	int eigentype,
	std::vector<double>& eigenvalues,
	std::vector<double>& eigenvectors);

//! Compute the eigensystem for a symmetric sparse matrix by ARPACK library
/*!
*	\param A A sparse matrix.
*	\param num_of_eigenvalue The number of required eigenvalues.
*	\param eigentype The order type for eigenvalues, \sa EIGENVALUE_TYPE .
*	\param eigenvalues The array for eigenvalues
*	\param eigenvectors The array for eigenvectors

*/
bool solve_sym_eigensystem_by_ARPACK(Sparse_Matrix *A,
	int num_of_eigenvalue,
	int eigentype,
	std::vector<double>& eigenvalues,
	std::vector<double>& eigenvectors);

#endif

#ifdef USE_SUPERLU
bool solve_by_SUPERLU(Sparse_Matrix *m_sparse_matrix);
#endif

//@}
