#pragma once

#include <vector>
#include <cstdlib>
#include <ostream>

// \addtogroup MathSuite
//@{
//! Sparse Entry class \ingroup MathSuite
class Sparse_Entry
{
public:

	//! Index ID
	size_t index;

	//! Real value
	double value;

public:
	//! constructor
	Sparse_Entry(size_t ind = 0, double v = 0);

	//! The compare function for sorting
	bool operator<(const Sparse_Entry & m_r) const;
};

//! Symmetric status
enum SPARSE_SYMMETRIC_TYPE
{
	NOSYM, /*!< general case */
	SYM_UPPER, /*!< symmetric (store upper triangular part) */
	SYM_LOWER, /*!< symmetric (store lower triangular part) */
	SYM_BOTH   /*!< symmetric (store both upper and lower triangular part) */
};

//!     Storage
enum SPARSE_STORAGE_TYPE
{
	CCS, /*!< compress column format */
	CRS, /*!< compress row format */
	TRIPLE /*!< row-wise coordinate format */
};

//! Array type
enum ARRAY_TYPE
{
	FORTRAN_TYPE, /*!< the index starts from 1 */
	C_TYPE /*!< the index starts from 0 */
};

typedef std::vector<Sparse_Entry> Sparse_Vector;

//////////////////////////////////////////////////////////////////////////
//! Sparse Matrix Class
class Sparse_Matrix
{
private:

	SPARSE_SYMMETRIC_TYPE symmetric_type;
	SPARSE_STORAGE_TYPE storage_type;

	size_t nrows; //!< number of rows
	size_t ncols; //!< number of columns
	size_t m_num_column_of_RHS; //!< number of rhs
	std::vector< Sparse_Vector* > entryset;
	std::vector< double > right_hand_side, solution;

public:

	//! Sparse matrix constructor
	/*!
	* \param m row dimension
	* \param n column dimension
	* \param symmetrictype
	* \param symmetrictype
	*/
	Sparse_Matrix(const size_t m, const size_t n, const SPARSE_SYMMETRIC_TYPE symmetrictype = NOSYM, const SPARSE_STORAGE_TYPE storagetype = CCS, const size_t num_column_of_RHS = 1);

	Sparse_Matrix(const Sparse_Matrix* mat);
	//! Sparse matrix destructor
	~Sparse_Matrix();

	//! clear
	void reset2zero();

	//! clear_right_hand_side
	void clear_rhs();

	//! adjust rhs
	void adjust_num_of_rhs(const size_t num_column_of_RHS);

	//! fill the right hand side
	void fill_rhs_entry(const size_t row_index, const double val);
	//! set the right hand side
	void set_rhs_entry(const size_t row_index, const double val);
	//! fill the right hand side
	void fill_rhs_entry(const size_t row_index, const size_t which_col, const double val);
	//! set the right hand side
	void set_rhs_entry(const size_t row_index, const size_t which_col, const double val);
	//! fill matrix entry \f$  Mat_{row_index, col_index} += val \f$
	void fill_entry(const size_t row_index, const size_t col_index, const double val);
	//! set matrix entry \f$  Mat_{row_index, col_index} = val \f$
	void set_entry(const size_t row_index, const size_t col_index, const double val);
	//! get matrix entry $Mat_{row_index, col_index}$
	double get_entry(const size_t row_index, const size_t col_index);
	//! get the pointer of matrix entry $Mat_{row_index, col_index}$
	double* get_entry_pointer(const size_t row_index, const size_t col_index);
	//////////////////////////////////////////////////////////////////////////

	//! get the number of nonzeros
	const size_t get_nonzeros() const;
	//! get the row dimension
	const size_t rows() const;
	//! get the column dimension
	const size_t cols() const;
	//! return the number column of the right hand side
	const size_t get_num_rhs() const;
	//! return the symmetric state
	const bool issymmetric() const;

	//! tell whether the matrix is upper or lower symmetric
	const bool issym_store_upper_or_lower() const;

	//! return symmetric state
	const SPARSE_SYMMETRIC_TYPE get_symmetric_type() const;

	//! tell whether the matrix is square
	const bool issquare() const;

	//! return the storage format
	const SPARSE_STORAGE_TYPE get_storage_type() const;
	//! get the entryset
	const std::vector<Sparse_Vector*>& get_entryset() const;
	std::vector<Sparse_Vector*>& get_entryset();
	//! get the rhs
	const std::vector<double>& get_rhs() const;
	std::vector<double>& get_rhs();
	//! get the solution
	const std::vector<double>& get_solution() const;
	std::vector<double>& get_solution();
	//////////////////////////////////////////////////////////////////////////
	void get_compress_data(std::vector<size_t>& rowind, std::vector<size_t>& colptr, std::vector<double>& values, const ARRAY_TYPE array_type = C_TYPE,
		std::vector<double>* diag = 0);
	//////////////////////////////////////////////////////////////////////////
	void get_compress_data(std::vector<int>& rowind, std::vector<int>& colptr, std::vector<double>& values, const ARRAY_TYPE array_type = C_TYPE,
		std::vector<double>* diag = 0);
	//! remove sym info
	void remove_sym_info();
	//////////////////////////////////////////////////////////////////////////
private:
	//! fill_internal_entry
	void fill_internal_entry(const size_t row_index, const size_t col_index, const double val);
	//! set_internal_entry
	void set_internal_entry(const size_t row_index, const size_t col_index, const double val);
	//////////////////////////////////////////////////////////////////////////
};

//! print sparse matrix
std::ostream & operator<<(std::ostream & output_stream, const Sparse_Matrix * mat);
std::ostream & operator<<(std::ostream & output_stream, const Sparse_Matrix & mat);
//@}

//! multiplication \f$ Y = A X \f$
/*!
*	\param A sparse matrix
*	\param X input column vector
*	\param Y output column vector
*/
void multiply(const Sparse_Matrix *A, double *X, double *Y);
//! multiplication \f$ Y = A^T X \f$
/*!
*	\param A sparse matrix
*	\param X input column vector
*	\param Y output column vector
*/
void transpose_multiply(const Sparse_Matrix *A, double *X, double *Y);

//! multiplication \f$ Y = (A^T A) X \f$
/*!
*	\param A sparse matrix
*	\param X input column vector
*	\param Y output column vector
*/
void transpose_self_multiply(const Sparse_Matrix *A, double *X, double *Y, double *tmp_vec);
//! multiplication \f$ Y = (A A^T) X \f$
/*!
*	\param A sparse matrix
*	\param X input column vector
*	\param Y output column vector
*/
void self_transpose_multiply(const Sparse_Matrix *A, double *X, double *Y, double *tmp_vec);
//! convert matrix storage
/*!
*	\param A sparse matrix
*	\param store the specified storage type
*	\param sym how to store the converted symmetric matrix. only valid, when A is symmetric
*/
Sparse_Matrix* convert_sparse_matrix_storage(const Sparse_Matrix *A, SPARSE_STORAGE_TYPE store, SPARSE_SYMMETRIC_TYPE sym);
//! matrix transposition
/*!
*	\param A sparse matrix
*	\param store the specified storage type
*	\param sym how to store the converted symmetric matrix. only valid, when A is symmetric
*/
Sparse_Matrix* transpose(Sparse_Matrix* A, SPARSE_STORAGE_TYPE store, SPARSE_SYMMETRIC_TYPE sym);
//! compute \f$ A^T A \f$
/*!
*	\param A sparse matrix
*	\param positive indicate whether the result is symmetric definite positive.
*	\param store the specified storage type
*	\param sym how to store the converted symmetric matrix. only valid, when A is symmetric
*	\param apply_to_rhs  perform A^T B
*/
Sparse_Matrix* TransposeTimesSelf(Sparse_Matrix* mat, SPARSE_STORAGE_TYPE store = CCS, SPARSE_SYMMETRIC_TYPE sym = SYM_LOWER, bool apply_to_rhs = true);

//! compute \f$ A^T A \f$
Sparse_Matrix* CCS_TransposeTimesSelf(Sparse_Matrix* mat, SPARSE_SYMMETRIC_TYPE sym, bool apply_to_rhs);

//! utils: remove zero elements from a sparse vector
void remove_zero_element(Sparse_Vector& vec, double tiny_value);

//! utils: save sparse pattern to a BMP file
void save_sparse_pattern(const char filename[], const Sparse_Matrix* mat, unsigned int resolution = 400);
