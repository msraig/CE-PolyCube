/*! \file region2d.h
 \brief Two dimensional regions for arrays and matrices (Region2D, const_Region2D, Region2D_iterator).

 Contains definitions of Region2D and const_Region2D classes which allow
 working with a rectangular region of a matrix without creating a new
 matrix, and the definition of the class Region2D_iterator which defines an
 iterator over all elements of a Region2D.

 */

/*
 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

// Modified from region2d.h in:
/*
 *
 * Template Numerical Toolkit (TNT): Linear Algebra Module
 *
 * Mathematical and Computational Sciences Division
 * National Institute of Technology,
 * Gaithersburg, MD USA
 *
 *
 * This software was developed at the National Institute of Standards and
 * Technology (NIST) by employees of the Federal Government in the course
 * of their official duties. Pursuant to title 17 Section 105 of the
 * United States Code, this software is not subject to copyright protection
 * and is in the public domain.  The Template Numerical Toolkit (TNT) is
 * an experimental system.  NIST assumes no responsibility whatsoever for
 * its use by other parties, and makes no guarantees, expressed or implied,
 * about its quality, reliability, or any other characteristic.
 *
 * BETA VERSION INCOMPLETE AND SUBJECT TO CHANGE
 * see http://math.nist.gov/tnt for latest updates.
 *
 */

#ifndef SCPPNT_REGION2D_H
#define SCPPNT_REGION2D_H

#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt.h"
#include "slice_pointer.h"
#include "scppnt_error.h"
#include "matop.h"
#include "index.h"
#else
#include "scppnt/scppnt.h"
#include "scppnt/slice_pointer.h"
#include "scppnt/scppnt_error.h"
#include "scppnt/matop.h"
#include "scppnt/index.h"
#endif

#ifndef SCPPNT_NO_IO
#include <iostream>
#endif

namespace SCPPNT
{

  /*!
   \brief Foward iterator over all elements of a matrix region.

   Foward iterator over all elements of a matrix region. Can be used
   to define an iterator over all elements of a matrix in either
   row-major or column-major order.
   For instance, Region2D_iterator<double, Matrix<double>::columns_iterator>
   will be an iterator over all elements of a matrix in column-major order
   even if the elements of the matrix are stored in row-major order.

   \param T Type of elements iterator points to. This is not the same as
   iterator_traits<iterator_traits<IT2D>::value_type>::value_type
   if IT2D is a constant iterator (e.g., Matrix<T>::const_rows_iterator).

   \param IT2D A two dimensional iterator (iterator over iterators to row elements or
   iterator over iterators to column elements). Possible values of this parameter
   are: Matrix<T>::rows_iterator, Matrix<T>::columns_iterator, Matrix<T>::const_rows_iterator,
   Matrix<T>::const_columns_iterator. This parameter determines whether the order of
   iteration is row-major or column-major.

   The techinique used for logical comparison operators that allow
   constant types to be compared with non-constant types is described
   in:

   Austern, Matt (2001, January). Defining iterators and const iterators.
   C/C++ Users Journal (www.cuj.com), 74-79.

   */
  template<class T, class IT2D> class Region2D_iterator
  {
public:

    // Types needed for iterator_traits
    typedef T value_type; //!< Type of element pointed to
    typedef std::forward_iterator_tag iterator_category; //!< Only allows forward iteration
    typedef Subscript difference_type; //!< Type for number of elements separating two iterators

    typedef IT2D iterator_2D; //! Iterator over row or column iterators

    // "typename" keyword added on next three lines (ww, 12-2-2007)
    typedef typename std::iterator_traits<iterator_2D>::value_type iterator_1D; //! Iterator over row or column elements

    typedef typename std::iterator_traits<iterator_1D>::pointer pointer; //! Pointer to type

    typedef typename std::iterator_traits<iterator_1D>::reference reference; //! Reference to type pointed to

    //! Constructor. Second and third arguments are 1-based indices of initial element iterator points to.
    Region2D_iterator(IT2D it, Subscript n2d, Subscript n1d, Subscript initial2 = 1,
        Subscript initial1 = 1);

    //! Return reference to current element pointed to
    reference operator*() const
    {
      return *current_;
    }

    //! Increment to point to next element and return modified Region2D_iterator (preincrement)
    Region2D_iterator& operator++()
    {
      ++current_;
      if (current_ == end_)
      {
        ++i2D_;
        if (i2D_ != end2D_)
        {
          current_ = *i2D_;
          end_ = current_ + n1D_;
        }
      }
      return *this;
    }

    //! Increment to point to next element and return Region2D_iterator pointing to original element (postincrement)
    Region2D_iterator operator++(int)
    {
      Region2D_iterator t(*this);
      ++current_;
      if (current_ == end_)
      {
        ++i2D_;
        if (i2D_ != end2D_)
        {
          current_ = *i2D_;
          end_ = current_ + n1D_;
        }
      }
      return t;
    }

#ifdef SCPPNT_MEMBER_COMPARISONS

    //! Returns true if two Region2D_iterator objects point to the same element
    bool operator==(const Region2D_iterator<T, IT2D> &rhs)
    {
      return current_ == rhs.current_;
    }

    //! Returns true if two Region2D_iterator objects point to the different elements
    bool operator!=(const Region2D_iterator<T, IT2D> &rhs)
    {
      return current_ != rhs.current_;
    }

#else

    //! Returns true if two Region2D_iterator objects point to the same element
    friend bool operator==(const Region2D_iterator<T, IT2D> &lhs,
        const Region2D_iterator<T, IT2D> &rhs)
    {
      return lhs.current_ == rhs.current_;
    }

    //! Returns true if two Region2D_iterator objects point to the different elements
    friend bool operator!=(const Region2D_iterator<T, IT2D> &lhs,
        const Region2D_iterator<T, IT2D> &rhs)
    {
      return lhs.current_ != rhs.current_;
    }

#endif

private:

    //! Iterator over second-order elements (row iterators or column iterators)
    iterator_2D i2D_;

    //! Points to current element pointed to by current first-order iterator (element of current row or column)
    iterator_1D current_;

    //! Points to one past last element pointed to by current first-order iterator (last element of current row or column)
    iterator_1D end_;

    //! Points to one past last element pointed to by second-order iterator (last row or column iterator)
    iterator_2D end2D_;

    //! Number of elements pointed to by each first-order iterator (number of rows or columns)
    Subscript n1D_;

  };

  /*!
   \brief Constructor

   \param it Iterator over row iterators or column iterators.
   \param n2d Number of valid elements pointed to by it.
   \param n1d Number of valid elements pointed to by each of the iterators it[0], it[1], etc.
   \param initial2 Initial element for 2D iterator (one-offset)
   \param initial1 Initial element  for each 1D iterator (one-offset)
   */
  template<class T, class IT2D> Region2D_iterator<T, IT2D>::Region2D_iterator(IT2D it,
      Subscript n2d, Subscript n1d, Subscript initial2, Subscript initial1) :
    i2D_(it), n1D_(n1d)
  {
    if (initial2 >= n2d || initial1 >= n1d)
      throw BadDimension("SCPPNT::Region2D_iterator::region2D_iterator");
    end2D_ = it + n2d;
    i2D_ += initial2 - 1;
    end_ = *i2D_ + n1D_;
    current_ = *i2D_ + initial1 - 1;
  }

  /*!
   \brief A rectangular subset of a matrix.

   View of a rectangular subset of a matrix. All operations
   that can be done with a SCPPNT::Matrix can be done with
   a Region2D. Note this is a view of the original matrix
   so all operations take place on elements of the original matrix,
   not copies.

   \param Array2D Matrix type that forms the basis of the view.
   Requires interface consistent with SCPPNT::Matrix.

   */
  template<class Array2D> class Region2D
  {
public:

    typedef Array2D array_type; //!< Type of matrix underlying region
    typedef typename Array2D::size_type size_type; //!< Subscript type
    typedef typename Array2D::value_type value_type; //!< Type of elements stored in matrix
    typedef typename Array2D::element_type element_type; //!< Type of elements stored in matrix
    typedef typename Array2D::pointer pointer; //!< Pointer to type stored in matrix
    typedef typename Array2D::reference reference; //!< Reference to type stored in matrix
    typedef typename Array2D::const_reference const_reference; //!< Reference to type stored in constant matrix

    /**** Iterators ****/

    //! Iterator over elements of a row
    typedef typename Array2D::row_iterator row_iterator;

    //! Iterator over constant elements of a row
    typedef typename Array2D::const_row_iterator const_row_iterator;

    //! Iterator over elements of a column
    typedef typename Array2D::column_iterator column_iterator;

    //! Iterator over constant elements of a column
    typedef typename Array2D::const_column_iterator const_column_iterator;

    /*
     rows_iterator - Iterator over row_iterators for rows of matrix.
     columns_iterator - Iterator over column_iterators for columns of matrix.
     */
#ifdef SCPPNT_BOUNDS_CHECK
    typedef slice_pointer_base<row_iterator, row_iterator*, row_iterator&> rows_iterator;
    typedef slice_pointer_base<column_iterator, column_iterator*, column_iterator&>
        columns_iterator;
    typedef
    slice_pointer_base<const_row_iterator, const_row_iterator*, const_row_iterator&>
        const_rows_iterator;
    typedef
    slice_pointer_base<const_column_iterator, const_column_iterator*, const_column_iterator&>
        const_columns_iterator;
#else
    //! Iterator over row iterators (points to row_iterator object for a row)
    typedef row_iterator* rows_iterator;

    //! Iterator over row iterators (points to const_row_iterator object for a row)
    typedef const_row_iterator* const_rows_iterator;

    //! Iterator over column iterators (points to column_iterator object for a column)
    typedef column_iterator* columns_iterator;

    //! Iterator over column iterators (points to const_column_iterator object for a column)
    typedef const_column_iterator* const_columns_iterator;
#endif

    //! Iterator over all elements      ("typename" keyword added on next two lines, ww, 12-6-2007)
    typedef Region2D_iterator<value_type, typename Region2D<Array2D>::rows_iterator> iterator;
    //! Iterator over all constant elements
    typedef Region2D_iterator<value_type, typename Region2D<Array2D>::const_rows_iterator>
        const_iterator;

    //! Iterator over elements of a matrix diagonal
    typedef typename Array2D::diag_iterator diag_iterator;
    //! Iterator over constant elements of a matrix diagonal
    typedef typename Array2D::const_diag_iterator const_diag_iterator;

    //! Return matrix that is the basis of the view
    const array_type & array() const
    {
      return A_;
    }

    //! Return number of rows in matrix
    Subscript num_rows() const
    {
      return dim_[0];
    }

    //! Return number of columns in matrix
    Subscript num_columns() const
    {
      return dim_[1];
    } //Changed "num_cols" to "num_columns". ww, 12-6-2007

    //! Returns lower bound of subscript
    Subscript lbound() const
    {
      return A_.lbound();
    }

    //! Return number of rows if argument is 1, number of columns if argument is 2
    Subscript dim(Subscript i) const
    {
      if (i< 1 || i> 2)throw InvalidArgument("argument must be 1 or 2", "SCPPNT::Region2D::dim");
      return (i==1) ? dim_[0] : dim_[1];
    }

    //! Return total number of elements in matrix
    Subscript size() const
    {
      return dim_[0]*dim_[1];
    }

    /**** constructors ****/

    /*! \brief Constructor from lower and upper limits for rows and columns
     
     \param A Matrix to use as basis of region.
     \param i1 Index in original matrix of first row in region (1-based).
     \param i2 Index in original matrix of last row in region (1-based).
     \param j1 Index in original matrix of first column in region (1-based).
     \param j2 Index in original matrix of last column in region (1-based).
     
     */
    Region2D(Array2D &A, Subscript i1, Subscript i2, Subscript j1, Subscript j2) : A_(A)
    {
      initialize(i1, i2, j1, j2);
    }

    /*! \brief Constructor using row and column index ranges
     
     \param A Matrix to use as basis of region
     \param I Gives first and last rows of original matrix to use in region.
     \param J Gives first and last columns of original matrix to use in region.
     */
    Region2D(Array2D &A, const Index1D &I, const Index1D &J) : A_(A)
    {
      initialize(I.lbound(), I.ubound(), J.lbound(), J.ubound());
    }

    //! Create a region containing entire matrix (convertion of matrix to region)
    Region2D(Array2D &A) : A_(A)
    {
      initialize(1, A.num_rows(), 1, A.num_columns()); // changed num_cols to "num_columns". ww, 12-6-2007
    }

    /*! \brief Constructor from region
     
     \param A Region to use as basis of region.
     \param i1 Index in original region of first row in new region (1-based).
     \param i2 Index in original region of last row in new region (1-based).
     \param j1 Index in original region of first column in new region (1-based).
     \param j2 Index in original region of last column in new region (1-based).
     */
    Region2D(Region2D<Array2D> &A, Subscript i1, Subscript i2, Subscript j1, Subscript j2)
    : A_(A.A_)
    {
      initialize(i1 + A.offset_[0], i2 + A.offset_[0], j1 + A.offset_[1], j2 + A.offset_[1]);
    }

    //! Copy constructor
    Region2D(const Region2D<Array2D> &A) : A_(A.A_)
    {
      initialize(1 + A.offset_[0], A.dim_[0] + A.offset_[0],
          1 + A.offset_[1], A.dim_[1] + A.offset_[1]);
    }

    //! Destructor
    ~Region2D()
    {
      if (row_ != 0) delete [] (row_);
      if (column_ != 0) delete [] (column_); // changed "col_" to "columns_", ww, 12-6-2007.
      if (const_row_ != 0) delete [] (const_row_);
      if (const_column_ != 0) delete [] (const_column_);} // changed "const_col_" to "const_column_". ww, 12-6-2007. 

    /**** Subscripting ****/

    //! Return iterator to elements in row i+1 (0-based subscripting)
    row_iterator operator[](Subscript i)
    {
      return row_[i];
    }

    //! Return iterator to constant elements in row i+1 (0-based subscripting)
    const_row_iterator operator[](Subscript i) const
    {
      return row_[i];
    }

    //! Return element in row i and column j, where i and j are 1-based indices
    reference operator()(Subscript i, Subscript j)
    {
      return row_[i-1][j-1];
    }

    //! Return constant element in row i and column j, where i and j are 1-based indices
    const_reference operator()(Subscript i, Subscript j) const
    {
      return row_[i-1][j-1];
    }

    /*! \brief Return new region that is subregion of this region.
     
     \param i1 Smallest row in current region to be contained in subregion (1-based)
     \param i2 Largest row in current region to be contained in subregion (1-based)
     \param j1 Smallest column in current region to be contained in subregion (1-based)
     \param j2 Largest column in current region to be contained in subregion (1-based)
     */
    Region2D<Array2D> operator()(Subscript i1, Subscript i2,
        Subscript j1, Subscript j2)
    {
      return Region2D<Array2D>(A_,
          i1+offset_[0], offset_[0] + i2,
          j1+offset_[1], offset_[1] + j2);
    }

    /*! \brief Return new region that is subregion of this region.
     
     \param I Smallest and largest rows of current retion to be contained in subregion (1-based)
     \param J Smallest and largest columns of current retion to be contained in subregion (1-based)
     */
    Region2D<Array2D> operator()(const Index1D &I,
        const Index1D &J)
    {
      return Region2D<Array2D>(A_, I.lbound()+offset_[0],
          offset_[0] + I.ubound(), offset_[1]+J.lbound(),
          offset_[1] + J.ubound());
    }

    /**** Iterators ****/

    //! Return iterator pointing to row iterator for first row (iterator over row iterators)
    rows_iterator begin_rows();

    //! Return iterator pointing to row iterator for first row (constant version)
    const_rows_iterator begin_rows() const;

    //! Return iterator pointing to column iterator for first column (iterator over column iterators)
    columns_iterator begin_columns();

    //! Return iterator pointing to column iterator for first column (constant version)
    const_columns_iterator begin_columns() const;

    //! Return iterator pointing to row iterator for one past last row (iterator over row iterators)
    rows_iterator end_rows();

    //! Return iterator pointing to row iterator for one past last row (constant version)
    const_rows_iterator end_rows() const;

    //! Return iterator pointing to column iterator for one past last column (iterator over column iterators)
    columns_iterator end_columns();

    //! Return iterator pointing to column iterator for one past last column (constant version)
    const_columns_iterator end_columns() const;

    // Return row or column iterator to the first or one past the last element in a particular row or column.
    // Index is one-offset (first row or column is 1).

    //! Return iterator pointing to first element in row 'index' (1-offset)
    row_iterator begin_row(Subscript index)
    {
      return row_[index-1];
    }

    //! Return iterator pointing to first element in row 'index' (1-offset)
    const_row_iterator begin_row(Subscript index) const
    {
      return const_row_[index-1];
    }

    //! Return iterator pointing to first element in column 'index' (1-offset)
    column_iterator begin_column(Subscript index)
    {
      return column_[index-1];
    }

    //! Return iterator pointing to first element in column 'index' (1-offset)
    const_column_iterator begin_column(Subscript index) const
    {
      return const_column_[index-1];
    }

    //! Return iterator pointing to one past last element in row 'index' (1-offset)
    row_iterator end_row(Subscript index)
    {
      return row_[index];
    }

    //! Return iterator pointing to one past last element in row 'index' (1-offset)
    const_row_iterator end_row(Subscript index) const
    {
      return const_row_[index];
    }

    //! Return iterator pointing to one past last element in column 'index' (1-offset)
    column_iterator end_column(Subscript index)
    {
      return column_[index];
    }

    //! Return iterator pointing to one past last element in column 'index' (1-offset)
    const_column_iterator end_column(Subscript index) const
    {
      return const_column_[index];
    }

    // Iterators over all elements of matrix.
    // Should only be used when the order in which the iteration takes place
    // over all elements of the matrix does not matter

    //! Return iterator pointing to first element of matrix
    iterator begin()
    {
      return iterator(begin_rows(), num_rows(), num_columns());
    } // Changed "num_cols" to "num_columns". ww, 12-6-2007.

    //! Return iterator pointing to first element of matrix (constant version)
    const_iterator begin() const
    {
      return const_iterator(begin_rows(), num_rows(), num_columns());
    } // Changed "num_cols" to "num_columns". ww, 12-6-2007.

    //! Return iterator pointing to one past last element of matrix
    iterator end()
    {
      return iterator(begin_rows(), num_rows(), num_columns(), dim_[0]+1, 1);
    } // Changed "num_cols" to "num_columns". ww, 12-6-2007.

    //! Return iterator pointing to one past last element of matrix (constant version)
    const_iterator end() const
    {
      return const_iterator(begin_rows(), num_rows(), num_columns(), dim_[0]+1, 1);
    } // Changed "num_cols" to "num_columns". ww, 12-6-2007.


    // Iterators over diagonals.

    //! Returns iterator pointing to first element of a matrix diagonal
    diag_iterator begin_diagonal(Subscript row, Subscript column);
    // Changed "col" to "column". ww, 12-6-2207.

    //! Returns iterator pointing to first element of matrix diagonal
    const_diag_iterator begin_diagonal(Subscript row, Subscript column) const;
    // Changed "col" to "column". ww, 12-6-2207.

    //! Returns iterator pointing to one past the last element of a matrix diagonal
    diag_iterator end_diagonal(Subscript row, Subscript column);
    // Changed "col" to "column". ww, 12-6-2207.

    //! Returns iterator pointing to one past the last element of matrix diagonal
    const_diag_iterator end_diagonal(Subscript row, Subscript column) const;
    // Changed "col" to "column". ww, 12-6-2207.

    /**** Assignment operators ****/

    //! Scalar assignment
    Region2D<Array2D>& operator=(const value_type& scalar)
    {
      iterator it = begin();
      for (Subscript j = size(); j--; ++it)
      *it = scalar;
    }

    //! Assignment operator
    Region2D<Array2D> & operator=(const Region2D<Array2D> &R);

    /* The following three member templates must be defined within
     the class definition for some compilers (such as Visual C++ 6).
     The definitions of other members are given outside the class definition */

    //! Add matrix rhs to matrix
    template<class MAT>
    Region2D<Array2D>& operator+=(const MAT &rhs)
    {
      matadd(*this, rhs);
      return *this;
    }

    //! Subtract matrix rhs from matrix
    template<class MAT>
    Region2D<Array2D>& operator-=(const MAT &rhs)
    {
      matsub(*this, rhs);
      return *this;
    }

    //! Matrix multiply of matrix times rhs (rhs must be square)
    template<class MAT>
    Region2D<Array2D>& operator*=(const MAT &rhs)
    {
      matmult_assign(*this, rhs);
      return *this;
    }

    //! Add scalar to all elements of the matrix
    Region2D<Array2D>& operator+=(const value_type& value)
    { iterator it = begin();
      for (Subscript j = size(); j--; ++it)
      *it += value;
    }

    //! Subtract scalar from all elements of the matrix
    Region2D<Array2D>& operator-=(const value_type& value)
    {
      iterator it = begin();
      for (Subscript j = size(); j--; ++it)
      *it -= value;
    }

    //! Multiply each element of matrix by scalar value
    Region2D<Array2D>& operator*=(const value_type& value)
    { iterator it = begin(); for (Subscript j = size(); j--; ++it) *it *= value;}

    //! Divide each element of matrix by scalar value
    Region2D<Array2D>& operator/=(const value_type& value)
    {
      iterator it = begin();
      for (Subscript j = size(); j--; ++it)
      *it /= value;
    }

  private:

    //! Initialize data members
    void initialize(Subscript rowLow, Subscript rowHi, Subscript colLow, Subscript colHi);

    Array2D& A_; //!< Matrix that is the basis of the view.
    Subscript offset_[2]; //!< Zero-based offsets of first row and column from those in original matrix.
    Subscript dim_[2]; //!< Number of rows and columns in matrix view.

    row_iterator *row_; //!< points to array of pointers to rows of matrix
    const_row_iterator *const_row_;

    // Changed "col" to "column" in next two lines. ww, 12-6-2007
    column_iterator *column_; //!< points to array of pointers to columns of matrix
    const_column_iterator *const_column_;

    // Changed "col" to "column" in next two lines. ww, 12-6-2007
    //! Return the number of elements in the diagonal beginning at (row, column)
    Subscript diagonal_size(Subscript row, Subscript column) const;

  };

  /*!
   Initialize Region2D object. Allocate memory for row and column iterators.

   \param rowLow Lowest row in original matrix that is in region (1-based)
   \param rowHi Highest row in original matrix that is in region (1-based)
   \param colLow Lowest column in original matrix that is in region (1-based)
   \param colHi Highest column in original matrix that is in region (1-based)
   */
  template <class Array2D>
  void Region2D<Array2D>::initialize(Subscript rowLow, Subscript rowHi,
      Subscript colLow, Subscript colHi)
  {
    Subscript lowerBound = A_.lbound();

    if ( rowHi < rowLow || rowLow < lowerBound || rowHi > A_.num_rows() )
    throw RuntimeError("Bad row index range", "SCPPNT::Region2D<Array2D>::initialize()");
    if ( colHi < colLow || colLow < lowerBound || colHi > A_.num_columns() )
    throw RuntimeError("Bad column index range", "SCPPNT::Region2D<Array2D>::initialize()");

    offset_[0] = rowLow-lowerBound;
    offset_[1] = colLow-lowerBound;
    dim_[0] = rowHi-rowLow+1;
    dim_[1] = colHi-colLow+1;

    /* Allocate space for pointers to rows of matrix. */
    row_ = new row_iterator[dim_[0]+1];
    const_row_ = new const_row_iterator[dim_[0]+1];

    /* assign pointers to elements of row_ */
    Subscript i;
    for (i=0; i<dim_[0]; i++)
    {
      row_[i] = A_.begin_row(i+offset_[0]+1) + offset_[1];
#ifdef SCPPNT_BOUNDS_CHECK
      row_[i].set_size(dim_[1]);
#endif
      const_row_[i] = row_[i];
    }

    /* Assign iterator to point to one past last row iterator */
    row_[i] = row_[i-1] + offset_[1];
    const_row_[i] = row_[i];
#ifdef SCPPNT_BOUNDS_CHECK
    row_[i].set_size(1);
#endif

    /* Allocate space for pointers to columns of matrix. */
    column_ = new column_iterator[dim_[1]+1];
    const_column_ = new const_column_iterator[dim_[1]+1];

    /* assign pointers to elements of column_ */
    for (i=0; i<dim_[1]; ++i)
    {
      column_[i] = A_.begin_column(i+offset_[1]+1) + offset_[0];
#ifdef SCPPNT_BOUNDS_CHECK
      column_[i].set_size(dim_[0]);
#endif
      const_column_[i] = column_[i];
    }

    /* Assign iterator to point to one past last column iterator */
    column_[i] = column_[i-1] + offset_[0];
    const_column_[i] = column_[i];
#ifdef SCPPNT_BOUNDS_CHECK
    column_[i].set_size(1);
#endif

  }

  //! Return iterator pointing to row iterator for first row (iterator over row iterators)
  template <class Array2D>
  inline typename Region2D<Array2D>::rows_iterator Region2D<Array2D>::begin_rows()
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return rows_iterator(row_, 1, dim_[0]);
#else
    return row_;
#endif
  }

  //! Return iterator pointing to row iterator for first row (constant version)
  template <class Array2D>
  inline typename Region2D<Array2D>::const_rows_iterator Region2D<Array2D>::begin_rows()
  const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return const_rows_iterator(const_row_, 1, dim_[0]);
#else
    return const_row_;
#endif
  }

  //! Return iterator pointing to column iterator for first column (iterator over column iterators)
  template <class Array2D>
  inline typename Region2D<Array2D>::columns_iterator Region2D<Array2D>::begin_columns()
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return columns_iterator(column_, 1, dim_[1]);
#else
    return column_;
#endif
  }

  //! Return iterator pointing to column iterator for first column (constant version)
  template <class Array2D>
  inline typename Region2D<Array2D>::const_columns_iterator Region2D<Array2D>::begin_columns() const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return const_columns_iterator(const_column_, 1, dim_[1]);
#else
    return const_column_;
#endif
  }

  // Return iterator pointing to row iterator for one past last row (iterator over row iterators)
  template <class Array2D>
  inline typename Region2D<Array2D>::rows_iterator Region2D<Array2D>::end_rows()
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return rows_iterator(row_ + dim_[0], 1, 1);
#else
    return row_ + dim_[0];
#endif
  }

  // Return iterator pointing to row iterator for one past last row (iterator over row iterators)
  template <class Array2D>
  inline typename Region2D<Array2D>::const_rows_iterator Region2D<Array2D>::end_rows() const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return const_rows_iterator(const_row_ + dim_[0], 1, 1);
#else
    return (const_rows_iterator) row_ + dim_[0];
#endif
  }

  // Return iterator pointing to column iterator for one past last column (iterator over column iterators)
  template <class Array2D>
  inline typename Region2D<Array2D>::columns_iterator Region2D<Array2D>::end_columns()
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return columns_iterator(column_ + dim_[1], 1, 1);
#else
    return column_ + dim_[1];
#endif
  }

  // Return iterator pointing to column iterator for one past last column (iterator over column iterators)
  template <class Array2D>
  inline typename Region2D<Array2D>::const_columns_iterator Region2D<Array2D>::end_columns() const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return const_columns_iterator(const_column_ + dim_[1], 1, 1);
#else
    return const_column_ + dim_[1];
#endif
  }

  // Returns number of elements in diagonal with first element (row, column)
  template <class Array2D>
  inline Subscript Region2D<Array2D>::diagonal_size(Subscript row, Subscript column) const
  {
    Subscript ncolumn = dim_[0] - column + 1;
    Subscript nrow = dim_[1] - row + 1;
    return (ncolumn > nrow) ? nrow : ncolumn;

  }

  /*!

   Returns iterator pointing to the first element of the matrix diagonal given by
   (row, column), (row+1, column+1), (row+2, column+2), ...
   For example, if row=1 and column=1 then an iterator pointing to the first element of the
   main diagnonal is returned.

   \param row 1-based row number for initial element of diagonal.
   \param column 1-based column number for initial element of diagonal.
   */

  template <class Array2D>
  inline typename Region2D<Array2D>::diag_iterator Region2D<Array2D>::begin_diagonal(Subscript row, Subscript column)
  {

    diag_iterator it = A_.begin_diagonal(row + offset_[0], column + offset_[1]);

#ifdef SCPPNT_BOUNDS_CHECK
    it.set_size(diagonal_size(row, column));
#endif

    return it;
  }

  /*!

   Returns iterator pointing to the first element of the matrix diagonal given by
   (row, column), (row+1, column+1), (row+2, column+2), ...
   For example, if row=1 and column=1 then an iterator pointing to the first element of the
   main diagnonal is returned.

   \param row 1-based row number for initial element of diagonal.
   \param column 1-based column number for initial element of diagonal.
   */
  template <class Array2D>
  inline typename Region2D<Array2D>::const_diag_iterator Region2D<Array2D>::begin_diagonal(Subscript row, Subscript column) const
  {
    const_diag_iterator it = A_.begin_diagonal(row + offset_[0], column + offset_[1]);

#ifdef SCPPNT_BOUNDS_CHECK
    it.set_size(diagonal_size(row, column));
#endif

    return it;
  }

  /*!

   Returns iterator pointing to the one past last element of the matrix diagonal given by
   (row, column), (row+1, column+1), (row+2, column+2), ...
   For example, for a square matrix with n rows and columns end_diagonal(1, 1)
   returns an iterator pointing to the element (row+n, column+n), which is one past
   the last element of the main diagonal.

   \param row 1-based row number for initial element of diagonal.
   \param column 1-based column number for initial element of diagonal.
   */
  template <class Array2D>
  inline typename Region2D<Array2D>::diag_iterator Region2D<Array2D>::end_diagonal(Subscript row, Subscript column)
  {
    diag_iterator it = begin_diagonal(row, column);
    return it + diagonal_size(row, column);
  }

  /*!

   Returns iterator pointing to the one past last element of the matrix diagonal given by
   (row, column), (row+1, column+1), (row+2, column+2), ...
   For example, for a square matrix with n rows and columns end_diagonal(1, 1)
   returns an iterator pointing to the element (row+n, column+n), which is one past
   the last element of the main diagonal.

   \param row 1-based row number for initial element of diagonal.
   \param column 1-based column number for initial element of diagonal.
   */
  template <class Array2D>
  inline typename Region2D<Array2D>::const_diag_iterator Region2D<Array2D>::end_diagonal(Subscript row, Subscript column) const
  {
    const_diag_iterator it = begin_diagonal(row, column);
    return it + diagonal_size(row, column);
  }

  /*!
   \brief Assignment operator.

   Assign contents of this region to those of
   region R.
   */
  template <class Array2D>
  Region2D<Array2D> & Region2D<Array2D>::operator=(const Region2D<Array2D> &R)
  {
    Subscript M = num_rows();
    Subscript N = num_columns();

    if (M != R.num_rows() || N != R.num_columns())
    {
      throw BadDimension("Region2D<Array2D>::operator=(const Region2D)");
    }

    const_iterator ir = R.begin();
    iterator ithis = begin();
    for (int i = M*N; i--; ++ir, ++ithis)
    {
      *ithis = *ir;
    }
    return *this;
  }

#ifndef SCPPNT_NO_IO
  /*!
   \brief Write a matrix to an output stream.

   Writes the number of rows and columnumns, followed
   by a new line. The following lines contain the
   elements in the rows of the matrix.
   */
  template <class Array2D>
  std::ostream& operator<<(std::ostream &s, const Region2D<Array2D> &A)
  {
    Subscript M=A.num_rows();
    Subscript N=A.num_columns();

    Subscript start = A.lbound();
    Subscript Mend = M + A.lbound() - 1;
    Subscript Nend = N + A.lbound() - 1;

    s << M << "  " << N << std::endl;
    for (Subscript i=start; i<=Mend; i++)
    {
      for (Subscript j=start; j<=Nend; j++)
      {
        s << A(i,j) << " ";
      }
      s << std::endl;
    }
    return s;
  }
#endif

  // *******************[ basic matrix operators ]***************************

  //! Add two matrices and return new matrix giving sum
  template <class Array2D, class M>
  inline M operator+(const Region2D<Array2D> &A,
      const M &B)
  {
    M sum(A.num_rows(), A.num_columns());
    matadd(A, B, sum);
    return sum;
  }

  //! Subtract matrix B from matrix A and return matrix giving difference
  template <class Array2D, class M>
  inline M operator-(const Region2D<Array2D> &A,
      const M &B)
  {
    M diff(A.num_rows(), A.num_columns());
    matsub(A, B, diff);
    return diff;
  }

  //! Multiplication operator returning a new matrix containing the matrix product A*B
  template <class Array2D, class M>
  inline M operator*(const Region2D<Array2D> &A,
      const M &B)
  {
    M c(A.num_rows(), B.num_columns());
    matmult(c, A, B);
    return c;
  }

  //! Return a new vector containing the product of the matrix A and the vector x
  template <class T, class Array2D>
  inline Vector<T> operator*(const Region2D<Array2D> &A, const Vector<T> &x)
  {
    return matrix_times_vector< Region2D<Array2D>, Vector<T> >(A, x);
  }

  /*!
   \brief Constant version of Region2D

   View of a rectangular subset of a constant matrix. All operations
   that can be done with a SCPPNT::Matrix can be done with
   a Region2D. Note this is a view of the original matrix
   so all operations take place using elements of the original matrix,
   not copies.

   \param Array2D Matrix type that forms the basis of the view.
   Requires interface consistent with SCPPNT::Matrix.

   */
  template <class Array2D>
  class const_Region2D
  {
  public:

    typedef const Array2D array_type; //!< Type of matrix underlying region
    typedef typename Array2D::size_type size_type; //!< Subscript type
    typedef typename Array2D::value_type value_type; //!< Type of elements stored in matrix
    typedef typename Array2D::element_type element_type; //!< Type of elements stored in matrix
    typedef typename Array2D::pointer pointer; //!< Pointer to type stored in matrix
    typedef typename Array2D::reference reference; //!< Reference to type stored in matrix
    typedef typename Array2D::const_reference const_reference; //!< Reference to type stored in constant matrix

    /**** Iterators ****/

    //! Iterator over constant elements of a row
    typedef typename Array2D::const_row_iterator const_row_iterator;
    //! Iterator over elements of a column
    typedef typename Array2D::const_column_iterator const_column_iterator;

    /*
     rows_iterator - Iterator over row_iterators for rows of matrix.
     columns_iterator - Iterator over column_iterators for columns of matrix.
     */

#ifdef SCPPNT_BOUNDS_CHECK
    typedef slice_pointer_base<const_row_iterator, const_row_iterator*, const_row_iterator&>
    const_rows_iterator;
    typedef slice_pointer_base<const_column_iterator, const_column_iterator*, const_column_iterator&>
    const_columns_iterator;
#else
    //! Iterator over row iterators (points to const_row_iterator object for a row)
    typedef const_row_iterator* const_rows_iterator;

    //! Iterator over column iterators (points to const_column_iterator object for a column)
    typedef const_column_iterator* const_columns_iterator;
#endif

    //! Iterator over all constant elements     ("typename" keyword added. ww, 12-8-2007)
    typedef Region2D_iterator<value_type, typename const_Region2D<Array2D>::const_rows_iterator>
    const_iterator;

    //! Iterator over constant elements of a matrix diagonal
    typedef typename Array2D::const_diag_iterator const_diag_iterator;

    //! Return matrix that is the basis of the view
    const array_type & array() const
    {
      return A_;
    }

    //! Return number of rows in matrix
    Subscript num_rows() const
    {
      return dim_[0];
    }

    //! Return number of columns in matrix
    Subscript num_columns() const
    {
      return dim_[1];
    }

    //! Returns lower bound of subscript
    Subscript lbound() const
    {
      return A_.lbound();
    }

    //! Return number of rows if argument is 1, number of columns if argument is 2
    Subscript dim(Subscript i) const
    {
      if (i < 1 || i > 2) throw InvalidArgument("argument must be 1 or 2", "SCPPNT::const_Region2D::dim");
      return (i==1) ? dim_[0] : dim_[1];
    }

    //! Return total number of elements in matrix
    Subscript size() const
    {
      return dim_[0]*dim_[1];
    }

    /**** constructors ****/

    /*! \brief Constructor from lower and upper limits for rows and columns
     
     \param A Matrix to use as basis of region.
     \param i1 Index in original matrix of first row in region (1-based).
     \param i2 Index in original matrix of last row in region (1-based).
     \param j1 Index in original matrix of first column in region (1-based).
     \param j2 Index in original matrix of last column in region (1-based).
     
     */
    const_Region2D(const Array2D &A, Subscript i1, Subscript i2, Subscript j1,
        Subscript j2) : A_(A)
    {
      initialize(i1, i2, j1, j2);
    }

    /*! \brief Constructor using row and column index ranges
     
     \param A Matrix to use as basis of region
     \param I Gives first and last rows of original matrix to use in region.
     \param J Gives first and last columns of original matrix to use in region.
     */
    const_Region2D(const Array2D &A, const Index1D &I, const Index1D &J) : A_(A)
    {
      initialize(I.lbound(), I.ubound(), J.lbound(), J.ubound());
    }

    //! Create a region containing entire matrix (convertion of matrix to region)
    const_Region2D(const Array2D &A) : A_(A)
    {
      initialize(1, A.num_rows(), 1, A.num_columns());
    }

    /*! \brief Constructor from region
     
     \param A Region to use as basis of region.
     \param i1 Index in original region of first row in new region (1-based).
     \param i2 Index in original region of last row in new region (1-based).
     \param j1 Index in original region of first column in new region (1-based).
     \param j2 Index in original region of last column in new region (1-based).
     */
    const_Region2D(const const_Region2D<Array2D> &A, Subscript i1, Subscript i2, Subscript j1, Subscript j2)
    : A_(A.A_)
    {
      initialize(i1 + A.offset_[0], i2 + A.offset_[0], j1 + A.offset_[1], j2 + A.offset_[1]);
    }

    //! Copy constructor
    const_Region2D(const const_Region2D<Array2D> &A) : A_(A.A_)
    {
      initialize(1 + A.offset_[0], A.dim_[0] + A.offset_[0], 1 + A.offset_[1], A.dim_[1] + A.offset_[1]);
    }

    //! Destructor
    ~const_Region2D()
    {
      if (const_row_ != 0) delete [] (const_row_);
      if (const_column_ != 0) delete [] (const_column_);
    }

    /**** Subscripting ****/

    //! Return iterator to constant elements in row i+1 (0-based subscripting)
    const_row_iterator operator[](Subscript i) const
    {
      return const_row_[i];
    }

    //! Return constant element in row i and column j, where i and j are 1-based indices
    const_reference operator()(Subscript i, Subscript j) const
    {
      return const_row_[i-1][j-1];
    }

    /*! \brief Return new region that is subregion of this region.
     
     \param i1 Smallest row in current region to be contained in subregion (1-based)
     \param i2 Largest row in current region to be contained in subregion (1-based)
     \param j1 Smallest column in current region to be contained in subregion (1-based)
     \param j2 Largest column in current region to be contained in subregion (1-based)
     */
    const_Region2D<Array2D> operator()(Subscript i1, Subscript i2,
        Subscript j1, Subscript j2)
    {
      return const_Region2D<Array2D>(A_,
          i1+offset_[0], offset_[0] + i2,
          j1+offset_[1], offset_[1] + j2);
    }

    /*! \brief Return new region that is subregion of this region.
     
     \param I Smallest and largest rows of current retion to be contained in subregion (1-based)
     \param I Smallest and largest columns of current retion to be contained in subregion (1-based)
     */
    const_Region2D<Array2D> operator()(const Index1D &I,
        const Index1D &J)
    {
      return const_Region2D<Array2D>(A_, I.lbound()+offset_[0],
          offset_[0] + I.ubound(), offset_[1]+J.lbound(),
          offset_[1] + J.ubound());
    }

    /**** Iterators ****/

    //! Return iterator pointing to row iterator for first row (iterator over row iterators)
    const_rows_iterator begin_rows() const;

    //! Return iterator pointing to column iterator for first column (iterator over column iterators)
    const_columns_iterator begin_columns() const;

    //! Return iterator pointing to row iterator for one past last row (iterator over row iterators)
    const_rows_iterator end_rows() const;

    //! Return iterator pointing to column iterator for one past last column (iterator over column iterators)
    const_columns_iterator end_columns() const;

    // Return row or column iterator for particular row or column.
    // Index is one-offset (first row or column is 1).
    // row(N+1) returns a row_iterator to one past the last row, where there are N rows
    // column(M+1) returns a column_iterator to one past the last column, where there are M columns
    //

    //! Return iterator pointing to first element in row 'index' (1-offset)
    const_row_iterator begin_row(Subscript index) const
    {
      return const_row_[index-1];
    }

    //! Return iterator over constant elements in column 'index' (1-offset)
    const_column_iterator begin_column(Subscript index) const
    {
      return const_column_[index-1];
    }

    //! Return iterator pointing to one past last element in row 'index' (1-offset)
    const_row_iterator end_row(Subscript index) const
    {
      return const_row_[index];
    }

    //! Return iterator pointing to one past last element in column 'index' (1-offset)
    const_column_iterator end_column(Subscript index) const
    {
      return const_column_[index];
    }

    // Iterators over all elements of matrix.
    // Should only be used when the order in which the iteration takes place
    // over all elements of the matrix does not matter

    //! Return iterator pointing to first element of matrix (constant version)
    const_iterator begin() const
    {
      return const_iterator(begin_rows(), num_rows(), num_columns());
    }

    //! Return iterator pointing to one past last element of matrix
    const_iterator end() const
    {
      return const_iterator(begin_rows(), num_rows(), num_columns(), dim_[0]+1, 1);
    }

    // Iterators over diagonals.

    //! Returns iterator pointing to first element of matrix diagonal
    const_diag_iterator begin_diagonal(Subscript row, Subscript column) const;

    //! Returns iterator pointing to one past the last element of matrix diagonal
    const_diag_iterator end_diagonal(Subscript row, Subscript column) const;

  private:

    //! Initialize data members
    void initialize(Subscript rowLow, Subscript rowHi, Subscript colLow, Subscript colHi);

    const Array2D& A_; //!< Matrix that is the basis of the view.
    Subscript offset_[2]; //!< Zero-based offsets of first row and column from those in original matrix.
    Subscript dim_[2]; //!< Number of rows and columns in matrix view.

    const_row_iterator *const_row_; //!< points to array of pointers to rows of matrix

    const_column_iterator *const_column_; //!< points to array of pointers to columns of matrix

    //! Return the number of elements in the diagonal beginning at (row, column)
    Subscript diagonal_size(Subscript row, Subscript column) const;
  };

  /*!
   Initialize const_Region2D object. Allocate memory for row and column iterators.

   \param rowLow Lowest row in original matrix that is in region (1-based)
   \param rowHi Highest row in original matrix that is in region (1-based)
   \param colLow Lowest column in original matrix that \
  s in region (1-based)
   \param colHi Highest column in original matrix that is in region (1-based)
   */
  template <class Array2D>
  void const_Region2D<Array2D>::initialize(Subscript rowLow, Subscript rowHi,
      Subscript colLow, Subscript colHi)
  {
    Subscript lowerBound = A_.lbound();

    if ( rowHi < rowLow || rowLow < lowerBound || rowHi > A_.num_rows())
    throw RuntimeError("Bad row index range", "SCPPNT::const_Region2D<Array2D>::initialize()");
    if ( colHi < colLow || colLow < lowerBound || colHi > A_.num_columns())
    throw RuntimeError("Bad column index range", "SCPPNT::const_Region2D<Array2D>::initialize()");

    offset_[0] = rowLow-lowerBound;
    offset_[1] = colLow-lowerBound;
    dim_[0] = rowHi-rowLow+1;
    dim_[1] = colHi-colLow+1;

    /* Allocate space for pointers to rows of matrix. */
    const_row_ = new const_row_iterator[dim_[0]+1];

    /* assign pointers to elements of row_ */
    Subscript i;
    for (i=0; i<dim_[0]; i++)
    {
      const_row_[i] = A_.row(i+offset_[0]+1) + offset_[1];
#ifdef SCPPNT_BOUNDS_CHECK
      const_row_[i].set_size(dim_[1]);
#endif
    }

    /* Assign iterator to point to one past last row iterator */
    // Note: This code section was in Brad Hanson's last know revision,
    //       but it cannot run since neither row_ nor i are defined here.
    //       Since const_row_[] is already determined in the code above,
    //       the following section is commented out (ww, 12-09-2007)
    // row_[i] = row_[i-1] + offset_[1];
    // const_row_[i] = row_[i];
    // #ifdef SCPPNT_BOUNDS_CHECK
    // row_[i].set_size(1);
    // #endif

    /* Allocate space for pointers to columns of matrix. */
    const_column_ = new const_column_iterator[dim_[1]+1];

    /* assign pointers to elements of column_ */
    for (i=0; i<dim_[1]; ++i)
    {
      const_column_[i] = A_.column(i+offset_[1]+1) + offset_[0];
#ifdef SCPPNT_BOUNDS_CHECK
      const_column_[i].set_size(dim_[0]);
#endif
    }

    /* Assign iterator to point to one past last column iterator */
    // Note: This code section was in Brad Hanson's last known revision,
    //       but it cannot run since neither column_ nor i are defined here.
    //       Since const_row_[] is already determined in the code above,
    //       the following section is commented out (ww, 12-09-2007)
    // column_[i] = column_[i-1] + offset_[0];
    // const_column_[i] = column_[i];
    // #ifdef SCPPNT_BOUNDS_CHECK
    // column_[i].set_size(1);
    // #endif
  }

  //! Return iterator pointing to row iterator for first row (constant version)
  template <class Array2D>
  inline typename const_Region2D<Array2D>::const_rows_iterator const_Region2D<Array2D>::begin_rows() const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return const_rows_iterator(const_row_, 1, dim_[0]);
#else
    return const_row_;
#endif
  }

  //! Return iterator pointing to column iterator for first column (constant version)
  template <class Array2D>
  inline typename const_Region2D<Array2D>::const_columns_iterator const_Region2D<Array2D>::begin_columns() const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return const_columns_iterator(const_column_, 1, dim_[1]);
#else
    return const_column_;
#endif
  }

  // Return iterator pointing to row iterator for one past last row (iterator over row iterators)
  template <class Array2D>
  inline typename const_Region2D<Array2D>::const_rows_iterator const_Region2D<Array2D>::end_rows() const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return const_rows_iterator(const_row_ + dim_[0], 1, 1);
#else
    return (const_rows_iterator) row_ + dim_[0];
#endif
  }

  // Return iterator pointing to column iterator for one past last column (iterator over column iterators)
  template <class Array2D>
  inline typename const_Region2D<Array2D>::const_columns_iterator const_Region2D<Array2D>::end_columns() const
  {
#ifdef SCPPNT_BOUNDS_CHECK
    return const_columns_iterator(const_column_ + dim_[1], 1, 1);
#else
    return const_column_ + dim_[1];
#endif
  }

  // Returns number of elements in diagonal with first element (row, column)
  template <class Array2D>
  inline Subscript const_Region2D<Array2D>::diagonal_size(Subscript row, Subscript column) const
  {
    Subscript ncolumn = dim_[0] - column + 1;
    Subscript nrow = dim_[1] - row + 1;

    return (ncolumn > nrow) ? nrow : ncolumn;
  }

  /*!

   Returns iterator pointing to the first element of the matrix diagonal given by
   (row, column), (row+1, column+1), (row+2, column+2), ...
   For example, if row=1 and column=1 then an iterator pointing to the first element of the
   main diagnonal is returned.

   \param row 1-based row number for initial element of diagonal.
   \param column 1-based column number for initial element of diagonal.
   */
  template <class Array2D>
  inline typename const_Region2D<Array2D>::const_diag_iterator
  const_Region2D<Array2D>::begin_diagonal(Subscript row, Subscript column) const
  {
    const_diag_iterator it = A_.begin_diagonal(row + offset_[0], column + offset_[1]);

#ifdef SCPPNT_BOUNDS_CHECK
    it.set_size(diagonal_size(row, column));
#endif

    return it;
  }

  /*!

   Returns iterator pointing to the one past last element of the matrix diagonal given by
   (row, column), (row+1, column+1), (row+2, column+2), ...
   For example, for a square matrix with n rows and columns end_diagonal(1, 1)
   returns an iterator pointing to the element (row+n, column+n), which is one past
   the last element of the main diagonal.

   \param row 1-based row number for initial element of diagonal.
   \param column 1-based column number for initial element of diagonal.
   */
  template <class Array2D>
  inline typename const_Region2D<Array2D>::const_diag_iterator
  const_Region2D<Array2D>::end_diagonal(Subscript row, Subscript column) const
  {

    const_diag_iterator it = begin_diagonal(row, column);
    return it + diagonal_size(row, column);
  }

#ifndef SCPPNT_NO_IO
  /*!
   \brief Write a matrix to an output stream.

   Writes the number of rows and columns, followed
   by a new line. The following lines contain the
   elements in the rows of the matrix.
   */
  template <class Array2D>
  std::ostream& operator<<(std::ostream &s, const const_Region2D<Array2D> &A)
  {
    Subscript M=A.num_rows();
    Subscript N=A.num_columns();

    Subscript start = A.lbound();
    Subscript Mend = M + A.lbound() - 1;
    Subscript Nend = N + A.lbound() - 1;

    s << M << "  " << N << std::endl;
    for (Subscript i=start; i<=Mend; i++)
    {
      for (Subscript j=start; j<=Nend; j++)
      {
        s << A(i,j) << " ";
      }
      s << std::endl;
    }

    return s;
  }
#endif

  // *******************[ basic matrix operators ]***************************

  //! Add two matrices and return new matrix giving sum
  template <class Array2D, class M>
  inline M operator+(const const_Region2D<Array2D> &A,
      const M &B)
  {
    M sum(A.num_rows(), A.num_columns());
    matadd(A, B, sum);
    return sum;
  }

  //! Subtract matrix B from matrix A and return matrix giving difference
  template <class Array2D, class M>
  inline M operator-(const const_Region2D<Array2D> &A,
      const M &B)
  {
    M diff(A.num_rows(), A.num_columns());
    matsub(A, B, diff);
    return diff;
  }

  //! Multiplication operator returning a new matrix containing the matrix product A*B
  template <class Array2D, class M>
  inline M operator*(const const_Region2D<Array2D> &A,
      const M &B)
  {
    M c(A.num_rows(), B.num_columns());
    matmult(c, A, B);
    return c;
  }

  //! Return a new vector containing the product of the matrix A and the vector x
  template <class T, class Array2D>
  inline Vector<T> operator*(const const_Region2D<Array2D> &A, const Vector<T> &x)
  {
    return matrix_times_vector< const_Region2D<Array2D>, Vector<T> >(A, x);
  }

} // namespace SCPPNT

#endif
// SCPPNT_TRANSV_H
