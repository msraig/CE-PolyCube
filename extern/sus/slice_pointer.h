/*!	\file slice_pointer.h
 \brief Definition of the generalized pointer class slice_pointer_base.
 
 Contains definition of the slice_pointer_base class that is
 a generalization of a normal pointer allowing
 consecutive elements to be stored non-consecutively in memory
 (consecutive elements are separated from one another by a constant distance).

 */

/*
 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

#ifndef SCPPNT_SLICE_POINTER_H
#define SCPPNT_SLICE_POINTER_H

#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt.h"
#include "scppnt_error.h"
#else
#include "scppnt/scppnt.h"
#include "scppnt/scppnt_error.h"
#endif

#include <iterator> // for struct random_access_iterator_category
namespace SCPPNT
{

  /*! 
   \brief Acts like T* but allows consecutive elements to be stored non-consecutively in memory.
   
   This class is a generalization of a normal pointer that allows 
   consecutive elements to be stored non-consecutively in memory
   (consecutive elements are separated from one another by a constant distance).
   
   \param T Type of elements pointed to.
   \param PTR Type of pointer to element (T* or const T*)
   \param REF Type of reference to element (T& or const T&)
   
   The PTR and REF template parameters allow both constant and non-constant types 
   to be defined using one class. This technique is described in
   
   Austern, Matt (1999). Generic programming and the STL.
   Reading, MA: Addison-Wesley.
   
   A slice_pointer_base<T, T*, T&> with a distance of 1 separating consecutive elements
   should have the same behavior as type T*, and a slice_pointer_base<T, const T*, const T&> 
   with a distance of 1 separating consecutive elements should have the same behavior
   as type const T*.
   
   An example of the use of a slice_pointer_base is to allow iterators over column
   elements in row-major matrices stored in one contiguous block of memory,
   where the matrix has M rows and two consecutive elements in a column
   are separated by M-1 elements.
   
   This class does not manage the memory it points to. In other words,
   slice_pointer_base objects are initialized with pointers
   that point to allocated areas of memory, and the memory pointed to is not
   released when slice_pointer_base objects are destroyed.
   
   The techinique used for logical comparison operators that allow
   constant types to be compared with non-constant types is described
   in:
   
   Austern, Matt (2001, January). Defining iterators and const iterators. 
   C/C++ Users Journal (www.cuj.com), 74-79.
   
   */
  template<class T, class PTR, class REF> class slice_pointer_base
  {

public:

    // Types needed for iterator_traits
    typedef T value_type; //!< Type pointed to
    typedef PTR pointer; //!< Pointer to type
    typedef REF reference; //!< Reference to type pointed to
    typedef std::random_access_iterator_tag iterator_category; //!< This is a random access iterator
    typedef Subscript difference_type; //!< Type for number of elements separating two iterators

    //! default constructor
    slice_pointer_base() :
      current(0), span(0)
    {
    }

    /*! \brief Constructor that initializes the slice_pointer_base
     
     \param first Pointer to first element.
     
     \param span Distance between consecutive elements (i-th element is given
     by first + (i-1) * span. If span == 1 then the object behaves the same as a T*.
     
     \param size Number of elements pointed to. 
     This argument is only used when bounds checking is enabled.
     */
    slice_pointer_base(pointer first, Subscript span, Subscript size)
    {
      set(first, span, size);
    }

    //! Copy constructor for non-const type, and construction of const type from non-const type.
#ifdef SCPPNT_BOUNDS_CHECK
    slice_pointer_base(const slice_pointer_base<T, T*, T&> &p)
    {
      set(p.get_first(), p.get_span(), p.get_size());
      current = p.get_current();
    }
#else
    slice_pointer_base(const slice_pointer_base<T, T*, T&> &p)
    {
      set(p.get_current(), p.get_span(), 0);
    }
#endif

    /*! \brief Set values of slice_pointer_base members.
     
     \param first Pointer to first element.
     
     \param span Distance between consecutive elements (i-th element is given
     by first + (i-1) * span. If span == 1 then behavior is the same as T*.
     
     \param size Number of elements pointed to. 
     This argument is only used when bounds checking is enabled.
     */
#ifdef SCPPNT_BOUNDS_CHECK
    void set(pointer first, Subscript span, Subscript size)
#else
    void set(pointer first, Subscript span, Subscript)
    // Prevents warning message about unused argument when size not used
#endif
    {
      current = first;
      this->span = span;
#ifdef SCPPNT_BOUNDS_CHECK
      this->first = first;
      last = first + (size-1) * span;
#endif
    }

    // Access functions.

    //! Return pointer to current element
    pointer get_current() const
    {
      return current;
    }

    //! Return distance separating consecutive elements
    Subscript get_span() const
    {
      return span;
    }

#ifdef SCPPNT_BOUNDS_CHECK
    //! Set number of elements pointed from current element through last element
    void set_size(Subscript size)
    {
      last = current + (size-1) * span;
      first = current;
    }

    //! Return number of elements pointed to
    Subscript get_size() const
    {
      return last-first+1;
    }

    //! Return pointer to first element.
    pointer get_first() const
    {
      return first;
    }
#endif

    //! Return current element pointed to.
    reference operator*() const
    {
#ifdef SCPPNT_BOUNDS_CHECK
      if (!current || current < first || current > last)
throw        BoundsError("SCPPNT::slice_pointer_base::operator*");
#endif
        return *current;
      }

      //! Dereference
      pointer operator->() const
      {
#ifdef SCPPNT_BOUNDS_CHECK
        if (!current || current < first || current > last)
        throw BoundsError("SCPPNT::const_slice_pointer_base::operator->");
#endif
        return current;
      }

      //! zero-offset subscripting
      reference operator[](Subscript i) const
      {
#ifdef SCPPNT_BOUNDS_CHECK
        if (!current || (current+i) < first || (current+i) > last)
        throw BoundsError("SCPPNT::slice_pointer_base::operator[]");
#endif
        return current[i];
      }

      //! Increment to point to next element and return modified slice_pointer_base (preincrement)
      slice_pointer_base& operator++()
      {
        current += span;
        return *this;
      }

      //! Increment to point to next element and return new slice_pointer_base pointing to original element (postincrement)
      slice_pointer_base operator++(int)
      {
        slice_pointer_base t(*this);
        current += span;
        return t;
      }

      //! Decrement to point to previous element and return modified slice_pointer_base (predecrement)
      slice_pointer_base& operator--()
      {
        current -= span;
        return *this;
      }

      //! Decrement to point to previous element and return new slice_pointer_base pointing to original element (postdecrement)
      slice_pointer_base operator--(int)
      {
        slice_pointer_base t(*this);
        current -= span;
        return t;
      }

      //! Increment to point to element that is i elements greater than current element.
      slice_pointer_base& operator+=(Subscript i)
      {
        current += span*i;
        return *this;
      }

      //! Decrement to point to element that is i elements less than current element.
      slice_pointer_base& operator-=(Subscript i)
      {
        current -= span*i;
        return *this;
      }

      //! Return new slice_pointer_base that points to i elements greater than the current element
      slice_pointer_base operator+(Subscript i) const
      {
        return slice_pointer_base<T,PTR,REF>(*this) += i;
      }

      //! Return new slice_pointer_base that points to i elements less than the current element
      slice_pointer_base operator-(Subscript i) const
      {
        return slice_pointer_base<T,PTR,REF>(*this) -= i;
      }

#ifdef SCPPNT_MEMBER_COMPARISONS

      //! Return number of elements separating current elements pointed to by lhs and rhs
      Subscript operator-(const slice_pointer_base<T,PTR,REF> &rhs)
      {
        if (get_span() != rhs.get_span())
        throw InvalidArgument("Mismatching spans", "SCPPNT::slice_pointer_base::operator-()");
        return (get_current() - rhs.get_current()) / lhs.get_span();
      }

      //! Returns true if two slice_pointer_base objects point to the same element
      bool operator==(const slice_pointer_base<T,PTR,REF> &rhs)
      {
        return current == rhs.current;
      }

      //! Returns true if two slice_pointer_base objects point to the different elements
      bool operator!=(const slice_pointer_base<T,PTR,REF> &rhs)
      {
        return current != rhs.current;
      }

      //! Returns true if lhs points to an element greater than the element pointed to by rhs
      bool operator>(const slice_pointer_base<T,PTR,REF> &rhs)
      {
        return current > rhs.current;
      }

      //! Returns true if lhs points to an element greater than or equal to the element pointed to by rhs
      bool operator>=(const slice_pointer_base<T,PTR,REF> &rhs)
      {
        return current >= rhs.current;
      }

      //! Returns true if lhs points to an element less than the element pointed to by rhs
      bool operator<(const slice_pointer_base<T,PTR,REF> &rhs)
      {
        return current < rhs.current;
      }

      //! Returns true if lhs points to an element less than or equal to the element pointed to by rhs
      bool operator<=(const slice_pointer_base<T,PTR,REF> &rhs)
      {
        return current <= rhs.current;
      }

#else

      //! Return number of elements separating current elements pointed to by lhs and rhs
      friend Subscript operator-(const slice_pointer_base<T,PTR,REF> &lhs, const slice_pointer_base<T,PTR,REF> &rhs)
      {
        if (lhs.get_span() != rhs.get_span())
        throw InvalidArgument("Mismatching spans", "SCPPNT::slice_pointer_base::operator-()");
        return (lhs.get_current() - rhs.get_current()) / lhs.get_span();
      }

      //! Returns true if two slice_pointer_base objects point to the same element
      friend bool operator==(const slice_pointer_base<T,PTR,REF> &lhs, const slice_pointer_base<T,PTR,REF> &rhs)
      {
        return lhs.current == rhs.current;
      }

      //! Returns true if two slice_pointer_base objects point to the different elements
      friend bool operator!=(const slice_pointer_base<T,PTR,REF> &lhs, const slice_pointer_base<T,PTR,REF> &rhs)
      {
        return lhs.current != rhs.current;
      }

      //! Returns true if lhs points to an element greater than the element pointed to by rhs
      friend bool operator>(const slice_pointer_base<T,PTR,REF> &lhs, const slice_pointer_base<T,PTR,REF> &rhs)
      {
        return lhs.current > rhs.current;
      }

      //! Returns true if lhs points to an element greater than or equal to the element pointed to by rhs
      friend bool operator>=(const slice_pointer_base<T,PTR,REF> &lhs, const slice_pointer_base<T,PTR,REF> &rhs)
      {
        return lhs.current >= rhs.current;
      }

      //! Returns true if lhs points to an element less than the element pointed to by rhs
      friend bool operator<(const slice_pointer_base<T,PTR,REF> &lhs, const slice_pointer_base<T,PTR,REF> &rhs)
      {
        return lhs.current < rhs.current;
      }

      //! Returns true if lhs points to an element less than or equal to the element pointed to by rhs
      friend bool operator<=(const slice_pointer_base<T,PTR,REF> &lhs, const slice_pointer_base<T,PTR,REF> &rhs)
      {
        return lhs.current <= rhs.current;
      }

#endif

    private:

      //! Pointer to current element.
      pointer current;

      //! Distance separating consecutive elements. For a normal pointer the distance 1.
      Subscript span;

#ifdef SCPPNT_BOUNDS_CHECK
      pointer first; //!< pointer to first valid element
      pointer last; //!< pointer to last valid element
#endif

    };

  } // namespace SCPPNT

#endif // SCPPNT_SLICE_POINTER_H
