#ifndef Sparse_Entry_H
#define Sparse_Entry_H

//! Lite Sparse Entry class \ingroup MathSuite
class Sparse_Entry
{
public:

	//! Index ID
	size_t index;

	//! Real value
	double value;

public:
	//! constructor
	Sparse_Entry(size_t ind = 0, double v = 0) :
	  index(ind), value(v)
	  {
	  }

	  //! destructor
	  ~Sparse_Entry()
	  {
	  }

	  //! The compare function for sorting
	  inline bool operator<(const Sparse_Entry & m_r) const
	  {
		  return index < m_r.index;
	  }
};

#endif //Lite_Sparse_Entry_H
