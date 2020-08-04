#pragma once

// Implementation of SmallMatT
#include <cassert>
#include <cstring>
namespace ig
{
	namespace impl
	{
		template<typename T>
		T CalculateDeterminant2x2(const CMat<T, 2, 2> &m)
		{
			T det = m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
			return det;
		}

		template<typename T>
		void CalculateInverse2x2(const CMat<T, 2, 2> &m, CMat<T, 2, 2> &res)
		{
			// cannot to itself.
			assert(&m != &res);
			T det = CalculateDeterminant2x2(m);
			res(0, 0) = m(1, 1);
			res(0, 1) = -m(0, 1);
			res(1, 0) = -m(1, 0);
			res(1, 1) = m(0, 0);
			res *= 1 / det;
		}

		template<typename T>
		T CalculateDeterminant3x3(const CMat<T, 3, 3> &m)
		{
			T det = m(0, 0) * m(1, 1) * m(2, 2) + m(1, 0) * m(2, 1) * m(0, 2) + m(2, 0) * m(0, 1) * m(1, 2);
			det -= m(0, 0) * m(2, 1) * m(1, 2) + m(2, 0) * m(1, 1) * m(0, 2) + m(1, 0) * m(0, 1) * m(2, 2);
			return det;
		}

		template<typename T>
		void CalculateInverse3x3(const CMat<T, 3, 3> &m, CMat<T, 3, 3> &res)
		{
			assert(&m != &res);
			T det = CalculateDeterminant3x3(m);

			res(0, 0) = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1);
			res(0, 1) = m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2);
			res(0, 2) = m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1);
			
			res(1, 0) = m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2);
			res(1, 1) = m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0);
			res(1, 2) = m(0, 2) * m(1, 0) - m(0, 0) * m(1, 2);
			
			res(2, 0) = m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0);
			res(2, 1) = m(0, 1) * m(2, 0) - m(0, 0) * m(2, 1);
			res(2, 2) = m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);

			res *= 1 / det;
		}

		template<typename T>
		T CalculateDeterminant4x4(const CMat<T, 4, 4> &m)
		{
			T det = m(0, 0) * m(1, 1) * m(2, 2) * m(3, 3) + m(0, 0) * m(1, 2) * m(2, 3) * m(3, 1) + m(0, 0) * m(1, 3) * m(2, 1) * m(3, 2) +
				m(0, 1) * m(1, 0) * m(2, 3) * m(3, 2) + m(0, 1) * m(1, 2) * m(2, 0) * m(3, 3) + m(0, 1) * m(1, 3) * m(2, 2) * m(3, 0) +
				m(0, 2) * m(1, 0) * m(2, 1) * m(3, 3) + m(0, 2) * m(1, 1) * m(2, 3) * m(3, 0) + m(0, 2) * m(1, 3) * m(2, 0) * m(3, 1) +
				m(0, 3) * m(1, 0) * m(2, 2) * m(3, 1) + m(0, 3) * m(1, 1) * m(2, 0) * m(3, 2) + m(0, 3) * m(1, 2) * m(2, 1) * m(3, 0) -
				m(0, 0) * m(1, 1) * m(2, 3) * m(3, 2) - m(0, 0) * m(1, 2) * m(2, 1) * m(3, 3) - m(0, 0) * m(1, 3) * m(2, 2) * m(3, 1) -
				m(0, 1) * m(1, 0) * m(2, 2) * m(3, 3) - m(0, 1) * m(1, 2) * m(2, 3) * m(3, 0) - m(0, 1) * m(1, 3) * m(2, 0) * m(3, 2) -
				m(0, 2) * m(1, 0) * m(2, 3) * m(3, 1) - m(0, 2) * m(1, 1) * m(2, 0) * m(3, 3) - m(0, 2) * m(1, 3) * m(2, 1) * m(3, 0) -
				m(0, 3) * m(1, 0) * m(2, 1) * m(3, 2) - m(0, 3) * m(1, 1) * m(2, 2) * m(3, 0) - m(0, 3) * m(1, 2) * m(2, 0) * m(3, 1);
			return det;
		}

		template<typename T>
		void CalculateInverse4x4(const CMat<T, 4, 4> &m, CMat<T, 4, 4> &res)
		{
			assert(&m != &res);
			T det = CalculateDeterminant4x4(m);
			
			res(0, 0) = m(1, 1) * m(2, 2) * m(3, 3) + m(1, 2) * m(2, 3) * m(3, 1) + m(1, 3) * m(2, 1) * m(3, 2) - m(1, 1) * m(2, 3) * m(3, 2) - m(1, 2) * m(2, 1) * m(3, 3) - m(1, 3) * m(2, 2) * m(3, 1);
			res(0, 1) = m(0, 1) * m(2, 3) * m(3, 2) + m(0, 2) * m(2, 1) * m(3, 3) + m(0, 3) * m(2, 2) * m(3, 1) - m(0, 1) * m(2, 2) * m(3, 3) - m(0, 2) * m(2, 3) * m(3, 1) - m(0, 3) * m(2, 1) * m(3, 2);
			res(0, 2) = m(0, 1) * m(1, 2) * m(3, 3) + m(0, 2) * m(1, 3) * m(3, 1) + m(0, 3) * m(1, 1) * m(3, 2) - m(0, 1) * m(1, 3) * m(3, 2) - m(0, 2) * m(1, 1) * m(3, 3) - m(0, 3) * m(1, 2) * m(3, 1);
			res(0, 3) = m(0, 1) * m(1, 3) * m(2, 2) + m(0, 2) * m(1, 1) * m(2, 3) + m(0, 3) * m(1, 2) * m(2, 1) - m(0, 1) * m(1, 2) * m(2, 3) - m(0, 2) * m(1, 3) * m(2, 1) - m(0, 3) * m(1, 1) * m(2, 2);

			res(1, 0) = m(1, 0) * m(2, 3) * m(3, 2) + m(1, 2) * m(2, 0) * m(3, 3) + m(1, 3) * m(2, 2) * m(3, 0) - m(1, 0) * m(2, 2) * m(3, 3) - m(1, 2) * m(2, 3) * m(3, 0) - m(1, 3) * m(2, 0) * m(3, 2);
			res(1, 1) = m(0, 0) * m(2, 2) * m(3, 3) + m(0, 2) * m(2, 3) * m(3, 0) + m(0, 3) * m(2, 0) * m(3, 2) - m(0, 0) * m(2, 3) * m(3, 2) - m(0, 2) * m(2, 0) * m(3, 3) - m(0, 3) * m(2, 2) * m(3, 0);
			res(1, 2) = m(0, 0) * m(1, 3) * m(3, 2) + m(0, 2) * m(1, 0) * m(3, 3) + m(0, 3) * m(1, 2) * m(3, 0) - m(0, 0) * m(1, 2) * m(3, 3) - m(0, 2) * m(1, 3) * m(3, 0) - m(0, 3) * m(1, 0) * m(3, 2);
			res(1, 3) = m(0, 0) * m(1, 2) * m(2, 3) + m(0, 2) * m(1, 3) * m(2, 0) + m(0, 3) * m(1, 0) * m(2, 2) - m(0, 0) * m(1, 3) * m(2, 2) - m(0, 2) * m(1, 0) * m(2, 3) - m(0, 3) * m(1, 2) * m(2, 0);

			res(2, 0) = m(1, 0) * m(2, 1) * m(3, 3) + m(1, 1) * m(2, 3) * m(3, 0) + m(1, 3) * m(2, 0) * m(3, 1) - m(1, 0) * m(2, 3) * m(3, 1) - m(1, 1) * m(2, 0) * m(3, 3) - m(1, 3) * m(2, 1) * m(3, 0);
			res(2, 1) = m(0, 0) * m(2, 3) * m(3, 1) + m(0, 1) * m(2, 0) * m(3, 3) + m(0, 3) * m(2, 1) * m(3, 0) - m(0, 0) * m(2, 1) * m(3, 3) - m(0, 1) * m(2, 3) * m(3, 0) - m(0, 3) * m(2, 0) * m(3, 1);
			res(2, 2) = m(0, 0) * m(1, 1) * m(3, 3) + m(0, 1) * m(1, 3) * m(3, 0) + m(0, 3) * m(1, 0) * m(3, 1) - m(0, 0) * m(1, 3) * m(3, 1) - m(0, 1) * m(1, 0) * m(3, 3) - m(0, 3) * m(1, 1) * m(3, 0);
			res(2, 3) = m(0, 0) * m(1, 3) * m(2, 1) + m(0, 1) * m(1, 0) * m(2, 3) + m(0, 3) * m(1, 1) * m(2, 0) - m(0, 0) * m(1, 1) * m(2, 3) - m(0, 1) * m(1, 3) * m(2, 0) - m(0, 3) * m(1, 0) * m(2, 1);

			res(3, 0) = m(1, 0) * m(2, 2) * m(3, 1) + m(1, 1) * m(2, 0) * m(3, 2) + m(1, 2) * m(2, 1) * m(3, 0) - m(1, 0) * m(2, 1) * m(3, 2) - m(1, 1) * m(2, 2) * m(3, 0) - m(1, 2) * m(2, 0) * m(3, 1);
			res(3, 1) = m(0, 0) * m(2, 1) * m(3, 2) + m(0, 1) * m(2, 2) * m(3, 0) + m(0, 2) * m(2, 0) * m(3, 1) - m(0, 0) * m(2, 2) * m(3, 1) - m(0, 1) * m(2, 0) * m(3, 2) - m(0, 2) * m(2, 1) * m(3, 0);
			res(3, 2) = m(0, 0) * m(1, 2) * m(3, 1) + m(0, 1) * m(1, 0) * m(3, 2) + m(0, 2) * m(1, 1) * m(3, 0) - m(0, 0) * m(1, 1) * m(3, 2) - m(0, 1) * m(1, 2) * m(3, 0) - m(0, 2) * m(1, 0) * m(3, 1);
			res(3, 3) = m(0, 0) * m(1, 1) * m(2, 2) + m(0, 1) * m(1, 2) * m(2, 0) + m(0, 2) * m(1, 0) * m(2, 1) - m(0, 0) * m(1, 2) * m(2, 1) - m(0, 1) * m(1, 0) * m(2, 2) - m(0, 2) * m(1, 1) * m(2, 0);

			res *= 1 / det;
		}
	}

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols>::CMat() {}

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols>::CMat(T rhs[])
	{
		memcpy_s(d, Rows * Cols * sizeof(T), rhs, Rows * Cols * sizeof(T));
	}

	// only for 2x2
	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols>::CMat(T v11, T v12,
		T v21, T v22)
	{
		static_assert(Rows == 2 && Cols == 2, "This constructor can only be applied to 2x2 matrix");
		m[0][0] = v11; m[0][1] = v12;
		m[1][0] = v21; m[1][1] = v22;
	}

	// only for 3x3
	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols>::CMat(T v11, T v12, T v13,
		T v21, T v22, T v23,
		T v31, T v32, T v33)
	{
		static_assert(Rows == 3 && Cols == 3, "This constructor can only be applied to 3x3 matrix");
		m[0][0] = v11; m[0][1] = v12; m[0][2] = v13; 
		m[1][0] = v21; m[1][1] = v22; m[1][2] = v23;
		m[2][0] = v31; m[2][1] = v32; m[2][2] = v33; 
	}

	// only for 4x4
	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols>::CMat(T v11, T v12, T v13, T v14,
		T v21, T v22, T v23, T v24,
		T v31, T v32, T v33, T v34,
		T v41, T v42, T v43, T v44)
	{
		static_assert(Rows == 4 && Cols == 4, "This constructor can only be applied to 4x4 matrix");
		m[0][0] = v11; m[0][1] = v12; m[0][2] = v13; m[0][3] = v14;
		m[1][0] = v21; m[1][1] = v22; m[1][2] = v23; m[1][3] = v24;
		m[2][0] = v31; m[2][1] = v32; m[2][2] = v33; m[2][3] = v34;
		m[3][0] = v41; m[3][1] = v42; m[3][2] = v43; m[3][3] = v44;
	}

	template<typename T, int Rows, int Cols>
	inline T& CMat<T, Rows, Cols>::operator()(int32_t row, int32_t col)
	{
		return m[row][col];
	}

	template<typename T, int Rows, int Cols>
	inline const T& CMat<T, Rows, Cols>::operator()(int32_t row, int32_t col) const
	{
		return m[row][col];
	}

	template<typename T, int Rows, int Cols>
	inline T& CMat<T, Rows, Cols>::operator[](int32_t idx)
	{
		return d[idx];
	}

	template<typename T, int Rows, int Cols>
	inline const T& CMat<T, Rows, Cols>::operator[](int32_t idx) const
	{
		return d[idx];
	}

	// operators
	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols>& CMat<T, Rows, Cols>::operator+=(const CMat& rhs)
	{
		for (int32_t i = 0; i < Rows * Cols; ++i)
		{
			d[i] += rhs.d[i];
		}
		return *this;
	}

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols>& CMat<T, Rows, Cols>::operator-=(const CMat& rhs)
	{
		for (int32_t i = 0; i < Rows * Cols; ++i)
		{
			d[i] -= rhs.d[i];
		}
		return *this;
	}

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols>& CMat<T, Rows, Cols>::operator*=(const CMat& rhs)
	{
		static_assert(Rows == Cols, "Operator *=() is only supported on square matrix.");
		CMat lhs(*this);
		for (int32_t r = 0; r<Rows; r++)
		for (int32_t c = 0; c<Cols; c++)
		{
			m[r][c] = 0;
			for (int32_t i = 0; i < Rows; ++i)
			{
				m[r][c] += lhs.m[r][i] * rhs.m[i][c];
			}
		}
		return *this;
	}

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols>& CMat<T, Rows, Cols>::operator*=(const T& rhs)
	{
		for (int32_t i = 0; i < Rows * Cols; ++i)
		{
			d[i] *= rhs;
		}
		return *this;
	}

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols>& CMat<T, Rows, Cols>::operator/=(const T& rhs)
	{
		for (int32_t i = 0; i < Rows * Cols; ++i)
		{
			d[i] /= rhs;
		}
		return *this;
	}

	// conversion operator
	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols>::operator T*()
	{
		return d;
	}

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols>::operator const T*() const
	{
		return d;
	}

	template<typename T, int Rows, int Cols>
	inline void CMat<T, Rows, Cols>::Identity()
	{
		static_assert(Rows == Cols, "Identity() is only supported on square matrix.");
		memset(d, 0, sizeof(T)* Rows * Cols);
		for (int32_t i = 0; i < Rows; ++i)
		{
			m[i][i] = 1;
		}
	}

	template<typename T, int Rows, int Cols>
	inline T CMat<T, Rows, Cols>::Determinant() const
	{
		static_assert(0, "Determinant() is not implemented for this matrix type.");
		return 0;
	}

	template<>
	inline float CMat<float, 2, 2>::Determinant() const
	{
		return impl::CalculateDeterminant2x2(*this);
	}

	// Specialization for 2x2, 3x3, 4x4 float & double
	template<>
	inline double CMat<double, 2, 2>::Determinant() const
	{
		return impl::CalculateDeterminant2x2(*this);
	}

	template<>
	inline float CMat<float, 3, 3>::Determinant() const
	{
		return impl::CalculateDeterminant3x3(*this);
	}
	template<>
	inline double CMat<double, 3, 3>::Determinant() const
	{
		return impl::CalculateDeterminant3x3(*this);
	}

	template<>
	inline float CMat<float, 4, 4>::Determinant() const
	{
		return impl::CalculateDeterminant4x4(*this);
	}
	template<>
	inline double CMat<double, 4, 4>::Determinant() const
	{
		return impl::CalculateDeterminant4x4(*this);
	}

	template<typename T, int Rows, int Cols>
	inline void CMat<T, Rows, Cols>::Inverse()
	{
		CMat<T, Rows, Cols> tmp;
		InverseTo(tmp);
		*this = tmp;
	}

	template<typename T, int Rows, int Cols>
	inline void CMat<T, Rows, Cols>::InverseTo(CMat &rhs) const
	{
		static_assert(0, "Inverse() is not implemented for this matrix type.");
	}

	// Specialization for 2x2, 3x3, 4x4 float & double
	template<>
	inline void CMat<float, 2, 2>::InverseTo(CMat<float, 2, 2> &rhs) const
	{
		impl::CalculateInverse2x2(*this, rhs);
	}

	template<>
	inline void CMat<double, 2, 2>::InverseTo(CMat<double, 2, 2> &rhs) const
	{
		impl::CalculateInverse2x2(*this, rhs);
	}
	
	template<>
	inline void CMat<float, 3, 3>::InverseTo(CMat<float, 3, 3> &rhs) const
	{
		impl::CalculateInverse3x3(*this, rhs);
	}
	
	template<>
	inline void CMat<double, 3, 3>::InverseTo(CMat<double, 3, 3> &rhs) const
	{
		impl::CalculateInverse3x3(*this, rhs);
	}

	template<>
	inline void CMat<float, 4, 4>::InverseTo(CMat<float, 4, 4> &rhs) const
	{
		impl::CalculateInverse4x4(*this, rhs);
	}
	
	template<>
	inline void CMat<double, 4, 4>::InverseTo(CMat<double, 4, 4> &rhs) const
	{
		impl::CalculateInverse4x4(*this, rhs);
	}

	template<typename T, int Rows, int Cols>
	inline void CMat<T, Rows, Cols>::TransposeTo(CMat<T, Cols, Rows>& rhs) const
	{
		// cannot to itself.
		assert(&rhs != this);
		for (int32_t r = 0; r<Rows; r++)
		for (int32_t c = 0; c < Cols; c++)
		{
			rhs.m[r][c] = m[c][r];
		}
	}

	template<typename T, int Rows, int Cols>
	inline void CMat<T, Rows, Cols>::Transpose()
	{
		static_assert(Rows == Cols, "Transpose() can only be applied to square matrix.");
		CMat<T, Rows, Cols> tmp;
		TransposeTo(tmp);
		*this = tmp;
	}


	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols> operator+(const CMat<T, Rows, Cols>& lhs, const CMat<T, Rows, Cols>& rhs)
	{
		CMat<T, Rows, Cols> ret(lhs);
		ret += rhs;
		return ret;
	}

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols> operator-(const CMat<T, Rows, Cols>& lhs, const CMat<T, Rows, Cols>& rhs)
	{
		CMat<T, Rows, Cols> ret(lhs);
		ret -= rhs;
		return ret;
	}

	template<typename T, int M, int N, int P>
	inline CMat<T, M, N> operator*(const CMat<T, M, P>& lhs, const CMat<T, P, N>& rhs)
	{
		CMat<T, M, N> ret;
		for (int32_t r = 0; r<M; r++)
		for (int32_t c = 0; c<N; c++)
		{
			ret(r, c) = 0;
			for (int32_t i = 0; i < P; ++i)
			{
				ret(r, c) += lhs(r, i) * rhs(i, c);
			}
		}
		return ret;
	}

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols> operator*(const CMat<T, Rows, Cols>& lhs, const T& rhs)
	{
		CMat<T, Rows, Cols> ret(lhs);
		ret *= rhs;
		return ret;
	}

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols> operator/(CMat<T, Rows, Cols>& lhs, const T& rhs)
	{
		CMat<T, Rows, Cols> ret(lhs);
		ret /= rhs;
		return ret;
	}

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols> operator*(const T& lhs, const CMat<T, Rows, Cols>& rhs)
	{
		CMat<T, Rows, Cols> ret(rhs);
		ret *= lhs;
		return ret;
	}

	template<typename T, int Rows, int Cols>
	inline CVec<T, Cols> operator*(const typename CVec<T, Rows>& lhs, const CMat<T, Rows, Cols>& rhs)
	{
		CVec<T, Cols> ret;
		for (int32_t i = 0; i < Cols; ++i)
		{
			ret[i] = 0;
			for (int32_t j = 0; j < Rows; ++j)
			{
				ret[i] += lhs[j] * rhs(j, i);
			}
		}
		return ret;
	}

	template<typename T, int Rows, int Cols>
	inline CVec<T, Rows> operator*(const CMat<T, Rows, Cols>& lhs, const CVec<T, Cols>& rhs)
	{
		CVec<T, Rows> ret;
		for (int32_t i = 0; i < Rows; ++i)
		{
			ret[i] = 0;
			for (int32_t j = 0; j < Cols; ++j)
			{
				ret[i] += lhs(i, j) * rhs[j];
			}
		}
		return ret;
	}
}

