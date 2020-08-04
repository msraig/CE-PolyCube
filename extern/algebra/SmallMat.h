#pragma once

// Usage: Small Matrix template type.
// Created: 2015.8.24
// Last Modified: 2015.8.24
// 2015.8.24: Initial Version.							
// Version: 0.1.150824.1550

#include "SmallVec.h"

namespace ig
{
	template<typename T, int Rows, int Cols>
	class CMat
	{
	private:
		union
		{
			T d[Rows * Cols];
			T m[Rows][Cols];
		};

	public:
		CMat();

		// general
		explicit CMat(T rhs[]);
		
		// only for 2x2
		CMat(T v11, T v12,
			T v21, T v22);

		// only for 3x3
		CMat(T v11, T v12, T v13, 
			T v21, T v22, T v23, 
			T v31, T v32, T v33);

		// only for 4x4
		CMat(T v11, T v12, T v13, T v14,
			T v21, T v22, T v23, T v24,
			T v31, T v32, T v33, T v34,
			T v41, T v42, T v43, T v44);

	public:
		// accessors
		T& operator()(int32_t row, int32_t col);
		
		const T& operator()(int32_t row, int32_t col) const;
		

		T& operator[](int32_t idx);
		

		const T& operator[](int32_t idx) const;
		

		// operators
		CMat& operator+=(const CMat& rhs);
		

		CMat& operator-=(const CMat& rhs);
		

		CMat& operator*=(const CMat& rhs);
		

		CMat& operator*=(const T& rhs);
		

		CMat& operator/=(const T& rhs);
		

		// conversion operator
		operator T*();
		
		operator const T*() const;
		

	public: // Commonly used operations
		void Identity();
		
		T Determinant() const;

		void Inverse();
		void InverseTo(CMat& rhs) const;

		void Transpose();
		void TransposeTo(CMat<T, Cols, Rows>& rhs) const;

	};

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols> operator+(const CMat<T, Rows, Cols>& lhs, const CMat<T, Rows, Cols>& rhs);

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols> operator-(const CMat<T, Rows, Cols>& lhs, const CMat<T, Rows, Cols>& rhs);

	template<typename T, int M, int N, int P>
	inline CMat<T, M, N> operator*(const CMat<T, M, P>& lhs, const CMat<T, P, N>& rhs);

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols> operator*(const CMat<T, Rows, Cols>& lhs, const T& rhs);

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols> operator/(CMat<T, Rows, Cols>& lhs, const T& rhs);

	template<typename T, int Rows, int Cols>
	inline CMat<T, Rows, Cols> operator*(const T& lhs, const CMat<T, Rows, Cols>& rhs);

	template<typename T, int Rows, int Cols>
	inline CVec<T, Cols> operator*(const typename CVec<T, Rows>& lhs, const CMat<T, Rows, Cols>& rhs);

	template<typename T, int Rows, int Cols>
	inline CVec<T, Rows> operator*(const CMat<T, Rows, Cols>& lhs, const CVec<T, Cols>& rhs);
}

#include "SmallMat.inl"
