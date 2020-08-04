#pragma once

// Usage: Small Vector type.
// Created: 2015.7.24
// Last Modified: 2015.9.24
// 2015.7.24: Initial Version.
// 2015.8.4: Added L2NormSqr(), L2Norm(), Normalize().
// 2015.8.17: Added type traits.
// 2015.8.17: Moved typedef to SmallVecType.h
// 2015.8.26: Removed CVec(const T*rhs), and added CVec(const T (&rhs)[N]).
// 2015.9.24: merged with yang's vector class
// Version: 0.1.150826.1441

#include <cstdint>
#include <istream>
#include <ostream>

namespace ig
{
	template<typename T, int N>
	class CVec
	{
	private:
		T v[N]; // data

	public:
		typedef T ChannelType;
		static const int32_t ByteSize = sizeof(T)* N;
		static const int32_t ChannelNum = N;
		static const int32_t ChannelSize = sizeof(T);

	public: // constructors
		CVec();

		// copy constructor & assignment operator
		// use the default one

		// other constructors

		explicit CVec(const T& rhs);

		explicit CVec(const T& x, const T& y);  // 2D only

		explicit CVec(const T& x, const T& y, const T& z);  // 3D only

		explicit CVec(const T& x, const T& y, const T& z, const T& w);  // 4D only

		//explicit CVec(const T* rhs), removed, since CVec(0) is ambiguous.
		explicit CVec(const T(&rhs)[N]); // general

		//////////////////////////////////////////////////////////////////////////
		//access
		inline operator const T *() const
		{
			return v;
		}
		//////////////////////////////////////////////////////////////////////////
		inline operator T *()
		{
			return v;
		}

	public: // accessors
		T x() const;
		T y() const;
		T z() const;
		T w() const;

		T& x();
		T& y();
		T& z();
		T& w();

	public: // operators
		const T& operator[](int32_t idx) const;

		T& operator[](int32_t idx);

		const T& operator()(int32_t idx) const;

		T& operator()(int32_t idx);

		// Unary operators
		CVec operator+() const;

		CVec operator-() const;

		// op= Arithmetic operations: only support vectors with same type and dimension
		const CVec& operator+=(const CVec& rhs);

		const CVec& operator+=(const T& rhs);

		const CVec& operator-=(const CVec& rhs);

		const CVec& operator-=(const T& rhs);

		const CVec& operator*=(const CVec& rhs);

		const CVec& operator*=(const T& rhs);

		const CVec& operator/=(const CVec& rhs);

		const CVec& operator/=(const T& rhs);

		// Assign scalar
		CVec& operator=(const T& rhs);

		// Conversion operators
		operator T() const;

		// == operator
		bool operator==(const CVec& rhs) const;
		bool operator!=(const CVec& rhs) const;
		bool operator<(const CVec& rhs) const;
		bool operator<=(const CVec& rhs) const;
		bool operator>(const CVec& rhs) const;
		bool operator>=(const CVec& rhs) const;

	public:
		// Commonly used operations
		T Dot(const CVec &rhs) const;

		CVec Cross(const CVec &TV) const;

		CVec UnitCross(const CVec &TV) const;

		T L2NormSqr() const;

		T L2Norm() const;

		T Length() const;
		T SquaredLength() const;

		T Normalize();

		T NormalizeTo(CVec &rhs) const;

		T DotPerp(const CVec &rhs) const;

	private:
		int comparearrays(const CVec &rkV) const;

		void set(const T &scalar);
	};

	template<typename T, int N>
	CVec<T, N> operator+(const CVec<T, N>& lhs, const CVec<T, N>& rhs);
	template<typename T, int N, typename U>
	CVec<T, N> operator+(const CVec<T, N>& lhs, const U& rhs);  // scalar
	template<typename T, int N, typename U>
	CVec<T, N> operator+(const U& lhs, const CVec<T, N>& rhs);  // scalar

	template<typename T, int N>
	CVec<T, N> operator-(const CVec<T, N>& lhs, const CVec<T, N>& rhs);
	template<typename T, int N, typename U>
	CVec<T, N> operator-(const CVec<T, N>& lhs, const U& rhs);  // scalar
	template<typename T, int N, typename U>
	CVec<T, N> operator-(const U& lhs, const CVec<T, N>& rhs);  // scalar

	template<typename T, int N>
	CVec<T, N> operator*(const CVec<T, N>& lhs, const CVec<T, N>& rhs);
	template<typename T, int N, typename U>
	CVec<T, N> operator*(const CVec<T, N>& lhs, const U& rhs);  // scalar
	template<typename T, int N, typename U>
	CVec<T, N> operator*(const U& lhs, const CVec<T, N>& rhs);  // scalar

	template<typename T, int N>
	CVec<T, N> operator/(const CVec<T, N>& lhs, const CVec<T, N>& rhs);
	template<typename T, int N, typename U>
	CVec<T, N> operator/(const CVec<T, N>& lhs, const U& rhs);  // scalar
	template<typename T, int N, typename U>
	CVec<T, N> operator/(const U& lhs, const CVec<T, N>& rhs);  // scalar

	template <typename T, int N>
	std::ostream &operator << (std::ostream &s, const CVec<T, N> &A);

	template <typename T, int N>
	std::istream &operator >> (std::istream &s, const CVec<T, N> &A);
}

// type traits for small-vec
namespace ig
{
	template<typename T>
	class TypeTraits;

	template<typename T, int N>
	class TypeTraits<CVec<T, N> >
	{
	public:
		typedef typename CVec<T, N>::ChannelType ChannelType;
		typedef CVec<T, N> Type;
		//static const bool IsAggregate = std::is_arithmetic<T>::value;
		static const size_t ByteSize = CVec<T, N>::ByteSize;
		static const size_t ChannelNum = CVec<T, N>::ChannelNum;
		static const size_t ChannelSize = CVec<T, N>::ChannelSize;
	};
	typedef CVec<double, 2> Vector2d;
	typedef CVec<float, 2> Vector2f;
	typedef CVec<double, 3> Vector3d;
	typedef CVec<float, 3> Vector3f;
	typedef CVec<double, 4> Vector4d;
	typedef CVec<float, 4> Vector4f;
}

#include "SmallVec.inl"
