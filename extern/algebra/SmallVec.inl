#pragma once

namespace ig
{
	template<typename T, int N>
	inline CVec<T, N>::CVec()
	{
		//memset(v, 0, sizeof(T)*N);
	}
	// copy constructor & assignment operator
	// use the default one

	// other constructors
	template<typename T, int N>
	inline CVec<T, N>::CVec(const T& rhs)
	{
		for (int i = 0; i < N; i++)
			v[i] = rhs;
	}
	template<typename T, int N>
	inline CVec<T, N>::CVec(const T& x, const T& y)  // 2D only
	{
		static_assert(2 == N, "CVec( const T& x,  const T& y ) can only be applied to 2D vector");
		v[0] = x;
		v[1] = y;
	}
	template<typename T, int N>
	inline CVec<T, N>::CVec(const T& x, const T& y, const T& z)  // 3D only
	{
		static_assert(3 == N, "CVec( const T& x,  const T& y, const T& z ) can only be applied to 3D vector");
		v[0] = x;
		v[1] = y;
		v[2] = z;
	}
	template<typename T, int N>
	inline CVec<T, N>::CVec(const T& x, const T& y, const T& z, const T& w)  // 4D only
	{
		static_assert(4 == N, "CVec( const T& x,  const T& y, const T& z, const T& w ) can only be applied to 4D vector");
		v[0] = x;
		v[1] = y;
		v[2] = z;
		v[3] = w;
	}
	template<typename T, int N>
	inline CVec<T, N>::CVec(const T(&rhs)[N])
	{
		// haven't tested whether memcpy is faster.
		for (int i = 0; i < N; i++)
			v[i] = rhs[i];
	}

	template<typename T, int N>
	void CVec<T, N>::set(const T &scalar)
	{
		for (int i = 0; i < N; i++)
		{
			v[i] = scalar;
		}
	}
	template<typename T, int N>
	int CVec<T, N>::comparearrays(const CVec &rkV) const
	{
		return memcmp(v, rkV.v, N * sizeof(T));
	}

	// accessors
	template<typename T, int N>
	inline T CVec<T, N>::x() const
	{
		return v[0];
	}
	template<typename T, int N>
	inline T CVec<T, N>::y() const
	{
		static_assert(2 <= N, "y() can only be applied to vectors with dim >=2");
		return v[1];
	}
	template<typename T, int N>
	inline T CVec<T, N>::z() const
	{
		static_assert(3 <= N, "z() can only be applied to vectors with dim >=3");
		return v[2];
	}
	template<typename T, int N>
	inline T CVec<T, N>::w() const
	{
		static_assert(4 <= N, "w() can only be applied to vectors with dim >=4");
		return v[3];
	}

	template<typename T, int N>
	inline T& CVec<T, N>::x()
	{
		return v[0];
	}
	template<typename T, int N>
	inline T& CVec<T, N>::y()
	{
		static_assert(2 <= N, "y() can only be applied to vectors with dim >=2");
		return v[1];
	}
	template<typename T, int N>
	inline T& CVec<T, N>::z()
	{
		static_assert(3 <= N, "z() can only be applied to vectors with dim >=3");
		return v[2];
	}
	template<typename T, int N>
	inline T& CVec<T, N>::w()
	{
		static_assert(4 <= N, "w() can only be applied to vectors with dim >=4");
		return v[3];
	}

	template<typename T, int N>
	inline const T& CVec<T, N>::operator[](int32_t idx) const
	{
		return v[idx];
	}
	template<typename T, int N>
	inline T& CVec<T, N>::operator[](int32_t idx)
	{
		return v[idx];
	}

	template<typename T, int N>
	inline const T& CVec<T, N>::operator()(int32_t idx) const
	{
		return v[idx - 1];
	}
	template<typename T, int N>
	inline T& CVec<T, N>::operator()(int32_t idx)
	{
		return v[idx - 1];
	}

	// Unary operators
	template<typename T, int N>
	inline CVec<T, N> CVec<T, N>::operator+() const
	{
		return *this;
	}
	template<typename T, int N>
	inline CVec<T, N> CVec<T, N>::operator-() const
	{
		CVec<T, N> tmp(*this);
		for (int i = 0; i < N; i++)
			tmp.v[i] = -tmp.v[i];
		return tmp;
	}

	// op= Arithmetic operations: only support vectors with same type and dimension
	template<typename T, int N>
	inline const CVec<T, N>& CVec<T, N>::operator+=(const CVec& rhs)
	{
		for (int i = 0; i < N; i++)
			v[i] += rhs.v[i];
		return *this;
	}

	template<typename T, int N>
	inline const CVec<T, N>& CVec<T, N>::operator+=(const T& rhs)
	{
		for (int i = 0; i < N; i++)
			v[i] += rhs;
		return *this;
	}

	template<typename T, int N>
	inline const CVec<T, N>& CVec<T, N>::operator-=(const CVec& rhs)
	{
		for (int i = 0; i < N; i++)
			v[i] -= rhs.v[i];
		return *this;
	}

	template<typename T, int N>
	inline const CVec<T, N>& CVec<T, N>::operator-=(const T& rhs)
	{
		for (int i = 0; i < N; i++)
			v[i] -= rhs;
		return *this;
	}

	template<typename T, int N>
	inline const CVec<T, N>& CVec<T, N>::operator*=(const CVec& rhs)
	{
		for (int i = 0; i < N; i++)
			v[i] *= rhs.v[i];
		return *this;
	}

	template<typename T, int N>
	inline const CVec<T, N>& CVec<T, N>::operator*=(const T& rhs)
	{
		for (int i = 0; i < N; i++)
			v[i] *= rhs;
		return *this;
	}

	template<typename T, int N>
	inline const CVec<T, N>& CVec<T, N>::operator/=(const CVec& rhs)
	{
		for (int i = 0; i < N; i++)
		{
			if (rhs.v[i] != 0)
				v[i] /= rhs.v[i];
		}
		return *this;
	}

	template<typename T, int N>
	inline const CVec<T, N>& CVec<T, N>::operator/=(const T& rhs)
	{
		if (rhs != 0)
		{
			for (int i = 0; i < N; i++)
				v[i] /= rhs;
		}
		return *this;
	}

	// Assign scalar
	template<typename T, int N>
	inline CVec<T, N>& CVec<T, N>::operator=(const T& rhs)
	{
		//static_assert(1 == N, "operator =(scalar) only supports 1D-Vector");
		//v[0] = rhs;
		for (int i = 0; i < N; i++)
			v[i] = rhs;
		return *this;
	}

	// Conversion operators
	template<typename T, int N>
	inline CVec<T, N>::operator T() const
	{
		static_assert(1 == N, "T operator T() can only be applied to 1D vector");
		return v[0];
	}

	// operator==
	//template<typename T, int N>
	//inline bool CVec<T, N>::operator==(const CVec<T, N>& rhs) const
	//{
	//	bool ret = (v[0] == rhs.v[0]);
	//	for (int i = 1; i < N; i++)
	//		ret = ret && (v[i] == rhs.v[i]);
	//	return ret;
	//}
	template<typename T, int N>
	inline bool CVec<T, N>::operator== (const CVec<T, N>& rhs) const
	{
		return comparearrays(rhs) == 0;
	}
	template<typename T, int N>
	inline bool CVec<T, N>::operator!= (const CVec<T, N>& rhs) const
	{
		return  comparearrays(rhs) != 0;
	}
	template<typename T, int N>
	inline bool CVec<T, N>::operator< (const CVec<T, N>& rhs) const
	{
		return comparearrays(rhs) < 0;
	}
	template<typename T, int N>
	inline bool CVec<T, N>::operator<= (const CVec<T, N>& rhs) const
	{
		return  comparearrays(rhs) <= 0;
	}
	template<typename T, int N>
	inline bool CVec<T, N>::operator> (const CVec<T, N>& rhs) const
	{
		return  comparearrays(rhs) > 0;
	}
	template<typename T, int N>
	inline bool CVec<T, N>::operator>= (const CVec<T, N>& rhs) const
	{
		return  comparearrays(rhs) >= 0;
	}

	// Commonly used operations
	template<typename T, int N>
	inline T CVec<T, N>::Dot(const CVec &rhs) const
	{
		T ret = v[0] * rhs.v[0];
		for (int i = 1; i < N; i++)
			ret += v[i] * rhs.v[i];
		return ret;
	}

	template<typename T, int N>
	inline CVec<T, N> CVec<T, N>::Cross(const CVec &TV) const
	{
		static_assert(3 == N, "CVec Cross( const CVec &v ) only supports 3D-Vector");
		return CVec<T, N>(v[1] * TV.v[2] - v[2] * TV.v[1], v[2] * TV.v[0] - v[0] * TV.v[2], v[0] * TV.v[1] - v[1] * TV.v[0]);

		//static_assert(3 <= N, "CVec Cross( const CVec &v ) only supports N(>=3)D-Vector");

		//if (N == 3)
		//{
		//	return CVec<T, N>(v[1] * TV.v[2] - v[2] * TV.v[1], v[2] * TV.v[0] - v[0] * TV.v[2], v[0] * TV.v[1] - v[1] * TV.v[0]);
		//}
		//else
		//{
		//	CVec<T, N> tmp;
		//	for (int i = 0; i < N; i++)
		//	{
		//		int id1 = (i + 1) % N;
		//		int id2 = (i + 2) % N;
		//		tmp[i] = v[id1] * TV.v[id2] - v[id2] * TV.v[id1];
		//	}
		//	return tmp;
		//}
	}

	template<typename T, int N>
	inline CVec<T, N> CVec<T, N>::UnitCross(const CVec &TV) const
	{
		CVec V = Cross(TV);
		V.Normalize();
		return V;
	}

	template<typename T, int N>
	inline T CVec<T, N>::L2NormSqr() const
	{
		return this->Dot(*this);
	}

	template<typename T, int N>
	inline T CVec<T, N>::L2Norm() const
	{
		return std::sqrt(this->Dot(*this));
	}

	template<typename T, int N>
	inline T CVec<T, N>::Length() const
	{
		return std::sqrt(this->Dot(*this));
	}
	template<typename T, int N>
	inline T CVec<T, N>::SquaredLength() const
	{
		return this->Dot(*this);
	}

	template<typename T, int N>
	inline T CVec<T, N>::Normalize()
	{
		T len = L2Norm();

		if (len != (T)0)
		{
			(*this) *= (T)1.0 / len;
		}
		else
		{
			set((T) 0.0);
		}
		return len;
	}

	template<typename T, int N>
	inline T CVec<T, N>::NormalizeTo(CVec &rhs) const
	{
		rhs = *this;
		return rhs.Normalize();
	}

	template<typename T, int N>
	inline T CVec<T, N>::DotPerp(const CVec &rhs) const
	{
		static_assert(N >= 2, "T DotPerp( const CVec &v ) only supports N(>=3)D-Vector");
		return v[0] * rhs.v[1] - v[1] * rhs.v[0];
	}

	// Operators
	template<typename T, int N>
	inline CVec<T, N> operator+(const CVec<T, N>& lhs, const CVec<T, N>& rhs)
	{
		CVec<T, N> tmp(lhs); tmp += rhs; return tmp;
	}
	template<typename T, int N, typename U>
	inline CVec<T, N> operator+(const CVec<T, N>& lhs, const U& rhs)  // scalar
	{
		CVec<T, N> tmp(lhs); tmp += (T)rhs; return tmp;
	}
	template<typename T, int N, typename U>
	inline CVec<T, N> operator+(const U& lhs, const CVec<T, N>& rhs)  // scalar
	{
		CVec<T, N> tmp(rhs); tmp += lhs; return tmp;
	}

	template<typename T, int N>
	inline CVec<T, N> operator*(const CVec<T, N>& lhs, const CVec<T, N>& rhs)
	{
		CVec<T, N> tmp(lhs); tmp *= rhs; return tmp;
	}
	template<typename T, int N, typename U>
	inline CVec<T, N> operator*(const CVec<T, N>& lhs, const U& rhs)  // scalar
	{
		CVec<T, N> tmp(lhs); tmp *= (T)rhs; return tmp;
	}
	template<typename T, int N, typename U>
	inline CVec<T, N> operator*(const U& lhs, const CVec<T, N>& rhs)  // scalar
	{
		CVec<T, N> tmp(rhs); tmp *= lhs; return tmp;
	}

	template<typename T, int N>
	inline CVec<T, N> operator-(const CVec<T, N>& lhs, const CVec<T, N>& rhs)
	{
		CVec<T, N> tmp(lhs); tmp -= rhs; return tmp;
	}
	template<typename T, int N, typename U>
	inline CVec<T, N> operator-(const CVec<T, N>& lhs, const U& rhs)  // scalar
	{
		CVec<T, N> tmp(lhs); tmp -= (T)rhs; return tmp;
	}
	template<typename T, int N, typename U>
	inline CVec<T, N> operator-(const U& lhs, const CVec<T, N>& rhs)  // scalar
	{
		CVec<T, N> tmp((T)lhs); tmp -= rhs; return tmp;
	}

	template<typename T, int N>
	inline CVec<T, N> operator/(const CVec<T, N>& lhs, const CVec<T, N>& rhs)
	{
		CVec<T, N> tmp(lhs); tmp /= rhs; return tmp;
	}
	template<typename T, int N, typename U>
	inline CVec<T, N> operator/(const CVec<T, N>& lhs, const U& rhs)  // scalar
	{
		CVec<T, N> tmp(lhs); tmp /= (T)rhs; return tmp;
	}
	template<typename T, int N, typename U>
	inline CVec<T, N> operator/(const U& lhs, const CVec<T, N>& rhs)  // scalar
	{
		CVec<T, N> tmp((T)lhs); tmp /= rhs; return tmp;
	}

	template <typename T, int N>
	std::ostream &operator << (std::ostream &s, const CVec<T, N> &A)
	{
		for (int i = 0; i < N - 1; i++)
		{
			s << A[i] << " ";
		}

		s << A[N - 1];
		return s;
	}

	template <typename T, int N>
	std::istream &operator >> (std::istream &s, CVec<T, N> &A)
	{
		for (int i = 0; i < N; i++)
		{
			s >> A[i];
		}

		return s;
	}
}