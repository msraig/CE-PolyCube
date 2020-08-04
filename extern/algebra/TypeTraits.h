#pragma once
#include <cstdint>
// Usage: type traits.
// Created: 2015.8.17
// Last Modified: 2015.8.17
// 2015.8.17: Initial Version.	
// Version: 0.1.150817.1800

namespace ig
{
	// type traits for intrinsic arithmetic types.
	// e.g. float, int, double, byte, etc.
	template<typename T>
	class TypeTraits
	{
	public:
		typedef T ChannelType;
		typedef T Type;
		//static const bool IsAggregate = std::is_arithmetic<T>::value;
		static const int32_t ByteSize = sizeof(T);
		static const int32_t ChannelNum = 1;
		static const int32_t ChannelSize = sizeof(T);
	};
}