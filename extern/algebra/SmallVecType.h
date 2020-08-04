#pragma once

// Usage: SmallVec type-def.
// Created: 2015.8.17
// Last Modified:2015.8.17
// 2015.8.17: Initial Version.	
// Version: 0.1.150817.2143


#include "SmallVec.h"

namespace ig
{
	/*typedef CVec<uint8_t, 1> byte1;
	typedef CVec<uint8_t, 2> byte2;
	typedef CVec<uint8_t, 3> byte3;
	typedef CVec<uint8_t, 4> byte4;

	typedef CVec<uint32_t, 1> uint1;
	typedef CVec<uint32_t, 2> uint2;
	typedef CVec<uint32_t, 3> uint3;
	typedef CVec<uint32_t, 4> uint4;

	typedef CVec<int, 1> int1;
	typedef CVec<int, 2> int2;
	typedef CVec<int, 3> int3;
	typedef CVec<int, 4> int4;

	typedef CVec<float, 1> float1;
	typedef CVec<float, 2> float2;
	typedef CVec<float, 3> float3;
	typedef CVec<float, 4> float4;

	typedef CVec<double, 1> double1;
	typedef CVec<double, 2> double2;
	typedef CVec<double, 3> double3;
	typedef CVec<double, 4> double4;*/

	typedef CVec<uint8_t, 1> Vec1b;
	typedef CVec<uint8_t, 2> Vec2b;
	typedef CVec<uint8_t, 3> Vec3b;
	typedef CVec<uint8_t, 4> Vec4b;

	typedef CVec<uint32_t, 1> Vec1u;
	typedef CVec<uint32_t, 2> Vec2u;
	typedef CVec<uint32_t, 3> Vec3u;
	typedef CVec<uint32_t, 4> Vec4u;

	typedef CVec<int, 1> Vec1i;
	typedef CVec<int, 2> Vec2i;
	typedef CVec<int, 3> Vec3i;
	typedef CVec<int, 4> Vec4i;

	typedef CVec<float, 1> Vec1f;
	typedef CVec<float, 2> Vec2f;
	typedef CVec<float, 3> Vec3f;
	typedef CVec<float, 4> Vec4f;

	typedef CVec<double, 1> Vec1d;
	typedef CVec<double, 2> Vec2d;
	typedef CVec<double, 3> Vec3d;
	typedef CVec<double, 4> Vec4d;
}