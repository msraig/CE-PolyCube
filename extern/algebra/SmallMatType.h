#pragma once

// Usage: SmallMat type-def.
// Created: 2015.8.24
// Last Modified:2015.8.24
// 2015.8.24: Initial Version.	
// Version: 0.1.150824.1550


#include "SmallMat.h"

namespace ig
{
	typedef CMat<float, 2, 2> Mat2x2f;
	typedef CMat<double, 2, 2> Mat2x2d;

	typedef CMat<float, 3, 3> Mat3x3f;
	typedef CMat<double, 3, 3> Mat3x3d;

	typedef CMat<float, 4, 4> Mat4x4f;
	typedef CMat<double, 4, 4> Mat4x4d;
}