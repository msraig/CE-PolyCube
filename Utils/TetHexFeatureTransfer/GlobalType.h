#pragma once
//define type for HexStructure

#define Interior_RegularE 4
#define Boundary_RegularE 2

namespace HexStructure
{
	typedef unsigned Int;
	typedef double Real;
	//const int hex_face_table[6][4] =
	//{
	//	{ 0,1,2,3 },
	//	{ 4,5,6,7 },
	//	{ 0,1,5,4 },
	//	{ 0,4,7,3 },
	//	{ 3,2,6,7 },
	//	{ 1,5,6,2 },
	//};

	//share the same order as cell_faces list
	const int hex_face_table[6][4] =
	{
		{ 3,2,1,0 },
		{ 7,6,5,4 },
		//{ 7,4,5,6 },
		{ 0,1,5,4 },
		{ 3,2,6,7 },
		{ 0,3,7,4 },
		{ 1,2,6,5 },
	};

	//outside shares the same order as hex face table
	const int hex_face_outside[6][4] =
	{
		{ 3,2,1,0 },
		{ 7,4,5,6 },
		{ 0,1,5,4 },
		{ 7,6,2,3 },
		{ 0,4,7,3 },
		{ 2,6,5,1 },
	};



	const int hex_tetra_table[8][4] =
	{
		{ 0,3,4,1 },
		{ 1,0,5,2 },
		{ 2,1,6,3 },
		{ 3,2,7,0 },
		{ 4,7,5,0 },
		{ 5,4,6,1 },
		{ 6,5,7,2 },
		{ 7,6,4,3 },
	};
}

