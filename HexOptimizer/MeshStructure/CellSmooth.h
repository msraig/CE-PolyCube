#pragma once
#include <cstddef>

namespace MeshLib
{
	class CellCorner
	{
	public:
		CellCorner(size_t cid = 0, int c = 0, int ctag = true)
			:cell_id(cid), corner(c), cornertag(ctag)
		{}
	public:
		size_t cell_id;
		int corner;
		int cornertag;  // 0: corner; 1: 8 - 20 types; 2: 21 - 32 types
	};

	class CellFace
	{
	public:
		CellFace(size_t cid = 0, int c = 0)
			:cell_id(cid), face(c)
		{}
	public:
		size_t cell_id;
		int face;
	};
}