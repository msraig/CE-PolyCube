/*===========================================================================*\
 *                                                                           *
 *                            OpenVolumeMesh                                 *
 *        Copyright (C) 2011 by Computer Graphics Group, RWTH Aachen         *
 *                        www.openvolumemesh.org                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of OpenVolumeMesh.                                     *
 *                                                                           *
 *  OpenVolumeMesh is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenVolumeMesh is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenVolumeMesh.  If not,                              *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *
 *   $Revision$                                                         *
 *   $Date$                    *
 *   $LastChangedBy$                                                *
 *                                                                           *
\*===========================================================================*/

#ifndef HEXAHEDRALMESHITERATORS_HH
#define HEXAHEDRALMESHITERATORS_HH

#include "../Core/Iterators.hh"

namespace OpenVolumeMesh {

class HexahedralMeshTopologyKernel;


class CellSheetCellIter : public BaseIterator<CellHandle,
	CellHandle> {
private:
	typedef BaseIterator<CellHandle,
			CellHandle>    BaseIter;
public:
	CellSheetCellIter(const CellHandle& _ref_h, const unsigned char _orthDir,
			const HexahedralMeshTopologyKernel* _mesh);

	CellSheetCellIter& operator=(const CellSheetCellIter& _c) {
		BaseIter::operator=(_c);
		neighb_sheet_cell_hs_ = _c.neighb_sheet_cell_hs_;
		cur_it_ = neighb_sheet_cell_hs_.begin();
		return *this;
	}

	// Post increment/decrement operator
	CellSheetCellIter operator++(int) {
		CellSheetCellIter cpy = *this;
		++(*this);
		return cpy;
	}
	CellSheetCellIter operator--(int) {
		CellSheetCellIter cpy = *this;
		--(*this);
		return cpy;
	}
	CellSheetCellIter operator+(int _n) {
		CellSheetCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	CellSheetCellIter operator-(int _n) {
		CellSheetCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	CellSheetCellIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	CellSheetCellIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	CellSheetCellIter& operator++();
	CellSheetCellIter& operator--();

private:
	std::set<CellHandle> neighb_sheet_cell_hs_;
	std::set<CellHandle>::const_iterator cur_it_;
};


class HalfFaceSheetHalfFaceIter : public BaseIterator<HalfFaceHandle,
	HalfFaceHandle> {
private:
	typedef BaseIterator<HalfFaceHandle,
			HalfFaceHandle> BaseIter;
public:
	HalfFaceSheetHalfFaceIter(const HalfFaceHandle& _ref_h,
			const HexahedralMeshTopologyKernel* _mesh);
	HalfFaceSheetHalfFaceIter& operator=(const HalfFaceSheetHalfFaceIter& _c) {
		BaseIter::operator=(_c);
		adjacent_halffaces_ = _c.adjacent_halffaces_;
		cur_it_ = adjacent_halffaces_.begin();
		return *this;
	}

	// Post increment/decrement operator
	HalfFaceSheetHalfFaceIter operator++(int) {
		HalfFaceSheetHalfFaceIter cpy = *this;
		++(*this);
		return cpy;
	}
	HalfFaceSheetHalfFaceIter operator--(int) {
		HalfFaceSheetHalfFaceIter cpy = *this;
		--(*this);
		return cpy;
	}
	HalfFaceSheetHalfFaceIter operator+(int _n) {
		HalfFaceSheetHalfFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	HalfFaceSheetHalfFaceIter operator-(int _n) {
		HalfFaceSheetHalfFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	HalfFaceSheetHalfFaceIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	HalfFaceSheetHalfFaceIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	HalfFaceSheetHalfFaceIter& operator++();
	HalfFaceSheetHalfFaceIter& operator--();

	const EdgeHandle& common_edge() const { return *edge_it_; }

private:
	std::vector<HalfFaceHandle> adjacent_halffaces_;
	std::vector<HalfFaceHandle>::const_iterator cur_it_;
	std::vector<EdgeHandle> common_edges_;
    std::vector<EdgeHandle>::const_iterator edge_it_;
};

/** \brief Iterate over all vertices of a hexahedron in a specific order
 *
 * Vertices are addressed in the following order:
 *
 *      5-------6
 *     /|      /|
 *    / |     / |
 *   3-------2  |
 *   |  4----|--7
 *   | /     | /
 *   |/      |/
 *   0-------1
 */

class HexVertexIter : public BaseIterator<CellHandle,
    VertexHandle> {
private:
    typedef BaseIterator<CellHandle,
            VertexHandle> BaseIter;
public:
    HexVertexIter(const CellHandle& _ref_h,
            const HexahedralMeshTopologyKernel* _mesh);
    HexVertexIter& operator=(const HexVertexIter& _c) {
        BaseIter::operator=(_c);
        vertices_ = _c.vertices_;
        cur_it_ = vertices_.begin();
        return *this;
    }

    // Post increment/decrement operator
    HexVertexIter operator++(int) {
        HexVertexIter cpy = *this;
        ++(*this);
        return cpy;
    }
    HexVertexIter operator--(int) {
        HexVertexIter cpy = *this;
        --(*this);
        return cpy;
    }
    HexVertexIter operator+(int _n) {
        HexVertexIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    HexVertexIter operator-(int _n) {
        HexVertexIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    HexVertexIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    HexVertexIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    HexVertexIter& operator++();
    HexVertexIter& operator--();

private:
    std::vector<VertexHandle> vertices_;
    std::vector<VertexHandle>::const_iterator cur_it_;
};

} // Namespace OpenVolumeMesh

#endif /* HEXAHEDRALMESHITERATORS_HH */
