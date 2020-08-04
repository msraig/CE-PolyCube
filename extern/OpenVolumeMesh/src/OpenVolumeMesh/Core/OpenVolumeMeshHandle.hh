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

#ifndef OPENVolumeMESHHANDLE_HH_
#define OPENVolumeMESHHANDLE_HH_

#include <iostream>
#include <vector>
#include <OpenVolumeMesh/System/FunctionalInclude.hh>
#include <algorithm>

namespace OpenVolumeMesh {

// Define handle types in order to distinguish different entities by their indices
class OpenVolumeMeshHandle {
public:
    // Default constructor
	explicit OpenVolumeMeshHandle(int _idx) : idx_(_idx) {};

	OpenVolumeMeshHandle& operator=(int _idx) {
		idx_ = _idx;
		return *this;
	}

	OpenVolumeMeshHandle& operator=(const OpenVolumeMeshHandle& _idx) {
		idx_ = _idx.idx_;
		return *this;
	}

	inline bool is_valid() const { return idx_ != -1; }

	inline bool operator<(const OpenVolumeMeshHandle& _idx) const { return (this->idx_ < _idx.idx_); }

	inline bool operator<(int _idx) const { return idx_ < _idx; }

	inline bool operator>(const OpenVolumeMeshHandle& _idx) const { return (this->idx_ > _idx.idx_); }

    inline bool operator>(int _idx) const { return idx_ > _idx; }

	inline bool operator==(const OpenVolumeMeshHandle& _h) const { return _h.idx_ == this->idx_; }

	inline bool operator!=(const OpenVolumeMeshHandle& _h) const { return _h.idx_ != this->idx_; }

	inline const int& idx() const { return idx_; }

	inline void idx(const int& _idx) { idx_ = _idx; }

	inline operator int() const { return idx_; }

	void reset() { idx_ = -1; }

	inline void cmp_idx(const int& id, bool& r)
	{
		r = id < idx_ ? true : false;
	}
	inline void decrease_one() { idx_ -= 1;}
	inline void decrease_two() { idx_ -= 2;}

private:
	int idx_;
};

// Default entity handles

class VertexHandle   : public OpenVolumeMeshHandle { public: VertexHandle(int _idx = -1)   : OpenVolumeMeshHandle(_idx) {} };
class EdgeHandle     : public OpenVolumeMeshHandle { public: EdgeHandle(int _idx = -1)     : OpenVolumeMeshHandle(_idx) {} };
class FaceHandle     : public OpenVolumeMeshHandle { public: FaceHandle(int _idx = -1)     : OpenVolumeMeshHandle(_idx) {} };
class CellHandle     : public OpenVolumeMeshHandle { public: CellHandle(int _idx = -1)     : OpenVolumeMeshHandle(_idx) {} };
class HalfEdgeHandle : public OpenVolumeMeshHandle { public: HalfEdgeHandle(int _idx = -1) : OpenVolumeMeshHandle(_idx) {} };
class HalfFaceHandle : public OpenVolumeMeshHandle { public: HalfFaceHandle(int _idx = -1) : OpenVolumeMeshHandle(_idx) {} };

// Helper class that is used to decrease all handles
// exceeding a certain threshold

class VHandleCorrection {
public:
    VHandleCorrection(VertexHandle _thld) : thld_(_thld) {}
    void correctValue(VertexHandle& _h) {
        if(_h > thld_) _h.idx(_h.idx() - 1);
    }
private:
    VertexHandle thld_;
};
class HEHandleCorrection {
public:
    HEHandleCorrection(HalfEdgeHandle _thld) : thld_(_thld) {}
    void correctVecValue(std::vector<HalfEdgeHandle>& _vec) {
        std::for_each(_vec.begin(), _vec.end(), fun::bind(&HEHandleCorrection::correctValue, this, fun::placeholders::_1));
    }
    void correctValue(HalfEdgeHandle& _h) {
        if(_h > thld_) _h.idx(_h.idx() - 2);
    }
private:
    HalfEdgeHandle thld_;
};
class HFHandleCorrection {
public:
    HFHandleCorrection(HalfFaceHandle _thld) : thld_(_thld) {}
    void correctVecValue(std::vector<HalfFaceHandle>& _vec) {
        std::for_each(_vec.begin(), _vec.end(), fun::bind(&HFHandleCorrection::correctValue, this, fun::placeholders::_1));
    }
    void correctValue(HalfFaceHandle& _h) {
        if(_h > thld_) _h.idx(_h.idx() - 2);
    }
private:
    HalfFaceHandle thld_;
};
class CHandleCorrection {
public:
    CHandleCorrection(CellHandle _thld) : thld_(_thld) {}
    void correctValue(CellHandle& _h) {
        if(_h > thld_) _h.idx(_h.idx() - 1);
    }
private:
    CellHandle thld_;
};

bool operator==(const int& _lhs, const OpenVolumeMeshHandle& _rhs);

bool operator==(const unsigned int& _lhs, const OpenVolumeMeshHandle& _rhs);

bool operator!=(const int& _lhs, const OpenVolumeMeshHandle& _rhs);

bool operator!=(const unsigned int& _lhs, const OpenVolumeMeshHandle& _rhs);

std::ostream& operator<<(std::ostream& _ostr, const OpenVolumeMeshHandle& _handle);

std::istream& operator>>(std::istream& _istr, OpenVolumeMeshHandle& _handle);

} // Namespace OpenVolumeMesh

#endif /* OPENVolumeMESHHANDLE_HH_ */
