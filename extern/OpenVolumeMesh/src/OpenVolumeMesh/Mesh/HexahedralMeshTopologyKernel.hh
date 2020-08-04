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

#ifndef HEXAHEDRALMESHTOPOLOGYKERNEL_HH
#define HEXAHEDRALMESHTOPOLOGYKERNEL_HH

#include <set>

#include "../Core/TopologyKernel.hh"
#include "HexahedralMeshIterators.hh"

namespace OpenVolumeMesh {

/**
 * \class HexahedralMeshTopologyKernel
 *
 * \brief A data structure basing on PolyhedralMesh with specializations for hexahedra.
 *
 * The hexahedron has an induced "virtual" coordinate system. This supposes
 * the incident half-faces to be given in a specific order.
 * See the following figure for an illustration of the induced
 * coordinate system.
 *
 * \image html induced_coordsys.png
 *
 * The abbreviations XF, XB, etc. are short for
 *
 * \li \c XF: X-axis front face
 * \li \c XB: X-axis back face
 * \li \c YF: Y-axis front face
 * \li \c ...
 *
 * The axes refer to the intrinsic "virtual" axes of the hexahedron.
 * The incident half-faces have to be defined in the following order:
 *
 * \li \c 1. XF
 * \li \c 2. XB
 * \li \c 3. YF
 * \li \c 4. YB
 * \li \c 5. ZF
 * \li \c 6. ZB
 */

class HexahedralMeshTopologyKernel : public TopologyKernel {
public:

    // Orientation constants
    static const unsigned char XF = 0;
    static const unsigned char XB = 1;
    static const unsigned char YF = 2;
    static const unsigned char YB = 3;
    static const unsigned char ZF = 4;
    static const unsigned char ZB = 5;
    static const unsigned char INVALID = 6;

    static inline unsigned char opposite_orientation(const unsigned char _d) {
        return (_d % 2 == 0 ? _d + 1 : _d - 1);
    }

    // Constructor
    HexahedralMeshTopologyKernel();

    // Destructor
    ~HexahedralMeshTopologyKernel();

    // Overridden function
    virtual FaceHandle add_face(const std::vector<HalfEdgeHandle>& _halfedges, bool _topologyCheck = false);

    // Overridden function
    virtual FaceHandle add_face(const std::vector<VertexHandle>& _vertices);

    /// Overridden function
    virtual CellHandle add_cell(const std::vector<HalfFaceHandle>& _halffaces, bool _topologyCheck = false);

private:

    bool check_halfface_ordering(const std::vector<HalfFaceHandle>& _hfs) const;

public:

    /** \brief Add cell via incident vertices
     *
     * Test whether all required faces are already defined
     * and, if not, create them.
     * Give vertices in the following order:
     *
     *      5-------6
     *     /|      /|
     *    / |     / |
     *   3-------2  |
     *   |  4----|--7
     *   | /     | /
     *   |/      |/
     *   0-------1
     *
     * @param _vertices A list of vertices in the correct order
     *
     * @return The new hexahedron's cell handle
     */
    CellHandle add_cell(const std::vector<VertexHandle>& _vertices);

    // ======================= Specialized Iterators =============================

    friend class CellSheetCellIter;
    friend class HalfFaceSheetHalfFaceIter;
    friend class HexVertexIter;

    typedef class CellSheetCellIter CellSheetCellIter;
    typedef class HalfFaceSheetHalfFaceIter HalfFaceSheetHalfFaceIter;
    typedef class HexVertexIter HexVertexIter;

    CellSheetCellIter csc_iter(const CellHandle& _ref_h, const unsigned char _orthDir) const {
        return CellSheetCellIter(_ref_h, _orthDir, this);
    }

    HalfFaceSheetHalfFaceIter hfshf_iter(const HalfFaceHandle& _ref_h) const {
        return HalfFaceSheetHalfFaceIter(_ref_h, this);
    }

    HexVertexIter hv_iter(const CellHandle& _ref_h) const {
        return HexVertexIter(_ref_h, this);
    }

    // ======================= Connectivity functions =============================

    inline HalfFaceHandle opposite_halfface_handle_in_cell(const HalfFaceHandle& _hfh, const CellHandle& _ch) const {

        assert((unsigned int)_ch.idx() < TopologyKernel::cells_.size());

        if(orientation(_hfh, _ch) == XF) return xback_halfface(_ch);
        if(orientation(_hfh, _ch) == XB) return xfront_halfface(_ch);
        if(orientation(_hfh, _ch) == YF) return yback_halfface(_ch);
        if(orientation(_hfh, _ch) == YB) return yfront_halfface(_ch);
        if(orientation(_hfh, _ch) == ZF) return zback_halfface(_ch);
        if(orientation(_hfh, _ch) == ZB) return zfront_halfface(_ch);

        return TopologyKernel::InvalidHalfFaceHandle;
    }

    inline HalfFaceHandle xfront_halfface(const CellHandle& _ch) const {

        assert((unsigned int)_ch.idx() < TopologyKernel::cells_.size());

        return TopologyKernel::cell(_ch).halffaces()[XF];
    }

    inline HalfFaceHandle xback_halfface(const CellHandle& _ch) const {

        assert((unsigned int)_ch.idx() < TopologyKernel::cells_.size());

        return TopologyKernel::cell(_ch).halffaces()[XB];
    }

    inline HalfFaceHandle yfront_halfface(const CellHandle& _ch) const {

        assert((unsigned int)_ch.idx() < TopologyKernel::cells_.size());

        return TopologyKernel::cell(_ch).halffaces()[YF];
    }

    inline HalfFaceHandle yback_halfface(const CellHandle& _ch) const {

        assert((unsigned int)_ch.idx() < TopologyKernel::cells_.size());

        return TopologyKernel::cell(_ch).halffaces()[YB];
    }

    inline HalfFaceHandle zfront_halfface(const CellHandle& _ch) const {

        assert((unsigned int)_ch.idx() < TopologyKernel::cells_.size());

        return TopologyKernel::cell(_ch).halffaces()[ZF];
    }

    inline HalfFaceHandle zback_halfface(const CellHandle& _ch) const {

        assert((unsigned int)_ch.idx() < TopologyKernel::cells_.size());

        return TopologyKernel::cell(_ch).halffaces()[ZB];
    }

    unsigned char orientation(const HalfFaceHandle& _hfh, const CellHandle& _ch) const {

        assert((unsigned int)_ch.idx() < TopologyKernel::cells_.size());

        std::vector<HalfFaceHandle> halffaces = TopologyKernel::cell(_ch).halffaces();
        for(unsigned int i = 0; i < halffaces.size(); ++i) {
            if(halffaces[i] == _hfh) return (unsigned char)i;
        }

        return INVALID;
    }

    static inline unsigned char orthogonal_orientation(const unsigned char _o1, const unsigned char _o2) {

        if(_o1 == XF && _o2 == YF) return ZF;
        if(_o1 == XF && _o2 == YB) return ZB;
        if(_o1 == XF && _o2 == ZF) return YB;
        if(_o1 == XF && _o2 == ZB) return YF;
        if(_o1 == XB && _o2 == YF) return ZB;
        if(_o1 == XB && _o2 == YB) return ZF;
        if(_o1 == XB && _o2 == ZF) return YF;
        if(_o1 == XB && _o2 == ZB) return YB;

        if(_o1 == YF && _o2 == XF) return ZB;
        if(_o1 == YF && _o2 == XB) return ZF;
        if(_o1 == YF && _o2 == ZF) return XF;
        if(_o1 == YF && _o2 == ZB) return XB;
        if(_o1 == YB && _o2 == XF) return ZF;
        if(_o1 == YB && _o2 == XB) return ZB;
        if(_o1 == YB && _o2 == ZF) return XB;
        if(_o1 == YB && _o2 == ZB) return XF;

        if(_o1 == ZF && _o2 == YF) return XB;
        if(_o1 == ZF && _o2 == YB) return XF;
        if(_o1 == ZF && _o2 == XF) return YF;
        if(_o1 == ZF && _o2 == XB) return YB;
        if(_o1 == ZB && _o2 == YF) return XF;
        if(_o1 == ZB && _o2 == YB) return XB;
        if(_o1 == ZB && _o2 == XF) return YB;
        if(_o1 == ZB && _o2 == XB) return YF;

        return INVALID;

    }

    inline HalfFaceHandle get_oriented_halfface(const unsigned char _o, const CellHandle& _ch) const {

        if(_o == XF) return xfront_halfface(_ch);
        if(_o == XB) return xback_halfface(_ch);
        if(_o == YF) return yfront_halfface(_ch);
        if(_o == YB) return yback_halfface(_ch);
        if(_o == ZF) return zfront_halfface(_ch);
        if(_o == ZB) return zback_halfface(_ch);
        return TopologyKernel::InvalidHalfFaceHandle;
    }

    HalfFaceHandle adjacent_halfface_on_sheet(const HalfFaceHandle& _hfh, const HalfEdgeHandle& _heh) const {

        if(!TopologyKernel::has_face_bottom_up_incidences()) {
            std::cerr << "No bottom-up incidences computed so far, could not get adjacent halfface on sheet!" << std::endl;
            return TopologyKernel::InvalidHalfFaceHandle;
        }

        HalfFaceHandle n_hf = _hfh;
        HalfEdgeHandle n_he = _heh;

        // Try the 1st way
        while(true) {
            n_hf = TopologyKernel::adjacent_halfface_in_cell(n_hf, n_he);
            if(n_hf == TopologyKernel::InvalidHalfFaceHandle) break;
            n_hf = TopologyKernel::opposite_halfface_handle(n_hf);
            if(n_hf == TopologyKernel::InvalidHalfFaceHandle) break;
            HalfEdgeHandle o_he = TopologyKernel::opposite_halfedge_handle(n_he);
            if(o_he == TopologyKernel::InvalidHalfEdgeHandle) break;
            n_hf = TopologyKernel::adjacent_halfface_in_cell(n_hf, o_he);
            if(n_hf == TopologyKernel::InvalidHalfFaceHandle) break;
            else return n_hf;
        }

        n_hf = TopologyKernel::opposite_halfface_handle(_hfh);
        n_he = TopologyKernel::opposite_halfedge_handle(_heh);

        // Try the 2nd way
        while(true) {
            n_hf = TopologyKernel::adjacent_halfface_in_cell(n_hf, n_he);
            if(n_hf == TopologyKernel::InvalidHalfFaceHandle) break;
            n_hf = TopologyKernel::opposite_halfface_handle(n_hf);
            if(n_hf == TopologyKernel::InvalidHalfFaceHandle) break;
            HalfEdgeHandle o_he = TopologyKernel::opposite_halfedge_handle(n_he);
            if(o_he == TopologyKernel::InvalidHalfEdgeHandle) break;
            n_hf = TopologyKernel::adjacent_halfface_in_cell(n_hf, o_he);
            if(n_hf == TopologyKernel::InvalidHalfFaceHandle) break;
            else return TopologyKernel::opposite_halfface_handle(n_hf);
        }

        return TopologyKernel::InvalidHalfFaceHandle;
    }

    HalfFaceHandle adjacent_halfface_on_surface(const HalfFaceHandle& _hfh, const HalfEdgeHandle& _heh) const {

        for(OpenVolumeMesh::HalfEdgeHalfFaceIter hehf_it = TopologyKernel::hehf_iter(_heh);
                hehf_it.valid(); ++hehf_it) {
            if(*hehf_it == _hfh) continue;
            if(TopologyKernel::is_boundary(*hehf_it)) {
                return *hehf_it;
            }
            if(TopologyKernel::is_boundary(TopologyKernel::opposite_halfface_handle(*hehf_it))) {
                return TopologyKernel::opposite_halfface_handle(*hehf_it);
            }
        }
        return TopologyKernel::InvalidHalfFaceHandle;
    }

    HalfFaceHandle neighboring_outside_halfface(const HalfFaceHandle& _hfh, const HalfEdgeHandle& _heh) const {

        if(!TopologyKernel::has_face_bottom_up_incidences()) {
            std::cerr << "No bottom-up incidences computed so far, could not get neighboring outside halfface!" << std::endl;
            return TopologyKernel::InvalidHalfFaceHandle;
        }

        for(OpenVolumeMesh::HalfEdgeHalfFaceIter hehf_it = TopologyKernel::hehf_iter(_heh);
                hehf_it; ++hehf_it) {
            if(*hehf_it == _hfh) continue;
            if(TopologyKernel::is_boundary(*hehf_it)) return *hehf_it;
            if(TopologyKernel::is_boundary(TopologyKernel::opposite_halfface_handle(*hehf_it)))
                return TopologyKernel::opposite_halfface_handle(*hehf_it);
        }

        return TopologyKernel::InvalidHalfFaceHandle;
    }

private:

    const HalfFaceHandle& get_adjacent_halfface(const HalfFaceHandle& _hfh, const HalfEdgeHandle& _heh,
            const std::vector<HalfFaceHandle>& _halffaces) const;

};

} // Namespace OpenVolumeMesh

#endif /* HEXAHEDRALMESHTOPOLOGYKERNEL_HH */
