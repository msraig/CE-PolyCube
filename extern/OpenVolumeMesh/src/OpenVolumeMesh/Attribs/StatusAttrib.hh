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

#ifndef STATUSATTRIB_HH_
#define STATUSATTRIB_HH_

#include <cassert>

#include "../Core/OpenVolumeMeshHandle.hh"
#include "OpenVolumeMeshStatus.hh"
#include "../Core/PropertyDefines.hh"

namespace OpenVolumeMesh {

// Forward declaration
class TopologyKernel;

class StatusAttrib {
public:
    explicit StatusAttrib(TopologyKernel& _kernel);
    ~StatusAttrib();

    const OpenVolumeMeshStatus& operator[](const VertexHandle& _h) const {
        return v_status_[_h.idx()];
    }

    OpenVolumeMeshStatus& operator[](const VertexHandle& _h) {
        return v_status_[_h.idx()];
    }

    const OpenVolumeMeshStatus& operator[](const EdgeHandle& _h) const {
        return e_status_[_h.idx()];
    }

    OpenVolumeMeshStatus& operator[](const EdgeHandle& _h) {
        return e_status_[_h.idx()];
    }

    const OpenVolumeMeshStatus& operator[](const HalfEdgeHandle& _h) const {
        return he_status_[_h.idx()];
    }

    OpenVolumeMeshStatus& operator[](const HalfEdgeHandle& _h) {
        return he_status_[_h.idx()];
    }

    const OpenVolumeMeshStatus& operator[](const FaceHandle& _h) const {
        return f_status_[_h.idx()];
    }

    OpenVolumeMeshStatus& operator[](const FaceHandle& _h) {
        return f_status_[_h.idx()];
    }

    const OpenVolumeMeshStatus& operator[](const HalfFaceHandle& _h) const {
        return hf_status_[_h.idx()];
    }

    OpenVolumeMeshStatus& operator[](const HalfFaceHandle& _h) {
        return hf_status_[_h.idx()];
    }

    const OpenVolumeMeshStatus& operator[](const CellHandle& _h) const {
        return c_status_[_h.idx()];
    }

    OpenVolumeMeshStatus& operator[](const CellHandle& _h) {
        return c_status_[_h.idx()];
    }

    const OpenVolumeMeshStatus& mesh_status() const {
        OpenVolumeMeshHandle h(0);
        return m_status_[h.idx()];
    }

    OpenVolumeMeshStatus& mesh_status() {
        OpenVolumeMeshHandle h(0);
        return m_status_[h.idx()];
    }

    typedef VertexPropertyT<OpenVolumeMeshStatus>::const_iterator   const_vstatus_iterator;
    typedef VertexPropertyT<OpenVolumeMeshStatus>::iterator         vstatus_iterator;
    typedef EdgePropertyT<OpenVolumeMeshStatus>::const_iterator     const_estatus_iterator;
    typedef EdgePropertyT<OpenVolumeMeshStatus>::iterator           estatus_iterator;
    typedef HalfEdgePropertyT<OpenVolumeMeshStatus>::const_iterator const_hestatus_iterator;
    typedef HalfEdgePropertyT<OpenVolumeMeshStatus>::iterator       hestatus_iterator;
    typedef FacePropertyT<OpenVolumeMeshStatus>::const_iterator     const_fstatus_iterator;
    typedef FacePropertyT<OpenVolumeMeshStatus>::iterator           fstatus_iterator;
    typedef HalfFacePropertyT<OpenVolumeMeshStatus>::const_iterator const_hfstatus_iterator;
    typedef HalfFacePropertyT<OpenVolumeMeshStatus>::iterator       hfstatus_iterator;
    typedef CellPropertyT<OpenVolumeMeshStatus>::const_iterator     const_cstatus_iterator;
    typedef CellPropertyT<OpenVolumeMeshStatus>::iterator           cstatus_iterator;

    // Iterator access
    VertexPropertyT<OpenVolumeMeshStatus>::const_iterator vstatus_begin() const {
        return v_status_.begin();
    }
    VertexPropertyT<OpenVolumeMeshStatus>::iterator vstatus_begin() {
        return v_status_.begin();
    }
    VertexPropertyT<OpenVolumeMeshStatus>::const_iterator vstatus_end() const {
        return v_status_.end();
    }
    VertexPropertyT<OpenVolumeMeshStatus>::iterator vstatus_end() {
        return v_status_.end();
    }

    EdgePropertyT<OpenVolumeMeshStatus>::const_iterator estatus_begin() const {
        return e_status_.begin();
    }
    EdgePropertyT<OpenVolumeMeshStatus>::iterator estatus_begin() {
        return e_status_.begin();
    }
    EdgePropertyT<OpenVolumeMeshStatus>::const_iterator estatus_end() const {
        return e_status_.end();
    }
    EdgePropertyT<OpenVolumeMeshStatus>::iterator estatus_end() {
        return e_status_.end();
    }

    HalfEdgePropertyT<OpenVolumeMeshStatus>::const_iterator hestatus_begin() const {
        return he_status_.begin();
    }
    HalfEdgePropertyT<OpenVolumeMeshStatus>::iterator hestatus_begin() {
        return he_status_.begin();
    }
    HalfEdgePropertyT<OpenVolumeMeshStatus>::const_iterator hestatus_end() const {
        return he_status_.end();
    }
    HalfEdgePropertyT<OpenVolumeMeshStatus>::iterator hestatus_end() {
        return he_status_.end();
    }

    FacePropertyT<OpenVolumeMeshStatus>::const_iterator fstatus_begin() const {
        return f_status_.begin();
    }
    FacePropertyT<OpenVolumeMeshStatus>::iterator fstatus_begin() {
        return f_status_.begin();
    }
    FacePropertyT<OpenVolumeMeshStatus>::const_iterator fstatus_end() const {
        return f_status_.end();
    }
    FacePropertyT<OpenVolumeMeshStatus>::iterator fstatus_end() {
        return f_status_.end();
    }

    HalfFacePropertyT<OpenVolumeMeshStatus>::const_iterator hfstatus_begin() const {
        return hf_status_.begin();
    }
    HalfFacePropertyT<OpenVolumeMeshStatus>::iterator hfstatus_begin() {
        return hf_status_.begin();
    }
    HalfFacePropertyT<OpenVolumeMeshStatus>::const_iterator hfstatus_end() const {
        return hf_status_.end();
    }
    HalfFacePropertyT<OpenVolumeMeshStatus>::iterator hfstatus_end() {
        return hf_status_.end();
    }

    CellPropertyT<OpenVolumeMeshStatus>::const_iterator cstatus_begin() const {
        return c_status_.begin();
    }
    CellPropertyT<OpenVolumeMeshStatus>::iterator cstatus_begin() {
        return c_status_.begin();
    }
    CellPropertyT<OpenVolumeMeshStatus>::const_iterator cstatus_end() const {
        return c_status_.end();
    }
    CellPropertyT<OpenVolumeMeshStatus>::iterator cstatus_end() {
        return c_status_.end();
    }

    /**
     * \brief Delete all entities that have been marked as deleted
     *
     * This function deletes all entities that have been marked as deleted.
     * It proceeds bottom-up, starting with the vertices. All higher
     * dimensional entities that are incident to a deleted entity are
     * automatically marked deleted, too. Once this first pass is through,
     * one can additionally delete all resulting non-manifold configurations
     * in a second pass (triggered by the parameter of this function).
     * This step proceeds as follows: Delete all n-dimensional entities
     * (starting with n = 2), that are not incident to at least one
     * entity of dimension n + 1. Note that the second pass requires bottom-up
     * incidences to be available. Compute them by calling update_incidences().
     *
     * @param _preserveManifoldness Pass true if the mesh is required to stay three-manifold
     */
    void garbage_collection(bool _preserveManifoldness = false);

private:

    void mark_higher_dim_entities();

    TopologyKernel& kernel_;

    VertexPropertyT<OpenVolumeMeshStatus> v_status_;
    EdgePropertyT<OpenVolumeMeshStatus> e_status_;
    HalfEdgePropertyT<OpenVolumeMeshStatus> he_status_;
    FacePropertyT<OpenVolumeMeshStatus> f_status_;
    HalfFacePropertyT<OpenVolumeMeshStatus> hf_status_;
    CellPropertyT<OpenVolumeMeshStatus> c_status_;
    MeshPropertyT<OpenVolumeMeshStatus> m_status_;
};

} // Namespace OpenVolumeMesh

#endif /* STATUSATTRIB_HH_ */
