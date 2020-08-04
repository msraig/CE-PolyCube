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

#include "ResourceManager.hh"
#include "BaseProperty.hh"

namespace OpenVolumeMesh {

ResourceManager::ResourceManager() {
}

ResourceManager::~ResourceManager() {

    // Delete persistent props
    clearVec(vertex_props_);
    clearVec(edge_props_);
    clearVec(halfedge_props_);
    clearVec(face_props_);
    clearVec(halfface_props_);
    clearVec(cell_props_);
    clearVec(mesh_props_);
}

void ResourceManager::resize_vprops(unsigned int _nv) {

    resize_props(vertex_props_, _nv);
}

void ResourceManager::resize_eprops(unsigned int _ne) {

    resize_props(edge_props_, _ne);
    resize_props(halfedge_props_, _ne*2u);
}

void ResourceManager::resize_fprops(unsigned int _nf) {

    resize_props(face_props_, _nf);
    resize_props(halfface_props_, _nf*2u);
}

void ResourceManager::resize_cprops(unsigned int _nc) {

    resize_props(cell_props_, _nc);
}

void ResourceManager::vertex_deleted(const VertexHandle& _h) {

    entity_deleted(vertex_props_, _h);
}

void ResourceManager::edge_deleted(const EdgeHandle& _h) {

    entity_deleted(edge_props_, _h);
    entity_deleted(halfedge_props_, OpenVolumeMeshHandle(_h.idx()*2 + 1));
    entity_deleted(halfedge_props_, OpenVolumeMeshHandle(_h.idx()*2));
}

void ResourceManager::face_deleted(const FaceHandle& _h) {

    entity_deleted(face_props_, _h);
    entity_deleted(halfface_props_, OpenVolumeMeshHandle(_h.idx()*2 + 1));
    entity_deleted(halfface_props_, OpenVolumeMeshHandle(_h.idx()*2));
}

void ResourceManager::cell_deleted(const CellHandle& _h) {

    entity_deleted(cell_props_, _h);
}

void ResourceManager::release_property(VertexPropHandle _handle) {

    remove_property(vertex_props_, _handle.idx());
}

void ResourceManager::release_property(EdgePropHandle _handle) {

    remove_property(edge_props_, _handle.idx());
}

void ResourceManager::release_property(HalfEdgePropHandle _handle) {

    remove_property(halfedge_props_, _handle.idx());
}

void ResourceManager::release_property(FacePropHandle _handle) {

    remove_property(face_props_, _handle.idx());
}

void ResourceManager::release_property(HalfFacePropHandle _handle) {

    remove_property(halfface_props_, _handle.idx());
}

void ResourceManager::release_property(CellPropHandle _handle) {

    remove_property(cell_props_, _handle.idx());
}

void ResourceManager::release_property(MeshPropHandle _handle) {

    remove_property(mesh_props_, _handle.idx());
}

void ResourceManager::delete_multiple_vertex_props(const std::vector<bool>& _tags) {

    Properties::iterator vp_it = vertex_props_.begin();
    Properties::iterator vp_end = vertex_props_.end();
    for(; vp_it != vp_end; ++vp_it) {
        (*vp_it)->delete_multiple_entries(_tags);
    }
}

void ResourceManager::delete_multiple_edge_props(const std::vector<bool>& _tags) {

    Properties::iterator ep_it = edge_props_.begin();
    Properties::iterator ep_end = edge_props_.end();
    for(; ep_it != ep_end; ++ep_it) {
        (*ep_it)->delete_multiple_entries(_tags);
    }
    // Create tags vector for halfedges
    std::vector<bool> hetags;
    for(std::vector<bool>::const_iterator t_it = _tags.begin(),
            t_end = _tags.end(); t_it != t_end; ++t_it) {
        hetags.push_back(*t_it);
        hetags.push_back(*t_it);
    }
    Properties::iterator hep_it = halfedge_props_.begin();
    Properties::iterator hep_end = halfedge_props_.end();
    for(; hep_it != hep_end; ++hep_it) {
        (*hep_it)->delete_multiple_entries(hetags);
    }
}

void ResourceManager::delete_multiple_face_props(const std::vector<bool>& _tags) {

    Properties::iterator fp_it = face_props_.begin();
    Properties::iterator fp_end = face_props_.end();
    for(; fp_it != fp_end; ++fp_it) {
        (*fp_it)->delete_multiple_entries(_tags);
    }
    // Create tags vector for halffaces
    std::vector<bool> hftags;
    for(std::vector<bool>::const_iterator t_it = _tags.begin(),
            t_end = _tags.end(); t_it != t_end; ++t_it) {
        hftags.push_back(*t_it);
        hftags.push_back(*t_it);
    }
    Properties::iterator hfp_it = halfface_props_.begin();
    Properties::iterator hfp_end = halfface_props_.end();
    for(; hfp_it != hfp_end; ++hfp_it) {
        (*hfp_it)->delete_multiple_entries(hftags);
    }
}

void ResourceManager::delete_multiple_cell_props(const std::vector<bool>& _tags) {

    Properties::iterator cp_it = cell_props_.begin();
    Properties::iterator cp_end = cell_props_.end();
    for(; cp_it != cp_end; ++cp_it) {
        (*cp_it)->delete_multiple_entries(_tags);
    }
}

} // Namespace OpenVolumeMesh
