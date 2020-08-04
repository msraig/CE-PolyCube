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

#ifndef RESOURCEMANAGER_HH_
#define RESOURCEMANAGER_HH_

#include <string>
#include <vector>

#include "OpenVolumeMeshProperty.hh"
#include "PropertyHandles.hh"

namespace OpenVolumeMesh {

// Forward declarations
class BaseProperty;
template <class T>
class VertexPropertyT;
template <class T>
class EdgePropertyT;
template <class T>
class HalfEdgePropertyT;
template <class T>
class FacePropertyT;
template <class T>
class HalfFacePropertyT;
template <class T>
class CellPropertyT;
template <class T>
class MeshPropertyT;
template <class PropT, class HandleT>
class PropertyPtr;

class ResourceManager {
public:
    ResourceManager();
    virtual ~ResourceManager();

    template <class PropT, class HandleT> friend class PropertyPtr;

    /// Change size of stored vertex properties
    void resize_vprops(unsigned int _nv);

    /// Change size of stored edge properties
    void resize_eprops(unsigned int _ne);

    /// Change size of stored face properties
    void resize_fprops(unsigned int _nf);

    /// Change size of stored cell properties
    void resize_cprops(unsigned int _nc);

protected:

    void vertex_deleted(const VertexHandle& _h);

    void edge_deleted(const EdgeHandle& _h);

    void face_deleted(const FaceHandle& _h);

    void cell_deleted(const CellHandle& _h);

public:

    void clear_vertex_props() { clearVec(vertex_props_); }

    void clear_edge_props() { clearVec(edge_props_); }

    void clear_halfedge_props() { clearVec(halfedge_props_); }

    void clear_face_props() { clearVec(face_props_); }

    void clear_halfface_props() { clearVec(halfface_props_); }

    void clear_cell_props() { clearVec(cell_props_); }

    void clear_mesh_props() { clearVec(mesh_props_); }

    /// Get number of vertices in mesh
    virtual unsigned int n_vertices() const = 0;
    /// Get number of edges in mesh
    virtual unsigned int n_edges() const = 0;
    /// Get number of halfedges in mesh
    virtual unsigned int n_halfedges() const = 0;
    /// Get number of faces in mesh
    virtual unsigned int n_faces() const = 0;
    /// Get number of halffaces in mesh
    virtual unsigned int n_halffaces() const = 0;
    /// Get number of cells in mesh
    virtual unsigned int n_cells() const = 0;

    template<class T> VertexPropertyT<T> request_vertex_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> EdgePropertyT<T> request_edge_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> HalfEdgePropertyT<T> request_halfedge_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> FacePropertyT<T> request_face_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> HalfFacePropertyT<T> request_halfface_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> CellPropertyT<T> request_cell_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> MeshPropertyT<T> request_mesh_property(const std::string& _name = std::string(), const T _def = T());

private:

    void release_property(VertexPropHandle _handle);

    void release_property(EdgePropHandle _handle);

    void release_property(HalfEdgePropHandle _handle);

    void release_property(FacePropHandle _handle);

    void release_property(HalfFacePropHandle _handle);

    void release_property(CellPropHandle _handle);

    void release_property(MeshPropHandle _handle);

public:

    unsigned int n_vertex_props() const { return vertex_props_.size(); }

    unsigned int n_edge_props() const { return edge_props_.size(); }

    unsigned int n_halfedge_props() const { return halfedge_props_.size(); }

    unsigned int n_face_props() const { return face_props_.size(); }

    unsigned int n_halfface_props() const { return halfface_props_.size(); }

    unsigned int n_cell_props() const { return cell_props_.size(); }

    unsigned int n_mesh_props() const { return mesh_props_.size(); }

    template<class T> void set_persistent(VertexPropertyT<T>& _prop, bool _flag = true);

    template<class T> void set_persistent(EdgePropertyT<T>& _prop, bool _flag = true);

    template<class T> void set_persistent(HalfEdgePropertyT<T>& _prop, bool _flag = true);

    template<class T> void set_persistent(FacePropertyT<T>& _prop, bool _flag = true);

    template<class T> void set_persistent(HalfFacePropertyT<T>& _prop, bool _flag = true);

    template<class T> void set_persistent(CellPropertyT<T>& _prop, bool _flag = true);

    template<class T> void set_persistent(MeshPropertyT<T>& _prop, bool _flag = true);

    typedef std::vector<BaseProperty*> Properties;

    Properties::const_iterator vertex_props_begin() const { return vertex_props_.begin(); }

    Properties::const_iterator vertex_props_end() const { return vertex_props_.end(); }

    Properties::const_iterator edge_props_begin() const { return edge_props_.begin(); }

    Properties::const_iterator edge_props_end() const { return edge_props_.end(); }

    Properties::const_iterator halfedge_props_begin() const { return halfedge_props_.begin(); }

    Properties::const_iterator halfedge_props_end() const { return halfedge_props_.end(); }

    Properties::const_iterator face_props_begin() const { return face_props_.begin(); }

    Properties::const_iterator face_props_end() const { return face_props_.end(); }

    Properties::const_iterator halfface_props_begin() const { return halfface_props_.begin(); }

    Properties::const_iterator halfface_props_end() const { return halfface_props_.end(); }

    Properties::const_iterator cell_props_begin() const { return cell_props_.begin(); }

    Properties::const_iterator cell_props_end() const { return cell_props_.end(); }

    Properties::const_iterator mesh_props_begin() const { return mesh_props_.begin(); }

    Properties::const_iterator mesh_props_end() const { return mesh_props_.end(); }

    template <class FullPropT, class PropIterT>
    bool property_exists(const PropIterT& _begin, const PropIterT& _end, const std::string& _name) const {

        if(_name.length() == 0) {
            std::cerr << "Checking for the existance of anonymous properties is ambiguous!" << std::endl;
            return false;
        }

        PropIterT it = _begin;
        for(; it != _end; ++it) {
            if((*it)->name() == _name && dynamic_cast<FullPropT*>(*it) != NULL) {
                return true;
            }
        }
        return false;
    }

    template <class PropT>
    bool vertex_property_exists(const std::string& _name) const {
        return property_exists<VertexPropertyT<PropT> >(vertex_props_begin(), vertex_props_end(), _name);
    }

    template <class PropT>
    bool edge_property_exists(const std::string& _name) const {
        return property_exists<EdgePropertyT<PropT> >(edge_props_begin(), edge_props_end(), _name);
    }

    template <class PropT>
    bool halfedge_property_exists(const std::string& _name) const {
        return property_exists<HalfEdgePropertyT<PropT> >(halfedge_props_begin(), halfedge_props_end(), _name);
    }

    template <class PropT>
    bool face_property_exists(const std::string& _name) const {
        return property_exists<FacePropertyT<PropT> >(face_props_begin(), face_props_end(), _name);
    }

    template <class PropT>
    bool halfface_property_exists(const std::string& _name) const {
        return property_exists<HalfFacePropertyT<PropT> >(halfface_props_begin(), halfface_props_end(), _name);
    }

    template <class PropT>
    bool cell_property_exists(const std::string& _name) const {
        return property_exists<CellPropertyT<PropT> >(cell_props_begin(), cell_props_end(), _name);
    }

    template <class PropT>
    bool mesh_property_exists(const std::string& _name) const {
        return property_exists<MeshPropertyT<PropT> >(mesh_props_begin(), mesh_props_end(), _name);
    }

protected:

    void delete_multiple_vertex_props(const std::vector<bool>& _tags);

    void delete_multiple_edge_props(const std::vector<bool>& _tags);

    void delete_multiple_face_props(const std::vector<bool>& _tags);

    void delete_multiple_cell_props(const std::vector<bool>& _tags);

private:

    template<class StdVecT>
    void resize_props(StdVecT& _vec, unsigned int _n);

    template<class StdVecT>
    void entity_deleted(StdVecT& _vec, const OpenVolumeMeshHandle& _h);

    template<class StdVecT>
    void remove_property(StdVecT& _vec, size_t _idx);

    template<class StdVecT, class PropT, class HandleT, class T>
    PropT request_property(StdVecT& _vec, const std::string& _name, size_t _size, const T _def = T());

    template<class PropT>
    void set_persistentT(PropT& _prop, bool _flag);

    template<class StdVecT>
    void clearVec(StdVecT& _vec);

    Properties vertex_props_;

    Properties edge_props_;

    Properties halfedge_props_;

    Properties face_props_;

    Properties halfface_props_;

    Properties cell_props_;

    Properties mesh_props_;
};

}

#if defined(INCLUDE_TEMPLATES) && !defined(RESOURCEMANAGERT_CC)
#include "ResourceManagerT.cc"
#endif

#endif /* RESOURCEMANAGER_HH_ */
