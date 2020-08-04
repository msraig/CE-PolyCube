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

#define FILEMANAGERT_CC

#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <typeinfo>

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

#include "FileManager.hh"

namespace OpenVolumeMesh {

namespace IO {

//==================================================

template <class MeshT>
bool FileManager::readFile(const std::string& _filename, MeshT& _mesh,
    bool _topologyCheck, bool _computeBottomUpIncidences) const {

    std::ifstream iff(_filename.c_str(), std::ios::in);

    if(!iff.good()) {
        std::cerr << "Error: Could not open file " << _filename << " for reading!" << std::endl;
        iff.close();
        return false;
    }

    std::stringstream sstr;
    std::string line;
    std::string s_tmp;
    unsigned int c = 0u;
    typedef typename MeshT::PointT Point;
    Point v = Point(0.0, 0.0, 0.0);
    unsigned int v1 = 0; unsigned int v2 = 0;

    _mesh.clear();
    // Temporarily disable bottom-up incidences
    // since it's way faster to first add all the
    // geometry and compute them in one pass afterwards
    _mesh.enable_bottom_up_incidences(false);

    /*
     * Header
     */

    bool header_found = true;

    // Get first line
    getCleanLine(iff, line);
    sstr.str(line);

    // Check header
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "OVM") {
        //iff.close();
        header_found = false;
        std::cerr << "The specified file might not be in OpenVolumeMesh format!" << std::endl;
        //return false;
    }

    // Get ASCII/BINARY string
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp == "BINARY") {
        iff.close();
        std::cerr << "Binary files are not supported at the moment!" << std::endl;
        return false;
    }

    /*
     * Vertices
     */
    if(!header_found) {
        sstr.clear();
        sstr.str(line);
    } else {
        getCleanLine(iff, line);
        sstr.clear();
        sstr.str(line);
    }

    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "VERTICES") {
        iff.close();
        std::cerr << "No vertex section defined!" << std::endl;
        return false;
    } else {

        // Read in number of vertices
        getCleanLine(iff, line);
        sstr.clear();
        sstr.str(line);
        sstr >> c;

        // Read in vertices
        for(unsigned int i = 0u; i < c; ++i) {

            getCleanLine(iff, line);
            sstr.clear();
            sstr.str(line);
            sstr >> v[0];
            sstr >> v[1];
            sstr >> v[2];
            _mesh.add_vertex(v);
        }
    }

    /*
     * Edges
     */
    getCleanLine(iff, line);
    sstr.clear();
    sstr.str(line);
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "EDGES") {
        iff.close();
        std::cerr << "No edge section defined!" << std::endl;
        return false;
    } else {

        // Read in number of edges
        getCleanLine(iff, line);
        sstr.clear();
        sstr.str(line);
        sstr >> c;

        // Read in edges
        for(unsigned int i = 0u; i < c; ++i) {

            getCleanLine(iff, line);
            sstr.clear();
            sstr.str(line);
            sstr >> v1;
            sstr >> v2;
            _mesh.add_edge(VertexHandle(v1), VertexHandle(v2), true);
        }
    }

    /*
     * Faces
     */
    getCleanLine(iff, line);
    sstr.clear();
    sstr.str(line);
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "FACES") {
        iff.close();
        std::cerr << "No face section defined!" << std::endl;
        return false;
    } else {

        // Read in number of faces
        getCleanLine(iff, line);
        sstr.clear();
        sstr.str(line);
        sstr >> c;

        // Read in faces
        for(unsigned int i = 0u; i < c; ++i) {

            getCleanLine(iff, line);
            sstr.clear();
            sstr.str(line);

            std::vector<HalfEdgeHandle> hes;

            // Get face valence
            unsigned int val = 0u;
            sstr >> val;

            // Read half-edge indices
            for(unsigned int e = 0; e < val; ++e) {
                sstr >> v1;
                hes.push_back(HalfEdgeHandle(v1));
            }

            _mesh.add_face(hes, _topologyCheck);
        }
    }

    /*
     * Cells
     */
    getCleanLine(iff, line);
    sstr.clear();
    sstr.str(line);
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "POLYHEDRA") {
        iff.close();
        std::cerr << "No polyhedra section defined!" << std::endl;
        return false;
    } else {

        // Read in number of cells
        getCleanLine(iff, line);
        sstr.clear();
        sstr.str(line);
        sstr >> c;

        // Read in cells
        for(unsigned int i = 0u; i < c; ++i) {

            getCleanLine(iff, line);
            sstr.clear();
            sstr.str(line);

            std::vector<HalfFaceHandle> hfs;

            // Get cell valence
            unsigned int val = 0u;
            sstr >> val;

            // Read half-face indices
            for(unsigned int f = 0; f < val; ++f) {
                sstr >> v1;
                hfs.push_back(HalfFaceHandle(v1));
            }

            _mesh.add_cell(hfs, _topologyCheck);
        }
    }

    while(!iff.eof()) {
        // "End of file reached while searching for input!"
        // is thrown here. \TODO Fix it!

        // Read property
        readProperty(iff, _mesh);
    }

    iff.close();

    if(_computeBottomUpIncidences) {
        // Compute bottom-up incidences
        _mesh.enable_bottom_up_incidences(true);
    }

    std::cerr << "######## openvolumemesh info #########" << std::endl;
    std::cerr << "#vertices: " << _mesh.n_vertices() << std::endl;
    std::cerr << "#edges:    " << _mesh.n_edges() << std::endl;
    std::cerr << "#faces:    " << _mesh.n_faces() << std::endl;
    std::cerr << "#cells:    " << _mesh.n_cells() << std::endl;
    std::cerr << "######################################" << std::endl;

    return true;
}

//==================================================

template <class MeshT>
void FileManager::readProperty(std::istream& _iff, MeshT& _mesh) const {

    std::string line, entity_t, prop_t, name;
    std::stringstream sstr;

    getCleanLine(_iff, line);

    if(line.empty()) return;

    sstr.clear();
    sstr.str(line);
    sstr >> entity_t;
    std::transform(entity_t.begin(), entity_t.end(), entity_t.begin(), ::tolower);
    sstr >> prop_t;
    std::transform(prop_t.begin(), prop_t.end(), prop_t.begin(), ::tolower);
    name = line;
    extractQuotedText(name);

    if(prop_t == typeName<int>()) generateGenericProperty<int, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<unsigned int>()) generateGenericProperty<unsigned int, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<short>()) generateGenericProperty<short, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<long>()) generateGenericProperty<long, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<unsigned long>()) generateGenericProperty<unsigned long, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<char>()) generateGenericProperty<char, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<unsigned char>()) generateGenericProperty<unsigned char, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<bool>()) generateGenericProperty<bool, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<float>()) generateGenericProperty<float, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<double>()) generateGenericProperty<double, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<std::string>()) generateGenericProperty<std::string, MeshT>(entity_t, name, _iff, _mesh);
}

//==================================================

template <class PropT, class MeshT>
void FileManager::generateGenericProperty(const std::string& _entity_t, const std::string& _name,
                                          std::istream& _iff, MeshT& _mesh) const {

    if(_entity_t == "vprop") {
        VertexPropertyT<PropT> prop = _mesh.template request_vertex_property<PropT>(_name);
        prop.deserialize(_iff);
        _mesh.set_persistent(prop);
    } else if(_entity_t == "eprop") {
        EdgePropertyT<PropT> prop = _mesh.template request_edge_property<PropT>(_name);
        prop.deserialize(_iff);
        _mesh.set_persistent(prop);
    } else if(_entity_t == "heprop") {
        HalfEdgePropertyT<PropT> prop = _mesh.template request_halfedge_property<PropT>(_name);
        prop.deserialize(_iff);
        _mesh.set_persistent(prop);
    } else if(_entity_t == "fprop") {
        FacePropertyT<PropT> prop = _mesh.template request_face_property<PropT>(_name);
        prop.deserialize(_iff);
        _mesh.set_persistent(prop);
    } else if(_entity_t == "hfprop") {
        HalfFacePropertyT<PropT> prop = _mesh.template request_halfface_property<PropT>(_name);
        prop.deserialize(_iff);
        _mesh.set_persistent(prop);
    } else if(_entity_t == "cprop") {
        CellPropertyT<PropT> prop = _mesh.template request_cell_property<PropT>(_name);
        prop.deserialize(_iff);
        _mesh.set_persistent(prop);
    } else if(_entity_t == "mprop") {
        MeshPropertyT<PropT> prop = _mesh.template request_mesh_property<PropT>(_name);
        prop.deserialize(_iff);
        _mesh.set_persistent(prop);
    }
}

//==================================================

template<class MeshT>
bool FileManager::writeFile(const std::string& _filename, const MeshT& _mesh) const {

    typedef typename MeshT::Face Face;
    std::ofstream off(_filename.c_str(), std::ios::out);

    if(!off.good()) {
        std::cerr << "Error: Could not open file " << _filename << " for writing!" << std::endl;
        off.close();
        return false;
    }

    // Write header
    off << "OVM ASCII" << std::endl;

    unsigned int n_vertices(_mesh.n_vertices());
    off << "Vertices" << std::endl;
    off << n_vertices << std::endl;

    typedef typename MeshT::PointT Point;

    // write vertices
    for(VertexIter v_it = _mesh.v_iter(); v_it; ++v_it) {

        Point v = _mesh.vertex(*v_it);
        off << v[0] << " " << v[1] << " " << v[2] << std::endl;
    }

    unsigned int n_edges(_mesh.n_edges());
    off << "Edges" << std::endl;
    off << n_edges << std::endl;

    // write edges
    for(EdgeIter e_it = _mesh.e_iter(); e_it; ++e_it) {

        VertexHandle from_vertex = _mesh.edge(*e_it).from_vertex();
        VertexHandle to_vertex = _mesh.edge(*e_it).to_vertex();
        off << from_vertex << " " << to_vertex << std::endl;
    }

    unsigned int n_faces(_mesh.n_faces());
    off << "Faces" << std::endl;
    off << n_faces << std::endl;

    // write faces
    for(FaceIter f_it = _mesh.f_iter(); f_it; ++f_it) {

        off << _mesh.face(*f_it).halfedges().size() << " ";

        std::vector<HalfEdgeHandle> halfedges = _mesh.face(*f_it).halfedges();

        for(typename std::vector<HalfEdgeHandle>::const_iterator it = halfedges.begin(); it
                != halfedges.end(); ++it) {

            off << it->idx();

            if((it + 1) != halfedges.end())
                off << " ";
        }

        off << std::endl;
    }

    unsigned int n_cells(_mesh.n_cells());
    off << "Polyhedra" << std::endl;
    off << n_cells << std::endl;

    for(CellIter c_it = _mesh.c_iter(); c_it; ++c_it) {

        off << _mesh.cell(*c_it).halffaces().size() << " ";

        std::vector<HalfFaceHandle> halffaces = _mesh.cell(*c_it).halffaces();

        for(typename std::vector<HalfFaceHandle>::const_iterator it = halffaces.begin(); it
                != halffaces.end(); ++it) {

            off << it->idx();

            if((it + 1) != halffaces.end())
                off << " ";
        }

        off << std::endl;
    }

    // write vertex props
    writeProps(off, _mesh.vertex_props_begin(), _mesh.vertex_props_end());
    // write edge props
    writeProps(off, _mesh.edge_props_begin(), _mesh.edge_props_end());
    // write halfedge props
    writeProps(off, _mesh.halfedge_props_begin(), _mesh.halfedge_props_end());
    // write face props
    writeProps(off, _mesh.face_props_begin(), _mesh.face_props_end());
    // write halfface props
    writeProps(off, _mesh.halfface_props_begin(), _mesh.halfface_props_end());
    // write cell props
    writeProps(off, _mesh.cell_props_begin(), _mesh.cell_props_end());
    // write mesh props
    writeProps(off, _mesh.mesh_props_begin(), _mesh.mesh_props_end());

    off.close();

    return true;
}

//==================================================

template<class IteratorT>
void FileManager::writeProps(std::ostream& _ostr, const IteratorT& _begin, const IteratorT& _end) const {

    // write props
    for(IteratorT p_it = _begin;
            p_it != _end; ++p_it) {
        if(!(*p_it)->persistent()) continue;
        if((*p_it)->anonymous()) {
            std::cerr << "Serialization of anonymous properties is not supported!" << std::endl;
            continue;
        }
        (*p_it)->serialize(_ostr);
    }
}

//==================================================

} // Namespace IO

} // Namespace FileManager
