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

#define PROPERTYDEFINEST_CC

#include "PropertyDefines.hh"
#include "PropertyPtr.hh"

namespace OpenVolumeMesh {

/// Property classes for the different entity types
template<class T>
VertexPropertyT<T>::VertexPropertyT(const std::string& _name, ResourceManager& _resMan, VertexPropHandle _handle, const T _def) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, VertexPropHandle>(new OpenVolumeMeshPropertyT<T>(_name, _def), _resMan, _handle) {

}

template<class T>
void VertexPropertyT<T>::serialize(std::ostream& _ostr) const {
    _ostr << "VProp ";
    _ostr << typeName<T>() << " ";
    PropertyPtr<OpenVolumeMeshPropertyT<T>, VertexPropHandle>::get()->serialize(_ostr);
}

template<class T>
void VertexPropertyT<T>::deserialize(std::istream& _istr) {
    PropertyPtr<OpenVolumeMeshPropertyT<T>, VertexPropHandle>::get()->deserialize(_istr);
}

template<class T>
EdgePropertyT<T>::EdgePropertyT(const std::string& _name, ResourceManager& _resMan, EdgePropHandle _handle, const T _def) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, EdgePropHandle>(new OpenVolumeMeshPropertyT<T>(_name, _def), _resMan, _handle) {

}

template<class T>
void EdgePropertyT<T>::serialize(std::ostream& _ostr) const {
    _ostr << "EProp ";
    _ostr << typeName<T>() << " ";
    PropertyPtr<OpenVolumeMeshPropertyT<T>, EdgePropHandle>::get()->serialize(_ostr);
}

template<class T>
void EdgePropertyT<T>::deserialize(std::istream& _istr) {
    PropertyPtr<OpenVolumeMeshPropertyT<T>, EdgePropHandle>::get()->deserialize(_istr);
}

template<class T>
HalfEdgePropertyT<T>::HalfEdgePropertyT(const std::string& _name, ResourceManager& _resMan, HalfEdgePropHandle _handle, const T _def) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, HalfEdgePropHandle>(new OpenVolumeMeshPropertyT<T>(_name, _def), _resMan, _handle) {

}

template<class T>
void HalfEdgePropertyT<T>::serialize(std::ostream& _ostr) const {
    _ostr << "HEProp ";
    _ostr << typeName<T>() << " ";
    PropertyPtr<OpenVolumeMeshPropertyT<T>, HalfEdgePropHandle>::get()->serialize(_ostr);
}

template<class T>
void HalfEdgePropertyT<T>::deserialize(std::istream& _istr) {
    PropertyPtr<OpenVolumeMeshPropertyT<T>, HalfEdgePropHandle>::get()->deserialize(_istr);
}

template<class T>
FacePropertyT<T>::FacePropertyT(const std::string& _name, ResourceManager& _resMan, FacePropHandle _handle, const T _def) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, FacePropHandle>(new OpenVolumeMeshPropertyT<T>(_name, _def), _resMan, _handle) {

}

template<class T>
void FacePropertyT<T>::serialize(std::ostream& _ostr) const {
    _ostr << "FProp ";
    _ostr << typeName<T>() << " ";
    PropertyPtr<OpenVolumeMeshPropertyT<T>, FacePropHandle>::get()->serialize(_ostr);
}

template<class T>
void FacePropertyT<T>::deserialize(std::istream& _istr) {
    PropertyPtr<OpenVolumeMeshPropertyT<T>, FacePropHandle>::get()->deserialize(_istr);
}

template<class T>
HalfFacePropertyT<T>::HalfFacePropertyT(const std::string& _name, ResourceManager& _resMan, HalfFacePropHandle _handle, const T _def) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, HalfFacePropHandle>(new OpenVolumeMeshPropertyT<T>(_name, _def), _resMan, _handle) {

}

template<class T>
void HalfFacePropertyT<T>::serialize(std::ostream& _ostr) const {
    _ostr << "HFProp ";
    _ostr << typeName<T>() << " ";
    PropertyPtr<OpenVolumeMeshPropertyT<T>, HalfFacePropHandle>::get()->serialize(_ostr);
}

template<class T>
void HalfFacePropertyT<T>::deserialize(std::istream& _istr) {
    PropertyPtr<OpenVolumeMeshPropertyT<T>, HalfFacePropHandle>::get()->deserialize(_istr);
}

template<class T>
CellPropertyT<T>::CellPropertyT(const std::string& _name, ResourceManager& _resMan, CellPropHandle _handle, const T _def) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, CellPropHandle>(new OpenVolumeMeshPropertyT<T>(_name, _def), _resMan, _handle) {

}

template<class T>
void CellPropertyT<T>::serialize(std::ostream& _ostr) const {
    _ostr << "CProp ";
    _ostr << typeName<T>() << " ";
    PropertyPtr<OpenVolumeMeshPropertyT<T>, CellPropHandle>::get()->serialize(_ostr);
}

template<class T>
void CellPropertyT<T>::deserialize(std::istream& _istr) {
    PropertyPtr<OpenVolumeMeshPropertyT<T>, CellPropHandle>::get()->deserialize(_istr);
}

template<class T>
MeshPropertyT<T>::MeshPropertyT(const std::string& _name, ResourceManager& _resMan, MeshPropHandle _handle, const T _def) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, MeshPropHandle>(new OpenVolumeMeshPropertyT<T>(_name, _def), _resMan, _handle) {

}

template<class T>
void MeshPropertyT<T>::serialize(std::ostream& _ostr) const {
    _ostr << "MProp ";
    _ostr << typeName<T>() << " ";
    PropertyPtr<OpenVolumeMeshPropertyT<T>, MeshPropHandle>::get()->serialize(_ostr);
}

template<class T>
void MeshPropertyT<T>::deserialize(std::istream& _istr) {
    PropertyPtr<OpenVolumeMeshPropertyT<T>, MeshPropHandle>::get()->deserialize(_istr);
}

template <class T>
const std::string typeName() {
    throw std::runtime_error("Serialization is not supported for these data types!");
}

} // Namespace OpenVolumeMesh
