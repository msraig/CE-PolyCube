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

#include <iostream>

#include "BaseEntities.hh"

namespace OpenVolumeMesh {

std::ostream& operator<<(std::ostream& _os, const OpenVolumeMeshEdge& _edge) {
    return _os << "(" << _edge.from_vertex() << ", " << _edge.to_vertex() << ")";
}

std::ostream& operator<<(std::ostream& _os, const OpenVolumeMeshFace& _face) {
    _os << "(";
    for(std::vector<HalfEdgeHandle>::const_iterator it =
            _face.halfedges().begin(); it < _face.halfedges().end(); ++it) {
        _os << *it;
        if(it + 1 < _face.halfedges().end())
            _os << ", ";
    }
    _os << ")";
    return _os;
}

std::ostream& operator<<(std::ostream& _os, const OpenVolumeMeshCell& _cell) {
    _os << "(";
    for(std::vector<HalfFaceHandle>::const_iterator it =
            _cell.halffaces().begin(); it < _cell.halffaces().end(); ++it) {
        _os << *it;
        if(it + 1 < _cell.halffaces().end())
            _os << ", ";
    }
    _os << ")";
    return _os;
}

} // Namespace OpenVolumeMesh
