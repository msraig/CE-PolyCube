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

#ifndef PROPERTYHANDLES_HH_
#define PROPERTYHANDLES_HH_

#include "OpenVolumeMeshHandle.hh"

namespace OpenVolumeMesh {

// Defines for property handles
class VertexPropHandle      : public OpenVolumeMeshHandle { public: VertexPropHandle(int _idx = -1)     : OpenVolumeMeshHandle(_idx) {} };
class EdgePropHandle        : public OpenVolumeMeshHandle { public: EdgePropHandle(int _idx = -1)       : OpenVolumeMeshHandle(_idx) {} };
class HalfEdgePropHandle    : public OpenVolumeMeshHandle { public: HalfEdgePropHandle(int _idx = -1)   : OpenVolumeMeshHandle(_idx) {} };
class FacePropHandle        : public OpenVolumeMeshHandle { public: FacePropHandle(int _idx = -1)       : OpenVolumeMeshHandle(_idx) {} };
class HalfFacePropHandle    : public OpenVolumeMeshHandle { public: HalfFacePropHandle(int _idx = -1)   : OpenVolumeMeshHandle(_idx) {} };
class CellPropHandle        : public OpenVolumeMeshHandle { public: CellPropHandle(int _idx = -1)       : OpenVolumeMeshHandle(_idx) {} };
class MeshPropHandle        : public OpenVolumeMeshHandle { public: MeshPropHandle(int _idx = -1)       : OpenVolumeMeshHandle(_idx) {} };

} // Namespace OpenVolumeMesh

#endif /* PROPERTYHANDLES_HH_ */
