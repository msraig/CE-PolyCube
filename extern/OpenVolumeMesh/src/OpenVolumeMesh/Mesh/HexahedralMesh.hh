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

#ifndef HEXAHEDRALMESH_HH_
#define HEXAHEDRALMESH_HH_

#include "HexahedralMeshTopologyKernel.hh"
#include "../Core/GeometryKernel.hh"

namespace OpenVolumeMesh {

/*
 * Predefines for most common mesh types
 */
typedef GeometryKernel<Geometry::Vec2i, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV2i;
typedef GeometryKernel<Geometry::Vec2ui, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV2ui;
typedef GeometryKernel<Geometry::Vec2f, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV2f;
typedef GeometryKernel<Geometry::Vec2d, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV2d;
typedef GeometryKernel<Geometry::Vec2c, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV2c;
typedef GeometryKernel<Geometry::Vec2uc, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV2uc;
typedef GeometryKernel<Geometry::Vec3i, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV3i;
typedef GeometryKernel<Geometry::Vec3ui, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV3ui;
typedef GeometryKernel<Geometry::Vec3f, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV3f;
typedef GeometryKernel<Geometry::Vec3d, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV3d;
typedef GeometryKernel<Geometry::Vec3c, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV3c;
typedef GeometryKernel<Geometry::Vec3uc, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV3uc;
typedef GeometryKernel<Geometry::Vec4i, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV4i;
typedef GeometryKernel<Geometry::Vec4ui, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV4ui;
typedef GeometryKernel<Geometry::Vec4f, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV4f;
typedef GeometryKernel<Geometry::Vec4d, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV4d;
typedef GeometryKernel<Geometry::Vec4c, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV4c;
typedef GeometryKernel<Geometry::Vec4uc, HexahedralMeshTopologyKernel> GeometricHexahedralMeshV4uc;

typedef HexahedralMeshTopologyKernel TopologicHexahedralMesh;

} // Namespace OpenVolumeMesh

#endif /* HEXAHEDRALMESH_HH_ */
