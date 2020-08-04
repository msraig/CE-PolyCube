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

#ifndef BASEENTITIES_HH_
#define BASEENTITIES_HH_

#include <vector>
#include <cassert>
#include "OpenVolumeMeshHandle.hh"

namespace OpenVolumeMesh {

class OpenVolumeMeshEdge {
friend class TopologyKernel;
public:
    OpenVolumeMeshEdge(const VertexHandle& _fromVertex,
                       const VertexHandle& _toVertex) :
        fromVertex_(_fromVertex),
        toVertex_(_toVertex) {
    }

    virtual ~OpenVolumeMeshEdge() {
    }

    const VertexHandle& from_vertex() const {
        return fromVertex_;
    }
    const VertexHandle& to_vertex() const {
        return toVertex_;
    }

protected:

    void set_from_vertex(const VertexHandle& _vertex) {
        fromVertex_ = _vertex;
    }
    void set_to_vertex(const VertexHandle& _vertex) {
        toVertex_ = _vertex;
    }

private:
    VertexHandle fromVertex_;
    VertexHandle toVertex_;
};

// Stream operator for edges
std::ostream& operator<<(std::ostream& _os, const OpenVolumeMeshEdge& _edge);

//***************************************************************************

class OpenVolumeMeshFace {
friend class TopologyKernel;
public:
    OpenVolumeMeshFace(const std::vector<HalfEdgeHandle>& _halfedges) :
        halfedges_(_halfedges) {
    }

    virtual ~OpenVolumeMeshFace() {
    }

    const std::vector<HalfEdgeHandle>& halfedges() const {
        return halfedges_;
    }

	const std::vector<VertexHandle>& vertices() const {
		return vertices_;
	}

	std::vector<HalfEdgeHandle>& get_halfedges()
	{
		return halfedges_;
	}

	const unsigned get_halfedges_size() const { return halfedges_.size(); }
	void get_halfedges_id(unsigned i, HalfEdgeHandle& he) const
	{
		assert( i>=0 && i < halfedges_.size());
		he = halfedges_[i];
	}

protected:

    void set_halfedges(const std::vector<HalfEdgeHandle>& _halfedges)
	{
        halfedges_ = _halfedges;
    }

	void set_vertices(const std::vector<VertexHandle>& _vertices) 
	{
		vertices_ = _vertices;
	}

private:
    std::vector<HalfEdgeHandle> halfedges_;
	std::vector<VertexHandle> vertices_;
};

// Stream operator for faces
std::ostream& operator<<(std::ostream& _os, const OpenVolumeMeshFace& _face);

//***************************************************************************

class OpenVolumeMeshCell {
friend class TopologyKernel;
public:
    OpenVolumeMeshCell(const std::vector<HalfFaceHandle>& _halffaces) :
        halffaces_(_halffaces) {
    }

    virtual ~OpenVolumeMeshCell() {
    }

    const std::vector<HalfFaceHandle>& halffaces() const {
        return halffaces_;
    }

	std::vector<HalfFaceHandle>& get_halffaces()
	{
		return halffaces_;
	}

	void update_cell_halffacs(int hfh_id0, int hfh_id1 , int hf_id_th)
	{
		unsigned hes_size = halffaces_.size(); unsigned ii = 0;
		int it_id = 0;
		
		while( ii < hes_size )
		{
			it_id = halffaces_[ii].idx();
			if( it_id == hfh_id0 || it_id == hfh_id1 )
			{
				halffaces_.erase( halffaces_.begin() + ii );
				--hes_size;
			}
			else
			{
				if(it_id > hf_id_th)
				{ halffaces_[ii].idx( it_id - 2 );}
			}
			++ii;
		}
	}

	const std::vector<VertexHandle>& vertices() const
	{
		return vertices_;
	}
	void set_vertices(const std::vector<VertexHandle>& _vertices)
	{
		vertices_ = _vertices;
	}

protected:

    void set_halffaces(const std::vector<HalfFaceHandle>& _halffaces) {
        halffaces_ = _halffaces;
    }

private:
	std::vector<HalfFaceHandle> halffaces_;
	std::vector<VertexHandle> vertices_;
};

// Stream operator for cells
std::ostream& operator<<(std::ostream& _os, const OpenVolumeMeshCell& _cell);

} // Namespace OpenVolumeMesh

#endif /* BASEENTITIES_HH_ */
