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

#ifndef GEOMETRYKERNEL_HH_
#define GEOMETRYKERNEL_HH_

#include <cassert>

#include "../Geometry/VectorT.hh"
#include "TopologyKernel.hh"

namespace OpenVolumeMesh {

template <class VecT, class TopologyKernelT = TopologyKernel>
class GeometryKernel : public TopologyKernelT {
public:

    typedef VecT PointT;
    typedef TopologyKernelT KernelT;

    /// Constructor
    GeometryKernel() {}

    /// Destructor
    ~GeometryKernel() {}

    /// Override of empty add_vertex function
    virtual VertexHandle add_vertex() { return add_vertex(VecT()); }

    /// Add a geometric point to the mesh
    VertexHandle add_vertex(const VecT& _p) {

        // Store vertex in list
        vertices_.push_back(_p);

        // Get handle of recently created vertex
        return KernelT::add_vertex();
    }

    /// Set the coordinates of point _vh
    void set_vertex(const VertexHandle& _vh, const VecT& _p) {

        assert(_vh.idx() < (int)vertices_.size());

        vertices_[_vh.idx()] = _p;
    }

    /// Get point _vh's coordinates
    const VecT& vertex(const VertexHandle& _vh) const {
        return vertices_[_vh.idx()];
    }

    virtual VertexIter delete_vertex(const VertexHandle& _h) {
        assert(_h.idx() < (int)TopologyKernel::n_vertices());

        VertexIter nV = TopologyKernelT::delete_vertex(_h);

        vertices_.erase(vertices_.begin() + _h.idx());

        return nV;
    }

protected:

    virtual void delete_multiple_vertices(const std::vector<bool>& _tag) {

        assert(_tag.size() == TopologyKernelT::n_vertices());

        std::vector<VecT> newVertices;

        typename std::vector<VecT>::const_iterator v_it = vertices_.begin();

        for(std::vector<bool>::const_iterator t_it = _tag.begin(),
                t_end = _tag.end(); t_it != t_end; ++t_it, ++v_it) {

            if(!(*t_it)) {
                // Not marked as deleted

                newVertices.push_back(*v_it);
            }
        }

        // Swap vertices
        vertices_.swap(newVertices);

        TopologyKernelT::delete_multiple_vertices(_tag);
    }

public:

    virtual void clear(bool _clearProps = true) {

        vertices_.clear();
        TopologyKernelT::clear(_clearProps);
    }

    typename PointT::value_type length(const EdgeHandle& _eh) const {

        const typename TopologyKernelT::Edge& e = TopologyKernelT::edge(_eh);
        return (vertex(e.to_vertex()) - vertex(e.from_vertex())).length();
    }

    PointT vector(const EdgeHandle& _eh) const {

        const typename TopologyKernelT::Edge& e = TopologyKernelT::edge(_eh);
        return (vertex(e.to_vertex()) - vertex(e.from_vertex()));
    }

    PointT barycenter(const EdgeHandle& _eh) const {
        return PointT(0.5 * vertex(TopologyKernelT::edge(_eh).from_vertex()) +
                      0.5 * vertex(TopologyKernelT::edge(_eh).to_vertex()));
    }

    PointT barycenter(const FaceHandle& _fh) const {
        PointT p(typename PointT::value_type(0));
        typename PointT::value_type valence = 0;
        HalfFaceVertexIter hfv_it =
                TopologyKernelT::hfv_iter(TopologyKernelT::halfface_handle(_fh, 0));
        for(; hfv_it.valid(); ++hfv_it, valence += 1) {
            p += vertex(*hfv_it);
        }
        p /= valence;
        return p;
    }

    PointT barycenter(const CellHandle& _ch) const {
        PointT p(typename PointT::value_type(0));
        typename PointT::value_type valence = 0;
        CellVertexIter cv_it = TopologyKernelT::cv_iter(_ch);
        for(; cv_it.valid(); ++cv_it, valence += 1) {
            p += vertex(*cv_it);
        }
        p /= valence;
        return p;
    }

    void clone_vertices(std::vector<VecT>& _copy) const {
        _copy.clear();
        _copy.reserve(vertices_.size());
        std::copy(vertices_.begin(), vertices_.end(), std::back_inserter(_copy));
    }

    void swap_vertices(std::vector<VecT>& _copy) {
        if(_copy.size() != vertices_.size()) {
            std::cerr << "Vertex vectors differ in size!" << std::endl;
            return;
        }
        std::swap(vertices_, _copy);
    }

private:

    std::vector<VecT> vertices_;
};

} // Namespace OpenVolumeMesh

#endif /* GEOMETRYKERNEL_HH_ */
