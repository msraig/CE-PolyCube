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

#ifndef TOPOLOGYKERNEL_HH_
#define TOPOLOGYKERNEL_HH_

#include <set>
#include <vector>
#include <iostream>
#include <cassert>

#include "BaseEntities.hh"
#include "OpenVolumeMeshHandle.hh"
#include "ResourceManager.hh"
#include "Iterators.hh"

namespace OpenVolumeMesh {

class TopologyKernel : public ResourceManager {
public:

    TopologyKernel();
    virtual ~TopologyKernel();

    /*
     * Defines and constants
     */

    static const VertexHandle   InvalidVertexHandle;
    static const EdgeHandle     InvalidEdgeHandle;
    static const FaceHandle     InvalidFaceHandle;
    static const CellHandle     InvalidCellHandle;
    static const HalfEdgeHandle InvalidHalfEdgeHandle;
    static const HalfFaceHandle InvalidHalfFaceHandle;

    typedef OpenVolumeMeshEdge Edge;
    typedef OpenVolumeMeshFace Face;
    typedef OpenVolumeMeshCell Cell;

    // Add StatusAttrib to list of friend classes
    // since it provides a garbage collection
    // that needs access to some protected methods
    friend class StatusAttrib;

    //=====================================================================
    // Iterators
    //=====================================================================

    friend class VertexOHalfEdgeIter;
    friend class HalfEdgeHalfFaceIter;
    friend class VertexCellIter;
    friend class HalfEdgeCellIter;
    friend class CellVertexIter;
    friend class CellCellIter;
    friend class HalfFaceVertexIter;
    friend class BoundaryHalfFaceHalfFaceIter;
    friend class BoundaryFaceIter;
    friend class VertexIter;
    friend class EdgeIter;
    friend class HalfEdgeIter;
    friend class FaceIter;
    friend class HalfFaceIter;
    friend class CellIter;

    VertexOHalfEdgeIter voh_iter(const VertexHandle& _h) const {
        return VertexOHalfEdgeIter(_h, this);
    }

    HalfEdgeHalfFaceIter hehf_iter(const HalfEdgeHandle& _h) const {
        return HalfEdgeHalfFaceIter(_h, this);
    }

    VertexCellIter vc_iter(const VertexHandle& _h) const {
        return VertexCellIter(_h, this);
    }

    HalfEdgeCellIter hec_iter(const HalfEdgeHandle& _h) const {
        return HalfEdgeCellIter(_h, this);
    }

    CellVertexIter cv_iter(const CellHandle& _h) const {
        return CellVertexIter(_h, this);
    }

    CellCellIter cc_iter(const CellHandle& _h) const {
        return CellCellIter(_h, this);
    }

    HalfFaceVertexIter hfv_iter(const HalfFaceHandle& _h) const {
	// HalfFaceVertexIter hfv_iter(const HalfFaceHandle& _h) {
		//hfv_it_.construct_iter(_h, this);
		//return hfv_it_;
        return HalfFaceVertexIter(_h, this);
    }

    BoundaryHalfFaceHalfFaceIter bhfhf_iter(const HalfFaceHandle& _ref_h) const {
        return BoundaryHalfFaceHalfFaceIter(_ref_h, this);
    }

    BoundaryFaceIter bf_iter() const {
        return BoundaryFaceIter(this);
    }

    VertexIter v_iter() const {
        return VertexIter(this);
    }

    VertexIter vertices_begin() const {
        return VertexIter(this, VertexHandle(0));
    }

    VertexIter vertices_end() const {
        return VertexIter(this, VertexHandle(n_vertices()));
    }

    EdgeIter e_iter() const {
        return EdgeIter(this);
    }

    EdgeIter edges_begin() const {
        return EdgeIter(this, EdgeHandle(0));
    }

    EdgeIter edges_end() const {
        return EdgeIter(this, EdgeHandle(edges_.size()));
    }

    HalfEdgeIter he_iter() const {
        return HalfEdgeIter(this);
    }

    HalfEdgeIter halfedges_begin() const {
        return HalfEdgeIter(this, HalfEdgeHandle(0));
    }

    HalfEdgeIter halfedges_end() const {
        return HalfEdgeIter(this, HalfEdgeHandle(edges_.size() * 2));
    }

    FaceIter f_iter() const {
        return FaceIter(this);
    }

    FaceIter faces_begin() const {
        return FaceIter(this, FaceHandle(0));
    }

    FaceIter faces_end() const {
        return FaceIter(this, FaceHandle(faces_.size()));
    }

    HalfFaceIter hf_iter() const {
        return HalfFaceIter(this);
    }

    HalfFaceIter halffaces_begin() const {
        return HalfFaceIter(this, HalfFaceHandle(0));
    }

    HalfFaceIter halffaces_end() const {
        return HalfFaceIter(this, HalfFaceHandle(faces_.size() * 2));
    }

    CellIter c_iter() const {
        return CellIter(this);
    }

    CellIter cells_begin() const {
        return CellIter(this, CellHandle(0));
    }

    CellIter cells_end() const {
        return CellIter(this, CellHandle(cells_.size()));
    }

    /*
     * Virtual functions with implementation
     */

    /// Get number of vertices in mesh
    virtual unsigned int n_vertices()   const { return n_vertices_; }
    /// Get number of edges in mesh
    virtual unsigned int n_edges()      const { return edges_.size(); }
    /// Get number of halfedges in mesh
    virtual unsigned int n_halfedges()  const { return edges_.size() * 2u; }
    /// Get number of faces in mesh
    virtual unsigned int n_faces()      const { return faces_.size(); }
    /// Get number of halffaces in mesh
    virtual unsigned int n_halffaces()  const { return faces_.size() * 2u; }
    /// Get number of cells in mesh
    virtual unsigned int n_cells()      const { return cells_.size(); }

protected:
	HalfFaceVertexIter hfv_it_;

private:

    // Cache total vertex number
    unsigned int n_vertices_;

public:

    /// Add abstract vertex
    virtual VertexHandle add_vertex();

    //=======================================================================

    /// Add edge
    virtual EdgeHandle add_edge(const VertexHandle& _fromVertex, const VertexHandle& _toHandle, bool _allowDuplicates = false);

    /// Add face via incident edges
    virtual FaceHandle add_face(const std::vector<HalfEdgeHandle>& _halfedges, bool _topologyCheck = false);

    /// Add face via incident vertices
    virtual FaceHandle add_face(const std::vector<VertexHandle>& _vertices);

    /// Add cell via incident halffaces
    virtual CellHandle add_cell(const std::vector<HalfFaceHandle>& _halffaces, bool _topologyCheck = false);

    /// Set the vertices of an edge
    void set_edge(const EdgeHandle& _eh, const VertexHandle& _fromVertex, const VertexHandle& _toVertex);

    /// Set the half-edges of a face
    void set_face(const FaceHandle& _fh, const std::vector<HalfEdgeHandle>& _hes);

    /// Set the half-faces of a cell
    void set_cell(const CellHandle& _ch, const std::vector<HalfFaceHandle>& _hfs);

    /*
     * Non-virtual functions
     */

    /// Get edge with handle _edgeHandle
    const Edge& edge(const EdgeHandle& _edgeHandle) const;

    /// Get face with handle _faceHandle
    const Face& face(const FaceHandle& _faceHandle) const;

    /// Get cell with handle _cellHandle
    const Cell& cell(const CellHandle& _cellHandle) const;

    /// Get edge with handle _edgeHandle
    Edge& edge(const EdgeHandle& _edgeHandle);

    /// Get face with handle _faceHandle
    Face& face(const FaceHandle& _faceHandle);

    /// Get cell with handle _cellHandle
    Cell& cell(const CellHandle& _cellHandle);

    /// Get edge that corresponds to halfedge with handle _halfEdgeHandle
    const Edge halfedge(const HalfEdgeHandle& _halfEdgeHandle) const;

    /// Get face that corresponds to halfface with handle _halfFaceHandle
    const Face halfface(const HalfFaceHandle& _halfFaceHandle) const;

	//added by Xiaoming Fu
	void get_halfedges_from_halfface(const HalfFaceHandle& _halfFaceHandle, std::vector<HalfEdgeHandle>& hes) const;
	void get_vertices_from_halfface(const HalfFaceHandle& _halfFaceHandle, std::vector<VertexHandle>& vs) const;
	void get_vertices_from_halfface_n(const HalfFaceHandle& _halfFaceHandle, std::vector<VertexHandle>& vs, const int& n) const;
	void get_halffaces_from_cell(const CellHandle& _cellHandle, std::vector<HalfFaceHandle>& hfs) const;

    /// Get opposite halfedge that corresponds to halfedge with handle _halfEdgeHandle
    const Edge opposite_halfedge(const HalfEdgeHandle& _halfEdgeHandle) const;

    /// Get opposite halfface that corresponds to halfface with handle _halfFaceHandle
    const Face opposite_halfface(const HalfFaceHandle& _halfFaceHandle) const;

    // Get halfedge from vertex _vh1 to _vh2
    const HalfEdgeHandle halfedge(const VertexHandle& _vh1, const VertexHandle& _vh2) const;

    // Get half-face from list of incident vertices (in connected order)
    // Note: Only the first three vertices are checked
    const HalfFaceHandle halfface(const std::vector<VertexHandle>& _vs) const;

    // Get half-face from list of incident half-edges
    // Note: Only the first two half-edges are checked
    const HalfFaceHandle halfface(const std::vector<HalfEdgeHandle>& _hes) const;

    /// Get next halfedge within a halfface
    const HalfEdgeHandle next_halfedge_in_halfface(const HalfEdgeHandle& _heh, const HalfFaceHandle& _hfh) const;

    /// Get previous halfedge within a halfface
    const HalfEdgeHandle prev_halfedge_in_halfface(const HalfEdgeHandle& _heh, const HalfFaceHandle& _hfh) const;

    /// Get valence of vertex (number of incident edges)
    inline unsigned int valence(const VertexHandle& _vh) const {
        if(!v_bottom_up_) {
            std::cerr << "Could not get vertex valence: No bottom-up incidences for vertices available!" << std::endl;
            return 0u;
        }
        assert((unsigned int)_vh.idx() < outgoing_hes_per_vertex_.size());
        return outgoing_hes_per_vertex_[_vh.idx()].size();
    }

    /// Get valence of edge (number of incident faces)
    inline unsigned int valence(const EdgeHandle& _eh) const {
        if(!e_bottom_up_) {
            std::cerr << "Could not get edge valence: No bottom-up incidences for edges available!" << std::endl;
            return 0u;
        }
        assert((unsigned int)halfedge_handle(_eh, 0).idx() < incident_hfs_per_he_.size());
        return incident_hfs_per_he_[halfedge_handle(_eh, 0).idx()].size();
    }

    /// Get valence of face (number of incident edges)
    inline unsigned int valence(const FaceHandle& _fh) const {

        assert((unsigned int)_fh.idx() < faces_.size());
        return face(_fh).halfedges().size();
    }

    /// Get valence of cell (number of incident faces)
    inline unsigned int valence(const CellHandle& _ch) const {

        assert((unsigned int)_ch.idx() < cells_.size());
        return cell(_ch).halffaces().size();
    }

    //=====================================================================
    // Delete entities
    //=====================================================================

public:

    virtual VertexIter delete_vertex(const VertexHandle& _h);

    virtual EdgeIter delete_edge(const EdgeHandle& _h);

    virtual FaceIter delete_face(const FaceHandle& _h);

    virtual CellIter delete_cell(const CellHandle& _h);

	//added_by_Xiaoming Fu
	virtual VertexIter delete_vertex_swap(const VertexHandle& _h);
	virtual void delete_edge_swap(const EdgeHandle& _h);
	virtual void delete_face_swap(const FaceHandle& _h);
	virtual void delete_cell_swap(const CellHandle& _h);

protected:

    virtual void delete_multiple_vertices(const std::vector<bool>& _tag);

    virtual void delete_multiple_edges(const std::vector<bool>& _tag);

    virtual void delete_multiple_faces(const std::vector<bool>& _tag);

    virtual void delete_multiple_cells(const std::vector<bool>& _tag);

    class EdgeCorrector {
    public:
        EdgeCorrector(const std::vector<int>& _newIndices) :
            newIndices_(_newIndices) {}

        void operator()(Edge& _edge) {
            _edge.set_from_vertex(VertexHandle(newIndices_[_edge.from_vertex().idx()]));
            _edge.set_to_vertex(VertexHandle(newIndices_[_edge.to_vertex().idx()]));
        }
    private:
        const std::vector<int>& newIndices_;
    };

    class FaceCorrector {
    public:
        FaceCorrector(const std::vector<int>& _newIndices) :
            newIndices_(_newIndices) {}

        void operator()(Face& _face) {
            std::vector<HalfEdgeHandle> hes = _face.halfedges();
            for(std::vector<HalfEdgeHandle>::iterator he_it = hes.begin(),
                    he_end = hes.end(); he_it != he_end; ++he_it) {

                EdgeHandle eh = edge_handle(*he_it);
                unsigned char opp = (he_it->idx() - halfedge_handle(eh, 0).idx());
                *he_it = halfedge_handle(newIndices_[eh.idx()], opp);
            }
            _face.set_halfedges(hes);
        }
    private:
        const std::vector<int>& newIndices_;
    };

    class CellCorrector {
    public:
        CellCorrector(const std::vector<int>& _newIndices) :
            newIndices_(_newIndices) {}

        void operator()(Cell& _cell) {
            std::vector<HalfFaceHandle> hfs = _cell.halffaces();
            for(std::vector<HalfFaceHandle>::iterator hf_it = hfs.begin(),
                    hf_end = hfs.end(); hf_it != hf_end; ++hf_it) {

                FaceHandle fh = face_handle(*hf_it);
                unsigned char opp = (hf_it->idx() - halfface_handle(fh, 0).idx());
                *hf_it = halfface_handle(newIndices_[fh.idx()], opp);
            }
            _cell.set_halffaces(hfs);
        }
    private:
        const std::vector<int>& newIndices_;
    };

public:

    /** \brief Delete range of cells
     *
     * Deletes all cells in range [_first, _last].
     *
     * @param _first Iterator to first cell that is to be deleted
     * @param _last Iterator to last cell that is to be deleted
     * @return An iterator to the first cell after the deleted range
     */
    CellIter delete_cell_range(const CellIter& _first, const CellIter& _last);

public:

    /// Clear whole mesh
    virtual void clear(bool _clearProps = true) {

        edges_.clear();
        faces_.clear();
        cells_.clear();
        outgoing_hes_per_vertex_.clear();
        incident_hfs_per_he_.clear();
        incident_cell_per_hf_.clear();
        n_vertices_ = 0;

        if(_clearProps) {

            // Delete all property data
            clear_vertex_props();
            clear_edge_props();
            clear_halfedge_props();
            clear_face_props();
            clear_halfface_props();
            clear_cell_props();
            clear_mesh_props();

        } else {
            // Resize props
            resize_vprops(0u);
            resize_eprops(0u);
            resize_fprops(0u);
            resize_cprops(0u);
        }
    }

    //=====================================================================
    // Bottom-up Incidences
    //=====================================================================

public:

    void enable_bottom_up_incidences(bool _enable = true) {

        enable_vertex_bottom_up_incidences(_enable);
        enable_edge_bottom_up_incidences(_enable);
        enable_face_bottom_up_incidences(_enable);
    }

    void enable_vertex_bottom_up_incidences(bool _enable = true) {

        if(_enable && !v_bottom_up_) {
            // Vertex bottom-up incidences have to be
            // recomputed for the whole mesh
            compute_vertex_bottom_up_incidences();
        }

        if(!_enable) {
            outgoing_hes_per_vertex_.clear();
        }

        v_bottom_up_ = _enable;
    }

    void enable_edge_bottom_up_incidences(bool _enable = true) {

        if(_enable && !e_bottom_up_) {
            // Edge bottom-up incidences have to be
            // recomputed for the whole mesh
            compute_edge_bottom_up_incidences();

            if(f_bottom_up_) {
                std::for_each(edges_begin(), edges_end(),
                              fun::bind(&TopologyKernel::reorder_incident_halffaces, this, fun::placeholders::_1));
            }
        }

        if(!_enable) {
            incident_hfs_per_he_.clear();
        }

        e_bottom_up_ = _enable;
    }
	void write_incident_hfs_per_he()
	{
		FILE* f_hf_he = fopen("A:\\Code\\VolumeMeshProcessing\\Models\\Tetrahedrization\\sphere\\hehf.txt", "w");
		for(unsigned ii=0;ii<incident_hfs_per_he_.size();++ii)
		{
			for(unsigned jj=0;jj<incident_hfs_per_he_[ii].size(); ++jj)
			{
				fprintf(f_hf_he, "%d ", incident_hfs_per_he_[ii][jj].idx() );
			}
			fprintf(f_hf_he, "\n");
		}
		fclose(f_hf_he);
	}
    void enable_face_bottom_up_incidences(bool _enable = true) {

        bool updateOrder = false;
        if(_enable && !f_bottom_up_) {
            // Face bottom-up incidences have to be
            // recomputed for the whole mesh
            compute_face_bottom_up_incidences();

            updateOrder = true;
        }

        if(!_enable) {
            incident_cell_per_hf_.clear();
        }

        f_bottom_up_ = _enable;

        if(updateOrder) {
            if(e_bottom_up_) {
                std::for_each(edges_begin(), edges_end(),
                              fun::bind(&TopologyKernel::reorder_incident_halffaces, this, fun::placeholders::_1));
            }
        }
    }

    bool has_full_bottom_up_incidences() const {
        return (has_vertex_bottom_up_incidences() &&
                has_edge_bottom_up_incidences() &&
                has_face_bottom_up_incidences());
    }

    bool has_vertex_bottom_up_incidences() const { return v_bottom_up_; }

    bool has_edge_bottom_up_incidences() const { return e_bottom_up_; }

    bool has_face_bottom_up_incidences() const { return f_bottom_up_; }

	void set_add_cell_edge_flag();

private:

    void compute_vertex_bottom_up_incidences();

    void compute_edge_bottom_up_incidences();

    void compute_face_bottom_up_incidences();

    void reorder_incident_halffaces(const EdgeHandle& _eh);

    // Outgoing halfedges per vertex
    std::vector<std::vector<HalfEdgeHandle> > outgoing_hes_per_vertex_;

    // Incident halffaces per (directed) halfedge
    std::vector<std::vector<HalfFaceHandle> > incident_hfs_per_he_;

    // Incident cell (at most one) per halfface
    std::vector<CellHandle> incident_cell_per_hf_;

    bool v_bottom_up_;

    bool e_bottom_up_;

    bool f_bottom_up_;

	std::vector<int> add_cell_edge_flag;

    //=====================================================================
    // Connectivity
    //=====================================================================

public:

    /// Get halfface that is adjacent (w.r.t. a common halfedge) within the same cell
    HalfFaceHandle adjacent_halfface_in_cell(const HalfFaceHandle& _halfFaceHandle, const HalfEdgeHandle& _halfEdgeHandle) const;

    /// Get cell that is incident to the given halfface
    CellHandle incident_cell(const HalfFaceHandle& _halfFaceHandle) const;

    bool is_boundary(const HalfFaceHandle& _halfFaceHandle) const {
        return _halfFaceHandle.idx() >= 0 && (unsigned int)_halfFaceHandle.idx() < incident_cell_per_hf_.size() &&
                incident_cell_per_hf_[_halfFaceHandle.idx()] == InvalidCellHandle;
    }

    bool is_boundary(const FaceHandle& _faceHandle) const {
        return  is_boundary(halfface_handle(_faceHandle, 0)) ||
                is_boundary(halfface_handle(_faceHandle, 1));
    }

    bool is_boundary(const EdgeHandle& _edgeHandle) const {
        if(!e_bottom_up_) {
            std::cerr << "Error: Function is_boundary() needs bottom-up incidences for edges!" << std::endl;
            return false;
        }
        for(HalfEdgeHalfFaceIter hehf_it = hehf_iter(halfedge_handle(_edgeHandle, 0));
                hehf_it.valid(); ++hehf_it) {
            if(is_boundary(face_handle(*hehf_it))) {
                return true;
            }
        }
        return false;
    }

    bool is_boundary(const HalfEdgeHandle& _halfedgeHandle) const {
        if(!e_bottom_up_) {
            std::cerr << "Error: Function is_boundary() needs bottom-up incidences for edges!" << std::endl;
            return false;
        }
        for(HalfEdgeHalfFaceIter hehf_it = hehf_iter(_halfedgeHandle);
                hehf_it.valid(); ++hehf_it) {
            if(is_boundary(face_handle(*hehf_it))) {
                return true;
            }
        }
        return false;
    }

    bool is_boundary(const VertexHandle& _vertexHandle) const {
        if(!v_bottom_up_) {
            std::cerr << "Error: Function is_boundary() needs bottom-up incidences for vertices!" << std::endl;
            return false;
        }
        for(VertexOHalfEdgeIter voh_it = voh_iter(_vertexHandle); voh_it.valid(); ++voh_it) {
            if(is_boundary(*voh_it)) return true;
        }
        return false;
    }

    unsigned int n_vertices_in_cell(const CellHandle& _ch) const {

        std::set<VertexHandle> vertices;
        std::vector<HalfFaceHandle> hfs = cell(_ch).halffaces();
        for(std::vector<HalfFaceHandle>::const_iterator hf_it = hfs.begin();
                hf_it != hfs.end(); ++hf_it) {
            std::vector<HalfEdgeHandle> hes = halfface(*hf_it).halfedges();
            for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin();
                he_it != hes.end(); ++he_it) {
                vertices.insert(halfedge(*he_it).to_vertex());
            }
        }
        return vertices.size();
    }

    //=========================================================================

    /*
     * Non-virtual functions
     */

    const Edge opposite_halfedge(const Edge& _edge) const {
        return Edge(_edge.to_vertex(), _edge.from_vertex());
    }

    const Face opposite_halfface(const Face& _face) const 
	{
		std::vector<HalfEdgeHandle> temp_opp_halfedges;
		Face return_face( temp_opp_halfedges );
		std::vector<HalfEdgeHandle>& opp_halfedges = return_face.get_halfedges();
		unsigned n = _face.halfedges().size(); opp_halfedges.resize(n);
		for(unsigned i=0;i<n;++i)
		{
			opp_halfedges[n - 1 - i] = opposite_halfedge_handle(_face.halfedges()[i]);
		}
		return return_face;
        /*for(std::vector<HalfEdgeHandle>::const_iterator it = _face.halfedges().begin(); it
                != _face.halfedges().end(); ++it) 
		{
            opp_halfedges.insert(opp_halfedges.begin(), opposite_halfedge_handle(*it));
        }*/
       //return Face(opp_halfedges);
    }

    /*
     * Static functions
     */

    /// Conversion function
    static inline HalfEdgeHandle halfedge_handle(const EdgeHandle& _h, const unsigned char _subIdx) {
        // Is handle in range?
        if(_h.idx() < 0 || _subIdx > 1) return InvalidHalfEdgeHandle;
        return HalfEdgeHandle((2 * _h.idx()) + (_subIdx ? 1 : 0));
    }

    /// Conversion function
    static inline HalfFaceHandle halfface_handle(const FaceHandle& _h, const unsigned char _subIdx) {
        // Is handle in range?
        if(_h.idx() < 0 || _subIdx > 1) return InvalidHalfFaceHandle;
        return HalfFaceHandle((2 * _h.idx()) + (_subIdx ? 1 : 0));
    }

    /// Handle conversion
    static inline EdgeHandle edge_handle(const HalfEdgeHandle& _h) {
        // Is handle in range?
        if(_h.idx() < 0) return InvalidEdgeHandle;
        return EdgeHandle((int)(_h.idx() / 2));
    }

    static inline FaceHandle face_handle(const HalfFaceHandle& _h) {
        // Is handle in range?
        if(_h.idx() < 0) return InvalidFaceHandle;
        return FaceHandle((int)(_h.idx() / 2));
    }

    static inline HalfEdgeHandle opposite_halfedge_handle(const HalfEdgeHandle& _h) {
        // Is handle in range?
        if(_h.idx() < 0) return InvalidHalfEdgeHandle;

        // Is handle even?
        if(_h.idx() % 2 == 0) {
            return HalfEdgeHandle(_h.idx() + 1);
        }
        return HalfEdgeHandle(_h.idx() - 1);
    }

    static inline HalfFaceHandle opposite_halfface_handle(const HalfFaceHandle& _h) {
        // Is handle in range?
        if(_h.idx() < 0) return InvalidHalfFaceHandle;

        // Is handle even?
        if(_h.idx() % 2 == 0) {
            return HalfFaceHandle(_h.idx() + 1);
        }
        return HalfFaceHandle(_h.idx() - 1);
    }

protected:

    // List of edges
    std::vector<Edge> edges_;

    // List of faces
    std::vector<Face> faces_;

    // List of cells
    std::vector<Cell> cells_;
};

}

#endif /* TOPOLOGYKERNEL_HH_ */
