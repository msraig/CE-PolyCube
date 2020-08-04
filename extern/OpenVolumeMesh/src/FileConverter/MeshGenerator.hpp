/*
 * MeshGenerator.hh
 *
 *  Created on: Mar 15, 2012
 *      Author: kremer
 */

#ifndef MESHGENERATOR_HH_
#define MESHGENERATOR_HH_

#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include <boost/shared_ptr.hpp>
#include <boost/progress.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>

class MeshGenerator {
private:

    typedef OpenVolumeMesh::VertexHandle    VertexHandle;
    typedef OpenVolumeMesh::EdgeHandle      EdgeHandle;
    typedef OpenVolumeMesh::HalfEdgeHandle  HalfEdgeHandle;
    typedef OpenVolumeMesh::FaceHandle      FaceHandle;
    typedef OpenVolumeMesh::HalfFaceHandle  HalfFaceHandle;
    typedef OpenVolumeMesh::CellHandle      CellHandle;

    typedef boost::tuple<VertexHandle, VertexHandle, VertexHandle> FaceTuple;

public:

    typedef OpenVolumeMesh::GeometricPolyhedralMeshV3d PolyhedralMesh;

    typedef OpenVolumeMesh::Geometry::Vec3d Vec3d;

    MeshGenerator(PolyhedralMesh& _mesh) : v_component_(0), mesh_(_mesh), progress_() {}
    MeshGenerator(const MeshGenerator& _cpy) :
        v_component_(_cpy.v_component_),
        vertex_(0.0, 0.0, 0.0),
        c_vertices_(),
        faceMap_(),
        mesh_(_cpy.mesh_),
        progress_() {}

    void add_vertex_component(double _comp) {

        if(v_component_ > 2) {
            std::cerr << "Vertices of dimension higher than three not supported!" << std::endl;
            return;
        }
        vertex_[v_component_] = _comp;
        ++v_component_;
        if(v_component_ == 3) {
            add_vertex();
        }
    }

    void add_vertex() {

        OpenVolumeMesh::VertexHandle vh = mesh_.add_vertex(vertex_);
        //std::cerr << "Added vertex " << mesh_.vertex(vh) << std::endl;
        v_component_ = 0;
    }

    void add_cell_vertex(unsigned int _idx) {

        assert(_idx > 0);

        c_vertices_.push_back(OpenVolumeMesh::VertexHandle((int)_idx - 1));
        if(c_vertices_.size() == 4) {

            add_tetrahedral_cell();
//            std::cerr << "Adding cell (" << c_vertices_[0] << ", " << c_vertices_[1] <<
//                    ", " << c_vertices_[2] << ", " << c_vertices_[3] << ")" << std::endl;
            c_vertices_.clear();
        }
    }

    void set_num_cells(unsigned int _n) {

        if(progress_.get() == NULL) {
            progress_.reset(new boost::progress_display(_n));
        }
    }

    void add_tetrahedral_cell() {

        if(c_vertices_.size() != 4) {
            std::cerr << "The specified cell is not incident to four vertices!" << std::endl;
            return;
        }

        // Get cell's mid-point
        Vec3d midP(0.0, 0.0, 0.0);
        double valence = 0.0;
        for(std::vector<OpenVolumeMesh::VertexHandle>::const_iterator it = c_vertices_.begin();
                it != c_vertices_.end(); ++it) {
            midP += mesh_.vertex(*it);
            valence += 1.0;
        }
        midP /= valence;

        // Sort vertex vector
        std::sort(c_vertices_.begin(), c_vertices_.end());

        std::vector<FaceTuple> tuples;

        // Create face tuple for all vertex combinations
        tuples.push_back(FaceTuple(c_vertices_[0], c_vertices_[1], c_vertices_[2]));
        tuples.push_back(FaceTuple(c_vertices_[1], c_vertices_[2], c_vertices_[3]));
        tuples.push_back(FaceTuple(c_vertices_[0], c_vertices_[2], c_vertices_[3]));
        tuples.push_back(FaceTuple(c_vertices_[0], c_vertices_[1], c_vertices_[3]));

        // Collect cell's half-faces in here
        std::vector<HalfFaceHandle> cell_halffaces;

        for(std::vector<FaceTuple>::const_iterator it = tuples.begin();
                it != tuples.end(); ++it) {

            // Check if face exists for current tuple
            FaceMap::iterator f = faceMap_.find(*it);
            if(f == faceMap_.end()) {
                // Face does not exist, create it

                // Find right orientation, s.t. normal
                // points inside the cell

                Vec3d e1 = mesh_.vertex(it->get<1>()) - mesh_.vertex(it->get<0>());
                Vec3d e2 = mesh_.vertex(it->get<2>()) - mesh_.vertex(it->get<1>());

                // Get face normal (cross product)
                Vec3d n = (e1 % e2).normalize();

                std::vector<VertexHandle> v_vec;
                v_vec.push_back(it->get<0>());
                v_vec.push_back(it->get<1>());
                v_vec.push_back(it->get<2>());
                FaceHandle fh = mesh_.add_face(v_vec);

                // Add face to face map
                faceMap_[*it] = fh;

                // Check whether normal points inside cell
                if(((midP - mesh_.vertex(it->get<0>())) | n) > 0.0) {

                    // Normal points inside cell, just add half-face 0
                    // Add corresponding half-face to cell definition
                    cell_halffaces.push_back(mesh_.halfface_handle(fh, 0));

                } else {

                    // Normal points outside cell, just add half-face 1
                    // Add corresponding half-face to cell definition
                    cell_halffaces.push_back(mesh_.halfface_handle(fh, 1));
                }

            } else {

                // Face exists, find right orientation
                FaceHandle fh = f->second;

                std::vector<HalfEdgeHandle> hes = mesh_.face(fh).halfedges();

                assert(hes.size() == 3);

                Vec3d e1 = mesh_.vertex(mesh_.halfedge(hes[0]).to_vertex()) -
                        mesh_.vertex(mesh_.halfedge(hes[0]).from_vertex());
                Vec3d e2 = mesh_.vertex(mesh_.halfedge(hes[1]).to_vertex()) -
                        mesh_.vertex(mesh_.halfedge(hes[1]).from_vertex());

                Vec3d n = (e1 % e2).normalize();

                if(((midP - mesh_.vertex(mesh_.halfedge(hes[0]).from_vertex())) | n) > 0.0) {
                    // Normal points inside cell
                    cell_halffaces.push_back(mesh_.halfface_handle(fh, 0));
                } else {
                    // Normal points outisde cell
                    cell_halffaces.push_back(mesh_.halfface_handle(fh, 1));
                }
            }
        }

        // Check whether cell definition contains four half-faces
        assert(cell_halffaces.size() == 4);

        // Finally, add cell
#ifndef NDEBUG
        mesh_.add_cell(cell_halffaces, true);
#else
        mesh_.add_cell(cell_halffaces, false);
#endif

        // Increase progress counter
        if((progress_.get() != NULL) && (progress_->expected_count() != 0))
            ++(*progress_);
    }

private:

    typedef std::map<FaceTuple, OpenVolumeMesh::FaceHandle> FaceMap;

    unsigned int v_component_;
    OpenVolumeMesh::Geometry::Vec3d vertex_;

    std::vector<VertexHandle> c_vertices_;

    FaceMap faceMap_;

    PolyhedralMesh& mesh_;

    boost::shared_ptr<boost::progress_display> progress_;
};

#endif /* MESHGENERATOR_HH_ */
