#pragma warning( disable :4018 4267 4244 4838)
#include "sus_lib.h"
#include <array>
#include <limits>
#include "Polycube_Deformation.h"
#include "omp.h"
#include "Sparse_Matrix.h"
#include "Sparse_Solver.h"
#include "ConvexQuadOptimization.h"
#include "cholmod.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "mosek.h"
#include <Queue>
#include <set>
#include "fmath.hpp"
#include "SmallMat.h"
#include "deformation_LBFGS.h"
#include "Helper.h"
#include "adjust_orientation_LBFGS.h"
#include "permutohedral.h"
#include "Polycube_Boundary_Map_IVF.h"
#define FEATURE_COEFF 100.0
using ig::CVec;
using OpenVolumeMesh::VertexHandle;
using OpenVolumeMesh::EdgeHandle;
using Eigen::Vector3d;
using Eigen::Matrix3d;
polycube_deformation_interface::polycube_deformation_interface()
{
    prepare_ok = false;
    sigma_s = 1.0; sigma_r = 0.3; RGNF_iter_count = 1;
    fixed_chart.clear();
    fixed_label.clear();
    feature_edge_change_flag = false;
    deformation_prepare_OK = false;
    polyline_ns_flag = false;
}
polycube_deformation_interface::~polycube_deformation_interface()
{
    ;
}
void polycube_deformation_interface::prepare_for_deformation(VolumeMesh* mesh_)
{
    std::vector<OpenVolumeMesh::VertexHandle> hfv_vec(3);
    std::vector<int> hfv_vec_id(3);
    std::vector<OpenVolumeMesh::HalfFaceHandle> cell_hfh_vec(4);
    unsigned nv = mesh_->n_vertices(); unsigned nc = mesh_->n_cells();
    vertex_cell.clear(); cell_vertex.clear();
    vertex_cell.resize(nv); cell_vertex.resize(nc);
    cell_vertex_right_order.clear();
    cell_vertex_right_order.resize(nc);
    cell_vertex_vertex.clear(); cell_S.resize(nc);
    cell_vertex_vertex.resize(nc);
    vertex_cell_vertex.clear(); vertex_cell_vertex.resize(nv);
    vcv_S.clear(); vcv_S.resize(nv);
    cell_volume.clear(); cell_volume.resize(nc, 0.0);
    Eigen::Matrix3d VS, IS; std::vector<double> S_(9);
    for (OpenVolumeMesh::CellIter c_it = mesh_->cells_begin(); c_it != mesh_->cells_end(); ++c_it)
    {
        double cv_count = 0.0; int c_id = c_it->idx();
        const std::vector<OpenVolumeMesh::VertexHandle>& vertices_ = mesh_->cell(*c_it).vertices();
        cell_vertex[c_id].clear();
        for (unsigned i = 0; i < vertices_.size(); ++i)
        {
            int v_id = vertices_[i];
            cell_vertex[c_id].push_back(v_id);
        }
        OpenVolumeMesh::Geometry::Vec3d cp = mesh_->vertex(vertices_[0]);
        OpenVolumeMesh::Geometry::Vec3d cr = mesh_->vertex(vertices_[1]) - cp;
        OpenVolumeMesh::Geometry::Vec3d cs = mesh_->vertex(vertices_[2]) - cp;
        OpenVolumeMesh::Geometry::Vec3d ct = mesh_->vertex(vertices_[3]) - cp;
        VS(0, 0) = cr[0]; VS(1, 0) = cr[1]; VS(2, 0) = cr[2];
        VS(0, 1) = cs[0]; VS(1, 1) = cs[1]; VS(2, 1) = cs[2];
        VS(0, 2) = ct[0]; VS(1, 2) = ct[1]; VS(2, 2) = ct[2];
        IS = VS.inverse();
        cell_S[c_id] = IS;
        if (VS.determinant() >= 0)
        {
            cell_vertex_right_order[c_id] = cell_vertex[c_id];
        }
        else
        {
            cell_vertex_right_order[c_id] = cell_vertex[c_id];
            cell_vertex_right_order[c_id][0] = cell_vertex[c_id][1];
            cell_vertex_right_order[c_id][1] = cell_vertex[c_id][0];
        }
        mesh_->get_halffaces_from_cell(*c_it, cell_hfh_vec);
        for (unsigned i = 0; i < vertices_.size(); ++i)
        {
            int v_id = vertices_[i];
            if (v_id < 0 || v_id > nv)
            {
                printf("%d\n", v_id);
            }
            if (vertex_cell[v_id].size() == 0) vertex_cell[v_id].reserve(15);
            vertex_cell[v_id].push_back(c_it->idx());
            for (unsigned j = 0; j < cell_hfh_vec.size(); ++j)
            {
                bool find_hf = true;
                mesh_->get_vertices_from_halfface(cell_hfh_vec[j], hfv_vec);
                for (unsigned jj = 0; jj < hfv_vec.size(); ++jj)
                {
                    if (hfv_vec[jj] == v_id)
                    {
                        find_hf = false;
                        break;
                    }
                }
                if (find_hf)
                {
                    hfv_vec_id[0] = hfv_vec[0].idx(); hfv_vec_id[1] = hfv_vec[1].idx(); hfv_vec_id[2] = hfv_vec[2].idx();
                    cell_vertex_vertex[c_id].push_back(hfv_vec_id);
                    vertex_cell_vertex[v_id].push_back(hfv_vec_id);
                    OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(vertices_[i]);
                    OpenVolumeMesh::Geometry::Vec3d r = mesh_->vertex(hfv_vec[0]) - p;
                    OpenVolumeMesh::Geometry::Vec3d s = mesh_->vertex(hfv_vec[1]) - p;
                    OpenVolumeMesh::Geometry::Vec3d t = mesh_->vertex(hfv_vec[2]) - p;
                    VS(0, 0) = r[0]; VS(1, 0) = r[1]; VS(2, 0) = r[2];
                    VS(0, 1) = s[0]; VS(1, 1) = s[1]; VS(2, 1) = s[2];
                    VS(0, 2) = t[0]; VS(1, 2) = t[1]; VS(2, 2) = t[2];
                    cell_volume[c_id] = -VS.determinant();
                    IS = VS.inverse();
                    S_[0] = IS(0, 0); S_[1] = IS(0, 1); S_[2] = IS(0, 2);
                    S_[3] = IS(1, 0); S_[4] = IS(1, 1); S_[5] = IS(1, 2);
                    S_[6] = IS(2, 0); S_[7] = IS(2, 1); S_[8] = IS(2, 2);
                    vcv_S[v_id].push_back(S_);
                    break;
                }
            }
        }
    }
    dpx.resize(nv); dpy.resize(nv); dpz.resize(nv); max_vc_size = 0; src_pos.resize(nv);
    is_bv.clear(); is_bv.resize(nv, -1);
    for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
    {
        int v_id = v_it->idx();
        OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(*v_it);
        dpx[v_id] = p[0]; dpy[v_id] = p[1]; dpz[v_id] = p[2];
        src_pos[v_id] = p;
        int vc_size = vertex_cell[v_id].size();
        if (vc_size > max_vc_size)
        {
            max_vc_size = vc_size;
        }
        if (mesh_->is_boundary(*v_it))
        {
            is_bv[v_id] = 1;
        }
    }
    int max_num_t = omp_get_max_threads();
    vc_pos_x_omp.resize(max_num_t); vc_pos_y_omp.resize(max_num_t); vc_pos_z_omp.resize(max_num_t);
    vc_pos_x2_omp.resize(max_num_t); vc_pos_y2_omp.resize(max_num_t); vc_pos_z2_omp.resize(max_num_t);
    vc_pos_x3_omp.resize(max_num_t); vc_pos_y3_omp.resize(max_num_t); vc_pos_z3_omp.resize(max_num_t);
    vc_n_cross_x_omp.resize(max_num_t); vc_n_cross_y_omp.resize(max_num_t); vc_n_cross_z_omp.resize(max_num_t);
    vc_S_omp.resize(max_num_t); exp_vec_omp.resize(max_num_t); large_dis_flag_omp.resize(max_num_t);
    gx_vec_omp.resize(max_num_t); gy_vec_omp.resize(max_num_t); gz_vec_omp.resize(max_num_t);
    mu_vec_omp.resize(max_num_t);
    for (int i = 0; i < max_num_t; ++i)
    {
        vc_pos_x_omp[i].resize(max_vc_size); vc_pos_y_omp[i].resize(max_vc_size); vc_pos_z_omp[i].resize(max_vc_size);
        vc_pos_x2_omp[i].resize(4 * max_vc_size); vc_pos_y2_omp[i].resize(4 * max_vc_size); vc_pos_z2_omp[i].resize(4 * max_vc_size);
        vc_pos_x3_omp[i].resize(4 * max_vc_size); vc_pos_y3_omp[i].resize(4 * max_vc_size); vc_pos_z3_omp[i].resize(4 * max_vc_size);
        vc_n_cross_x_omp[i].resize(max_vc_size); vc_n_cross_y_omp[i].resize(max_vc_size); vc_n_cross_z_omp[i].resize(max_vc_size);
        vc_S_omp[i].resize(max_vc_size); exp_vec_omp[i].resize(max_vc_size); large_dis_flag_omp[i].resize(max_vc_size);
        gx_vec_omp[i].resize(max_vc_size); gy_vec_omp[i].resize(max_vc_size); gz_vec_omp[i].resize(max_vc_size);
        mu_vec_omp[i].resize(max_vc_size);
        for (int j = 0; j < max_vc_size; ++j)
        {
            vc_S_omp[i][j].resize(9);
        }
    }
    bfv_id.clear(); bfv_id.resize(nc, OpenVolumeMesh::Geometry::Vec3i(-1, -1, -1));
    bvf_id.clear(); bvf_id.resize(nv); boundary_face_number = 0; avg_boundary_edge_length = 0.0;
    bef_id.clear(); bef_id.resize(mesh_->n_edges());
    for (OpenVolumeMesh::FaceIter f_it = mesh_->faces_begin(); f_it != mesh_->faces_end(); ++f_it)
    {
        OpenVolumeMesh::HalfFaceHandle hfh0 = mesh_->halfface_handle(*f_it, 0);
        OpenVolumeMesh::HalfFaceHandle hfh1 = mesh_->halfface_handle(*f_it, 1);
        OpenVolumeMesh::CellHandle ch0 = mesh_->incident_cell(hfh0);
        OpenVolumeMesh::CellHandle ch1 = mesh_->incident_cell(hfh1);
        if (ch0 == VolumeMesh::InvalidCellHandle && ch1 != VolumeMesh::InvalidCellHandle)
        {
            int hfv_count = 0; int c_id = ch1.idx();
            for (OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh0); hfv_it; ++hfv_it)
            {
                int hfv_id = hfv_it->idx();
                bfv_id[c_id][hfv_count] = hfv_id;
                bvf_id[hfv_id].push_back(c_id);
                ++hfv_count;
            }
            OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face(*f_it);
            std::vector<OpenVolumeMesh::HalfEdgeHandle>& heh_vec = face.get_halfedges();
            for (int i = 0; i < heh_vec.size(); ++i)
            {
                int e_id = mesh_->edge_handle(heh_vec[i]);
                bef_id[e_id].push_back(c_id);
            }
            ++boundary_face_number;
        }
        else if (ch1 == VolumeMesh::InvalidCellHandle && ch0 != VolumeMesh::InvalidCellHandle)
        {
            int hfv_count = 0; int c_id = ch0.idx();
            for (OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh1); hfv_it; ++hfv_it)
            {
                int hfv_id = hfv_it->idx();
                bfv_id[c_id][hfv_count] = hfv_id;
                bvf_id[hfv_id].push_back(c_id);
                ++hfv_count;
            }
            OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face(*f_it);
            std::vector<OpenVolumeMesh::HalfEdgeHandle>& heh_vec = face.get_halfedges();
            for (int i = 0; i < heh_vec.size(); ++i)
            {
                int e_id = mesh_->edge_handle(heh_vec[i]);
                bef_id[e_id].push_back(c_id);
            }
            ++boundary_face_number;
        }
        else if (ch1 == VolumeMesh::InvalidCellHandle && ch0 == VolumeMesh::InvalidCellHandle)
        {
            printf("Error : Both two halffaces have no cells!\n");
        }
    }
    bvv_id.clear(); bvv_id.resize(nv);
    for (int i = 0; i < nv; ++i)
    {
        if (bvf_id[i].size() > 0)
        {
            std::vector<int>& one_ring_f_id = bvf_id[i]; int one_ring_size = one_ring_f_id.size();
            std::vector<int> one_ring_v_id(one_ring_size); std::vector<int> one_ring_v_f_id(one_ring_size);
            std::vector<int> one_ring_order(one_ring_size);
            for (int j = 0; j < one_ring_size; ++j)
            {
                OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[one_ring_f_id[j]];
                for (int k = 0; k < 3; ++k)
                {
                    if (i == fv_id[k])
                    {
                        one_ring_order[j] = k;
                        break;
                    }
                }
            }
            OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[one_ring_f_id[0]];
            int start_v = fv_id[(one_ring_order[0] + 1) % 3];
            int v2 = start_v; one_ring_v_id[0] = v2; one_ring_v_f_id[0] = one_ring_f_id[0];
            int v3 = fv_id[(one_ring_order[0] + 2) % 3];
            for (int j = 1; j < one_ring_size; ++j)
            {
                v2 = v3;
                if (v2 == start_v) break;
                one_ring_v_id[j] = v2;
                for (int k = 1; k < one_ring_size; ++k)
                {
                    OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[one_ring_f_id[k]];
                    if (fv_id[(one_ring_order[k] + 1) % 3] == v2)
                    {
                        v3 = fv_id[(one_ring_order[k] + 2) % 3];
                        one_ring_v_f_id[j] = one_ring_f_id[k];
                        break;
                    }
                }
            }
            bvv_id[i] = one_ring_v_id;
            bvf_id[i] = one_ring_v_f_id;
        }
    }
    a2b.clear();
    a2b.resize(nv, -1);
    boundary_verts.clear();
    for (size_t i = 0; i < nv; i++)
    {
        if (bvf_id[i].size() > 0)
        {
            a2b[i] = boundary_verts.size();
            boundary_verts.push_back(i);
        }
    }
    avg_boundary_edge_length = 0.0;
    min_boundary_edge_length = -1;
    OpenVolumeMesh::Geometry::Vec3d p0, p1, p2; double count = 0.0;
    for (int i = 0; i < bfv_id.size(); ++i)
    {
        OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
        if (one_bfv_id[0] < 0) continue;
        p0[0] = dpx[one_bfv_id[0]]; p0[1] = dpy[one_bfv_id[0]]; p0[2] = dpz[one_bfv_id[0]];
        p1[0] = dpx[one_bfv_id[1]]; p1[1] = dpy[one_bfv_id[1]]; p1[2] = dpz[one_bfv_id[1]];
        p2[0] = dpx[one_bfv_id[2]]; p2[1] = dpy[one_bfv_id[2]]; p2[2] = dpz[one_bfv_id[2]];
        avg_boundary_edge_length += (p0 - p1).norm();
        avg_boundary_edge_length += (p1 - p2).norm();
        avg_boundary_edge_length += (p2 - p0).norm();
        if (i == 0)
        {
            min_boundary_edge_length = (p0 - p1).norm();
        }
        if (min_boundary_edge_length > (p0 - p1).norm())
            min_boundary_edge_length = (p0 - p1).norm();
        if (min_boundary_edge_length > (p2 - p1).norm())
            min_boundary_edge_length = (p2 - p1).norm();
        if (min_boundary_edge_length > (p2 - p0).norm())
            min_boundary_edge_length = (p2 - p0).norm();
        count += 3.0;
    }
    avg_boundary_edge_length /= count;
    boundary_vertface_number = 0;
    std::vector<int> cell_flag(bfv_id.size(), -1);
    for (size_t i = 0; i < nv; i++)
    {
        for (size_t j = 0; j < vertex_cell[i].size(); j++)
        {
            cell_flag[vertex_cell[i][j]] = 1;
        }
    }
    for (size_t i = 0; i < nc; i++)
    {
        if (cell_flag[i] == 1)
            boundary_vertface_number++;
    }
    printf("Average Boundary Edge length : %e\n", avg_boundary_edge_length);
    printf("Min boundary edge length: %e\n", min_boundary_edge_length);
    change_big_flag.clear(); last_bfn.clear();
    change_big_flag.resize(nc, -1); last_bfn.resize(nc);
    bound_K = 1.0;
    prepare_ok = true;
}
void polycube_deformation_interface::prepare_for_deformation(TetStructure<double>* tet_mesh_)
{
    if (deformation_prepare_OK) return;
    deformation_prepare_OK = true;
    std::vector<OpenVolumeMesh::VertexHandle> hfv_vec(3);
    std::vector<int> hfv_vec_id(3);
    std::vector<OpenVolumeMesh::HalfFaceHandle> cell_hfh_vec(4);
    unsigned nv = tet_mesh_->tetra_vertices.size(); unsigned nc = tet_mesh_->tetras.size();
    vert_update_type.resize(nv, 0);
    vertex_cell.clear(); cell_vertex.clear();
    vertex_cell.resize(nv); cell_vertex.resize(nc);
    cell_vertex_right_order.clear();
    cell_vertex_right_order.resize(nc);
    cell_vertex_vertex.clear(); cell_S.resize(nc);
    cell_vertex_vertex.resize(nc);
    vertex_cell_vertex.clear(); vertex_cell_vertex.resize(nv);
    vcv_S.clear(); vcv_S.resize(nv);
    cell_volume.clear(); cell_volume.resize(nc, 0.0);
    Eigen::Matrix3d VS, IS; std::vector<double> S_(9);
    static int s_tet_id[4][3] = { { 2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2} };
    const std::vector<Tetrahedron<double>*> &tetras = tet_mesh_->tetras;
    const std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
    for (size_t c_id = 0; c_id < nc; ++c_id)
    {
        double cv_count = 0.0;
        TetVertex<double> **tet_vert_ = tet_mesh_->tetras[c_id]->vertex;
        cell_vertex[c_id].clear();
        for (unsigned i = 0; i < 4; ++i)
        {
            int v_id = tet_vert_[i]->id;
            cell_vertex[c_id].push_back(v_id);
        }
        CVec<double, 3> cp = tet_vert_[0]->pos;
        CVec<double, 3> cr = tet_vert_[1]->pos - cp;
        CVec<double, 3> cs = tet_vert_[2]->pos - cp;
        CVec<double, 3> ct = tet_vert_[3]->pos - cp;
        VS(0, 0) = cr[0]; VS(1, 0) = cr[1]; VS(2, 0) = cr[2];
        VS(0, 1) = cs[0]; VS(1, 1) = cs[1]; VS(2, 1) = cs[2];
        VS(0, 2) = ct[0]; VS(1, 2) = ct[1]; VS(2, 2) = ct[2];
        IS = VS.inverse();
        cell_S[c_id] = IS;
        if (VS.determinant() >= 0)
        {
            cell_vertex_right_order[c_id] = cell_vertex[c_id];
        }
        else
        {
            cell_vertex_right_order[c_id] = cell_vertex[c_id];
            cell_vertex_right_order[c_id][0] = cell_vertex[c_id][1];
            cell_vertex_right_order[c_id][1] = cell_vertex[c_id][0];
        }
        for (unsigned i = 0; i < 4; ++i)
        {
            int v_id = tet_vert_[i]->id;
            if (v_id < 0 || v_id > nv)
            {
                printf("%d\n", v_id);
            }
            if (vertex_cell[v_id].size() == 0) vertex_cell[v_id].reserve(15);
            vertex_cell[v_id].push_back(c_id);
            hfv_vec_id[0] = tetras[c_id]->vertex[s_tet_id[i][0]]->id;
            hfv_vec_id[1] = tetras[c_id]->vertex[s_tet_id[i][1]]->id;
            hfv_vec_id[2] = tetras[c_id]->vertex[s_tet_id[i][2]]->id;
            cell_vertex_vertex[c_id].push_back(hfv_vec_id);
            vertex_cell_vertex[v_id].push_back(hfv_vec_id);
            CVec<double, 3> p = tetras[c_id]->vertex[i]->pos;
            CVec<double, 3> r = tetras[c_id]->vertex[s_tet_id[i][0]]->pos - p;
            CVec<double, 3> s = tetras[c_id]->vertex[s_tet_id[i][1]]->pos - p;
            CVec<double, 3> t = tetras[c_id]->vertex[s_tet_id[i][2]]->pos - p;
            VS(0, 0) = r[0]; VS(1, 0) = r[1]; VS(2, 0) = r[2];
            VS(0, 1) = s[0]; VS(1, 1) = s[1]; VS(2, 1) = s[2];
            VS(0, 2) = t[0]; VS(1, 2) = t[1]; VS(2, 2) = t[2];
            cell_volume[c_id] = -VS.determinant();
            IS = VS.inverse();
            S_[0] = IS(0, 0); S_[1] = IS(0, 1); S_[2] = IS(0, 2);
            S_[3] = IS(1, 0); S_[4] = IS(1, 1); S_[5] = IS(1, 2);
            S_[6] = IS(2, 0); S_[7] = IS(2, 1); S_[8] = IS(2, 2);
            vcv_S[v_id].push_back(S_);
        }
    }
    dpx.clear(); dpy.clear(); dpz.clear();
    dpx.resize(nv); dpy.resize(nv); dpz.resize(nv); max_vc_size = 0; src_pos.resize(nv);
    is_bv.clear(); is_bv.resize(nv, -1);
    for (int i = 0; i < nv; i++)
    {
        int v_id = i;
        CVec<double, 3> p = tetra_vertices[i]->pos;
        dpx[v_id] = p[0]; dpy[v_id] = p[1]; dpz[v_id] = p[2];
        OpenVolumeMesh::Geometry::Vec3d p_temp;
        p_temp[0] = p[0]; p_temp[1] = p[1]; p_temp[2] = p[2];
        src_pos[v_id] = p_temp;
        int vc_size = vertex_cell[v_id].size();
        if (vc_size > max_vc_size)
        {
            max_vc_size = vc_size;
        }
        if (tetra_vertices[i]->boundary)
        {
            is_bv[v_id] = 1;
        }
    }
    int max_num_t = omp_get_max_threads();
    vc_pos_x_omp.resize(max_num_t); vc_pos_y_omp.resize(max_num_t); vc_pos_z_omp.resize(max_num_t);
    vc_pos_x2_omp.resize(max_num_t); vc_pos_y2_omp.resize(max_num_t); vc_pos_z2_omp.resize(max_num_t);
    vc_pos_x3_omp.resize(max_num_t); vc_pos_y3_omp.resize(max_num_t); vc_pos_z3_omp.resize(max_num_t);
    vc_n_cross_x_omp.resize(max_num_t); vc_n_cross_y_omp.resize(max_num_t); vc_n_cross_z_omp.resize(max_num_t);
    vc_S_omp.resize(max_num_t); exp_vec_omp.resize(max_num_t); large_dis_flag_omp.resize(max_num_t);
    gx_vec_omp.resize(max_num_t); gy_vec_omp.resize(max_num_t); gz_vec_omp.resize(max_num_t);
    mu_vec_omp.resize(max_num_t);
    for (int i = 0; i < max_num_t; ++i)
    {
        vc_pos_x_omp[i].resize(max_vc_size); vc_pos_y_omp[i].resize(max_vc_size); vc_pos_z_omp[i].resize(max_vc_size);
        vc_pos_x2_omp[i].resize(4 * max_vc_size); vc_pos_y2_omp[i].resize(4 * max_vc_size); vc_pos_z2_omp[i].resize(4 * max_vc_size);
        vc_pos_x3_omp[i].resize(4 * max_vc_size); vc_pos_y3_omp[i].resize(4 * max_vc_size); vc_pos_z3_omp[i].resize(4 * max_vc_size);
        vc_n_cross_x_omp[i].resize(max_vc_size); vc_n_cross_y_omp[i].resize(max_vc_size); vc_n_cross_z_omp[i].resize(max_vc_size);
        vc_S_omp[i].resize(max_vc_size); exp_vec_omp[i].resize(max_vc_size); large_dis_flag_omp[i].resize(max_vc_size);
        gx_vec_omp[i].resize(max_vc_size); gy_vec_omp[i].resize(max_vc_size); gz_vec_omp[i].resize(max_vc_size);
        mu_vec_omp[i].resize(max_vc_size);
        for (int j = 0; j < max_vc_size; ++j)
        {
            vc_S_omp[i][j].resize(9);
        }
    }
    bfv_id.clear(); bfv_id.resize(nc, OpenVolumeMesh::Geometry::Vec3i(-1, -1, -1));
    bvf_id.clear(); bvf_id.resize(nv); boundary_face_number = 0; avg_boundary_edge_length = 0.0;
    int edge_count = 0;
    std::map<std::pair<int, int>, int> edge2id;
    std::map<std::pair<int, int>, int>::iterator edge_it;
    for (size_t i = 0; i < nc; i++)
    {
        if (!tetras[i]->is_onboundary())
        {
            continue;
        }
        int bf_id = -1;
        for (size_t j = 0; j < 4; j++)
        {
            if (!tetras[i]->vertex[j]->boundary)
            {
                bf_id = j;
                break;
            }
        }
        if (bf_id == -1)
        {
            std::cout << "four verts on boundary!" << std::endl;
            return;
        }
        int idx0, idx1, idx2;
        idx0 = tetras[i]->vertex[s_tet_id[bf_id][0]]->id;
        idx1 = tetras[i]->vertex[s_tet_id[bf_id][1]]->id;
        idx2 = tetras[i]->vertex[s_tet_id[bf_id][2]]->id;
        int min_id, mid_id, max_id;
        min_id = std::min(std::min(idx0, idx1), idx2);
        max_id = std::max(std::max(idx0, idx1), idx2);
        mid_id = idx0 + idx1 + idx2 - min_id - max_id;
        std::pair<int, int> edge0(min_id, mid_id), edge1(min_id, max_id), edge2(mid_id, max_id);
        edge_it = edge2id.find(edge0);
        if (edge_it == edge2id.end())
        {
            edge2id.insert(std::pair<std::pair<int, int>, int>(edge0, edge_count));
            edge_count++;
        }
        edge_it = edge2id.find(edge1);
        if (edge_it == edge2id.end())
        {
            edge2id.insert(std::pair<std::pair<int, int>, int>(edge1, edge_count));
            edge_count++;
        }
        edge_it = edge2id.find(edge2);
        if (edge_it == edge2id.end())
        {
            edge2id.insert(std::pair<std::pair<int, int>, int>(edge2, edge_count));
            edge_count++;
        }
    }
    bef_id.clear(); bef_id.resize(edge_count);
    for (size_t i = 0; i < nc; i++)
    {
        if (!tetras[i]->is_onboundary())
        {
            continue;
        }
        int bf_id = -1;
        for (size_t j = 0; j < 4; j++)
        {
            if (!tetras[i]->vertex[j]->boundary)
            {
                bf_id = j;
                break;
            }
        }
        if (bf_id == -1)
        {
            std::cout << "four verts on boundary!" << std::endl;
            return;
        }
        ++boundary_face_number;
        int idx0, idx1, idx2;
        idx0 = tetras[i]->vertex[s_tet_id[bf_id][0]]->id;
        idx1 = tetras[i]->vertex[s_tet_id[bf_id][1]]->id;
        idx2 = tetras[i]->vertex[s_tet_id[bf_id][2]]->id;
        bfv_id[i][0] = idx0;
        bfv_id[i][1] = idx1;
        bfv_id[i][2] = idx2;
        bvf_id[idx0].push_back(i);
        bvf_id[idx1].push_back(i);
        bvf_id[idx2].push_back(i);
        int min_id, mid_id, max_id;
        min_id = std::min(std::min(idx0, idx1), idx2);
        max_id = std::max(std::max(idx0, idx1), idx2);
        mid_id = idx0 + idx1 + idx2 - min_id - max_id;
        std::pair<int, int> edge0(min_id, mid_id), edge1(min_id, max_id), edge2(mid_id, max_id);
        int e_id0 = edge2id[edge0];
        bef_id[e_id0].push_back(i);
        int e_id1 = edge2id[edge1];
        bef_id[e_id1].push_back(i);
        int e_id2 = edge2id[edge2];
        bef_id[e_id2].push_back(i);
    }
    bvv_id.clear(); bvv_id.resize(nv);
    for (int i = 0; i < nv; ++i)
    {
        if (bvf_id[i].size() > 0)
        {
            std::vector<int>& one_ring_f_id = bvf_id[i]; int one_ring_size = one_ring_f_id.size();
            std::vector<int> one_ring_v_id(one_ring_size); std::vector<int> one_ring_v_f_id(one_ring_size);
            std::vector<int> one_ring_order(one_ring_size);
            for (int j = 0; j < one_ring_size; ++j)
            {
                OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[one_ring_f_id[j]];
                for (int k = 0; k < 3; ++k)
                {
                    if (i == fv_id[k])
                    {
                        one_ring_order[j] = k;
                        break;
                    }
                }
            }
            OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[one_ring_f_id[0]];
            int start_v = fv_id[(one_ring_order[0] + 1) % 3];
            int v2 = start_v; one_ring_v_id[0] = v2; one_ring_v_f_id[0] = one_ring_f_id[0];
            int v3 = fv_id[(one_ring_order[0] + 2) % 3];
            for (int j = 1; j < one_ring_size; ++j)
            {
                v2 = v3;
                if (v2 == start_v) break;
                one_ring_v_id[j] = v2;
                for (int k = 1; k < one_ring_size; ++k)
                {
                    OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[one_ring_f_id[k]];
                    if (fv_id[(one_ring_order[k] + 1) % 3] == v2)
                    {
                        v3 = fv_id[(one_ring_order[k] + 2) % 3];
                        one_ring_v_f_id[j] = one_ring_f_id[k];
                        break;
                    }
                }
            }
            bvv_id[i] = one_ring_v_id;
            bvf_id[i] = one_ring_v_f_id;
        }
    }
    avg_boundary_edge_length = 0.0;
    min_boundary_edge_length = -1.0;
    OpenVolumeMesh::Geometry::Vec3d p0, p1, p2; double count = 0.0;
    for (int i = 0; i < bfv_id.size(); ++i)
    {
        OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
        if (one_bfv_id[0] < 0) continue;
        p0[0] = dpx[one_bfv_id[0]]; p0[1] = dpy[one_bfv_id[0]]; p0[2] = dpz[one_bfv_id[0]];
        p1[0] = dpx[one_bfv_id[1]]; p1[1] = dpy[one_bfv_id[1]]; p1[2] = dpz[one_bfv_id[1]];
        p2[0] = dpx[one_bfv_id[2]]; p2[1] = dpy[one_bfv_id[2]]; p2[2] = dpz[one_bfv_id[2]];
        avg_boundary_edge_length += (p0 - p1).norm();
        avg_boundary_edge_length += (p1 - p2).norm();
        avg_boundary_edge_length += (p2 - p0).norm();
        if (min_boundary_edge_length < 0.0)
        {
            min_boundary_edge_length = (p0 - p1).norm();
        }
        if (min_boundary_edge_length > (p0 - p1).norm())
            min_boundary_edge_length = (p0 - p1).norm();
        if (min_boundary_edge_length > (p2 - p1).norm())
            min_boundary_edge_length = (p2 - p1).norm();
        if (min_boundary_edge_length > (p2 - p0).norm())
            min_boundary_edge_length = (p2 - p0).norm();
        count += 3.0;
    }
    avg_boundary_edge_length /= count;
    printf("Average Boundary Edge length : %e\n", avg_boundary_edge_length);
    printf("Min Boundary Edge length : %e\n", min_boundary_edge_length);
    change_big_flag.clear(); last_bfn.clear();
    change_big_flag.resize(nc, -1); last_bfn.resize(nc);
    bound_K = 1.0;
    prepare_ok = true;
}
void polycube_deformation_interface::load_ori_tet(TetStructure<double>* tet_mesh_)
{
    unsigned nv = tet_mesh_->tetra_vertices.size(); unsigned nc = tet_mesh_->tetras.size();
    Eigen::Matrix3d VS, IS;
    cell_S_ori.clear();
    cell_S_ori.resize(nc);
    for (size_t c_id = 0; c_id < nc; ++c_id)
    {
        double cv_count = 0.0;
        TetVertex<double> **tet_vert_ = tet_mesh_->tetras[c_id]->vertex;
        cell_vertex[c_id].clear();
        for (unsigned i = 0; i < 4; ++i)
        {
            int v_id = tet_vert_[i]->id;
            cell_vertex[c_id].push_back(v_id);
        }
        CVec<double, 3> cp = tet_vert_[0]->pos;
        CVec<double, 3> cr = tet_vert_[1]->pos - cp;
        CVec<double, 3> cs = tet_vert_[2]->pos - cp;
        CVec<double, 3> ct = tet_vert_[3]->pos - cp;
        VS(0, 0) = cr[0]; VS(1, 0) = cr[1]; VS(2, 0) = cr[2];
        VS(0, 1) = cs[0]; VS(1, 1) = cs[1]; VS(2, 1) = cs[2];
        VS(0, 2) = ct[0]; VS(1, 2) = ct[1]; VS(2, 2) = ct[2];
        IS = VS.inverse();
        cell_S_ori[c_id] = IS;
    }
}
void polycube_deformation_interface::assign_pos_mesh(VolumeMesh* mesh_, bool r_order)
{
    if (r_order)
    {
        for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
        {
            int v_id = v_it->idx();
            OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(*v_it);
            dpx[v_id] = p[0]; dpy[v_id] = p[1]; dpz[v_id] = p[2];
        }
    }
    else
    {
        for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
        {
            int v_id = v_it->idx();
            OpenVolumeMesh::Geometry::Vec3d p(dpx[v_id], dpy[v_id], dpz[v_id]);
            mesh_->set_vertex(*v_it, p);
        }
    }
}
void polycube_deformation_interface::assign_pos_feature(std::vector<double> &fx, std::vector<double> &fy, std::vector<double> &fz, bool r_order)
{
    assert(!feature_v2e.empty());
    if (r_order)
    {
        if (dpx_feature.empty())
        {
            dpx_feature = dpx;
            dpy_feature = dpy;
            dpz_feature = dpz;
        }
        for (size_t i = 0; i < dpx_feature.size(); i++)
        {
            if (feature_v2e[i].size() != 0)
            {
                dpx_feature[i] = fx[i];
                dpy_feature[i] = fy[i];
                dpz_feature[i] = fz[i];
            }
        }
    }
    else
    {
        for (size_t i = 0; i < dpx_feature.size(); i++)
        {
            if (feature_v2e[i].size() != 0)
            {
                fx[i] = dpx_feature[i];
                fy[i] = dpy_feature[i];
                fz[i] = dpz_feature[i];
            }
        }
    }
}
void polycube_deformation_interface::assign_pos_mesh(TetStructure<double>* tet_mesh_, bool r_order)
{
    if (r_order)
    {
        int nv = tet_mesh_->tetra_vertices.size();
        const std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
        if (dpx.size() == 0)
        {
            dpx.resize(nv);
            dpy.resize(nv);
            dpz.resize(nv);
        }
        for (size_t i = 0; i < nv; i++)
        {
            CVec<double, 3> p = tetra_vertices[i]->pos;
            dpx[i] = p[0]; dpy[i] = p[1]; dpz[i] = p[2];
        }
    }
    else
    {
        int nv = tet_mesh_->tetra_vertices.size();
        std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
        for (size_t i = 0; i < nv; i++)
        {
            CVec<double, 3> p(dpx[i], dpy[i], dpz[i]);
            tetra_vertices[i]->pos = p;
        }
    }
}
void polycube_deformation_interface::compute_distortion(VolumeMesh* mesh_)
{
    Eigen::Matrix3d VS, A; Eigen::Vector3d a;
    double max_cd = 0.0; double min_cd = 1e30; double avg_cd = 0.0; int max_cd_cell_id = -1;
    double max_iso = 0.0; double min_iso = 1e30; double avg_iso = 0.0; int max_iso_cell_id = -1;
    double max_vol = 0.0; double min_vol = 1e30; double avg_vol = 0.0; int max_vol_cell_id = -1;
    double cd_count = 0.0; int flip_count = 0; flipped_cell.clear();
    std::vector<double> all_iso_d; std::vector<double> all_con_d; std::vector<double> all_vol_d;
    int nc = mesh_->n_cells();
    for (int c_id = 0; c_id < nc; ++c_id)
    {
        std::vector<int>& cv_id = cell_vertex[c_id];
        double cpx = dpx[cv_id[0]]; double cpy = dpy[cv_id[0]]; double cpz = dpz[cv_id[0]];
        VS(0, 0) = dpx[cv_id[1]] - cpx; VS(1, 0) = dpy[cv_id[1]] - cpy; VS(2, 0) = dpz[cv_id[1]] - cpz;
        VS(0, 1) = dpx[cv_id[2]] - cpx; VS(1, 1) = dpy[cv_id[2]] - cpy; VS(2, 1) = dpz[cv_id[2]] - cpz;
        VS(0, 2) = dpx[cv_id[3]] - cpx; VS(1, 2) = dpy[cv_id[3]] - cpy; VS(2, 2) = dpz[cv_id[3]] - cpz;
        A = VS * cell_S[c_id];
        if (A.determinant() < 0)
        {
            ++flip_count;
            flipped_cell.push_back(OpenVolumeMesh::CellHandle(c_id));
            all_iso_d.push_back(-1);
            all_con_d.push_back(-1);
            all_vol_d.push_back(-1);
        }
        else
        {
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(A);
            a = svd.singularValues();
            double cd = a(0) / a(2);
            if (cd > max_cd) { max_cd = cd; max_cd_cell_id = c_id; }
            if (cd < min_cd) min_cd = cd;
            avg_cd += cd; cd_count += 1.0;
            all_con_d.push_back(cd);
            double iso_d = a(0) > 1.0 / a(2) ? a(0) : 1.0 / a(2);
            if (iso_d > max_iso) { max_iso = iso_d; max_iso_cell_id = c_id; }
            if (iso_d < min_iso) min_iso = iso_d;
            avg_iso += iso_d;
            all_iso_d.push_back(iso_d);
            double dJ = a(0)*a(1)*a(2);
            double vol_d = dJ > 1.0 / dJ ? dJ : 1.0 / dJ;
            if (vol_d > max_vol) { max_vol = vol_d; max_vol_cell_id = c_id; }
            if (vol_d < min_vol) min_vol = vol_d;
            avg_vol += vol_d;
            all_vol_d.push_back(vol_d);
        }
    }
    avg_cd /= cd_count; avg_iso /= cd_count; avg_vol /= cd_count;
    double std_id = 0.0; double std_cd = 0.0; double std_vd = 0.0;
    for (int i = 0; i < all_iso_d.size(); ++i)
    {
        if (all_iso_d[i] > 0)
        {
            std_cd += (avg_cd - all_con_d[i])*(avg_cd - all_con_d[i]);
            std_id += (avg_iso - all_iso_d[i])*(avg_iso - all_iso_d[i]);
            std_vd += (avg_vol - all_vol_d[i])*(avg_vol - all_vol_d[i]);
        }
    }
    std_id = std::sqrt(std_id / (cd_count - 1)); std_cd = std::sqrt(std_cd / (cd_count - 1)); std_vd = std::sqrt(std_vd / (cd_count - 1));
    printf("---------------------------------------------------------\n");
    printf("Flip Count : %d\n", flip_count);
    printf("Isometric Distortion: %f/%f/%f/%f\n", max_iso, min_iso, avg_iso, std_id);
    printf("conformal Distortion: %f/%f/%f/%f\n", max_cd, min_cd, avg_cd, std_cd);
    printf("volume Distortion: %f/%f/%f/%f\n", max_vol, min_vol, avg_vol, std_vd);
}
void polycube_deformation_interface::update_face_center()
{
    int count = 0; int nc = bfv_id.size(); OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
    face_center.resize(nc);
    for (int i = 0; i < nc; ++i)
    {
        OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
        if (one_bfv_id[0] < 0) continue;
        p0[0] = dpx[one_bfv_id[0]]; p0[1] = dpy[one_bfv_id[0]]; p0[2] = dpz[one_bfv_id[0]];
        p1[0] = dpx[one_bfv_id[1]]; p1[1] = dpy[one_bfv_id[1]]; p1[2] = dpz[one_bfv_id[1]];
        p2[0] = dpx[one_bfv_id[2]]; p2[1] = dpy[one_bfv_id[2]]; p2[2] = dpz[one_bfv_id[2]];
        OpenVolumeMesh::Geometry::Vec3d pt = (p0 + p1 + p2) / 3.0;
        OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0)*0.5;
        face_center[i] = pt;
        ++count;
    }
}
ig::CVec<double, 3> polycube_deformation_interface::compute_distortion(TetStructure<double>* tet_mesh_)
{
    Eigen::Matrix3d VS, A; Eigen::Vector3d a;
    double max_cd = 0.0; double min_cd = 1e30; double avg_cd = 0.0; int max_cd_cell_id = -1;
    double max_iso = 0.0; double min_iso = 1e30; double avg_iso = 0.0; int max_iso_cell_id = -1;
    double max_vol = 0.0; double min_vol = 1e30; double avg_vol = 0.0; int max_vol_cell_id = -1;
    double cd_count = 0.0; int flip_count = 0; flipped_cell.clear();
    std::vector<double> all_iso_d; std::vector<double> all_con_d; std::vector<double> all_vol_d;
    int nc = tet_mesh_->tetras.size();
    for (int c_id = 0; c_id < nc; ++c_id)
    {
        std::vector<int>& cv_id = cell_vertex[c_id];
        double cpx = dpx[cv_id[0]]; double cpy = dpy[cv_id[0]]; double cpz = dpz[cv_id[0]];
        VS(0, 0) = dpx[cv_id[1]] - cpx; VS(1, 0) = dpy[cv_id[1]] - cpy; VS(2, 0) = dpz[cv_id[1]] - cpz;
        VS(0, 1) = dpx[cv_id[2]] - cpx; VS(1, 1) = dpy[cv_id[2]] - cpy; VS(2, 1) = dpz[cv_id[2]] - cpz;
        VS(0, 2) = dpx[cv_id[3]] - cpx; VS(1, 2) = dpy[cv_id[3]] - cpy; VS(2, 2) = dpz[cv_id[3]] - cpz;
        if (cell_S_ori.size() != 0)
            A = VS * cell_S_ori[c_id];
        else
            A = VS * cell_S[c_id];
        if (A.determinant() < 0)
        {
            ++flip_count;
            flipped_cell.push_back(OpenVolumeMesh::CellHandle(c_id));
            all_iso_d.push_back(-1);
            all_con_d.push_back(-1);
            all_vol_d.push_back(-1);
        }
        else
        {
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(A);
            a = svd.singularValues();
            double cd = a(0) / a(2);
            if (cd > max_cd) { max_cd = cd; max_cd_cell_id = c_id; }
            if (cd < min_cd) min_cd = cd;
            avg_cd += cd; cd_count += 1.0;
            all_con_d.push_back(cd);
            double iso_d = a(0) > 1.0 / a(2) ? a(0) : 1.0 / a(2);
            if (iso_d > max_iso) { max_iso = iso_d; max_iso_cell_id = c_id; }
            if (iso_d < min_iso) min_iso = iso_d;
            iso_d = std::min(iso_d, 10.0);
            avg_iso += iso_d;
            all_iso_d.push_back(iso_d);
            double dJ = a(0)*a(1)*a(2);
            double vol_d = dJ > 1.0 / dJ ? dJ : 1.0 / dJ;
            if (vol_d > max_vol) { max_vol = vol_d; max_vol_cell_id = c_id; }
            if (vol_d < min_vol) min_vol = vol_d;
            avg_vol += vol_d;
            all_vol_d.push_back(vol_d);
        }
    }
    avg_cd /= cd_count; avg_iso /= cd_count; avg_vol /= cd_count;
    double std_id = 0.0; double std_cd = 0.0; double std_vd = 0.0;
    for (int i = 0; i < all_iso_d.size(); ++i)
    {
        if (all_iso_d[i] > 0)
        {
            std_cd += (avg_cd - all_con_d[i])*(avg_cd - all_con_d[i]);
            std_id += (avg_iso - all_iso_d[i])*(avg_iso - all_iso_d[i]);
            std_vd += (avg_vol - all_vol_d[i])*(avg_vol - all_vol_d[i]);
        }
    }
    std_id = std::sqrt(std_id / (cd_count - 1)); std_cd = std::sqrt(std_cd / (cd_count - 1)); std_vd = std::sqrt(std_vd / (cd_count - 1));
    printf("---------------------------------------------------------\n");
    printf("Flip Count : %d\n", flip_count);
    printf("Isometric Distortion: %f/%f/%f/%f\n", max_iso, min_iso, avg_iso, std_id);
    printf("conformal Distortion: %f/%f/%f/%f\n", max_cd, min_cd, avg_cd, std_cd);
    printf("volume Distortion: %f/%f/%f/%f\n", max_vol, min_vol, avg_vol, std_vd);
    return ig::CVec<double, 3>(avg_iso, avg_cd, avg_vol);
}
ig::CVec<double, 3> polycube_deformation_interface::compute_volumetric_distortion(TetStructure<double>* tet_mesh_)
{
    Eigen::Matrix3d VS, A; Eigen::Vector3d a;
    double max_cd = 0.0; double min_cd = 1e30; double avg_cd = 0.0; int max_cd_cell_id = -1;
    double max_iso = 0.0; double min_iso = 1e30; double avg_iso = 0.0; int max_iso_cell_id = -1;
    double max_vol = 0.0; double min_vol = 1e30; double avg_vol = 0.0; int max_vol_cell_id = -1;
    double cd_count = 0.0; int flip_count = 0; flipped_cell.clear();
    double volume_count = 0.0;
    std::vector<double> all_iso_d; std::vector<double> all_con_d; std::vector<double> all_vol_d;
    int nc = tet_mesh_->tetras.size();
    for (int c_id = 0; c_id < nc; ++c_id)
    {
        TetVertex<double> **tet_vert_ = tet_mesh_->tetras[c_id]->vertex;
        bool boundary_flag = false;
        for (unsigned i = 0; i < 4; ++i)
        {
            if (tet_vert_[i]->boundary)
                boundary_flag = true;
        }
        std::vector<int>& cv_id = cell_vertex[c_id];
        double cpx = dpx[cv_id[0]]; double cpy = dpy[cv_id[0]]; double cpz = dpz[cv_id[0]];
        VS(0, 0) = dpx[cv_id[1]] - cpx; VS(1, 0) = dpy[cv_id[1]] - cpy; VS(2, 0) = dpz[cv_id[1]] - cpz;
        VS(0, 1) = dpx[cv_id[2]] - cpx; VS(1, 1) = dpy[cv_id[2]] - cpy; VS(2, 1) = dpz[cv_id[2]] - cpz;
        VS(0, 2) = dpx[cv_id[3]] - cpx; VS(1, 2) = dpy[cv_id[3]] - cpy; VS(2, 2) = dpz[cv_id[3]] - cpz;
        if (cell_S_ori.size() != 0)
            A = VS * cell_S_ori[c_id];
        else
            A = VS * cell_S[c_id];
        if (A.determinant() < 0)
        {
            ++flip_count;
            flipped_cell.push_back(OpenVolumeMesh::CellHandle(c_id));
            all_iso_d.push_back(-1);
            all_con_d.push_back(-1);
            all_vol_d.push_back(-1);
        }
        else
        {
            double volume = tet_mesh_->tetras[c_id]->compute_tet_volume();
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(A);
            a = svd.singularValues();
            double cd = a(0) / a(2);
            if (cd > max_cd) { max_cd = cd; max_cd_cell_id = c_id; }
            if (cd < min_cd) min_cd = cd;
            avg_cd += volume * cd;
            cd_count += 1.0;
            all_con_d.push_back(volume * cd);
            double iso_d = a(0) > 1.0 / a(2) ? a(0) : 1.0 / a(2);
            if (iso_d > max_iso) { max_iso = iso_d; max_iso_cell_id = c_id; }
            if (iso_d < min_iso) min_iso = iso_d;
            avg_iso += volume * iso_d;
            all_iso_d.push_back(volume * iso_d);
            volume_count += volume;
            double dJ = a(0)*a(1)*a(2);
            double vol_d = dJ > 1.0 / dJ ? dJ : 1.0 / dJ;
            if (vol_d > max_vol) { max_vol = vol_d; max_vol_cell_id = c_id; }
            if (vol_d < min_vol) min_vol = vol_d;
            avg_vol += volume * vol_d;
            all_vol_d.push_back(volume * vol_d);
        }
    }
    avg_cd /= volume_count;
    avg_iso /= volume_count;
    avg_vol /= volume_count;
    double std_id = 0.0; double std_cd = 0.0; double std_vd = 0.0;
    for (int i = 0; i < all_iso_d.size(); ++i)
    {
        if (all_iso_d[i] > 0)
        {
            std_cd += (avg_cd - all_con_d[i])*(avg_cd - all_con_d[i]);
            std_id += (avg_iso - all_iso_d[i])*(avg_iso - all_iso_d[i]);
            std_vd += (avg_vol - all_vol_d[i])*(avg_vol - all_vol_d[i]);
        }
    }
    std_id = std::sqrt(std_id / (cd_count - 1)); std_cd = std::sqrt(std_cd / (cd_count - 1)); std_vd = std::sqrt(std_vd / (cd_count - 1));
    printf("---------------------------------------------------------\n");
    printf("Volumetric weighted Isometric Distortion: %f/%f/%f/%f\n", max_iso, min_iso, avg_iso, std_id);
    printf("Volumetric weighted conformal Distortion: %f/%f/%f/%f\n", max_cd, min_cd, avg_cd, std_cd);
    printf("Volumetric weighted volume Distortion: %f/%f/%f/%f\n", max_vol, min_vol, avg_vol, std_vd);
    return ig::CVec<double, 3>(avg_iso, avg_cd, avg_vol);
}
double polycube_deformation_interface::compute_energy_ratio_filter_deform()
{
    int nc = cell_S.size(); int nv = bvv_id.size();
    int max_num_t = omp_get_max_threads();
    std::vector<Eigen::Matrix3d> VS(max_num_t), A(max_num_t), A2(max_num_t);
    std::vector<double> distortion_v(nc, -1.0);
#pragma omp parallel for
    for (int c_id = 0; c_id < nc; ++c_id)
    {
        int id = omp_get_thread_num();
        const int* cv_id = cell_vertex[c_id].data();
        int v_id_0 = cv_id[0]; int v_id_1 = cv_id[1]; int v_id_2 = cv_id[2]; int v_id_3 = cv_id[3];
        double x0 = dpx[v_id_0]; double y0 = dpy[v_id_0]; double z0 = dpz[v_id_0];
        VS[id](0, 0) = dpx[v_id_1] - x0; VS[id](1, 0) = dpy[v_id_1] - y0; VS[id](2, 0) = dpz[v_id_1] - z0;
        VS[id](0, 1) = dpx[v_id_2] - x0; VS[id](1, 1) = dpy[v_id_2] - y0; VS[id](2, 1) = dpz[v_id_2] - z0;
        VS[id](0, 2) = dpx[v_id_3] - x0; VS[id](1, 2) = dpy[v_id_3] - y0; VS[id](2, 2) = dpz[v_id_3] - z0;
        A[id] = VS[id] * cell_S[c_id];
        A2[id] = A[id].transpose()*A[id];
        double AF = A2[id](0, 0) + A2[id](1, 1) + A2[id](2, 2);
        double AF_05 = std::sqrt(AF);
        double AF2 = A2[id].squaredNorm();
        double AF_I = (AF*AF - AF2)*0.5; double AF_I_05 = std::sqrt(AF_I);
        double det_A = A[id].determinant();
        double i_det_A = 1.0 / det_A;
        double g = AF_05 * AF_I_05;
        double e = g * i_det_A;
        double k = ((e*e - 1.0)*0.125 + 0.5*(det_A + i_det_A)) * 0.5;
        if (k > 60) k = 60;
        distortion_v[c_id] = std::exp(k);
    }
    if (distortion_big_count.size() != nv)
    {
        distortion_big_count.clear();
        distortion_big_count.resize(nv, 0);
    }
    if (distortion_big_cell.size() != nc)
    {
        distortion_big_cell.clear();
        distortion_big_cell.resize(nc, -1);
    }
    double large_dis_th = std::exp(bound_K);
    double e_sum = 0.0;
    double n_diff = 0.0; double count_diff = 0.0; local_energy_ratio.clear(); local_energy_ratio.resize(nc, -1);
    OpenVolumeMesh::Geometry::Vec3d p0, p1, p2; double max_ratio = 0.0; double ratio = 0.0;
    std::vector<double> n_diff_v(nc, 10.0);
    for (int i = 0; i < nc; ++i)
    {
        OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[i];
        if (fv_id[0] >= 0)
        {
            OpenVolumeMesh::Geometry::Vec3d& ln = target_bfn[i];
            p0[0] = dpx[fv_id[0]]; p0[1] = dpy[fv_id[0]]; p0[2] = dpz[fv_id[0]];
            p1[0] = dpx[fv_id[1]]; p1[1] = dpy[fv_id[1]]; p1[2] = dpz[fv_id[1]];
            p2[0] = dpx[fv_id[2]]; p2[1] = dpy[fv_id[2]]; p2[2] = dpz[fv_id[2]];
            OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0).normalize();
            n_diff += 0.5*(n - ln).sqrnorm();
            count_diff += 1.0;
            n_diff_v[i] = 0.5*(n - ln).sqrnorm();
            e_sum += distortion_v[i];
            local_energy_ratio[i] = distortion_v[i] / (1e-10 + 0.5*(n - ln).sqrnorm());
            if (local_energy_ratio[i] > 1e16) local_energy_ratio[i] = 1e16;
        }
        if (distortion_v[i] > large_dis_th)
        {
            distortion_big_cell[i] = 1;
            const int* cv_id = cell_vertex[i].data();
            for (int j = 0; j < 4; ++j)
            {
                distortion_big_count[cv_id[j]] += 1;
            }
        }
    }
    n_diff /= count_diff;
    e_sum /= count_diff;
    double s = e_sum / (n_diff + 1e-8);
    printf("1 %f %f %f\n", e_sum, n_diff, s);
    double up_diff_th = 2.0*n_diff; double down_diff_th = 0.1*n_diff;
    double avg_n_diff = 0.0; count_diff = 0.0; e_sum = 0.0;
    for (int i = 0; i < nc; ++i)
    {
        if (n_diff_v[i] < up_diff_th && n_diff_v[i] > down_diff_th)
        {
            avg_n_diff += n_diff_v[i]; count_diff += 1.0;
            e_sum += distortion_v[i];
        }
    }
    if (std::abs(count_diff) < 1e-8)
    {
        s = 1.0;
    }
    else
    {
        avg_n_diff /= count_diff;
        e_sum /= count_diff;
        s = e_sum / (avg_n_diff + 1e-8);
    }
    printf("2 %f %f %f\n", e_sum, avg_n_diff, s);
    printf("--------------------------------------------\n");
    if (s > 1e16) s = 1e16;
    for (int i = 0; i < nc; ++i)
    {
        if (local_energy_ratio[i] > 2 * s)
        {
            local_energy_ratio[i] = 2 * s;
        }
        else if (local_energy_ratio[i] < 0.1* s)
        {
            local_energy_ratio[i] = s * 0.1;
        }
    }
    return s;
}
std::pair<double, double> polycube_deformation_interface::compute_chart_normal_diff()
{
    OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
    double max_max, max_avg;
    max_max = -1;
    max_avg = -1;
    for (size_t i = 0; i < chart_tet_idx.size(); i++)
    {
        double max_diff = -1.0;
        double avg_diff = 0.0;
        for (size_t j = 0; j < chart_tet_idx[i].size(); j++)
        {
            OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[chart_tet_idx[i][j]];
            assert(fv_id[0] >= 0);
            if (fv_id[0] >= 0)
            {
                OpenVolumeMesh::Geometry::Vec3d& ln = target_bfn[chart_tet_idx[i][j]];
                p0[0] = dpx[fv_id[0]]; p0[1] = dpy[fv_id[0]]; p0[2] = dpz[fv_id[0]];
                p1[0] = dpx[fv_id[1]]; p1[1] = dpy[fv_id[1]]; p1[2] = dpz[fv_id[1]];
                p2[0] = dpx[fv_id[2]]; p2[1] = dpy[fv_id[2]]; p2[2] = dpz[fv_id[2]];
                OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0).normalize();
                avg_diff += (n - ln).sqrnorm();
                if (max_diff < (n - ln).sqrnorm())
                    max_diff = (n - ln).sqrnorm();
            }
        }
        avg_diff = avg_diff / chart_tet_idx[i].size();
        if (max_avg < avg_diff)
            max_avg = avg_diff;
        if (max_max < max_diff)
            max_max = max_diff;
    }
    return (std::pair<double, double>(max_avg, max_max));
}
void polycube_deformation_interface::set_vertex_color(int n_color, const std::vector<int>& v_color, int nv)
{
    number_of_color = n_color;
    vertex_diff_color.clear();
    if (number_of_color == 0)
    {
        vertex_diff_color.resize(nv);
    }
    else
    {
        vertex_diff_color.resize(number_of_color);
    }
    for (unsigned i = 0; i < nv; ++i)
    {
        if (number_of_color == 0)
        {
            vertex_diff_color[i].push_back(i);
        }
        else
        {
            vertex_diff_color[v_color[i] - 1].push_back(i);
        }
    }
}
void polycube_deformation_interface::set_feature_edges_ovm(const std::vector<std::pair<int, int>> &edge_pairs, VolumeMesh * mesh_)
{
    if (mesh_ == NULL) return;
    feature_edge_flag.clear();
    feature_edge_flag.resize(mesh_->n_edges(), false);
    feature_e2v.clear();
    corner_vert_flag.clear();
    for (size_t i = 0; i < edge_pairs.size(); i++)
    {
        int id0 = edge_pairs[i].first;
        int id1 = edge_pairs[i].second;
        int edge_id = -1;
        OpenVolumeMesh::VertexHandle vh = OpenVolumeMesh::VertexHandle(id0);
        for (auto it = mesh_->voh_iter(vh); it; ++it)
        {
            OpenVolumeMesh::EdgeHandle eh = mesh_->edge_handle(*it);
            OpenVolumeMesh::OpenVolumeMeshEdge ovme = mesh_->edge(eh);
            if (ovme.from_vertex().idx() == id1 || ovme.to_vertex().idx() == id1)
            {
                edge_id = eh.idx();
                break;
            }
        }
        assert(edge_id != -1);
        feature_edge_flag[edge_id] = true;
    }
    feature_edge_change_flag = true;
}
void polycube_deformation_interface::build_feature_edges_connectivity(VolumeMesh * mesh_, bool refine_long_edge, int min_edge_size)
{
    if (feature_edge_change_flag == false) return;
    if (feature_edge_flag.empty()) return;
    set_feature_edge_array();
    feature_edge_length.clear();
    corner_vert_flag.resize(mesh_->n_vertices(), false);
    feature_v2e.clear();
    feature_v2e.resize(mesh_->n_vertices());
    feature_v2v.clear();
    feature_v2v.resize(mesh_->n_vertices());
    feature_neighbor_vert2cellpair.clear();
    feature_neighbor_vert2cellpair.resize(mesh_->n_vertices());
    for (size_t i = 0; i < feature_edge_flag.size(); i++)
    {
        if (feature_edge_flag[i])
        {
            OpenVolumeMesh::EdgeHandle eh = OpenVolumeMesh::EdgeHandle(i);
            OpenVolumeMesh::OpenVolumeMeshEdge ovme = mesh_->edge(eh);
            int id0 = ovme.from_vertex().idx();
            int id1 = ovme.to_vertex().idx();
            feature_v2e[id0].push_back(i);
            feature_v2e[id1].push_back(i);
            feature_v2v[id0].push_back(id1);
            feature_v2v[id1].push_back(id0);
            std::set<int> twoneighborface;
            OpenVolumeMesh::HalfEdgeHandle heh = mesh_->halfedge_handle(eh, 0);
            for (auto it = mesh_->hehf_iter(heh); it; ++it)
            {
                OpenVolumeMesh::FaceHandle fh = mesh_->face_handle(*it);
                if (mesh_->is_boundary(fh))
                {
                    twoneighborface.insert(fh.idx());
                }
            }
            assert(twoneighborface.size() == 2);
            std::vector<int> face;
            for (auto it : twoneighborface)
            {
                face.push_back(it);
            }
            std::vector<int> twocell;
            int neighbor_vert[2] = { -1, -1 };
            for (size_t j = 0; j < 2; j++)
            {
                OpenVolumeMesh::FaceHandle fh = OpenVolumeMesh::FaceHandle(face[j]);
                OpenVolumeMesh::HalfFaceHandle hfh = mesh_->halfface_handle(fh, 0);
                OpenVolumeMesh::CellHandle ch = mesh_->incident_cell(hfh);
                if (ch.is_valid())
                    twocell.push_back(ch.idx());
                else
                {
                    OpenVolumeMesh::HalfFaceHandle hfh1 = mesh_->halfface_handle(fh, 1);
                    OpenVolumeMesh::CellHandle ch1 = mesh_->incident_cell(hfh1);
                    assert(ch1.is_valid());
                    twocell.push_back(ch1.idx());
                }
                for (auto it = mesh_->hfv_iter(hfh); it; ++it)
                {
                    int idx = (*it).idx();
                    if (idx != id0 && idx != id1)
                    {
                        neighbor_vert[j] = idx;
                        break;
                    }
                }
            }
            assert(neighbor_vert[0] != -1 && neighbor_vert[1] != -1);
            assert(twocell.size() == 2);
            feature_neighbor_vert2cellpair[neighbor_vert[0]].push_back(std::pair<int, int>(twocell[0], twocell[1]));
            feature_neighbor_vert2cellpair[neighbor_vert[1]].push_back(std::pair<int, int>(twocell[1], twocell[0]));
        }
    }
    std::vector<int> feature_edges;
    std::vector<std::pair<int, int>> feature_e2v_backup = feature_e2v;
    feature_e2v.clear();
    feature_e2v.resize(mesh_->n_edges(), std::pair<int, int>(-1, -1));
    feature_neighbor.clear();
    feature_neighbor.resize(mesh_->n_edges());
    for (size_t i = 0; i < feature_edge_flag.size(); i++)
    {
        if (feature_edge_flag[i] == 1)
        {
            feature_edges.push_back(i);
            OpenVolumeMesh::EdgeHandle eh = OpenVolumeMesh::EdgeHandle(i);
            OpenVolumeMesh::OpenVolumeMeshEdge ovme = mesh_->edge(eh);
            feature_e2v[i].first = ovme.from_vertex().idx();
            feature_e2v[i].second = ovme.to_vertex().idx();
        }
    }
    for (size_t i = 0; i < feature_edges.size(); i++)
    {
        int eid = feature_edges[i];
        int v[2];
        v[0] = feature_e2v[eid].first;
        v[1] = feature_e2v[eid].second;
        for (size_t j = 0; j < 2; j++)
        {
            if (feature_v2e[v[j]].size() == 2)
            {
                int othere = feature_v2e[v[j]][0];
                if (othere == eid) othere = feature_v2e[v[j]][1];
                feature_neighbor[eid].push_back(othere);
            }
        }
    }
    long_feature_edge.clear();
    feature_edge2long.clear();
    feature_edge2long.resize(mesh_->n_edges(), -1);
    std::set<int> feature_edge_set(feature_edges.begin(), feature_edges.end());
    while (!feature_edge_set.empty())
    {
        std::vector<int> one_long_edge;
        auto it = feature_edge_set.begin();
        std::queue<int> q;
        q.push(*it);
        one_long_edge.push_back(*it);
        std::vector<int> edge_color(mesh_->n_edges(), 0);
        edge_color[*it] = 1;
        while (!q.empty())
        {
            int front = q.front();
            q.pop();
            for (auto n : feature_neighbor[front])
            {
                if (edge_color[n] != 1)
                {
                    edge_color[n] = 1;
                    one_long_edge.push_back(n);
                    q.push(n);
                }
            }
        }
        long_feature_edge.push_back(one_long_edge);
        for (auto e : one_long_edge)
        {
            feature_edge_set.erase(e);
            feature_edge2long[e] = long_feature_edge.size() - 1;
        }
    }
    if (feature_e2v_backup.empty())
    {
        for (size_t i = 0; i < long_feature_edge.size(); i++)
        {
            std::vector<int> edge_color(mesh_->n_edges(), 0);
            std::queue<int> q;
            q.push(long_feature_edge[i][0]);
            edge_color[long_feature_edge[i][0]] = 1;
            while (!q.empty())
            {
                int front = q.front(); q.pop();
                int first = feature_e2v[front].first, second = feature_e2v[front].second;
                if (feature_v2e[first].size() == 2)
                {
                    int othere = feature_v2e[first][0];
                    if (othere == front) othere = feature_v2e[first][1];
                    if (edge_color[othere] == 0)
                    {
                        edge_color[othere] = 1;
                        if (feature_e2v[othere].first == first)
                        {
                            int tmp = feature_e2v[othere].first;
                            feature_e2v[othere].first = feature_e2v[othere].second;
                            feature_e2v[othere].second = tmp;
                        }
                        q.push(othere);
                    }
                }
                if (feature_v2e[second].size() == 2)
                {
                    int othere = feature_v2e[second][0];
                    if (othere == front) othere = feature_v2e[second][1];
                    if (edge_color[othere] == 0)
                    {
                        edge_color[othere] = 1;
                        if (feature_e2v[othere].second == second)
                        {
                            int tmp = feature_e2v[othere].first;
                            feature_e2v[othere].first = feature_e2v[othere].second;
                            feature_e2v[othere].second = tmp;
                        }
                        q.push(othere);
                    }
                }
            }
        }
    }
    else
    {
        feature_e2v = feature_e2v_backup;
    }
    for (size_t i = 0; i < feature_v2e.size(); i++)
    {
        if (feature_v2e[i].size() == 1 || feature_v2e[i].size() >= 3)
            corner_vert_flag[i] = true;
    }
    if (refine_long_edge)
    {
        refine_corner_and_long_edge(min_edge_size);
    }
    feature_edge_change_flag = false;
}
void polycube_deformation_interface::refine_corner_and_long_edge(int min_edge_size)
{
    filter_feature_edge_ring(2, 1.0, false);
    const double dihedral_angle = 100 * 3.1415926535897932384626433832795 / 180;
    for (int i = 0; i < feature_v2e.size(); i++)
    {
        if (feature_v2e[i].size() == 2)
        {
            double center_x = dpx[i];
            double center_y = dpy[i];
            double center_z = dpz[i];
            int eid0 = feature_v2e[i][0];
            int eid1 = feature_v2e[i][1];
            int otherid0 = feature_e2v[eid0].first;
            if (otherid0 == i) otherid0 = feature_e2v[eid0].second;
            int otherid1 = feature_e2v[eid1].first;
            if (otherid1 == i) otherid1 = feature_e2v[eid1].second;
            double x0 = dpx[otherid0];
            double y0 = dpy[otherid0];
            double z0 = dpz[otherid0];
            double x1 = dpx[otherid1];
            double y1 = dpy[otherid1];
            double z1 = dpz[otherid1];
            OpenVolumeMesh::Geometry::Vec3d n0(x0 - center_x, y0 - center_y, z0 - center_z);
            OpenVolumeMesh::Geometry::Vec3d n1(x1 - center_x, y1 - center_y, z1 - center_z);
            n0.normalize();
            n1.normalize();
            double tmp_dot = OpenVolumeMesh::Geometry::dot(n0, n1);
            if (acos(tmp_dot) < dihedral_angle) corner_vert_flag[i] = true;
        }
    }
    feature_long_edge_redefinition();
}
void polycube_deformation_interface::feature_long_edge_redefinition(int min_edge_size)
{
    long_feature_edge.clear();
    feature_edge2long.clear();
    feature_edge2long.resize(bef_id.size(), -1);
    std::set<int> feature_edge_set(feature_edge_array.begin(), feature_edge_array.end());
    while (!feature_edge_set.empty())
    {
        std::vector<int> one_long_edge;
        auto it = feature_edge_set.begin();
        std::queue<int> q;
        q.push(*it);
        one_long_edge.push_back(*it);
        std::vector<int> edge_color(bef_id.size(), 0);
        edge_color[*it] = 1;
        while (!q.empty())
        {
            int front = q.front();
            q.pop();
            int v0 = feature_e2v[front].first;
            int v1 = feature_e2v[front].second;
            if (corner_vert_flag[v0] != true)
            {
                assert(feature_v2e[v0].size() == 2);
                int othere = feature_v2e[v0][0];
                if (othere == front) othere = feature_v2e[v0][1];
                if (edge_color[othere] != 1)
                {
                    edge_color[othere] = 1;
                    one_long_edge.push_back(othere);
                    q.push(othere);
                }
            }
            if (corner_vert_flag[v1] != true)
            {
                assert(feature_v2e[v1].size() == 2);
                int othere = feature_v2e[v1][0];
                if (othere == front) othere = feature_v2e[v1][1];
                if (edge_color[othere] != 1)
                {
                    edge_color[othere] = 1;
                    one_long_edge.push_back(othere);
                    q.push(othere);
                }
            }
        }
        long_feature_edge.push_back(one_long_edge);
        for (auto e : one_long_edge)
        {
            feature_edge_set.erase(e);
            feature_edge2long[e] = long_feature_edge.size() - 1;
        }
    }
    for (size_t i = 0; i < long_feature_edge.size(); i++)
    {
        std::vector<std::vector<int>> local_vv(bvf_id.size());
        int start = -1;
        for (size_t j = 0; j < long_feature_edge[i].size(); j++)
        {
            int eid = long_feature_edge[i][j];
            int v0 = feature_e2v[eid].first;
            int v1 = feature_e2v[eid].second;
            local_vv[v0].push_back(v1);
            local_vv[v1].push_back(v0);
            if (corner_vert_flag[v0]) start = v0;
            if (corner_vert_flag[v1]) start = v1;
        }
        if (start == -1)
        {
            for (size_t j = 0; j < long_feature_edge[i].size(); j++)
            {
                int eid = long_feature_edge[i][j];
                int v0 = feature_e2v[eid].first;
                int v1 = feature_e2v[eid].second;
                if (local_vv[v0].size() == 2)
                {
                    int eid0 = feature_v2e[v0][0];
                    int eid1 = feature_v2e[v0][1];
                    OpenVolumeMesh::Geometry::Vec3d tn0 = target_fea_n[eid0];
                    OpenVolumeMesh::Geometry::Vec3d tn1 = target_fea_n[eid1];
                    if (tn0 != tn1)
                    {
                        start = v0;
                        corner_vert_flag[start] = true;
                        break;
                    }
                }
                if (local_vv[v1].size() == 2)
                {
                    int eid0 = feature_v2e[v1][0];
                    int eid1 = feature_v2e[v1][1];
                    OpenVolumeMesh::Geometry::Vec3d tn0 = target_fea_n[eid0];
                    OpenVolumeMesh::Geometry::Vec3d tn1 = target_fea_n[eid1];
                    if (tn0 != tn1)
                    {
                        start = v1;
                        corner_vert_flag[start] = true;
                        break;
                    }
                }
            }
        }
        if (start == -1)
        {
            start = feature_e2v[long_feature_edge[i][0]].first;
        }
        int count = 1;
        std::vector<int> order2v;
        order2v.push_back(start);
        std::vector<int> v_color(bvf_id.size(), -1);
        v_color[start] = 1;
        bool find_flag = true;
        int nextid = -1;
        while (find_flag)
        {
            if (nextid != -1)
            {
                start = nextid;
                nextid = -1;
                v_color[start] = 1;
                order2v.push_back(start);
                count++;
            }
            for (size_t j = 0; j < local_vv[start].size(); j++)
            {
                if (v_color[local_vv[start][j]] != 1)
                {
                    nextid = local_vv[start][j];
                    break;
                }
            }
            if (nextid == -1) find_flag = false;
        }
        int prev_id = 0;
        for (int j = 1; j < count - 1; j++)
        {
            int vid = order2v[j];
            int dist0 = abs(j - prev_id);
            int dist1 = abs(count - 1 - j);
            if (dist0 >= min_edge_size && dist1 >= min_edge_size)
            {
                int eid0 = feature_v2e[vid][0];
                int eid1 = feature_v2e[vid][1];
                OpenVolumeMesh::Geometry::Vec3d tn0 = target_fea_n[eid0];
                OpenVolumeMesh::Geometry::Vec3d tn1 = target_fea_n[eid1];
                if (tn0 != tn1)
                {
                    corner_vert_flag[vid] = true;
                    prev_id = j;
                }
            }
        }
    }
    long_feature_edge.clear();
    feature_edge2long.clear();
    feature_edge2long.resize(bef_id.size(), -1);
    feature_edge_set.clear();
    for (auto it : feature_edge_array) feature_edge_set.insert(it);
    while (!feature_edge_set.empty())
    {
        std::vector<int> one_long_edge;
        auto it = feature_edge_set.begin();
        std::queue<int> q;
        q.push(*it);
        one_long_edge.push_back(*it);
        std::vector<int> edge_color(bef_id.size(), 0);
        edge_color[*it] = 1;
        while (!q.empty())
        {
            int front = q.front();
            q.pop();
            int v0 = feature_e2v[front].first;
            int v1 = feature_e2v[front].second;
            if (corner_vert_flag[v0] != true)
            {
                assert(feature_v2e[v0].size() == 2);
                int othere = feature_v2e[v0][0];
                if (othere == front) othere = feature_v2e[v0][1];
                if (edge_color[othere] != 1)
                {
                    edge_color[othere] = 1;
                    one_long_edge.push_back(othere);
                    q.push(othere);
                }
            }
            if (corner_vert_flag[v1] != true)
            {
                assert(feature_v2e[v1].size() == 2);
                int othere = feature_v2e[v1][0];
                if (othere == front) othere = feature_v2e[v1][1];
                if (edge_color[othere] != 1)
                {
                    edge_color[othere] = 1;
                    one_long_edge.push_back(othere);
                    q.push(othere);
                }
            }
        }
        long_feature_edge.push_back(one_long_edge);
        for (auto e : one_long_edge)
        {
            feature_edge_set.erase(e);
            feature_edge2long[e] = long_feature_edge.size() - 1;
        }
    }
}
void polycube_deformation_interface::corner_redefinition_polyline()
{
    get_polyline_normal();
    corner_vert_flag.clear();
    corner_vert_flag.resize(bvf_id.size(), false);
    for (size_t i = 0; i < feature_v2e.size(); i++)
    {
        if (feature_v2e[i].size() == 1 || feature_v2e[i].size() >= 3)
            corner_vert_flag[i] = true;
        else if (feature_v2e[i].size() == 2)
        {
            OpenVolumeMesh::Geometry::Vec3d n0 = polyline_target_n[feature_v2e[i][0]];
            OpenVolumeMesh::Geometry::Vec3d n1 = polyline_target_n[feature_v2e[i][1]];
            double dot = OpenVolumeMesh::Geometry::dot(n0, n1);
            if (dot < 0.1)
            {
                corner_vert_flag[i] = true;
            }
        }
    }
}
void polycube_deformation_interface::repair_features(VolumeMesh * mesh_, bool polyline_guidance)
{
    int n_level = 2;
    double ss = 1.0;
    double sr = 0.3;
    if (!polyline_guidance)
    {
        bool feature_OK = false;
        build_feature_edges_connectivity(mesh_);
        filter_boundary_normal_feature(ss, sr, 3);
        filter_feature_edge_ring(n_level, ss);
        feature_OK = check_and_repair_feature(false);
        while (!feature_OK)
        {
            check_and_repair_feature(true);
            build_feature_edges_connectivity(mesh_);
            filter_boundary_normal_feature(ss, sr, 3);
            filter_feature_edge_ring(n_level, ss);
            feature_OK = check_and_repair_feature(false);
        }
    }
    else
    {
        bool feature_OK = false;
        build_feature_edges_connectivity(mesh_);
        feature_long_edge_redefinition();
        if (polyline_target_n.empty()) get_polyline_normal();
        target_fea_n = polyline_target_n;
        feature_OK = check_and_repair_feature(false);
        while (!feature_OK)
        {
            check_and_repair_feature(true);
            build_feature_edges_connectivity(mesh_);
            feature_long_edge_redefinition();
            feature_OK = check_and_repair_feature(false);
        }
    }
}
bool polycube_deformation_interface::check_and_repair_feature(bool repair_status, bool polyline_guidance)
{
    bool feature_OK = true;
    std::set<int> deleted_long_id;
    std::map<OpenVolumeMesh::Geometry::Vec3d, int> normal2idx;
    normal2idx[OpenVolumeMesh::Geometry::Vec3d(1.0, 0.0, 0.0)] = 0;
    normal2idx[OpenVolumeMesh::Geometry::Vec3d(-1.0, 0.0, 0.0)] = 1;
    normal2idx[OpenVolumeMesh::Geometry::Vec3d(0.0, 1.0, 0.0)] = 2;
    normal2idx[OpenVolumeMesh::Geometry::Vec3d(0.0, -1.0, 0.0)] = 3;
    normal2idx[OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, 1.0)] = 4;
    normal2idx[OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, -1.0)] = 5;
    std::vector<double> long_edge_distortion(long_feature_edge.size(), 0.0);
    for (size_t i = 0; i < long_feature_edge.size(); i++)
    {
        int one_edge_size = long_feature_edge[i].size();
        double sum = 0;
        for (size_t j = 0; j < one_edge_size; j++)
        {
            int eid = long_feature_edge[i][j];
            int v0 = feature_e2v[eid].first, v1 = feature_e2v[eid].second;
            OpenVolumeMesh::Geometry::Vec3d p0(dpx[v0], dpy[v0], dpz[v0]);
            OpenVolumeMesh::Geometry::Vec3d p1(dpx[v1], dpy[v1], dpz[v1]);
            OpenVolumeMesh::Geometry::Vec3d tmpn = (p1 - p0).normalize();
            sum += (tmpn - target_fea_n[eid]).length();
        }
        sum = sum / one_edge_size;
        long_edge_distortion[i] = sum;
    }
    std::set<int> allaxis;
    std::set<int> alln;
    for (int i = 0; i < feature_v2e.size(); i++)
    {
        switch (feature_v2e[i].size())
        {
        case 6:
        {
            const std::vector<int> &eid6 = feature_v2e[i];
            alln.clear();
            for (size_t j = 0; j < eid6.size(); j++)
            {
                OpenVolumeMesh::Geometry::Vec3d tar_n = target_fea_n[eid6[j]];
                int first = feature_e2v[eid6[j]].first;
                if (first != i) tar_n = -tar_n;
                int tmpn = normal2idx[tar_n];
                alln.insert(tmpn);
            }
            if (alln.size() == 6)
            {
                break;
            }
            ;
        }
        case 5:
            ;
        case 4:
        {
            feature_OK = false;
            int maxid = -1;
            double maxdist = 0.0;
            bool already_delete_flag = false;
            for (size_t j = 0; j < feature_v2e[i].size(); j++)
            {
                int eid = feature_v2e[i][j];
                int long_edge_id = feature_edge2long[eid];
                auto it = deleted_long_id.find(long_edge_id);
                if (it != deleted_long_id.end())
                    already_delete_flag = true;
            }
            if (already_delete_flag)
                continue;
            for (size_t j = 0; j < feature_v2e[i].size(); j++)
            {
                int eid = feature_v2e[i][j];
                double tmpdist = long_edge_distortion[feature_edge2long[eid]];
                if (tmpdist > maxdist)
                {
                    maxdist = tmpdist;
                    maxid = j;
                }
            }
            assert(maxid != -1);
            deleted_long_id.insert(feature_edge2long[feature_v2e[i][maxid]]);
            break;
        }
        case 3:
        {
            const std::vector<int> &eid3 = feature_v2e[i];
            int n[3], axis[3];
            allaxis.clear();
            alln.clear();
            for (size_t j = 0; j < 3; j++)
            {
                OpenVolumeMesh::Geometry::Vec3d tar_n = target_fea_n[eid3[j]];
                int first = feature_e2v[eid3[j]].first;
                if (first != i) tar_n = -tar_n;
                n[j] = normal2idx[tar_n];
                axis[j] = n[j] / 2;
                allaxis.insert(axis[j]);
                alln.insert(n[j]);
            }
            if (allaxis.size() == 1)
            {
                feature_OK = false;
                assert(alln.size() <= 2);
                if (alln.size() == 1)
                {
                    bool already_delete_flag = false;
                    for (size_t j = 0; j < feature_v2e[i].size(); j++)
                    {
                        int eid = feature_v2e[i][j];
                        int long_edge_id = feature_edge2long[eid];
                        auto it = deleted_long_id.find(long_edge_id);
                        if (it != deleted_long_id.end())
                            already_delete_flag = true;
                    }
                    if (already_delete_flag)
                        continue;
                    int minid = 0;
                    double mindist = long_edge_distortion[feature_edge2long[eid3[0]]];
                    for (size_t j = 1; j < 3; j++)
                    {
                        double tmp_dist = long_edge_distortion[feature_edge2long[eid3[j]]];
                        if (tmp_dist < mindist)
                        {
                            mindist = tmp_dist;
                            minid = j;
                        }
                    }
                    for (int j = 0; j < 3; j++)
                    {
                        if (j == minid) continue;
                        deleted_long_id.insert(feature_edge2long[eid3[j]]);
                    }
                }
                else
                {
                    std::vector<int> twosamedirid;
                    for (size_t j = 0; j < 3; j++)
                    {
                        for (size_t k = j + 1; k < 3; k++)
                        {
                            if (n[j] == n[k])
                            {
                                twosamedirid.push_back(j);
                                twosamedirid.push_back(k);
                                break;
                            }
                        }
                    }
                    assert(twosamedirid.size() == 2);
                    int le0 = feature_edge2long[eid3[twosamedirid[0]]];
                    int le1 = feature_edge2long[eid3[twosamedirid[1]]];
                    bool already_delete_flag = false;
                    auto it = deleted_long_id.find(le0);
                    if (it != deleted_long_id.end())
                        continue;
                    it = deleted_long_id.find(le1);
                    if (it != deleted_long_id.end())
                        continue;
                    double dist0 = long_edge_distortion[feature_edge2long[eid3[twosamedirid[0]]]];
                    double dist1 = long_edge_distortion[feature_edge2long[eid3[twosamedirid[1]]]];
                    if (dist0 > dist1)
                    {
                        deleted_long_id.insert(feature_edge2long[eid3[twosamedirid[0]]]);
                    }
                    else
                    {
                        deleted_long_id.insert(feature_edge2long[eid3[twosamedirid[1]]]);
                    }
                }
            }
            else if (allaxis.size() == 2)
            {
                feature_OK = false;
                assert(alln.size() != 1);
                if (alln.size() == 2)
                {
                    std::vector<int> twosamedirid;
                    for (size_t j = 0; j < 3; j++)
                    {
                        for (size_t k = j + 1; k < 3; k++)
                        {
                            if (n[j] == n[k])
                            {
                                twosamedirid.push_back(j);
                                twosamedirid.push_back(k);
                                break;
                            }
                        }
                    }
                    assert(twosamedirid.size() == 2);
                    int le0 = feature_edge2long[eid3[twosamedirid[0]]];
                    int le1 = feature_edge2long[eid3[twosamedirid[1]]];
                    bool already_delete_flag = false;
                    auto it = deleted_long_id.find(le0);
                    if (it != deleted_long_id.end())
                        continue;
                    it = deleted_long_id.find(le1);
                    if (it != deleted_long_id.end())
                        continue;
                    double dist0 = long_edge_distortion[feature_edge2long[eid3[twosamedirid[0]]]];
                    double dist1 = long_edge_distortion[feature_edge2long[eid3[twosamedirid[1]]]];
                    if (dist0 > dist1)
                    {
                        deleted_long_id.insert(feature_edge2long[eid3[twosamedirid[0]]]);
                    }
                    else
                    {
                        deleted_long_id.insert(feature_edge2long[eid3[twosamedirid[1]]]);
                    }
                }
                else
                {
                    std::vector<int> twosamediraxis;
                    for (size_t j = 0; j < 3; j++)
                    {
                        for (size_t k = j + 1; k < 3; k++)
                        {
                            if (axis[j] == axis[k])
                            {
                                twosamediraxis.push_back(j);
                                twosamediraxis.push_back(k);
                                break;
                            }
                        }
                    }
                    assert(twosamediraxis.size() == 2);
                    int del_id = 3 - twosamediraxis[0] - twosamediraxis[1];
                    deleted_long_id.insert(feature_edge2long[eid3[del_id]]);
                }
            }
            break;
        }
        case 2:
        {
            int eid0 = feature_v2e[i][0];
            int eid1 = feature_v2e[i][1];
            OpenVolumeMesh::Geometry::Vec3d tn0 = target_fea_n[eid0];
            OpenVolumeMesh::Geometry::Vec3d tn1 = target_fea_n[eid1];
            int first0 = feature_e2v[eid0].first;
            int first1 = feature_e2v[eid1].first;
            if (first0 != i) tn0 = -tn0;
            if (first1 != i) tn1 = -tn1;
            int n0 = normal2idx[tn0];
            int n1 = normal2idx[tn1];
            int axis0 = n0 / 2;
            int axis1 = n1 / 2;
            if (axis0 == axis1)
            {
                if (abs(n0 - n1) == 0)
                {
                    feature_OK = false;
                    deleted_long_id.insert(feature_edge2long[eid0]);
                }
            }
            break;
        }
        case 1:
        {
            int eid = feature_v2e[i][0];
            int other_id = feature_e2v[eid].first;
            if (other_id == i) other_id = feature_e2v[eid].second;
            if (feature_v2e[other_id].size() == 1)
            {
                feature_OK = false;
                deleted_long_id.insert(feature_edge2long[eid]);
            }
            break;
        }
        default:
            break;
        }
    }
    assert(!target_bfn.empty());
    double delete_th_ratio = 0.6;
    for (size_t i = 0; i < long_feature_edge.size(); i++)
    {
        double n_wrong_edge = 0.0;
        for (size_t j = 0; j < long_feature_edge[i].size(); j++)
        {
            int eid = long_feature_edge[i][j];
            int nc0 = bef_id[eid][0];
            int nc1 = bef_id[eid][1];
            OpenVolumeMesh::Geometry::Vec3d tn0 = target_bfn[nc0], tn1 = target_bfn[nc1];
            if (tn0 == tn1)
            {
                n_wrong_edge = n_wrong_edge + 1.0;
            }
        }
        if (n_wrong_edge / long_feature_edge[i].size() > delete_th_ratio)
        {
            feature_OK = false;
            deleted_long_id.insert(i);
        }
    }
    if (repair_status)
    {
        if (!deleted_long_id.empty())
        {
            for (auto i : deleted_long_id)
            {
                for (size_t j = 0; j < long_feature_edge[i].size(); ++j)
                {
                    feature_edge_flag[long_feature_edge[i][j]] = false;
                }
            }
            feature_edge_change_flag = true;
            feature_edge_array.clear();
        }
    }
    if (repair_status)
    {
    }
    else
    {
        std::cout << "Check Feature Edge Correct: " << feature_OK << std::endl;
    }
    return feature_OK;
}
void polycube_deformation_interface::filter_boundary_normal_RGNF(double sigma_s_, double sigma_r_, int iter_count_, bool xyz_flag)
{
    face_center.resize(bfv_id.size());
    std::vector<double> six_pos(6 * boundary_face_number, 0.0);
    std::vector<double> N0(3 * boundary_face_number);
    int count = 0; OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
    for (int i = 0; i < bfv_id.size(); ++i)
    {
        OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
        if (one_bfv_id[0] < 0) continue;
        p0[0] = dpx[one_bfv_id[0]]; p0[1] = dpy[one_bfv_id[0]]; p0[2] = dpz[one_bfv_id[0]];
        p1[0] = dpx[one_bfv_id[1]]; p1[1] = dpy[one_bfv_id[1]]; p1[2] = dpz[one_bfv_id[1]];
        p2[0] = dpx[one_bfv_id[2]]; p2[1] = dpy[one_bfv_id[2]]; p2[2] = dpz[one_bfv_id[2]];
        OpenVolumeMesh::Geometry::Vec3d pt = (p0 + p1 + p2) / 3.0;
        OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0)*0.5;
        face_center[i] = pt;
        if (xyz_flag)
        {
            double fa = n.norm()*0.5;
            double max_component = std::abs(n[0]); int max_id = 0;
            OpenVolumeMesh::Geometry::Vec3d axis(0.0, 0.0, 0.0);
            for (int j = 1; j < 3; j++)
            {
                if (std::abs(n[j]) > max_component)
                {
                    max_component = std::abs(n[j]);
                    max_id = j;
                }
            }
            axis[max_id] = n[max_id] / max_component;
            n = axis * fa;
        }
        six_pos[6 * count + 0] = pt[0]; six_pos[6 * count + 1] = pt[1]; six_pos[6 * count + 2] = pt[2];
        N0[3 * count + 0] = n[0]; N0[3 * count + 1] = n[1]; N0[3 * count + 2] = n[2];
        ++count;
    }
    double inv_pos = 1.0 / (sigma_s_ *avg_boundary_edge_length);
    double inv_nor = 1.0 / sigma_r_;
    for (int i = 0; i < count; ++i)
    {
        six_pos[6 * i + 0] *= inv_pos; six_pos[6 * i + 1] *= inv_pos; six_pos[6 * i + 2] *= inv_pos;
    }
    std::vector<double> N1(3 * boundary_face_number, 0.0);
    for (int i = 0; i < iter_count_; i++)
    {
        PermutohedralLattice::filter(six_pos.data(), 6, N0.data(), 3, boundary_face_number, N1.data());
        for (int j = 0; j < boundary_face_number; j++)
        {
            OpenVolumeMesh::Geometry::Vec3d tn(N1[3 * j + 0], N1[3 * j + 1], N1[3 * j + 2]);
            double d = tn.norm();
            if (d < 1.0e-16)
            {
                OpenVolumeMesh::Geometry::Vec3d tn2(N0[3 * j + 0], N0[3 * j + 1], N0[3 * j + 2]);
                tn2.normalize();
                N1[3 * j + 0] = tn2[0]; N1[3 * j + 1] = tn2[1]; N1[3 * j + 2] = tn2[2];
            }
            else
            {
                N1[3 * j + 0] = tn[0] / d; N1[3 * j + 1] = tn[1] / d; N1[3 * j + 2] = tn[2] / d;
            }
            six_pos[6 * j + 3] = N1[3 * j + 0] * inv_nor;
            six_pos[6 * j + 4] = N1[3 * j + 1] * inv_nor;
            six_pos[6 * j + 5] = N1[3 * j + 2] * inv_nor;
        }
    }
    count = 0;
    if (target_bfn.size() != bfv_id.size())
    {
        target_bfn.resize(bfv_id.size());
    }
    for (int i = 0; i < bfv_id.size(); ++i)
    {
        if (bfv_id[i][0] < 0) continue;
        {
            double max_component = std::abs(N1[count * 3 + 0]); int max_id = 0;
            OpenVolumeMesh::Geometry::Vec3d axis(0.0, 0.0, 0.0);
            for (int j = 1; j < 3; j++)
            {
                if (std::abs(N1[count * 3 + j]) > max_component)
                {
                    max_component = std::abs(N1[count * 3 + j]);
                    max_id = j;
                }
            }
            axis[max_id] = N1[count * 3 + max_id] / max_component;
            target_bfn[i] = axis;
        }
        ++count;
    }
}
void polycube_deformation_interface::filter_boundary_normal_ring(double sigma_s_, int iter_count_, bool xyz_flag )
{
    avg_boundary_edge_length = 0.0;
    std::vector<OpenVolumeMesh::Geometry::Vec3d> N0(boundary_face_number);
    std::vector<OpenVolumeMesh::Geometry::Vec3d> fc(boundary_face_number);
    int count = 0; int nc = bfv_id.size(); OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
    std::vector<int> map_face(nc, -1);
    for (int i = 0; i < nc; ++i)
    {
        OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
        if (one_bfv_id[0] < 0) continue;
        p0[0] = dpx[one_bfv_id[0]]; p0[1] = dpy[one_bfv_id[0]]; p0[2] = dpz[one_bfv_id[0]];
        p1[0] = dpx[one_bfv_id[1]]; p1[1] = dpy[one_bfv_id[1]]; p1[2] = dpz[one_bfv_id[1]];
        p2[0] = dpx[one_bfv_id[2]]; p2[1] = dpy[one_bfv_id[2]]; p2[2] = dpz[one_bfv_id[2]];
        avg_boundary_edge_length += (p0 - p1).norm();
        avg_boundary_edge_length += (p1 - p2).norm();
        avg_boundary_edge_length += (p2 - p0).norm();
        OpenVolumeMesh::Geometry::Vec3d pt = (p0 + p1 + p2) / 3.0;
        OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0)*0.5;
        if (xyz_flag)
        {
            double fa = n.norm()*0.5;
            double max_component = std::abs(n[0]); int max_id = 0;
            OpenVolumeMesh::Geometry::Vec3d axis(0.0, 0.0, 0.0);
            for (int j = 1; j < 3; j++)
            {
                if (std::abs(n[j]) > max_component)
                {
                    max_component = std::abs(n[j]);
                    max_id = j;
                }
            }
            axis[max_id] = n[max_id] / max_component;
            n = axis * fa;
        }
        fc[count] = pt; N0[count] = n;
        map_face[i] = count;
        ++count;
    }
    avg_boundary_edge_length /= (3.0*count);
    double inv_pos = 1.0 / (sigma_s_ *avg_boundary_edge_length);
    double inv_pos2 = inv_pos * inv_pos*0.5;
    std::vector<int> visited_f(nc, -1); std::vector<int> ring_face; ring_face.reserve(25);
    std::vector<OpenVolumeMesh::Geometry::Vec3d> N1(boundary_face_number);
    for (int i = 0; i < iter_count_; i++)
    {
        for (int f_id = 0; f_id < nc; ++f_id)
        {
            OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[f_id];
            if (one_bfv_id[0] < 0) continue;
#if 1
            for (int j = 0; j < 3; ++j)
            {
                std::vector<int>& one_bvf_id = bvf_id[one_bfv_id[j]];
                for (int k = 0; k < one_bvf_id.size(); ++k)
                {
                    if (visited_f[one_bvf_id[k]] < 0)
                    {
                        visited_f[one_bvf_id[k]] = 1;
                        ring_face.push_back(one_bvf_id[k]);
                    }
                }
            }
#else
            for (int j = 0; j < 3; ++j)
            {
                std::vector<int>& one_bvf_id = bvf_id[one_bfv_id[j]];
                for (int k = 0; k < one_bvf_id.size(); ++k)
                {
                    if (visited_f[one_bvf_id[k]] < 0)
                    {
                        visited_f[one_bvf_id[k]] = 1;
                        ring_face.push_back(one_bvf_id[k]);
                    }
                }
            }
            for (int ijk = 0; ijk < ring_face.size(); ++ijk)
            {
                OpenVolumeMesh::Geometry::Vec3i& one_bfv_id_ = bfv_id[ring_face[ijk]];
                for (int j = 0; j < 3; ++j)
                {
                    std::vector<int>& one_bvf_id = bvf_id[one_bfv_id_[j]];
                    for (int k = 0; k < one_bvf_id.size(); ++k)
                    {
                        if (visited_f[one_bvf_id[k]] < 0)
                        {
                            visited_f[one_bvf_id[k]] = 1;
                            ring_face.push_back(one_bvf_id[k]);
                        }
                    }
                }
            }
#endif
            OpenVolumeMesh::Geometry::Vec3d& fc_0 = fc[map_face[f_id]];
            OpenVolumeMesh::Geometry::Vec3d n(0, 0, 0);
            for (int j = 0; j < ring_face.size(); ++j)
            {
                double len2 = (fc[map_face[ring_face[j]]] - fc_0).sqrnorm();
                double w = std::exp(-len2 * inv_pos2);
                n += w * N0[map_face[ring_face[j]]];
            }
            double d = n.norm();
            if (d < 1.0e-16)
            {
                N1[map_face[f_id]] = N0[map_face[f_id]];
                N1[map_face[f_id]].normalize();
            }
            else
            {
                N1[map_face[f_id]] = n / d;
            }
            for (int j = 0; j < ring_face.size(); ++j)
            {
                visited_f[ring_face[j]] = -1;
            }
            ring_face.clear(); ring_face.reserve(25);
        }
    }
    count = 0;
    if (target_bfn.size() != bfv_id.size())
    {
        target_bfn.resize(bfv_id.size());
    }
    for (int i = 0; i < bfv_id.size(); ++i)
    {
        if (bfv_id[i][0] < 0) continue;
        {
            target_bfn[i] = N1[count];
        }
        ++count;
    }
}
void polycube_deformation_interface::filter_feature_edge_ring(int n_level, double ss, bool consider_face_normal)
{
    if (feature_edge_array.empty()) set_feature_edge_array();
    if (feature_edge_array.empty()) return;
    target_fea_n.clear();
    target_fea_n.resize(feature_edge_flag.size(), OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, 0.0));
    std::vector<OpenVolumeMesh::Geometry::Vec3d> fe_center, fe_normal;
    fe_center.resize(bef_id.size());
    fe_normal.resize(bef_id.size());
    double avg_feature_edge_length = 0.0;
    OpenVolumeMesh::Geometry::Vec3d sixnormal[6] = {
        OpenVolumeMesh::Geometry::Vec3d(1.0, 0.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(-1.0, 0.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, 1.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, -1.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, 1.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, -1.0),
    };
    std::map<OpenVolumeMesh::Geometry::Vec3d, int> normal2idx;
    normal2idx[OpenVolumeMesh::Geometry::Vec3d(1.0, 0.0, 0.0)] = 0;
    normal2idx[OpenVolumeMesh::Geometry::Vec3d(-1.0, 0.0, 0.0)] = 1;
    normal2idx[OpenVolumeMesh::Geometry::Vec3d(0.0, 1.0, 0.0)] = 2;
    normal2idx[OpenVolumeMesh::Geometry::Vec3d(0.0, -1.0, 0.0)] = 3;
    normal2idx[OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, 1.0)] = 4;
    normal2idx[OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, -1.0)] = 5;
    for (size_t i = 0; i < feature_edge_array.size(); i++)
    {
        int eid = feature_edge_array[i];
        int v0 = feature_e2v[eid].first, v1 = feature_e2v[eid].second;
        OpenVolumeMesh::Geometry::Vec3d p0(0.0, 0.0, 0.0), p1(0.0, 0.0, 0.0);
        p0[0] = dpx[v0];
        p0[1] = dpy[v0];
        p0[2] = dpz[v0];
        p1[0] = dpx[v1];
        p1[1] = dpy[v1];
        p1[2] = dpz[v1];
        fe_center[eid] = (p0 + p1) / 2.0;
        fe_normal[eid] = (p1 - p0).normalize();
        avg_feature_edge_length += (p1 - p0).norm();
    }
    avg_feature_edge_length = avg_feature_edge_length / (1.0 * feature_edge_array.size());
    double inv_pos = 1.0 / (ss *avg_feature_edge_length);
    double inv_pos2 = inv_pos * inv_pos* 0.5;
    std::vector<OpenVolumeMesh::Geometry::Vec3d> weighted_normal;
    weighted_normal.resize(bef_id.size());
    for (size_t i = 0; i < feature_edge_array.size(); i++)
    {
        int eid = feature_edge_array[i];
        std::vector<int> edge_color(bef_id.size(), -1);
        weighted_normal[eid] = fe_normal[eid];
        edge_color[eid] = 1;
        OpenVolumeMesh::Geometry::Vec3d c = fe_center[eid];
        std::vector<int> nb;
        nb.push_back(eid);
        for (size_t j = 0; j < n_level; j++)
        {
            std::vector<int> nb_backup = nb;
            for (size_t k = 0; k < nb_backup.size(); k++)
            {
                int tmp_eid = nb_backup[k];
                for (size_t p = 0; p < feature_neighbor[tmp_eid].size(); p++)
                {
                    int tmp_nb = feature_neighbor[tmp_eid][p];
                    if (edge_color[tmp_nb] != 1)
                    {
                        edge_color[tmp_nb] = 1;
                        nb.push_back(tmp_nb);
                    }
                }
            }
        }
        assert(nb.size() <= 2 * n_level + 1);
        for (size_t j = 0; j < nb.size(); j++)
        {
            if (nb[j] != eid)
            {
                double len2 = (fe_center[nb[j]] - c).sqrnorm();
                double w = std::exp(-len2 * inv_pos2);
                weighted_normal[eid] += w * fe_normal[nb[j]];
            }
        }
        weighted_normal[eid] = weighted_normal[eid].normalize();
        int max_id = -1;
        double max_dot = -1.0;
        if (consider_face_normal)
        {
            int nc0 = bef_id[eid][0];
            int nc1 = bef_id[eid][1];
            auto it0 = normal2idx.find(target_bfn[nc0]);
            auto it1 = normal2idx.find(target_bfn[nc1]);
            assert(it0 != normal2idx.end() && it1 != normal2idx.end());
            int nid0 = it0->second, nid1 = it1->second;
            int nid0_pair = 2 * (nid0 / 2) + 1 - nid0 % 2;
            int nid1_pair = 2 * (nid1 / 2) + 1 - nid1 % 2;
            std::vector<int> normal_color(6, -1);
            normal_color[nid0] = 1;
            normal_color[nid1] = 1;
            normal_color[nid0_pair] = 1;
            normal_color[nid1_pair] = 1;
            for (int j = 0; j < 6; j++)
            {
                if (normal_color[j] != 1)
                {
                    double tmp_dot = OpenVolumeMesh::Geometry::dot(weighted_normal[eid], sixnormal[j]);
                    if (tmp_dot > max_dot)
                    {
                        max_dot = tmp_dot;
                        max_id = j;
                    }
                }
            }
        }
        else
        {
            for (int j = 0; j < 6; j++)
            {
                {
                    double tmp_dot = OpenVolumeMesh::Geometry::dot(weighted_normal[eid], sixnormal[j]);
                    if (tmp_dot > max_dot)
                    {
                        max_dot = tmp_dot;
                        max_id = j;
                    }
                }
            }
        }
        assert(max_id != -1);
        target_fea_n[eid] = sixnormal[max_id];
    }
}
void polycube_deformation_interface::get_polyline_normal()
{
    OpenVolumeMesh::Geometry::Vec3d sixnormal[6] = {
        OpenVolumeMesh::Geometry::Vec3d(1.0, 0.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(-1.0, 0.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, 1.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, -1.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, 1.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, -1.0),
    };
    polyline_target_n.resize(feature_edge_flag.size(), OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, 0.0));
    for (size_t i = 0; i < feature_edge_flag.size(); i++)
    {
        if (feature_edge_flag[i] == false) continue;
        int eid = i;
        int v0 = feature_e2v[eid].first, v1 = feature_e2v[eid].second;
        OpenVolumeMesh::Geometry::Vec3d p0(0.0, 0.0, 0.0), p1(0.0, 0.0, 0.0);
        p0[0] = dpx_feature[v0];
        p0[1] = dpy_feature[v0];
        p0[2] = dpz_feature[v0];
        p1[0] = dpx_feature[v1];
        p1[1] = dpy_feature[v1];
        p1[2] = dpz_feature[v1];
        OpenVolumeMesh::Geometry::Vec3d tmpnormal = (p1 - p0).normalize();
        int max_id = -1;
        double max_dot = -1.0;
        for (int j = 0; j < 6; j++)
        {
            double tmp_dot = OpenVolumeMesh::Geometry::dot(tmpnormal, sixnormal[j]);
            if (tmp_dot > max_dot)
            {
                max_dot = tmp_dot;
                max_id = j;
            }
        }
        assert(max_id != -1);
        polyline_target_n[eid] = sixnormal[max_id];
    }
}
void polycube_deformation_interface::set_feature_edge_array()
{
    feature_edge_array.clear();
    for (int i = 0; i < feature_edge_flag.size(); i++)
    {
        if (feature_edge_flag[i] == true)
        {
            feature_edge_array.push_back(i);
        }
    }
}
void polycube_deformation_interface::set_feature_edge_length()
{
    feature_edge_length.clear();
    for (int i = 0; i < feature_edge_flag.size(); i++)
    {
        if (feature_edge_flag[i] == true)
        {
            int v0 = feature_e2v[i].first;
            int v1 = feature_e2v[i].second;
            OpenVolumeMesh::Geometry::Vec3d p0(dpx[v0], dpy[v0], dpz[v0]);
            OpenVolumeMesh::Geometry::Vec3d p1(dpx[v1], dpy[v1], dpz[v1]);
            feature_edge_length.push_back((p0 - p1).length());
        }
    }
    int count = 0;
    feature_g2l.resize(bvf_id.size(), -1);
    for (size_t i = 0; i < feature_v2v.size(); i++)
    {
        if (!feature_v2v[i].empty())
        {
            feature_g2l[i] = count++;
            feature_l2g.push_back(i);
        }
    }
}
void polycube_deformation_interface::filter_boundary_normal_feature(double ss, double sr, int level)
{
    bool feature_settled = true;
    if (feature_edge_flag.empty()) feature_settled = false;
    if (target_bfn.size() != bfv_id.size())
    {
        target_bfn.resize(bfv_id.size());
    }
    if (bfe_id.empty())
    {
        bfe_id.resize(bfv_id.size());
        for (size_t i = 0; i < bef_id.size(); i++)
        {
            for (size_t j = 0; j < bef_id[i].size(); j++)
            {
                int fid = bef_id[i][j];
                bfe_id[fid].push_back(i);
            }
        }
        for (size_t i = 0; i < bfe_id.size(); i++)
        {
            assert(bfe_id[i].size() == 3 || bfe_id[i].size() == 0);
        }
    }
    OpenVolumeMesh::Geometry::Vec3d sixnormal[6] = {
        OpenVolumeMesh::Geometry::Vec3d(1.0, 0.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(-1.0, 0.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, 1.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, -1.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, 1.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, -1.0),
    };
    avg_boundary_edge_length = 0.0;
    std::vector<OpenVolumeMesh::Geometry::Vec3d> N0(boundary_face_number);
    std::vector<OpenVolumeMesh::Geometry::Vec3d> fc(boundary_face_number);
    int count = 0; int nc = bfv_id.size(); OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
    std::vector<int> facea2b(nc, -1);
    face_center.resize(nc);
    for (int i = 0; i < nc; ++i)
    {
        OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
        if (one_bfv_id[0] < 0) continue;
        p0[0] = dpx[one_bfv_id[0]]; p0[1] = dpy[one_bfv_id[0]]; p0[2] = dpz[one_bfv_id[0]];
        p1[0] = dpx[one_bfv_id[1]]; p1[1] = dpy[one_bfv_id[1]]; p1[2] = dpz[one_bfv_id[1]];
        p2[0] = dpx[one_bfv_id[2]]; p2[1] = dpy[one_bfv_id[2]]; p2[2] = dpz[one_bfv_id[2]];
        avg_boundary_edge_length += (p0 - p1).norm();
        avg_boundary_edge_length += (p1 - p2).norm();
        avg_boundary_edge_length += (p2 - p0).norm();
        OpenVolumeMesh::Geometry::Vec3d pt = (p0 + p1 + p2) / 3.0;
        OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0)*0.5;
        fc[count] = pt; N0[count] = n;
        face_center[i] = pt;
        facea2b[i] = count;
        ++count;
    }
    avg_boundary_edge_length /= (3.0*count);
    double inv_pos = 1.0 / (ss *avg_boundary_edge_length);
    double inv_pos2 = inv_pos * inv_pos*0.5;
    double inv_nor2 = 1.0 / (2.0 * sr * sr);
    double range = 1.0 * level * avg_boundary_edge_length;
    int avg_nb_size = 0;
    for (size_t i = 0; i < nc; i++)
    {
        OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
        if (one_bfv_id[0] < 0) continue;
        OpenVolumeMesh::Geometry::Vec3d center = fc[facea2b[i]];
        OpenVolumeMesh::Geometry::Vec3d cnormal = N0[facea2b[i]];
        std::vector<int> level_flag(nc, -1);
        std::vector<std::vector<int>> level_cells(2 * level);
        level_cells[0].push_back(i);
        level_flag[i] = 1;
        for (size_t j = 1; j < 2 * level; j++)
        {
            for (size_t i1 = 0; i1 < level_cells[j - 1].size(); i1++)
            {
                int tmpc = level_cells[j - 1][i1];
                for (size_t i2 = 0; i2 < 3; i2++)
                {
                    int tmpv = bfv_id[tmpc][i2];
                    for (size_t i3 = 0; i3 < bvf_id[tmpv].size(); i3++)
                    {
                        int tmpcc = bvf_id[tmpv][i3];
                        if (level_flag[tmpcc] == -1)
                        {
                            level_flag[tmpcc] = 1;
                            level_cells[j].push_back(tmpcc);
                        }
                    }
                }
            }
        }
        std::vector<int> nb;
        std::vector<int> cell_color(nc, -1);
        nb.push_back(i);
        cell_color[i] = 1;
        std::queue<int> q;
        q.push(i);
        while (!q.empty())
        {
            int front = q.front(); q.pop();
            for (size_t j = 0; j < bfe_id[front].size(); j++)
            {
                int eid = bfe_id[front][j];
                int otherid = bef_id[eid][0];
                if (otherid == i) otherid = bef_id[eid][1];
                if (feature_settled && feature_edge_flag[eid] == 1)
                    continue;
                if (cell_color[otherid] != 1 && level_flag[otherid] == 1)
                {
                    if ((fc[facea2b[otherid]] - center).length() < range)
                    {
                        cell_color[otherid] = 1;
                        q.push(otherid);
                        nb.push_back(otherid);
                    }
                }
            }
        }
        avg_nb_size += nb.size();
        OpenVolumeMesh::Geometry::Vec3d wn(0.0, 0.0, 0.0);
        for (size_t j = 0; j < nb.size(); j++)
        {
            int tmp = nb[j];
            double len2 = (fc[facea2b[tmp]] - center).sqrnorm();
            double normal2 = (N0[facea2b[tmp]] - cnormal).sqrnorm();
            double w = std::exp(-len2 * inv_pos2 - normal2 * inv_nor2);
            w = std::exp(-len2 * inv_pos2);
            wn += w * N0[facea2b[tmp]];
        }
        wn.normalize();
        int max_id = -1;
        double max_dot = -1.0;
        for (size_t j = 0; j < 6; j++)
        {
            double tmp_dot = OpenVolumeMesh::Geometry::dot(wn, sixnormal[j]);
            if (tmp_dot > max_dot)
            {
                max_dot = tmp_dot;
                max_id = j;
            }
        }
        assert(max_id != -1);
        target_bfn[i] = sixnormal[max_id];
    }
}
bool polycube_deformation_interface::target_normal_check_repair(bool repair_flag)
{
    if (bfe_id.empty())
    {
        bfe_id.resize(bfv_id.size());
        for (size_t i = 0; i < bef_id.size(); i++)
        {
            for (size_t j = 0; j < bef_id[i].size(); j++)
            {
                int fid = bef_id[i][j];
                bfe_id[fid].push_back(i);
            }
        }
        for (size_t i = 0; i < bfe_id.size(); i++)
        {
            assert(bfe_id[i].size() == 3 || bfe_id[i].size() == 0);
        }
    }
    error_edges.clear();
    for (size_t i = 0; i < bef_id.size(); i++)
    {
        std::vector<OpenVolumeMesh::Geometry::Vec3d> twonormal;
        for (size_t j = 0; j < bef_id[i].size(); j++)
        {
            int fid = bef_id[i][j];
            twonormal.push_back(target_bfn[fid]);
        }
        if (twonormal.size() == 2)
        {
            double tmpdot = OpenVolumeMesh::Geometry::dot(twonormal[0], twonormal[1]);
            if (tmpdot < -0.5)
            {
                std::cout << "!!!!!Opposite Target Direction along One Edge!!!!! " << i << std::endl;
                error_edges.push_back(i);
            }
        }
    }
    return true;
}
void polycube_deformation_interface::exp_mips_deformation_refine_polycube_omp_feature_nf(int max_iter, int iter2, double energy_power, double angle_area_ratio, VolumeMesh *mesh_, bool use_xyz, bool use_RGNF, bool use_half_sigma_s, bool feature_flag, bool boundary_global_opt)
{
    double ga = 1; double lamda = 0.1;
    if (feature_edge_change_flag) build_feature_edges_connectivity(mesh_);
    if (feature_edge_length.empty()) set_feature_edge_length();
    if (!polyline_target_n.empty())
    {
        target_fea_n = polyline_target_n;
    }
    for (int j = 0; j < max_iter; ++j)
    {
        if (feature_flag)
        {
            if (lambda_array.empty())
            {
                lambda_array.push_back(1.0);
                lambda_array.push_back(1000);
                lambda_array.push_back(1.0);
            }
            else
            {
                lambda_array[0] = 10 * lambda_array[0];
                lambda_array[2] = lambda_array[2] / 10;
                if (lambda_array[0] > 1e16) lambda_array[0] = 1e16;
            }
            for (int k = 0; k < iter2; ++k)
            {
                optimize_all_verts_nf();
            }
            double sr = 0.3;
        }
    }
}
bool polycube_deformation_interface::check_polyline_orthogonality(double threshold)
{
    double ave_f(0), max_f(0);
    double poorcount(0.0);
    const static OpenVolumeMesh::Geometry::Vec3d sixnormal[6] = {
        OpenVolumeMesh::Geometry::Vec3d(1.0, 0.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(-1.0, 0.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, 1.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, -1.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, 1.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, -1.0),
    };
    int n_fea = (int)feature_edge_array.size();
    for (size_t i = 0; i < n_fea; i++)
    {
        int eid = feature_edge_array[i];
        int aid0 = feature_e2v[eid].first, aid1 = feature_e2v[eid].second;
        double x1, x2, y1, y2, z1, z2;
        x1 = dpx_feature[aid0];
        y1 = dpy_feature[aid0];
        z1 = dpz_feature[aid0];
        x2 = dpx_feature[aid1];
        y2 = dpy_feature[aid1];
        z2 = dpz_feature[aid1];
        OpenVolumeMesh::Geometry::Vec3d v1(x1, y1, z1);
        OpenVolumeMesh::Geometry::Vec3d v2(x2, y2, z2);
        OpenVolumeMesh::Geometry::Vec3d tmpn = v2 - v1;
        tmpn.normalize();
        double dot = 0;
        for (size_t j = 0; j < 6; j++)
        {
            dot = std::max(dot, OpenVolumeMesh::Geometry::dot(tmpn, sixnormal[j]));
        }
        if (dot > 1) dot = 1;
        dot = std::acos(dot) / 3.1415926 * 180;
        ave_f += dot;
        if (max_f < dot)
            max_f = dot;
        if (dot > 5)
            poorcount = poorcount + 1.0;
    }
    ave_f = ave_f / n_fea;
    double ratio = poorcount / n_fea;
    std::cout << "polyline difference angle > 5бу: " << ratio << std::endl;
    std::cout << "average: " << ave_f << std::endl;
    return threshold > ratio;
}
bool polycube_deformation_interface::check_polycube_planarity(double threshold)
{
    double ave_f(0), max_f(0);
    double facecount(0.0);
    double maxcount(0.0);
    const static OpenVolumeMesh::Geometry::Vec3d sixnormal[6] = {
        OpenVolumeMesh::Geometry::Vec3d(1.0, 0.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(-1.0, 0.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, 1.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, -1.0, 0.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, 1.0),
        OpenVolumeMesh::Geometry::Vec3d(0.0, 0.0, -1.0),
    };
    for (size_t i = 0; i < bfv_id.size(); i++)
    {
        if (bfv_id[i][0] == -1) continue;
        facecount = facecount + 1.0;
        int v_id_0 = bfv_id[i][0];
        int v_id_1 = bfv_id[i][1];
        int v_id_2 = bfv_id[i][2];
        double x0, y0, z0;
        double x1, x2, y1, y2, z1, z2;
        x0 = dpx[v_id_0];
        y0 = dpy[v_id_0];
        z0 = dpz[v_id_0];
        x1 = dpx[v_id_1];
        y1 = dpy[v_id_1];
        z1 = dpz[v_id_1];
        x2 = dpx[v_id_2];
        y2 = dpy[v_id_2];
        z2 = dpz[v_id_2];
        OpenVolumeMesh::Geometry::Vec3d v0(x0, y0, z0);
        OpenVolumeMesh::Geometry::Vec3d v1(x1, y1, z1);
        OpenVolumeMesh::Geometry::Vec3d v2(x2, y2, z2);
        OpenVolumeMesh::Geometry::Vec3d tmpn = OpenVolumeMesh::Geometry::cross(v1 - v0, v2 - v0);
        tmpn.normalize();
        double dot = 0;
        for (size_t j = 0; j < 6; j++)
        {
            dot = std::max(dot, OpenVolumeMesh::Geometry::dot(tmpn, sixnormal[j]));
        }
        if (dot > 1) dot = 1;
        dot = std::acos(dot) / 3.1415926 * 180;
        ave_f += dot;
        if (max_f < dot)
            max_f = dot;
        if (dot > 5)
            maxcount = maxcount + 1.0;
    }
    ave_f = ave_f / facecount;
    double ratio = maxcount / facecount;
    std::cout << "difference angle > 5бу: " << ratio << std::endl;
    std::cout << "average: " << ave_f << std::endl;
    return threshold > ratio;
}
void polycube_deformation_interface::deformation_polylines(VolumeMesh* mesh_, bool assign_current, bool ns_flag, double smooth_factor)
{
    if (feature_edge_change_flag) build_feature_edges_connectivity(mesh_);
    if (feature_edge_length.empty()) set_feature_edge_length();
    polyline_ns_flag = ns_flag;
    {
        if (lambda_array_polyline.empty())
        {
            if (polyline_ns_flag)
            {
                lambda_array_polyline.push_back(2);
            }
            else
            {
                lambda_array_polyline.push_back(0.01);
            }
            lambda_array_polyline.push_back(smooth_factor);
        }
        else
        {
            lambda_array_polyline[0] = 10 * lambda_array_polyline[0];
            if (lambda_array_polyline[0] > 1e16) lambda_array_polyline[0] = 1e16;
            lambda_array_polyline[1] = lambda_array_polyline[1] / 10;
        }
        optimize_all_polylines(assign_current);
        filter_feature_edge_ring(2, sigma_s, false);
    }
}
void polycube_deformation_interface::exp_mips_deformation_refine_polycube_omp_constrained(int max_iter, int iter2, double energy_power, double angle_area_ratio, int nv, bool use_xyz, bool use_RGNF, bool use_half_sigma_s, bool change_mu)
{
    double ga = 1; double lamda = 0.1;
    double temp_ss = sigma_s;
    if (use_half_sigma_s) sigma_s *= 0.5;
    for (int a = 0; a < max_iter; ++a)
    {
        if (fixed_chart.size() == 0)
        {
            if (use_RGNF) filter_boundary_normal_RGNF(sigma_s, sigma_r, RGNF_iter_count, use_xyz);
            else filter_boundary_normal_ring(2 * sigma_s, RGNF_iter_count, true);
        }
        use_half_sigma_s = true;
        ga = compute_energy_ratio_filter_deform();
        if (change_mu)
        {
            int cut_face_tet_size = cut_face_tet_idx.size();
            for (size_t i = 0; i < cut_face_tet_size; i++)
            {
                local_energy_ratio[cut_face_tet_idx[i]] = 5 * local_energy_ratio[cut_face_tet_idx[i]];
            }
        }
        else
        {
            int cut_face_tet_size = cut_face_tet_idx.size();
            for (size_t i = 0; i < cut_face_tet_size; i++)
            {
                local_energy_ratio[cut_face_tet_idx[i]] = 1.2 * local_energy_ratio[cut_face_tet_idx[i]];
            }
        }
        int non_zero = 0;
        for (size_t i = 0; i < distortion_big_count.size(); i++)
        {
            if (distortion_big_count[i] != 0)
                non_zero++;
        }
        int max_num_t = omp_get_max_threads();
        distortion_big_count_omp.resize(max_num_t);
        for (size_t k = 0; k < max_num_t; k++)
        {
            distortion_big_count_omp[k].resize(distortion_big_count.size(), 0);
            std::fill(distortion_big_count_omp[k].begin(), distortion_big_count_omp[k].end(), 0);
        }
        vc_count_omp.resize(max_num_t);
        for (int k = 0; k < iter2; ++k)
        {
            int n_color = nv;
            if (number_of_color != 0) n_color = number_of_color;
            for (unsigned i = 0; i < n_color; ++i)
            {
                int one_color_size = vertex_diff_color[i].size();
#pragma omp parallel for schedule(dynamic, 1)
                for (int j = 0; j < one_color_size; ++j)
                {
                    int id = omp_get_thread_num();
                    int temp_vert = vertex_diff_color[i][j];
                    if (vert_update_type[temp_vert] == 3)
                    {
                        int common_vert_idx = three_cut_common_vertex_idx[vert_belong_three_cut_array[temp_vert]];
                        int mat_type = three_cut_types[vert_belong_three_cut_array[temp_vert]];
                        exp_mips_deformation_refine_one_polycube_normal_constrained_three_cut(temp_vert, angle_area_ratio, energy_power, id, three_cut_adjacent_one_cut_type[vert_belong_three_cut_array[temp_vert]], common_vert_idx, ga, use_half_sigma_s);
                    }
                    else if (vert_update_type[temp_vert] == 0)
                        exp_mips_deformation_refine_one_polycube_normal_constrained(temp_vert, angle_area_ratio, energy_power, id, ga, use_half_sigma_s);
                    else if (vert_update_type[temp_vert] == 1)
                    {
                        int common_vert_idx = cut_common_verts_idx[vert_belong_cut_array[temp_vert] << 1];
                        int mat_type = cut_types[vert_belong_cut_array[temp_vert]];
                        exp_mips_deformation_refine_one_polycube_normal_constrained_equal_face(temp_vert, angle_area_ratio, energy_power, id, type_matrix[mat_type], common_vert_idx, ga, use_half_sigma_s);
                    }
                    else
                    {
                        continue;
                    }
                }
                for (size_t j = 0; j < cut_common_verts_idx.size(); j++)
                {
                    int temp_vert = cut_common_verts_idx[j];
                    if (vert_update_type[temp_vert] == 4)
                    {
                        int mat_type = cut_types[vert_belong_cut_array[temp_vert]];
                        exp_mips_deformation_refine_one_polycube_normal_constrained_single_line(temp_vert, angle_area_ratio, energy_power, 0, type_four_matrix[mat_type], ga, use_half_sigma_s);
                    }
                }
                for (size_t a = 0; a < max_num_t; a++)
                {
                    for (size_t b = 0; b < distortion_big_count.size(); b++)
                    {
                        distortion_big_count[b] += distortion_big_count_omp[a][b];
                    }
                }
            }
        }
        if (bound_K < 4.0)
        {
            bound_K += 0.1;
        }
    }
    sigma_s = temp_ss;
}
void polycube_deformation_interface::optimize_all_verts_nf()
{
    std::cout << "__________________________________" << std::endl;
    std::cout << "begin optimize all verts" << std::endl;
    std::cout << "lambda: " << lambda_array[0] << std::endl;
    deform_pointer dp;
    dp.pdi = this;
    std::vector<double> solution(3 * bvf_id.size(), 0.0);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < bvf_id.size(); i++)
    {
        int vid = i;
        solution[3 * i] = dpx[vid];
        solution[3 * i + 1] = dpy[vid];
        solution[3 * i + 2] = dpz[vid];
    }
    de_LBFGS_all_nf(&dp, solution);
    assert(solution.size() == 3 * bvf_id.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < bvf_id.size(); i++)
    {
        int vid = i;
        dpx[vid] = solution[3 * i];
        dpy[vid] = solution[3 * i + 1];
        dpz[vid] = solution[3 * i + 2];
    }
    return;
}
void polycube_deformation_interface::optimize_all_polylines(bool assign_current)
{
    std::cout << "__________________________________" << std::endl;
    std::cout << "begin optimize all verts" << std::endl;
    deform_pointer dp;
    dp.pdi = this;
    std::vector<double> solution(3 * feature_l2g.size(), 0.0);
    if (!assign_current)
    {
        if (dpx_feature.empty())
        {
            dpx_feature = dpx;
            dpy_feature = dpy;
            dpz_feature = dpz;
        }
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < feature_l2g.size(); i++)
        {
            int lid = i;
            int gid = feature_l2g[lid];
            solution[3 * i] = dpx_feature[gid];
            solution[3 * i + 1] = dpy_feature[gid];
            solution[3 * i + 2] = dpz_feature[gid];
        }
    }
    else
    {
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < feature_l2g.size(); i++)
        {
            int lid = i;
            int gid = feature_l2g[lid];
            solution[3 * i] = dpx[gid];
            solution[3 * i + 1] = dpy[gid];
            solution[3 * i + 2] = dpz[gid];
        }
    }
    de_LBFGS_polylines(&dp, solution);
    assert(solution.size() == 3 * feature_l2g.size());
    if (assign_current)
    {
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < feature_l2g.size(); i++)
        {
            int lid = i;
            int gid = feature_l2g[lid];
            dpx[gid] = solution[3 * i];
            dpy[gid] = solution[3 * i + 1];
            dpz[gid] = solution[3 * i + 2];
        }
    }
    else
    {
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < feature_l2g.size(); i++)
        {
            int lid = i;
            int gid = feature_l2g[lid];
            dpx_feature[gid] = solution[3 * i];
            dpy_feature[gid] = solution[3 * i + 1];
            dpz_feature[gid] = solution[3 * i + 2];
        }
    }
    return;
}
void polycube_deformation_interface::compute_face_normal(std::vector< OpenVolumeMesh::Geometry::Vec3d >& bfn)
{
    int nc = bfv_id.size();
    OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
    for (int i = 0; i < nc; ++i)
    {
        OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[i];
        if (fv_id[0] >= 0)
        {
            p0[0] = dpx[fv_id[0]]; p0[1] = dpy[fv_id[0]]; p0[2] = dpz[fv_id[0]];
            p1[0] = dpx[fv_id[1]]; p1[1] = dpy[fv_id[1]]; p1[2] = dpz[fv_id[1]];
            p2[0] = dpx[fv_id[2]]; p2[1] = dpy[fv_id[2]]; p2[2] = dpz[fv_id[2]];
            bfn[i] = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0).normalize();
        }
    }
}
bool polycube_deformation_interface::check_normal_difference(double th)
{
    double cos_th = std::cos(th); int nc = bfv_id.size();
    OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
    if (change_big_flag.size() != nc) change_big_flag.resize(nc, -1);
    for (int i = 0; i < nc; ++i)
    {
        OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[i];
        change_big_flag[i] = -1;
        if (fv_id[0] >= 0)
        {
            OpenVolumeMesh::Geometry::Vec3d& ln = last_bfn[i];
            p0[0] = dpx[fv_id[0]]; p0[1] = dpy[fv_id[0]]; p0[2] = dpz[fv_id[0]];
            p1[0] = dpx[fv_id[1]]; p1[1] = dpy[fv_id[1]]; p1[2] = dpz[fv_id[1]];
            p2[0] = dpx[fv_id[2]]; p2[1] = dpy[fv_id[2]]; p2[2] = dpz[fv_id[2]];
            OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0).normalize();
            double cos_ = OpenVolumeMesh::Geometry::dot(n, ln);
            if (cos_ < cos_th) change_big_flag[i] = 1;
            ln = n;
        }
    }
    return true;
}
void polycube_deformation_interface::set_fixed_chart_label(std::vector<int> &chart, std::vector<int> &label, bool is_ovm)
{
    fixed_chart = chart;
    fixed_label = label;
    if (fixed_chart.size() == 0 || fixed_label.size() == 0)
        return;
    target_bfn.clear();
    target_bfn.resize(chart.size(), OpenVolumeMesh::Geometry::Vec3d(0, 0, 0));
    int nc = chart.size();
    int max_chart = -1;
    for (size_t i = 0; i < nc; i++)
    {
        if (chart[i] > max_chart)
            max_chart = chart[i];
    }
    chart_tet_idx.clear();
    chart_tet_idx.resize(max_chart + 1);
    for (size_t i = 0; i < nc; i++)
    {
        if (chart[i] >= 0)
            chart_tet_idx[chart[i]].push_back(i);
    }
    OpenVolumeMesh::Geometry::Vec3d n[6];
    for (size_t i = 0; i < 6; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            n[i][j] = 0;
        }
    }
    if (is_ovm)
    {
        n[0][0] = 1;
        n[1][0] = -1;
        n[2][1] = 1;
        n[3][1] = -1;
        n[4][2] = 1;
        n[5][2] = -1;
    }
    else
    {
        n[0][0] = -1;
        n[1][0] = 1;
        n[2][1] = -1;
        n[3][1] = 1;
        n[4][2] = -1;
        n[5][2] = 1;
    }
    for (size_t i = 0; i < nc; i++)
    {
        if (label[i] != -1)
            target_bfn[i] = n[label[i]];
    }
}
void polycube_deformation_interface::set_equal_faces(
    const std::vector<std::vector<std::vector<int>>> &faces_array,
    const std::vector<int>& common_verts_idx,
    const std::vector<int>& cut_types_array,
    std::vector<std::pair<int, int>> &chart_pair,
    std::vector<int> &three_cut_types_array,
    std::vector<int> &three_cut_common_vert,
    std::vector<std::vector<std::array<unsigned int, 3>>> &three_cut_vert,
    const std::vector<std::array<int, 3>> &temp_three_cut_adjacent_one_cut_type
)
{
    assert(faces_array.size() == (common_verts_idx.size() >> 1) && faces_array.size() == cut_types_array.size());
    assert(three_cut_types_array.size() == three_cut_common_vert.size() && three_cut_types_array.size() == three_cut_vert.size());
    if (faces_array.size() == 0)
        return;
    vert_update_type.clear();
    vert_belong_cut_array.clear();
    cut_types.clear();
    cut_common_verts_idx.clear();
    vert_pairs_map.clear();
    vert_update_type.resize(bvf_id.size(), 0);
    vert_belong_cut_array.resize(bvf_id.size(), -1);
    cut_types = cut_types_array;
    cut_common_verts_idx = common_verts_idx;
    equal_triangles.clear();
    for (size_t i = 0; i < faces_array.size(); i++)
    {
        equal_triangles.insert(equal_triangles.end(), faces_array[i].begin(), faces_array[i].end());
    }
    for (size_t i = 0; i < faces_array.size(); i++)
    {
        for (size_t j = 0; j < faces_array[i].size(); j++)
        {
            int id1, id2, id3, id4, id5, id0;
            id0 = faces_array[i][j][0];
            id1 = faces_array[i][j][1];
            id2 = faces_array[i][j][2];
            id3 = faces_array[i][j][3];
            id4 = faces_array[i][j][4];
            id5 = faces_array[i][j][5];
            vert_update_type[id0] = 1;
            vert_update_type[id1] = 1;
            vert_update_type[id2] = 1;
            vert_update_type[id3] = 2;
            vert_update_type[id4] = 2;
            vert_update_type[id5] = 2;
            vert_belong_cut_array[id0] = i;
            vert_belong_cut_array[id1] = i;
            vert_belong_cut_array[id2] = i;
            vert_belong_cut_array[id3] = i;
            vert_belong_cut_array[id4] = i;
            vert_belong_cut_array[id5] = i;
            vert_pairs_map[id0] = id3;
            vert_pairs_map[id1] = id4;
            vert_pairs_map[id2] = id5;
            vert_pairs_map[id3] = id0;
            vert_pairs_map[id4] = id1;
            vert_pairs_map[id5] = id2;
        }
    }
    for (size_t i = 0; i < cut_common_verts_idx.size(); i++)
    {
        int id = cut_common_verts_idx[i];
        vert_update_type[id] = 4;
    }
    int min_chart = chart_pair[0].first;
    for (size_t i = 0; i < chart_pair.size(); i++)
    {
        if (chart_pair[i].first < min_chart)
            min_chart = chart_pair[i].first;
        if (chart_pair[i].second < min_chart)
            min_chart = chart_pair[i].second;
    }
    cut_face_tet_idx.clear();
    for (size_t i = 0; i < fixed_chart.size(); i++)
    {
        if (fixed_chart[i] >= min_chart)
            cut_face_tet_idx.push_back(i);
    }
    vert_belong_three_cut_array.clear();
    three_cut_types.clear();
    three_cut_common_vertex_idx.clear();
    three_cut_vert_to_pair.clear();
    vert_belong_three_cut_array.resize(bvf_id.size(), 0);
    three_cut_types = three_cut_types_array;
    three_cut_common_vertex_idx = three_cut_common_vert;
    three_cut_adjacent_one_cut_type = temp_three_cut_adjacent_one_cut_type;
    for (size_t i = 0; i < three_cut_vert.size(); i++)
    {
        for (size_t j = 0; j < three_cut_vert[i].size(); j++)
        {
            int id0, id1, id2;
            id0 = three_cut_vert[i][j][0];
            id1 = three_cut_vert[i][j][1];
            id2 = three_cut_vert[i][j][2];
            vert_update_type[id0] = 3;
            vert_update_type[id1] = 2;
            vert_update_type[id2] = 2;
            vert_belong_three_cut_array[id0] = i;
            vert_belong_three_cut_array[id1] = i;
            vert_belong_three_cut_array[id2] = i;
            std::pair<int, int> tmp_pair(id1, id2);
            three_cut_vert_to_pair.insert(std::pair<int, std::pair<int, int>>(id0, tmp_pair));
        }
        vert_update_type[three_cut_common_vert[i]] = 2;
    }
}
void polycube_deformation_interface::adjust_orientation(VolumeMesh* mesh_, bool feature_flag)
{
    OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
    std::vector<OpenVolumeMesh::Geometry::Vec3d> all_bfn;
    for (int i = 0; i < bfv_id.size(); ++i)
    {
        OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
        if (one_bfv_id[0] < 0) continue;
        p0[0] = dpx[one_bfv_id[0]]; p0[1] = dpy[one_bfv_id[0]]; p0[2] = dpz[one_bfv_id[0]];
        p1[0] = dpx[one_bfv_id[1]]; p1[1] = dpy[one_bfv_id[1]]; p1[2] = dpz[one_bfv_id[1]];
        p2[0] = dpx[one_bfv_id[2]]; p2[1] = dpy[one_bfv_id[2]]; p2[2] = dpz[one_bfv_id[2]];
        OpenVolumeMesh::Geometry::Vec3d pt = (p0 + p1 + p2) / 3.0;
        OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0);
        all_bfn.push_back(n.normalize());
    }
    if (feature_flag)
    {
        for (size_t i = 0; i < feature_edge_flag.size(); i++)
        {
            if (feature_edge_flag[i] == true)
            {
                OpenVolumeMesh::EdgeHandle eh = OpenVolumeMesh::EdgeHandle(i);
                OpenVolumeMesh::OpenVolumeMeshEdge ovme = mesh_->edge(eh);
                OpenVolumeMesh::Geometry::Vec3d p0 = mesh_->vertex(ovme.from_vertex());
                OpenVolumeMesh::Geometry::Vec3d p1 = mesh_->vertex(ovme.to_vertex());
                OpenVolumeMesh::Geometry::Vec3d n = p1 - p0;
                all_bfn.push_back(n.normalize());
            }
        }
    }
    std::vector<double> X(3, 0.0);
    for (int i = 0; i < 1; ++i)
    {
        ad_LBFGS(all_bfn, X);
    }
    double alpha = X[0]; double beta = X[1]; double gamma = X[2];
    printf("%f %f %f\n", alpha, beta, gamma);
    double cos_phi = std::cos(X[0]), sin_phi = std::sin(X[0]);
    double cos_theta = std::cos(X[1]), sin_theta = std::sin(X[1]);
    double cos_xi = std::cos(X[2]), sin_xi = std::sin(X[2]);
    Matrix3d Mtheta, Mphi, Msi;
    Mtheta.setZero(); Mphi.setZero(); Msi.setZero();
    Mphi(0, 0) = 1; Mphi(1, 1) = Mphi(2, 2) = cos_phi, Mphi(1, 2) = sin_phi, Mphi(2, 1) = -sin_phi;
    Mtheta(1, 1) = 1; Mtheta(0, 0) = Mtheta(2, 2) = cos_theta, Mtheta(2, 0) = sin_theta, Mtheta(0, 2) = -sin_theta;
    Msi(2, 2) = 1; Msi(0, 0) = Msi(1, 1) = cos_xi, Msi(0, 1) = sin_xi, Msi(1, 0) = -sin_xi;
    Matrix3d M = Msi * Mphi * Mtheta;
    OpenVolumeMesh::Geometry::Vec3d p;
    for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
    {
        OpenVolumeMesh::Geometry::Vec3d q = mesh_->vertex(*v_it);
        Vector3d qvec(q[0], q[1], q[2]);
        Vector3d pvec = M * qvec;
        p[0] = pvec[0];
        p[1] = pvec[1];
        p[2] = pvec[2];
        mesh_->set_vertex(*v_it, p);
        int v_id = v_it->idx();
        dpx[v_id] = p[0] ; dpy[v_id] = p[1] ; dpz[v_id] = p[2] ;
    }
}
void polycube_deformation_interface::adjust_orientation_feature(VolumeMesh* mesh_, bool polyline_flag)
{
    OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
    std::vector<OpenVolumeMesh::Geometry::Vec3d> all_bfn;
    if (!polyline_flag)
    {
        assert(feature_edge_flag.size() != 0);
        for (size_t i = 0; i < feature_edge_flag.size(); i++)
        {
            if (feature_edge_flag[i] == true)
            {
                OpenVolumeMesh::EdgeHandle eh = OpenVolumeMesh::EdgeHandle(i);
                OpenVolumeMesh::OpenVolumeMeshEdge ovme = mesh_->edge(eh);
                OpenVolumeMesh::Geometry::Vec3d p0 = mesh_->vertex(ovme.from_vertex());
                OpenVolumeMesh::Geometry::Vec3d p1 = mesh_->vertex(ovme.to_vertex());
                OpenVolumeMesh::Geometry::Vec3d n = p1 - p0;
                all_bfn.push_back(n.normalize());
            }
        }
    }
    else
    {
        if (!dpx_feature.empty())
        {
            assert(feature_edge_flag.size() != 0);
            for (size_t i = 0; i < feature_edge_flag.size(); i++)
            {
                if (feature_edge_flag[i] == true)
                {
                    OpenVolumeMesh::EdgeHandle eh = OpenVolumeMesh::EdgeHandle(i);
                    OpenVolumeMesh::OpenVolumeMeshEdge ovme = mesh_->edge(eh);
                    int id0 = ovme.from_vertex().idx();
                    int id1 = ovme.to_vertex().idx();
                    OpenVolumeMesh::Geometry::Vec3d p0;
                    OpenVolumeMesh::Geometry::Vec3d p1;
                    p0[0] = dpx_feature[id0];
                    p0[1] = dpy_feature[id0];
                    p0[2] = dpz_feature[id0];
                    p1[0] = dpx_feature[id1];
                    p1[1] = dpy_feature[id1];
                    p1[2] = dpz_feature[id1];
                    OpenVolumeMesh::Geometry::Vec3d n = p1 - p0;
                    all_bfn.push_back(n.normalize());
                }
            }
        }
    }
    std::vector<double> X(3, 0.0);
    for (int i = 0; i < 1; ++i)
    {
        ad_LBFGS(all_bfn, X);
    }
    double alpha = X[0]; double beta = X[1]; double gamma = X[2];
    printf("%f %f %f\n", alpha, beta, gamma);
    if (alpha > std::numeric_limits<double>::min() && alpha < std::numeric_limits<double>::max() && gamma > std::numeric_limits<double>::min() && gamma < std::numeric_limits<double>::max() && beta > std::numeric_limits<double>::min() && beta < std::numeric_limits<double>::max())
    {
        double cos_phi = std::cos(X[0]), sin_phi = std::sin(X[0]);
        double cos_theta = std::cos(X[1]), sin_theta = std::sin(X[1]);
        double cos_xi = std::cos(X[2]), sin_xi = std::sin(X[2]);
        Matrix3d Mtheta, Mphi, Msi;
        Mtheta.setZero(); Mphi.setZero(); Msi.setZero();
        Mphi(0, 0) = 1; Mphi(1, 1) = Mphi(2, 2) = cos_phi, Mphi(1, 2) = sin_phi, Mphi(2, 1) = -sin_phi;
        Mtheta(1, 1) = 1; Mtheta(0, 0) = Mtheta(2, 2) = cos_theta, Mtheta(2, 0) = sin_theta, Mtheta(0, 2) = -sin_theta;
        Msi(2, 2) = 1; Msi(0, 0) = Msi(1, 1) = cos_xi, Msi(0, 1) = sin_xi, Msi(1, 0) = -sin_xi;
        Matrix3d M = Msi * Mphi * Mtheta;
        OpenVolumeMesh::Geometry::Vec3d p;
        for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
        {
            OpenVolumeMesh::Geometry::Vec3d q = mesh_->vertex(*v_it);
            Vector3d qvec(q[0], q[1], q[2]);
            Vector3d pvec = M * qvec;
            p[0] = pvec[0];
            p[1] = pvec[1];
            p[2] = pvec[2];
            mesh_->set_vertex(*v_it, p);
            int v_id = v_it->idx();
            dpx[v_id] = p[0] ; dpy[v_id] = p[1] ; dpz[v_id] = p[2] ;
        }
        if (polyline_flag)
        {
            for (size_t i = 0; i < dpx_feature.size(); i++)
            {
                if (feature_v2e[i].size() > 0)
                {
                    Vector3d q;
                    q[0] = dpx_feature[i];
                    q[1] = dpy_feature[i];
                    q[2] = dpz_feature[i];
                    Vector3d pvec = M * q;
                    dpx_feature[i] = pvec[0];
                    dpy_feature[i] = pvec[1];
                    dpz_feature[i] = pvec[2];
                }
            }
        }
    }
}
void polycube_deformation_interface::adjust_orientation(TetStructure<double>* tet_mesh_)
{
    OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
    std::vector<OpenVolumeMesh::Geometry::Vec3d> all_bfn;
    for (int i = 0; i < bfv_id.size(); ++i)
    {
        OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
        if (one_bfv_id[0] < 0) continue;
        p0[0] = dpx[one_bfv_id[0]]; p0[1] = dpy[one_bfv_id[0]]; p0[2] = dpz[one_bfv_id[0]];
        p1[0] = dpx[one_bfv_id[1]]; p1[1] = dpy[one_bfv_id[1]]; p1[2] = dpz[one_bfv_id[1]];
        p2[0] = dpx[one_bfv_id[2]]; p2[1] = dpy[one_bfv_id[2]]; p2[2] = dpz[one_bfv_id[2]];
        OpenVolumeMesh::Geometry::Vec3d pt = (p0 + p1 + p2) / 3.0;
        OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0);
        all_bfn.push_back(n.normalize());
    }
    std::vector<double> X(3, 0.0);
    for (int i = 0; i < 1; ++i)
    {
        ad_LBFGS(all_bfn, X);
    }
    double alpha = X[0]; double beta = X[1]; double gamma = X[2];
    printf("%f %f %f\n", alpha, beta, gamma);
    double cos_phi = std::cos(X[0]), sin_phi = std::sin(X[0]);
    double cos_theta = std::cos(X[1]), sin_theta = std::sin(X[1]);
    double cos_xi = std::cos(X[2]), sin_xi = std::sin(X[2]);
    Matrix3d Mtheta, Mphi, Msi;
    Mtheta.setZero(); Mphi.setZero(); Msi.setZero();
    Mphi(0, 0) = 1; Mphi(1, 1) = Mphi(2, 2) = cos_phi, Mphi(1, 2) = sin_phi, Mphi(2, 1) = -sin_phi;
    Mtheta(1, 1) = 1; Mtheta(0, 0) = Mtheta(2, 2) = cos_theta, Mtheta(2, 0) = sin_theta, Mtheta(0, 2) = -sin_theta;
    Msi(2, 2) = 1; Msi(0, 0) = Msi(1, 1) = cos_xi, Msi(0, 1) = sin_xi, Msi(1, 0) = -sin_xi;
    Matrix3d M = Msi * Mphi * Mtheta;
    CVec<double, 3> p;
    int nv = tet_mesh_->tetra_vertices.size();
    std::vector<TetVertex<double> *> &tetra_vertices = tet_mesh_->tetra_vertices;
    for (size_t i = 0; i < nv; i++)
    {
        CVec<double, 3> q = tetra_vertices[i]->pos;
        Vector3d qvec(q[0], q[1], q[2]);
        Vector3d pvec = M * qvec;
        p[0] = pvec[0];
        p[1] = pvec[1];
        p[2] = pvec[2];
        tetra_vertices[i]->pos = p;
        int v_id = i;
        dpx[v_id] = p[0] ; dpy[v_id] = p[1] ; dpz[v_id] = p[2] ;
    }
}
void polycube_deformation_interface::exp_mips_deformation_refine_one_polycube(int v_id, double angle_area_ratio, double energy_power, const int& omp_id, double ga, bool update_all, bool feature_flag)
{
    bool bv_flag = false;
    if (bvf_id[v_id].size() > 0)
    {
        bv_flag = true;
    }
    else
    {
        if (distortion_big_count[v_id] == 0)
            return;
    }
    const double* p_dpx = dpx.data(); const double* p_dpy = dpy.data(); const double* p_dpz = dpz.data();
    double x0 = p_dpx[v_id]; double y0 = p_dpy[v_id]; double z0 = p_dpz[v_id];
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    double local_energy = 0.0; double gx = 0.0; double gy = 0.0; double gz = 0.0;
    double min_radius = 1e30; double len;
    unsigned jj = 0; double alpha = (1.0 - angle_area_ratio)*energy_power; double beta = energy_power - alpha;
    int vc_size = vertex_cell[v_id].size(); const int* vc_id = vertex_cell[v_id].data();
    int v_id_1, v_id_2, v_id_3, k;
    double D00, D10, D20, D01, D11, D21, D02, D12, D22;
    double C00, C01, C02, C10, C11, C12, C20, C21, C22;
    double A00, A01, A02, A10, A11, A12, A20, A21, A22;
    double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
    double AF, AF_05, I_AF_05, AF2, AF_I, AF_I_05, I_AF_I_05, exp_k, exp_e;
    double D_AF_x, D_AF_y, D_AF_z, D_AF2_x, D_AF2_y, D_AF2_z, D_AF_I_x, D_AF_I_y, D_AF_I_z, D_det_A_x, D_det_A_y, D_det_A_z, inv_det_A_2_05;
    double d_A00_x, d_A10_y, d_A20_z, d_A01_x, d_A11_y, d_A21_z, d_A02_x, d_A12_y, d_A22_z;
    double dvex, dvey, dvez, g, dgx, dgy, dgz, dex, dey, dez, e, mips_e, volume_e, d_mips_e_x, d_mips_e_y, d_mips_e_z;
    double* p_vc_pos_x = vc_pos_x2_omp[omp_id].data();
    double* p_vc_pos_y = vc_pos_y2_omp[omp_id].data();
    double* p_vc_pos_z = vc_pos_z2_omp[omp_id].data();
    double* p_vc_n_cross_x = vc_n_cross_x_omp[omp_id].data();
    double* p_vc_n_cross_y = vc_n_cross_y_omp[omp_id].data();
    double* p_vc_n_cross_z = vc_n_cross_z_omp[omp_id].data();
    double* exp_vec = exp_vec_omp[omp_id].data();
    double* gx_vec = gx_vec_omp[omp_id].data();
    double* gy_vec = gy_vec_omp[omp_id].data();
    double* gz_vec = gz_vec_omp[omp_id].data();
    double* mu_vec = mu_vec_omp[omp_id].data();
    std::vector<int> feature_edge_face_idx;
    std::vector<std::pair<int, int>> feature_neighber_vert_localidx_pairidx;
    std::vector<std::vector<double>>& vc_S = vc_S_omp[omp_id];
    for (unsigned i = 0; i < vc_size; ++i)
    {
        const int* vv_id = vertex_cell_vertex[v_id][i].data();
        memcpy(&vc_S[i][0], &vcv_S[v_id][i][0], 9 * sizeof(double));
        double* s_data = vc_S[i].data();
        k = 3 * i;
        v_id_1 = vv_id[0];
        v_id_2 = vv_id[1];
        v_id_3 = vv_id[2];
        x1 = p_dpx[v_id_1]; x2 = p_dpx[v_id_2]; x3 = p_dpx[v_id_3];
        y1 = p_dpy[v_id_1]; y2 = p_dpy[v_id_2]; y3 = p_dpy[v_id_3];
        z1 = p_dpz[v_id_1]; z2 = p_dpz[v_id_2]; z3 = p_dpz[v_id_3];
        p_vc_pos_x[k + 0] = x1; p_vc_pos_x[k + 1] = x2; p_vc_pos_x[k + 2] = x3;
        p_vc_pos_y[k + 0] = y1; p_vc_pos_y[k + 1] = y2; p_vc_pos_y[k + 2] = y3;
        p_vc_pos_z[k + 0] = z1; p_vc_pos_z[k + 1] = z2; p_vc_pos_z[k + 2] = z3;
        p_vc_n_cross_x[i] = (y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1);
        p_vc_n_cross_y[i] = (x3 - x1)*(z2 - z1) - (x2 - x1)*(z3 - z1);
        p_vc_n_cross_z[i] = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
        D00 = x1 - x0; D10 = y1 - y0; D20 = z1 - z0;
        D01 = x2 - x0; D11 = y2 - y0; D21 = z2 - z0;
        D02 = x3 - x0; D12 = y3 - y0; D22 = z3 - z0;
        len = D00 * D00 + D10 * D10 + D20 * D20; if (len < min_radius) min_radius = len;
        len = D01 * D01 + D11 * D11 + D21 * D21; if (len < min_radius) min_radius = len;
        len = D02 * D02 + D12 * D12 + D22 * D22; if (len < min_radius) min_radius = len;
        C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
        C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
        C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];
        A00 = C00 * D00 + C10 * D01 + C20 * D02; d_A00_x = -C00 - C10 - C20;
        A10 = C00 * D10 + C10 * D11 + C20 * D12; d_A10_y = -C00 - C10 - C20;
        A20 = C00 * D20 + C10 * D21 + C20 * D22; d_A20_z = -C00 - C10 - C20;
        A01 = C01 * D00 + C11 * D01 + C21 * D02; d_A01_x = -C01 - C11 - C21;
        A11 = C01 * D10 + C11 * D11 + C21 * D12; d_A11_y = -C01 - C11 - C21;
        A21 = C01 * D20 + C11 * D21 + C21 * D22; d_A21_z = -C01 - C11 - C21;
        A02 = C02 * D00 + C12 * D01 + C22 * D02; d_A02_x = -C02 - C12 - C22;
        A12 = C02 * D10 + C12 * D11 + C22 * D12; d_A12_y = -C02 - C12 - C22;
        A22 = C02 * D20 + C12 * D21 + C22 * D22; d_A22_z = -C02 - C12 - C22;
        A2_00 = A00 * A00 + A10 * A10 + A20 * A20;
        A2_01 = A00 * A01 + A10 * A11 + A20 * A21;
        A2_02 = A00 * A02 + A10 * A12 + A20 * A22;
        A2_11 = A01 * A01 + A11 * A11 + A21 * A21;
        A2_12 = A01 * A02 + A11 * A12 + A21 * A22;
        A2_22 = A02 * A02 + A12 * A12 + A22 * A22;
        AF = A2_00 + A2_11 + A2_22; AF_05 = std::sqrt(AF); I_AF_05 = 0.5 / AF_05;
        AF2 = A2_00 * A2_00 + A2_01 * A2_01 * 2 + A2_02 * A2_02 * 2 + A2_11 * A2_11 + A2_12 * A2_12 * 2 + A2_22 * A2_22;
        AF_I = (AF*AF - AF2)*0.5; AF_I_05 = std::sqrt(AF_I); I_AF_I_05 = 0.5 / AF_I_05;
        det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
        i_det_A = 1.0 / det_A;
        D_AF_x = (2.0*A00*d_A00_x + 2.0*A01*d_A01_x + 2.0*A02*d_A02_x);
        D_AF_y = (2.0*A10*d_A10_y + 2.0*A11*d_A11_y + 2.0*A12*d_A12_y);
        D_AF_z = (2.0*A20*d_A20_z + 2.0*A21*d_A21_z + 2.0*A22*d_A22_z);
        D_AF2_x = (2.0*A2_00*2.0*A00*d_A00_x + 4.0*A2_01*(A01*d_A00_x + A00 * d_A01_x) + 4.0*A2_02*(A02*d_A00_x + A00 * d_A02_x)
            + 2.0*A2_11*2.0*A01*d_A01_x + 2.0*A2_22*2.0*A02*d_A02_x + 4.0*A2_12*(A01*d_A02_x + A02 * d_A01_x));
        D_AF2_y = (2.0*A2_00*2.0*A10*d_A10_y + 4.0*A2_01*(A11*d_A10_y + A10 * d_A11_y) + 4.0*A2_02*(A12*d_A10_y + A10 * d_A12_y)
            + 2.0*A2_11*2.0*A11*d_A11_y + 2.0*A2_22*2.0*A12*d_A12_y + 4.0*A2_12*(A11*d_A12_y + A12 * d_A11_y));
        D_AF2_z = (2.0*A2_00*2.0*A20*d_A20_z + 4.0*A2_01*(A21*d_A20_z + A20 * d_A21_z) + 4.0*A2_02*(A22*d_A20_z + A20 * d_A22_z)
            + 2.0*A2_11*2.0*A21*d_A21_z + 2.0*A2_22*2.0*A22*d_A22_z + 4.0*A2_12*(A21*d_A22_z + A22 * d_A21_z));
        D_AF_I_x = (AF*D_AF_x - 0.5*D_AF2_x)*I_AF_I_05; D_AF_I_y = (AF*D_AF_y - 0.5*D_AF2_y)*I_AF_I_05; D_AF_I_z = (AF*D_AF_z - 0.5*D_AF2_z)*I_AF_I_05;
        D_AF_x *= I_AF_05; D_AF_y *= I_AF_05; D_AF_z *= I_AF_05;
        D_det_A_x = d_A00_x * A11 * A22 + d_A01_x * A12 * A20 + d_A02_x * A10 * A21 - d_A02_x * A11 * A20 - d_A01_x * A10 * A22 - d_A00_x * A12 * A21;
        D_det_A_y = A00 * d_A11_y * A22 + A01 * d_A12_y * A20 + A02 * d_A10_y * A21 - A02 * d_A11_y * A20 - A01 * d_A10_y * A22 - A00 * d_A12_y * A21;
        D_det_A_z = A00 * A11 * d_A22_z + A01 * A12 * d_A20_z + A02 * A10 * d_A21_z - A02 * A11 * d_A20_z - A01 * A10 * d_A22_z - A00 * A12 * d_A21_z;
        inv_det_A_2_05 = 0.5 *(1.0 - i_det_A * i_det_A);
        dvex = D_det_A_x * inv_det_A_2_05;
        dvey = D_det_A_y * inv_det_A_2_05;
        dvez = D_det_A_z * inv_det_A_2_05;
        g = AF_05 * AF_I_05;
        dgx = D_AF_x * AF_I_05 + AF_05 * D_AF_I_x;
        dgy = D_AF_y * AF_I_05 + AF_05 * D_AF_I_y;
        dgz = D_AF_z * AF_I_05 + AF_05 * D_AF_I_z;
        dex = (dgx *i_det_A - g * D_det_A_x *i_det_A*i_det_A);
        dey = (dgy *i_det_A - g * D_det_A_y *i_det_A*i_det_A);
        dez = (dgz *i_det_A - g * D_det_A_z *i_det_A*i_det_A);
        e = g * i_det_A;
        mips_e = (e*e - 1.0)*0.125;
        volume_e = 0.5*(det_A + i_det_A);
        exp_k = (alpha*mips_e + beta * volume_e);
        if (exp_k > 60) exp_k = 60;
        exp_vec[i] = exp_k;
        d_mips_e_x = e * 0.25 * dex;
        d_mips_e_y = e * 0.25 * dey;
        d_mips_e_z = e * 0.25 * dez;
        gx_vec[i] = alpha * d_mips_e_x + beta * dvex;
        gy_vec[i] = alpha * d_mips_e_y + beta * dvey;
        gz_vec[i] = alpha * d_mips_e_z + beta * dvez;
    }
#if 0
    fmath::expd_v(exp_vec, vc_size);
    for (int i = 0; i < vc_size; ++i)
    {
        exp_e = exp_vec[i];
        local_energy += exp_e;
        gx += gx_vec[i] * exp_e;
        gy += gy_vec[i] * exp_e;
        gz += gz_vec[i] * exp_e;
    }
#else
    for (int i = 0; i < vc_size; ++i)
    {
        exp_e = std::exp(exp_vec[i]);
        local_energy += exp_e;
        gx += gx_vec[i] * exp_e;
        gy += gy_vec[i] * exp_e;
        gz += gz_vec[i] * exp_e;
    }
#endif
    double mu = 1.0; int bvf_size = 0; int bvv_size = 0;
    double* p_vf_pos_x = NULL; double* p_vf_pos_y = NULL; double* p_vf_pos_z = NULL;
    if (bv_flag)
    {
        double normal_e = 0.0; double smooth_e = 0;
        double ne_gx = 0.0; double ne_gy = 0.0; double ne_gz = 0.0;
        double se_gx = 0.0; double se_gy = 0.0; double se_gz = 0.0;
        p_vf_pos_x = vc_pos_x3_omp[omp_id].data();
        p_vf_pos_y = vc_pos_y3_omp[omp_id].data();
        p_vf_pos_z = vc_pos_z3_omp[omp_id].data();
        std::vector<int>& one_vv_id = bvv_id[v_id];
        std::vector<int>& one_vf_id = bvf_id[v_id];
        bvv_size = one_vv_id.size(); bvf_size = bvv_size;
        for (int i = 0; i < bvv_size; ++i)
        {
            int j = (i + 1) % bvv_size; int k = (i + 2) % bvv_size;
            x1 = p_dpx[one_vv_id[i]]; y1 = p_dpy[one_vv_id[i]]; z1 = p_dpz[one_vv_id[i]];
            x2 = p_dpx[one_vv_id[j]]; y2 = p_dpy[one_vv_id[j]]; z2 = p_dpz[one_vv_id[j]];
            x3 = p_dpx[one_vv_id[k]]; y3 = p_dpy[one_vv_id[k]]; z3 = p_dpz[one_vv_id[k]];
            p_vf_pos_x[4 * i + 0] = x1; p_vf_pos_x[4 * i + 1] = x2; p_vf_pos_x[4 * i + 2] = x3;
            p_vf_pos_y[4 * i + 0] = y1; p_vf_pos_y[4 * i + 1] = y2; p_vf_pos_y[4 * i + 2] = y3;
            p_vf_pos_z[4 * i + 0] = z1; p_vf_pos_z[4 * i + 1] = z2; p_vf_pos_z[4 * i + 2] = z3;
            OpenVolumeMesh::Geometry::Vec3d& tn = target_bfn[one_vf_id[i]];
            p_vf_pos_x[4 * i + 3] = tn[0]; p_vf_pos_y[4 * i + 3] = tn[1]; p_vf_pos_z[4 * i + 3] = tn[2];
            double nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
            double ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
            double nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);
            double len2 = nx * nx + ny * ny + nz * nz; double len = std::sqrt(len2);
            double E_ne = 0.5*((nx / len - tn[0])*(nx / len - tn[0]) + (ny / len - tn[1])*(ny / len - tn[1]) + (nz / len - tn[2])*(nz / len - tn[2]));
            if (update_all)
            {
                mu_vec[i] = local_energy_ratio[one_vf_id[i]];
            }
            else
            {
                mu_vec[i] = 1.0;
            }
            normal_e += mu_vec[i] * E_ne;
            double d_nx_x0 = 0.0; double d_nx_y0 = z1 - z2; double d_nx_z0 = y2 - y1;
            double d_ny_x0 = z2 - z1; double d_ny_y0 = 0.0; double d_ny_z0 = x1 - x2;
            double d_nz_x0 = y1 - y2; double d_nz_y0 = x2 - x1; double d_nz_z0 = 0.0;
            double d_len2_x0 = (nx*d_nx_x0 + ny * d_ny_x0 + nz * d_nz_x0);
            double d_len2_y0 = (nx*d_nx_y0 + ny * d_ny_y0 + nz * d_nz_y0);
            double d_len2_z0 = (nx*d_nx_z0 + ny * d_ny_z0 + nz * d_nz_z0);
            double d_len_x0 = d_len2_x0 / len; double d_len_y0 = d_len2_y0 / len; double d_len_z0 = d_len2_z0 / len;
            double d_nx_len_x0 = d_nx_x0 / len - nx * d_len_x0 / (len*len);
            double d_nx_len_y0 = d_nx_y0 / len - nx * d_len_y0 / (len*len);
            double d_nx_len_z0 = d_nx_z0 / len - nx * d_len_z0 / (len*len);
            double d_ny_len_x0 = d_ny_x0 / len - ny * d_len_x0 / (len*len);
            double d_ny_len_y0 = d_ny_y0 / len - ny * d_len_y0 / (len*len);
            double d_ny_len_z0 = d_ny_z0 / len - ny * d_len_z0 / (len*len);
            double d_nz_len_x0 = d_nz_x0 / len - nz * d_len_x0 / (len*len);
            double d_nz_len_y0 = d_nz_y0 / len - nz * d_len_y0 / (len*len);
            double d_nz_len_z0 = d_nz_z0 / len - nz * d_len_z0 / (len*len);
            ne_gx += mu_vec[i] * (d_nx_len_x0*(nx / len - tn[0]) + d_ny_len_x0 * (ny / len - tn[1]) + d_nz_len_x0 * (nz / len - tn[2]));
            ne_gy += mu_vec[i] * (d_nx_len_y0*(nx / len - tn[0]) + d_ny_len_y0 * (ny / len - tn[1]) + d_nz_len_y0 * (nz / len - tn[2]));
            ne_gz += mu_vec[i] * (d_nx_len_z0*(nx / len - tn[0]) + d_ny_len_z0 * (ny / len - tn[1]) + d_nz_len_z0 * (nz / len - tn[2]));
        }
        if (feature_flag)
        {
            if (feature_v2e[v_id].size() != 0)
            {
                for (size_t i = 0; i < feature_v2e[v_id].size(); i++)
                {
                    int bf0 = bef_id[feature_v2e[v_id][i]][0];
                    int bf1 = bef_id[feature_v2e[v_id][i]][1];
                    double tmp_mu = 0.0;
                    if (update_all)
                    {
                        tmp_mu = FEATURE_COEFF;
                    }
                    else
                    {
                        tmp_mu = 1.0;
                    }
                    int bf0_id(-1), bf1_id(-1);
                    for (size_t j = 0; j < bvf_size; j++)
                    {
                        if (one_vf_id[j] == bf0)
                            bf0_id = j;
                        if (one_vf_id[j] == bf1)
                            bf1_id = j;
                    }
                    assert(bf0_id != -1 && bf1_id != -1 && bf0_id != bf1_id);
                    feature_edge_face_idx.push_back(bf0_id);
                    feature_edge_face_idx.push_back(bf1_id);
                    int v1_id[2] = { bvv_id[v_id][bf0_id], bvv_id[v_id][bf1_id] };
                    int v2_id[2] = { bvv_id[v_id][(bf0_id + 1) % bvv_size], bvv_id[v_id][(bf1_id + 1) % bvv_size] };
                    double nx[2], ny[2], nz[2];
                    double nxn[2], nyn[2], nzn[2];
                    double x1[2], x2[2], y1[2], y2[2], z1[2], z2[2];
                    double len2[2], len[2];
                    for (size_t j = 0; j < 2; j++)
                    {
                        x1[j] = p_dpx[v1_id[j]];
                        y1[j] = p_dpy[v1_id[j]];
                        z1[j] = p_dpz[v1_id[j]];
                        x2[j] = p_dpx[v2_id[j]];
                        y2[j] = p_dpy[v2_id[j]];
                        z2[j] = p_dpz[v2_id[j]];
                        nx[j] = (y0 - y1[j])*(z0 - z2[j]) - (y0 - y2[j])*(z0 - z1[j]);
                        ny[j] = (x0 - x2[j])*(z0 - z1[j]) - (x0 - x1[j])*(z0 - z2[j]);
                        nz[j] = (x0 - x1[j])*(y0 - y2[j]) - (x0 - x2[j])*(y0 - y1[j]);
                        len2[j] = nx[j] * nx[j] + ny[j] * ny[j] + nz[j] * nz[j];
                        len[j] = std::sqrt(len2[j]);
                        nxn[j] = nx[j] / len[j];
                        nyn[j] = ny[j] / len[j];
                        nzn[j] = nz[j] / len[j];
                    }
                    double n1n2 = nxn[0] * nxn[1] + nyn[0] * nyn[1] + nzn[0] * nzn[1];
                    double n1n22 = n1n2 * n1n2;
                    normal_e += tmp_mu * n1n22;
                    double dx_nxn[2], dy_nyn[2], dz_nzn[2], dy_nxn[2], dz_nxn[2], dx_nyn[2], dz_nyn[2], dx_nzn[2], dy_nzn[2];
                    for (size_t j = 0; j < 2; j++)
                    {
                        double tmp_coord[3][3] = {
                            {x0, y0, z0},
                            {x1[j], y1[j], z1[j]},
                            {x2[j], y2[j], z2[j]},
                        };
                        double normal_coord[3] = { nx[j], ny[j], nz[j] };
                        double der_normal[3][3] = {
                            {0.0, z1[j] - z2[j], y2[j] - y1[j] },
                            {z2[j] - z1[j], 0.0, x1[j] - x2[j] },
                            {y1[j] - y2[j], x2[j] - x1[j], 0.0 },
                        };
                        double der_len[3];
                        for (size_t a = 0; a < 3; a++)
                        {
                            der_len[a] = (der_normal[0][a] * normal_coord[0] + der_normal[1][a] * normal_coord[1] + der_normal[2][a] * normal_coord[2]) / len[j];
                        }
                        double der_normaln[3][3];
                        for (size_t a = 0; a < 3; a++)
                        {
                            for (size_t b = 0; b < 3; b++)
                            {
                                der_normaln[a][b] = (der_normal[a][b] * len[j] - normal_coord[a] * der_len[b]) / len2[j];
                            }
                        }
                        dx_nxn[j] = der_normaln[0][0]; dy_nxn[j] = der_normaln[0][1]; dz_nxn[j] = der_normaln[0][2];
                        dx_nyn[j] = der_normaln[1][0]; dy_nyn[j] = der_normaln[1][1]; dz_nyn[j] = der_normaln[1][2];
                        dx_nzn[j] = der_normaln[2][0]; dy_nzn[j] = der_normaln[2][1]; dz_nzn[j] = der_normaln[2][2];
                    }
                    double dx_n1n2, dy_n1n2, dz_n1n2;
                    dx_n1n2 = nyn[0] * dx_nyn[1] + dx_nyn[0] * nyn[1] + nzn[0] * dx_nzn[1] + dx_nzn[0] * nzn[1] + nxn[0] * dx_nxn[1] + dx_nxn[0] * nxn[1];
                    dy_n1n2 = nxn[0] * dy_nxn[1] + dy_nxn[0] * nxn[1] + nzn[0] * dy_nzn[1] + dy_nzn[0] * nzn[1] + nyn[0] * dy_nyn[1] + dy_nyn[0] * nyn[1];
                    dz_n1n2 = nxn[0] * dz_nxn[1] + dz_nxn[0] * nxn[1] + nyn[0] * dz_nyn[1] + dz_nyn[0] * nyn[1] + nzn[0] * dz_nzn[1] + dz_nzn[0] * nzn[1];
                    ne_gx += 2 * tmp_mu * n1n2 * dx_n1n2;
                    ne_gy += 2 * tmp_mu * n1n2 * dy_n1n2;
                    ne_gz += 2 * tmp_mu * n1n2 * dz_n1n2;
                }
                assert(feature_edge_face_idx.size() >= 4);
            }
            if (feature_neighbor_vert2cellpair[v_id].size() != 0)
            {
                for (size_t i = 0; i < feature_neighbor_vert2cellpair[v_id].size(); i++)
                {
                    int bf0 = feature_neighbor_vert2cellpair[v_id][i].first;
                    int bf1 = feature_neighbor_vert2cellpair[v_id][i].second;
                    double tmp_mu = FEATURE_COEFF;
                    int bf0_id(-1);
                    for (size_t j = 0; j < bvf_size; j++)
                    {
                        if (one_vf_id[j] == bf0)
                            bf0_id = j;
                    }
                    assert(bf0_id != -1);
                    feature_neighber_vert_localidx_pairidx.push_back(std::pair<int, int>(bf0_id, bf1));
                    int v1_id = bvv_id[v_id][bf0_id];
                    int v2_id = bvv_id[v_id][(bf0_id + 1) % bvv_size];
                    double nx, ny, nz, nxn, nyn, nzn;
                    double len2, len;
                    x1 = p_dpx[v1_id];
                    y1 = p_dpy[v1_id];
                    z1 = p_dpz[v1_id];
                    x2 = p_dpx[v2_id];
                    y2 = p_dpy[v2_id];
                    z2 = p_dpz[v2_id];
                    nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
                    ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
                    nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);
                    len2 = nx * nx + ny * ny + nz * nz;
                    len = std::sqrt(len2);
                    nxn = nx / len;
                    nyn = ny / len;
                    nzn = nz / len;
                    double nxn2, nyn2, nzn2;
                    nxn2 = last_bfn[bf1][0];
                    nyn2 = last_bfn[bf1][1];
                    nzn2 = last_bfn[bf1][2];
                    double n1n2 = nxn * nxn2 + nyn * nyn2 + nzn * nzn2;
                    normal_e += tmp_mu * n1n2 * n1n2;
                    double dx_nxn[2], dy_nyn[2], dz_nzn[2], dy_nxn[2], dz_nxn[2], dx_nyn[2], dz_nyn[2], dx_nzn[2], dy_nzn[2];
                    for (size_t j = 0; j < 1; j++)
                    {
                        double tmp_coord[3][3] = {
                            {x0, y0, z0},
                            {x1, y1, z1},
                            {x2, y2, z2},
                        };
                        double normal_coord[3] = { nx, ny, nz };
                        double der_normal[3][3] = {
                            {0.0, z1 - z2, y2 - y1 },
                            {z2 - z1, 0.0, x1 - x2 },
                            {y1 - y2, x2 - x1, 0.0 },
                        };
                        double der_len[3];
                        for (size_t a = 0; a < 3; a++)
                        {
                            der_len[a] = der_len[a] = (der_normal[0][a] * normal_coord[0] + der_normal[1][a] * normal_coord[1] + der_normal[2][a] * normal_coord[2]) / len;
                        }
                        double der_normaln[3][3];
                        for (size_t a = 0; a < 3; a++)
                        {
                            for (size_t b = 0; b < 3; b++)
                            {
                                der_normaln[a][b] = (der_normal[a][b] * len - normal_coord[a] * der_len[b]) / len2;
                            }
                        }
                        dx_nxn[j] = der_normaln[0][0]; dy_nxn[j] = der_normaln[0][1]; dz_nxn[j] = der_normaln[0][2];
                        dx_nyn[j] = der_normaln[1][0]; dy_nyn[j] = der_normaln[1][1]; dz_nyn[j] = der_normaln[1][2];
                        dx_nzn[j] = der_normaln[2][0]; dy_nzn[j] = der_normaln[2][1]; dz_nzn[j] = der_normaln[2][2];
                    }
                    double dx_n1n2, dy_n1n2, dz_n1n2;
                    dx_n1n2 = dx_nyn[0] * nyn2 + dx_nzn[0] * nzn2 + dx_nxn[0] * nxn2;
                    dy_n1n2 = dy_nxn[0] * nxn2 + dy_nzn[0] * nzn2 + dy_nyn[0] * nyn2;
                    dz_n1n2 = dz_nxn[0] * nxn2 + dz_nyn[0] * nyn2 + dz_nzn[0] * nzn2;
                    ne_gx += 2 * tmp_mu * n1n2 * dx_n1n2;
                    ne_gy += 2 * tmp_mu * n1n2 * dy_n1n2;
                    ne_gz += 2 * tmp_mu * n1n2 * dz_n1n2;
                }
            }
        }
        if (update_all)
        {
            local_energy += normal_e;
            gx += ne_gx; gy += ne_gy; gz += ne_gz;
        }
        else
        {
            double s = local_energy / (normal_e + 1e-8);
            {
                ga = s;
            }
            if (s > 1e16) ga = 1e16;
            for (int i = 0; i < bvf_size; ++i)
            {
                mu_vec[i] = ga;
            }
            local_energy += ga * normal_e;
            gx += ga * ne_gx; gy += ga * ne_gy; gz += ga * ne_gz;
        }
    }
    min_radius = std::sqrt(min_radius);
    double dx = gx, dy = gy, dz = gz;
    double move_d = sqrt(dx*dx + dy * dy + dz * dz);
    if (move_d > min_radius)
    {
        double t = min_radius / move_d;
        dx *= t; dy *= t; dz *= t;
        move_d = min_radius;
    }
    double npx = x0 - dx; double npy = y0 - dy; double npz = z0 - dz;
    while (!local_check_negative_volume4(vc_size, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, p_vc_n_cross_x, p_vc_n_cross_y, p_vc_n_cross_z, npx, npy, npz))
    {
        dx *= 0.8; dy *= 0.8; dz *= 0.8; move_d *= 0.8;
        if (move_d < 1e-8)
        {
            npx = x0; npy = y0; npz = z0;
            move_d = 1e-8; break;
        }
        else
        {
            npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
        }
    }
#if 1
    double new_e = 0;
    bool e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec, feature_edge_face_idx, feature_neighber_vert_localidx_pairidx);
    if (!e_ok)
    {
        while (!e_ok)
        {
            dx *= 0.2; dy *= 0.2; dz *= 0.2; move_d *= 0.2;
            if (move_d < 1e-8)
            {
                npx = x0; npy = y0; npz = z0;
                move_d = 1e-8; break;
            }
            else
            {
                npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
            }
            new_e = 0;
            e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec, feature_edge_face_idx, feature_neighber_vert_localidx_pairidx);
        }
        if (e_ok)
        {
            dx *= 5; dy *= 5; dz *= 5; move_d *= 5; e_ok = false;
            while (!e_ok)
            {
                dx *= 0.5; dy *= 0.5; dz *= 0.5; move_d *= 0.5;
                if (move_d < 1e-8)
                {
                    npx = x0; npy = y0; npz = z0;
                    move_d = 1e-8; break;
                }
                else
                {
                    npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
                }
                new_e = 0;
                e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec, feature_edge_face_idx, feature_neighber_vert_localidx_pairidx);
            }
            if (e_ok)
            {
                dx *= 2; dy *= 2; dz *= 2; move_d *= 2; e_ok = false;
                while (!e_ok)
                {
                    dx *= 0.875; dy *= 0.875; dz *= 0.875; move_d *= 0.875;
                    if (move_d < 1e-8)
                    {
                        npx = x0; npy = y0; npz = z0;
                        move_d = 1e-8; break;
                    }
                    else
                    {
                        npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
                    }
                    new_e = 0;
                    e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec, feature_edge_face_idx, feature_neighber_vert_localidx_pairidx);
                }
            }
        }
    }
#else
    double new_e = 0;
    bool e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
    while (!e_ok)
    {
        dx *= 0.8; dy *= 0.8; dz *= 0.8; move_d *= 0.8;
        if (move_d < 1e-8)
        {
            npx = x0; npy = y0; npz = z0;
            move_d = 1e-8; break;
        }
        else
        {
            npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
        }
        new_e = 0;
        e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
    }
#endif
    dpx[v_id] = npx; dpy[v_id] = npy; dpz[v_id] = npz;
    check_big_distortion_cell(vc_size, large_dis_flag_omp[omp_id], p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, npx, npy, npz);
    for (int i = 0; i < vc_size; ++i)
    {
        int c_id = vertex_cell[v_id][i];
        if (large_dis_flag_omp[omp_id][i] == 1)
        {
            if (distortion_big_cell[c_id] == 1)
            {
            }
            else
            {
                distortion_big_cell[c_id] = 1;
                for (int j = 0; j < 4; ++j)
                {
                    distortion_big_count_omp[omp_id][cell_vertex[c_id][j]] += 1;
                }
            }
        }
        else
        {
            if (distortion_big_cell[c_id] == 1)
            {
                distortion_big_cell[c_id] = -1;
                for (int j = 0; j < 4; ++j)
                {
                    distortion_big_count_omp[omp_id][cell_vertex[c_id][j]] -= 1;
                }
            }
            else
            {
            }
        }
    }
}
void polycube_deformation_interface::exp_mips_deformation_refine_one_polycube_normal_constrained(int v_id, double angle_area_ratio, double energy_power, const int& omp_id, double ga, bool update_all)
{
    bool bv_flag = false;
    if (bvf_id[v_id].size() > 0)
    {
        bv_flag = true;
    }
    else
    {
        if (distortion_big_count[v_id] == 0)
            return;
    }
    double gx_all = 0.0; double gy_all = 0.0; double gz_all = 0.0;
    double min_radius = 1e30;
    const double* p_dpx = dpx.data(); const double* p_dpy = dpy.data(); const double* p_dpz = dpz.data();
    double x0 = p_dpx[v_id]; double y0 = p_dpy[v_id]; double z0 = p_dpz[v_id];
    int vc_size = vertex_cell[v_id].size();
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    double local_energy = 0.0; double gx = 0.0; double gy = 0.0; double gz = 0.0;
    double len;
    unsigned jj = 0; double alpha = (1.0 - angle_area_ratio)*energy_power; double beta = energy_power - alpha;
    const int* vc_id = vertex_cell[v_id].data();
    int v_id_1, v_id_2, v_id_3, k;
    double D00, D10, D20, D01, D11, D21, D02, D12, D22;
    double C00, C01, C02, C10, C11, C12, C20, C21, C22;
    double A00, A01, A02, A10, A11, A12, A20, A21, A22;
    double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
    double AF, AF_05, I_AF_05, AF2, AF_I, AF_I_05, I_AF_I_05, exp_k, exp_e;
    double D_AF_x, D_AF_y, D_AF_z, D_AF2_x, D_AF2_y, D_AF2_z, D_AF_I_x, D_AF_I_y, D_AF_I_z, D_det_A_x, D_det_A_y, D_det_A_z, inv_det_A_2_05;
    double d_A00_x, d_A10_y, d_A20_z, d_A01_x, d_A11_y, d_A21_z, d_A02_x, d_A12_y, d_A22_z;
    double dvex, dvey, dvez, g, dgx, dgy, dgz, dex, dey, dez, e, mips_e, volume_e, d_mips_e_x, d_mips_e_y, d_mips_e_z;
    double* p_vc_pos_x = vc_pos_x2_omp[omp_id].data();
    double* p_vc_pos_y = vc_pos_y2_omp[omp_id].data();
    double* p_vc_pos_z = vc_pos_z2_omp[omp_id].data();
    double* p_vc_n_cross_x = vc_n_cross_x_omp[omp_id].data();
    double* p_vc_n_cross_y = vc_n_cross_y_omp[omp_id].data();
    double* p_vc_n_cross_z = vc_n_cross_z_omp[omp_id].data();
    double* exp_vec = exp_vec_omp[omp_id].data();
    double* gx_vec = gx_vec_omp[omp_id].data();
    double* gy_vec = gy_vec_omp[omp_id].data();
    double* gz_vec = gz_vec_omp[omp_id].data();
    double* mu_vec = mu_vec_omp[omp_id].data();
    std::vector<std::vector<double>>& vc_S = vc_S_omp[omp_id];
    for (unsigned i = 0; i < vc_size; ++i)
    {
        const int* vv_id = vertex_cell_vertex[v_id][i].data();
        memcpy(&vc_S[i][0], &vcv_S[v_id][i][0], 9 * sizeof(double));
        double* s_data = vc_S[i].data();
        k = 3 * i;
        v_id_1 = vv_id[0];
        v_id_2 = vv_id[1];
        v_id_3 = vv_id[2];
        x1 = p_dpx[v_id_1]; x2 = p_dpx[v_id_2]; x3 = p_dpx[v_id_3];
        y1 = p_dpy[v_id_1]; y2 = p_dpy[v_id_2]; y3 = p_dpy[v_id_3];
        z1 = p_dpz[v_id_1]; z2 = p_dpz[v_id_2]; z3 = p_dpz[v_id_3];
        p_vc_pos_x[k + 0] = x1; p_vc_pos_x[k + 1] = x2; p_vc_pos_x[k + 2] = x3;
        p_vc_pos_y[k + 0] = y1; p_vc_pos_y[k + 1] = y2; p_vc_pos_y[k + 2] = y3;
        p_vc_pos_z[k + 0] = z1; p_vc_pos_z[k + 1] = z2; p_vc_pos_z[k + 2] = z3;
        p_vc_n_cross_x[i] = (y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1);
        p_vc_n_cross_y[i] = (x3 - x1)*(z2 - z1) - (x2 - x1)*(z3 - z1);
        p_vc_n_cross_z[i] = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
        D00 = x1 - x0; D10 = y1 - y0; D20 = z1 - z0;
        D01 = x2 - x0; D11 = y2 - y0; D21 = z2 - z0;
        D02 = x3 - x0; D12 = y3 - y0; D22 = z3 - z0;
        len = D00 * D00 + D10 * D10 + D20 * D20; if (len < min_radius) min_radius = len;
        len = D01 * D01 + D11 * D11 + D21 * D21; if (len < min_radius) min_radius = len;
        len = D02 * D02 + D12 * D12 + D22 * D22; if (len < min_radius) min_radius = len;
        C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
        C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
        C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];
        A00 = C00 * D00 + C10 * D01 + C20 * D02; d_A00_x = -C00 - C10 - C20;
        A10 = C00 * D10 + C10 * D11 + C20 * D12; d_A10_y = -C00 - C10 - C20;
        A20 = C00 * D20 + C10 * D21 + C20 * D22; d_A20_z = -C00 - C10 - C20;
        A01 = C01 * D00 + C11 * D01 + C21 * D02; d_A01_x = -C01 - C11 - C21;
        A11 = C01 * D10 + C11 * D11 + C21 * D12; d_A11_y = -C01 - C11 - C21;
        A21 = C01 * D20 + C11 * D21 + C21 * D22; d_A21_z = -C01 - C11 - C21;
        A02 = C02 * D00 + C12 * D01 + C22 * D02; d_A02_x = -C02 - C12 - C22;
        A12 = C02 * D10 + C12 * D11 + C22 * D12; d_A12_y = -C02 - C12 - C22;
        A22 = C02 * D20 + C12 * D21 + C22 * D22; d_A22_z = -C02 - C12 - C22;
        A2_00 = A00 * A00 + A10 * A10 + A20 * A20;
        A2_01 = A00 * A01 + A10 * A11 + A20 * A21;
        A2_02 = A00 * A02 + A10 * A12 + A20 * A22;
        A2_11 = A01 * A01 + A11 * A11 + A21 * A21;
        A2_12 = A01 * A02 + A11 * A12 + A21 * A22;
        A2_22 = A02 * A02 + A12 * A12 + A22 * A22;
        AF = A2_00 + A2_11 + A2_22; AF_05 = std::sqrt(AF); I_AF_05 = 0.5 / AF_05;
        AF2 = A2_00 * A2_00 + A2_01 * A2_01 * 2 + A2_02 * A2_02 * 2 + A2_11 * A2_11 + A2_12 * A2_12 * 2 + A2_22 * A2_22;
        AF_I = (AF*AF - AF2)*0.5; AF_I_05 = std::sqrt(AF_I); I_AF_I_05 = 0.5 / AF_I_05;
        det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
        i_det_A = 1.0 / det_A;
        D_AF_x = (2.0*A00*d_A00_x + 2.0*A01*d_A01_x + 2.0*A02*d_A02_x);
        D_AF_y = (2.0*A10*d_A10_y + 2.0*A11*d_A11_y + 2.0*A12*d_A12_y);
        D_AF_z = (2.0*A20*d_A20_z + 2.0*A21*d_A21_z + 2.0*A22*d_A22_z);
        D_AF2_x = (2.0*A2_00*2.0*A00*d_A00_x + 4.0*A2_01*(A01*d_A00_x + A00 * d_A01_x) + 4.0*A2_02*(A02*d_A00_x + A00 * d_A02_x)
            + 2.0*A2_11*2.0*A01*d_A01_x + 2.0*A2_22*2.0*A02*d_A02_x + 4.0*A2_12*(A01*d_A02_x + A02 * d_A01_x));
        D_AF2_y = (2.0*A2_00*2.0*A10*d_A10_y + 4.0*A2_01*(A11*d_A10_y + A10 * d_A11_y) + 4.0*A2_02*(A12*d_A10_y + A10 * d_A12_y)
            + 2.0*A2_11*2.0*A11*d_A11_y + 2.0*A2_22*2.0*A12*d_A12_y + 4.0*A2_12*(A11*d_A12_y + A12 * d_A11_y));
        D_AF2_z = (2.0*A2_00*2.0*A20*d_A20_z + 4.0*A2_01*(A21*d_A20_z + A20 * d_A21_z) + 4.0*A2_02*(A22*d_A20_z + A20 * d_A22_z)
            + 2.0*A2_11*2.0*A21*d_A21_z + 2.0*A2_22*2.0*A22*d_A22_z + 4.0*A2_12*(A21*d_A22_z + A22 * d_A21_z));
        D_AF_I_x = (AF*D_AF_x - 0.5*D_AF2_x)*I_AF_I_05; D_AF_I_y = (AF*D_AF_y - 0.5*D_AF2_y)*I_AF_I_05; D_AF_I_z = (AF*D_AF_z - 0.5*D_AF2_z)*I_AF_I_05;
        D_AF_x *= I_AF_05; D_AF_y *= I_AF_05; D_AF_z *= I_AF_05;
        D_det_A_x = d_A00_x * A11 * A22 + d_A01_x * A12 * A20 + d_A02_x * A10 * A21 - d_A02_x * A11 * A20 - d_A01_x * A10 * A22 - d_A00_x * A12 * A21;
        D_det_A_y = A00 * d_A11_y * A22 + A01 * d_A12_y * A20 + A02 * d_A10_y * A21 - A02 * d_A11_y * A20 - A01 * d_A10_y * A22 - A00 * d_A12_y * A21;
        D_det_A_z = A00 * A11 * d_A22_z + A01 * A12 * d_A20_z + A02 * A10 * d_A21_z - A02 * A11 * d_A20_z - A01 * A10 * d_A22_z - A00 * A12 * d_A21_z;
        inv_det_A_2_05 = 0.5 *(1.0 - i_det_A * i_det_A);
        dvex = D_det_A_x * inv_det_A_2_05;
        dvey = D_det_A_y * inv_det_A_2_05;
        dvez = D_det_A_z * inv_det_A_2_05;
        g = AF_05 * AF_I_05;
        dgx = D_AF_x * AF_I_05 + AF_05 * D_AF_I_x;
        dgy = D_AF_y * AF_I_05 + AF_05 * D_AF_I_y;
        dgz = D_AF_z * AF_I_05 + AF_05 * D_AF_I_z;
        dex = (dgx *i_det_A - g * D_det_A_x *i_det_A*i_det_A);
        dey = (dgy *i_det_A - g * D_det_A_y *i_det_A*i_det_A);
        dez = (dgz *i_det_A - g * D_det_A_z *i_det_A*i_det_A);
        e = g * i_det_A;
        mips_e = (e*e - 1.0)*0.125;
        volume_e = 0.5*(det_A + i_det_A);
        exp_k = (alpha*mips_e + beta * volume_e);
        if (exp_k > 60) exp_k = 60;
        exp_vec[i] = exp_k;
        d_mips_e_x = e * 0.25 * dex;
        d_mips_e_y = e * 0.25 * dey;
        d_mips_e_z = e * 0.25 * dez;
        gx_vec[i] = alpha * d_mips_e_x + beta * dvex;
        gy_vec[i] = alpha * d_mips_e_y + beta * dvey;
        gz_vec[i] = alpha * d_mips_e_z + beta * dvez;
    }
#if 0
    fmath::expd_v(exp_vec, vc_size);
    for (int i = 0; i < vc_size; ++i)
    {
        exp_e = exp_vec[i];
        local_energy += exp_e;
        gx += gx_vec[i] * exp_e;
        gy += gy_vec[i] * exp_e;
        gz += gz_vec[i] * exp_e;
    }
#else
    for (int i = 0; i < vc_size; ++i)
    {
        exp_e = std::exp(exp_vec[i]);
        local_energy += exp_e;
        gx += gx_vec[i] * exp_e;
        gy += gy_vec[i] * exp_e;
        gz += gz_vec[i] * exp_e;
    }
#endif
    double mu = 1.0; int bvf_size = 0; int bvv_size = 0;
    double* p_vf_pos_x = NULL; double* p_vf_pos_y = NULL; double* p_vf_pos_z = NULL;
    if (bv_flag)
    {
        double normal_e = 0.0; double smooth_e = 0;
        double ne_gx = 0.0; double ne_gy = 0.0; double ne_gz = 0.0;
        double se_gx = 0.0; double se_gy = 0.0; double se_gz = 0.0;
        p_vf_pos_x = vc_pos_x3_omp[omp_id].data();
        p_vf_pos_y = vc_pos_y3_omp[omp_id].data();
        p_vf_pos_z = vc_pos_z3_omp[omp_id].data();
        std::vector<int>& one_vv_id = bvv_id[v_id];
        std::vector<int>& one_vf_id = bvf_id[v_id];
        bvv_size = one_vv_id.size(); bvf_size = bvv_size;
        for (int i = 0; i < bvv_size; ++i)
        {
            int j = (i + 1) % bvv_size; int k = (i + 2) % bvv_size;
            x1 = p_dpx[one_vv_id[i]]; y1 = p_dpy[one_vv_id[i]]; z1 = p_dpz[one_vv_id[i]];
            x2 = p_dpx[one_vv_id[j]]; y2 = p_dpy[one_vv_id[j]]; z2 = p_dpz[one_vv_id[j]];
            x3 = p_dpx[one_vv_id[k]]; y3 = p_dpy[one_vv_id[k]]; z3 = p_dpz[one_vv_id[k]];
            p_vf_pos_x[4 * i + 0] = x1; p_vf_pos_x[4 * i + 1] = x2; p_vf_pos_x[4 * i + 2] = x3;
            p_vf_pos_y[4 * i + 0] = y1; p_vf_pos_y[4 * i + 1] = y2; p_vf_pos_y[4 * i + 2] = y3;
            p_vf_pos_z[4 * i + 0] = z1; p_vf_pos_z[4 * i + 1] = z2; p_vf_pos_z[4 * i + 2] = z3;
            OpenVolumeMesh::Geometry::Vec3d& tn = target_bfn[one_vf_id[i]];
            p_vf_pos_x[4 * i + 3] = tn[0]; p_vf_pos_y[4 * i + 3] = tn[1]; p_vf_pos_z[4 * i + 3] = tn[2];
            double nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
            double ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
            double nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);
            double len2 = nx * nx + ny * ny + nz * nz; double len = std::sqrt(len2);
            double E_ne = 0.5*((nx / len - tn[0])*(nx / len - tn[0]) + (ny / len - tn[1])*(ny / len - tn[1]) + (nz / len - tn[2])*(nz / len - tn[2]));
            if (update_all)
            {
                mu_vec[i] = local_energy_ratio[one_vf_id[i]];
            }
            else
            {
                mu_vec[i] = 1.0;
            }
            normal_e += mu_vec[i] * E_ne;
            double d_nx_x0 = 0.0; double d_nx_y0 = z1 - z2; double d_nx_z0 = y2 - y1;
            double d_ny_x0 = z2 - z1; double d_ny_y0 = 0.0; double d_ny_z0 = x1 - x2;
            double d_nz_x0 = y1 - y2; double d_nz_y0 = x2 - x1; double d_nz_z0 = 0.0;
            double d_len2_x0 = (nx*d_nx_x0 + ny * d_ny_x0 + nz * d_nz_x0);
            double d_len2_y0 = (nx*d_nx_y0 + ny * d_ny_y0 + nz * d_nz_y0);
            double d_len2_z0 = (nx*d_nx_z0 + ny * d_ny_z0 + nz * d_nz_z0);
            double d_len_x0 = d_len2_x0 / len; double d_len_y0 = d_len2_y0 / len; double d_len_z0 = d_len2_z0 / len;
            double d_nx_len_x0 = d_nx_x0 / len - nx * d_len_x0 / (len*len);
            double d_nx_len_y0 = d_nx_y0 / len - nx * d_len_y0 / (len*len);
            double d_nx_len_z0 = d_nx_z0 / len - nx * d_len_z0 / (len*len);
            double d_ny_len_x0 = d_ny_x0 / len - ny * d_len_x0 / (len*len);
            double d_ny_len_y0 = d_ny_y0 / len - ny * d_len_y0 / (len*len);
            double d_ny_len_z0 = d_ny_z0 / len - ny * d_len_z0 / (len*len);
            double d_nz_len_x0 = d_nz_x0 / len - nz * d_len_x0 / (len*len);
            double d_nz_len_y0 = d_nz_y0 / len - nz * d_len_y0 / (len*len);
            double d_nz_len_z0 = d_nz_z0 / len - nz * d_len_z0 / (len*len);
            ne_gx += mu_vec[i] * (d_nx_len_x0*(nx / len - tn[0]) + d_ny_len_x0 * (ny / len - tn[1]) + d_nz_len_x0 * (nz / len - tn[2]));
            ne_gy += mu_vec[i] * (d_nx_len_y0*(nx / len - tn[0]) + d_ny_len_y0 * (ny / len - tn[1]) + d_nz_len_y0 * (nz / len - tn[2]));
            ne_gz += mu_vec[i] * (d_nx_len_z0*(nx / len - tn[0]) + d_ny_len_z0 * (ny / len - tn[1]) + d_nz_len_z0 * (nz / len - tn[2]));
        }
        if (update_all)
        {
            local_energy += normal_e;
            gx += ne_gx; gy += ne_gy; gz += ne_gz;
        }
        else
        {
            double s = local_energy / (normal_e + 1e-8);
            {
                ga = s;
            }
            if (s > 1e16) ga = 1e16;
            for (int i = 0; i < bvf_size; ++i)
            {
                mu_vec[i] = ga;
            }
            local_energy += ga * normal_e;
            gx += ga * ne_gx; gy += ga * ne_gy; gz += ga * ne_gz;
        }
    }
    min_radius = std::sqrt(min_radius);
    double dx = gx, dy = gy, dz = gz;
    double move_d = sqrt(dx*dx + dy * dy + dz * dz);
    if (move_d > min_radius)
    {
        double t = min_radius / move_d;
        dx *= t; dy *= t; dz *= t;
        move_d = min_radius;
    }
    double npx = x0 - dx; double npy = y0 - dy; double npz = z0 - dz;
    while (!local_check_negative_volume4(vc_size, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, p_vc_n_cross_x, p_vc_n_cross_y, p_vc_n_cross_z, npx, npy, npz))
    {
        dx *= 0.8; dy *= 0.8; dz *= 0.8; move_d *= 0.8;
        if (move_d < 1e-8)
        {
            npx = x0; npy = y0; npz = z0;
            move_d = 1e-8; break;
        }
        else
        {
            npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
        }
    }
#if 1
    double new_e = 0;
    bool e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec);
    if (!e_ok)
    {
        while (!e_ok)
        {
            dx *= 0.2; dy *= 0.2; dz *= 0.2; move_d *= 0.2;
            if (move_d < 1e-8)
            {
                npx = x0; npy = y0; npz = z0;
                move_d = 1e-8; break;
            }
            else
            {
                npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
            }
            new_e = 0;
            e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec);
        }
        if (e_ok)
        {
            dx *= 5; dy *= 5; dz *= 5; move_d *= 5; e_ok = false;
            while (!e_ok)
            {
                dx *= 0.5; dy *= 0.5; dz *= 0.5; move_d *= 0.5;
                if (move_d < 1e-8)
                {
                    npx = x0; npy = y0; npz = z0;
                    move_d = 1e-8; break;
                }
                else
                {
                    npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
                }
                new_e = 0;
                e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec);
            }
            if (e_ok)
            {
                dx *= 2; dy *= 2; dz *= 2; move_d *= 2; e_ok = false;
                while (!e_ok)
                {
                    dx *= 0.875; dy *= 0.875; dz *= 0.875; move_d *= 0.875;
                    if (move_d < 1e-8)
                    {
                        npx = x0; npy = y0; npz = z0;
                        move_d = 1e-8; break;
                    }
                    else
                    {
                        npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
                    }
                    new_e = 0;
                    e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec);
                }
            }
        }
    }
#else
    double new_e = 0;
    bool e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
    while (!e_ok)
    {
        dx *= 0.8; dy *= 0.8; dz *= 0.8; move_d *= 0.8;
        if (move_d < 1e-8)
        {
            npx = x0; npy = y0; npz = z0;
            move_d = 1e-8; break;
        }
        else
        {
            npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
        }
        new_e = 0;
        e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
    }
#endif
    dpx[v_id] = npx; dpy[v_id] = npy; dpz[v_id] = npz;
    check_big_distortion_cell(vc_size, large_dis_flag_omp[omp_id], p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, npx, npy, npz);
    for (int i = 0; i < vc_size; ++i)
    {
        int c_id = vertex_cell[v_id][i];
        if (large_dis_flag_omp[omp_id][i] == 1)
        {
            if (distortion_big_cell[c_id] == 1)
            {
            }
            else
            {
                distortion_big_cell[c_id] = 1;
                for (int j = 0; j < 4; ++j)
                {
                    distortion_big_count_omp[omp_id][cell_vertex[c_id][j]] += 1;
                }
            }
        }
        else
        {
            if (distortion_big_cell[c_id] == 1)
            {
                distortion_big_cell[c_id] = -1;
                for (int j = 0; j < 4; ++j)
                {
                    distortion_big_count_omp[omp_id][cell_vertex[c_id][j]] -= 1;
                }
            }
            else
            {
            }
        }
    }
}
void polycube_deformation_interface::exp_mips_deformation_refine_one_polycube_normal_constrained_single_line(int v_id, double angle_area_ratio, double energy_power, const int& omp_id, const double move_dir[3], double ga, bool update_all)
{
    bool bv_flag = false;
    if (bvf_id[v_id].size() > 0)
    {
        bv_flag = true;
    }
    else
    {
        if (distortion_big_count[v_id] == 0)
            return;
    }
    double gx_all = 0.0; double gy_all = 0.0; double gz_all = 0.0;
    double min_radius = 1e30;
    const double* p_dpx = dpx.data(); const double* p_dpy = dpy.data(); const double* p_dpz = dpz.data();
    double x0 = p_dpx[v_id]; double y0 = p_dpy[v_id]; double z0 = p_dpz[v_id];
    int vc_size = vertex_cell[v_id].size();
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    double local_energy = 0.0; double gx = 0.0; double gy = 0.0; double gz = 0.0;
    double len;
    unsigned jj = 0; double alpha = (1.0 - angle_area_ratio)*energy_power; double beta = energy_power - alpha;
    const int* vc_id = vertex_cell[v_id].data();
    int v_id_1, v_id_2, v_id_3, k;
    double D00, D10, D20, D01, D11, D21, D02, D12, D22;
    double C00, C01, C02, C10, C11, C12, C20, C21, C22;
    double A00, A01, A02, A10, A11, A12, A20, A21, A22;
    double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
    double AF, AF_05, I_AF_05, AF2, AF_I, AF_I_05, I_AF_I_05, exp_k, exp_e;
    double D_AF_x, D_AF_y, D_AF_z, D_AF2_x, D_AF2_y, D_AF2_z, D_AF_I_x, D_AF_I_y, D_AF_I_z, D_det_A_x, D_det_A_y, D_det_A_z, inv_det_A_2_05;
    double d_A00_x, d_A10_y, d_A20_z, d_A01_x, d_A11_y, d_A21_z, d_A02_x, d_A12_y, d_A22_z;
    double dvex, dvey, dvez, g, dgx, dgy, dgz, dex, dey, dez, e, mips_e, volume_e, d_mips_e_x, d_mips_e_y, d_mips_e_z;
    double* p_vc_pos_x = vc_pos_x2_omp[omp_id].data();
    double* p_vc_pos_y = vc_pos_y2_omp[omp_id].data();
    double* p_vc_pos_z = vc_pos_z2_omp[omp_id].data();
    double* p_vc_n_cross_x = vc_n_cross_x_omp[omp_id].data();
    double* p_vc_n_cross_y = vc_n_cross_y_omp[omp_id].data();
    double* p_vc_n_cross_z = vc_n_cross_z_omp[omp_id].data();
    double* exp_vec = exp_vec_omp[omp_id].data();
    double* gx_vec = gx_vec_omp[omp_id].data();
    double* gy_vec = gy_vec_omp[omp_id].data();
    double* gz_vec = gz_vec_omp[omp_id].data();
    double* mu_vec = mu_vec_omp[omp_id].data();
    std::vector<std::vector<double>>& vc_S = vc_S_omp[omp_id];
    for (unsigned i = 0; i < vc_size; ++i)
    {
        const int* vv_id = vertex_cell_vertex[v_id][i].data();
        memcpy(&vc_S[i][0], &vcv_S[v_id][i][0], 9 * sizeof(double));
        double* s_data = vc_S[i].data();
        k = 3 * i;
        v_id_1 = vv_id[0];
        v_id_2 = vv_id[1];
        v_id_3 = vv_id[2];
        x1 = p_dpx[v_id_1]; x2 = p_dpx[v_id_2]; x3 = p_dpx[v_id_3];
        y1 = p_dpy[v_id_1]; y2 = p_dpy[v_id_2]; y3 = p_dpy[v_id_3];
        z1 = p_dpz[v_id_1]; z2 = p_dpz[v_id_2]; z3 = p_dpz[v_id_3];
        p_vc_pos_x[k + 0] = x1; p_vc_pos_x[k + 1] = x2; p_vc_pos_x[k + 2] = x3;
        p_vc_pos_y[k + 0] = y1; p_vc_pos_y[k + 1] = y2; p_vc_pos_y[k + 2] = y3;
        p_vc_pos_z[k + 0] = z1; p_vc_pos_z[k + 1] = z2; p_vc_pos_z[k + 2] = z3;
        p_vc_n_cross_x[i] = (y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1);
        p_vc_n_cross_y[i] = (x3 - x1)*(z2 - z1) - (x2 - x1)*(z3 - z1);
        p_vc_n_cross_z[i] = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
        D00 = x1 - x0; D10 = y1 - y0; D20 = z1 - z0;
        D01 = x2 - x0; D11 = y2 - y0; D21 = z2 - z0;
        D02 = x3 - x0; D12 = y3 - y0; D22 = z3 - z0;
        len = D00 * D00 + D10 * D10 + D20 * D20; if (len < min_radius) min_radius = len;
        len = D01 * D01 + D11 * D11 + D21 * D21; if (len < min_radius) min_radius = len;
        len = D02 * D02 + D12 * D12 + D22 * D22; if (len < min_radius) min_radius = len;
        C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
        C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
        C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];
        A00 = C00 * D00 + C10 * D01 + C20 * D02; d_A00_x = -C00 - C10 - C20;
        A10 = C00 * D10 + C10 * D11 + C20 * D12; d_A10_y = -C00 - C10 - C20;
        A20 = C00 * D20 + C10 * D21 + C20 * D22; d_A20_z = -C00 - C10 - C20;
        A01 = C01 * D00 + C11 * D01 + C21 * D02; d_A01_x = -C01 - C11 - C21;
        A11 = C01 * D10 + C11 * D11 + C21 * D12; d_A11_y = -C01 - C11 - C21;
        A21 = C01 * D20 + C11 * D21 + C21 * D22; d_A21_z = -C01 - C11 - C21;
        A02 = C02 * D00 + C12 * D01 + C22 * D02; d_A02_x = -C02 - C12 - C22;
        A12 = C02 * D10 + C12 * D11 + C22 * D12; d_A12_y = -C02 - C12 - C22;
        A22 = C02 * D20 + C12 * D21 + C22 * D22; d_A22_z = -C02 - C12 - C22;
        A2_00 = A00 * A00 + A10 * A10 + A20 * A20;
        A2_01 = A00 * A01 + A10 * A11 + A20 * A21;
        A2_02 = A00 * A02 + A10 * A12 + A20 * A22;
        A2_11 = A01 * A01 + A11 * A11 + A21 * A21;
        A2_12 = A01 * A02 + A11 * A12 + A21 * A22;
        A2_22 = A02 * A02 + A12 * A12 + A22 * A22;
        AF = A2_00 + A2_11 + A2_22; AF_05 = std::sqrt(AF); I_AF_05 = 0.5 / AF_05;
        AF2 = A2_00 * A2_00 + A2_01 * A2_01 * 2 + A2_02 * A2_02 * 2 + A2_11 * A2_11 + A2_12 * A2_12 * 2 + A2_22 * A2_22;
        AF_I = (AF*AF - AF2)*0.5; AF_I_05 = std::sqrt(AF_I); I_AF_I_05 = 0.5 / AF_I_05;
        det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
        i_det_A = 1.0 / det_A;
        D_AF_x = (2.0*A00*d_A00_x + 2.0*A01*d_A01_x + 2.0*A02*d_A02_x);
        D_AF_y = (2.0*A10*d_A10_y + 2.0*A11*d_A11_y + 2.0*A12*d_A12_y);
        D_AF_z = (2.0*A20*d_A20_z + 2.0*A21*d_A21_z + 2.0*A22*d_A22_z);
        D_AF2_x = (2.0*A2_00*2.0*A00*d_A00_x + 4.0*A2_01*(A01*d_A00_x + A00 * d_A01_x) + 4.0*A2_02*(A02*d_A00_x + A00 * d_A02_x)
            + 2.0*A2_11*2.0*A01*d_A01_x + 2.0*A2_22*2.0*A02*d_A02_x + 4.0*A2_12*(A01*d_A02_x + A02 * d_A01_x));
        D_AF2_y = (2.0*A2_00*2.0*A10*d_A10_y + 4.0*A2_01*(A11*d_A10_y + A10 * d_A11_y) + 4.0*A2_02*(A12*d_A10_y + A10 * d_A12_y)
            + 2.0*A2_11*2.0*A11*d_A11_y + 2.0*A2_22*2.0*A12*d_A12_y + 4.0*A2_12*(A11*d_A12_y + A12 * d_A11_y));
        D_AF2_z = (2.0*A2_00*2.0*A20*d_A20_z + 4.0*A2_01*(A21*d_A20_z + A20 * d_A21_z) + 4.0*A2_02*(A22*d_A20_z + A20 * d_A22_z)
            + 2.0*A2_11*2.0*A21*d_A21_z + 2.0*A2_22*2.0*A22*d_A22_z + 4.0*A2_12*(A21*d_A22_z + A22 * d_A21_z));
        D_AF_I_x = (AF*D_AF_x - 0.5*D_AF2_x)*I_AF_I_05; D_AF_I_y = (AF*D_AF_y - 0.5*D_AF2_y)*I_AF_I_05; D_AF_I_z = (AF*D_AF_z - 0.5*D_AF2_z)*I_AF_I_05;
        D_AF_x *= I_AF_05; D_AF_y *= I_AF_05; D_AF_z *= I_AF_05;
        D_det_A_x = d_A00_x * A11 * A22 + d_A01_x * A12 * A20 + d_A02_x * A10 * A21 - d_A02_x * A11 * A20 - d_A01_x * A10 * A22 - d_A00_x * A12 * A21;
        D_det_A_y = A00 * d_A11_y * A22 + A01 * d_A12_y * A20 + A02 * d_A10_y * A21 - A02 * d_A11_y * A20 - A01 * d_A10_y * A22 - A00 * d_A12_y * A21;
        D_det_A_z = A00 * A11 * d_A22_z + A01 * A12 * d_A20_z + A02 * A10 * d_A21_z - A02 * A11 * d_A20_z - A01 * A10 * d_A22_z - A00 * A12 * d_A21_z;
        inv_det_A_2_05 = 0.5 *(1.0 - i_det_A * i_det_A);
        dvex = D_det_A_x * inv_det_A_2_05;
        dvey = D_det_A_y * inv_det_A_2_05;
        dvez = D_det_A_z * inv_det_A_2_05;
        g = AF_05 * AF_I_05;
        dgx = D_AF_x * AF_I_05 + AF_05 * D_AF_I_x;
        dgy = D_AF_y * AF_I_05 + AF_05 * D_AF_I_y;
        dgz = D_AF_z * AF_I_05 + AF_05 * D_AF_I_z;
        dex = (dgx *i_det_A - g * D_det_A_x *i_det_A*i_det_A);
        dey = (dgy *i_det_A - g * D_det_A_y *i_det_A*i_det_A);
        dez = (dgz *i_det_A - g * D_det_A_z *i_det_A*i_det_A);
        e = g * i_det_A;
        mips_e = (e*e - 1.0)*0.125;
        volume_e = 0.5*(det_A + i_det_A);
        exp_k = (alpha*mips_e + beta * volume_e);
        if (exp_k > 60) exp_k = 60;
        exp_vec[i] = exp_k;
        d_mips_e_x = e * 0.25 * dex;
        d_mips_e_y = e * 0.25 * dey;
        d_mips_e_z = e * 0.25 * dez;
        gx_vec[i] = alpha * d_mips_e_x + beta * dvex;
        gy_vec[i] = alpha * d_mips_e_y + beta * dvey;
        gz_vec[i] = alpha * d_mips_e_z + beta * dvez;
    }
#if 0
    fmath::expd_v(exp_vec, vc_size);
    for (int i = 0; i < vc_size; ++i)
    {
        exp_e = exp_vec[i];
        local_energy += exp_e;
        gx += gx_vec[i] * exp_e;
        gy += gy_vec[i] * exp_e;
        gz += gz_vec[i] * exp_e;
    }
#else
    for (int i = 0; i < vc_size; ++i)
    {
        exp_e = std::exp(exp_vec[i]);
        local_energy += exp_e;
        gx += gx_vec[i] * exp_e;
        gy += gy_vec[i] * exp_e;
        gz += gz_vec[i] * exp_e;
    }
#endif
    double mu = 1.0; int bvf_size = 0; int bvv_size = 0;
    double* p_vf_pos_x = NULL; double* p_vf_pos_y = NULL; double* p_vf_pos_z = NULL;
    if (bv_flag)
    {
        double normal_e = 0.0; double smooth_e = 0;
        double ne_gx = 0.0; double ne_gy = 0.0; double ne_gz = 0.0;
        double se_gx = 0.0; double se_gy = 0.0; double se_gz = 0.0;
        p_vf_pos_x = vc_pos_x3_omp[omp_id].data();
        p_vf_pos_y = vc_pos_y3_omp[omp_id].data();
        p_vf_pos_z = vc_pos_z3_omp[omp_id].data();
        std::vector<int>& one_vv_id = bvv_id[v_id];
        std::vector<int>& one_vf_id = bvf_id[v_id];
        bvv_size = one_vv_id.size(); bvf_size = bvv_size;
        for (int i = 0; i < bvv_size; ++i)
        {
            int j = (i + 1) % bvv_size; int k = (i + 2) % bvv_size;
            x1 = p_dpx[one_vv_id[i]]; y1 = p_dpy[one_vv_id[i]]; z1 = p_dpz[one_vv_id[i]];
            x2 = p_dpx[one_vv_id[j]]; y2 = p_dpy[one_vv_id[j]]; z2 = p_dpz[one_vv_id[j]];
            x3 = p_dpx[one_vv_id[k]]; y3 = p_dpy[one_vv_id[k]]; z3 = p_dpz[one_vv_id[k]];
            p_vf_pos_x[4 * i + 0] = x1; p_vf_pos_x[4 * i + 1] = x2; p_vf_pos_x[4 * i + 2] = x3;
            p_vf_pos_y[4 * i + 0] = y1; p_vf_pos_y[4 * i + 1] = y2; p_vf_pos_y[4 * i + 2] = y3;
            p_vf_pos_z[4 * i + 0] = z1; p_vf_pos_z[4 * i + 1] = z2; p_vf_pos_z[4 * i + 2] = z3;
            OpenVolumeMesh::Geometry::Vec3d& tn = target_bfn[one_vf_id[i]];
            p_vf_pos_x[4 * i + 3] = tn[0]; p_vf_pos_y[4 * i + 3] = tn[1]; p_vf_pos_z[4 * i + 3] = tn[2];
            double nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
            double ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
            double nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);
            double len2 = nx * nx + ny * ny + nz * nz; double len = std::sqrt(len2);
            double E_ne = 0.5*((nx / len - tn[0])*(nx / len - tn[0]) + (ny / len - tn[1])*(ny / len - tn[1]) + (nz / len - tn[2])*(nz / len - tn[2]));
            if (update_all)
            {
                mu_vec[i] = local_energy_ratio[one_vf_id[i]];
            }
            else
            {
                mu_vec[i] = 1.0;
            }
            normal_e += mu_vec[i] * E_ne;
            double d_nx_x0 = 0.0; double d_nx_y0 = z1 - z2; double d_nx_z0 = y2 - y1;
            double d_ny_x0 = z2 - z1; double d_ny_y0 = 0.0; double d_ny_z0 = x1 - x2;
            double d_nz_x0 = y1 - y2; double d_nz_y0 = x2 - x1; double d_nz_z0 = 0.0;
            double d_len2_x0 = (nx*d_nx_x0 + ny * d_ny_x0 + nz * d_nz_x0);
            double d_len2_y0 = (nx*d_nx_y0 + ny * d_ny_y0 + nz * d_nz_y0);
            double d_len2_z0 = (nx*d_nx_z0 + ny * d_ny_z0 + nz * d_nz_z0);
            double d_len_x0 = d_len2_x0 / len; double d_len_y0 = d_len2_y0 / len; double d_len_z0 = d_len2_z0 / len;
            double d_nx_len_x0 = d_nx_x0 / len - nx * d_len_x0 / (len*len);
            double d_nx_len_y0 = d_nx_y0 / len - nx * d_len_y0 / (len*len);
            double d_nx_len_z0 = d_nx_z0 / len - nx * d_len_z0 / (len*len);
            double d_ny_len_x0 = d_ny_x0 / len - ny * d_len_x0 / (len*len);
            double d_ny_len_y0 = d_ny_y0 / len - ny * d_len_y0 / (len*len);
            double d_ny_len_z0 = d_ny_z0 / len - ny * d_len_z0 / (len*len);
            double d_nz_len_x0 = d_nz_x0 / len - nz * d_len_x0 / (len*len);
            double d_nz_len_y0 = d_nz_y0 / len - nz * d_len_y0 / (len*len);
            double d_nz_len_z0 = d_nz_z0 / len - nz * d_len_z0 / (len*len);
            ne_gx += mu_vec[i] * (d_nx_len_x0*(nx / len - tn[0]) + d_ny_len_x0 * (ny / len - tn[1]) + d_nz_len_x0 * (nz / len - tn[2]));
            ne_gy += mu_vec[i] * (d_nx_len_y0*(nx / len - tn[0]) + d_ny_len_y0 * (ny / len - tn[1]) + d_nz_len_y0 * (nz / len - tn[2]));
            ne_gz += mu_vec[i] * (d_nx_len_z0*(nx / len - tn[0]) + d_ny_len_z0 * (ny / len - tn[1]) + d_nz_len_z0 * (nz / len - tn[2]));
        }
        if (update_all)
        {
            local_energy += normal_e;
            gx += ne_gx; gy += ne_gy; gz += ne_gz;
        }
        else
        {
            double s = local_energy / (normal_e + 1e-8);
            {
                ga = s;
            }
            if (s > 1e16) ga = 1e16;
            for (int i = 0; i < bvf_size; ++i)
            {
                mu_vec[i] = ga;
            }
            local_energy += ga * normal_e;
            gx += ga * ne_gx; gy += ga * ne_gy; gz += ga * ne_gz;
        }
    }
    min_radius = std::sqrt(min_radius);
    double dx = gx, dy = gy, dz = gz;
    dx = dx * move_dir[0];
    dy = dy * move_dir[1];
    dz = dz * move_dir[2];
    double move_d = sqrt(dx*dx + dy * dy + dz * dz);
    if (move_d > min_radius)
    {
        double t = min_radius / move_d;
        dx *= t; dy *= t; dz *= t;
        move_d = min_radius;
    }
    double npx = x0 - dx; double npy = y0 - dy; double npz = z0 - dz;
    while (!local_check_negative_volume4(vc_size, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, p_vc_n_cross_x, p_vc_n_cross_y, p_vc_n_cross_z, npx, npy, npz))
    {
        dx *= 0.8; dy *= 0.8; dz *= 0.8; move_d *= 0.8;
        if (move_d < 1e-8)
        {
            npx = x0; npy = y0; npz = z0;
            move_d = 1e-8; break;
        }
        else
        {
            npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
        }
    }
#if 1
    double new_e = 0;
    bool e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec);
    if (!e_ok)
    {
        while (!e_ok)
        {
            dx *= 0.2; dy *= 0.2; dz *= 0.2; move_d *= 0.2;
            if (move_d < 1e-8)
            {
                npx = x0; npy = y0; npz = z0;
                move_d = 1e-8; break;
            }
            else
            {
                npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
            }
            new_e = 0;
            e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec);
        }
        if (e_ok)
        {
            dx *= 5; dy *= 5; dz *= 5; move_d *= 5; e_ok = false;
            while (!e_ok)
            {
                dx *= 0.5; dy *= 0.5; dz *= 0.5; move_d *= 0.5;
                if (move_d < 1e-8)
                {
                    npx = x0; npy = y0; npz = z0;
                    move_d = 1e-8; break;
                }
                else
                {
                    npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
                }
                new_e = 0;
                e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec);
            }
            if (e_ok)
            {
                dx *= 2; dy *= 2; dz *= 2; move_d *= 2; e_ok = false;
                while (!e_ok)
                {
                    dx *= 0.875; dy *= 0.875; dz *= 0.875; move_d *= 0.875;
                    if (move_d < 1e-8)
                    {
                        npx = x0; npy = y0; npz = z0;
                        move_d = 1e-8; break;
                    }
                    else
                    {
                        npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
                    }
                    new_e = 0;
                    e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec);
                }
            }
        }
    }
#else
    double new_e = 0;
    bool e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
    while (!e_ok)
    {
        dx *= 0.8; dy *= 0.8; dz *= 0.8; move_d *= 0.8;
        if (move_d < 1e-8)
        {
            npx = x0; npy = y0; npz = z0;
            move_d = 1e-8; break;
        }
        else
        {
            npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
        }
        new_e = 0;
        e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
    }
#endif
    dpx[v_id] = npx; dpy[v_id] = npy; dpz[v_id] = npz;
    check_big_distortion_cell(vc_size, large_dis_flag_omp[omp_id], p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, npx, npy, npz);
    for (int i = 0; i < vc_size; ++i)
    {
        int c_id = vertex_cell[v_id][i];
        if (large_dis_flag_omp[omp_id][i] == 1)
        {
            if (distortion_big_cell[c_id] == 1)
            {
            }
            else
            {
                distortion_big_cell[c_id] = 1;
                for (int j = 0; j < 4; ++j)
                {
                    distortion_big_count_omp[omp_id][cell_vertex[c_id][j]] += 1;
                }
            }
        }
        else
        {
            if (distortion_big_cell[c_id] == 1)
            {
                distortion_big_cell[c_id] = -1;
                for (int j = 0; j < 4; ++j)
                {
                    distortion_big_count_omp[omp_id][cell_vertex[c_id][j]] -= 1;
                }
            }
            else
            {
            }
        }
    }
}
const std::array<Eigen::Vector3d, 8> polycube_deformation_interface::a_three_cut_opposite_dir{
    {
        { -1, -1, -1 },
        { -1, -1, 1 },
        { -1, 1, -1 },
        { -1, 1, 1 },
        { 1, -1, -1 },
        { 1, -1, 1 },
        { -1, -1, 1 },
        { 1, 1, 1 }
    }
};
void polycube_deformation_interface::exp_mips_deformation_refine_one_polycube_normal_constrained_three_cut(
    const int &v_id,
    const double &angle_area_ratio,
    const double &energy_power,
    const int &omp_id,
    const std::array<int, 3> &adjacent_one_cut_type,
    const int &common_vert,
    const double &ga,
    const bool update_all
)
{
    static const std::array<Eigen::Vector3d, 12> s_a_cut_idx2normal{
        {
           { 1, -1, 0 },
           { 1, 1, 0 },
           { 1, 1, 0 },
           { 1, -1, 0 },
           {-1, 0, 1 },
           { 1, 0, 1 },
           { 1, 0, 1 },
           {-1, 0, 1 },
           { 0, -1, 1 },
           { 0, 1, 1 },
           { 0, 1, 1 },
           { 0, -1, 1 }
        }
    };
    bool bv_flag = false;
    if (bvf_id[v_id].size() > 0)
    {
        bv_flag = true;
    }
    double local_energy_all = 0.0;
    const double* p_dpx = dpx.data(); const double* p_dpy = dpy.data(); const double* p_dpz = dpz.data();
    std::array<Eigen::Vector3d, 3> a_new_coord{
        {
            { dpx[v_id], dpy[v_id], dpz[v_id] }
        }
    };
    std::array<int, 3> a_equal_vertex_id{
        v_id,
        three_cut_vert_to_pair[v_id].first,
        three_cut_vert_to_pair[v_id].second
    };
    std::array<double, 3> a_min_radius{
        1e30, 1e30, 1e30
    };
    std::array<double, 3> a_vc_size{
        vertex_cell[a_equal_vertex_id[0]].size(),
        vertex_cell[a_equal_vertex_id[1]].size(),
        vertex_cell[a_equal_vertex_id[2]].size()
    };
    std::vector<double> vc_pos_x2_omp_1 = vc_pos_x2_omp[omp_id];
    std::vector<double> vc_pos_y2_omp_1 = vc_pos_y2_omp[omp_id];
    std::vector<double> vc_pos_z2_omp_1 = vc_pos_z2_omp[omp_id];
    std::vector<double> vc_pos_x2_omp_2 = vc_pos_x2_omp[omp_id];
    std::vector<double> vc_pos_y2_omp_2 = vc_pos_y2_omp[omp_id];
    std::vector<double> vc_pos_z2_omp_2 = vc_pos_z2_omp[omp_id];
    std::array<double *, 3> a_p_vc_pos_x{
        vc_pos_x2_omp[omp_id].data(),
        vc_pos_x2_omp_1.data(),
        vc_pos_x2_omp_2.data()
    };
    std::array<double *, 3> a_p_vc_pos_y{
        vc_pos_y2_omp[omp_id].data(),
        vc_pos_y2_omp_1.data(),
        vc_pos_y2_omp_2.data()
    };
    std::array<double *, 3> a_p_vc_pos_z{
        vc_pos_z2_omp[omp_id].data(),
        vc_pos_z2_omp_1.data(),
        vc_pos_z2_omp_2.data()
    };
    std::vector<double> vc_n_cross_x_omp_1 = vc_n_cross_x_omp[omp_id];
    std::vector<double> vc_n_cross_y_omp_1 = vc_n_cross_y_omp[omp_id];
    std::vector<double> vc_n_cross_z_omp_1 = vc_n_cross_z_omp[omp_id];
    std::vector<double> vc_n_cross_x_omp_2 = vc_n_cross_x_omp[omp_id];
    std::vector<double> vc_n_cross_y_omp_2 = vc_n_cross_y_omp[omp_id];
    std::vector<double> vc_n_cross_z_omp_2 = vc_n_cross_z_omp[omp_id];
    std::array<double *, 3> a_p_vc_n_cross_x{
        vc_n_cross_x_omp[omp_id].data(),
        vc_n_cross_x_omp_1.data(),
        vc_n_cross_x_omp_2.data()
    };
    std::array<double *, 3> a_p_vc_n_cross_y{
        vc_n_cross_y_omp[omp_id].data(),
        vc_n_cross_y_omp_1.data(),
        vc_n_cross_y_omp_2.data()
    };
    std::array<double *, 3> a_p_vc_n_cross_z{
        vc_n_cross_z_omp[omp_id].data(),
        vc_n_cross_z_omp_1.data(),
        vc_n_cross_z_omp_2.data()
    };
    std::vector<double> exp_vec_omp_1 = exp_vec_omp[omp_id];
    std::vector<double> exp_vec_omp_2 = exp_vec_omp[omp_id];
    std::array<double *, 3> a_exp_vec{
        exp_vec_omp[omp_id].data(),
        exp_vec_omp_1.data(),
        exp_vec_omp_2.data()
    };
    std::vector<double> gx_vec_omp_1 = gx_vec_omp[omp_id];
    std::vector<double> gy_vec_omp_1 = gy_vec_omp[omp_id];
    std::vector<double> gz_vec_omp_1 = gz_vec_omp[omp_id];
    std::vector<double> mu_vec_omp_1 = mu_vec_omp[omp_id];
    std::vector<double> gx_vec_omp_2 = gx_vec_omp[omp_id];
    std::vector<double> gy_vec_omp_2 = gy_vec_omp[omp_id];
    std::vector<double> gz_vec_omp_2 = gz_vec_omp[omp_id];
    std::vector<double> mu_vec_omp_2 = mu_vec_omp[omp_id];
    std::array<double *, 3> a_gx_vec{
        gx_vec_omp[omp_id].data(),
        gx_vec_omp_1.data(),
        gx_vec_omp_2.data()
    };
    std::array<double *, 3> a_gy_vec{
        gy_vec_omp[omp_id].data(),
        gy_vec_omp_1.data(),
        gy_vec_omp_2.data()
    };
    std::array<double *, 3> a_gz_vec{
        gz_vec_omp[omp_id].data(),
        gz_vec_omp_1.data(),
        gz_vec_omp_2.data()
    };
    std::array<double *, 3> a_mu_vec{
        mu_vec_omp[omp_id].data(),
        mu_vec_omp_1.data(),
        mu_vec_omp_2.data()
    };
    std::vector<double> vc_pos_x3_omp_1 = vc_pos_x3_omp[omp_id];
    std::vector<double> vc_pos_y3_omp_1 = vc_pos_y3_omp[omp_id];
    std::vector<double> vc_pos_z3_omp_1 = vc_pos_z3_omp[omp_id];
    std::vector<double> vc_pos_x3_omp_2 = vc_pos_x3_omp[omp_id];
    std::vector<double> vc_pos_y3_omp_2 = vc_pos_y3_omp[omp_id];
    std::vector<double> vc_pos_z3_omp_2 = vc_pos_z3_omp[omp_id];
    std::array<double *, 3> a_p_vf_pos_x{
        vc_pos_x3_omp[omp_id].data(),
        vc_pos_x3_omp_1.data(),
        vc_pos_x3_omp_2.data()
    };
    std::array<double *, 3> a_p_vf_pos_y{
        vc_pos_y3_omp[omp_id].data(),
        vc_pos_y3_omp_1.data(),
        vc_pos_y3_omp_2.data()
    };
    std::array<double *, 3> a_p_vf_pos_z{
        vc_pos_z3_omp[omp_id].data(),
        vc_pos_z3_omp_1.data(),
        vc_pos_z3_omp_2.data()
    };
    std::array<std::vector<std::vector<double>>, 3> avv_vc_S{ vc_S_omp[omp_id], vc_S_omp[omp_id], vc_S_omp[omp_id] };
    std::array<int, 3> a_bvf_size{ 0, 0, 0 };
    double alpha = (1.0 - angle_area_ratio) * energy_power;
    double beta = energy_power - alpha;
    double normal_e_all = 0.0;
    double gx_new, gy_new, gz_new;
    compute_vertex_update_info(
        update_all, bv_flag, a_equal_vertex_id[0],
        p_dpx, p_dpy, p_dpz,
        alpha, beta, avv_vc_S[0],
        a_p_vf_pos_x[0], a_p_vf_pos_y[0], a_p_vf_pos_z[0],
        a_p_vc_pos_x[0], a_p_vc_pos_y[0], a_p_vc_pos_z[0],
        a_p_vc_n_cross_x[0], a_p_vc_n_cross_y[0], a_p_vc_n_cross_z[0],
        a_exp_vec[0], a_mu_vec[0], a_gx_vec[0], a_gy_vec[0], a_gz_vec[0],
        local_energy_all, normal_e_all, gx_new, gy_new, gz_new, a_min_radius[0], a_bvf_size[0]
    );
    Eigen::Vector3d g_all(gx_new, gy_new, gz_new);
    compute_vertex_update_info(
        update_all, bv_flag, a_equal_vertex_id[1],
        p_dpx, p_dpy, p_dpz,
        alpha, beta, avv_vc_S[1],
        a_p_vf_pos_x[1], a_p_vf_pos_y[1], a_p_vf_pos_z[1],
        a_p_vc_pos_x[1], a_p_vc_pos_y[1], a_p_vc_pos_z[1],
        a_p_vc_n_cross_x[1], a_p_vc_n_cross_y[1], a_p_vc_n_cross_z[1],
        a_exp_vec[1], a_mu_vec[1], a_gx_vec[1], a_gy_vec[1], a_gz_vec[1],
        local_energy_all, normal_e_all, gx_new, gy_new, gz_new, a_min_radius[1], a_bvf_size[1]
    );
    g_all[0] += type_matrix[adjacent_one_cut_type[0]][0][0] * gx_new + type_matrix[adjacent_one_cut_type[0]][0][1] * gy_new + type_matrix[adjacent_one_cut_type[0]][0][2] * gz_new;
    g_all[1] += type_matrix[adjacent_one_cut_type[0]][1][0] * gx_new + type_matrix[adjacent_one_cut_type[0]][1][1] * gy_new + type_matrix[adjacent_one_cut_type[0]][1][2] * gz_new;
    g_all[2] += type_matrix[adjacent_one_cut_type[0]][2][0] * gx_new + type_matrix[adjacent_one_cut_type[0]][2][1] * gy_new + type_matrix[adjacent_one_cut_type[0]][2][2] * gz_new;
    compute_vertex_update_info(
        update_all, bv_flag, a_equal_vertex_id[2],
        p_dpx, p_dpy, p_dpz,
        alpha, beta, avv_vc_S[2],
        a_p_vf_pos_x[2], a_p_vf_pos_y[2], a_p_vf_pos_z[2],
        a_p_vc_pos_x[2], a_p_vc_pos_y[2], a_p_vc_pos_z[2],
        a_p_vc_n_cross_x[2], a_p_vc_n_cross_y[2], a_p_vc_n_cross_z[2],
        a_exp_vec[2], a_mu_vec[2], a_gx_vec[2], a_gy_vec[2], a_gz_vec[2],
        local_energy_all, normal_e_all, gx_new, gy_new, gz_new, a_min_radius[2], a_bvf_size[2]
    );
    g_all[0] += type_matrix[adjacent_one_cut_type[1]][0][0] * gx_new + type_matrix[adjacent_one_cut_type[1]][0][1] * gy_new + type_matrix[adjacent_one_cut_type[1]][0][2] * gz_new;
    g_all[1] += type_matrix[adjacent_one_cut_type[1]][1][0] * gx_new + type_matrix[adjacent_one_cut_type[1]][1][1] * gy_new + type_matrix[adjacent_one_cut_type[1]][1][2] * gz_new;
    g_all[2] += type_matrix[adjacent_one_cut_type[1]][2][0] * gx_new + type_matrix[adjacent_one_cut_type[1]][2][1] * gy_new + type_matrix[adjacent_one_cut_type[1]][2][2] * gz_new;
    double min_radius = std::min(std::min(a_min_radius[0], a_min_radius[1]), a_min_radius[2]);
    double move_d = g_all.norm();
    if (move_d > min_radius)
    {
        double t = min_radius / move_d;
        g_all *= t;
        move_d = min_radius;
    }
    Eigen::Vector3d common_vertex_coord(dpx[common_vert], dpy[common_vert], dpz[common_vert]);
    g_all = project2plane(s_a_cut_idx2normal[adjacent_one_cut_type[2]], g_all);
    a_new_coord[0] -= g_all;
    get_pair_coordinate(a_new_coord[0][0], a_new_coord[0][1], a_new_coord[0][2], type_matrix[adjacent_one_cut_type[0]], common_vertex_coord[0], common_vertex_coord[1], common_vertex_coord[2], a_new_coord[1][0], a_new_coord[1][1], a_new_coord[1][2]);
    get_pair_coordinate(a_new_coord[0][0], a_new_coord[0][1], a_new_coord[0][2], type_matrix[adjacent_one_cut_type[1]], common_vertex_coord[0], common_vertex_coord[1], common_vertex_coord[2], a_new_coord[2][0], a_new_coord[2][1], a_new_coord[2][2]);
    while (!local_check_negative_volume4(
        a_vc_size[0], a_p_vc_pos_x[0], a_p_vc_pos_y[0], a_p_vc_pos_z[0],
        a_p_vc_n_cross_x[0], a_p_vc_n_cross_y[0], a_p_vc_n_cross_z[0],
        a_new_coord[0].x(), a_new_coord[0].y(), a_new_coord[0].z()))
    {
        g_all *= 0.8; move_d *= 0.8;
        a_new_coord[0][0] = dpx[v_id];
        a_new_coord[0][1] = dpy[v_id];
        a_new_coord[0][2] = dpz[v_id];
        if (move_d < 1e-8)
        {
            move_d = 1e-8; break;
        }
        else
        {
            a_new_coord[0] -= g_all;
        }
    }
    while (!local_check_negative_volume4(
        a_vc_size[1], a_p_vc_pos_x[1], a_p_vc_pos_y[1], a_p_vc_pos_z[1],
        a_p_vc_n_cross_x[1], a_p_vc_n_cross_y[1], a_p_vc_n_cross_z[1],
        a_new_coord[1].x(), a_new_coord[1].y(), a_new_coord[1].z()))
    {
        g_all *= 0.8; move_d *= 0.8;
        a_new_coord[0][0] = dpx[v_id];
        a_new_coord[0][1] = dpy[v_id];
        a_new_coord[0][2] = dpz[v_id];
        if (move_d < 1e-8)
        {
            get_pair_coordinate(a_new_coord[0][0], a_new_coord[0][1], a_new_coord[0][2], type_matrix[adjacent_one_cut_type[0]], common_vertex_coord[0], common_vertex_coord[1], common_vertex_coord[2], a_new_coord[1][0], a_new_coord[1][1], a_new_coord[1][2]);
            move_d = 1e-8; break;
        }
        else
        {
            a_new_coord[0] -= g_all;
            get_pair_coordinate(a_new_coord[0][0], a_new_coord[0][1], a_new_coord[0][2], type_matrix[adjacent_one_cut_type[0]], common_vertex_coord[0], common_vertex_coord[1], common_vertex_coord[2], a_new_coord[1][0], a_new_coord[1][1], a_new_coord[1][2]);
        }
    }
    while (!local_check_negative_volume4(
        a_vc_size[2], a_p_vc_pos_x[2], a_p_vc_pos_y[2], a_p_vc_pos_z[2],
        a_p_vc_n_cross_x[2], a_p_vc_n_cross_y[2], a_p_vc_n_cross_z[2],
        a_new_coord[2].x(), a_new_coord[2].y(), a_new_coord[2].z()))
    {
        g_all *= 0.8; move_d *= 0.8;
        a_new_coord[0][0] = dpx[v_id];
        a_new_coord[0][1] = dpy[v_id];
        a_new_coord[0][2] = dpz[v_id];
        if (move_d < 1e-8)
        {
            get_pair_coordinate(a_new_coord[0][0], a_new_coord[0][1], a_new_coord[0][2], type_matrix[adjacent_one_cut_type[1]], common_vertex_coord[0], common_vertex_coord[1], common_vertex_coord[2], a_new_coord[2][0], a_new_coord[2][1], a_new_coord[2][2]);
            move_d = 1e-8; break;
        }
        else
        {
            a_new_coord[0] -= g_all;
            get_pair_coordinate(a_new_coord[0][0], a_new_coord[0][1], a_new_coord[0][2], type_matrix[adjacent_one_cut_type[1]], common_vertex_coord[0], common_vertex_coord[1], common_vertex_coord[2], a_new_coord[2][0], a_new_coord[2][1], a_new_coord[2][2]);
        }
    }
    std::array<double, 3> a_new_e{ 0, 0, 0 };
    bool e_ok = compute_exp_misp_energy_refine_polycube_three_cut(
        a_vc_size, local_energy_all, a_new_e,
        a_p_vc_pos_x, a_p_vc_pos_y, a_p_vc_pos_z,
        avv_vc_S, a_exp_vec, a_new_coord, alpha, beta, ga,
        a_bvf_size, a_p_vf_pos_x, a_p_vf_pos_y, a_p_vf_pos_z, a_mu_vec
    );
    if (!e_ok)
    {
        while (!e_ok)
        {
            g_all *= 0.2; move_d *= 0.2;
            a_new_coord[0][0] = dpx[v_id];
            a_new_coord[0][1] = dpy[v_id];
            a_new_coord[0][2] = dpz[v_id];
            if (move_d < 1e-8)
            {
                get_pair_coordinate(a_new_coord[0][0], a_new_coord[0][1], a_new_coord[0][2], type_matrix[adjacent_one_cut_type[0]], common_vertex_coord[0], common_vertex_coord[1], common_vertex_coord[2], a_new_coord[1][0], a_new_coord[1][1], a_new_coord[1][2]);
                get_pair_coordinate(a_new_coord[0][0], a_new_coord[0][1], a_new_coord[0][2], type_matrix[adjacent_one_cut_type[1]], common_vertex_coord[0], common_vertex_coord[1], common_vertex_coord[2], a_new_coord[2][0], a_new_coord[2][1], a_new_coord[2][2]);
                move_d = 1e-8; break;
            }
            else
            {
                a_new_coord[0] -= g_all;
                get_pair_coordinate(a_new_coord[0][0], a_new_coord[0][1], a_new_coord[0][2], type_matrix[adjacent_one_cut_type[0]], common_vertex_coord[0], common_vertex_coord[1], common_vertex_coord[2], a_new_coord[1][0], a_new_coord[1][1], a_new_coord[1][2]);
                get_pair_coordinate(a_new_coord[0][0], a_new_coord[0][1], a_new_coord[0][2], type_matrix[adjacent_one_cut_type[1]], common_vertex_coord[0], common_vertex_coord[1], common_vertex_coord[2], a_new_coord[2][0], a_new_coord[2][1], a_new_coord[2][2]);
            }
            a_new_e[0] = 0;
            e_ok = compute_exp_misp_energy_refine_polycube_three_cut(
                a_vc_size, local_energy_all, a_new_e,
                a_p_vc_pos_x, a_p_vc_pos_y, a_p_vc_pos_z,
                avv_vc_S, a_exp_vec, a_new_coord, alpha, beta, ga,
                a_bvf_size, a_p_vf_pos_x, a_p_vf_pos_y, a_p_vf_pos_z, a_mu_vec
            );
        }
    }
    dpx[v_id] = a_new_coord[0].x(); dpy[v_id] = a_new_coord[0].y(); dpz[v_id] = a_new_coord[0].z();
    get_pair_coordinate(a_new_coord[0][0], a_new_coord[0][1], a_new_coord[0][2], type_matrix[adjacent_one_cut_type[0]], common_vertex_coord[0], common_vertex_coord[1], common_vertex_coord[2], a_new_coord[1][0], a_new_coord[1][1], a_new_coord[1][2]);
    dpx[a_equal_vertex_id[1]] = a_new_coord[1].x(); dpy[a_equal_vertex_id[1]] = a_new_coord[1].y(); dpz[a_equal_vertex_id[1]] = a_new_coord[1].z();
    get_pair_coordinate(a_new_coord[0][0], a_new_coord[0][1], a_new_coord[0][2], type_matrix[adjacent_one_cut_type[1]], common_vertex_coord[0], common_vertex_coord[1], common_vertex_coord[2], a_new_coord[2][0], a_new_coord[2][1], a_new_coord[2][2]);
    dpx[a_equal_vertex_id[2]] = a_new_coord[2].x(); dpy[a_equal_vertex_id[2]] = a_new_coord[2].y(); dpz[a_equal_vertex_id[2]] = a_new_coord[2].z();
}
void polycube_deformation_interface::exp_mips_deformation_refine_one_polycube_normal_constrained_equal_face(int v_id, double angle_area_ratio, double energy_power, const int& omp_id, const double sym_M[3][3], int common_vert, double ga, bool update_all)
{
    bool bv_flag = false;
    if (bvf_id[v_id].size() > 0)
    {
        bv_flag = true;
    }
    int pair_id = vert_pairs_map[v_id];
    double local_energy_all = 0.0;
    double gx_all = 0.0; double gy_all = 0.0; double gz_all = 0.0;
    double min_radius = 1e30;
    double min_radius_pair = 1e30;
    const double* p_dpx = dpx.data(); const double* p_dpy = dpy.data(); const double* p_dpz = dpz.data();
    double x0_all = p_dpx[v_id]; double y0_all = p_dpy[v_id]; double z0_all = p_dpz[v_id];
    double center_x, center_y, center_z;
    center_x = p_dpx[common_vert];
    center_y = p_dpy[common_vert];
    center_z = p_dpz[common_vert];
    int vc_size = vertex_cell[v_id].size();
    int vc_size_pair = vertex_cell[pair_id].size();
    double* p_vc_pos_x = vc_pos_x2_omp[omp_id].data();
    double* p_vc_pos_y = vc_pos_y2_omp[omp_id].data();
    double* p_vc_pos_z = vc_pos_z2_omp[omp_id].data();
    std::vector<double> vc_pos_x2_omp_pair = vc_pos_x2_omp[omp_id];
    std::vector<double> vc_pos_y2_omp_pair = vc_pos_y2_omp[omp_id];
    std::vector<double> vc_pos_z2_omp_pair = vc_pos_z2_omp[omp_id];
    double* p_vc_pos_x_pair = vc_pos_x2_omp_pair.data();
    double* p_vc_pos_y_pair = vc_pos_y2_omp_pair.data();
    double* p_vc_pos_z_pair = vc_pos_z2_omp_pair.data();
    double* p_vc_n_cross_x = vc_n_cross_x_omp[omp_id].data();
    double* p_vc_n_cross_y = vc_n_cross_y_omp[omp_id].data();
    double* p_vc_n_cross_z = vc_n_cross_z_omp[omp_id].data();
    std::vector<double> vc_n_cross_x_omp_pair = vc_n_cross_x_omp[omp_id];
    std::vector<double> vc_n_cross_y_omp_pair = vc_n_cross_y_omp[omp_id];
    std::vector<double> vc_n_cross_z_omp_pair = vc_n_cross_z_omp[omp_id];
    double* p_vc_n_cross_x_pair = vc_n_cross_x_omp_pair.data();
    double* p_vc_n_cross_y_pair = vc_n_cross_y_omp_pair.data();
    double* p_vc_n_cross_z_pair = vc_n_cross_z_omp_pair.data();
    double* exp_vec = exp_vec_omp[omp_id].data();
    std::vector<double> exp_vec_omp_pair = exp_vec_omp[omp_id];
    double* exp_vec_pair = exp_vec_omp_pair.data();
    double* gx_vec = gx_vec_omp[omp_id].data();
    double* gy_vec = gy_vec_omp[omp_id].data();
    double* gz_vec = gz_vec_omp[omp_id].data();
    double* mu_vec = mu_vec_omp[omp_id].data();
    std::vector<double> gx_vec_omp_pair = gx_vec_omp[omp_id];
    std::vector<double> gy_vec_omp_pair = gy_vec_omp[omp_id];
    std::vector<double> gz_vec_omp_pair = gz_vec_omp[omp_id];
    std::vector<double> mu_vec_omp_pair = mu_vec_omp[omp_id];
    double* gx_vec_pair = gx_vec_omp_pair.data();
    double* gy_vec_pair = gy_vec_omp_pair.data();
    double* gz_vec_pair = gz_vec_omp_pair.data();
    double* mu_vec_pair = mu_vec_omp_pair.data();
    std::vector<double> vc_pos_x3_omp_pair = vc_pos_x3_omp[omp_id];
    std::vector<double> vc_pos_y3_omp_pair = vc_pos_y3_omp[omp_id];
    std::vector<double> vc_pos_z3_omp_pair = vc_pos_z3_omp[omp_id];
    double* p_vf_pos_x = vc_pos_x3_omp[omp_id].data(); double* p_vf_pos_y = vc_pos_y3_omp[omp_id].data(); double* p_vf_pos_z = vc_pos_z3_omp[omp_id].data();
    double* p_vf_pos_x_pair = vc_pos_x3_omp_pair.data(); double* p_vf_pos_y_pair = vc_pos_y3_omp_pair.data(); double* p_vf_pos_z_pair = vc_pos_z3_omp_pair.data();
    int bvf_size = 0;
    int bvf_size_pair = 0;
    std::vector<std::vector<double>>& vc_S = vc_S_omp[omp_id];
    std::vector<std::vector<double>> vc_S_pair = vc_S_omp[omp_id];
    double alpha = (1.0 - angle_area_ratio) * energy_power;
    double beta = energy_power - alpha;
    double normal_e_all = 0.0;
    double gx_new, gy_new, gz_new;
    compute_vertex_update_info(
        update_all, bv_flag, v_id,
        p_dpx, p_dpy, p_dpz,
        alpha, beta, vc_S,
        p_vf_pos_x, p_vf_pos_y, p_vf_pos_z,
        p_vc_pos_x, p_vc_pos_y, p_vc_pos_z,
        p_vc_n_cross_x, p_vc_n_cross_y, p_vc_n_cross_z,
        exp_vec, mu_vec, gx_vec, gy_vec, gz_vec,
        local_energy_all, normal_e_all, gx_new, gy_new, gz_new, min_radius, bvf_size
    );
    gx_all = gx_new; gy_all = gy_new; gz_all = gz_new;
    compute_vertex_update_info(
        update_all, bv_flag, pair_id,
        p_dpx, p_dpy, p_dpz,
        alpha, beta, vc_S_pair,
        p_vf_pos_x_pair, p_vf_pos_y_pair, p_vf_pos_z_pair,
        p_vc_pos_x_pair, p_vc_pos_y_pair, p_vc_pos_z_pair,
        p_vc_n_cross_x_pair, p_vc_n_cross_y_pair, p_vc_n_cross_z_pair,
        exp_vec_pair, mu_vec_pair, gx_vec_pair, gy_vec_pair, gz_vec_pair,
        local_energy_all, normal_e_all, gx_new, gy_new, gz_new, min_radius_pair, bvf_size_pair
    );
    gx_all += sym_M[0][0] * gx_new + sym_M[0][1] * gy_new + sym_M[0][2] * gz_new;
    gy_all += sym_M[1][0] * gx_new + sym_M[1][1] * gy_new + sym_M[1][2] * gz_new;
    gz_all += sym_M[2][0] * gx_new + sym_M[2][1] * gy_new + sym_M[2][2] * gz_new;
    min_radius = std::min(min_radius, min_radius_pair);
    double dx = gx_all, dy = gy_all, dz = gz_all;
    double move_d = sqrt(dx*dx + dy * dy + dz * dz);
    if (move_d > min_radius)
    {
        double t = min_radius / move_d;
        dx *= t; dy *= t; dz *= t;
        move_d = min_radius;
    }
    double npx = x0_all - dx; double npy = y0_all - dy; double npz = z0_all - dz;
    double npx_pair, npy_pair, npz_pair;
    get_pair_coordinate(npx, npy, npz, sym_M, center_x, center_y, center_z, npx_pair, npy_pair, npz_pair);
    if (1)
    {
        while (!local_check_negative_volume4(vc_size, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, p_vc_n_cross_x, p_vc_n_cross_y, p_vc_n_cross_z, npx, npy, npz))
        {
            dx *= 0.8; dy *= 0.8; dz *= 0.8; move_d *= 0.8;
            if (move_d < 1e-8)
            {
                npx = x0_all; npy = y0_all; npz = z0_all;
                move_d = 1e-8; break;
            }
            else
            {
                npx = x0_all - dx; npy = y0_all - dy; npz = z0_all - dz;
            }
        }
    }
    if (1)
    {
        while (!local_check_negative_volume4(vc_size_pair, p_vc_pos_x_pair, p_vc_pos_y_pair, p_vc_pos_z_pair, p_vc_n_cross_x_pair, p_vc_n_cross_y_pair, p_vc_n_cross_z_pair, npx_pair, npy_pair, npz_pair))
        {
            dx *= 0.8; dy *= 0.8; dz *= 0.8; move_d *= 0.8;
            if (move_d < 1e-8)
            {
                npx = x0_all; npy = y0_all; npz = z0_all;
                get_pair_coordinate(npx, npy, npz, sym_M, center_x, center_y, center_z, npx_pair, npy_pair, npz_pair);
                move_d = 1e-8; break;
            }
            else
            {
                npx = x0_all - dx; npy = y0_all - dy; npz = z0_all - dz;
                get_pair_coordinate(npx, npy, npz, sym_M, center_x, center_y, center_z, npx_pair, npy_pair, npz_pair);
            }
        }
    }
#if 1
    double new_e = 0;
    double new_e_pair = 0;
    bool e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy_all, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec, vc_size_pair, 0.0, new_e_pair, p_vc_pos_x_pair, p_vc_pos_y_pair, p_vc_pos_z_pair, vc_S_pair, exp_vec_pair, npx_pair, npy_pair, npz_pair, ga, bvf_size_pair, p_vf_pos_x_pair, p_vf_pos_y_pair, p_vf_pos_z_pair, mu_vec_pair);
    if (!e_ok)
    {
        while (!e_ok)
        {
            dx *= 0.2; dy *= 0.2; dz *= 0.2; move_d *= 0.2;
            if (move_d < 1e-8)
            {
                npx = x0_all; npy = y0_all; npz = z0_all;
                get_pair_coordinate(npx, npy, npz, sym_M, center_x, center_y, center_z, npx_pair, npy_pair, npz_pair);
                move_d = 1e-8; break;
            }
            else
            {
                npx = x0_all - dx; npy = y0_all - dy; npz = z0_all - dz;
                get_pair_coordinate(npx, npy, npz, sym_M, center_x, center_y, center_z, npx_pair, npy_pair, npz_pair);
            }
            new_e = 0;
            e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy_all, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec, vc_size_pair, 0.0, new_e_pair, p_vc_pos_x_pair, p_vc_pos_y_pair, p_vc_pos_z_pair, vc_S_pair, exp_vec_pair, npx_pair, npy_pair, npz_pair, ga, bvf_size_pair, p_vf_pos_x_pair, p_vf_pos_y_pair, p_vf_pos_z_pair, mu_vec_pair);
        }
        if (0)
        {
            dx *= 5; dy *= 5; dz *= 5; move_d *= 5; e_ok = false;
            while (!e_ok)
            {
                dx *= 0.5; dy *= 0.5; dz *= 0.5; move_d *= 0.5;
                if (move_d < 1e-8)
                {
                    npx = x0_all; npy = y0_all; npz = z0_all;
                    get_pair_coordinate(npx, npy, npz, sym_M, center_x, center_y, center_z, npx_pair, npy_pair, npz_pair);
                    move_d = 1e-8; break;
                }
                else
                {
                    npx = x0_all - dx; npy = y0_all - dy; npz = z0_all - dz;
                    get_pair_coordinate(npx, npy, npz, sym_M, center_x, center_y, center_z, npx_pair, npy_pair, npz_pair);
                }
                new_e = 0;
                e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy_all, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec, vc_size_pair, 0.0, new_e_pair, p_vc_pos_x_pair, p_vc_pos_y_pair, p_vc_pos_z_pair, vc_S_pair, exp_vec_pair, npx_pair, npy_pair, npz_pair, ga, bvf_size_pair, p_vf_pos_x_pair, p_vf_pos_y_pair, p_vf_pos_z_pair, mu_vec_pair);
            }
            if (0)
            {
                dx *= 2; dy *= 2; dz *= 2; move_d *= 2; e_ok = false;
                while (!e_ok)
                {
                    dx *= 0.875; dy *= 0.875; dz *= 0.875; move_d *= 0.875;
                    if (move_d < 1e-8)
                    {
                        npx = x0_all; npy = y0_all; npz = z0_all;
                        get_pair_coordinate(npx, npy, npz, sym_M, center_x, center_y, center_z, npx_pair, npy_pair, npz_pair);
                        move_d = 1e-8; break;
                    }
                    else
                    {
                        npx = x0_all - dx; npy = y0_all - dy; npz = z0_all - dz;
                        get_pair_coordinate(npx, npy, npz, sym_M, center_x, center_y, center_z, npx_pair, npy_pair, npz_pair);
                    }
                    new_e = 0;
                    e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy_all, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec, vc_size_pair, 0.0, new_e_pair, p_vc_pos_x_pair, p_vc_pos_y_pair, p_vc_pos_z_pair, vc_S_pair, exp_vec_pair, npx_pair, npy_pair, npz_pair, ga, bvf_size_pair, p_vf_pos_x_pair, p_vf_pos_y_pair, p_vf_pos_z_pair, mu_vec_pair);
                }
            }
        }
    }
#else
    double new_e = 0;
    bool e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
    while (!e_ok)
    {
        dx *= 0.8; dy *= 0.8; dz *= 0.8; move_d *= 0.8;
        if (move_d < 1e-8)
        {
            npx = x0; npy = y0; npz = z0;
            move_d = 1e-8; break;
        }
        else
        {
            npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
        }
        new_e = 0;
        e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
    }
#endif
    dpx[v_id] = npx; dpy[v_id] = npy; dpz[v_id] = npz;
    get_pair_coordinate(npx, npy, npz, sym_M, center_x, center_y, center_z, npx_pair, npy_pair, npz_pair);
    dpx[pair_id] = npx_pair; dpy[pair_id] = npy_pair; dpz[pair_id] = npz_pair;
}
void polycube_deformation_interface::get_pair_coordinate(double x, double y, double z, const double M[3][3], double center_x, double center_y, double center_z, double &x_pair, double & y_pair, double &z_pair)
{
    x_pair = 0.0;
    y_pair = 0.0;
    z_pair = 0.0;
    double rot_vec[3];
    double ori_vec[3] = { x - center_x, y - center_y, z - center_z };
    for (size_t b = 0; b < 3; b++)
    {
        rot_vec[b] = 0;
        for (size_t c = 0; c < 3; c++)
        {
            rot_vec[b] += M[b][c] * ori_vec[c];
        }
    }
    x_pair = center_x + rot_vec[0];
    y_pair = center_y + rot_vec[1];
    z_pair = center_z + rot_vec[2];
}
void polycube_deformation_interface::check_big_distortion_cell(
    const int& vc_size, std::vector<int>& large_distortion_flag,
    const double* posx, const double* posy, const double* posz,
    const std::vector<std::vector<double> >& vc_S,
    const double& npx, const double& npy, const double& npz)
{
    double D00, D10, D20, D01, D11, D21, D02, D12, D22;
    double C00, C01, C02, C10, C11, C12, C20, C21, C22;
    double A00, A01, A02, A10, A11, A12, A20, A21, A22;
    double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
    double AF, AF2, AF_I, e2, mips_e, volume_e, k;
    for (unsigned i = 0; i < vc_size; ++i)
    {
        const double* s_data = vc_S[i].data();
        int j = 3 * i;
        D00 = posx[j + 0] - npx; D10 = posy[j + 0] - npy; D20 = posz[j + 0] - npz;
        D01 = posx[j + 1] - npx; D11 = posy[j + 1] - npy; D21 = posz[j + 1] - npz;
        D02 = posx[j + 2] - npx; D12 = posy[j + 2] - npy; D22 = posz[j + 2] - npz;
        C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
        C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
        C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];
        A00 = C00 * D00 + C10 * D01 + C20 * D02;
        A10 = C00 * D10 + C10 * D11 + C20 * D12;
        A20 = C00 * D20 + C10 * D21 + C20 * D22;
        A01 = C01 * D00 + C11 * D01 + C21 * D02;
        A11 = C01 * D10 + C11 * D11 + C21 * D12;
        A21 = C01 * D20 + C11 * D21 + C21 * D22;
        A02 = C02 * D00 + C12 * D01 + C22 * D02;
        A12 = C02 * D10 + C12 * D11 + C22 * D12;
        A22 = C02 * D20 + C12 * D21 + C22 * D22;
        A2_00 = A00 * A00 + A10 * A10 + A20 * A20;
        A2_01 = A00 * A01 + A10 * A11 + A20 * A21;
        A2_02 = A00 * A02 + A10 * A12 + A20 * A22;
        A2_11 = A01 * A01 + A11 * A11 + A21 * A21;
        A2_12 = A01 * A02 + A11 * A12 + A21 * A22;
        A2_22 = A02 * A02 + A12 * A12 + A22 * A22;
        AF = A2_00 + A2_11 + A2_22;
        AF2 = A2_00 * A2_00 + A2_01 * A2_01 * 2 + A2_02 * A2_02 * 2 + A2_11 * A2_11 + A2_12 * A2_12 * 2 + A2_22 * A2_22;
        AF_I = (AF*AF - AF2)*0.5;
        det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
        i_det_A = 1 / det_A;
        e2 = AF * AF_I * (i_det_A*i_det_A);
        mips_e = (e2 - 1.0)*0.125;
        volume_e = 0.5*(det_A + i_det_A);
        k = (mips_e + volume_e)*0.5;
        if (k > bound_K) large_distortion_flag[i] = 1;
        else large_distortion_flag[i] = -1;
    }
}
bool polycube_deformation_interface::local_check_negative_volume4(const int& vc_size, const double* posx, const double* posy, const double* posz,
    const double* vc_n_cross_x, const double* vc_n_cross_y, const double* vc_n_cross_z,
    const double& npx, const double& npy, const double& npz)
{
    int j = 0;
    for (unsigned i = 0; i < vc_size; ++i)
    {
        j = 3 * i;
        double c_v = vc_n_cross_x[i] * (npx - posx[j]) + vc_n_cross_y[i] * (npy - posy[j]) + vc_n_cross_z[i] * (npz - posz[j]);
        if (c_v < 1e-20) return false;
    }
    return true;
}
Eigen::Vector3d polycube_deformation_interface::rotate_vector(const double &theta, const Eigen::Vector3d &axis, const Eigen::Vector3d &x)
{
    Eigen::Vector3d axis_normalized = axis / axis.norm();
    Eigen::Matrix3d rotation_matrix;
    double c = cos(theta), s = sin(theta);
    rotation_matrix <<
        axis_normalized[0] * axis_normalized[0] * (1 - c) + c, axis_normalized[0] * axis_normalized[1] * (1 - c) - axis_normalized[2] * s, axis_normalized[0] * axis_normalized[2] * (1 - c) + axis_normalized[1] * s,
        axis_normalized[0] * axis_normalized[1] * (1 - c) + axis_normalized[2] * s, axis_normalized[1] * axis_normalized[1] * (1 - c) + c, axis_normalized[1] * axis_normalized[2] * (1 - c) - axis_normalized[0] * s,
        axis_normalized[0] * axis_normalized[2] * (1 - c) - axis_normalized[1] * s, axis_normalized[1] * axis_normalized[2] * (1 - c) + axis_normalized[0] * s, axis_normalized[2] * axis_normalized[2] * (1 - c) + c;
    return rotation_matrix * x;
}
Eigen::Vector3d polycube_deformation_interface::project2plane(const Eigen::Vector3d &normal, const Eigen::Vector3d &d)
{
    Eigen::Vector3d n_normalized = normal.normalized();
    return d - d.dot(n_normalized) * n_normalized;
}
void polycube_deformation_interface::compute_vertex_update_info(
    const bool update_all, const bool bv_flag,
    const int &v_id,
    const double *&p_dpx, const double *&p_dpy, const double *&p_dpz,
    const double &alpha, const double &beta,
    std::vector<std::vector<double>> &vc_S,
    double *&p_vf_pos_x, double *&p_vf_pos_y, double *&p_vf_pos_z,
    double *&p_vc_pos_x, double *&p_vc_pos_y, double *&p_vc_pos_z,
    double *&p_vc_n_cross_x, double *&p_vc_n_cross_y, double *&p_vc_n_cross_z,
    double *&exp_vec, double *&mu_vec,
    double *&gx_vec, double *&gy_vec, double *&gz_vec,
    double &local_energy_all, double &normal_e_all,
    double &gx_new, double &gy_new, double &gz_new,
    double &min_radius, int &bvf_size
)
{
    double x0 = p_dpx[v_id]; double y0 = p_dpy[v_id]; double z0 = p_dpz[v_id];
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    double gx = 0.0; double gy = 0.0; double gz = 0.0;
    double len;
    unsigned jj = 0;
    const int* vc_id = vertex_cell[v_id].data();
    int vc_size = vertex_cell[v_id].size();
    int v_id_1, v_id_2, v_id_3, k;
    double D00, D10, D20, D01, D11, D21, D02, D12, D22;
    double C00, C01, C02, C10, C11, C12, C20, C21, C22;
    double A00, A01, A02, A10, A11, A12, A20, A21, A22;
    double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
    double AF, AF_05, I_AF_05, AF2, AF_I, AF_I_05, I_AF_I_05, exp_k, exp_e;
    double D_AF_x, D_AF_y, D_AF_z, D_AF2_x, D_AF2_y, D_AF2_z, D_AF_I_x, D_AF_I_y, D_AF_I_z, D_det_A_x, D_det_A_y, D_det_A_z, inv_det_A_2_05;
    double d_A00_x, d_A10_y, d_A20_z, d_A01_x, d_A11_y, d_A21_z, d_A02_x, d_A12_y, d_A22_z;
    double dvex, dvey, dvez, g, dgx, dgy, dgz, dex, dey, dez, e, mips_e, volume_e, d_mips_e_x, d_mips_e_y, d_mips_e_z;
    for (unsigned i = 0; i < vc_size; ++i)
    {
        const int* vv_id = vertex_cell_vertex[v_id][i].data();
        memcpy(&vc_S[i][0], &vcv_S[v_id][i][0], 9 * sizeof(double));
        double* s_data = vc_S[i].data();
        k = 3 * i;
        v_id_1 = vv_id[0];
        v_id_2 = vv_id[1];
        v_id_3 = vv_id[2];
        x1 = p_dpx[v_id_1]; x2 = p_dpx[v_id_2]; x3 = p_dpx[v_id_3];
        y1 = p_dpy[v_id_1]; y2 = p_dpy[v_id_2]; y3 = p_dpy[v_id_3];
        z1 = p_dpz[v_id_1]; z2 = p_dpz[v_id_2]; z3 = p_dpz[v_id_3];
        p_vc_pos_x[k + 0] = x1; p_vc_pos_x[k + 1] = x2; p_vc_pos_x[k + 2] = x3;
        p_vc_pos_y[k + 0] = y1; p_vc_pos_y[k + 1] = y2; p_vc_pos_y[k + 2] = y3;
        p_vc_pos_z[k + 0] = z1; p_vc_pos_z[k + 1] = z2; p_vc_pos_z[k + 2] = z3;
        p_vc_n_cross_x[i] = (y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1);
        p_vc_n_cross_y[i] = (x3 - x1)*(z2 - z1) - (x2 - x1)*(z3 - z1);
        p_vc_n_cross_z[i] = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
        D00 = x1 - x0; D10 = y1 - y0; D20 = z1 - z0;
        D01 = x2 - x0; D11 = y2 - y0; D21 = z2 - z0;
        D02 = x3 - x0; D12 = y3 - y0; D22 = z3 - z0;
        len = D00 * D00 + D10 * D10 + D20 * D20; if (len < min_radius) min_radius = len;
        len = D01 * D01 + D11 * D11 + D21 * D21; if (len < min_radius) min_radius = len;
        len = D02 * D02 + D12 * D12 + D22 * D22; if (len < min_radius) min_radius = len;
        C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
        C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
        C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];
        A00 = C00 * D00 + C10 * D01 + C20 * D02; d_A00_x = -C00 - C10 - C20;
        A10 = C00 * D10 + C10 * D11 + C20 * D12; d_A10_y = -C00 - C10 - C20;
        A20 = C00 * D20 + C10 * D21 + C20 * D22; d_A20_z = -C00 - C10 - C20;
        A01 = C01 * D00 + C11 * D01 + C21 * D02; d_A01_x = -C01 - C11 - C21;
        A11 = C01 * D10 + C11 * D11 + C21 * D12; d_A11_y = -C01 - C11 - C21;
        A21 = C01 * D20 + C11 * D21 + C21 * D22; d_A21_z = -C01 - C11 - C21;
        A02 = C02 * D00 + C12 * D01 + C22 * D02; d_A02_x = -C02 - C12 - C22;
        A12 = C02 * D10 + C12 * D11 + C22 * D12; d_A12_y = -C02 - C12 - C22;
        A22 = C02 * D20 + C12 * D21 + C22 * D22; d_A22_z = -C02 - C12 - C22;
        A2_00 = A00 * A00 + A10 * A10 + A20 * A20;
        A2_01 = A00 * A01 + A10 * A11 + A20 * A21;
        A2_02 = A00 * A02 + A10 * A12 + A20 * A22;
        A2_11 = A01 * A01 + A11 * A11 + A21 * A21;
        A2_12 = A01 * A02 + A11 * A12 + A21 * A22;
        A2_22 = A02 * A02 + A12 * A12 + A22 * A22;
        AF = A2_00 + A2_11 + A2_22; AF_05 = std::sqrt(AF); I_AF_05 = 0.5 / AF_05;
        AF2 = A2_00 * A2_00 + A2_01 * A2_01 * 2 + A2_02 * A2_02 * 2 + A2_11 * A2_11 + A2_12 * A2_12 * 2 + A2_22 * A2_22;
        AF_I = (AF*AF - AF2)*0.5; AF_I_05 = std::sqrt(AF_I); I_AF_I_05 = 0.5 / AF_I_05;
        det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
        i_det_A = 1.0 / det_A;
        D_AF_x = (2.0*A00*d_A00_x + 2.0*A01*d_A01_x + 2.0*A02*d_A02_x);
        D_AF_y = (2.0*A10*d_A10_y + 2.0*A11*d_A11_y + 2.0*A12*d_A12_y);
        D_AF_z = (2.0*A20*d_A20_z + 2.0*A21*d_A21_z + 2.0*A22*d_A22_z);
        D_AF2_x = (2.0*A2_00*2.0*A00*d_A00_x + 4.0*A2_01*(A01*d_A00_x + A00 * d_A01_x) + 4.0*A2_02*(A02*d_A00_x + A00 * d_A02_x)
            + 2.0*A2_11*2.0*A01*d_A01_x + 2.0*A2_22*2.0*A02*d_A02_x + 4.0*A2_12*(A01*d_A02_x + A02 * d_A01_x));
        D_AF2_y = (2.0*A2_00*2.0*A10*d_A10_y + 4.0*A2_01*(A11*d_A10_y + A10 * d_A11_y) + 4.0*A2_02*(A12*d_A10_y + A10 * d_A12_y)
            + 2.0*A2_11*2.0*A11*d_A11_y + 2.0*A2_22*2.0*A12*d_A12_y + 4.0*A2_12*(A11*d_A12_y + A12 * d_A11_y));
        D_AF2_z = (2.0*A2_00*2.0*A20*d_A20_z + 4.0*A2_01*(A21*d_A20_z + A20 * d_A21_z) + 4.0*A2_02*(A22*d_A20_z + A20 * d_A22_z)
            + 2.0*A2_11*2.0*A21*d_A21_z + 2.0*A2_22*2.0*A22*d_A22_z + 4.0*A2_12*(A21*d_A22_z + A22 * d_A21_z));
        D_AF_I_x = (AF*D_AF_x - 0.5*D_AF2_x)*I_AF_I_05; D_AF_I_y = (AF*D_AF_y - 0.5*D_AF2_y)*I_AF_I_05; D_AF_I_z = (AF*D_AF_z - 0.5*D_AF2_z)*I_AF_I_05;
        D_AF_x *= I_AF_05; D_AF_y *= I_AF_05; D_AF_z *= I_AF_05;
        D_det_A_x = d_A00_x * A11 * A22 + d_A01_x * A12 * A20 + d_A02_x * A10 * A21 - d_A02_x * A11 * A20 - d_A01_x * A10 * A22 - d_A00_x * A12 * A21;
        D_det_A_y = A00 * d_A11_y * A22 + A01 * d_A12_y * A20 + A02 * d_A10_y * A21 - A02 * d_A11_y * A20 - A01 * d_A10_y * A22 - A00 * d_A12_y * A21;
        D_det_A_z = A00 * A11 * d_A22_z + A01 * A12 * d_A20_z + A02 * A10 * d_A21_z - A02 * A11 * d_A20_z - A01 * A10 * d_A22_z - A00 * A12 * d_A21_z;
        inv_det_A_2_05 = 0.5 *(1.0 - i_det_A * i_det_A);
        dvex = D_det_A_x * inv_det_A_2_05;
        dvey = D_det_A_y * inv_det_A_2_05;
        dvez = D_det_A_z * inv_det_A_2_05;
        g = AF_05 * AF_I_05;
        dgx = D_AF_x * AF_I_05 + AF_05 * D_AF_I_x;
        dgy = D_AF_y * AF_I_05 + AF_05 * D_AF_I_y;
        dgz = D_AF_z * AF_I_05 + AF_05 * D_AF_I_z;
        dex = (dgx *i_det_A - g * D_det_A_x *i_det_A*i_det_A);
        dey = (dgy *i_det_A - g * D_det_A_y *i_det_A*i_det_A);
        dez = (dgz *i_det_A - g * D_det_A_z *i_det_A*i_det_A);
        e = g * i_det_A;
        mips_e = (e*e - 1.0)*0.125;
        volume_e = 0.5*(det_A + i_det_A);
        exp_k = (alpha*mips_e + beta * volume_e);
        if (exp_k > 60) exp_k = 60;
        exp_vec[i] = exp_k;
        d_mips_e_x = e * 0.25 * dex;
        d_mips_e_y = e * 0.25 * dey;
        d_mips_e_z = e * 0.25 * dez;
        gx_vec[i] = alpha * d_mips_e_x + beta * dvex;
        gy_vec[i] = alpha * d_mips_e_y + beta * dvey;
        gz_vec[i] = alpha * d_mips_e_z + beta * dvez;
    }
    for (int i = 0; i < vc_size; ++i)
    {
        exp_e = std::exp(exp_vec[i]);
        local_energy_all += exp_e;
        gx += gx_vec[i] * exp_e;
        gy += gy_vec[i] * exp_e;
        gz += gz_vec[i] * exp_e;
    }
    double mu = 1.0; int bvv_size = 0;
    if (bv_flag)
    {
        double normal_e = 0.0; double smooth_e = 0;
        double ne_gx = 0.0; double ne_gy = 0.0; double ne_gz = 0.0;
        double se_gx = 0.0; double se_gy = 0.0; double se_gz = 0.0;
        std::vector<int>& one_vv_id = bvv_id[v_id];
        std::vector<int>& one_vf_id = bvf_id[v_id];
        bvv_size = one_vv_id.size(); bvf_size = bvv_size;
        for (int i = 0; i < bvv_size; ++i)
        {
            int j = (i + 1) % bvv_size; int k = (i + 2) % bvv_size;
            x1 = p_dpx[one_vv_id[i]]; y1 = p_dpy[one_vv_id[i]]; z1 = p_dpz[one_vv_id[i]];
            x2 = p_dpx[one_vv_id[j]]; y2 = p_dpy[one_vv_id[j]]; z2 = p_dpz[one_vv_id[j]];
            x3 = p_dpx[one_vv_id[k]]; y3 = p_dpy[one_vv_id[k]]; z3 = p_dpz[one_vv_id[k]];
            p_vf_pos_x[4 * i + 0] = x1; p_vf_pos_x[4 * i + 1] = x2; p_vf_pos_x[4 * i + 2] = x3;
            p_vf_pos_y[4 * i + 0] = y1; p_vf_pos_y[4 * i + 1] = y2; p_vf_pos_y[4 * i + 2] = y3;
            p_vf_pos_z[4 * i + 0] = z1; p_vf_pos_z[4 * i + 1] = z2; p_vf_pos_z[4 * i + 2] = z3;
            OpenVolumeMesh::Geometry::Vec3d& tn = target_bfn[one_vf_id[i]];
            p_vf_pos_x[4 * i + 3] = tn[0]; p_vf_pos_y[4 * i + 3] = tn[1]; p_vf_pos_z[4 * i + 3] = tn[2];
            double nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
            double ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
            double nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);
            double len2 = nx * nx + ny * ny + nz * nz; double len = std::sqrt(len2);
            double E_ne = 0.5*((nx / len - tn[0])*(nx / len - tn[0]) + (ny / len - tn[1])*(ny / len - tn[1]) + (nz / len - tn[2])*(nz / len - tn[2]));
            if (update_all)
            {
                mu_vec[i] = local_energy_ratio[one_vf_id[i]];
            }
            else
            {
                mu_vec[i] = 1.0;
            }
            normal_e += mu_vec[i] * E_ne;
            normal_e_all += mu_vec[i] * E_ne;
            double d_nx_x0 = 0.0; double d_nx_y0 = z1 - z2; double d_nx_z0 = y2 - y1;
            double d_ny_x0 = z2 - z1; double d_ny_y0 = 0.0; double d_ny_z0 = x1 - x2;
            double d_nz_x0 = y1 - y2; double d_nz_y0 = x2 - x1; double d_nz_z0 = 0.0;
            double d_len2_x0 = (nx*d_nx_x0 + ny * d_ny_x0 + nz * d_nz_x0);
            double d_len2_y0 = (nx*d_nx_y0 + ny * d_ny_y0 + nz * d_nz_y0);
            double d_len2_z0 = (nx*d_nx_z0 + ny * d_ny_z0 + nz * d_nz_z0);
            double d_len_x0 = d_len2_x0 / len; double d_len_y0 = d_len2_y0 / len; double d_len_z0 = d_len2_z0 / len;
            double d_nx_len_x0 = d_nx_x0 / len - nx * d_len_x0 / (len*len);
            double d_nx_len_y0 = d_nx_y0 / len - nx * d_len_y0 / (len*len);
            double d_nx_len_z0 = d_nx_z0 / len - nx * d_len_z0 / (len*len);
            double d_ny_len_x0 = d_ny_x0 / len - ny * d_len_x0 / (len*len);
            double d_ny_len_y0 = d_ny_y0 / len - ny * d_len_y0 / (len*len);
            double d_ny_len_z0 = d_ny_z0 / len - ny * d_len_z0 / (len*len);
            double d_nz_len_x0 = d_nz_x0 / len - nz * d_len_x0 / (len*len);
            double d_nz_len_y0 = d_nz_y0 / len - nz * d_len_y0 / (len*len);
            double d_nz_len_z0 = d_nz_z0 / len - nz * d_len_z0 / (len*len);
            ne_gx += mu_vec[i] * (d_nx_len_x0*(nx / len - tn[0]) + d_ny_len_x0 * (ny / len - tn[1]) + d_nz_len_x0 * (nz / len - tn[2]));
            ne_gy += mu_vec[i] * (d_nx_len_y0*(nx / len - tn[0]) + d_ny_len_y0 * (ny / len - tn[1]) + d_nz_len_y0 * (nz / len - tn[2]));
            ne_gz += mu_vec[i] * (d_nx_len_z0*(nx / len - tn[0]) + d_ny_len_z0 * (ny / len - tn[1]) + d_nz_len_z0 * (nz / len - tn[2]));
        }
        local_energy_all += normal_e_all;
        gx_new = gx + ne_gx; gy_new = gy + ne_gy; gz_new = gz + ne_gz;
    }
    min_radius = std::sqrt(min_radius);
}
bool polycube_deformation_interface::compute_exp_misp_energy_refine_polycube(
    const int& vc_size, const double& old_e, double& new_e,
    const double* posx, const double* posy, const double* posz,
    const std::vector<std::vector<double> >& vc_S, double* exp_vec,
    const double& npx, const double& npy, const double& npz,
    double alpha, double beta,
    const double& ga, const int& vf_size,
    const double* vf_px, const double* vf_py, const double* vf_pz, const double* mu)
{
    new_e = 0.0;
    double D00, D10, D20, D01, D11, D21, D02, D12, D22;
    double C00, C01, C02, C10, C11, C12, C20, C21, C22;
    double A00, A01, A02, A10, A11, A12, A20, A21, A22;
    double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
    double AF, AF2, AF_I, e2, mips_e, volume_e, k, exp_e;
    for (unsigned i = 0; i < vc_size; ++i)
    {
        const double* s_data = vc_S[i].data();
        int j = 3 * i;
        D00 = posx[j + 0] - npx; D10 = posy[j + 0] - npy; D20 = posz[j + 0] - npz;
        D01 = posx[j + 1] - npx; D11 = posy[j + 1] - npy; D21 = posz[j + 1] - npz;
        D02 = posx[j + 2] - npx; D12 = posy[j + 2] - npy; D22 = posz[j + 2] - npz;
        C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
        C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
        C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];
        A00 = C00 * D00 + C10 * D01 + C20 * D02;
        A10 = C00 * D10 + C10 * D11 + C20 * D12;
        A20 = C00 * D20 + C10 * D21 + C20 * D22;
        A01 = C01 * D00 + C11 * D01 + C21 * D02;
        A11 = C01 * D10 + C11 * D11 + C21 * D12;
        A21 = C01 * D20 + C11 * D21 + C21 * D22;
        A02 = C02 * D00 + C12 * D01 + C22 * D02;
        A12 = C02 * D10 + C12 * D11 + C22 * D12;
        A22 = C02 * D20 + C12 * D21 + C22 * D22;
        det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
        A2_00 = A00 * A00 + A10 * A10 + A20 * A20;
        A2_01 = A00 * A01 + A10 * A11 + A20 * A21;
        A2_02 = A00 * A02 + A10 * A12 + A20 * A22;
        A2_11 = A01 * A01 + A11 * A11 + A21 * A21;
        A2_12 = A01 * A02 + A11 * A12 + A21 * A22;
        A2_22 = A02 * A02 + A12 * A12 + A22 * A22;
        i_det_A = 1.0 / det_A;
        AF = A2_00 + A2_11 + A2_22;
        AF2 = A2_00 * A2_00 + A2_01 * A2_01 * 2 + A2_02 * A2_02 * 2 + A2_11 * A2_11 + A2_12 * A2_12 * 2 + A2_22 * A2_22;
        AF_I = (AF*AF - AF2)*0.5;
        e2 = AF * AF_I * (i_det_A*i_det_A);
        mips_e = (e2 - 1.0)*0.125;
        volume_e = 0.5*(det_A + i_det_A);
        k = alpha * mips_e + beta * volume_e;
        if (k > 60) k = 60;
        exp_vec[i] = k;
    }
#if 0
    fmath::expd_v(exp_vec, vc_size);
    for (int i = 0; i < vc_size; ++i)
    {
        new_e += exp_vec[i];
    }
#else
    for (int i = 0; i < vc_size; ++i)
    {
        exp_e = std::exp(exp_vec[i]);
        new_e += exp_e;
    }
#endif
    if (new_e > old_e) return false;
    double tnx, tny, tnz; double x1, x2, y1, y2, z1, z2, x3, y3, z3;
    for (int i = 0; i < vf_size; ++i)
    {
        int j = 4 * i;
        x1 = vf_px[j + 0]; x2 = vf_px[j + 1]; x3 = vf_px[j + 2]; tnx = vf_px[j + 3];
        y1 = vf_py[j + 0]; y2 = vf_py[j + 1]; y3 = vf_py[j + 2]; tny = vf_py[j + 3];
        z1 = vf_pz[j + 0]; z2 = vf_pz[j + 1]; z3 = vf_pz[j + 2]; tnz = vf_pz[j + 3];
        double nx = (npy - y1)*(npz - z2) - (npy - y2)*(npz - z1);
        double ny = (npx - x2)*(npz - z1) - (npx - x1)*(npz - z2);
        double nz = (npx - x1)*(npy - y2) - (npx - x2)*(npy - y1);
        double len = std::sqrt(nx*nx + ny * ny + nz * nz); double inv_len = 1 / len;
        double nx_d = nx * inv_len - tnx; double ny_d = ny * inv_len - tny; double nz_d = nz * inv_len - tnz;
        double E_ne = nx_d * nx_d + ny_d * ny_d + nz_d * nz_d;
        new_e += 0.5*E_ne*mu[i];
        if (new_e > old_e) return false;
    }
    return true;
}
bool polycube_deformation_interface::compute_exp_misp_energy_refine_polycube(
    const int& vc_size, const double& old_e, double& new_e,
    const double* posx, const double* posy, const double* posz,
    const std::vector<std::vector<double> >& vc_S, double* exp_vec,
    const double& npx, const double& npy, const double& npz,
    double alpha, double beta,
    const double& ga, const int& vf_size,
    const double* vf_px, const double* vf_py, const double* vf_pz, const double* mu, const std::vector<int> &feature_face_pair, const std::vector<std::pair<int, int>> &feature_neighber_vert_localidx_pairidx)
{
    new_e = 0.0;
    double D00, D10, D20, D01, D11, D21, D02, D12, D22;
    double C00, C01, C02, C10, C11, C12, C20, C21, C22;
    double A00, A01, A02, A10, A11, A12, A20, A21, A22;
    double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
    double AF, AF2, AF_I, e2, mips_e, volume_e, k, exp_e;
    for (unsigned i = 0; i < vc_size; ++i)
    {
        const double* s_data = vc_S[i].data();
        int j = 3 * i;
        D00 = posx[j + 0] - npx; D10 = posy[j + 0] - npy; D20 = posz[j + 0] - npz;
        D01 = posx[j + 1] - npx; D11 = posy[j + 1] - npy; D21 = posz[j + 1] - npz;
        D02 = posx[j + 2] - npx; D12 = posy[j + 2] - npy; D22 = posz[j + 2] - npz;
        C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
        C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
        C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];
        A00 = C00 * D00 + C10 * D01 + C20 * D02;
        A10 = C00 * D10 + C10 * D11 + C20 * D12;
        A20 = C00 * D20 + C10 * D21 + C20 * D22;
        A01 = C01 * D00 + C11 * D01 + C21 * D02;
        A11 = C01 * D10 + C11 * D11 + C21 * D12;
        A21 = C01 * D20 + C11 * D21 + C21 * D22;
        A02 = C02 * D00 + C12 * D01 + C22 * D02;
        A12 = C02 * D10 + C12 * D11 + C22 * D12;
        A22 = C02 * D20 + C12 * D21 + C22 * D22;
        det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
        A2_00 = A00 * A00 + A10 * A10 + A20 * A20;
        A2_01 = A00 * A01 + A10 * A11 + A20 * A21;
        A2_02 = A00 * A02 + A10 * A12 + A20 * A22;
        A2_11 = A01 * A01 + A11 * A11 + A21 * A21;
        A2_12 = A01 * A02 + A11 * A12 + A21 * A22;
        A2_22 = A02 * A02 + A12 * A12 + A22 * A22;
        i_det_A = 1.0 / det_A;
        AF = A2_00 + A2_11 + A2_22;
        AF2 = A2_00 * A2_00 + A2_01 * A2_01 * 2 + A2_02 * A2_02 * 2 + A2_11 * A2_11 + A2_12 * A2_12 * 2 + A2_22 * A2_22;
        AF_I = (AF*AF - AF2)*0.5;
        e2 = AF * AF_I * (i_det_A*i_det_A);
        mips_e = (e2 - 1.0)*0.125;
        volume_e = 0.5*(det_A + i_det_A);
        k = alpha * mips_e + beta * volume_e;
        if (k > 60) k = 60;
        exp_vec[i] = k;
    }
#if 0
    fmath::expd_v(exp_vec, vc_size);
    for (int i = 0; i < vc_size; ++i)
    {
        new_e += exp_vec[i];
    }
#else
    for (int i = 0; i < vc_size; ++i)
    {
        exp_e = std::exp(exp_vec[i]);
        new_e += exp_e;
    }
#endif
    if (new_e > old_e) return false;
    double tnx, tny, tnz; double x1, x2, y1, y2, z1, z2, x3, y3, z3;
    for (int i = 0; i < vf_size; ++i)
    {
        int j = 4 * i;
        x1 = vf_px[j + 0]; x2 = vf_px[j + 1]; x3 = vf_px[j + 2]; tnx = vf_px[j + 3];
        y1 = vf_py[j + 0]; y2 = vf_py[j + 1]; y3 = vf_py[j + 2]; tny = vf_py[j + 3];
        z1 = vf_pz[j + 0]; z2 = vf_pz[j + 1]; z3 = vf_pz[j + 2]; tnz = vf_pz[j + 3];
        double nx = (npy - y1)*(npz - z2) - (npy - y2)*(npz - z1);
        double ny = (npx - x2)*(npz - z1) - (npx - x1)*(npz - z2);
        double nz = (npx - x1)*(npy - y2) - (npx - x2)*(npy - y1);
        double len = std::sqrt(nx*nx + ny * ny + nz * nz); double inv_len = 1 / len;
        double nx_d = nx * inv_len - tnx; double ny_d = ny * inv_len - tny; double nz_d = nz * inv_len - tnz;
        double E_ne = nx_d * nx_d + ny_d * ny_d + nz_d * nz_d;
        new_e += 0.5*E_ne*mu[i];
        if (new_e > old_e) return false;
    }
    if (feature_face_pair.size() != 0)
    {
        int pair_size = feature_face_pair.size() / 2;
        for (size_t i = 0; i < pair_size; i++)
        {
            int fid[2] = { feature_face_pair[2 * i], feature_face_pair[2 * i + 1] };
            double tmp_mu = FEATURE_COEFF;
            double nx[2], ny[2], nz[2];
            double nxn[2], nyn[2], nzn[2];
            double x1[2], x2[2], y1[2], y2[2], z1[2], z2[2];
            double len2[2], len[2];
            for (size_t j = 0; j < 2; j++)
            {
                int k1 = 4 * fid[j];
                x1[j] = vf_px[k1 + 0];
                y1[j] = vf_py[k1 + 0];
                z1[j] = vf_pz[k1 + 0];
                x2[j] = vf_px[k1 + 1];
                y2[j] = vf_py[k1 + 1];
                z2[j] = vf_pz[k1 + 1];
                nx[j] = (npy - y1[j])*(npz - z2[j]) - (npy - y2[j])*(npz - z1[j]);
                ny[j] = (npx - x2[j])*(npz - z1[j]) - (npx - x1[j])*(npz - z2[j]);
                nz[j] = (npx - x1[j])*(npy - y2[j]) - (npx - x2[j])*(npy - y1[j]);
                len2[j] = nx[j] * nx[j] + ny[j] * ny[j] + nz[j] * nz[j];
                len[j] = std::sqrt(len2[j]);
                nxn[j] = nx[j] / len[j];
                nyn[j] = ny[j] / len[j];
                nzn[j] = nz[j] / len[j];
            }
            double n1n2 = nxn[0] * nxn[1] + nyn[0] * nyn[1] + nzn[0] * nzn[1];
            double n1n22 = n1n2 * n1n2;
            new_e += tmp_mu * n1n22;
            if (new_e > old_e) return false;
        }
    }
    if (feature_neighber_vert_localidx_pairidx.size() != 0)
    {
        for (size_t i = 0; i < feature_neighber_vert_localidx_pairidx.size(); i++)
        {
            int local_fid = feature_neighber_vert_localidx_pairidx[i].first;
            int fid2 = feature_neighber_vert_localidx_pairidx[i].second;
            double tmp_mu = FEATURE_COEFF;
            double nx, ny, nz;
            double nxn[2], nyn[2], nzn[2];
            double x1, x2, y1, y2, z1, z2;
            double len2, len;
            for (size_t j = 0; j < 1; j++)
            {
                int k1 = 4 * local_fid;
                x1 = vf_px[k1 + 0];
                y1 = vf_py[k1 + 0];
                z1 = vf_pz[k1 + 0];
                x2 = vf_px[k1 + 1];
                y2 = vf_py[k1 + 1];
                z2 = vf_pz[k1 + 1];
                nx = (npy - y1)*(npz - z2) - (npy - y2)*(npz - z1);
                ny = (npx - x2)*(npz - z1) - (npx - x1)*(npz - z2);
                nz = (npx - x1)*(npy - y2) - (npx - x2)*(npy - y1);
                len2 = nx * nx + ny * ny + nz * nz;
                len = std::sqrt(len2);
                nxn[j] = nx / len;
                nyn[j] = ny / len;
                nzn[j] = nz / len;
            }
            nxn[1] = last_bfn[fid2][0];
            nyn[1] = last_bfn[fid2][1];
            nzn[1] = last_bfn[fid2][2];
            double n1n2 = nxn[0] * nxn[1] + nyn[0] * nyn[1] + nzn[0] * nzn[1];
            double n1n22 = n1n2 * n1n2;
            new_e += tmp_mu * n1n22;
            if (new_e > old_e) return false;
        }
    }
    return true;
}
bool polycube_deformation_interface::compute_exp_misp_energy_refine_polycube(const int& vc_size, const double& old_e, double& new_e,
    const double* posx, const double* posy, const double* posz,
    const std::vector<std::vector<double> >& vc_S, double* exp_vec,
    const double& npx, const double& npy, const double& npz,
    double alpha, double beta,
    const double& ga, const int& vf_size,
    const double* vf_px, const double* vf_py, const double* vf_pz, const double* mu,
    const int& vc_size_pair, const double& old_e_pair, double& new_e_pair,
    const double* posx_pair, const double* posy_pair, const double* posz_pair,
    const std::vector<std::vector<double> >& vc_S_pair, double* exp_vec_pair,
    const double& npx_pair, const double& npy_pair, const double& npz_pair,
    const double& ga_pair, const int& vf_size_pair,
    const double* vf_px_pair, const double* vf_py_pair, const double* vf_pz_pair, const double* mu_pair)
{
    new_e = 0.0;
    new_e_pair = 0.0;
    double new_e_total = 0.0;
    double old_e_total = old_e + old_e_pair;
    compute_exp_misp_energy_refine_polycube(vc_size, old_e, new_e, posx, posy, posz, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, vf_size, vf_px, vf_py, vf_pz, mu);
    compute_exp_misp_energy_refine_polycube(vc_size_pair, old_e_pair, new_e_pair, posx_pair, posy_pair, posz_pair, vc_S_pair, exp_vec_pair, npx_pair, npy_pair, npz_pair, alpha, beta, ga_pair, vf_size_pair, vf_px_pair, vf_py_pair, vf_pz_pair, mu_pair);
    new_e_total = new_e + new_e_pair;
    if (new_e_total > old_e_total) return false;
    return true;
}
bool polycube_deformation_interface::compute_exp_misp_energy_refine_polycube_three_cut(
    const std::array<double, 3> &a_vc_size, const double &old_e, std::array<double, 3> &a_new_e,
    const std::array<double *, 3> &a_posx, const std::array<double *, 3> &a_posy, const std::array<double *, 3> &a_posz,
    const std::array<std::vector<std::vector<double>>, 3> &avv_vc_S, std::array<double *, 3> &a_exp_vec,
    const std::array<Eigen::Vector3d, 3> &a_new_coord,
    const double &alpha, const double &beta, const double &ga, const std::array<int, 3> &a_vf_size,
    const std::array<double *, 3> &a_vf_px, const std::array<double *, 3> a_vf_py, const std::array<double *, 3> &a_vf_pz, const std::array<double *, 3> a_mu
)
{
    double sum_new_e = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        a_new_e[i] = 0.0;
        compute_exp_misp_energy_refine_polycube(
            a_vc_size[i], 0.0, a_new_e[i], a_posx[i], a_posy[i], a_posz[i],
            avv_vc_S[i], a_exp_vec[i], a_new_coord[i].x(), a_new_coord[i].y(), a_new_coord[i].z(),
            alpha, beta, ga, a_vf_size[i], a_vf_px[i], a_vf_py[i], a_vf_pz[i], a_mu[i]
        );
        sum_new_e += a_new_e[i];
    }
    if (sum_new_e > old_e)
        return false;
    return true;
}
