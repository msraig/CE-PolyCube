#pragma warning( disable : 4477 4018 4267 4244 4838)
#include "Polycube_Flattening.h"
#include "omp.h"
#include "Sparse_Matrix.h"
#include "Sparse_Solver.h"
#include "ConvexQuadOptimization.h"
#include "cholmod.h"
#include "SmallMat.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "mosek.h"
#include <Queue>
#include <fstream>
#include "permutohedral.h"
#include "adjust_orientation_LBFGS.h"
#include "Polycube_Boundary_Map_IVF.h"
#include "fmath.hpp"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include "Helper.h"
#define MIN_EDGE_DIST 1e-6
using namespace boost;
using ig::CVec;
polycube_flattening_interface::polycube_flattening_interface()
{
    prepare_ok = false;
    is_auto_deformation = true;
    deformation_v_id.clear(); deformation_new_p.clear(); deformation_ok.clear();
    handles_region_id.clear(); region_handles.clear();
    dis_th = 5.0;
    sigma_s = 1.0; sigma_r = 10; 
    asimptotic_epsilon = 1e-8;
    EpsSafetyFactor = 1.0e6;
    meps = 2.220446e-13;
    EpsilonEfectivo = meps * EpsSafetyFactor;
    umbral_factor = 1e-5;
    AABB_Tree = NULL;
    cut_to_chart_pair.clear();
    cut_to_chart_pair_neighbor.clear();
    cut_types.clear();
    cut_common_verts_idx.clear();
    vert_pairs_map.clear();
    equal_triangles.clear();
    three_cut_common_vert.clear();
    three_cut_vert.clear();
    three_cut_vert_flag.clear();
    three_cut_adjacent_one_cut_index.clear();
    hex_meshing_flag = true;
    update_updownchart_flag = false;
}
polycube_flattening_interface::~polycube_flattening_interface()
{
    if (AABB_Tree)
    {
        delete AABB_Tree;
    }
}
void polycube_flattening_interface::load_boundary_face_label(const std::vector<int> &label, const std::vector<int> &chart, TetStructure<double>* tet_mesh_)
{
    polycube_edges.clear();
    polycube_edge_distortion_iso.clear();
    polycube_edge_distortion_conf.clear();
    polycube_edge_distortion_vol.clear();
    polycube_edge_face_distortion.clear();
    polycube_edge_layer_distortion.clear();
    polycube_edge_length.clear();
    polycube_short_edges.clear();
    int nc = tet_mesh_->tetras.size();
    const std::vector<Tetrahedron<double>*> &tetras = tet_mesh_->tetras;
    const std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
    target_bfn.resize(nc, OpenVolumeMesh::Geometry::Vec3d(0, 0, 0));
    bf_chart.clear(); bf_chart.resize(nc, -1);
    std::vector<int> chart_xyz(nc, -1);
    int max_chart_id = 0;
    assert(label.size() == chart.size());
    if (label.size() == 0 || chart.size() == 0)
        return;
    for (size_t i = 0; i < label.size(); i++)
    {
        int l = label[i]; int cid = i;
        switch (l)
        {
        case 0:
            target_bfn[cid][0] = 1.0;
            chart_xyz[cid] = 0;
            break;
        case 1:
            target_bfn[cid][1] = 1.0;
            chart_xyz[cid] = 2;
            break;
        case 2:
            target_bfn[cid][2] = 1.0;
            chart_xyz[cid] = 4;
            break;
        case 3:
            target_bfn[cid][0] = -1.0;
            chart_xyz[cid] = 1;
            break;
        case 4:
            target_bfn[cid][1] = -1.0;
            chart_xyz[cid] = 3;
            break;
        case 5:
            target_bfn[cid][2] = -1.0;
            chart_xyz[cid] = 5;
            break;
        default:
            break;
        }
        bf_chart[cid] = chart[i];
        if (bf_chart[cid] > max_chart_id)
        {
            max_chart_id = bf_chart[cid];
        }
    }
    polycube_chart.clear(); polycube_chart.resize(max_chart_id + 1);
    polycube_chart_label.clear(); polycube_chart_label.resize(max_chart_id + 1);
    chart_mean_value.clear(); chart_mean_value.resize(max_chart_id + 1, 0.0);
    for (int i = 0; i < nc; ++i)
    {
        int chart_id = bf_chart[i];
        if (chart_id < 0) continue;
        polycube_chart[chart_id].push_back(i);
        polycube_chart_label[chart_id] = chart_xyz[i];
    }
    int nv = bvf_id.size();
    vertex_type.clear(); vertex_type.resize(nv, OpenVolumeMesh::Geometry::Vec3i(-1, -1, -1));
    for (int i = 0; i < nv; ++i)
    {
        std::vector<int>& one_vf_id = bvf_id[i];
        int bvf_size = one_vf_id.size();
        if (bvf_size > 0)
        {
            for (int j = 0; j < bvf_size; ++j)
            {
                vertex_type[i][chart_xyz[one_vf_id[j]] / 2] = bf_chart[one_vf_id[j]];
            }
        }
    }
    int ne = bef_id.size();
    std::vector<OpenVolumeMesh::Geometry::Vec3i> edge_two_chart;
    for (int i = 0; i < ne; ++i)
    {
        std::vector<int>& one_ef_id = bef_id[i];
        int bef_size = one_ef_id.size();
        if (bef_size == 2)
        {
            if (bf_chart[one_ef_id[0]] != bf_chart[one_ef_id[1]] && bf_chart[one_ef_id[1]] >= 0 && bf_chart[one_ef_id[0]] >= 0)
            {
                int l0 = bf_chart[one_ef_id[0]]; int l1 = bf_chart[one_ef_id[1]];
                if (l0 < l1)
                {
                    edge_two_chart.push_back(OpenVolumeMesh::Geometry::Vec3i(l0, l1, i));
                }
                else
                {
                    edge_two_chart.push_back(OpenVolumeMesh::Geometry::Vec3i(l1, l0, i));
                }
            }
        }
        else if (bef_size != 0)
        {
            printf("Error boundary edge face %d\n", bef_size);
        }
    }
    std::vector<int> empty_(1);
    std::vector<OpenVolumeMesh::Geometry::Vec2i> two_label;
    edge_with_same_label.clear();
    polycube_edge_idx2same_label_idx.clear();
    for (int i = 0; i < edge_two_chart.size(); ++i)
    {
        OpenVolumeMesh::Geometry::Vec3i& two_chart = edge_two_chart[i];
        int inserted_id = -1;
        for (int j = 0; j < two_label.size(); ++j)
        {
            if (two_chart[0] == two_label[j][0] && two_chart[1] == two_label[j][1])
            {
                inserted_id = j;
                edge_with_same_label[j].push_back(two_chart[2]);
            }
        }
        if (inserted_id < 0)
        {
            two_label.push_back(OpenVolumeMesh::Geometry::Vec2i(two_chart[0], two_chart[1]));
            int label1 = polycube_chart_label[two_chart[0]];
            int label2 = polycube_chart_label[two_chart[1]];
            if (label1 / 2 == label2 / 2)
            {
                std::cout << "!!!!!!!!!!!!!!!!!Error: illegal label!!!!!!!!!!!!!!!!!!!" << std::endl;
                std::cout << "Chart info: " << two_chart[0] << " " << two_chart[1] << std::endl;
                std::cout << "Label info: " << label1 << " " << label2 << std::endl;
                return;
            }
            empty_[0] = two_chart[2];
            edge_with_same_label.push_back(empty_);
        }
    }
    up_chart_id.clear(); down_chart_id.clear(); diff_up_down.clear();
    key_edges.clear();
    key_edges = edge_with_same_label;
    get_chart_distortion();
    for (int i = 0; i < edge_with_same_label.size(); ++i)
    {
        std::vector<int> two_corner;
        double edge_distortion = 0.0;
        double edge_distortion_conf = 0.0;
        double edge_distortion_vol = 0.0;
        for (int j = 0; j < edge_with_same_label[i].size(); ++j)
        {
            int edge_id = edge_with_same_label[i][j];
            std::pair<int, int> temp_pair = id2edge[edge_id];
            int fv = temp_pair.first; int tv = temp_pair.second;
            if (vertex_type[fv][0] >= 0 && vertex_type[fv][1] >= 0 && vertex_type[fv][2] >= 0)
            {
                two_corner.push_back(fv);
            }
            if (vertex_type[tv][0] >= 0 && vertex_type[tv][1] >= 0 && vertex_type[tv][2] >= 0)
            {
                two_corner.push_back(tv);
            }
            int face1 = bef_id[edge_id][0];
            int face2 = bef_id[edge_id][1];
            edge_distortion += all_iso_d[face1];
            edge_distortion += all_iso_d[face2];
            edge_distortion_vol += all_vol_d[face1];
            edge_distortion_vol += all_vol_d[face2];
            edge_distortion_conf += all_con_d[face1];
            edge_distortion_conf += all_con_d[face2];
        }
        edge_distortion = edge_distortion / edge_with_same_label[i].size();
        edge_distortion_vol = edge_distortion_vol / edge_with_same_label[i].size();
        edge_distortion_conf = edge_distortion_conf / edge_with_same_label[i].size();
        std::vector<OpenVolumeMesh::Geometry::Vec2i> final_tc;
        std::vector<std::pair<int, int>> final_tc_chart;
        if (two_corner.size() != 2)
        {
            std::vector<int> visited(two_corner.size(), -1);
            for (int ijk = 0; ijk < two_corner.size(); ++ijk)
            {
                if (visited[ijk] == 1) continue;
                visited[ijk] = 1; int seed_v = two_corner[ijk]; int last_v = -1;
                std::vector<int>& one_vv = bvv_id[seed_v]; int label_count = 0;
                for (int j = 0; j < one_vv.size(); ++j)
                {
                    label_count = 0;
                    for (int k = 0; k < 3; ++k)
                    {
                        if (vertex_type[one_vv[j]][k] == two_label[i][0] || vertex_type[one_vv[j]][k] == two_label[i][1])
                        {
                            ++label_count;
                        }
                    }
                    if (label_count == 2)
                    {
                        last_v = seed_v;
                        seed_v = one_vv[j];
                        break;
                    }
                }
                if (vertex_type[seed_v][0] >= 0 && vertex_type[seed_v][1] >= 0 && vertex_type[seed_v][2] >= 0)
                {
                    std::pair<int, int> temp_edge;
                    temp_edge.first = two_corner[ijk] < seed_v ? two_corner[ijk] : seed_v;
                    temp_edge.second = two_corner[ijk] + seed_v - temp_edge.first;
                    int edge_id = edge2id[temp_edge];
                    if (edge_id >= 0)
                    {
                        int l0 = bf_chart[bef_id[edge_id][0]]; int l1 = bf_chart[bef_id[edge_id][1]];
                        if (l0 > l1)
                        {
                            int temp = l0;
                            l0 = l1;
                            l1 = temp;
                        }
                        if (l0 == two_label[i][0] && l1 == two_label[i][1])
                        {
                            final_tc.push_back(OpenVolumeMesh::Geometry::Vec2i(two_corner[ijk], seed_v));
                            for (int kj = 0; kj < two_corner.size(); ++kj)
                            {
                                if (two_corner[kj] == seed_v) visited[kj] = 1;
                            }
                        }
                    }
                }
                else
                {
                    while (1)
                    {
                        one_vv = bvv_id[seed_v];
                        for (int j = 0; j < one_vv.size(); ++j)
                        {
                            label_count = 0;
                            for (int k = 0; k < 3; ++k)
                            {
                                if (vertex_type[one_vv[j]][k] == two_label[i][0] || vertex_type[one_vv[j]][k] == two_label[i][1])
                                {
                                    ++label_count;
                                }
                            }
                            if (label_count == 2 && one_vv[j] != last_v)
                            {
                                last_v = seed_v;
                                seed_v = one_vv[j];
                                break;
                            }
                        }
                        if (vertex_type[seed_v][0] >= 0 && vertex_type[seed_v][1] >= 0 && vertex_type[seed_v][2] >= 0)
                        {
                            final_tc.push_back(OpenVolumeMesh::Geometry::Vec2i(two_corner[ijk], seed_v));
                            for (int kj = 0; kj < two_corner.size(); ++kj)
                            {
                                if (two_corner[kj] == seed_v) visited[kj] = 1;
                            }
                            break;
                        }
                    }
                }
            }
        }
        else
        {
            final_tc.push_back(OpenVolumeMesh::Geometry::Vec2i(two_corner[0], two_corner[1]));
        }
        for (size_t j = 0; j < final_tc.size(); j++)
        {
            std::vector<int> one_edge;
            one_edge.push_back(final_tc[j][0]);
            one_edge.push_back(final_tc[j][1]);
            one_edge.push_back(two_label[i][0]);
            one_edge.push_back(two_label[i][1]);
            polycube_edges.push_back(one_edge);
            polycube_edge_idx2same_label_idx.push_back(i);
            polycube_edge_distortion_iso.push_back(edge_distortion);
            polycube_edge_distortion_conf.push_back(edge_distortion_conf);
            polycube_edge_distortion_vol.push_back(edge_distortion_vol);
            polycube_edge_face_distortion.push_back(chart_iso_d[two_label[i][0]] + chart_iso_d[two_label[i][1]]);
            double length = (dpx[final_tc[j][0]] - dpx[final_tc[j][1]]) * (dpx[final_tc[j][0]] - dpx[final_tc[j][1]]);
            length += (dpy[final_tc[j][0]] - dpy[final_tc[j][1]]) * (dpy[final_tc[j][0]] - dpy[final_tc[j][1]]);
            length += (dpz[final_tc[j][0]] - dpz[final_tc[j][1]]) * (dpz[final_tc[j][0]] - dpz[final_tc[j][1]]);
            length = sqrt(length);
            polycube_edge_length.push_back(length);
        }
        if (final_tc.size() == 1)
        {
            polycube_short_edges.push_back(edge_with_same_label[i]);
        }
        else
        {
            std::vector<std::vector<int>> local_bve(bvf_id.size());
            for (size_t j = 0; j < edge_with_same_label[i].size(); j++)
            {
                int edgeid = edge_with_same_label[i][j];
                int v0 = id2edge[edgeid].first;
                int v1 = id2edge[edgeid].second;
                local_bve[v0].push_back(edgeid);
                local_bve[v1].push_back(edgeid);
            }
            for (size_t j = 0; j < final_tc.size(); j++)
            {
                std::vector<int> one_long_edge;
                int beginp = final_tc[j][0];
                int endp = final_tc[j][1];
                assert(local_bve[beginp].size() == 1 && local_bve[endp].size() == 1);
                int prev_edge(-1), next_edge(local_bve[beginp][0]);
                int cur_p(beginp), next_p(-1);
                while (next_edge != -1)
                {
                    one_long_edge.push_back(next_edge);
                    std::pair<int, int> tmp_pair = id2edge[next_edge];
                    if (tmp_pair.first == cur_p)
                    {
                        next_p = tmp_pair.second;
                    }
                    else
                    {
                        next_p = tmp_pair.first;
                    }
                    cur_p = next_p;
                    prev_edge = next_edge;
                    if (cur_p == endp)
                    {
                        next_edge = -1;
                    }
                    else
                    {
                        int e0 = local_bve[cur_p][0];
                        int e1 = local_bve[cur_p][1];
                        if (prev_edge == e0)
                        {
                            next_edge = e1;
                        }
                        else
                        {
                            next_edge = e0;
                        }
                    }
                }
                polycube_short_edges.push_back(one_long_edge);
            }
        }
        for (int j = 0; j < final_tc.size(); ++j)
        {
            int v0 = final_tc[j][0]; int v1 = final_tc[j][1];
            int l0 = -1; int l1 = -1; int k0 = -1; int k1 = -1;
            int two_label_axis0 = polycube_chart_label[two_label[i][0]] / 2;
            int two_label_axis1 = polycube_chart_label[two_label[i][1]] / 2;
            for (int j = 0; j < 3; ++j)
            {
                int v0_axis = polycube_chart_label[vertex_type[v0][j]] / 2;
                int v1_axis = polycube_chart_label[vertex_type[v1][j]] / 2;
                if (two_label_axis0 != v0_axis && two_label_axis1 != v0_axis)
                {
                    l0 = vertex_type[v0][j]; k0 = j;
                }
                if (two_label_axis0 != v1_axis && two_label_axis1 != v1_axis)
                {
                    l1 = vertex_type[v1][j]; k1 = j;
                }
            }
            if (k0 != k1)
            {
                printf("Error corner label %d %d %d %d\n", v0, k0, v1, k1);
            }
            CVec<double, 3> p0 = tetra_vertices[v0]->pos;
            CVec<double, 3> p1 = tetra_vertices[v1]->pos;
            if (p0[k0] > p1[k1])
            {
                up_chart_id.push_back(l0); down_chart_id.push_back(l1); diff_up_down.push_back(p0[k0] - p1[k1]);
            }
            else
            {
                up_chart_id.push_back(l1); down_chart_id.push_back(l0); diff_up_down.push_back(p1[k1] - p0[k0]);
            }
        }
    }
    assert(polycube_edges.size() == polycube_short_edges.size());
}
void polycube_flattening_interface::load_feature_edges_vtk(const std::vector<std::pair<int, int>> &feature_edge_array)
{
    if (feature_edge_array.empty()) return;
    feature_edge_flag_ori.clear();
    feature_edge_flag_ori.resize(id2edge.size(), false);
    feature_v2v.resize(bvf_id.size());
    feature_v2e.resize(bvf_id.size());
    for (size_t i = 0; i < feature_edge_array.size(); i++)
    {
        int v0 = feature_edge_array[i].first;
        int v1 = feature_edge_array[i].second;
        if (v0 > v1)
        {
            int tmp = v0;
            v0 = v1;
            v1 = tmp;
        }
        auto it = edge2id.find(std::pair<int, int>(v0, v1));
        assert(it != edge2id.end());
        feature_edge_flag_ori[it->second] = true;
        feature_v2e[v0].push_back(it->second);
        feature_v2e[v1].push_back(it->second);
        feature_v2v[v0].push_back(v1);
        feature_v2v[v1].push_back(v0);
    }
    pqedge_v2e.clear();
    pqedge_v2e.resize(bvf_id.size());
    pqedge_v2v.clear();
    pqedge_v2v.resize(bvf_id.size());
    for (size_t i = 0; i < polycube_short_edges.size(); i++)
    {
        for (size_t j = 0; j < polycube_short_edges[i].size(); j++)
        {
            int eid = polycube_short_edges[i][j];
            int v0 = id2edge[eid].first;
            int v1 = id2edge[eid].second;
            pqedge_v2e[v0].push_back(eid);
            pqedge_v2e[v1].push_back(eid);
            pqedge_v2v[v0].push_back(v1);
            pqedge_v2v[v1].push_back(v0);
        }
    }
}
bool polycube_flattening_interface::feature_extraction(double select_ratio, int min_feature_length)
{
    if (feature_edge_flag_ori.empty()) return false;
    feature_polycube_edge.clear();
    group_feature_edge.clear();
    feature_edge_flag_final.clear();
    feature_edge_flag_final.resize(id2edge.size(), false);
    std::vector<int> one_feature_edge;
    bool success_flag = true;
    for (size_t i = 0; i < polycube_short_edges.size(); i++)
    {
        one_feature_edge.clear();
        for (size_t j = 0; j < polycube_short_edges[i].size(); j++)
        {
            int eid = polycube_short_edges[i][j];
            if (feature_edge_flag_ori[eid] == true)
                one_feature_edge.push_back(eid);
        }
        double ratio = one_feature_edge.size() * 1.0 / polycube_short_edges[i].size();
        if (one_feature_edge.size() == polycube_short_edges[i].size())
        {
            feature_polycube_edge.push_back(i);
            feature_polycube_edge_segm.push_back(std::pair<double, double>(0.0, 1.0));
            group_feature_edge.push_back(one_feature_edge);
            group_feature_edge_start.push_back(polycube_edges[i][0]);
            for (size_t j = 0; j < one_feature_edge.size(); j++)
            {
                feature_edge_flag_final[one_feature_edge[j]] = true;
            }
        }
        else
        {
            if (ratio > select_ratio)
                success_flag = false;
            std::set<int> one_feature_edge_set;
            for (size_t j = 0; j < one_feature_edge.size(); j++)
            {
                one_feature_edge_set.insert(one_feature_edge[j]);
            }
            std::vector<int> leftpart_eid, rightpart_eid;
            std::vector<std::pair<int, int>> pq_left, pq_right;
            for (size_t j = 0; j < polycube_short_edges[i].size(); j++)
            {
                int eid = polycube_short_edges[i][j];
                pq_left.push_back(id2edge[eid]);
            }
            pq_right = pq_left;
            reorder_edge(pq_left, polycube_edges[i][0]);
            reorder_edge(pq_right, polycube_edges[i][1]);
            int right_end_id = -1;
            for (size_t j = 0; j < pq_left.size(); j++)
            {
                int minid = std::min(pq_left[j].first, pq_left[j].second);
                int maxid = std::max(pq_left[j].first, pq_left[j].second);
                auto it = edge2id.find(std::pair<int, int>(minid, maxid));
                assert(it != edge2id.end());
                int eid = it->second;
                auto feait = one_feature_edge_set.find(eid);
                if (feait != one_feature_edge_set.end())
                {
                    leftpart_eid.push_back(eid);
                }
                else
                {
                    break;
                }
            }
            for (size_t j = 0; j < pq_right.size(); j++)
            {
                int minid = std::min(pq_right[j].first, pq_right[j].second);
                int maxid = std::max(pq_right[j].first, pq_right[j].second);
                auto it = edge2id.find(std::pair<int, int>(minid, maxid));
                assert(it != edge2id.end());
                int eid = it->second;
                auto feait = one_feature_edge_set.find(eid);
                if (feait != one_feature_edge_set.end())
                {
                    rightpart_eid.push_back(eid);
                    right_end_id = pq_right[j].second;
                }
                else
                {
                    break;
                }
            }
            if (leftpart_eid.size() > min_feature_length)
            {
                group_feature_edge.push_back(leftpart_eid);
                group_feature_edge_start.push_back(polycube_edges[i][0]);
                feature_polycube_edge.push_back(i);
                std::vector<std::pair<int, int>> shortlongedge;
                shortlongedge.push_back(std::pair<int, int>(polycube_edges[i][0], pq_left[leftpart_eid.size() - 1].second));
                shortlongedge.push_back(std::pair<int, int>(polycube_edges[i][0], polycube_edges[i][1]));
                std::vector<double> shortlonglength;
                get_edge_length(shortlongedge, shortlonglength);
                double ratio = shortlonglength[0] / shortlonglength[1];
                feature_polycube_edge_segm.push_back(std::pair<double, double>(0.0, ratio));
                for (size_t j = 0; j < leftpart_eid.size(); j++)
                {
                    feature_edge_flag_final[leftpart_eid[j]] = true;
                }
            }
            if (rightpart_eid.size() > min_feature_length)
            {
                group_feature_edge.push_back(rightpart_eid);
                assert(right_end_id != -1);
                group_feature_edge_start.push_back(right_end_id);
                feature_polycube_edge.push_back(i);
                std::vector<std::pair<int, int>> shortlongedge;
                shortlongedge.push_back(std::pair<int, int>(polycube_edges[i][1], pq_right[rightpart_eid.size() - 1].second));
                shortlongedge.push_back(std::pair<int, int>(polycube_edges[i][1], polycube_edges[i][0]));
                std::vector<double> shortlonglength;
                get_edge_length(shortlongedge, shortlonglength);
                double ratio = shortlonglength[0] / shortlonglength[1];
                feature_polycube_edge_segm.push_back(std::pair<double, double>(1.0 - ratio, 1.0));
                for (size_t j = 0; j < rightpart_eid.size(); j++)
                {
                    feature_edge_flag_final[rightpart_eid[j]] = true;
                }
            }
        }
    }
    return success_flag;
}
void polycube_flattening_interface::get_edge_length(const std::vector<std::pair<int, int>> &edgearray, std::vector<double> &edge_length)
{
    assert(dpx.size() != 0);
    edge_length.clear();
    for (size_t i = 0; i < edgearray.size(); i++)
    {
        int id0 = edgearray[i].first, id1 = edgearray[i].second;
        double len2 = (dpx[id0] - dpx[id1]) * (dpx[id0] - dpx[id1]) + (dpy[id0] - dpy[id1]) * (dpy[id0] - dpy[id1]) + (dpz[id0] - dpz[id1]) * (dpz[id0] - dpz[id1]);
        edge_length.push_back(std::sqrt(len2));
    }
}
bool polycube_flattening_interface::repair_chartlabel_feature(int max_line_length, int max_area_size)
{
    std::vector<double> edgelength;
    get_edge_length(id2edge, edgelength);
    std::vector<bool> pqe_flag(bef_id.size(), false);
    for (size_t i = 0; i < polycube_short_edges.size(); i++)
    {
        for (size_t j = 0; j < polycube_short_edges[i].size(); j++)
        {
            pqe_flag[polycube_short_edges[i][j]] = true;
        }
    }
    bool corner_correct = true;
    std::set<int> iso_edges;
    for (int i = 0; i < feature_edge_flag_ori.size(); i++)
    {
        if (feature_edge_flag_ori[i])
        {
            if (!pqe_flag[i]) iso_edges.insert(i);
        }
    }
    std::vector<std::vector<int>> bfe(bfv_id.size());
    for (size_t i = 0; i < bef_id.size(); i++)
    {
        if (bef_id[i].empty()) continue;
        assert(bef_id[i].size() == 2);
        int f0 = bef_id[i][0];
        int f1 = bef_id[i][1];
        bfe[f0].push_back(i);
        bfe[f1].push_back(i);
    }
    std::vector<std::vector<int>> grouped_iso_edges;
    while (!iso_edges.empty())
    {
        auto first = iso_edges.begin();
        std::vector<int> one_edge;
        std::vector<bool> color(bef_id.size(), false);
        color[*first] = true;
        std::queue<int> q;
        q.push(*first);
        while (!q.empty())
        {
            int front = q.front();
            one_edge.push_back(front);
            q.pop();
            int twovert[] = { id2edge[front].first, id2edge[front].second };
            for (size_t i = 0; i < 2; i++)
            {
                for (size_t j = 0; j < feature_v2e[twovert[i]].size(); j++)
                {
                    int eid = feature_v2e[twovert[i]][j];
                    if (eid == front) continue;
                    if (iso_edges.find(eid) != iso_edges.end() && color[eid] == false)
                    {
                        color[eid] = true;
                        q.push(eid);
                    }
                }
            }
        }
        grouped_iso_edges.push_back(one_edge);
        for (size_t i = 0; i < one_edge.size(); i++)
        {
            iso_edges.erase(one_edge[i]);
        }
    }
    std::vector<std::vector<std::pair<int, int>>> grouped_edges_reordered;
    std::vector<int> valence_three_case_id;
    for (size_t i = 0; i < grouped_iso_edges.size(); i++)
    {
        if (grouped_iso_edges[i].size() > max_line_length) continue;
        std::map<int, std::vector<int>> local_v2e;
        for (size_t j = 0; j < grouped_iso_edges[i].size(); j++)
        {
            int eid = grouped_iso_edges[i][j];
            local_v2e[id2edge[eid].first].push_back(eid);
            local_v2e[id2edge[eid].second].push_back(eid);
        }
        std::vector<int> endpoints;
        bool valence_three_flag = false;
        for (auto it : local_v2e)
        {
            if (it.second.size() >= 3)
            {
                valence_three_flag = true;
                corner_correct = false;
                valence_three_error_feature_pt_id.push_back(it.first);
                valence_three_case_id.push_back(i);
            }
            if (it.second.size() == 1) endpoints.push_back(it.first);
        }
        if (valence_three_flag || endpoints.size() != 2) continue;
        bool legal_endpoint_flag = true;
        for (size_t j = 0; j < endpoints.size(); j++)
        {
            int vid = endpoints[j];
            if (feature_v2e[vid].size() == 1) legal_endpoint_flag = false;
            if (pqedge_v2e[vid].empty()) legal_endpoint_flag = false;
        }
        if (!legal_endpoint_flag) continue;
        std::vector<std::pair<int, int>> one_order_edge;
        for (size_t j = 0; j < grouped_iso_edges[i].size(); j++)
        {
            int eid = grouped_iso_edges[i][j];
            one_order_edge.push_back(id2edge[eid]);
        }
        reorder_edge(one_order_edge, endpoints[0]);
        grouped_edges_reordered.push_back(one_order_edge);
    }
    std::vector<std::vector<std::pair<int, int>>> grouped_edges_reordered_split;
    for (size_t i = 0; i < grouped_edges_reordered.size(); i++)
    {
        std::vector<int> split_id;
        for (int j = 1; j < grouped_edges_reordered[i].size(); j++)
        {
            int vid = grouped_edges_reordered[i][j].first;
            if (pqedge_v2e[vid].size() != 0) split_id.push_back(j);
        }
        split_id.push_back(grouped_edges_reordered[i].size());
        if (split_id.size() == 1)
        {
            grouped_edges_reordered_split.push_back(grouped_edges_reordered[i]);
        }
        else
        {
            int start = 0;
            for (size_t j = 0; j < split_id.size(); j++)
            {
                int end = split_id[j];
                std::vector<std::pair<int, int>> one_edge;
                for (size_t k = start; k < end; k++)
                {
                    one_edge.push_back(grouped_edges_reordered[i][k]);
                }
                grouped_edges_reordered_split.push_back(one_edge);
                start = end;
            }
        }
    }
    std::vector<std::vector<std::pair<int, int>>> grouped_loop;
    typedef adjacency_list < listS, vecS, undirectedS,
        no_property, property < edge_weight_t, double > > graph_t;
    typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
    typedef std::pair<int, int> Edge;
    int n_pqv = 0;
    for (size_t i = 0; i < pqedge_v2e.size(); i++)
    {
        if (!pqedge_v2e[i].empty()) n_pqv++;
    }
    std::vector<Edge> edge_array;
    std::vector<double> weights;
    for (size_t i = 0; i < polycube_short_edges.size(); i++)
    {
        for (size_t j = 0; j < polycube_short_edges[i].size(); j++)
        {
            int first = id2edge[polycube_short_edges[i][j]].first;
            int second = id2edge[polycube_short_edges[i][j]].second;
            edge_array.push_back(Edge(first, second));
            double d2 = (dpx[first] - dpx[second]) * (dpx[first] - dpx[second]) + (dpy[first] - dpy[second]) * (dpy[first] - dpy[second]) + (dpz[first] - dpz[second]) * (dpz[first] - dpz[second]);
            weights.push_back(std::sqrt(d2));
        }
    }
    graph_t g(edge_array.begin(), edge_array.end(), weights.begin(), n_pqv);
    property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
    for (size_t i = 0; i < grouped_edges_reordered_split.size(); i++)
    {
        grouped_loop.push_back(grouped_edges_reordered_split[i]);
        int end = grouped_edges_reordered_split[i].back().second;
        int start = grouped_edges_reordered_split[i].front().first;
        std::vector<vertex_descriptor> p(num_vertices(g));
        std::vector<double> d(num_vertices(g));
        vertex_descriptor s = vertex(start, g);
        dijkstra_shortest_paths(g, s,
            predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, g))).
            distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, g))));
        vertex_descriptor t = vertex(end, g);
        vertex_descriptor par = p[t];
        if (d[t] < DBL_MAX)
        {
            while (t != s)
            {
                grouped_loop.back().push_back(std::pair<int, int>(t, par));
                t = par;
                par = p[t];
            }
            assert(grouped_loop.back().front().first == grouped_loop.back().back().second);
        }
        else
        {
            grouped_loop.pop_back();
            continue;
        }
        if (grouped_loop.back().size() > 3 * grouped_edges_reordered_split[i].size())
        {
            grouped_loop.pop_back();
            continue;
        }
    }
    for (size_t i = 0; i < grouped_loop.size(); i++)
    {
        std::vector<bool> sel_face_flag(bfv_id.size(), false);
        std::set<int> one_loop_edge;
        for (size_t j = 0; j < grouped_loop[i].size(); j++)
        {
            int first = grouped_loop[i][j].first;
            int second = grouped_loop[i][j].second;
            if (first > second) std::swap(first, second);
            auto it = edge2id.find(std::pair<int, int>(first, second));
            assert(it != edge2id.end());
            one_loop_edge.insert(it->second);
        }
        bool all_feature_flag = true;
        for (auto e : one_loop_edge)
        {
            if (!feature_edge_flag_ori[e])
                all_feature_flag = false;
        }
        if (all_feature_flag) continue;
        std::vector<int> sel_faces;
        int seedface = bef_id[*one_loop_edge.begin()][0];
        sel_face_flag[seedface] = true;
        std::queue<int> q;
        q.push(seedface);
        while (!q.empty())
        {
            int front = q.front();
            sel_faces.push_back(q.front());
            q.pop();
            for (size_t j = 0; j < bfe[front].size(); j++)
            {
                int eid = bfe[front][j];
                if (one_loop_edge.find(eid) == one_loop_edge.end())
                {
                    int otherface = bef_id[eid][0];
                    if (otherface == front) otherface = bef_id[eid][1];
                    if (!sel_face_flag[otherface])
                    {
                        sel_face_flag[otherface] = true;
                        q.push(otherface);
                    }
                }
            }
        }
        int sel_face_count = sel_faces.size();
        if (sel_face_count > boundary_face_number - sel_face_count)
        {
            sel_faces.clear();
            for (size_t j = 0; j < sel_face_flag.size(); j++)
            {
                if (!sel_face_flag[j] && bfv_id[j][0] != -1)
                {
                    sel_faces.push_back(j);
                }
            }
            assert(sel_faces.size() == boundary_face_number - sel_face_count);
            sel_face_count = sel_faces.size();
        }
        if (sel_face_count > max_area_size || sel_face_count == 0) continue;
        std::set<int> sel_face_charts;
        for (size_t j = 0; j < sel_faces.size(); j++)
        {
            sel_face_charts.insert(bf_chart[sel_faces[j]]);
        }
        if (sel_face_charts.size() >= 3) continue;
        std::vector<double> loop_chart_count(polycube_chart_label.size(), 0.0);
        for (auto e : one_loop_edge)
        {
            int f0 = bef_id[e][0];
            int f1 = bef_id[e][1];
            loop_chart_count[bf_chart[f0]] += edgelength[e];
            loop_chart_count[bf_chart[f1]] += edgelength[e];
        }
        int final_chart = -1;
        double max_chart_count = -1;
        for (size_t j = 0; j < polycube_chart_label.size(); j++)
        {
            int chart = j;
            if (sel_face_charts.find(chart) == sel_face_charts.end())
            {
                if (max_chart_count < loop_chart_count[chart])
                {
                    max_chart_count = loop_chart_count[chart];
                    final_chart = chart;
                }
            }
        }
        assert(final_chart != -1 && max_chart_count != 0);
        for (size_t j = 0; j < sel_faces.size(); j++)
        {
            bf_chart[sel_faces[j]] = final_chart;
        }
    }
    for (size_t i = 0; i < valence_three_case_id.size(); i++)
    {
        std::map<int, std::vector<int>> local_v2e;
        int id = valence_three_case_id[i];
        for (size_t j = 0; j < grouped_iso_edges[id].size(); j++)
        {
            int eid = grouped_iso_edges[id][j];
            local_v2e[id2edge[eid].first].push_back(eid);
            local_v2e[id2edge[eid].second].push_back(eid);
        }
        for (auto it : local_v2e)
        {
            if (it.second.size() >= 3)
            {
                std::vector<std::set<int>> split_face_id;
                std::set<int> one_ring_face(bvf_id[it.first].begin(), bvf_id[it.first].end());
                std::map<int, bool> visited;
                for (auto f : one_ring_face)
                    visited[f] = false;
                while (!one_ring_face.empty())
                {
                    std::set<int> selected;
                    std::map<int, bool> cur_visited = visited;
                    auto first = one_ring_face.begin();
                    cur_visited[*first] = true;
                    std::queue<int> q;
                    q.push(*first);
                    while (!q.empty())
                    {
                        int front = q.front();
                        q.pop();
                        selected.insert(front);
                        for (size_t j = 0; j < bfe[front].size(); j++)
                        {
                            int eid = bfe[front][j];
                            if (feature_edge_flag_ori[eid]) continue;
                            int otherface = bef_id[eid][0];
                            if (otherface == front) otherface = bef_id[eid][1];
                            if (one_ring_face.find(otherface) != one_ring_face.end() && visited[otherface] != true)
                            {
                                q.push(otherface);
                                visited[otherface] = true;
                            }
                        }
                    }
                    split_face_id.push_back(selected);
                    for (auto f : selected)
                    {
                        one_ring_face.erase(f);
                    }
                }
                assert(split_face_id.size() == it.second.size());
                for (size_t j = 0; j < split_face_id.size(); j++)
                {
                    std::set<int> outsideface = split_face_id[j];
                    int n_search_layer = 6;
                    for (size_t k = 0; k < n_search_layer; k++)
                    {
                        std::set<int> added_face;
                        for (auto f : outsideface)
                        {
                            for (auto e : bfe[f])
                            {
                                if (!feature_edge_flag_ori[e])
                                {
                                    int otherface = bef_id[e][0];
                                    if (otherface == f) otherface = bef_id[e][1];
                                    if (outsideface.find(otherface) == outsideface.end())
                                    {
                                        added_face.insert(otherface);
                                    }
                                }
                            }
                        }
                        for (auto f : added_face)
                            split_face_id[j].insert(f);
                        outsideface = added_face;
                    }
                }
                for (size_t j = 0; j < split_face_id.size(); j++)
                {
                    int max_chart = -1;
                    int max_chart_count = -1;
                    std::vector<int> chart_count(polycube_chart_label.size(), 0);
                    for (auto f : split_face_id[j])
                    {
                        chart_count[bf_chart[f]]++;
                    }
                    for (size_t k = 0; k < polycube_chart_label.size(); k++)
                    {
                        if (max_chart_count == -1 || max_chart_count < chart_count[k])
                        {
                            max_chart_count = chart_count[k];
                            max_chart = k;
                        }
                    }
                    assert(max_chart != -1);
                    for (auto f : split_face_id[j])
                    {
                        bf_chart[f] = max_chart;
                    }
                }
            }
        }
    }
    return corner_correct;
}
bool polycube_flattening_interface::save_chartlabel(const char* filename)
{
    static int normalmap[] = { 0, 1, 2, 3, 4, 5 };
    std::ofstream ofs(filename);
    for (size_t i = 0; i < bf_chart.size(); i++)
    {
        if (bf_chart[i] < 0)
            ofs << i << " -1 -1\n";
        else
        {
            ofs << i << " " << bf_chart[i] << " " << normalmap[polycube_chart_label[bf_chart[i]]] << std::endl;
        }
    }
    ofs.close();
    return true;
}
void polycube_flattening_interface::reorder_edge(std::vector<std::pair<int, int>> &one_edge, int start)
{
    std::map<int, std::vector<int>> vert2edgeid;
    for (int i = 0; i < one_edge.size(); i++)
    {
        vert2edgeid[one_edge[i].first].push_back(i);
        vert2edgeid[one_edge[i].second].push_back(i);
    }
    auto findstart = vert2edgeid.find(start);
    assert(findstart != vert2edgeid.end());
    assert(findstart->second.size() == 1);
    if (findstart == vert2edgeid.end() || findstart->second.size() != 1)
    {
        return;
    }
    int start_eid = findstart->second[0];
    std::vector<std::pair<int, int>> one_edge_new;
    int cur = start, prev = -1, next = one_edge[start_eid].first;
    next = next == cur ? one_edge[start_eid].second : next;
    while (next != -1)
    {
        one_edge_new.push_back(std::pair<int, int>(cur, next));
        prev = cur;
        cur = next;
        next = -1;
        auto it = vert2edgeid.find(cur);
        assert(it != vert2edgeid.end());
        for (size_t i = 0; i < it->second.size(); i++)
        {
            int eid = it->second[i];
            int diffv = one_edge[eid].first == cur ? one_edge[eid].second : one_edge[eid].first;
            if (diffv != prev)
            {
                next = diffv;
                break;
            }
        }
    }
    assert(one_edge_new.size() == one_edge.size());
    one_edge = one_edge_new;
}
void polycube_flattening_interface::save_feature_edges_vtk_feaformat(const char* filename)
{
    if (group_feature_edge.empty()) return;
    std::ofstream ofs(filename);
    std::vector<std::pair<int, int>> long_feature;
    for (size_t i = 0; i < group_feature_edge.size(); i++)
    {
        std::vector<std::pair<int, int>> one_short_edge_pair;
        for (size_t j = 0; j < group_feature_edge[i].size(); j++)
        {
            int eid = group_feature_edge[i][j];
            one_short_edge_pair.push_back(std::pair<int, int>(id2edge[eid].first, id2edge[eid].second));
        }
        long_feature.insert(long_feature.end(), one_short_edge_pair.begin(), one_short_edge_pair.end());
    }
    ofs << long_feature.size() << std::endl;
    for (size_t i = 0; i < long_feature.size(); i++)
    {
        ofs << long_feature[i].first << " " << long_feature[i].second << std::endl;
    }
    assert(!feature_polycube_edge_segm.empty());
    ofs << "Feature PolyCube Edge" << std::endl;
    int useful_edge_count = 0;
    for (size_t i = 0; i < feature_polycube_edge.size(); i++)
    {
        int eid = feature_polycube_edge[i];
        if (feature_polycube_edge_segm[i].first < DBL_EPSILON && 1.0 - feature_polycube_edge_segm[i].second < DBL_EPSILON)
        {
            useful_edge_count++;
        }
    }
    ofs << useful_edge_count << std::endl;
    for (size_t i = 0; i < feature_polycube_edge.size(); i++)
    {
        int eid = feature_polycube_edge[i];
        if (feature_polycube_edge_segm[i].first < DBL_EPSILON && 1.0 - feature_polycube_edge_segm[i].second < DBL_EPSILON)
        {
            ofs << polycube_edges[eid][0] << " " << polycube_edges[eid][1] << std::endl;
        }
    }
    ofs.close();
}
void polycube_flattening_interface::save_feature_edges_vtk_tfeformat(const char* filename, TetStructure<double>* tet_mesh_ori)
{
    if (group_feature_edge.empty()) return;
    std::ofstream ofs(filename);
    ofs << "Tet Feature Line" << std::endl;
    ofs << group_feature_edge.size() << std::endl;
    for (size_t i = 0; i < group_feature_edge.size(); i++)
    {
        ofs << group_feature_edge[i].size() + 1 << std::endl;
        std::vector<std::pair<int, int>> one_short_edge_pair;
        for (size_t j = 0; j < group_feature_edge[i].size(); j++)
        {
            int eid = group_feature_edge[i][j];
            one_short_edge_pair.push_back(std::pair<int, int>(id2edge[eid].first, id2edge[eid].second));
        }
        reorder_edge(one_short_edge_pair, group_feature_edge_start[i]);
        ofs << tet_mesh_ori->tetra_vertices[one_short_edge_pair[0].first]->pos << " ";
        for (size_t j = 0; j < one_short_edge_pair.size(); j++)
        {
            int vid = one_short_edge_pair[j].second;
            ofs << tet_mesh_ori->tetra_vertices[vid]->pos << " ";
        }
        ofs << std::endl;
    }
    ofs.close();
}
void polycube_flattening_interface::save_error_feature_ptid_ptsformat(const char* filename)
{
    std::ofstream ofs(filename);
    ofs << valence_three_error_feature_pt_id.size() << std::endl;
    for (size_t i = 0; i < valence_three_error_feature_pt_id.size(); i++)
    {
        ofs << valence_three_error_feature_pt_id[i] << std::endl;
    }
    ofs.close();
}
void polycube_flattening_interface::save_feature_edges_vtk_vtkformat(const char* filename)
{
    if ((feature_edge_flag_final).empty()) return;
    std::vector<int> feature_edge_array;
    for (size_t i = 0; i < (feature_edge_flag_final).size(); i++)
    {
        if ((feature_edge_flag_final)[i] == true) feature_edge_array.push_back(i);
    }
    std::ofstream outputfile(filename);
    outputfile << "# vtk DataFile Version 3.0\n"
        << "mesh vtk data\n"
        << "ASCII\n"
        << "DATASET POLYDATA\n";
    int n_point = dpx.size();
    outputfile << "POINTS " << n_point << " double" << std::endl;
    for (size_t i = 0; i < n_point; i++)
    {
        outputfile << dpx[i] << " " << dpy[i] << " " << dpz[i] << std::endl;
    }
    std::vector<int> feature_point_color(n_point, 0);
    outputfile << "LINES " << feature_edge_array.size() << " " << 3 * feature_edge_array.size() << std::endl;
    for (size_t i = 0; i < feature_edge_array.size(); i++)
    {
        int id0 = id2edge[feature_edge_array[i]].first;
        int id1 = id2edge[feature_edge_array[i]].second;
        feature_point_color[id0] = 1;
        feature_point_color[id1] = 1;
        outputfile << "2 " << id0 << " " << id1 << std::endl;
    }
    outputfile << "POINT_DATA " << n_point << "\n"
        << "SCALARS V_Scalars int\nLOOKUP_TABLE V_Table" << std::endl;
    for (size_t i = 0; i < n_point; i++)
    {
        outputfile << feature_point_color[i] << std::endl;
    }
    outputfile.close();
}
void polycube_flattening_interface::save_feature_long_edges_vtk_vtkformat(const char* filename)
{
    std::vector<bool> polycube_vertex_flag(dpx.size(), false);
    for (size_t i = 0; i < polycube_edges.size(); i++)
    {
        polycube_vertex_flag[polycube_edges[i][0]] = true;
        polycube_vertex_flag[polycube_edges[i][1]] = true;
    }
    std::vector<int> vert_o2n(dpx.size(), -1);
    int count = 0;
    for (size_t i = 0; i < dpx.size(); i++)
    {
        if (polycube_vertex_flag[i])
            vert_o2n[i] = count++;
    }
    std::ofstream outputfile(filename);
    outputfile << "# vtk DataFile Version 3.0\n"
        << "mesh vtk data\n"
        << "ASCII\n"
        << "DATASET POLYDATA\n";
    outputfile << "POINTS " << count << " double" << std::endl;
    for (size_t i = 0; i < dpx.size(); i++)
    {
        if (polycube_vertex_flag[i])
            outputfile << dpx[i] << " " << dpy[i] << " " << dpz[i] << std::endl;
    }
    std::vector<int> feature_point_color(count, 0);
    outputfile << "LINES " << feature_polycube_edge.size() << " " << 3 * feature_polycube_edge.size() << std::endl;
    for (size_t i = 0; i < feature_polycube_edge.size(); i++)
    {
        int id0 = vert_o2n[polycube_edges[feature_polycube_edge[i]][0]];
        int id1 = vert_o2n[polycube_edges[feature_polycube_edge[i]][1]];
        feature_point_color[id0] = 1;
        feature_point_color[id1] = 1;
        outputfile << "2 " << id0 << " " << id1 << std::endl;
    }
    outputfile << "POINT_DATA " << count << "\n"
        << "SCALARS V_Scalars int\nLOOKUP_TABLE V_Table" << std::endl;
    for (size_t i = 0; i < count; i++)
    {
        outputfile << feature_point_color[i] << std::endl;
    }
    outputfile.close();
}
void polycube_flattening_interface::save_feature_long_edges_vtk_psfeformat(const char* filename)
{
    std::vector<bool> polycube_vertex_flag(dpx.size(), false);
    for (size_t i = 0; i < polycube_edges.size(); i++)
    {
        polycube_vertex_flag[polycube_edges[i][0]] = true;
        polycube_vertex_flag[polycube_edges[i][1]] = true;
    }
    std::vector<int> vert_o2n(dpx.size(), -1);
    int count = 0;
    for (size_t i = 0; i < dpx.size(); i++)
    {
        if (polycube_vertex_flag[i])
            vert_o2n[i] = count++;
    }
    std::ofstream outputfile(filename);
    outputfile << "POINTS " << count << std::endl;
    for (size_t i = 0; i < dpx.size(); i++)
    {
        if (polycube_vertex_flag[i])
            outputfile << dpx[i] << " " << dpy[i] << " " << dpz[i] << std::endl;
    }
    assert(feature_polycube_edge_segm.size() == feature_polycube_edge.size());
    outputfile << "EDGES " << feature_polycube_edge.size() << " " << std::endl;
    for (size_t i = 0; i < feature_polycube_edge.size(); i++)
    {
        int id0 = vert_o2n[polycube_edges[feature_polycube_edge[i]][0]];
        int id1 = vert_o2n[polycube_edges[feature_polycube_edge[i]][1]];
        outputfile << id0 << " " << id1 << " ";
        outputfile << feature_polycube_edge_segm[i].first << " " << feature_polycube_edge_segm[i].second << std::endl;
    }
    outputfile.close();
}
void polycube_flattening_interface::prepare_for_deformation(TetStructure<double> *tet_mesh_, const std::vector<double> &x_ori, const std::vector<double> &y_ori, const std::vector<double> &z_ori)
{
    std::vector<OpenVolumeMesh::VertexHandle> hfv_vec(3);
    std::vector<int> hfv_vec_id(3);
    std::vector<OpenVolumeMesh::HalfFaceHandle> cell_hfh_vec(4);
    unsigned nv = tet_mesh_->tetra_vertices.size(); unsigned nc = tet_mesh_->tetras.size();
    vertex_cell.clear(); cell_vertex.clear();
    vertex_cell.resize(nv); cell_vertex.resize(nc);
    cell_vertex_vertex.clear(); S.clear(); S_v.clear(); cell_S.resize(nc);
    cell_vertex_vertex.resize(nc); S.resize(nc); S_v.resize(nc);
    vertex_cell_vertex.clear(); vertex_cell_vertex.resize(nv);
    vcv_S.clear(); vcv_S.resize(nv);
    cell_volume.clear(); cell_volume.resize(nc, 0.0);
    Eigen::Matrix3d VS, IS; std::vector<double> S_(9);
    vertex_move_d.clear(); vertex_move_d.resize(nv, 1.0e30);
    static int s_tet_id[4][3] = { { 1, 2, 3}, {2, 0, 3}, {0, 1, 3}, {1, 0, 2} };
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
        int id0, id1, id2, id3;
        id0 = tet_vert_[0]->id;
        id1 = tet_vert_[1]->id;
        id2 = tet_vert_[2]->id;
        id3 = tet_vert_[3]->id;
        CVec<double, 3> cp;
        CVec<double, 3> cr;
        CVec<double, 3> cs;
        CVec<double, 3> ct;
        cp[0] = x_ori[id0]; cp[1] = y_ori[id0]; cp[2] = z_ori[id0];
        cr[0] = x_ori[id1]; cr[1] = y_ori[id1]; cr[2] = z_ori[id1];
        cs[0] = x_ori[id2]; cs[1] = y_ori[id2]; cs[2] = z_ori[id2];
        ct[0] = x_ori[id3]; ct[1] = y_ori[id3]; ct[2] = z_ori[id3];
        cr = cr - cp;
        cs = cs - cp;
        ct = ct - cp;
        VS(0, 0) = cr[0]; VS(1, 0) = cr[1]; VS(2, 0) = cr[2];
        VS(0, 1) = cs[0]; VS(1, 1) = cs[1]; VS(2, 1) = cs[2];
        VS(0, 2) = ct[0]; VS(1, 2) = ct[1]; VS(2, 2) = ct[2];
        IS = VS.inverse();
        cell_S[c_id] = IS;
        S[c_id].resize(4); S_v[c_id].resize(4);
        for (unsigned i = 0; i < 4; ++i)
        {
            S_v[c_id][i].reserve(9);
            int v_id = tet_vert_[i]->id;
            if (vertex_cell[v_id].size() == 0) vertex_cell[v_id].reserve(15);
            vertex_cell[v_id].push_back(c_id);
            hfv_vec_id[0] = tetras[c_id]->vertex[s_tet_id[i][0]]->id;
            hfv_vec_id[1] = tetras[c_id]->vertex[s_tet_id[i][1]]->id;
            hfv_vec_id[2] = tetras[c_id]->vertex[s_tet_id[i][2]]->id;
            cell_vertex_vertex[c_id].push_back(hfv_vec_id);
            vertex_cell_vertex[v_id].push_back(hfv_vec_id);
            CVec<double, 3> p;
            CVec<double, 3> r;
            CVec<double, 3> s;
            CVec<double, 3> t;
            id0 = tetras[c_id]->vertex[i]->id;
            id1 = hfv_vec_id[0];
            id2 = hfv_vec_id[1];
            id3 = hfv_vec_id[2];
            p[0] = x_ori[id0]; p[1] = y_ori[id0]; p[2] = z_ori[id0];
            r[0] = x_ori[id1]; r[1] = y_ori[id1]; r[2] = z_ori[id1];
            s[0] = x_ori[id2]; s[1] = y_ori[id2]; s[2] = z_ori[id2];
            t[0] = x_ori[id3]; t[1] = y_ori[id3]; t[2] = z_ori[id3];
            r = r - p;
            s = s - p;
            t = t - p;
            VS(0, 0) = r[0]; VS(1, 0) = r[1]; VS(2, 0) = r[2];
            VS(0, 1) = s[0]; VS(1, 1) = s[1]; VS(2, 1) = s[2];
            VS(0, 2) = t[0]; VS(1, 2) = t[1]; VS(2, 2) = t[2];
            cell_volume[c_id] = -VS.determinant();
            IS = VS.inverse();
            S_[0] = IS(0, 0); S_[1] = IS(0, 1); S_[2] = IS(0, 2);
            S_[3] = IS(1, 0); S_[4] = IS(1, 1); S_[5] = IS(1, 2);
            S_[6] = IS(2, 0); S_[7] = IS(2, 1); S_[8] = IS(2, 2);
            S[c_id][i] = IS;
            S_v[c_id][i] = S_;
            vcv_S[v_id].push_back(S_);
        }
    }
    alpha = 100.0; mu_A = 1e3;
    is_deformation_constraint.clear(); is_deformation_constraint.resize(nv, 0);
    deformation_handle_id.clear(); deformation_handle_id.resize(nv, -1);
    for (unsigned i = 0; i < deformation_v_id.size(); ++i)
    {
        is_deformation_constraint[deformation_v_id[i]] = 1;
        deformation_handle_id[deformation_v_id[i]] = i;
    }
    dpx.resize(nv); dpy.resize(nv); dpz.resize(nv); max_vc_size = 0; src_pos.resize(nv);
    for (int i = 0; i < nv; i++)
    {
        int v_id = i;
        CVec<double, 3> p = tetra_vertices[i]->pos;
        dpx[v_id] = p[0]; dpy[v_id] = p[1]; dpz[v_id] = p[2];
        OpenVolumeMesh::Geometry::Vec3d p_temp;
        p_temp[0] = x_ori[i]; p_temp[1] = y_ori[i]; p_temp[2] = z_ori[i];
        src_pos[v_id] = p_temp;
        int vc_size = vertex_cell[v_id].size();
        if (vc_size > max_vc_size)
        {
            max_vc_size = vc_size;
        }
    }
    bfv_id.clear(); bfv_id.resize(nc, OpenVolumeMesh::Geometry::Vec3i(-1, -1, -1));
    bvf_id.clear(); bvf_id.resize(nv); boundary_face_number = 0; avg_boundary_edge_length = 0.0;
    int edge_count = 0;
    edge2id.clear();
    id2edge.clear();
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
        std::pair<int, int> edge0(idx0, idx1), edge1(idx1, idx2), edge2(idx2, idx0);
        if (idx0 < idx1)
        {
            edge2id.insert(std::pair<std::pair<int, int>, int>(edge0, edge_count));
            id2edge.push_back(edge0);
            edge_count++;
        }
        if (idx1 < idx2)
        {
            edge2id.insert(std::pair<std::pair<int, int>, int>(edge1, edge_count));
            id2edge.push_back(edge1);
            edge_count++;
        }
        if (idx2 < idx0)
        {
            edge2id.insert(std::pair<std::pair<int, int>, int>(edge2, edge_count));
            id2edge.push_back(edge2);
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
    avg_boundary_edge_length = 0.0; OpenVolumeMesh::Geometry::Vec3d p0, p1, p2; double count = 0.0;
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
        count += 3.0;
    }
    avg_boundary_edge_length /= count;
    rest_v_id.clear(); rotate_v_id.clear(); rotate_radius.clear(); rotate_angle.clear();
#if 0
    deformation_v_id.clear(); deformation_new_p.clear(); deformation_ok.clear();
    is_deformation_constraint.clear(); is_deformation_constraint.resize(nv, 0);
    for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
    {
        int v_id = v_it->idx();
        OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(*v_it);
#if 1
        if (std::abs(p[2] + 19.5) < 1e-6)
        {
            rest_v_id.push_back(v_id);
            deformation_v_id.push_back(v_id); deformation_new_p.push_back(p); deformation_ok.push_back(1);
            is_deformation_constraint[v_id] = 1;
        }
#else
        if (std::abs(p[2] + 19.5) < 1e-6 && std::abs(std::abs(p[1]) - 5) < 1e-6 && std::abs(std::abs(p[0]) - 5) < 1e-6)
        {
            rest_v_id.push_back(v_id);
            deformation_v_id.push_back(v_id); deformation_new_p.push_back(p); deformation_ok.push_back(1);
            is_deformation_constraint[v_id] = 1;
        }
        else if (std::abs(p[2] + 19.5) < 1e-6 && std::abs(p[1]) < 1e-6 && std::abs(p[0]) < 1e-6)
        {
            rest_v_id.push_back(v_id);
            deformation_v_id.push_back(v_id); deformation_new_p.push_back(p); deformation_ok.push_back(1);
            is_deformation_constraint[v_id] = 1;
        }
#endif
        else if (std::abs(p[2] - 19.5) < 1e-6)
        {
            if (std::abs(p[1]) < 1e-6 && std::abs(p[0]) < 1e-6)
            {
                rest_v_id.push_back(v_id);
                deformation_v_id.push_back(v_id); deformation_new_p.push_back(p); deformation_ok.push_back(1);
            }
#if 1
            else
#else
            else if (std::abs(std::abs(p[1]) - 5) < 1e-6 && std::abs(std::abs(p[0]) - 5) < 1e-6)
#endif
            {
                rotate_v_id.push_back(v_id);
                rotate_radius.push_back(std::sqrt(p[1] * p[1] + p[0] * p[0]));
                rotate_angle.push_back(std::atan2(p[1], p[0]));
                deformation_v_id.push_back(v_id); deformation_new_p.push_back(p); deformation_ok.push_back(1);
            }
            is_deformation_constraint[v_id] = 1;
        }
    }
#endif
    prepare_ok = true;
}
void polycube_flattening_interface::build_AABB_Tree()
{
    if (AABB_Tree) { delete AABB_Tree; AABB_Tree = NULL; }
    unsigned nf = bfv_id.size();
    if (nf == 0) return;
    triangle_vectors.clear(); triangle_vectors.reserve(nf);
    OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
    for (int i = 0; i < nf; ++i)
    {
        OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[i];
        if (fv_id[0] >= 0)
        {
            p0 = OpenVolumeMesh::Geometry::Vec3d(dpx[fv_id[0]], dpy[fv_id[0]], dpz[fv_id[0]]);
            p1 = OpenVolumeMesh::Geometry::Vec3d(dpx[fv_id[1]], dpy[fv_id[1]], dpz[fv_id[1]]);
            p2 = OpenVolumeMesh::Geometry::Vec3d(dpx[fv_id[2]], dpy[fv_id[2]], dpz[fv_id[2]]);
            CGAL_double_3_Point q0(p0[0], p0[1], p0[2]);
            CGAL_double_3_Point q1(p1[0], p1[1], p1[2]);
            CGAL_double_3_Point q2(p2[0], p2[1], p2[2]);
            triangle_vectors.push_back(CGAL_3_Triangle(q0, q1, q2));
        }
    }
    AABB_Tree = new CGAL_AABB_Tree(triangle_vectors.begin(), triangle_vectors.end());
    AABB_Tree->accelerate_distance_queries();
    printf("Finish Constructing AABB Tree.\n");
    tan_smooth = true;
}
void polycube_flattening_interface::project_on_ref_mesh(OpenVolumeMesh::Geometry::Vec3d& p)
{
    CGAL_double_3_Point pos = AABB_Tree->closest_point(CGAL_double_3_Point(p[0], p[1], p[2]));
    p = OpenVolumeMesh::Geometry::Vec3d(pos.x(), pos.y(), pos.z());
}
void polycube_flattening_interface::project_on_ref_mesh(double& px, double& py, double& pz)
{
    CGAL_double_3_Point pos = AABB_Tree->closest_point(CGAL_double_3_Point(px, py, pz));
    px = pos.x(); py = pos.y(); pz = pos.z();
}
void polycube_flattening_interface::assign_pos_mesh(TetStructure<double>* tet_mesh_, bool r_order)
{
    if (r_order)
    {
        int nv = tet_mesh_->tetra_vertices.size();
        const std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
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
void polycube_flattening_interface::get_coord(const std::vector<int>& all_vert, std::vector<double> &coord_x, std::vector<double> &coord_y, std::vector<double> &coord_z)
{
    coord_x.clear();
    coord_y.clear();
    coord_z.clear();
    for (size_t i = 0; i < all_vert.size(); i++)
    {
        int v_id = all_vert[i];
        coord_x.push_back(dpx[v_id]);
        coord_y.push_back(dpy[v_id]);
        coord_z.push_back(dpz[v_id]);
    }
}
ig::CVec<double, 3> polycube_flattening_interface::compute_volumetric_distortion(TetStructure<double>* tet_mesh_)
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
    printf("Volumetric Weighted Isometric Distortion: %f/%f/%f/%f\n", max_iso, min_iso, avg_iso, std_id);
    printf("Volumetric Weighted conformal Distortion: %f/%f/%f/%f\n", max_cd, min_cd, avg_cd, std_cd);
    printf("Volumetric Weighted volume Distortion: %f/%f/%f/%f\n", max_vol, min_vol, avg_vol, std_vd);
    return ig::CVec<double, 3>(avg_iso, avg_cd, avg_vol);
}
ig::CVec<double, 3> polycube_flattening_interface::compute_distortion(TetStructure<double>* tet_mesh_)
{
    Eigen::Matrix3d VS, A; Eigen::Vector3d a;
    double max_cd = 0.0; double min_cd = 1e30; double avg_cd = 0.0; int max_cd_cell_id = -1;
    double max_iso = 0.0; double min_iso = 1e30; double avg_iso = 0.0; int max_iso_cell_id = -1;
    double max_vol = 0.0; double min_vol = 1e30; double avg_vol = 0.0; int max_vol_cell_id = -1;
    double cd_count = 0.0; int flip_count = 0; flipped_cell.clear();
    int nc = tet_mesh_->tetras.size();
    all_iso_d.clear(); all_con_d.clear(); all_vol_d.clear();
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
    return ig::CVec<double, 3>(avg_iso, avg_cd, avg_vol);
}
void polycube_flattening_interface::get_chart_distortion()
{
    chart_iso_d.clear();
    chart_iso_d.resize(polycube_chart.size());
    chart_con_d.clear();
    chart_con_d.resize(polycube_chart.size());
    chart_vol_d.clear();
    chart_vol_d.resize(polycube_chart.size());
    if (all_iso_d.size() == 0)
    {
        std::cout << "distortion not computed yet!" << std::endl;
        return;
    }
    for (size_t i = 0; i < polycube_chart.size(); i++)
    {
        double count = 0.0;
        double all_distortion = 0.0;
        double con_distortion = 0.0;
        double vol_distortion = 0.0;
        for (size_t j = 0; j < polycube_chart[i].size(); j++)
        {
            if (all_iso_d[polycube_chart[i][j]] > 0)
            {
                all_distortion += all_iso_d[polycube_chart[i][j]];
                con_distortion += all_con_d[polycube_chart[i][j]];
                vol_distortion += all_vol_d[polycube_chart[i][j]];
                count = count + 1.0;
            }
        }
        chart_iso_d[i] = all_distortion / count;
        chart_con_d[i] = con_distortion / count;
        chart_vol_d[i] = vol_distortion / count;
    }
}
void polycube_flattening_interface::get_polycube_edge_layer_distortion(int n_layer, TetStructure<double> *tet_mesh_)
{
    if (edge_with_same_label.size() == 0)
        return;
    std::vector<double> tmp_distortion(edge_with_same_label.size(), 0);
    polycube_edge_layer_distortion.clear();
    polycube_edge_layer_distortion.resize(polycube_edges.size());
    const std::vector<Tetrahedron<double>*> &tetras = tet_mesh_->tetras;
    const std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
    std::map<Tetrahedron<double>*, int> tet2id;
    int nc = tet_mesh_->tetras.size();
    for (size_t i = 0; i < nc; i++)
    {
        tet2id[tetras[i]] = i;
    }
    for (size_t i = 0; i < edge_with_same_label.size(); i++)
    {
        std::vector<int> edge_layer_flag;
        edge_layer_flag.resize(nc, -1);
        std::vector<int> edge_layer_idx;
        for (int j = 0; j < edge_with_same_label[i].size(); ++j)
        {
            int edge_id = edge_with_same_label[i][j];
            int face1 = bef_id[edge_id][0];
            int face2 = bef_id[edge_id][1];
            edge_layer_flag[face1] = 1;
            edge_layer_flag[face2] = 1;
            edge_layer_idx.push_back(face1);
            edge_layer_idx.push_back(face2);
        }
        for (size_t iter = 0; iter < n_layer; iter++)
        {
            std::vector<int> added_idx;
            for (size_t j = 0; j < edge_layer_idx.size(); j++)
            {
                const Tetrahedron<double>* temp_tet = tetras[edge_layer_idx[j]];
                for (size_t k = 0; k < 4; k++)
                {
                    Tetrahedron<double>* neighbor = temp_tet->neighborTet[k];
                    if (neighbor != NULL)
                    {
                        int idx = tet2id[neighbor];
                        if (edge_layer_flag[idx] < 0)
                        {
                            added_idx.push_back(idx);
                        }
                    }
                }
            }
            edge_layer_idx.insert(edge_layer_idx.begin(), added_idx.begin(), added_idx.end());
            for (size_t j = 0; j < added_idx.size(); j++)
            {
                edge_layer_flag[added_idx[j]] = 1;
            }
        }
        double chart_distortion = 0.0;
        double count = 0.0;
        for (size_t j = 0; j < edge_layer_idx.size(); j++)
        {
            int idx = edge_layer_idx[j];
            if (all_iso_d[idx] > 0)
            {
                chart_distortion += all_iso_d[idx];
                count = count + 1;
            }
        }
        tmp_distortion[i] = chart_distortion / count;
    }
    for (size_t i = 0; i < polycube_edges.size(); i++)
    {
        int tmp_idx = polycube_edge_idx2same_label_idx[i];
        polycube_edge_layer_distortion[i] = tmp_distortion[tmp_idx];
    }
}
void polycube_flattening_interface::compute_triangle_area()
{
    double max_diff = 0.0;
    double ave_diff = 0.0;
    for (size_t i = 0; i < equal_triangles.size(); i++)
    {
        double x[6], y[6], z[6];
        for (size_t j = 0; j < 6; j++)
        {
            x[j] = dpx[equal_triangles[i][j]];
            y[j] = dpy[equal_triangles[i][j]];
            z[j] = dpz[equal_triangles[i][j]];
        }
        Eigen::Vector3d vec1(x[1] - x[0], y[1] - y[0], z[1] - z[0]);
        Eigen::Vector3d vec2(x[2] - x[0], y[2] - y[0], z[2] - z[0]);
        double area1 = vec1.cross(vec2).norm();
        Eigen::Vector3d vec3(x[4] - x[3], y[4] - y[3], z[4] - z[3]);
        Eigen::Vector3d vec4(x[5] - x[3], y[5] - y[3], z[5] - z[3]);
        double area2 = vec3.cross(vec4).norm();
        double diff = abs(area1 - area2);
        if (max_diff < diff)
            max_diff = diff;
        ave_diff += diff;
    }
    ave_diff = ave_diff / (equal_triangles.size() + 1e-8);
}

void polycube_flattening_interface::find_all_corner(TetStructure<double>* tet_mesh_)
{
    int chart_size = polycube_chart_label.size();
    std::vector<int> with_big_diff_chart; std::vector<double> max_chart; std::vector<int> min_chart;
    const std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
    for (int i = 0; i < chart_size; ++i)
    {
        int chart_label = polycube_chart_label[i] / 2; double sum = 0.0;
        std::vector<int>& one_chart = polycube_chart[i]; double max_value = -1e30; double min_value = 1e30;
        for (int j = 0; j < one_chart.size(); ++j)
        {
            int bf_id = one_chart[j]; double one_face_sum = 0.0;
            for (int k = 0; k < 3; ++k)
            {
                one_face_sum += tetra_vertices[bfv_id[bf_id][k]]->pos[chart_label];
            }
            double xyz_value = one_face_sum / 3.0;
            sum += xyz_value;
            if (max_value < xyz_value) max_value = xyz_value;
            if (min_value > xyz_value) min_value = xyz_value;
        }
        double max_diff = max_value - min_value;
        chart_mean_value[i] = 0.5*(max_value + min_value);
        if (max_diff > 0.001)
        {
            with_big_diff_chart.push_back(i); max_chart.push_back(max_value); min_chart.push_back(min_value);
        }
    }
}
bool polycube_flattening_interface::find_all_chart_value_equal_face_mine(TetStructure<double> *tet_mesh_, double offset)
{
    find_all_corner(tet_mesh_);
    std::vector<MSKboundkeye> bkc; std::vector<double> blc; std::vector<double> buc;
    std::vector<MSKboundkeye> bkx; std::vector<double> blx; std::vector<double> bux;
    std::vector<MSKlidxt> aptrb; std::vector<MSKidxt> asub; std::vector<double> aval;
    std::vector<MSKidxt> qsubi; std::vector<MSKidxt> qsubj; std::vector<double> qval;
    double ratio_diff = 0.8;
    std::vector< Eigen::Triplet<double> > coef;
    int up_down_chart_count = up_chart_id.size();
    int chart_number = chart_mean_value.size();
    double diff_u_d; int up_chart, down_chart;
    for (int i = 0; i < up_down_chart_count; ++i)
    {
        up_chart = up_chart_id[i];
        down_chart = down_chart_id[i];
        diff_u_d = diff_up_down[i] * ratio_diff;
        if (hex_meshing_flag)
        {
            if (diff_u_d < cube_len)
            {
                diff_u_d = cube_len;
            }
            coef.push_back(Eigen::Triplet<double>(i, up_chart, +1.0));
            coef.push_back(Eigen::Triplet<double>(i, down_chart, -1.0));
            bkc.push_back(MSK_BK_LO); blc.push_back(diff_u_d / cube_len); buc.push_back(+MSK_INFINITY);
        }
        else
        {
            coef.push_back(Eigen::Triplet<double>(i, up_chart, +1.0));
            coef.push_back(Eigen::Triplet<double>(i, down_chart, -1.0));
            bkc.push_back(MSK_BK_LO); blc.push_back(diff_u_d); buc.push_back(+MSK_INFINITY);
        }
    }
    int cut_num = 0;
    cut_num = cut_to_chart_pair.size();
    for (size_t i = 0; i < cut_num; i++)
    {
        int chart1 = cut_to_chart_pair[i].first;
        int chart2 = cut_to_chart_pair[i].second;
        if (cut_to_chart_pair_neighbor.size() == 0)
        {
            std::vector<int> chart1_idx, chart2_idx;
            for (size_t j = 0; j < polycube_edges.size(); j++)
            {
                if ((polycube_edges[j][2] == chart1 && polycube_edges[j][3] != chart2)
                    || (polycube_edges[j][3] == chart1 && polycube_edges[j][2] != chart2))
                    chart1_idx.push_back(j);
                if ((polycube_edges[j][2] == chart2 && polycube_edges[j][3] != chart1)
                    || (polycube_edges[j][3] == chart2 && polycube_edges[j][2] != chart1))
                    chart2_idx.push_back(j);
            }
            int final_idx1(-1), final_idx2(-1);
            for (size_t j = 0; j < chart1_idx.size(); j++)
            {
                if (up_chart_id[chart1_idx[j]] == chart2 || down_chart_id[chart1_idx[j]] == chart2)
                {
                    final_idx1 = chart1_idx[j];
                    break;
                }
            }
            for (size_t j = 0; j < chart2_idx.size(); j++)
            {
                if (up_chart_id[chart2_idx[j]] == chart1 || down_chart_id[chart2_idx[j]] == chart1)
                {
                    final_idx2 = chart2_idx[j];
                    break;
                }
            }
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, up_chart_id[final_idx1], +1.0));
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, down_chart_id[final_idx1], -1.0));
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, up_chart_id[final_idx2], -1.0));
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, down_chart_id[final_idx2], +1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
        }
        else
        {
            int up1, up2, down1, down2;
            up1 = cut_to_chart_pair_neighbor[i].first;
            down1 = cut_to_chart_pair[i].second;
            up2 = cut_to_chart_pair[i].first;
            down2 = cut_to_chart_pair_neighbor[i].second;
            if (chart_mean_value[up1] < chart_mean_value[down1])
            {
                int temp = up1;
                up1 = down1;
                down1 = temp;
            }
            if (chart_mean_value[up2] < chart_mean_value[down2])
            {
                int temp = up2;
                up2 = down2;
                down2 = temp;
            }
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, up1, +1.0));
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, down1, -1.0));
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, up2, -1.0));
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, down2, +1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
        }
    }
    int three_cut_num = 0;
    three_cut_num = three_cut_adjacent_one_cut_index.size();
    for (size_t i = 0; i < three_cut_num; i++)
    {
        std::vector<int> six_chart, chart_x, chart_y, chart_z;
        for (size_t j = 0; j < 3; j++)
        {
            six_chart.push_back(cut_to_chart_pair[three_cut_adjacent_one_cut_index[i][j]].first);
            six_chart.push_back(cut_to_chart_pair[three_cut_adjacent_one_cut_index[i][j]].second);
        }
        for (size_t j = 0; j < 6; j++)
        {
            int temp_label = polycube_chart_label[six_chart[j]] / 2;
            if (temp_label == 0)
                chart_x.push_back(six_chart[j]);
            else if (temp_label == 1)
                chart_y.push_back(six_chart[j]);
            else if (temp_label == 2)
                chart_z.push_back(six_chart[j]);
        }
        assert(chart_x.size() == 2);
        assert(chart_y.size() == 2);
        assert(chart_z.size() == 2);
        coef.push_back(Eigen::Triplet<double>(3 * i + up_down_chart_count + cut_num, chart_x[0], +1.0));
        coef.push_back(Eigen::Triplet<double>(3 * i + up_down_chart_count + cut_num, chart_x[1], -1.0));
        bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
        coef.push_back(Eigen::Triplet<double>(3 * i + 1 + up_down_chart_count + cut_num, chart_y[0], +1.0));
        coef.push_back(Eigen::Triplet<double>(3 * i + 1 + up_down_chart_count + cut_num, chart_y[1], -1.0));
        bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
        coef.push_back(Eigen::Triplet<double>(3 * i + 2 + up_down_chart_count + cut_num, chart_z[0], +1.0));
        coef.push_back(Eigen::Triplet<double>(3 * i + 2 + up_down_chart_count + cut_num, chart_z[1], -1.0));
        bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
    }
    int n_valence_six_constraint = 0;
    std::map<int, std::set<int>> vert2chart;
    for (size_t i = 0; i < polycube_edges.size(); i++)
    {
        int id0 = polycube_edges[i][0];
        int id1 = polycube_edges[i][1];
        int chart0 = polycube_edges[i][2];
        int chart1 = polycube_edges[i][3];
        vert2chart[id0].insert(chart0);
        vert2chart[id0].insert(chart1);
        vert2chart[id1].insert(chart0);
        vert2chart[id1].insert(chart1);
    }
    for (auto vc : vert2chart)
    {
        if (vc.second.size() == 6)
        {
            std::vector<int> xchart, ychart, zchart;
            for (auto chart : vc.second)
            {
                int axis = polycube_chart_label[chart] / 2;
                switch (axis)
                {
                case 0:
                    xchart.push_back(chart);
                    break;
                case 1:
                    ychart.push_back(chart);
                    break;
                case 2:
                    zchart.push_back(chart);
                    break;
                default:
                    break;
                }
            }
            assert(xchart.size() == 2 && ychart.size() == 2 && zchart.size() == 2);
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint, xchart[0], +1.0));
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint, xchart[1], -1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 1, ychart[0], +1.0));
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 1, ychart[1], -1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 2, zchart[0], +1.0));
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 2, zchart[1], -1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
            n_valence_six_constraint += 3;
        }
    }
    Eigen::SparseMatrix<double> A(3 * three_cut_num + up_down_chart_count + cut_num + n_valence_six_constraint, chart_number);
    A.setFromTriplets(coef.begin(), coef.end());
    int size = A.nonZeros(); int cols = A.cols();
    aptrb.resize(cols + 1); asub.resize(size); aval.resize(size);
    aptrb[cols] = size;
    int ind = 0;
    for (int k = 0; k < cols; ++k)
    {
        aptrb[k] = ind;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
        {
            aval[ind] = it.value();
            asub[ind] = it.row();
            ++ind;
        }
    }
    std::vector<double> c(chart_number, 0.0);
    for (int i = 0; i < chart_number; ++i)
    {
        qsubi.push_back(i); qsubj.push_back(i); qval.push_back(2.0);
        if (hex_meshing_flag)
        {
            c[i] = -2.0 * chart_mean_value[i] / cube_len;
        }
        else
        {
            c[i] = -2.0 * chart_mean_value[i];
        }
    }
    bkx.resize(chart_number, MSK_BK_FR); blx.resize(chart_number, -MSK_INFINITY); bux.resize(chart_number, +MSK_INFINITY);
    std::vector<double> x(chart_number, 0.0); int mosek_count = 0;
    bool solve_success;
    if (hex_meshing_flag)
    {
        solve_success = solveConvexQuadPorgramming_mosek_integer(bkc, blc, buc, bkx, blx, bux, aptrb, asub, aval, qsubi, qsubj, qval, c, x);
    }
    else
    {
        solve_success = solveConvexQuadPorgramming_mosek(bkc, blc, buc, bkx, blx, bux, aptrb, asub, aval, qsubi, qsubj, qval, c, x);
    }
    while (!solve_success)
    {
        ratio_diff *= 0.8;
        for (int i = 0; i < up_down_chart_count; ++i)
        {
            diff_u_d = diff_up_down[i] * ratio_diff;
            if (hex_meshing_flag)
            {
                if (diff_u_d < cube_len)
                {
                    diff_u_d = cube_len;
                }
                blc[i] = diff_u_d / cube_len;
            }
            else
            {
                blc[i] = diff_u_d;
            }
        }
        if (hex_meshing_flag)
        {
            solve_success = solveConvexQuadPorgramming_mosek_integer(bkc, blc, buc, bkx, blx, bux, aptrb, asub, aval, qsubi, qsubj, qval, c, x);
        }
        else
        {
            solve_success = solveConvexQuadPorgramming_mosek(bkc, blc, buc, bkx, blx, bux, aptrb, asub, aval, qsubi, qsubj, qval, c, x);
        }
        ++mosek_count;
        if (mosek_count > 20) break;
    }
    if (mosek_count > 20)
    {
        printf("-------------------------------------------\nNo Solution!!!!!\n");
        return false;
    }
    else
    {
        for (int i = 0; i < chart_number; ++i)
        {
            if (hex_meshing_flag)
            {
                chart_mean_value[i] = cube_len * x[i];
                double offset_ratio = offset;
                if (polycube_chart_label[i] % 2 == 0)
                {
                    chart_mean_value[i] += cube_len * offset_ratio;
                }
                else
                {
                    chart_mean_value[i] -= cube_len * offset_ratio;
                }
            }
            else
            {
                chart_mean_value[i] = x[i];
            }
        }
        if (hex_meshing_flag)
        {
            for (size_t i = 0; i < cut_num; i++)
            {
                int chart1 = cut_to_chart_pair[i].first;
                int chart2 = cut_to_chart_pair[i].second;
                chart_mean_value[chart1] = cube_len * x[chart1];
                chart_mean_value[chart2] = cube_len * x[chart2];
            }
        }
    }
    return true;
}
bool polycube_flattening_interface::find_all_chart_value_equal_face_mine(TetStructure<double> *tet_mesh_, const std::vector<double>& chart_mean_value_ori, double lambda, bool add_constraint, bool set_min_diff, double offset)
{
    find_all_corner(tet_mesh_);
    if (add_constraint)
    {
        update_updownchart(tet_mesh_);
    }
    std::vector<MSKboundkeye> bkc; std::vector<double> blc; std::vector<double> buc;
    std::vector<MSKboundkeye> bkx; std::vector<double> blx; std::vector<double> bux;
    std::vector<MSKlidxt> aptrb; std::vector<MSKidxt> asub; std::vector<double> aval;
    std::vector<MSKidxt> qsubi; std::vector<MSKidxt> qsubj; std::vector<double> qval;
    double ratio_diff = 0.8;
    std::vector< Eigen::Triplet<double> > coef;
    int up_down_chart_count = up_chart_id.size();
    int chart_number = chart_mean_value.size();
    double diff_u_d; int up_chart, down_chart;
    for (int i = 0; i < up_down_chart_count; ++i)
    {
        up_chart = up_chart_id[i];
        down_chart = down_chart_id[i];
        diff_u_d = diff_up_down[i] * ratio_diff;
        if (set_min_diff)
        {
            if (diff_u_d < cube_len)
            {
                diff_u_d = cube_len;
            }
        }
        if (hex_meshing_flag)
        {
            if (diff_u_d < cube_len)
            {
                diff_u_d = cube_len;
            }
            coef.push_back(Eigen::Triplet<double>(i, up_chart, +1.0));
            coef.push_back(Eigen::Triplet<double>(i, down_chart, -1.0));
            bkc.push_back(MSK_BK_LO); blc.push_back(diff_u_d / cube_len); buc.push_back(+MSK_INFINITY);
        }
        else
        {
            coef.push_back(Eigen::Triplet<double>(i, up_chart, +1.0));
            coef.push_back(Eigen::Triplet<double>(i, down_chart, -1.0));
            bkc.push_back(MSK_BK_LO); blc.push_back(diff_u_d); buc.push_back(+MSK_INFINITY);
        }
    }
    int cut_num = 0;
    cut_num = cut_to_chart_pair.size();
    for (size_t i = 0; i < cut_num; i++)
    {
        int chart1 = cut_to_chart_pair[i].first;
        int chart2 = cut_to_chart_pair[i].second;
        if (cut_to_chart_pair_neighbor.size() == 0)
        {
            std::vector<int> chart1_idx, chart2_idx;
            for (size_t j = 0; j < polycube_edges.size(); j++)
            {
                if ((polycube_edges[j][2] == chart1 && polycube_edges[j][3] != chart2)
                    || (polycube_edges[j][3] == chart1 && polycube_edges[j][2] != chart2))
                    chart1_idx.push_back(j);
                if ((polycube_edges[j][2] == chart2 && polycube_edges[j][3] != chart1)
                    || (polycube_edges[j][3] == chart2 && polycube_edges[j][2] != chart1))
                    chart2_idx.push_back(j);
            }
            int final_idx1(-1), final_idx2(-1);
            for (size_t j = 0; j < chart1_idx.size(); j++)
            {
                if (up_chart_id[chart1_idx[j]] == chart2 || down_chart_id[chart1_idx[j]] == chart2)
                {
                    final_idx1 = chart1_idx[j];
                    break;
                }
            }
            for (size_t j = 0; j < chart2_idx.size(); j++)
            {
                if (up_chart_id[chart2_idx[j]] == chart1 || down_chart_id[chart2_idx[j]] == chart1)
                {
                    final_idx2 = chart2_idx[j];
                    break;
                }
            }
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, up_chart_id[final_idx1], +1.0));
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, down_chart_id[final_idx1], -1.0));
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, up_chart_id[final_idx2], -1.0));
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, down_chart_id[final_idx2], +1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
        }
        else
        {
            int up1, up2, down1, down2;
            up1 = cut_to_chart_pair_neighbor[i].first;
            down1 = cut_to_chart_pair[i].second;
            up2 = cut_to_chart_pair[i].first;
            down2 = cut_to_chart_pair_neighbor[i].second;
            if (chart_mean_value[up1] < chart_mean_value[down1])
            {
                int temp = up1;
                up1 = down1;
                down1 = temp;
            }
            if (chart_mean_value[up2] < chart_mean_value[down2])
            {
                int temp = up2;
                up2 = down2;
                down2 = temp;
            }
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, up1, +1.0));
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, down1, -1.0));
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, up2, -1.0));
            coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, down2, +1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
        }
    }
    int three_cut_num = 0;
    three_cut_num = three_cut_adjacent_one_cut_index.size();
    for (size_t i = 0; i < three_cut_num; i++)
    {
        std::vector<int> six_chart, chart_x, chart_y, chart_z;
        for (size_t j = 0; j < 3; j++)
        {
            six_chart.push_back(cut_to_chart_pair[three_cut_adjacent_one_cut_index[i][j]].first);
            six_chart.push_back(cut_to_chart_pair[three_cut_adjacent_one_cut_index[i][j]].second);
        }
        for (size_t j = 0; j < 6; j++)
        {
            int temp_label = polycube_chart_label[six_chart[j]] / 2;
            if (temp_label == 0)
                chart_x.push_back(six_chart[j]);
            else if (temp_label == 1)
                chart_y.push_back(six_chart[j]);
            else if (temp_label == 2)
                chart_z.push_back(six_chart[j]);
        }
        assert(chart_x.size() == 2);
        assert(chart_y.size() == 2);
        assert(chart_z.size() == 2);
        coef.push_back(Eigen::Triplet<double>(3 * i + up_down_chart_count + cut_num, chart_x[0], +1.0));
        coef.push_back(Eigen::Triplet<double>(3 * i + up_down_chart_count + cut_num, chart_x[1], -1.0));
        bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
        coef.push_back(Eigen::Triplet<double>(3 * i + 1 + up_down_chart_count + cut_num, chart_y[0], +1.0));
        coef.push_back(Eigen::Triplet<double>(3 * i + 1 + up_down_chart_count + cut_num, chart_y[1], -1.0));
        bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
        coef.push_back(Eigen::Triplet<double>(3 * i + 2 + up_down_chart_count + cut_num, chart_z[0], +1.0));
        coef.push_back(Eigen::Triplet<double>(3 * i + 2 + up_down_chart_count + cut_num, chart_z[1], -1.0));
        bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
    }
    int n_valence_six_constraint = 0;
    std::map<int, std::set<int>> vert2chart;
    for (size_t i = 0; i < polycube_edges.size(); i++)
    {
        int id0 = polycube_edges[i][0];
        int id1 = polycube_edges[i][1];
        int chart0 = polycube_edges[i][2];
        int chart1 = polycube_edges[i][3];
        vert2chart[id0].insert(chart0);
        vert2chart[id0].insert(chart1);
        vert2chart[id1].insert(chart0);
        vert2chart[id1].insert(chart1);
    }
    for (auto vc : vert2chart)
    {
        if (vc.second.size() == 6)
        {
            std::vector<int> xchart, ychart, zchart;
            for (auto chart : vc.second)
            {
                int axis = polycube_chart_label[chart] / 2;
                switch (axis)
                {
                case 0:
                    xchart.push_back(chart);
                    break;
                case 1:
                    ychart.push_back(chart);
                    break;
                case 2:
                    zchart.push_back(chart);
                    break;
                default:
                    break;
                }
            }
            assert(xchart.size() == 2 && ychart.size() == 2 && zchart.size() == 2);
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint, xchart[0], +1.0));
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint, xchart[1], -1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 1, ychart[0], +1.0));
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 1, ychart[1], -1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 2, zchart[0], +1.0));
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 2, zchart[1], -1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
            n_valence_six_constraint += 3;
        }
    }
    Eigen::SparseMatrix<double> A(3 * three_cut_num + up_down_chart_count + cut_num + n_valence_six_constraint, chart_number);
    A.setFromTriplets(coef.begin(), coef.end());
    int size = A.nonZeros(); int cols = A.cols();
    aptrb.resize(cols + 1); asub.resize(size); aval.resize(size);
    aptrb[cols] = size;
    int ind = 0;
    for (int k = 0; k < cols; ++k)
    {
        aptrb[k] = ind;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
        {
            aval[ind] = it.value();
            asub[ind] = it.row();
            ++ind;
        }
    }
    std::vector<double> c(chart_number, 0.0);
    for (int i = 0; i < chart_number; ++i)
    {
        qsubi.push_back(i); qsubj.push_back(i); qval.push_back(2.0);
        if (hex_meshing_flag)
        {
            c[i] = -2.0 * chart_mean_value[i] / cube_len;
        }
        else
        {
            c[i] = -2.0 * chart_mean_value[i];
        }
    }
    assert(chart_mean_value_ori.size() <= chart_mean_value.size());
    for (size_t i = 0; i < chart_mean_value_ori.size(); i++)
    {
        qval[i] = 2.0 + 2.0 * lambda;
        if (hex_meshing_flag)
        {
            c[i] += -2.0 * chart_mean_value_ori[i] * lambda / cube_len;
        }
        else
        {
            c[i] += -2.0 * chart_mean_value_ori[i] * lambda;
        }
    }
    bkx.resize(chart_number, MSK_BK_FR); blx.resize(chart_number, -MSK_INFINITY); bux.resize(chart_number, +MSK_INFINITY);
    std::vector<double> x(chart_number, 0.0); int mosek_count = 0;
    bool solve_success;
    if (hex_meshing_flag)
        solve_success = solveConvexQuadPorgramming_mosek_integer(bkc, blc, buc, bkx, blx, bux, aptrb, asub, aval, qsubi, qsubj, qval, c, x);
    else
        solve_success = solveConvexQuadPorgramming_mosek(bkc, blc, buc, bkx, blx, bux, aptrb, asub, aval, qsubi, qsubj, qval, c, x);
    while (!solve_success)
    {
        ratio_diff *= 0.8;
        for (int i = 0; i < up_down_chart_count; ++i)
        {
            diff_u_d = diff_up_down[i] * ratio_diff;
            if (set_min_diff)
            {
                if (diff_u_d < cube_len)
                {
                    diff_u_d = cube_len;
                }
            }
            if (hex_meshing_flag)
            {
                if (diff_u_d < cube_len)
                {
                    diff_u_d = cube_len;
                }
                blc[i] = diff_u_d / cube_len;
            }
            else
            {
                blc[i] = diff_u_d;
            }
        }
        if (hex_meshing_flag)
            solve_success = solveConvexQuadPorgramming_mosek_integer(bkc, blc, buc, bkx, blx, bux, aptrb, asub, aval, qsubi, qsubj, qval, c, x);
        else
            solve_success = solveConvexQuadPorgramming_mosek(bkc, blc, buc, bkx, blx, bux, aptrb, asub, aval, qsubi, qsubj, qval, c, x);
        ++mosek_count;
        if (mosek_count > 20) break;
    }
    if (mosek_count > 20)
    {
        printf("-------------------------------------------\nNo Solution!!!!!\n");
        return false;
    }
    else
    {
        for (int i = 0; i < chart_number; ++i)
        {
            if (hex_meshing_flag)
            {
                chart_mean_value[i] = cube_len * x[i];
                double offset_ratio = offset;
                if (polycube_chart_label[i] % 2 == 0)
                {
                    chart_mean_value[i] += cube_len * offset_ratio;
                }
                else
                {
                    chart_mean_value[i] -= cube_len * offset_ratio;
                }
            }
            else
            {
                chart_mean_value[i] = x[i];
            }
        }
        if (hex_meshing_flag)
            for (size_t i = 0; i < cut_num; i++)
            {
                int chart1 = cut_to_chart_pair[i].first;
                int chart2 = cut_to_chart_pair[i].second;
                chart_mean_value[chart1] = cube_len * x[chart1];
                chart_mean_value[chart2] = cube_len * x[chart2];
            }
    }
    return true;
}
bool polycube_flattening_interface::find_all_chart_value_preprocessing(TetStructure<double> *tet_mesh_)
{
    find_all_corner(tet_mesh_);
    std::vector<MSKboundkeye> bkc; std::vector<double> blc; std::vector<double> buc;
    std::vector<MSKboundkeye> bkx; std::vector<double> blx; std::vector<double> bux;
    std::vector<MSKlidxt> aptrb; std::vector<MSKidxt> asub; std::vector<double> aval;
    std::vector<MSKidxt> qsubi; std::vector<MSKidxt> qsubj; std::vector<double> qval;
    double ratio_diff = 0.8;
    std::vector< Eigen::Triplet<double> > coef;
    int up_down_chart_count = up_chart_id.size();
    int chart_number = chart_mean_value.size();
    double diff_u_d; int up_chart, down_chart;
    for (int i = 0; i < up_down_chart_count; ++i)
    {
        up_chart = up_chart_id[i];
        down_chart = down_chart_id[i];
        coef.push_back(Eigen::Triplet<double>(i, up_chart, +1.0));
        coef.push_back(Eigen::Triplet<double>(i, down_chart, -1.0));
        bkc.push_back(MSK_BK_LO); blc.push_back(cube_len); buc.push_back(+MSK_INFINITY);
    }
    int n_valence_six_constraint = 0;
    std::map<int, std::set<int>> vert2chart;
    for (size_t i = 0; i < polycube_edges.size(); i++)
    {
        int id0 = polycube_edges[i][0];
        int id1 = polycube_edges[i][1];
        int chart0 = polycube_edges[i][2];
        int chart1 = polycube_edges[i][3];
        vert2chart[id0].insert(chart0);
        vert2chart[id0].insert(chart1);
        vert2chart[id1].insert(chart0);
        vert2chart[id1].insert(chart1);
    }
    for (auto vc : vert2chart)
    {
        if (vc.second.size() == 6)
        {
            std::vector<int> xchart, ychart, zchart;
            for (auto chart : vc.second)
            {
                int axis = polycube_chart_label[chart] / 2;
                switch (axis)
                {
                case 0:
                    xchart.push_back(chart);
                    break;
                case 1:
                    ychart.push_back(chart);
                    break;
                case 2:
                    zchart.push_back(chart);
                    break;
                default:
                    break;
                }
            }
            assert(xchart.size() == 2 && ychart.size() == 2 && zchart.size() == 2);
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint, xchart[0], +1.0));
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint, xchart[1], -1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 1, ychart[0], +1.0));
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 1, ychart[1], -1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 2, zchart[0], +1.0));
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 2, zchart[1], -1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
            n_valence_six_constraint += 3;
        }
    }
    Eigen::SparseMatrix<double> A(up_down_chart_count + n_valence_six_constraint, chart_number);
    A.setFromTriplets(coef.begin(), coef.end());
    int size = A.nonZeros(); int cols = A.cols();
    aptrb.resize(cols + 1); asub.resize(size); aval.resize(size);
    aptrb[cols] = size;
    int ind = 0;
    for (int k = 0; k < cols; ++k)
    {
        aptrb[k] = ind;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
        {
            aval[ind] = it.value();
            asub[ind] = it.row();
            ++ind;
        }
    }
    std::vector<double> c(chart_number, 0.0);
    for (int i = 0; i < chart_number; ++i)
    {
        qsubi.push_back(i); qsubj.push_back(i); qval.push_back(2.0);
        c[i] = -2.0 * chart_mean_value[i];
    }
    bkx.resize(chart_number, MSK_BK_FR); blx.resize(chart_number, -MSK_INFINITY); bux.resize(chart_number, +MSK_INFINITY);
    std::vector<double> x(chart_number, 0.0); int mosek_count = 0;
    bool solve_success;
    solve_success = solveConvexQuadPorgramming_mosek(bkc, blc, buc, bkx, blx, bux, aptrb, asub, aval, qsubi, qsubj, qval, c, x);
    while (!solve_success)
    {
        ratio_diff *= 0.8;
        for (int i = 0; i < up_down_chart_count; ++i)
        {
            diff_u_d = diff_up_down[i] * ratio_diff;
            blc[i] = diff_u_d;
        }
        solve_success = solveConvexQuadPorgramming_mosek(bkc, blc, buc, bkx, blx, bux, aptrb, asub, aval, qsubi, qsubj, qval, c, x);
        ++mosek_count;
        if (mosek_count > 20) break;
    }
    if (mosek_count > 20)
    {
        printf("-------------------------------------------\nNo Solution!!!!!\n");
        return false;
    }
    else
    {
        for (int i = 0; i < chart_number; ++i)
        {
            chart_mean_value[i] = x[i];
        }
    }
    return true;
}
void polycube_flattening_interface::update_updownchart(TetStructure<double> *tet_mesh_)
{
    if (update_updownchart_flag)
        return;
    update_updownchart_flag = true;
    std::cout << "--------------------Adding Constraints for Flattening------------------------" << std::endl;
    if (bd_pts.empty()) get_bd_face_chart_label(tet_mesh_);
    std::set<std::pair<int, int>> updown_array, updown_array_ori;
    for (size_t i = 0; i < up_chart_id.size(); i++)
    {
        updown_array.insert(std::pair<int, int>(up_chart_id[i], down_chart_id[i]));
        updown_array_ori.insert(std::pair<int, int>(up_chart_id[i], down_chart_id[i]));
    }
    CVec<double, 3> axis_dir[3] = { CVec<double, 3>(1.0, 0.0, 0.0), CVec<double, 3>(0.0, 1.0, 0.0), CVec<double, 3>(0.0, 0.0, 1.0) };
    CVec<double, 3> leftdownpt = bd_pts[0];
    for (size_t i = 1; i < bd_pts.size(); i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            if (leftdownpt[j] > bd_pts[i][j])
                leftdownpt[j] = bd_pts[i][j];
        }
    }
    std::vector<std::vector<int>> axis2pqedgeid(3), axis2chart(3);
    std::vector<int> pqedge2axis(polycube_edges.size());
    std::vector<std::vector<int>> chart2bd_faceid(polycube_chart_label.size());
    assert(bd_chart.size() == bd_faces.size() / 3);
    for (size_t i = 0; i < bd_faces.size() / 3; i++)
    {
        chart2bd_faceid[bd_chart[i]].push_back(i);
    }
    for (size_t i = 0; i < polycube_edges.size(); i++)
    {
        int chart0 = polycube_edges[i][2];
        int chart1 = polycube_edges[i][3];
        int axis0 = polycube_chart_label[chart0] / 2;
        int axis1 = polycube_chart_label[chart1] / 2;
        int select_axis = 3 - axis0 - axis1;
        axis2pqedgeid[select_axis].push_back(i);
        pqedge2axis[i] = select_axis;
    }
    for (size_t i = 0; i < polycube_chart_label.size(); i++)
    {
        int axis = polycube_chart_label[i] / 2;
        axis2chart[axis].push_back(i);
    }
    for (int i = 0; i < 3; i++)
    {
        int curaxis = i;
        int otheraxis[2];
        otheraxis[0] = (i + 1) % 3;
        otheraxis[1] = (i + 2) % 3;
        std::vector<unsigned> cur_faces;
        std::vector<int> face_c2o;
        std::vector<double> coord[2];
        for (size_t j = 0; j < axis2chart[i].size(); j++)
        {
            int chartid = axis2chart[i][j];
            for (size_t k = 0; k < chart2bd_faceid[chartid].size(); k++)
            {
                int fid = chart2bd_faceid[chartid][k];
                cur_faces.push_back(bd_faces[3 * fid]);
                cur_faces.push_back(bd_faces[3 * fid + 1]);
                cur_faces.push_back(bd_faces[3 * fid + 2]);
                face_c2o.push_back(fid);
            }
        }
        for (size_t j = 0; j < 2; j++)
        {
            int axis = otheraxis[j];
            for (size_t k = 0; k < axis2pqedgeid[axis].size(); k++)
            {
                int pqid = axis2pqedgeid[axis][k];
                int vid[2];
                vid[0] = polycube_edges[pqid][0];
                vid[1] = polycube_edges[pqid][1];
                for (size_t iter = 0; iter < 2; iter++)
                {
                    switch (axis)
                    {
                    case 0:
                        coord[j].push_back(dpx[vid[iter]]);
                        break;
                    case 1:
                        coord[j].push_back(dpy[vid[iter]]);
                        break;
                    case 2:
                        coord[j].push_back(dpz[vid[iter]]);
                        break;
                    }
                }
            }
        }
        std::sort(coord[0].begin(), coord[0].end());
        std::sort(coord[1].begin(), coord[1].end());
        std::vector<double> coord_final[2];
        for (size_t j = 0; j < 2; j++)
        {
            for (size_t k = 0; k < coord[j].size(); k++)
            {
                if (k == 0 || std::abs(coord[j][k] - coord[j][k - 1]) > MIN_EDGE_DIST)
                {
                    coord_final[j].push_back(coord[j][k]);
                }
            }
        }
        std::vector<std::pair<CVec<double, 3>, CVec<double, 3>>> rays;
        std::vector<std::set<unsigned>> ignored_faces;
        std::vector<CVec<double, 3>> first_intersection_pts;
        std::vector<std::vector<int>> all_intersection_face_id_cur, all_intersection_face_id_ori;
        for (size_t j = 0; j < coord_final[0].size() - 1; j++)
        {
            for (size_t k = 0; k < coord_final[1].size() - 1; k++)
            {
                CVec<double, 3> begin_pt;
                begin_pt[curaxis] = leftdownpt[curaxis] - 1.0;
                begin_pt[otheraxis[0]] = (coord_final[0][j] + coord_final[0][j + 1]) / 2.0;
                begin_pt[otheraxis[1]] = (coord_final[1][k] + coord_final[1][k + 1]) / 2.0;
                rays.push_back(std::pair<CVec<double, 3>, CVec<double, 3>>(begin_pt, axis_dir[curaxis]));
            }
        }
        ig::CGALHelper::rays_triangles_intersection(rays, bd_pts, cur_faces, ignored_faces, first_intersection_pts, &all_intersection_face_id_cur);
        all_intersection_face_id_ori.resize(all_intersection_face_id_cur.size());
        for (size_t j = 0; j < all_intersection_face_id_cur.size(); j++)
        {
            for (size_t k = 0; k < all_intersection_face_id_cur[j].size(); k++)
            {
                all_intersection_face_id_ori[j].push_back(face_c2o[all_intersection_face_id_cur[j][k]]);
            }
            std::vector<std::pair<int, double>> intersection_chart_value_sorted_all, intersection_chart_value_sorted;
            std::vector<int> intersection_label_sorted;
            std::vector<bool> intersection_dir_sorted;
            for (size_t k = 0; k < all_intersection_face_id_ori[j].size(); k++)
            {
                if (k == 0 || bd_chart[all_intersection_face_id_ori[j][k]] != bd_chart[all_intersection_face_id_ori[j][k - 1]])
                {
                    int axis = bd_label[all_intersection_face_id_ori[j][k]] / 2;
                    if (axis == curaxis)
                    {
                        int chart = bd_chart[all_intersection_face_id_ori[j][k]];
                        intersection_chart_value_sorted_all.push_back(std::pair<int, double>(chart, chart_mean_value[chart]));
                    }
                }
            }
            std::sort(intersection_chart_value_sorted_all.begin(), intersection_chart_value_sorted_all.end(), [](std::pair<int, double>p1, std::pair<int, double>p2) {
                return p1.second < p2.second;
                });
            for (size_t k = 0; k < intersection_chart_value_sorted_all.size(); k++)
            {
                int chart = intersection_chart_value_sorted_all[k].first;
                int label = polycube_chart_label[chart];
                intersection_label_sorted.push_back(label);
                intersection_dir_sorted.push_back((bool)(1 - label % 2));
            }
            if (!intersection_dir_sorted.empty())
            {
                for (size_t k = 0; k < intersection_dir_sorted.size() - 1; k++)
                {
                    bool cur_dir = intersection_dir_sorted[k];
                    int next_id = -1;
                    for (int iter = k + 1; iter < intersection_dir_sorted.size(); iter++)
                    {
                        if (intersection_dir_sorted[iter] != cur_dir)
                        {
                            next_id = iter;
                            break;
                        }
                    }
                    if (next_id != -1)
                    {
                        updown_array.insert(std::pair<int, int>(intersection_chart_value_sorted_all[next_id].first, intersection_chart_value_sorted_all[k].first));
                    }
                }
            }
        }
    }
    std::vector<std::vector<int>> edges_with_same_chart;
    edges_with_same_chart.resize(polycube_chart_label.size());
    for (size_t i = 0; i < polycube_edges.size(); i++)
    {
        int chart_idx1 = polycube_edges[i][2];
        int chart_idx2 = polycube_edges[i][3];
        edges_with_same_chart[chart_idx1].push_back(i);
        edges_with_same_chart[chart_idx2].push_back(i);
    }
    for (size_t i = 0; i < polycube_chart_label.size(); i++)
    {
        if (i == 35)
        {
            int a = 1;
        }
        int curchart = i;
        int chartaxis = polycube_chart_label[curchart] / 2;
        std::map<int, std::vector<int>> local_a2p;
        for (size_t j = 0; j < edges_with_same_chart[curchart].size(); j++)
        {
            int eid = edges_with_same_chart[curchart][j];
            local_a2p[pqedge2axis[eid]].push_back(eid);
        }
        CVec<double, 3> center(0.0, 0.0, 0.0);
        center[chartaxis] = chart_mean_value[curchart];
        assert(local_a2p.size() == 2);
        for (auto a2p : local_a2p)
        {
            int axis = a2p.first;
            std::vector<std::pair<int, std::pair<double, double>>> pqidsegm;
            for (size_t j = 0; j < a2p.second.size(); j++)
            {
                int eid = a2p.second[j];
                int v0 = polycube_edges[eid][0];
                int v1 = polycube_edges[eid][1];
                CVec<double, 3> p0(dpx[v0], dpy[v0], dpz[v0]), p1(dpx[v1], dpy[v1], dpz[v1]);
                double s0 = (p0 - center).Dot(axis_dir[axis]);
                double s1 = (p1 - center).Dot(axis_dir[axis]);
                if (s0 > s1) std::swap(s0, s1);
                pqidsegm.push_back(std::pair<int, std::pair<double, double>>(eid, std::pair<double, double>(s0, s1)));
            }
            for (size_t i1 = 0; i1 < pqidsegm.size() - 1; i1++)
            {
                for (size_t i2 = i1 + 1; i2 < pqidsegm.size(); i2++)
                {
                    bool overlap = ig::CGALHelper::line_segm_intersection<double>(pqidsegm[i1].second, pqidsegm[i2].second);
                    if (overlap)
                    {
                        int eid1 = pqidsegm[i1].first, eid2 = pqidsegm[i2].first;
                        int otherchart1 = polycube_edges[eid1][2], otherchart2 = polycube_edges[eid2][2];
                        if (otherchart1 == curchart) otherchart1 = polycube_edges[eid1][3];
                        if (otherchart2 == curchart) otherchart2 = polycube_edges[eid2][3];
                        if (chart_mean_value[otherchart1] < chart_mean_value[otherchart2]) std::swap(otherchart1, otherchart2);
                        updown_array.insert(std::pair<int, int>(otherchart1, otherchart2));
                    }
                }
            }
        }
    }
    for (auto it : updown_array)
    {
        auto it_find = updown_array_ori.find(it);
        if (it_find == updown_array_ori.end())
        {
            up_chart_id.push_back(it.first);
            down_chart_id.push_back(it.second);
            diff_up_down.push_back(cube_len);
        }
    }
}
bool polycube_flattening_interface::find_all_chart_value_equal_face_rounding(TetStructure<double> *tet_mesh_, double offset)
{
    find_all_corner(tet_mesh_);
    std::vector<MSKboundkeye> bkc; std::vector<double> blc; std::vector<double> buc;
    std::vector<MSKboundkeye> bkx; std::vector<double> blx; std::vector<double> bux;
    std::vector<MSKlidxt> aptrb; std::vector<MSKidxt> asub; std::vector<double> aval;
    std::vector<MSKidxt> qsubi; std::vector<MSKidxt> qsubj; std::vector<double> qval;
    double ratio_diff = 0.8;
    std::vector< Eigen::Triplet<double> > coef;
    int up_down_chart_count = up_chart_id.size();
    int chart_number = chart_mean_value.size();
    double diff_u_d; int up_chart, down_chart;
    for (int i = 0; i < up_down_chart_count; ++i)
    {
        up_chart = up_chart_id[i];
        down_chart = down_chart_id[i];
        diff_u_d = diff_up_down[i] * ratio_diff;
        if (hex_meshing_flag)
        {
            if (diff_u_d < cube_len)
            {
                diff_u_d = cube_len;
            }
        }
        coef.push_back(Eigen::Triplet<double>(i, up_chart, +1.0));
        coef.push_back(Eigen::Triplet<double>(i, down_chart, -1.0));
        bkc.push_back(MSK_BK_LO); blc.push_back(diff_u_d); buc.push_back(+MSK_INFINITY);
    }
    int cut_num = 0;
    cut_num = cut_to_chart_pair.size();
    for (size_t i = 0; i < cut_to_chart_pair.size(); i++)
    {
        int chart1 = cut_to_chart_pair[i].first;
        int chart2 = cut_to_chart_pair[i].second;
        std::vector<int> chart1_idx, chart2_idx;
        for (size_t j = 0; j < polycube_edges.size(); j++)
        {
            if ((polycube_edges[j][2] == chart1 && polycube_edges[j][3] != chart2)
                || (polycube_edges[j][3] == chart1 && polycube_edges[j][2] != chart2))
                chart1_idx.push_back(j);
            if ((polycube_edges[j][2] == chart2 && polycube_edges[j][3] != chart1)
                || (polycube_edges[j][3] == chart2 && polycube_edges[j][2] != chart1))
                chart2_idx.push_back(j);
        }
        int final_idx1(-1), final_idx2(-1);
        for (size_t j = 0; j < chart1_idx.size(); j++)
        {
            if (up_chart_id[chart1_idx[j]] == chart2 || down_chart_id[chart1_idx[j]] == chart2)
            {
                final_idx1 = chart1_idx[j];
                break;
            }
        }
        for (size_t j = 0; j < chart2_idx.size(); j++)
        {
            if (up_chart_id[chart2_idx[j]] == chart1 || down_chart_id[chart2_idx[j]] == chart1)
            {
                final_idx2 = chart2_idx[j];
                break;
            }
        }
        coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, up_chart_id[final_idx1], +1.0));
        coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, down_chart_id[final_idx1], -1.0));
        coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, up_chart_id[final_idx2], -1.0));
        coef.push_back(Eigen::Triplet<double>(i + up_down_chart_count, down_chart_id[final_idx2], +1.0));
        bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
    }
    int n_valence_six_constraint = 0;
    std::map<int, std::set<int>> vert2chart;
    for (size_t i = 0; i < polycube_edges.size(); i++)
    {
        int id0 = polycube_edges[i][0];
        int id1 = polycube_edges[i][1];
        int chart0 = polycube_edges[i][2];
        int chart1 = polycube_edges[i][3];
        vert2chart[id0].insert(chart0);
        vert2chart[id0].insert(chart1);
        vert2chart[id1].insert(chart0);
        vert2chart[id1].insert(chart1);
    }
    for (auto vc : vert2chart)
    {
        if (vc.second.size() == 6)
        {
            std::vector<int> xchart, ychart, zchart;
            for (auto chart : vc.second)
            {
                int axis = polycube_chart_label[chart] / 2;
                switch (axis)
                {
                case 0:
                    xchart.push_back(chart);
                    break;
                case 1:
                    ychart.push_back(chart);
                    break;
                case 2:
                    zchart.push_back(chart);
                    break;
                default:
                    break;
                }
            }
            assert(xchart.size() == 2 && ychart.size() == 2 && zchart.size() == 2);
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint, xchart[0], +1.0));
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint, xchart[1], -1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 1, ychart[0], +1.0));
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 1, ychart[1], -1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 2, zchart[0], +1.0));
            coef.push_back(Eigen::Triplet<double>(up_down_chart_count + n_valence_six_constraint + 2, zchart[1], -1.0));
            bkc.push_back(MSK_BK_FX); blc.push_back(0.0); buc.push_back(0.0);
            n_valence_six_constraint += 3;
        }
    }
    Eigen::SparseMatrix<double> A(up_down_chart_count + cut_num + n_valence_six_constraint, chart_number);
    A.setFromTriplets(coef.begin(), coef.end());
    int size = A.nonZeros(); int cols = A.cols();
    aptrb.resize(cols + 1); asub.resize(size); aval.resize(size);
    aptrb[cols] = size;
    int ind = 0;
    for (int k = 0; k < cols; ++k)
    {
        aptrb[k] = ind;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
        {
            aval[ind] = it.value();
            asub[ind] = it.row();
            ++ind;
        }
    }
    std::vector<double> c(chart_number, 0.0);
    for (int i = 0; i < chart_number; ++i)
    {
        qsubi.push_back(i); qsubj.push_back(i); qval.push_back(2.0);
        c[i] = -2.0 * chart_mean_value[i];
    }
    bkx.resize(chart_number, MSK_BK_FR); blx.resize(chart_number, -MSK_INFINITY); bux.resize(chart_number, +MSK_INFINITY);
    std::vector<double> x(chart_number, 0.0); int mosek_count = 0;
    bool solve_success = solveConvexQuadPorgramming_mosek(bkc, blc, buc, bkx, blx, bux, aptrb, asub, aval, qsubi, qsubj, qval, c, x);
    while (!solve_success)
    {
        ratio_diff *= 0.8;
        for (int i = 0; i < up_down_chart_count; ++i)
        {
            diff_u_d = diff_up_down[i] * ratio_diff;
            if (hex_meshing_flag)
                if (diff_u_d < cube_len)
                {
                    diff_u_d = cube_len;
                }
            blc[i] = diff_u_d;
        }
        solve_success = solveConvexQuadPorgramming_mosek(bkc, blc, buc, bkx, blx, bux, aptrb, asub, aval, qsubi, qsubj, qval, c, x);
        ++mosek_count;
        if (mosek_count > 20) break;
    }
    if (mosek_count > 20)
    {
        printf("-------------------------------------------\nNo Solution!!!!!\n");
        return false;
    }
    else
    {
        for (int i = 0; i < chart_number; ++i)
        {
            if (hex_meshing_flag)
            {
                chart_mean_value[i] = cube_len * std::floor(x[i] / cube_len + 0.5);
                std::cout << "solution: " << chart_mean_value[i] / cube_len << std::endl;
                double offset_ratio = offset;
                if (polycube_chart_label[i] % 2 == 0)
                {
                    chart_mean_value[i] += cube_len * offset_ratio;
                }
                else
                {
                    chart_mean_value[i] -= cube_len * offset_ratio;
                }
            }
            else
            {
                chart_mean_value[i] = x[i];
                std::cout << "solution: " << chart_mean_value[i] << std::endl;
            }
        }
    }
    return true;
}
void polycube_flattening_interface::boundary_mapping_polycube_equal_face(TetStructure<double>* tet_mesh_)
{
    const std::vector<Tetrahedron<double>*> &tetras = tet_mesh_->tetras;
    std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
    int nv = tetra_vertices.size();
    int nc = tetras.size();
    std::vector<int> map_v(nv, -1); std::vector<int> map_c(nc, -1);
    int chart_size = polycube_chart_label.size();
    int chart_pair_size = cut_to_chart_pair.size();
    for (size_t i = 0; i < chart_pair_size; i++)
    {
        int matrix_type = cut_types[i];
        int chart1 = cut_to_chart_pair[i].first;
        int chart2 = cut_to_chart_pair[i].second;
        int var_v_count = 0; int var_c_count = 0;
        std::vector<int>& one_chart = polycube_chart[chart1];
        std::vector<int> chart_boundary_flag;
        chart_boundary_flag.resize(nv, -1);
        for (int j = 0; j < one_chart.size(); ++j)
        {
            int bf_id = one_chart[j];
            for (int k = 0; k < 3; ++k)
            {
                int v_id = bfv_id[bf_id][k];
                bool fixed = ((is_polycube_handles[3 * v_id + 0] == 1) && (is_polycube_handles[3 * v_id + 1] == 1) && (is_polycube_handles[3 * v_id + 2] == 1));
                if (!fixed)
                {
                    if (map_v[v_id] < 0)
                    {
                        map_v[v_id] = var_v_count; ++var_v_count;
                    }
                    if (map_c[bf_id] < 0)
                    {
                        map_c[bf_id] = var_c_count; ++var_c_count;
                    }
                }
                else
                {
                    chart_boundary_flag[v_id] = 1;
                }
            }
        }
        int chart_label = polycube_chart_label[chart1];
        int chart_label_pair = polycube_chart_label[chart2];
        polycube_boundary_mapping_IVF pbm;
        pbm.optimize_IVF(tet_mesh_, chart_label, src_pos, bfv_id, bvf_id, map_v, var_v_count, map_c, var_c_count);
        int common_vert_idx = cut_common_verts_idx[i];
        CVec<double, 3> np_common = tetra_vertices[common_vert_idx]->pos;
        for (size_t j = 0; j < nv; j++)
        {
            if (map_v[j] >= 0 || chart_boundary_flag[j] > 0)
            {
                if (three_cut_vert_flag.size() > 0)
                {
                    if (three_cut_vert_flag[j] > 0)
                        continue;
                }
                if (j == common_vert_idx)
                    continue;
                int pair_idx = vert_pairs_map[j];
                if (pair_idx < 0)
                {
                }
                CVec<double, 3> np = tetra_vertices[j]->pos;
                CVec<double, 3> np_pair;
                double rot_vec[3];
                double ori_vec[3] = { np[0] - np_common[0], np[1] - np_common[1], np[2] - np_common[2] };
                for (size_t b = 0; b < 3; b++)
                {
                    rot_vec[b] = 0.0;
                    for (size_t c = 0; c < 3; c++)
                    {
                        rot_vec[b] += type_matrix[matrix_type][b][c] * ori_vec[c];
                    }
                }
                np_pair[0] = np_common[0] + rot_vec[0];
                np_pair[1] = np_common[1] + rot_vec[1];
                np_pair[2] = np_common[2] + rot_vec[2];
                tetra_vertices[pair_idx]->pos = np_pair;
            }
        }
        for (int j = 0; j < nv; ++j)
        {
            if (map_v[j] >= 0 || chart_boundary_flag[j] > 0)
            {
                is_polycube_handles[3 * j + 0] = 1; is_polycube_handles[3 * j + 1] = 1; is_polycube_handles[3 * j + 2] = 1;
                int pair_idx = vert_pairs_map[j];
                is_polycube_handles[3 * pair_idx + 0] = 1; is_polycube_handles[3 * pair_idx + 1] = 1; is_polycube_handles[3 * pair_idx + 2] = 1;
            }
            map_v[j] = -1;
        }
        for (int j = 0; j < nc; ++j)
        {
            map_c[j] = -1;
        }
    }
    for (int i = 0; i < chart_size; ++i)
    {
        bool cut_flag = false;
        for (size_t j = 0; j < cut_to_chart_pair.size(); j++)
        {
            if (i == cut_to_chart_pair[j].first)
            {
                cut_flag = true;
                break;
            }
            if (i == cut_to_chart_pair[j].second)
            {
                cut_flag = true;
                break;
            }
        }
        if (cut_flag)
        {
            continue;
        }
        int var_v_count = 0; int var_c_count = 0;
        std::vector<int>& one_chart = polycube_chart[i];
        for (int j = 0; j < one_chart.size(); ++j)
        {
            int bf_id = one_chart[j];
            for (int k = 0; k < 3; ++k)
            {
                int v_id = bfv_id[bf_id][k];
                bool fixed = ((is_polycube_handles[3 * v_id + 0] == 1) && (is_polycube_handles[3 * v_id + 1] == 1) && (is_polycube_handles[3 * v_id + 2] == 1));
                if (!fixed)
                {
                    if (map_v[v_id] < 0)
                    {
                        map_v[v_id] = var_v_count; ++var_v_count;
                    }
                    if (map_c[bf_id] < 0)
                    {
                        map_c[bf_id] = var_c_count; ++var_c_count;
                    }
                }
            }
        }
        int chart_label = polycube_chart_label[i];
        polycube_boundary_mapping_IVF pbm;
        pbm.optimize_IVF(tet_mesh_, chart_label, src_pos, bfv_id, bvf_id, map_v, var_v_count, map_c, var_c_count);
        for (int j = 0; j < nv; ++j)
        {
            if (map_v[j] >= 0)
            {
                is_polycube_handles[3 * j + 0] = 1; is_polycube_handles[3 * j + 1] = 1; is_polycube_handles[3 * j + 2] = 1;
            }
            map_v[j] = -1;
        }
        for (int j = 0; j < nc; ++j)
        {
            map_c[j] = -1;
        }
    }
}
bool polycube_flattening_interface::deform_ARAP_polycube_equal_face(TetStructure<double>* tet_mesh_, bool int_solver_flag, bool add_constraint, bool set_min_diff, double offset, bool double_cube_length_flag)
{
    find_all_corner(tet_mesh_);
    if (add_constraint)
    {
        update_updownchart(tet_mesh_);
    }
    bool find_chart_flag;
    if (double_cube_length_flag)
        cube_len = cube_len * 2;
    if (!hex_meshing_flag)
    {
        if (set_min_diff)
        {
            find_chart_flag = find_all_chart_value_preprocessing(tet_mesh_);
        }
        else
        {
            find_chart_flag = find_all_chart_value_equal_face_mine(tet_mesh_, offset);
        }
    }
    else
    {
        if (int_solver_flag)
        {
            find_chart_flag = find_all_chart_value_equal_face_mine(tet_mesh_, offset);
        }
        else
        {
            find_chart_flag = find_all_chart_value_equal_face_rounding(tet_mesh_, offset);
        }
    }
    if (double_cube_length_flag)
        cube_len = cube_len / 2;
    if (!find_chart_flag)
        return false;
    int nv = vertex_type.size();
    is_polycube_handles.clear(); is_polycube_handles.resize(3 * nv, -1);
    std::vector<int> polycube_edge_v; std::vector<int> polycube_edge_v_neighbor;
    std::vector<int> polycube_edge_v_type;
    const std::vector<Tetrahedron<double>*> &tetras = tet_mesh_->tetras;
    std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
    std::vector<int> corner_pts_idx;
    std::vector<std::vector<int>> corner_faces;
    std::map<int, ig::CVec<double, 3>> idx2coord;
    corner_faces.resize(chart_mean_value.size());
    for (int i = 0; i < nv; i++)
    {
        OpenVolumeMesh::Geometry::Vec3i& v_type = vertex_type[i];
        if (v_type[0] >= 0 && v_type[1] >= 0 && v_type[2] >= 0)
        {
            corner_pts_idx.push_back(i);
            CVec<double, 3> np = tetra_vertices[i]->pos;
            idx2coord.insert(std::pair<int, ig::CVec<double, 3>>(i, np));
            corner_faces[v_type[0]].push_back(i);
            corner_faces[v_type[1]].push_back(i);
            corner_faces[v_type[2]].push_back(i);
        }
    }
    for (int i = 0; i < nv; ++i)
    {
        OpenVolumeMesh::Geometry::Vec3i& v_type = vertex_type[i];
        CVec<double, 3> np = tetra_vertices[i]->pos;
        if (v_type[0] >= 0)
        {
            np[0] = chart_mean_value[v_type[0]];
            dpx[i] = np[0];
            is_polycube_handles[3 * i + 0] = 1;
        }
        if (v_type[1] >= 0)
        {
            np[1] = chart_mean_value[v_type[1]];
            dpy[i] = np[1];
            is_polycube_handles[3 * i + 1] = 1;
        }
        if (v_type[2] >= 0)
        {
            np[2] = chart_mean_value[v_type[2]];
            dpz[i] = np[2];
            is_polycube_handles[3 * i + 2] = 1;
        }
        tetra_vertices[i]->pos = np;
        if (v_type[0] >= 0 && v_type[1] >= 0 && v_type[2] < 0)
        {
            polycube_edge_v.push_back(i);
            polycube_edge_v_type.push_back(2);
            is_polycube_handles[3 * i + 2] = 1;
        }
        else if (v_type[0] < 0 && v_type[1] >= 0 && v_type[2] >= 0)
        {
            polycube_edge_v.push_back(i);
            polycube_edge_v_type.push_back(0);
            is_polycube_handles[3 * i + 0] = 1;
        }
        else if (v_type[0] >= 0 && v_type[1] < 0 && v_type[2] >= 0)
        {
            polycube_edge_v.push_back(i);
            polycube_edge_v_type.push_back(1);
            is_polycube_handles[3 * i + 1] = 1;
        }
    }
    bool flag_translation = true;
    if (flag_translation)
    {
        std::vector<ig::CVec<double, 3>> chart_translation;
        for (size_t i = 0; i < chart_mean_value.size(); i++)
        {
            int face_size = corner_faces[i].size();
            ig::CVec<double, 3> tmp_translation(0.0, 0.0, 0.0);
            for (size_t j = 0; j < face_size; j++)
            {
                int tmp_idx = corner_faces[i][j];
                tmp_translation = tmp_translation + tetra_vertices[tmp_idx]->pos - idx2coord[tmp_idx];
            }
            tmp_translation = tmp_translation / (double)face_size;
            int ignore_label = polycube_chart_label[i];
            tmp_translation[ignore_label / 2] = 0.0;
            chart_translation.push_back(tmp_translation);
        }
        for (size_t i = 0; i < nv; i++)
        {
            std::vector<int> tmp_vert_belong_chart;
            OpenVolumeMesh::Geometry::Vec3i& v_type = vertex_type[i];
            if (v_type[0] >= 0)
                tmp_vert_belong_chart.push_back(v_type[0]);
            if (v_type[1] >= 0)
                tmp_vert_belong_chart.push_back(v_type[1]);
            if (v_type[2] >= 0)
                tmp_vert_belong_chart.push_back(v_type[2]);
            if (tmp_vert_belong_chart.size() == 1)
            {
                int tmp_chart = tmp_vert_belong_chart[0];
                CVec<double, 3> np = tetra_vertices[i]->pos;
                np = np + chart_translation[tmp_chart];
                tetra_vertices[i]->pos = np;
            }
        }
    }
    for (int i = 0; i < polycube_edge_v.size(); ++i)
    {
        int v_id = polycube_edge_v[i];
        CVec<double, 3> np = tetra_vertices[v_id]->pos;
        std::vector<int>& one_vv = bvv_id[v_id];
        int vv_size = one_vv.size(); int add_nv = 0;
        for (int j = 0; j < vv_size; ++j)
        {
            int vj_id = one_vv[j];
            OpenVolumeMesh::Geometry::Vec3i& vj_type = vertex_type[vj_id]; int chart_count = 0; bool can_add = false;
            if (vj_type[0] >= 0) ++chart_count;
            if (vj_type[1] >= 0) ++chart_count;
            if (vj_type[2] >= 0) ++chart_count;
            if (chart_count == 3)
            {
                can_add = true;
            }
            else if (chart_count == 2)
            {
                if (vj_type[polycube_edge_v_type[i]] < 0)
                {
                    can_add = true;
                }
            }
            if (can_add)
            {
                CVec<double, 3> npj = tetra_vertices[vj_id]->pos;
                bool nv_ok = true;
                for (int k = 0; k < 3; ++k)
                {
                    if (k != polycube_edge_v_type[i])
                    {
                        if (std::abs(npj[k] - np[k]) > 1e-10) nv_ok = false;
                    }
                }
                if (nv_ok)
                {
                    polycube_edge_v_neighbor.push_back(vj_id);
                    ++add_nv;
                }
            }
        }
        if (add_nv > 0 && add_nv != 2)
        {
            if (add_nv > 2)
            {
                int pevn_size = polycube_edge_v_neighbor.size();
                std::vector<int> delete_v;
                std::vector<int> visited_f(tetras.size(), 0);
                std::vector<int>& one_vf = bvf_id[v_id];
                for (int j = 0; j < one_vf.size(); ++j)
                {
                    visited_f[one_vf[j]] += 1;
                }
                for (int j = pevn_size - add_nv; j < pevn_size; ++j)
                {
                    int vj_id = polycube_edge_v_neighbor[j];
                    std::vector<int>& one_vj_f = bvf_id[vj_id];
                    int two_f[2]; int two_count = 0;
                    for (int k = 0; k < one_vj_f.size(); ++k)
                    {
                        if (visited_f[one_vj_f[k]] == 1)
                        {
                            two_f[two_count] = one_vj_f[k]; ++two_count;
                        }
                    }
                    if (bf_chart[two_f[0]] == bf_chart[two_f[1]])
                    {
                        delete_v.push_back(j);
                    }
                }
                for (int j = 0; j < delete_v.size(); ++j)
                {
                    polycube_edge_v_neighbor.erase(polycube_edge_v_neighbor.begin() + delete_v[delete_v.size() - 1 - j]);
                }
            }
        }
    }
    if (true)
    {
        if (polycube_edge_v_neighbor.size() == 2 * polycube_edge_v.size())
        {
            for (int i = 0; i < 20; ++i)
            {
                for (int j = 0; j < polycube_edge_v.size(); ++j)
                {
                    int v_id = polycube_edge_v[j];
                    CVec<double, 3> np = tetra_vertices[v_id]->pos;
                    CVec<double, 3> p0 = tetra_vertices[polycube_edge_v_neighbor[2 * j + 0]]->pos;
                    CVec<double, 3> p1 = tetra_vertices[polycube_edge_v_neighbor[2 * j + 1]]->pos;
                    int edge_type = polycube_edge_v_type[j];
                    np[edge_type] = 0.5*(p0[edge_type] + p1[edge_type]);
                    tetra_vertices[v_id]->pos = np;
                    dpx[v_id] = np[0]; dpy[v_id] = np[1]; dpz[v_id] = np[2];
                }
            }
        }
    }
    if (true)
    {
        assert(three_cut_common_vert.size() == three_cut_vert.size());
        for (size_t i = 0; i < three_cut_common_vert.size(); i++)
        {
            int common_vert = three_cut_common_vert[i];
            int n_segment = three_cut_vert[i].size();
            int end_idx[3];
            end_idx[0] = three_cut_vert[i][0][0];
            end_idx[1] = three_cut_vert[i][0][1];
            end_idx[2] = three_cut_vert[i][0][2];
            CVec<double, 3> common_pos = tetra_vertices[common_vert]->pos;
            CVec<double, 3> end_pos[3];
            end_pos[0] = tetra_vertices[end_idx[0]]->pos;
            end_pos[1] = tetra_vertices[end_idx[1]]->pos;
            end_pos[2] = tetra_vertices[end_idx[2]]->pos;
            CVec<double, 3> seg[3];
            seg[0] = (common_pos - end_pos[0]) / n_segment;
            seg[1] = (common_pos - end_pos[1]) / n_segment;
            seg[2] = (common_pos - end_pos[2]) / n_segment;
            for (int j = 1; j < three_cut_vert[i].size(); j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    int v_id = three_cut_vert[i][j][k];
                    CVec<double, 3> np = tetra_vertices[v_id]->pos;
                    np = end_pos[k] + j * seg[k];
                    tetra_vertices[v_id]->pos = np;
                    dpx[v_id] = np[0]; dpy[v_id] = np[1]; dpz[v_id] = np[2];
                }
            }
        }
    }
    boundary_mapping_polycube_equal_face(tet_mesh_);
    double min_singular_value = 0.001;
    int nc = tetras.size();
    Eigen::Matrix3d VS, Q, U, V, R;
    std::vector < std::vector < double >> Cell_R(nc);
    for (int c_id = 0; c_id < nc; ++c_id)
    {
        Cell_R[c_id].resize(9);
    }
    for (int iter = 0; iter < 3; ++iter)
    {
        for (int c_id = 0; c_id < nc; ++c_id)
        {
            std::vector<int>& cv_id = cell_vertex[c_id];
            double cpx = dpx[cv_id[0]]; double cpy = dpy[cv_id[0]]; double cpz = dpz[cv_id[0]];
            VS(0, 0) = dpx[cv_id[1]] - cpx; VS(1, 0) = dpy[cv_id[1]] - cpy; VS(2, 0) = dpz[cv_id[1]] - cpz;
            VS(0, 1) = dpx[cv_id[2]] - cpx; VS(1, 1) = dpy[cv_id[2]] - cpy; VS(2, 1) = dpz[cv_id[2]] - cpz;
            VS(0, 2) = dpx[cv_id[3]] - cpx; VS(1, 2) = dpy[cv_id[3]] - cpy; VS(2, 2) = dpz[cv_id[3]] - cpz;
            Q = VS * cell_S[c_id];
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(Q, Eigen::ComputeFullU | Eigen::ComputeFullV);
            U = svd.matrixU(); V = svd.matrixV();
            double sig0 = svd.singularValues()[0];
            double sig1 = svd.singularValues()[1];
            double sig2 = svd.singularValues()[2];
            if (Q.determinant() > 0)
            {
                R = V * U.transpose();
            }
            else
            {
                U(0, 2) = -U(0, 2); U(1, 2) = -U(1, 2); U(2, 2) = -U(2, 2);
                R = V * U.transpose();
            }
            Cell_R[c_id][0] = R(0, 0); Cell_R[c_id][1] = R(0, 1); Cell_R[c_id][2] = R(0, 2);
            Cell_R[c_id][3] = R(1, 0); Cell_R[c_id][4] = R(1, 1); Cell_R[c_id][5] = R(1, 2);
            Cell_R[c_id][6] = R(2, 0); Cell_R[c_id][7] = R(2, 1); Cell_R[c_id][8] = R(2, 2);
        }
        for (int i = 0; i < 3; ++i)
        {
            int var_count = 0; std::vector<int> map_v(nv, -1);
            for (int j = 0; j < nv; ++j)
            {
                if (is_polycube_handles[3 * j + i] != 1)
                {
                    map_v[j] = var_count;
                    ++var_count;
                }
            }
            Sparse_Matrix A = new Sparse_Matrix(nc * 3, var_count, NOSYM, CCS, 1);
            std::vector<double> cs(4);
            for (int c_id = 0; c_id < nc; ++c_id)
            {
                double cv = std::sqrt(std::abs(cell_volume[c_id])) * 100;
                std::vector<int>& cv_id = cell_vertex[c_id];
                Eigen::Matrix3d& CS = cell_S[c_id];
                std::vector<double>& CR = Cell_R[c_id];
                for (int j = 0; j < 3; ++j)
                {
                    cs[1] = CS(0, j)*cv; cs[2] = CS(1, j)*cv; cs[3] = CS(2, j)*cv;
                    cs[0] = -(cs[1] + cs[2] + cs[3]);
                    for (int k = 0; k < 4; ++k)
                    {
                        int var_id = map_v[cv_id[k]];
                        if (var_id < 0)
                        {
                            A.fill_rhs_entry(3 * c_id + j, -cs[k] * tetra_vertices[cv_id[k]]->pos[i]);
                        }
                        else
                        {
                            A.fill_entry(3 * c_id + j, var_id, cs[k]);
                        }
                    }
                    A.fill_rhs_entry(3 * c_id + j, cv*CR[i + 3 * j]);
                }
            }
            Sparse_Matrix* D = TransposeTimesSelf(&A, CCS, SYM_LOWER, true);
            solve_by_CHOLMOD(D);
            const std::vector<double>& xyz = D->get_solution();
            for (int j = 0; j < nv; ++j)
            {
                if (is_polycube_handles[3 * j + i] != 1)
                {
                    int var_id = map_v[j];
                    CVec<double, 3> np = tetra_vertices[j]->pos;
                    np[i] = xyz[var_id];
                    tetra_vertices[j]->pos = np;
                }
            }
            delete D;
        }
        for (int j = 0; j < nv; ++j)
        {
            CVec<double, 3> np = tetra_vertices[j]->pos;
            dpx[j] = np[0]; dpy[j] = np[1]; dpz[j] = np[2];
        }
    }
    for (int j = 0; j < nv; ++j)
    {
        CVec<double, 3> np = tetra_vertices[j]->pos;
        dpx[j] = np[0]; dpy[j] = np[1]; dpz[j] = np[2];
    }
    return true;
}
bool polycube_flattening_interface::deform_ARAP_polycube_equal_face(TetStructure<double>* tet_mesh_, const std::vector<double>& chart_mean_value_ori, double lambda, bool add_constraint, bool set_min_diff, double offset, bool double_cube_length_flag)
{
    if (add_constraint)
    {
        update_updownchart(tet_mesh_);
    }
    bool find_chart_flag;
    if (double_cube_length_flag)
        cube_len = cube_len * 2;
    if (!hex_meshing_flag)
    {
        if (set_min_diff)
        {
            find_chart_flag = find_all_chart_value_preprocessing(tet_mesh_);
        }
        else
        {
            find_chart_flag = find_all_chart_value_equal_face_mine(tet_mesh_, chart_mean_value_ori, lambda, add_constraint, set_min_diff, offset);
        }
    }
    else
    {
        find_chart_flag = find_all_chart_value_equal_face_mine(tet_mesh_, chart_mean_value_ori, lambda, add_constraint, set_min_diff, offset);
    }
    if (double_cube_length_flag)
        cube_len = cube_len / 2;
    if (!find_chart_flag)
        return false;
    int nv = vertex_type.size();
    is_polycube_handles.clear(); is_polycube_handles.resize(3 * nv, -1);
    std::vector<int> polycube_edge_v; std::vector<int> polycube_edge_v_neighbor;
    std::vector<int> polycube_edge_v_type;
    const std::vector<Tetrahedron<double>*> &tetras = tet_mesh_->tetras;
    std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
    std::vector<int> corner_pts_idx;
    std::vector<std::vector<int>> corner_faces;
    std::map<int, ig::CVec<double, 3>> idx2coord;
    corner_faces.resize(chart_mean_value.size());
    for (int i = 0; i < nv; i++)
    {
        OpenVolumeMesh::Geometry::Vec3i& v_type = vertex_type[i];
        if (v_type[0] >= 0 && v_type[1] >= 0 && v_type[2] >= 0)
        {
            corner_pts_idx.push_back(i);
            CVec<double, 3> np = tetra_vertices[i]->pos;
            idx2coord.insert(std::pair<int, ig::CVec<double, 3>>(i, np));
            corner_faces[v_type[0]].push_back(i);
            corner_faces[v_type[1]].push_back(i);
            corner_faces[v_type[2]].push_back(i);
        }
    }
    for (int i = 0; i < nv; ++i)
    {
        OpenVolumeMesh::Geometry::Vec3i& v_type = vertex_type[i];
        CVec<double, 3> np = tetra_vertices[i]->pos;
        if (v_type[0] >= 0)
        {
            np[0] = chart_mean_value[v_type[0]];
            dpx[i] = np[0];
            is_polycube_handles[3 * i + 0] = 1;
        }
        if (v_type[1] >= 0)
        {
            np[1] = chart_mean_value[v_type[1]];
            dpy[i] = np[1];
            is_polycube_handles[3 * i + 1] = 1;
        }
        if (v_type[2] >= 0)
        {
            np[2] = chart_mean_value[v_type[2]];
            dpz[i] = np[2];
            is_polycube_handles[3 * i + 2] = 1;
        }
        tetra_vertices[i]->pos = np;
        if (v_type[0] >= 0 && v_type[1] >= 0 && v_type[2] < 0)
        {
            polycube_edge_v.push_back(i);
            polycube_edge_v_type.push_back(2);
            is_polycube_handles[3 * i + 2] = 1;
        }
        else if (v_type[0] < 0 && v_type[1] >= 0 && v_type[2] >= 0)
        {
            polycube_edge_v.push_back(i);
            polycube_edge_v_type.push_back(0);
            is_polycube_handles[3 * i + 0] = 1;
        }
        else if (v_type[0] >= 0 && v_type[1] < 0 && v_type[2] >= 0)
        {
            polycube_edge_v.push_back(i);
            polycube_edge_v_type.push_back(1);
            is_polycube_handles[3 * i + 1] = 1;
        }
    }
    bool flag_translation = true;
    if (flag_translation)
    {
        std::vector<ig::CVec<double, 3>> chart_translation;
        for (size_t i = 0; i < chart_mean_value.size(); i++)
        {
            int face_size = corner_faces[i].size();
            ig::CVec<double, 3> tmp_translation(0.0, 0.0, 0.0);
            for (size_t j = 0; j < face_size; j++)
            {
                int tmp_idx = corner_faces[i][j];
                tmp_translation = tmp_translation + tetra_vertices[tmp_idx]->pos - idx2coord[tmp_idx];
            }
            tmp_translation = tmp_translation / (double)face_size;
            int ignore_label = polycube_chart_label[i];
            tmp_translation[ignore_label / 2] = 0.0;
            chart_translation.push_back(tmp_translation);
        }
        for (size_t i = 0; i < nv; i++)
        {
            std::vector<int> tmp_vert_belong_chart;
            OpenVolumeMesh::Geometry::Vec3i& v_type = vertex_type[i];
            if (v_type[0] >= 0)
                tmp_vert_belong_chart.push_back(v_type[0]);
            if (v_type[1] >= 0)
                tmp_vert_belong_chart.push_back(v_type[1]);
            if (v_type[2] >= 0)
                tmp_vert_belong_chart.push_back(v_type[2]);
            if (tmp_vert_belong_chart.size() == 1)
            {
                int tmp_chart = tmp_vert_belong_chart[0];
                CVec<double, 3> np = tetra_vertices[i]->pos;
                np = np + chart_translation[tmp_chart];
                tetra_vertices[i]->pos = np;
            }
        }
    }
    for (int i = 0; i < polycube_edge_v.size(); ++i)
    {
        int v_id = polycube_edge_v[i];
        CVec<double, 3> np = tetra_vertices[v_id]->pos;
        std::vector<int>& one_vv = bvv_id[v_id];
        int vv_size = one_vv.size(); int add_nv = 0;
        for (int j = 0; j < vv_size; ++j)
        {
            int vj_id = one_vv[j];
            OpenVolumeMesh::Geometry::Vec3i& vj_type = vertex_type[vj_id]; int chart_count = 0; bool can_add = false;
            if (vj_type[0] >= 0) ++chart_count;
            if (vj_type[1] >= 0) ++chart_count;
            if (vj_type[2] >= 0) ++chart_count;
            if (chart_count == 3)
            {
                can_add = true;
            }
            else if (chart_count == 2)
            {
                if (vj_type[polycube_edge_v_type[i]] < 0)
                {
                    can_add = true;
                }
            }
            if (can_add)
            {
                CVec<double, 3> npj = tetra_vertices[vj_id]->pos;
                bool nv_ok = true;
                for (int k = 0; k < 3; ++k)
                {
                    if (k != polycube_edge_v_type[i])
                    {
                        if (std::abs(npj[k] - np[k]) > 1e-10) nv_ok = false;
                    }
                }
                if (nv_ok)
                {
                    polycube_edge_v_neighbor.push_back(vj_id);
                    ++add_nv;
                }
            }
        }
        if (add_nv > 0 && add_nv != 2)
        {
            if (add_nv > 2)
            {
                int pevn_size = polycube_edge_v_neighbor.size();
                std::vector<int> delete_v;
                std::vector<int> visited_f(tetras.size(), 0);
                std::vector<int>& one_vf = bvf_id[v_id];
                for (int j = 0; j < one_vf.size(); ++j)
                {
                    visited_f[one_vf[j]] += 1;
                }
                for (int j = pevn_size - add_nv; j < pevn_size; ++j)
                {
                    int vj_id = polycube_edge_v_neighbor[j];
                    std::vector<int>& one_vj_f = bvf_id[vj_id];
                    int two_f[2]; int two_count = 0;
                    for (int k = 0; k < one_vj_f.size(); ++k)
                    {
                        if (visited_f[one_vj_f[k]] == 1)
                        {
                            two_f[two_count] = one_vj_f[k]; ++two_count;
                        }
                    }
                    if (bf_chart[two_f[0]] == bf_chart[two_f[1]])
                    {
                        delete_v.push_back(j);
                    }
                }
                for (int j = 0; j < delete_v.size(); ++j)
                {
                    polycube_edge_v_neighbor.erase(polycube_edge_v_neighbor.begin() + delete_v[delete_v.size() - 1 - j]);
                }
            }
        }
    }
    printf("polycube edge information: \n");
    printf("%d %d\n", polycube_edge_v.size(), polycube_edge_v_neighbor.size());
    if (true)
    {
        if (polycube_edge_v_neighbor.size() == 2 * polycube_edge_v.size())
        {
            for (int i = 0; i < 20; ++i)
            {
                for (int j = 0; j < polycube_edge_v.size(); ++j)
                {
                    int v_id = polycube_edge_v[j];
                    CVec<double, 3> np = tetra_vertices[v_id]->pos;
                    CVec<double, 3> p0 = tetra_vertices[polycube_edge_v_neighbor[2 * j + 0]]->pos;
                    CVec<double, 3> p1 = tetra_vertices[polycube_edge_v_neighbor[2 * j + 1]]->pos;
                    int edge_type = polycube_edge_v_type[j];
                    np[edge_type] = 0.5*(p0[edge_type] + p1[edge_type]);
                    tetra_vertices[v_id]->pos = np;
                    dpx[v_id] = np[0]; dpy[v_id] = np[1]; dpz[v_id] = np[2];
                }
            }
        }
    }
    if (true)
    {
        assert(three_cut_common_vert.size() == three_cut_vert.size());
        for (size_t i = 0; i < three_cut_common_vert.size(); i++)
        {
            int common_vert = three_cut_common_vert[i];
            int n_segment = three_cut_vert[i].size();
            int end_idx[3];
            end_idx[0] = three_cut_vert[i][0][0];
            end_idx[1] = three_cut_vert[i][0][1];
            end_idx[2] = three_cut_vert[i][0][2];
            CVec<double, 3> common_pos = tetra_vertices[common_vert]->pos;
            CVec<double, 3> end_pos[3];
            end_pos[0] = tetra_vertices[end_idx[0]]->pos;
            end_pos[1] = tetra_vertices[end_idx[1]]->pos;
            end_pos[2] = tetra_vertices[end_idx[2]]->pos;
            CVec<double, 3> seg[3];
            seg[0] = (common_pos - end_pos[0]) / n_segment;
            seg[1] = (common_pos - end_pos[1]) / n_segment;
            seg[2] = (common_pos - end_pos[2]) / n_segment;
            for (int j = 1; j < three_cut_vert[i].size(); j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    int v_id = three_cut_vert[i][j][k];
                    CVec<double, 3> np = tetra_vertices[v_id]->pos;
                    np = end_pos[k] + j * seg[k];
                    tetra_vertices[v_id]->pos = np;
                    dpx[v_id] = np[0]; dpy[v_id] = np[1]; dpz[v_id] = np[2];
                }
            }
        }
    }
    boundary_mapping_polycube_equal_face(tet_mesh_);
    if (false)
    {
        is_polycube_handles.clear(); is_polycube_handles.resize(3 * nv, -1);
        for (int i = 0; i < nv; ++i)
        {
            OpenVolumeMesh::Geometry::Vec3i& v_type = vertex_type[i];
            if (v_type[0] >= 0)
            {
                is_polycube_handles[3 * i + 0] = 1;
            }
            if (v_type[1] >= 0)
            {
                is_polycube_handles[3 * i + 1] = 1;
            }
            if (v_type[2] >= 0)
            {
                is_polycube_handles[3 * i + 2] = 1;
            }
        }
    }
    double min_singular_value = 0.001;
    int nc = tetras.size();
    Eigen::Matrix3d VS, Q, U, V, R;
    std::vector < std::vector < double >> Cell_R(nc);
    for (int c_id = 0; c_id < nc; ++c_id)
    {
        Cell_R[c_id].resize(9);
    }
    for (int iter = 0; iter < 3; ++iter)
    {
        for (int c_id = 0; c_id < nc; ++c_id)
        {
            std::vector<int>& cv_id = cell_vertex[c_id];
            double cpx = dpx[cv_id[0]]; double cpy = dpy[cv_id[0]]; double cpz = dpz[cv_id[0]];
            VS(0, 0) = dpx[cv_id[1]] - cpx; VS(1, 0) = dpy[cv_id[1]] - cpy; VS(2, 0) = dpz[cv_id[1]] - cpz;
            VS(0, 1) = dpx[cv_id[2]] - cpx; VS(1, 1) = dpy[cv_id[2]] - cpy; VS(2, 1) = dpz[cv_id[2]] - cpz;
            VS(0, 2) = dpx[cv_id[3]] - cpx; VS(1, 2) = dpy[cv_id[3]] - cpy; VS(2, 2) = dpz[cv_id[3]] - cpz;
            Q = VS * cell_S[c_id];
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(Q, Eigen::ComputeFullU | Eigen::ComputeFullV);
            U = svd.matrixU(); V = svd.matrixV();
            double sig0 = svd.singularValues()[0];
            double sig1 = svd.singularValues()[1];
            double sig2 = svd.singularValues()[2];
            if (Q.determinant() > 0)
            {
                R = V * U.transpose();
            }
            else
            {
                U(0, 2) = -U(0, 2); U(1, 2) = -U(1, 2); U(2, 2) = -U(2, 2);
                R = V * U.transpose();
            }
            Cell_R[c_id][0] = R(0, 0); Cell_R[c_id][1] = R(0, 1); Cell_R[c_id][2] = R(0, 2);
            Cell_R[c_id][3] = R(1, 0); Cell_R[c_id][4] = R(1, 1); Cell_R[c_id][5] = R(1, 2);
            Cell_R[c_id][6] = R(2, 0); Cell_R[c_id][7] = R(2, 1); Cell_R[c_id][8] = R(2, 2);
        }
        for (int i = 0; i < 3; ++i)
        {
            int var_count = 0; std::vector<int> map_v(nv, -1);
            for (int j = 0; j < nv; ++j)
            {
                if (is_polycube_handles[3 * j + i] != 1)
                {
                    map_v[j] = var_count;
                    ++var_count;
                }
            }
            Sparse_Matrix A = new Sparse_Matrix(nc * 3, var_count, NOSYM, CCS, 1);
            std::vector<double> cs(4);
            for (int c_id = 0; c_id < nc; ++c_id)
            {
                double cv = std::sqrt(std::abs(cell_volume[c_id])) * 100;
                std::vector<int>& cv_id = cell_vertex[c_id];
                Eigen::Matrix3d& CS = cell_S[c_id];
                std::vector<double>& CR = Cell_R[c_id];
                for (int j = 0; j < 3; ++j)
                {
                    cs[1] = CS(0, j)*cv; cs[2] = CS(1, j)*cv; cs[3] = CS(2, j)*cv;
                    cs[0] = -(cs[1] + cs[2] + cs[3]);
                    for (int k = 0; k < 4; ++k)
                    {
                        int var_id = map_v[cv_id[k]];
                        if (var_id < 0)
                        {
                            A.fill_rhs_entry(3 * c_id + j, -cs[k] * tetra_vertices[cv_id[k]]->pos[i]);
                        }
                        else
                        {
                            A.fill_entry(3 * c_id + j, var_id, cs[k]);
                        }
                    }
                    A.fill_rhs_entry(3 * c_id + j, cv*CR[i + 3 * j]);
                }
            }
            Sparse_Matrix* D = TransposeTimesSelf(&A, CCS, SYM_LOWER, true);
            solve_by_CHOLMOD(D);
            const std::vector<double>& xyz = D->get_solution();
            for (int j = 0; j < nv; ++j)
            {
                if (is_polycube_handles[3 * j + i] != 1)
                {
                    int var_id = map_v[j];
                    CVec<double, 3> np = tetra_vertices[j]->pos;
                    np[i] = xyz[var_id];
                    tetra_vertices[j]->pos = np;
                }
            }
            delete D;
        }
        for (int j = 0; j < nv; ++j)
        {
            CVec<double, 3> np = tetra_vertices[j]->pos;
            dpx[j] = np[0]; dpy[j] = np[1]; dpz[j] = np[2];
        }
    }
    for (int j = 0; j < nv; ++j)
    {
        CVec<double, 3> np = tetra_vertices[j]->pos;
        dpx[j] = np[0]; dpy[j] = np[1]; dpz[j] = np[2];
    }
    return true;
}
void polycube_flattening_interface::save_polycube_para(const char* filename, TetStructure<double>* tet_mesh_)
{
    const std::vector<Tetrahedron<double>*> &tetras = tet_mesh_->tetras;
    const std::vector<TetVertex<double>*> &tetra_vertices = tet_mesh_->tetra_vertices;
    FILE* f_der = fopen(filename, "w");
    fprintf(f_der, "%d\n", tetras.size() * 12);
    int nc = tetras.size();
    for (size_t i = 0; i < nc; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            CVec<double, 3> p = tetras[i]->vertex[j]->pos;
            CVec<double, 3> q = p / cube_len;
            for (size_t k = 0; k < 3; k++)
            {
                int tmp = (int)q[k];
                if (abs(tmp - q[k]) < ROUNDING_TH)
                    q[k] = (double)tmp;
            }
            fprintf(f_der, "%.19f ", (q[0]));
            fprintf(f_der, "%.19f ", (q[1]));
            fprintf(f_der, "%.19f ", (q[2]));
            fprintf(f_der, "\n");
        }
    }
    fclose(f_der);
}
bool polycube_flattening_interface::compute_eps(VolumeMesh* mesh_)
{
    Eigen::Matrix3d VS, A; Eigen::Vector3d a;
    eps = 1e-3;
    double min_det = 1e30; double avg_det = 0.0;
    int d_count = 0.0;
    Eigen::Matrix3d IS;
    IS << 1.0, -0.57735026918962576450914878050196, -0.40824829046386301636621401245098,
        0, 0, 1.2247448713915890490986420373529,
        0, 1.1547005383792515290182975610039, -0.4082482904638630163662140124509;
    for (OpenVolumeMesh::CellIter c_it = mesh_->cells_begin(); c_it != mesh_->cells_end(); ++c_it)
    {
        int c_id = c_it->idx(); OpenVolumeMesh::CellHandle ch = c_it.cur_handle();
        OpenVolumeMesh::VertexHandle pvh(cell_vertex[c_id][0]);
        OpenVolumeMesh::VertexHandle rvh(cell_vertex_vertex[c_id][0][0]);
        OpenVolumeMesh::VertexHandle svh(cell_vertex_vertex[c_id][0][1]);
        OpenVolumeMesh::VertexHandle tvh(cell_vertex_vertex[c_id][0][2]);
        OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(pvh);
        OpenVolumeMesh::Geometry::Vec3d r = mesh_->vertex(rvh) - p;
        OpenVolumeMesh::Geometry::Vec3d s = mesh_->vertex(svh) - p;
        OpenVolumeMesh::Geometry::Vec3d t = mesh_->vertex(tvh) - p;
        VS(0, 0) = r[0]; VS(1, 0) = r[1]; VS(2, 0) = r[2];
        VS(0, 1) = s[0]; VS(1, 1) = s[1]; VS(2, 1) = s[2];
        VS(0, 2) = t[0]; VS(1, 2) = t[1]; VS(2, 2) = t[2];
        A = VS * S[c_id][0];
        double det = A.determinant();
        if (det < min_det) min_det = det;
        avg_det += std::abs(det);
        if (det < 0) ++d_count;
    }
    avg_det /= mesh_->n_cells();
    double umbral = umbral_factor * avg_det;
    double radicando = EpsilonEfectivo * (EpsilonEfectivo - min_det);
    if (radicando <= 0.0)
        min_det = 0.0;
    else
        min_det = std::sqrt(radicando);
    eps = (min_det >= umbral) ? min_det : umbral;
    eps = eps * eps;
    printf("Flip Count : %d, %e\n", d_count, eps);
    return (d_count == 0);
}

void polycube_flattening_interface::set_bd_face_chart_label(const std::vector<CVec<double, 3>> &bd_pts_data, const std::vector<unsigned> &bd_faces_data, const std::vector<int> &bd_chart_data, const std::vector<int> &bd_label_data)
{
    bd_pts = bd_pts_data;
    bd_faces = bd_faces_data;
    bd_chart = bd_chart_data;
    bd_label = bd_label_data;
}
void polycube_flattening_interface::get_bd_face_chart_label(TetStructure<double> *tet_mesh_)
{
    if (!bd_pts.empty()) return;
    assert(!polycube_edges.empty());
    bd_pts.clear();
    bd_faces.clear();
    bd_chart.clear();
    bd_label.clear();
    bool cur_hex_flag = hex_meshing_flag;
    set_hex_meshing_flag(false);
    bool find_chart_flag = find_all_chart_value_equal_face_mine(tet_mesh_);
    set_hex_meshing_flag(cur_hex_flag);
    std::map<int, CVec<double, 3>> id2cornerpt;
    for (int i = 0; i < tet_mesh_->tetra_vertices.size(); i++)
    {
        OpenVolumeMesh::Geometry::Vec3i& v_type = vertex_type[i];
        if (v_type[0] >= 0 && v_type[1] >= 0 && v_type[2] >= 0)
        {
            CVec<double, 3> pt(chart_mean_value[v_type[0]], chart_mean_value[v_type[1]], chart_mean_value[v_type[2]]);
            id2cornerpt[i] = pt;
            dpx[i] = chart_mean_value[v_type[0]];
            dpy[i] = chart_mean_value[v_type[1]];
            dpz[i] = chart_mean_value[v_type[2]];
        }
    }
    std::map<std::pair<int, int>, int> pq_edge2idx;
    for (int i = 0; i < polycube_edges.size(); i++)
    {
        std::pair<int, int> point_pair1(polycube_edges[i][0], polycube_edges[i][1]);
        std::pair<int, int> point_pair2(polycube_edges[i][1], polycube_edges[i][0]);
        pq_edge2idx[point_pair1] = i;
        pq_edge2idx[point_pair2] = i;
    }
    std::vector<std::vector<std::pair<int, int>>> edges_with_same_chart;
    int max_chart_idx = 0;
    for (size_t i = 0; i < polycube_edges.size(); i++)
    {
        if (max_chart_idx < polycube_edges[i][2])
            max_chart_idx = polycube_edges[i][2];
        if (max_chart_idx < polycube_edges[i][3])
            max_chart_idx = polycube_edges[i][3];
    }
    edges_with_same_chart.resize(max_chart_idx + 1);
    for (size_t i = 0; i < polycube_edges.size(); i++)
    {
        int chart_idx1 = polycube_edges[i][2];
        int chart_idx2 = polycube_edges[i][3];
        std::pair<int, int> point_pair(polycube_edges[i][0], polycube_edges[i][1]);
        edges_with_same_chart[chart_idx1].push_back(point_pair);
        edges_with_same_chart[chart_idx2].push_back(point_pair);
    }
    std::vector<std::vector<std::vector<int>>> chart2loops;
    for (size_t i = 0; i < edges_with_same_chart.size(); i++)
    {
        std::vector<std::pair<int, int>> &one_chart = edges_with_same_chart[i];
        std::list<std::pair<int, int>> one_chart_list(one_chart.begin(), one_chart.end());
        std::vector<int> one_loop;
        std::vector<std::vector<int>> onechart_loops;
        int next_idx;
        while (!one_chart_list.empty())
        {
            one_loop.clear();
            one_loop.push_back(one_chart_list.front().first);
            one_loop.push_back(one_chart_list.front().second);
            one_chart_list.pop_front();
            next_idx = -1;
            int max_try_time = one_chart_list.size();
            int try_time = 0;
            while (!one_chart_list.empty() && try_time < max_try_time && next_idx != one_loop[0])
            {
                for (std::list<std::pair<int, int>>::iterator it = one_chart_list.begin(); it != one_chart_list.end(); ++it)
                {
                    if ((*it).first == one_loop.back())
                    {
                        next_idx = (*it).second;
                        one_chart_list.erase(it);
                        break;
                    }
                    if ((*it).second == one_loop.back())
                    {
                        next_idx = (*it).first;
                        one_chart_list.erase(it);
                        break;
                    }
                }
                try_time++;
                if (next_idx != one_loop.back() && next_idx != -1)
                    one_loop.push_back(next_idx);
            }
            assert(one_loop[0] == one_loop.back());
            if (one_loop[0] == one_loop.back())
            {
                onechart_loops.push_back(one_loop);
            }
        }
        chart2loops.push_back(onechart_loops);
    }
    for (size_t i = 0; i < chart2loops.size(); i++)
    {
        std::vector<double> face_coord_x, face_coord_y, face_coord_z;
        for (size_t k = 0; k < chart2loops[i].size(); k++)
        {
            chart2loops[i][k].erase(chart2loops[i][k].begin() + chart2loops[i][k].size() - 1);
            int face_size = chart2loops[i][k].size();
            for (size_t j = 0; j < face_size; j++)
            {
                auto it = id2cornerpt.find(chart2loops[i][k][j]);
                assert(it != id2cornerpt.end());
                face_coord_x.push_back(id2cornerpt[chart2loops[i][k][j]][0]);
                face_coord_y.push_back(id2cornerpt[chart2loops[i][k][j]][1]);
                face_coord_z.push_back(id2cornerpt[chart2loops[i][k][j]][2]);
            }
        }
        std::vector<std::vector<int>> new_faces;
        if (ig::SimpleTriangulation::sortface_area(face_coord_x, face_coord_y, face_coord_z, chart2loops[i]))
        {
            face_coord_x.clear();
            face_coord_y.clear();
            face_coord_z.clear();
            for (size_t k = 0; k < chart2loops[i].size(); k++)
            {
                int face_size = chart2loops[i][k].size();
                for (size_t j = 0; j < face_size; j++)
                {
                    face_coord_x.push_back(id2cornerpt[chart2loops[i][k][j]][0]);
                    face_coord_y.push_back(id2cornerpt[chart2loops[i][k][j]][1]);
                    face_coord_z.push_back(id2cornerpt[chart2loops[i][k][j]][2]);
                }
            }
        }
        bool triangulation_success = ig::SimpleTriangulation::triangulation(face_coord_x, face_coord_y, face_coord_z, chart2loops[i], new_faces);
        for (size_t j = 0; j < new_faces.size(); j++)
        {
            for (size_t k = 0; k < 3; k++)
            {
                bd_faces.push_back(new_faces[j][k]);
            }
            bd_chart.push_back(i);
            bd_label.push_back(polycube_chart_label[i]);
        }
    }
    std::map<int, int> ido2n;
    int count = 0;
    for (auto p : id2cornerpt)
    {
        bd_pts.push_back(p.second);
        ido2n[p.first] = count++;
    }
    for (size_t i = 0; i < bd_faces.size(); i++)
    {
        bd_faces[i] = ido2n[bd_faces[i]];
    }
}
void polycube_flattening_interface::load_deformation_result(const std::vector<double> &coord_ori, TetStructure<double> *tet_mesh_)
{
    int nv = tet_mesh_->tetra_vertices.size();
    dpx_ori.resize(nv);
    dpy_ori.resize(nv);
    dpz_ori.resize(nv);
    for (size_t i = 0; i < nv; i++)
    {
        dpx_ori[i] = coord_ori[3 * i];
        dpy_ori[i] = coord_ori[3 * i + 1];
        dpz_ori[i] = coord_ori[3 * i + 2];
    }
    prepare_for_deformation(tet_mesh_, dpx_ori, dpy_ori, dpz_ori);
    printf("Loading deformation result.............\n");
    build_AABB_Tree();
#if 1
    for (int i = 0; i < deformation_v_id.size(); ++i)
    {
        deformation_new_p[i][0] = dpx[deformation_v_id[i]];
        deformation_new_p[i][1] = dpy[deformation_v_id[i]];
        deformation_new_p[i][2] = dpz[deformation_v_id[i]];
    }
#endif
    compute_distortion(tet_mesh_);
}
void polycube_flattening_interface::set_equal_faces(const std::vector<std::vector<std::vector<int>>> &faces_array, const std::vector<int>& common_verts_idx, const std::vector<int>& cut_types_array, int n_vert, const std::vector<std::pair<int, int>> &chart_pair, const std::vector<std::pair<int, int>> &chart_pair_neighbor, const std::vector<int> &three_cut_common_vert_array, const std::vector<std::vector<std::array<unsigned int, 3>>> &three_cut_vert_array, const std::vector<std::array<int, 3>> &three_cut_adjacent_one_cut_index_array)
{
    assert(faces_array.size() == (common_verts_idx.size() >> 1) && faces_array.size() == cut_types_array.size());
    if (faces_array.size() == 0)
        return;
    vert_update_type.clear();
    vert_belong_cut_array.clear();
    cut_types.clear();
    cut_common_verts_idx.clear();
    cut_to_chart_pair.clear();
    vert_pairs_map.clear();
    equal_triangles.clear();
    vert_update_type.resize(n_vert, 0);
    vert_belong_cut_array.resize(n_vert, -1);
    cut_types = cut_types_array;
    for (size_t i = 0; i < common_verts_idx.size() / 2; i++)
    {
        cut_common_verts_idx.push_back(common_verts_idx[2 * i]);
    }
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
    cut_to_chart_pair = chart_pair;
    cut_to_chart_pair_neighbor = chart_pair_neighbor;
    three_cut_common_vert = three_cut_common_vert_array;
    three_cut_vert = three_cut_vert_array;
    three_cut_adjacent_one_cut_index = three_cut_adjacent_one_cut_index_array;
    three_cut_vert_flag.clear();
    three_cut_vert_flag.resize(n_vert, -1);
    assert(three_cut_common_vert.size() == three_cut_vert.size());
    for (size_t i = 0; i < three_cut_common_vert.size(); i++)
    {
        three_cut_vert_flag[three_cut_common_vert[i]] = 1;
        for (size_t j = 0; j < three_cut_vert[i].size(); j++)
        {
            for (size_t k = 0; k < 3; k++)
            {
                three_cut_vert_flag[three_cut_vert[i][j][k]] = 1;
            }
        }
    }
}
