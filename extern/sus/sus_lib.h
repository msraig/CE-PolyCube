#pragma once
#include <vector>

#define STANDARD_TET 0

//bool SUS_standard(const std::vector<double> &nodes_coord, const std::vector<int> &ref_nodes, const std::vector<unsigned int> &elems, std::vector<double> &output_coord, int max_smooth_iter = 20, int max_iter = 40);

bool SUS(const std::vector<double> &nodes_coord, const std::vector<int> &ref_nodes, const std::vector<unsigned int> &elems,std::vector<int> &node_normal_tag, const std::vector<double> &node_deform, std::vector<double> &output_coord, int max_smooth_iter = 20, int max_iter = 40);
//# of nodes_coord = 3 * # of nodes of ref_nodes

//bool test();
bool SUS_proj(const std::vector<double> &nodes_coord, const std::vector<int> &ref_nodes, const std::vector<unsigned int> &elems, std::vector<int> &node_normal_tag, const std::vector<double> &node_deform, std::vector<double> &proj_pts, const std::vector<int> &proj_faces, std::vector<double> &output_coord, int max_smooth_iter = 20, int max_iter = 40);
//only normal_x is useful to determine whether a vert is on boundary and should be optimize;


bool SUS_proj_laplacian(const std::vector<double> &nodes_coord, const std::vector<int> &ref_nodes, const std::vector<unsigned int> &elems, std::vector<int> &node_normal_tag, const std::vector<double> &node_deform, std::vector<double> &proj_pts, const std::vector<int> &proj_faces, std::vector<double> &output_coord, const std::vector<int> &boundary_idx, const std::vector<int> &boundary_faces, int max_smooth_iter = 20, int max_iter = 40);