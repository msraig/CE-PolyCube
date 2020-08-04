#pragma once

#include <vector>

//bool SUS_standard(std::vector<double> &nodes_coord, std::vector<int> &ref_nodes, std::vector<int> &elems, std::vector<double> &output_coord, int max_smooth_iter = 20, int max_iter = 40);

namespace SUSSpace {

	bool SUS_standard(const std::vector<double> &nodes_deform, const std::vector<int> &ref_nodes, const std::vector<int> &elems, std::vector<double> &output_coord, int max_smooth_iter = 20, int max_iter = 40);

	bool SUS_standard_polycube(const std::vector<double> &nodes_deform, const std::vector<int> &ref_nodes, const std::vector<int> &elems, const std::vector<int> &node_normal_tag, std::vector<double> &output_coord, int max_smooth_iter = 20, int max_iter = 40);
}



