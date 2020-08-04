#pragma once
#include <vector>

namespace SUSPolycube 
{
	//bool SUS(const std::vector<double> &nodes_coord, const std::vector<int> &ref_nodes, const std::vector<unsigned int> &elems, std::vector<int> &node_normal_tag, const std::vector<double> &node_deform, std::vector<double> &output_coord, int max_smooth_iter = 20, int max_iter = 40);
	//# of nodes_coord = 3 * # of nodes of ref_nodes
	bool SUS(const std::vector<double> &nodes_deform, const std::vector<int> &ref_nodes, const std::vector<unsigned int> &elems, std::vector<int> &node_normal_tag, const std::vector<double> &node_ori, std::vector<double> &output_coord, int max_smooth_iter = 20, int max_iter = 40);
	//bool test();
}

