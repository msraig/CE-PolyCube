#pragma once
#include <vector>
#include <set>

void greedy_graph_coloring(const size_t num_vertices, const std::vector<std::set<size_t>>& edges, std::vector<std::vector<size_t>>& colored_vertices);