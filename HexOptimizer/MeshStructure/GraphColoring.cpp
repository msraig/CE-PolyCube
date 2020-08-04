#include "GraphColoring.h"
#include <algorithm>

void greedy_graph_coloring(const size_t num_vertices, const std::vector<std::set<size_t>>& edges, std::vector<std::vector<size_t>>& colored_vertices)
{
	std::vector<int> result(num_vertices, -1);
	result[0] = 0;
	std::vector<bool> available(num_vertices, false);

	int max_color = 0;
	for (size_t u = 1; u < num_vertices; u++)
	{
		for (typename std::set<size_t>::iterator iter = edges[u].begin(); iter != edges[u].end(); iter++)
		{
			if (result[*iter] != -1)
				available[result[*iter]] = true;
		}
		int cr;
		for (cr = 0; cr < (int)num_vertices; cr++)
		{
			if (available[cr] == false)
				break;
		}

		result[u] = cr;
		max_color = std::max(cr, max_color);

		for (typename std::set<size_t>::iterator iter = edges[u].begin(); iter != edges[u].end(); iter++)
		{
			if (result[*iter] != -1)
				available[result[*iter]] = false;
		}
	}
	colored_vertices.resize(max_color + 1);
	for (size_t i = 0; i < num_vertices; i++)
	{
		colored_vertices[result[i]].push_back(i);
	}

}