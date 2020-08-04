#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <queue>
#include <assert.h>
#include "HexOperation.h"

#include <boost/config.hpp>
#include <iostream>
#include <fstream>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

using namespace boost;
using namespace HexStructure;
using ig::CVec;

bool load_tet_feature_vtk(const char* filename, std::vector<CVec<double, 3>> &tet_fea_point, std::vector<std::pair<int, int>> &tet_fea_edge)
{
	tet_fea_point.clear();
	tet_fea_edge.clear();
	std::ifstream ifs(filename);
	if (!ifs.is_open()) return false;
	
	bool find = false; int lines = 0;
	int vnum, edgenum;
	char s[1024], sread[1024], sread2[1024];
	while (!find)
	{
		ifs.getline(s, 1013);
		std::cout << s << std::endl;
		if (sscanf(s, "%s %d %s", &sread, &vnum, &sread2) == 3 && (strcmp(sread, "POINTS") == 0))
			find = true;
		if (++lines > 10)
		{
			throw std::runtime_error("cannot find head of VTK!");
			return false;
		}
	}
	
	tet_fea_point.resize(vnum);
	double x, y, z;
	for (size_t i = 0; i < vnum; i++)
	{
		ifs >> x >> y >> z;
		tet_fea_point[i][0] = x;
		tet_fea_point[i][1] = y;
		tet_fea_point[i][2] = z;
	}

	find = false;
	while (!find)
	{
		int tmp_int;
		ifs.getline(s, 1023);
		std::cout << s << std::endl;
		if (sscanf(s, "%s %d %d", &sread, &edgenum, &tmp_int) == 3 && (strcmp(sread, "LINES") == 0))
			find = true;
	}
	tet_fea_edge.resize(edgenum);
	for (size_t i = 0; i < edgenum; i++)
	{
		int tmp;
		ifs >> tmp >> tet_fea_edge[i].first >> tet_fea_edge[i].second;
	}

	ifs.close();
	return true;
}

bool load_tet_feature_psfe(const char* filename, std::vector<CVec<double, 3>> &tet_fea_point, std::vector<std::pair<int, int>> &tet_fea_edge, std::vector<std::pair<double, double>> &tet_fea_segm)
{
	tet_fea_point.clear();
	tet_fea_edge.clear();
	tet_fea_segm.clear();
	std::ifstream ifs(filename);
	if (!ifs.is_open()) return false;

	bool find = false; int lines = 0;
	int vnum, edgenum;
	char s[1024], sread[1024], sread2[1024];
	std::string str;
	ifs >> str;
	ifs >> vnum;
	std::cout << "Point Size: " << vnum << std::endl;

	tet_fea_point.resize(vnum);
	double x, y, z;
	for (size_t i = 0; i < vnum; i++)
	{
		ifs >> x >> y >> z;
		tet_fea_point[i][0] = x;
		tet_fea_point[i][1] = y;
		tet_fea_point[i][2] = z;
	}

	ifs >> str >> edgenum;
	std::cout << "Edge Size: " << edgenum << std::endl;
	
	tet_fea_edge.resize(edgenum);
	tet_fea_segm.resize(edgenum);
	for (size_t i = 0; i < edgenum; i++)
	{
		ifs >> tet_fea_edge[i].first >> tet_fea_edge[i].second >> tet_fea_segm[i].first >> tet_fea_segm[i].second;		
	}

	ifs.close();
	return true;
}

void write_hex_feature_vtk(const char* filename, const HexFrame &hf, const std::vector<int> &edge_array)
{
	std::ofstream ofs(filename);
	ofs << "# vtk DataFile Version 3.0\n"
		<< "mesh vtk data\n"
		<< "ASCII\n"
		<< "DATASET POLYDATA\n";
	ofs << "POINTS " << hf.mesh.Vs.size() << " double" << std::endl;

	for (size_t i = 0; i < hf.mesh.Vs.size(); i++)
	{
		ofs << hf.mesh.Vs[i].v << std::endl;
	}

	std::vector<int> feature_point_color(hf.mesh.Vs.size(), 0);
	ofs << "LINES " << edge_array.size() << " " << 3 * edge_array.size() << std::endl;
	for (size_t i = 0; i < edge_array.size(); i++)
	{
		int eid = edge_array[i];
		int id0 = hf.mesh.Es[eid].vs[0];
		int id1 = hf.mesh.Es[eid].vs[1];
		feature_point_color[id0] = 1;
		feature_point_color[id1] = 1;
		ofs << "2 " << id0 << " " << id1 << std::endl;
	}
	ofs << "POINT_DATA " << hf.mesh.Vs.size() << "\n"
		<< "SCALARS V_Scalars int\nLOOKUP_TABLE V_Table" << std::endl;

	for (size_t i = 0; i < hf.mesh.Vs.size(); i++)
	{
		ofs << feature_point_color[i] << std::endl;
	}

	ofs.close();
}

void write_hex_group_feature_feaformat(const char* filename, const std::vector<std::vector<std::pair<int, int>>> &edge_array_pair)
{
	//first not group then grouped info
	//grouped information does not contain merge information
	if (edge_array_pair.empty()) return;
	std::vector<std::pair<int, int>> all_short_edges;
	for (size_t i = 0; i < edge_array_pair.size(); i++)
	{
		all_short_edges.insert(all_short_edges.end(), edge_array_pair[i].begin(), edge_array_pair[i].end());
	}

	std::ofstream ofs(filename);
	if (!ofs.is_open()) return;
	ofs << all_short_edges.size() << std::endl;
	for (size_t i = 0; i < all_short_edges.size(); i++)
	{
		ofs << all_short_edges[i].first << " " << all_short_edges[i].second << std::endl;
	}

	ofs << "Grouped Hex Feature" << std::endl;
	ofs << edge_array_pair.size() << std::endl;
	for (size_t i = 0; i < edge_array_pair.size(); i++)
	{
		int one_edge_size = (int)edge_array_pair[i].size();
		ofs << one_edge_size << std::endl;
		for (size_t j = 0; j < one_edge_size; j++)
		{
			ofs << edge_array_pair[i][j].first << " " << edge_array_pair[i][j].second << std::endl;
		}
	}
	ofs.close();
}

void write_hex_group_feature_hfeformat(const char* filename, const HexFrame &hf, const std::vector<std::vector<std::pair<int, int>>> &edge_array_pair)
{
	std::ofstream ofs(filename);
	if (!ofs.is_open()) return;
	ofs << "Hex Feature Line" << std::endl;
	ofs << edge_array_pair.size() << std::endl;
	for (size_t i = 0; i < edge_array_pair.size(); i++)
	{
		int one_edge_size = (int)edge_array_pair[i].size();
		ofs << one_edge_size + 1 << std::endl;
		int vid = edge_array_pair[i][0].first;
		
		ofs << hf.mesh.Vs[vid].v << " ";

		for (size_t j = 0; j < one_edge_size; j++)
		{
			int vid = edge_array_pair[i][j].second;
			ofs << hf.mesh.Vs[vid].v << " ";
		}
		ofs << std::endl;
		
	}
	ofs.close();
}


void reorder_edge(std::vector<std::pair<int, int>> &one_edge, int start)
{
	std::map<int, std::vector<int>> vert2edgeid;
	for (int i = 0; i < one_edge.size(); i++)
	{
		vert2edgeid[one_edge[i].first].push_back(i);
		vert2edgeid[one_edge[i].second].push_back(i);
	}

	if (start == -1)
	{
		//no start
		for (auto vepair : vert2edgeid)
		{
			if (vepair.second.size() == 1)
			{
				start = vepair.first;
				break;
			}
		}
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
		//find next
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

int main(int argc, char** argv)
{
	if (argc < 5)
	{
		std::cout << "Usage1: " << argv[0] << " : InputPolycubeHex.vtk InputTetFeature.vtk InputOriHex.vtk OutputHexFeature.hfe" << std::endl;
		std::cout << "Usage2: " << argv[0] << " : InputCutPolycubeHex.vtk InputCutTetFeature.vtk InputCutOriHex.vtk  OutputHexFeature.hfe InputCut.pfe" << std::endl;
		return 0; 
	}
	//merging need to be take into consideration later, another .fea file is needed
	
	std::cout << "Input: ";
	for (size_t i = 0; i < argc; i++)
	{
		std::cout << argv[i] << " ";
	}
	std::cout << std::endl;

	HexFrame hf_polycube;
	HexFrame hf_ori;

	if (!read_hex_mesh_vtk(hf_polycube.mesh, argv[1]))
	{
		std::cout << "!!!InputPolycubeHex.vtk not FOUND!!!" << std::endl;
		return 0;
	}

	if (!read_hex_mesh_vtk(hf_ori.mesh, argv[3]))
	{
		std::cout << "!!!InputOriHex.vtk not FOUND!!!" << std::endl;
		return 0;
	}
	
	std::string outputname(argv[4]);
		
	assert(hf_ori.mesh.Vs.size() == hf_polycube.mesh.Vs.size());
	if (hf_ori.mesh.Vs.size() != hf_polycube.mesh.Vs.size())
	{
		std::cout << "Origin hex And PolyCube hex have different vert num." << std::endl;
		return 0;
	}

	build_connectivity(hf_ori.mesh);
	build_connectivity(hf_polycube.mesh);
	extract_singularity(hf_polycube.mesh, hf_polycube.singularity);
	//extract_metamesh(hf_polycube.mesh, hf_polycube.singularity, hf_polycube.meta_mesh);
	//write_metamesh_vtk(hf_polycube.mesh, hf_polycube.meta_mesh, (char *)"metamesh.vtk");

	//load tet feature
	std::vector<CVec<double, 3>> tet_fea_point;
	std::vector<std::pair<int, int>> tet_fea_edge;
	std::vector<std::pair<double, double>> tet_fea_segm;
	std::string input_fea_file(argv[2]);
	std::string input_fea_suffix = input_fea_file.substr(input_fea_file.find_last_of('.') + 1, input_fea_file.length() - input_fea_file.find_last_of('.') - 1);
	std::cout << "Input Fea Suffix: " << input_fea_suffix << std::endl;

	if (input_fea_suffix == "vtk")
	{
		if (!load_tet_feature_vtk(argv[2], tet_fea_point, tet_fea_edge))
		{
			std::cout << "!!!InputTetFeature.vtk not FOUND!!!" << std::endl;
			return 0;
		}
		tet_fea_segm.resize(tet_fea_edge.size(), std::pair<double, double>(0.0, 1.0));
	}
	else if (input_fea_suffix == "psfe")
	{
		if (!load_tet_feature_psfe(argv[2], tet_fea_point, tet_fea_edge, tet_fea_segm))
		{
			std::cout << "!!!InputTetFeature.psfe not FOUND!!!" << std::endl;
			return 0;
		}
	}
	
	

	std::set<int> corner_eid_set; //edge neighboring with 1 or 3 hex, in hex space
	std::vector<std::vector<int>> corner_v2e(hf_polycube.mesh.Vs.size()); //if size equals to 3, then it is corner
	std::vector<std::vector<int>> grouped_corner_eid;
	std::vector<std::vector<std::pair<int, int>>> grouped_corner_edge_pairs;
	int nhe = hf_polycube.mesh.Es.size();
	for (int i = 0; i < nhe; i++)
	{
		if (hf_polycube.mesh.Es[i].n_cs.size() == 1 || hf_polycube.mesh.Es[i].n_cs.size() == 3)
		{
			//corner e
			corner_eid_set.insert(i);
			int v0 = hf_polycube.mesh.Es[i].vs[0];
			int v1 = hf_polycube.mesh.Es[i].vs[1];
			corner_v2e[v0].push_back(i);
			corner_v2e[v1].push_back(i);
		}
	}
	
	//construct grouped_corner_eid
	while (!corner_eid_set.empty())
	{
		std::vector<int> one_long_edge_id;
		auto first = corner_eid_set.begin();
		std::vector<bool> color(hf_polycube.mesh.Es.size(), false);
		color[*first] = true;
		std::queue<int> q;
		q.push(*first);
		while (!q.empty())
		{
			int front = q.front();
			q.pop();
			one_long_edge_id.push_back(front);
			int v[2] = { hf_polycube.mesh.Es[front].vs[0], hf_polycube.mesh.Es[front].vs[1] };
			for (size_t i = 0; i < 2; i++)
			{
				int vid = v[i];
				if (corner_v2e[vid].size() >= 3) continue;
				assert(corner_v2e[vid].size() == 2);
				int otherid = corner_v2e[vid][0];
				if (otherid == front) otherid = corner_v2e[vid][1];
				if (color[otherid] == false)
				{
					color[otherid] = true;
					q.push(otherid);
				}
			}
		}

		for (size_t i = 0; i < one_long_edge_id.size(); i++)
		{
			int eid = one_long_edge_id[i];
			corner_eid_set.erase(eid);
		}
		grouped_corner_eid.push_back(one_long_edge_id);

	}
	
	std::vector<std::pair<int, int>> grouped_corner_epid;
	for (size_t i = 0; i < grouped_corner_eid.size(); i++)
	{
		std::vector<std::pair<int, int>> one_long_edge_pairs;
		std::map<std::pair<int, int>, int> local_edge_pair2id; //include reverse case
		for (size_t j = 0; j < grouped_corner_eid[i].size(); j++)
		{
			int eid = grouped_corner_eid[i][j];
			std::pair<int, int> tmppair(-1, -1);
			tmppair.first = hf_polycube.mesh.Es[eid].vs[0];
			tmppair.second = hf_polycube.mesh.Es[eid].vs[1];
			local_edge_pair2id[tmppair] = eid;
			one_long_edge_pairs.push_back(tmppair);
			std::swap(tmppair.first, tmppair.second);
			local_edge_pair2id[tmppair] = eid;
		}

		reorder_edge(one_long_edge_pairs, -1);
		grouped_corner_edge_pairs.push_back(one_long_edge_pairs);
		//write_hex_feature_vtk((std::to_string(i) + ".vtk").c_str(), hf_polycube, grouped_corner_eid[i]);
		grouped_corner_epid.push_back(std::pair<int, int>(one_long_edge_pairs.front().first, one_long_edge_pairs.back().second));
		
		//refine grouped_corner_eid[i]
		std::vector<int> new_edge_id;
		for (size_t j = 0; j < one_long_edge_pairs.size(); j++)
		{
			auto it = local_edge_pair2id.find(one_long_edge_pairs[j]);
			assert(it != local_edge_pair2id.end());
			new_edge_id.push_back(it->second);
		}
		grouped_corner_eid[i] = new_edge_id;
	}

	std::vector<std::vector<std::pair<int, int>>> grouped_short_edge_pair;
	std::vector<int> all_short_edge_id;
	std::set<int> used_corner_eid;
	std::set<std::pair<int, int>> tet_fea_edge_unique;
	for (size_t i = 0; i < tet_fea_edge.size(); i++)
	{
		int first = tet_fea_edge[i].first;
		int second = tet_fea_edge[i].second;
		std::pair<int, int> tmppair = tet_fea_edge[i];
		if (tmppair.first > tmppair.second) std::swap(tmppair.first, tmppair.second);
		tet_fea_edge_unique.insert(tmppair);

		int minid = -1;
		double mindist = -1.0;
		bool reverse_flag = true;
		CVec<double, 3> p0 = tet_fea_point[first];
		CVec<double, 3> p1 = tet_fea_point[second];

		for (int j = 0; j < grouped_corner_epid.size(); j++)
		{
			CVec<double, 3> ep0 = hf_polycube.mesh.Vs[grouped_corner_epid[j].first].v;
			CVec<double, 3> ep1 = hf_polycube.mesh.Vs[grouped_corner_epid[j].second].v;
			double dist = (p0 - ep0).Length() + (p1 - ep1).Length();
			if (mindist < 0.0 || mindist > dist)
			{
				mindist = dist;
				minid = j;
				reverse_flag = false;
			}
			
			dist = (p0 - ep1).Length() + (p1 - ep0).Length();
			if (mindist < 0.0 || mindist > dist)
			{
				mindist = dist;
				minid = j;
				reverse_flag = true;
			}
		}
		
		//std::cout << "min dist: " << mindist << std::endl;
		//std::cout << "min id: " << minid << std::endl;
		assert(minid != -1);

		auto findit = used_corner_eid.find(minid);
		if (findit != used_corner_eid.end())
		{
			//found duplicate one
			std::cout << "two point: " << first << " " << second << std::endl;
			std::cout << "!!!!!!!!Duplicate e id: " << minid << std::endl;
		}


		used_corner_eid.insert(minid);
		int segm_begin_id = 0, segm_end_id = grouped_corner_eid[minid].size();
		segm_begin_id = grouped_corner_eid[minid].size() * tet_fea_segm[i].first;
		segm_end_id = grouped_corner_eid[minid].size() * tet_fea_segm[i].second;

		//all_short_edge_id.insert(all_short_edge_id.end(), grouped_corner_eid[minid].begin(), grouped_corner_eid[minid].end());
		if (reverse_flag)
		{
			//reverse
			std::vector<std::pair<int, int>> reverse_one_edge;
			for (size_t j = 0; j < grouped_corner_edge_pairs[minid].size(); j++)
			{
				auto tmppair = grouped_corner_edge_pairs[minid][grouped_corner_edge_pairs[minid].size() - 1 - j];
				std::swap(tmppair.first, tmppair.second);
				reverse_one_edge.push_back(tmppair);
			}
			std::reverse(grouped_corner_eid[minid].begin(), grouped_corner_eid[minid].end());
			all_short_edge_id.insert(all_short_edge_id.end(), grouped_corner_eid[minid].begin() + segm_begin_id, grouped_corner_eid[minid].begin() + segm_end_id);

			std::vector<std::pair<int, int>> edge_segm(reverse_one_edge.begin() + segm_begin_id, reverse_one_edge.begin() + segm_end_id);
			grouped_short_edge_pair.push_back(edge_segm);
		}
		else
		{
			all_short_edge_id.insert(all_short_edge_id.end(), grouped_corner_eid[minid].begin() + segm_begin_id, grouped_corner_eid[minid].begin() + segm_end_id);
			std::vector<std::pair<int, int>> edge_segm(grouped_corner_edge_pairs[minid].begin() + segm_begin_id, grouped_corner_edge_pairs[minid].begin() + segm_end_id);
			grouped_short_edge_pair.push_back(edge_segm);
			//grouped_short_edge_pair.push_back(grouped_corner_edge_pairs[minid]);
		}
	}
	//there might be cases that both feature correspondents to the same polycube hex edge
	//assert(used_corner_eid.size() == tet_fea_edge.size());
	assert(used_corner_eid.size() == tet_fea_edge_unique.size());

	if (argc == 6)
	{
		std::ifstream input_fea_file(argv[5]);
		//redefine grouped_short_edges, all_short_edges not changed
		
		//first part: not group
		std::string tmp_str;
		std::getline(input_fea_file, tmp_str);
		if (tmp_str == "Grouped PolyCube Features")
		{
			std::vector<std::vector<std::pair<int, int>>> grouped_tet_feature;
			int n_group = 0;
			input_fea_file >> n_group;
			std::cout << "Group Size: " << n_group << std::endl;
			for (size_t i = 0; i < n_group; i++)
			{
				int one_group_size = 0;
				input_fea_file >> one_group_size;
				std::vector<std::pair<int, int>> one_tet_feature;
				for (size_t j = 0; j < one_group_size; j++)
				{
					std::pair<int, int> tmp_pair;
					input_fea_file >> tmp_pair.first >> tmp_pair.second;
					one_tet_feature.push_back(tmp_pair);
				}
				assert(one_tet_feature.size() == 1 || one_tet_feature.size() == 2);
				grouped_tet_feature.push_back(one_tet_feature);
			}

			assert(grouped_tet_feature.size() == n_group);
			//redefine grouped_short_edges
			std::vector<std::vector<std::pair<int, int>>> grouped_short_edges_new;
			for (size_t i = 0; i < n_group; i++)
			{
				int findid = -1;
				bool sameorder = true;
				std::pair<int, int> tmp_pair = grouped_tet_feature[i][0];
				for (int j = 0; j < tet_fea_edge.size(); j++)
				{
					if (tmp_pair.first == tet_fea_edge[j].first && tmp_pair.second == tet_fea_edge[j].second)
					{
						findid = j;
						break;
					}
					if (tmp_pair.second == tet_fea_edge[j].first && tmp_pair.first == tet_fea_edge[j].second)
					{
						findid = j;
						sameorder = false;
						break;
					}
				}

				assert(findid != -1);
				if (sameorder)
				{
					grouped_short_edges_new.push_back(grouped_short_edge_pair[findid]);
				}
				else
				{
					//reverse grouped_short_edges[findid]
					std::vector<std::pair<int, int>> reverse_feature;

					for (int j = grouped_short_edge_pair[findid].size() - 1; j >= 0 ; j--)
					{
						std::pair<int, int> tmp_pair = grouped_short_edge_pair[findid][j];
						std::swap(tmp_pair.first, tmp_pair.second);
						reverse_feature.push_back(tmp_pair);
					}
					grouped_short_edges_new.push_back(reverse_feature);
				}

			}
			grouped_short_edge_pair = grouped_short_edges_new;
		}

		//second part grouped part

		input_fea_file.close();
	}

	std::cout << "Group Size: " << grouped_short_edge_pair.size() << std::endl;
	//output short_edges
	//write_hex_feature_vtk((outputname.substr(0, outputname.length() - 4) + "_fea_polycube.vtk").c_str(), hf_polycube, all_short_edge_id);
	//write_hex_feature_vtk((outputname.substr(0, outputname.length() - 4) + "_fea_ori.vtk").c_str(), hf_ori, all_short_edge_id);
	//write_hex_group_feature_feaformat((outputname.substr(0, outputname.length() - 4) + ".fea").c_str(), grouped_short_edge_pair);
	write_hex_group_feature_hfeformat(outputname.c_str(), hf_ori, grouped_short_edge_pair);
	//system("pause");
	return 1;
}