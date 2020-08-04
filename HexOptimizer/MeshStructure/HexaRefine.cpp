#include "HexaRefine.h"
#include <set>
#include <omp.h>
#include <chrono>
#include "GraphColoring.h"

namespace MeshLib
{
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	HexaRefine<Real>::HexaRefine(VHexMesh<Real>* hexmesh, Mesh3D<Real>* bndmesh, int max_iter,
		bool fixboundary, bool surf_preserve, double threshold, bool use_ejac, bool opt_illegal_only)
		:HexaSmoothBase<Real>(hexmesh, bndmesh, max_iter, fixboundary, surf_preserve, threshold, use_ejac, opt_illegal_only)
	{
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	HexaRefine<Real>::~HexaRefine()
	{
		clear_feature();
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	bool HexaRefine<Real>::has_feature()
	{
		return !trimesh_featurelines.empty() && !hexmesh_featurelines.empty();
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaRefine<Real>::clear_feature()
	{
		for (size_t i = 0; i < trimesh_featurelines.size(); i++)
			delete trimesh_featurelines[i];
		for (size_t i = 0; i < hexmesh_featurelines.size(); i++)
			delete hexmesh_featurelines[i];
		trimesh_featurelines.resize(0);
		hexmesh_featurelines.resize(0);
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	bool HexaRefine<Real>::load_features(const char trifeaturename[], const char hexfeaturename[])
	{
		clear_feature();
		try
		{
			std::string FileName(trifeaturename);
			std::string exttri = FileName.substr(FileName.find_last_of(".") + 1);

			std::string FileName2(hexfeaturename);
			std::string exthex = FileName2.substr(FileName2.find_last_of(".") + 1);

			std::vector<TinyVector<Real, 3>> meshvertices;
			std::vector<ptrdiff_t> face_id_list;
			std::vector<std::vector<size_t>> trilines, hexlines, trinewlines, hexnewlines;

			if (exttri == "tfe")
			{
				meshvertices.resize(this->m_bndmesh->get_num_of_vertices());
				face_id_list.resize(this->m_bndmesh->get_num_of_vertices());
				for (ptrdiff_t i = 0; i < this->m_bndmesh->get_num_of_vertices(); i++)
				{
					meshvertices[i] = this->m_bndmesh->get_vertex(i)->pos;
					face_id_list[i] = i;
				}
				ig::KD_Tree<Real, 3> mybndkdtree(meshvertices, face_id_list);
				if (!this->load_feature_edge(trifeaturename, mybndkdtree, trilines)) throw 1;
			}
			else if (exttri == "polyline")
			{
				if (!this->load_feature_edge(trifeaturename, trilines)) throw 1;
			}
			else
				return false;

			if (exthex == "hfe")
			{
				meshvertices.resize(this->m_quadmesh->get_num_of_vertices());
				face_id_list.resize(this->m_quadmesh->get_num_of_vertices());

				for (ptrdiff_t i = 0; i < this->m_quadmesh->get_num_of_vertices(); i++)
				{
					meshvertices[i] = this->m_quadmesh->get_vertex(i)->pos;
					face_id_list[i] = i;
				}
				ig::KD_Tree<Real, 3> myquadkdtree(meshvertices, face_id_list);
				if (!this->load_feature_edge(hexfeaturename, myquadkdtree, hexlines)) throw 1;
			}
			else if (exthex == "polyline")
			{
				if (!this->load_feature_edge(hexfeaturename, hexlines)) throw 1;
			}
			else
				return false;

			if (trilines.size() != hexlines.size())
			{
				std::cout << "The input feauture sizes do not match" << std::endl;
				throw 1;
			}

			this->reorganize_all_features(trilines, trinewlines, hexlines, hexnewlines);
			prepare_feature_edges_for_trimesh(trinewlines);
			prepare_feature_edges_for_hexmesh(hexnewlines);

			this->feature_sensitive = true;
			return true;

			//if (this->load_feature_edge(trifeaturename, mybndkdtree, trilines)
			//	&& this->load_feature_edge(hexfeaturename, myquadkdtree, hexlines))
			//{
			//	if (trilines.size() != hexlines.size())
			//	{
			//		std::cout << "The input feauture sizes do not match" << std::endl;
			//		throw 1;
			//	}

			//	this->reorganize_all_features(trilines, trinewlines, hexlines, hexnewlines);
			//	prepare_feature_edges_for_trimesh(trinewlines);
			//	prepare_feature_edges_for_hexmesh(hexnewlines);

			//	this->feature_sensitive = true;
			//	return true;
			//}
			//else
			//	throw 2;
		}
		catch (...)
		{
			this->feature_sensitive = false;
			clear_feature();
			return false;
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaRefine<Real>::reorganize_linesegment(const std::vector < std::vector<size_t>>& shortlines, std::vector < std::vector<size_t>>& lines)
	{
		lines.resize(0);
		std::map<size_t, std::vector<size_t>> vertex2line;
		for (size_t i = 0; i < shortlines.size(); i++)
		{
			vertex2line[shortlines[i][0]].push_back(i);
			vertex2line[shortlines[i][1]].push_back(i);
		}
		std::set<size_t> usedline;
		std::vector<size_t> singleline;

		for (auto iter = vertex2line.begin(); iter != vertex2line.end(); iter++)
		{
			if (iter->second.size() == 1 || iter->second.size() > 2)
			{
				bool global_find = false;
				do
				{
					global_find = false;
					singleline.resize(0);
					singleline.push_back(iter->first);
					bool find = false;
					do
					{
						find = false;
						const std::vector<size_t>& myline = vertex2line[singleline.back()];
						for (size_t j = 0; j < myline.size(); j++)
						{
							if (usedline.count(myline[j]) == 0)
							{
								size_t candidate = shortlines[myline[j]][0] == singleline.back() ? shortlines[myline[j]][1] : shortlines[myline[j]][0];

								find = true;
								usedline.insert(myline[j]);
								singleline.push_back(candidate);
								break;
							}
						}
						if (find && vertex2line[singleline.back()].size() != 2)
							break;
					} while (find);
					if (singleline.size() > 1)
					{
						lines.push_back(singleline);
						global_find = true;
					}
				} while (global_find);
			}
		}

		for (auto iter = vertex2line.begin(); iter != vertex2line.end(); iter++)
		{
			if (iter->second.size() == 2)
			{
				singleline.resize(0);
				singleline.push_back(iter->first);
				bool find = false;
				do
				{
					find = false;
					const std::vector<size_t>& myline = vertex2line[singleline.back()];
					for (size_t j = 0; j < myline.size(); j++)
					{
						if (usedline.count(myline[j]) == 0)
						{
							size_t candidate = shortlines[myline[j]][0] == singleline.back() ? shortlines[myline[j]][1] : shortlines[myline[j]][0];
							if (vertex2line[candidate].size() == 2)
							{
								find = true;
								usedline.insert(myline[j]);
								singleline.push_back(candidate);
								break;
							}
						}
					}
					if (find && singleline.front() == singleline.back())
						break;
				} while (find);
				if (singleline.size() > 1)
				{
					lines.push_back(singleline);
				}
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaRefine<Real>::prepare_feature_edges_for_trimesh(const std::vector < std::vector<size_t>>& newlines)
	{
		std::vector < TinyVector<Real, 3>> linevertices;
		for (size_t i = 0; i < newlines.size(); i++)
		{
			linevertices.resize(newlines[i].size());
			for (size_t j = 0; j < newlines[i].size(); j++)
				linevertices[j] = this->m_bndmesh->get_vertex(newlines[i][j])->pos;
			FeatureLine<Real>* newfeature = new FeatureLine<Real>;
			newfeature->build_segment_tree(linevertices);
			newfeature->get_vert_indices().assign(newlines[i].begin(), newlines[i].end());
			trimesh_featurelines.push_back(newfeature);
		}
		std::cout << "# trimesh feature edges: " << trimesh_featurelines.size() << std::endl;

		this->trifeature_edges.resize(0);
		for (size_t i = 0; i < trimesh_featurelines.size(); i++)
		{
			const std::vector<size_t>& vertindice = trimesh_featurelines[i]->get_vert_indices();
			for (size_t j = 0; j < vertindice.size() - 1; j++)
			{
				this->trifeature_edges.push_back(vertindice[j]);
				this->trifeature_edges.push_back(vertindice[j + 1]);
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaRefine<Real>::prepare_feature_edges_for_hexmesh(const std::vector < std::vector<size_t>>& newlines)
	{
		for (size_t i = 0; i < newlines.size(); i++)
		{
			FeatureLine<Real>* newfeature = new FeatureLine<Real>;
			newfeature->get_vert_indices().assign(newlines[i].begin(), newlines[i].end());
			if (newlines[i].front() == newlines[i].back())
				newfeature->get_type() = LOOP;
			else
				newfeature->get_type() = SEGMENT;
			hexmesh_featurelines.push_back(newfeature);
		}
		std::cout << "# hexmesh feature edges: " << hexmesh_featurelines.size() << std::endl;

		this->feature_edges.resize(0);
		feature_edge_targetedgeID.resize(0);
		feature_point_targetID.assign(this->m_quadmesh->get_num_of_vertices(), -1);

		for (size_t i = 0; i < hexmesh_featurelines.size(); i++)
		{
			const std::vector<size_t>& edgevert = hexmesh_featurelines[i]->get_vert_indices();
			if (edgevert.front() != edgevert.back())
			{
				this->vertex_feature_tags[edgevert.front()] = 0;
				this->vertex_feature_tags[edgevert.back()] = 0;

				for (size_t j = 1; j < edgevert.size() - 1; j++)
				{
					if (this->vertex_feature_tags[edgevert[j]] == -1)
						this->vertex_feature_tags[edgevert[j]] = 1;
				}
				for (size_t j = 0; j < edgevert.size() - 1; j++)
				{
					feature_edge_targetedgeID.push_back(i);
					this->feature_edges.push_back(edgevert[j]);
					this->feature_edges.push_back(edgevert[j + 1]);
					feature_point_targetID[edgevert[j]] = i;
				}
				feature_point_targetID[edgevert.back()] = i;
			}
			else
			{
				for (size_t j = 0; j < edgevert.size() - 1; j++)
				{
					if (this->vertex_feature_tags[edgevert[j]] == -1)
						this->vertex_feature_tags[edgevert[j]] = 1;

					feature_edge_targetedgeID.push_back(i);
					this->feature_edges.push_back(edgevert[j]);
					this->feature_edges.push_back(edgevert[j + 1]);
					feature_point_targetID[edgevert[j]] = i;
				}
				feature_point_targetID[edgevert.back()] = i;
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaRefine<Real>::dump(const char meshfile[], const char trifeaturename[], const char quadname[], const char hexfeaturename[])
	{
		if (this->m_bndmesh)
			this->m_bndmesh->write_obj(meshfile);
		if (this->m_quadmesh)
			this->m_quadmesh->write_obj(quadname);
		std::ofstream trifeatureoutput(trifeaturename), hexfeatureoutput(hexfeaturename);
		if (!trimesh_featurelines.empty() && trifeatureoutput.is_open())
		{
			size_t len = 0;
			for (size_t i = 0; i < trimesh_featurelines.size(); i++)
			{
				len += trimesh_featurelines[i]->get_vert_indices().size() - 1;
			}
			trifeatureoutput << len << std::endl;
			for (size_t i = 0; i < trimesh_featurelines.size(); i++)
			{
				const std::vector<size_t>& vertindice = trimesh_featurelines[i]->get_vert_indices();
				for (size_t j = 0; j < vertindice.size() - 1; j++)
					trifeatureoutput << vertindice[j] << ' ' << vertindice[j + 1] << std::endl;
			}
			trifeatureoutput.close();
		}
		if (!hexmesh_featurelines.empty() && hexfeatureoutput.is_open())
		{
			size_t len = 0;
			for (size_t i = 0; i < hexmesh_featurelines.size(); i++)
			{
				len += hexmesh_featurelines[i]->get_vert_indices().size() - 1;
			}
			hexfeatureoutput << len << std::endl;
			for (size_t i = 0; i < hexmesh_featurelines.size(); i++)
			{
				const std::vector<size_t>& vertindice = hexmesh_featurelines[i]->get_vert_indices();
				for (size_t j = 0; j < vertindice.size() - 1; j++)
					hexfeatureoutput << vertindice[j] << ' ' << vertindice[j + 1] << std::endl;
			}
			hexfeatureoutput.close();
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaRefine<Real>::feature_projection(size_t surf_vertex_id, const TinyVector<Real, 3>& cur_point, TinyVector<Real, 3>& projection_point)
	{
		trimesh_featurelines[feature_point_targetID[surf_vertex_id]]->project(cur_point, projection_point);
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaRefine<Real>::corner_projection(size_t surf_vertex_id, const TinyVector<Real, 3>& cur_point, TinyVector<Real, 3>& projection_point)
	{
		FeatureLine<Real>* hline = hexmesh_featurelines[feature_point_targetID[surf_vertex_id]];
		FeatureLine<Real>* tline = trimesh_featurelines[feature_point_targetID[surf_vertex_id]];
		projection_point = hline->get_vert_indices().front() == surf_vertex_id ? this->m_bndmesh->get_vertex(tline->get_vert_indices().front())->pos :
			this->m_bndmesh->get_vertex(tline->get_vert_indices().back())->pos;
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaRefine<Real>::save_trifeature_polyline(const char trifeaturename[])
	{
		std::ofstream mout(trifeaturename);
		if (!mout.is_open()) return;

		if (!trimesh_featurelines.empty())
		{
			mout << trimesh_featurelines.size() << std::endl;
			for (size_t i = 0; i < trimesh_featurelines.size(); i++)
			{
				const std::vector<size_t>& vertindice = trimesh_featurelines[i]->get_vert_indices();
				mout << vertindice.size() << std::endl;
				for (size_t j = 0; j < vertindice.size() - 1; j++)
					mout << vertindice[j] << ' ';
				mout << vertindice.back() << std::endl;
			}
		}

		mout.close();
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaRefine<Real>::save_hexfeature_polyline(const char hexfeaturename[])
	{
		std::ofstream mout(hexfeaturename);
		if (!mout.is_open()) return;

		if (!hexmesh_featurelines.empty())
		{
			mout << hexmesh_featurelines.size() << std::endl;
			for (size_t i = 0; i < hexmesh_featurelines.size(); i++)
			{
				const std::vector<size_t>& vertindice = hexmesh_featurelines[i]->get_vert_indices();
				mout << vertindice.size() << std::endl;
				for (size_t j = 0; j < vertindice.size() - 1; j++)
					mout << vertindice[j] << ' ';
				mout << vertindice.back() << std::endl;
			}
		}
		mout.close();
	}
	//////////////////////////////////////////////////////////////////////////
	template class HexaRefine<double>;
	template class FeatureLine<double>;
}