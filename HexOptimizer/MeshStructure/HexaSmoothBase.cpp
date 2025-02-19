#include "HexaSmoothBase.h"
#include <set>
#include <omp.h>
#include <chrono>
#include "GraphColoring.h"

namespace MeshLib
{
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	HexaSmoothBase<Real>::HexaSmoothBase(VHexMesh<Real>* hexmesh, Mesh3D<Real>* bndmesh, int max_iter,
		bool fixboundary, bool surf_preserve, double threshold, bool use_ejac, bool opt_illegal_only)
		:m_hexmesh(hexmesh), m_bndmesh(bndmesh), max_iters(max_iter), fix_boundary(fixboundary), feature_sensitive(surf_preserve), bound(threshold),
		use_extend_jac(use_ejac), opt_illegal_vertex_only(opt_illegal_only), m_quadmesh(0), m_aabb(0)
	{
		if (m_hexmesh == 0 || m_hexmesh->get_elements().empty()) return;

		std::cout << "Hex mesh: " << m_hexmesh->get_vertices().size() << " vertices; " << m_hexmesh->get_elements().size() << " hexes." << std::endl;
		if (m_bndmesh)
			std::cout << "boundary mesh: " << m_bndmesh->get_num_of_vertices() << " vertices; " << m_bndmesh->get_num_of_faces() << " faces." << std::endl;

		////////////////////////
		initialize();
		////////////////////////
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	HexaSmoothBase<Real>::~HexaSmoothBase()
	{
		if (m_aabb)
			delete m_aabb;
		if (m_quadmesh)
			delete m_quadmesh;
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	bool HexaSmoothBase<Real>::load_features(const char trifeaturename[], const char hexfeaturename[])
	{
		std::cout << "please override this function" << std::endl;
		return false;
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::save_quadmesh(const char quadname[])
	{
		if (m_quadmesh)
		{
			const std::vector< TinyVector<Real, 3> >& vertices = m_hexmesh->get_vertices();

			for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
				m_quadmesh->get_vertex(i)->pos = vertices[reverse_ids[i]];
			m_quadmesh->write_obj(quadname);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::save_trimesh(const char triname[])
	{
		if (m_bndmesh)
		{
			m_bndmesh->write_obj(triname);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::save_trifeature(const char feaname[])
	{
		std::ofstream mfile(feaname);
		if (!mfile.is_open()) return;
		mfile << trifeature_edges.size() / 2 << std::endl;
		for (size_t i = 0; i < trifeature_edges.size(); i += 2)
			mfile << trifeature_edges[i] << ' ' << trifeature_edges[i + 1] << std::endl;
		mfile.close();
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::save_quadfeature(const char feaname[])
	{
		std::ofstream mfile(feaname);
		if (!mfile.is_open()) return;
		mfile << feature_edges.size() / 2 << std::endl;
		for (size_t i = 0; i < feature_edges.size(); i += 2)
			mfile << feature_edges[i] << ' ' << feature_edges[i + 1] << std::endl;
		mfile.close();
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	bool HexaSmoothBase<Real>::load_feature_edge(const char filename[], ig::KD_Tree<Real, 3>& kdtree, std::vector < std::vector<size_t>>& alllines)
	{
		std::string FileName(filename);
		std::string ext = FileName.substr(FileName.find_last_of(".") + 1);
		if (ext != "tfe" && ext != "hfe") return false;
		std::vector < std::vector<size_t>> newlines;
		alllines.resize(0);
		std::ifstream mfile(filename);
		if (mfile.is_open())
		{
			std::vector<size_t> singleline;
			std::string head;
			mfile >> head >> head >> head;
			int num_features = 0, npts = 0;
			mfile >> num_features;
			TinyVector<Real, 3> P, nP;
			ptrdiff_t id;

			for (int i = 0; i < num_features; i++)
			{
				mfile >> npts;
				singleline.resize(0);
				for (int j = 0; j < npts; j++)
				{
					mfile >> P;
					kdtree.find_Nearest_Point(P, nP, &id);
					singleline.push_back(id);
				}
				alllines.push_back(singleline);
			}
			return true;
		}
		else
			return false;
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	bool HexaSmoothBase<Real>::load_feature_edge(const char filename[], std::vector < std::vector<size_t>>& alllines)
	{
		std::string FileName(filename);
		std::string ext = FileName.substr(FileName.find_last_of(".") + 1);
		if (ext != "polyline") return false;
		std::vector < std::vector<size_t>> newlines;
		alllines.resize(0);
		std::ifstream mfile(filename);
		if (mfile.is_open())
		{
			std::vector<size_t> singleline;
			int num_features = 0, npts = 0, id = 0;
			mfile >> num_features;

			for (int i = 0; i < num_features; i++)
			{
				mfile >> npts;
				singleline.resize(0);
				for (int j = 0; j < npts; j++)
				{
					mfile >> id;
					singleline.push_back(id);
				}
				alllines.push_back(singleline);
			}
			return true;
		}
		else
			return false;
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::smooth(bool laplaciansmooth, bool verbose)
	{
		if (verbose)
		{
			print_scalejacobian();
			auto start = std::chrono::high_resolution_clock::now();
			std::cout << "perform laplacian smoothing first..." << std::endl;
			laplacian_smooth_init();
			//if (laplaciansmooth)
			fast_laplacian_smooth(20, !laplaciansmooth);
			//laplacian_smooth(10, false, !laplaciansmooth);
			optimize(max_iters, fix_boundary);
			auto finish = std::chrono::high_resolution_clock::now();
			print_scalejacobian();
			std::chrono::duration<double> elapsed = finish - start;
			std::cout << "Timing: " << elapsed.count() << " s\n";
		}
		else
		{
			laplacian_smooth_init();
			laplacian_smooth(10, false, !laplaciansmooth);
			optimize(max_iters, fix_boundary, verbose);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::feature_projection(size_t surf_vertex_id, const TinyVector<Real, 3>& cur_point, TinyVector<Real, 3>& projection_point)
	{
		std::cout << "This function should be overrided" << std::endl;
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::corner_projection(size_t surf_vertex_id, const TinyVector<Real, 3>& cur_point, TinyVector<Real, 3>& projection_point)
	{
		std::cout << "This function should be overrided" << std::endl;
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::initialize()
	{
		std::vector< TinyVector<Real, 3> >& vertices = m_hexmesh->get_vertices();
		std::vector< HexElement<Real>* >& elements = m_hexmesh->get_elements();

		//create aabb tree
		std::vector<TinyVector<Real, 3>> bnd_face_vertices_list;
		if (m_bndmesh)
		{
			for (ptrdiff_t i = 0; i < m_bndmesh->get_num_of_faces(); i++)
			{
				HE_face<Real>* hf = m_bndmesh->get_face(i);
				HE_edge<Real>* he = hf->edge; he = he->next;
				for (unsigned int k = 0; k < hf->valence - 2; k++)
				{
					bnd_face_vertices_list.push_back(hf->edge->pair->vert->pos);
					bnd_face_vertices_list.push_back(he->pair->vert->pos);
					bnd_face_vertices_list.push_back(he->vert->pos);
				}
			}
		}
		else
		{
			for (size_t i = 0; i < elements.size(); i++)
			{
				for (unsigned int j = 0; j < 6; j++)
				{
					if (elements[i]->neighborHex[j].first == 0)
					{
						bnd_face_vertices_list.push_back(vertices[elements[i]->indices[face_id[j][0]]]);
						bnd_face_vertices_list.push_back(vertices[elements[i]->indices[face_id[j][1]]]);
						bnd_face_vertices_list.push_back(vertices[elements[i]->indices[face_id[j][2]]]);
						bnd_face_vertices_list.push_back(vertices[elements[i]->indices[face_id[j][0]]]);
						bnd_face_vertices_list.push_back(vertices[elements[i]->indices[face_id[j][2]]]);
						bnd_face_vertices_list.push_back(vertices[elements[i]->indices[face_id[j][3]]]);
					}
				}
			}
		}
		if (m_aabb) delete m_aabb;
		m_aabb = new ig::AABB_Tree<Real>(bnd_face_vertices_list);

		//boundary tag
		boundaryvertexlabel.assign(vertices.size(), false);
		ave_edgelength = 0;
		size_t counter = 0;
		for (size_t i = 0; i < elements.size(); i++)
		{
			for (int j = 0; j < 6; j++)
			{
				if (elements[i]->neighborHex[j].first != 0) { continue; }
				boundaryfacets.push_back(CellFace(i, j));
				for (int k = 0; k < 4; k++)
				{
					const size_t& vertID = elements[i]->indices[face_id[j][k]];
					boundaryvertexlabel[vertID] = true;
					ave_edgelength += (vertices[elements[i]->indices[face_id[j][k]]] - vertices[elements[i]->indices[face_id[j][(k + 1) % 4]]]).SquaredLength();
					counter++;
				}
			}
		}

		vertex2cornerlist.resize(vertices.size());
		//initial corner list

		int numbound = use_extend_jac ? 32 : 8;
		for (size_t i = 0; i < elements.size(); i++)
		{
			for (int j = 0; j < numbound; j++)
			{
				const size_t& id = elements[i]->indices[corner_vert[j][0]];
				const size_t& v1 = elements[i]->indices[corner_vert[j][1]];
				const size_t& v2 = elements[i]->indices[corner_vert[j][2]];
				const size_t& v3 = elements[i]->indices[corner_vert[j][3]];

				CellCorner hc(i, j, j);
				vertex2cornerlist[id].push_back(hc);
				vertex2cornerlist[v1].push_back(hc);
				vertex2cornerlist[v2].push_back(hc);
				vertex2cornerlist[v3].push_back(hc);
			}
		}

		ave_edgelength = std::sqrt(ave_edgelength / counter);

		m_quadmesh = m_hexmesh->extract_boundary_quad_mesh(vert_map_ids);
		reverse_ids.resize(vert_map_ids.size());
		for (typename std::map<size_t, size_t>::iterator iter = vert_map_ids.begin();
			iter != vert_map_ids.end(); iter++)
		{
			reverse_ids[iter->second] = iter->first;
		}

		vertex_feature_tags.assign(m_quadmesh->get_num_of_vertices(), -1);
		graph_coloring();
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	Real HexaSmoothBase<Real>::compute_distortion_and_gradient(const size_t id, const TinyVector<Real, 3>& vertex, const OptData<Real>& data, bool& has_flip, TinyVector<Real, 3>* gradient/* = 0*/)
	{
		Real distortion = 0;
		std::vector< TinyVector<Real, 3> >& vertices = m_hexmesh->get_vertices();
		std::vector< HexElement<Real>* >& elements = m_hexmesh->get_elements();

		TinyVector<Real, 3> e[3];
		Real selen[3];
		Real max_e = 0;
		has_flip = false;
		TinyVector<Real, 3> L, R, tmp_gradient;
		if (gradient) gradient->reset_zero();

		for (size_t ic = 0; ic < vertex2cornerlist[id].size(); ic++)
		{
			const CellCorner& hc = vertex2cornerlist[id][ic];

			Real K = 1;
			if (hc.cornertag >= 8 && hc.cornertag < 20)
			{
				K = K1;
			}
			else if (hc.cornertag >= 20 && hc.cornertag < 32)
			{
				K = K2;
			}

			const size_t& v0 = elements[hc.cell_id]->indices[corner_vert[hc.corner][0]];
			const size_t& v1 = elements[hc.cell_id]->indices[corner_vert[hc.corner][1]];
			const size_t& v2 = elements[hc.cell_id]->indices[corner_vert[hc.corner][2]];
			const size_t& v3 = elements[hc.cell_id]->indices[corner_vert[hc.corner][3]];

			const TinyVector<Real, 3>& V0 = v0 == id ? vertex : vertices[v0];
			const TinyVector<Real, 3>& V1 = v1 == id ? vertex : vertices[v1];
			const TinyVector<Real, 3>& V2 = v2 == id ? vertex : vertices[v2];
			const TinyVector<Real, 3>& V3 = v3 == id ? vertex : vertices[v3];

			int index = -1;
			if (v0 == id)
				index = 0;
			else if (v1 == id)
				index = 1;
			else if (v2 == id)
				index = 2;
			else
				index = 3;

			bool is_boundary_edges = boundaryvertexlabel[v0] && boundaryvertexlabel[v1] &&
				boundaryvertexlabel[v2] && boundaryvertexlabel[v3];

			Real mybound = is_boundary_edges ? std::min(bound, 1.0e-2) : bound;

			e[0] = V1 - V0; e[1] = V2 - V0; e[2] = V3 - V0;
			selen[0] = e[0].SquaredLength(); selen[1] = e[1].SquaredLength(); selen[2] = e[2].SquaredLength();
			Real jacobian = e[2].Dot(e[0].Cross(e[1]));
			Real stelen = selen[0] * selen[1] * selen[2];
			Real telen = sqrt(stelen);

			if (telen < 1.0e-10)
				return DBL_MAX;

			Real scale_jacobian = telen == 0 ? DBL_MAX : jacobian / telen;

			Real tmp_corner_jacobian = scale_jacobian;
			if (scale_jacobian > K)
			{
				tmp_corner_jacobian = 1 + K - tmp_corner_jacobian;
			}
			else if (scale_jacobian < -K)
			{
				tmp_corner_jacobian = -(1 + K) - tmp_corner_jacobian;
			}
			else
			{
				tmp_corner_jacobian /= K;
			}

			tmp_corner_jacobian -= mybound;
			if (tmp_corner_jacobian <= 0)
			{
				has_flip = true;
			}
			else if (opt_illegal_vertex_only)
			{
				continue;
			}

			if (gradient)
			{
				switch (index)
				{
				case 0:
					L = -(e[0].Cross(e[1] - e[2]) + e[1].Cross(e[2]));
					R = -(e[0] / selen[0] + e[1] / selen[1] + e[2] / selen[2]);
					break;
				case 1:
					L = e[1].Cross(e[2]);
					R = e[0] / selen[0];
					break;
				case 2:
					L = e[2].Cross(e[0]);
					R = e[1] / selen[1];
					break;
				case 3:
					L = e[0].Cross(e[1]);
					R = e[2] / selen[2];
					break;
				default:
					break;
				}

				Real sdenom1 = tmp_corner_jacobian * tmp_corner_jacobian + threshold;
				Real denom1 = sqrt(sdenom1);
				Real denom = tmp_corner_jacobian + denom1;

				if (scale_jacobian > K)
				{
					tmp_gradient = ((1 + tmp_corner_jacobian / denom1) / (denom * denom * telen)) * (L - jacobian * R);
				}
				else if (scale_jacobian < -K)
				{
					tmp_gradient = ((1 + tmp_corner_jacobian / denom1) / (denom * denom * telen)) * (L - jacobian * R);
				}
				else
				{
					tmp_gradient = (-(1 + tmp_corner_jacobian / denom1) / (K * denom * denom * telen)) * (L - jacobian * R);
				}
				scale_jacobian = denom;
			}
			else
			{
				scale_jacobian = tmp_corner_jacobian + sqrt(tmp_corner_jacobian * tmp_corner_jacobian + threshold);
			}

			distortion += 1.0 / scale_jacobian;

			if (gradient)
				*gradient += tmp_gradient;

			max_e = std::max(max_e, std::max(std::max(selen[0], selen[1]), selen[2]));
		}

		if (boundaryvertexlabel[id])
		{
			TinyVector<Real, 3> Z = vertex - data.proj;
			distortion += data.weight * Z.SquaredLength();
			if (gradient)
				*gradient += 2.0 * data.weight * Z;
		}

		if (gradient)
		{
			gradient->Normalize(); *gradient *= sqrt(max_e);
		}

		return distortion;
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::optimize(int max_iter, bool fixboundary, bool verbose)
	{
		std::vector< TinyVector<Real, 3> >& vertices = m_hexmesh->get_vertices();
		std::vector< HexElement<Real>* >& elements = m_hexmesh->get_elements();

		if (verbose)
			std::cout << "optimizing";

		int global_iter = 0, local_iter = 10;

		std::vector<bool> vertex_flips(vertices.size(), false);

		std::vector< TinyVector<Real, 3>> backup(vertices);

		//std::vector<bool> vertexirregularlabel;
		//extract_irrugular_lines(vertexirregularlabel, false);

		int resolve = 0;

		do
		{
			int changed = 0;
			int has_flip = 0;

			for (size_t color = 0; color < colored_vertices.size(); color++)
			{
#pragma omp parallel for reduction(+:changed,has_flip) schedule(dynamic)
				for (ptrdiff_t vi = 0; vi < (ptrdiff_t)colored_vertices[color].size(); vi++)
				{
					size_t id = colored_vertices[color][vi];

					if (fixboundary && boundaryvertexlabel[id]) continue;

					OptData<Real> opt;
					int tag = -1; ptrdiff_t surf_id = -1;
					if (boundaryvertexlabel[id])
					{
						surf_id = vert_map_ids[id];
						tag = vertex_feature_tags[surf_id];
						//if (tag == 0) continue;
						opt.edgelen = 1;
						opt.proj = backup[id];
						opt.weight = 1.0e5;
						if (tag != -1)
							opt.weight = 1.0e6;
						//if (tag == -1)
						{
							//opt.weight = 1.0e3;
							/*
							if (vertexirregularlabel[surf_id] || vertex_flips[surf_id])
								opt.weight = 1.0e3;
							else
								opt.weight = 1.0e5;
								*/
						}
					}
					bool flip;
					TinyVector<Real, 3> gradient, oldgradient, tmpV, tmpV2;

					Real distortion = compute_distortion_and_gradient(id, vertices[id], opt, flip, &gradient);
					oldgradient = gradient;
					double l = 1;
					bool suc = false;

					for (int iter = 0; iter < local_iter; iter++)
					{
						tmpV = vertices[id] - l * gradient;

						if (boundaryvertexlabel[id])
						{
							if (feature_sensitive)
							{
								if (tag == 0)
								{
									corner_projection(surf_id, tmpV, tmpV);
								}
								else if (tag == 1)
								{
									// if (m_quadmesh->get_vertex(surf_id)->degree == 4)
										feature_projection(surf_id, tmpV, tmpV);
								}
								else
								{
									m_aabb->find_Nearest_Point(tmpV, tmpV);
								}
							}
							// else
							// 	m_aabb->find_Nearest_Point(tmpV, tmpV);
						}

						bool new_flip;
						Real new_distortion = compute_distortion_and_gradient(id, tmpV, opt, new_flip);

						if (new_distortion < distortion)
						{
							flip = new_flip;
							suc = true;
							break;
						}

						l *= 0.5;
					}
					vertex_flips[id] = flip;
					if (flip) has_flip++;

					if (suc)
					{
						vertices[id] = tmpV;
						changed++;
					}
				}
			}
			if (verbose)
			{
				std::cout << "-";
			}
			if (changed == 0) break;
			if (has_flip == 0 && opt_illegal_vertex_only) break;
			global_iter++;

			if (global_iter == max_iter && has_flip && resolve < 5)
			{
				if (verbose)
					std::cout << std::endl << "restart laplacian smoothing" << std::endl << "optimizing";
				fast_laplacian_smooth(5, false);
				backup.assign(vertices.begin(), vertices.end());
				global_iter = 0; resolve++;
			}
		} while (global_iter < max_iter);


// 		for (size_t color = 0; color < colored_vertices.size(); color++)
// 		{
// #pragma omp parallel for schedule(dynamic)
// 			for (ptrdiff_t vi = 0; vi < (ptrdiff_t)colored_vertices[color].size(); vi++)
// 			{
// 				size_t id = colored_vertices[color][vi];

// 				if (fixboundary && boundaryvertexlabel[id]) continue;

// 				OptData<Real> opt;
// 				int tag = -1; ptrdiff_t surf_id = -1;
// 				if (boundaryvertexlabel[id])
// 				{
// 					surf_id = vert_map_ids[id];
// 					tag = vertex_feature_tags[surf_id];

// 					if (feature_sensitive)
// 					{
// 						if (tag == 0)
// 						{
// 							// corner_projection(surf_id, tmpV, tmpV);
// 						}
// 						else if (tag == 1)
// 						{
// 							// if (m_quadmesh->get_vertex(surf_id)->degree == 4)
// 							// 	feature_projection(surf_id, tmpV, tmpV);
// 						}
// 						else
// 						{
// 							m_aabb->find_Nearest_Point(vertices[id], vertices[id]);
// 						}
// 					}
// 							// else
// 							// 	m_aabb->find_Nearest_Point(tmpV, tmpV);
// 				}
// 			}
// 		}


		if (verbose)
			std::cout << "done!" << std::endl;
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::extract_irrugular_lines(std::vector<bool>& vertexirregularlabel, bool detect_structure)
	{
		vertexirregularlabel.resize(m_quadmesh->get_num_of_vertices(), false);
		m_quadmesh->reset_edges_tag(false);

		for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
		{
			HE_vert<Real>* hv = m_quadmesh->get_vertex(i);
			if (vertex_feature_tags[i] == -1) continue;
			if (vertex_feature_tags[i] == 0 || hv->degree > 4)
			{
				HE_edge<Real>* he = hv->edge;
				do
				{
					HE_edge<Real>* next_edge = he;

					std::vector<HE_edge<Real>*> edgestore;
					bool suc = false;
					while (!next_edge->tag && featureedgelabels[next_edge->id] == false)
					{
						edgestore.push_back(next_edge);
						if (next_edge->vert->degree != 4 || vertex_feature_tags[next_edge->vert->id] != -1)
						{
							suc = true;
							break;
						}
						next_edge = next_edge->next->pair->next;
					}

					if (suc)
					{
						for (auto eiter = edgestore.begin(); eiter != edgestore.end(); eiter++)
						{
							(*eiter)->tag = (*eiter)->pair->tag = true;
							vertexirregularlabel[(*eiter)->vert->id] = vertexirregularlabel[(*eiter)->pair->vert->id] = true;
						}
					}
					he = he->pair->next;
				} while (he != hv->edge);
			}
		}
		//for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
		//{
		//	HE_vert<Real>* hv = m_quadmesh->get_vertex(i);
		//	//if (vertex_feature_tags[i] == 0)
		//	if (hv->degree != 4 && vertex_feature_tags[i] == 0)
		//		vertexirregularlabel[i] = true;
		//}

		return;

		for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
		{
			HE_vert<Real>* hv = m_quadmesh->get_vertex(i);
			if (!detect_structure)
			{
				if (hv->degree == 4 || vertex_feature_tags[i] != 1) continue;
			}
			else
			{
				if (hv->degree == 4) continue;
			}

			HE_edge<Real>* he = hv->edge;
			do
			{
				if (!detect_structure)
				{
					if (vertex_feature_tags[he->vert->id] != -1) break;
				}
				HE_edge<Real>* next_edge = he;
				while (!next_edge->tag)
				{
					next_edge->tag = next_edge->pair->tag = true;
					vertexirregularlabel[next_edge->vert->id] = vertexirregularlabel[next_edge->pair->vert->id] = true;
					if (next_edge->vert->degree != 4)
						break;
					next_edge = next_edge->next->pair->next;
				}
				he = he->pair->next;
			} while (he != hv->edge);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::print_scalejacobian()
	{
		std::vector< TinyVector<Real, 3> >& vertices = m_hexmesh->get_vertices();
		std::vector< HexElement<Real>* >& elements = m_hexmesh->get_elements();

		std::vector<Real> m_hex_jacobians(elements.size()), m_hex_jacobians_ext(elements.size());

		Real avg_jac = 0.0;
		Real min_jac = 1.0, min_jac_ext = 1.0;
		Real max_jac = -1.0;
		size_t neg_hexes = 0, neg_hexes_ext = 0;

#pragma omp parallel for reduction(+:neg_hexes,neg_hexes_ext,avg_jac)
		for (ptrdiff_t i = 0; i < (ptrdiff_t)elements.size(); i++)
		{
			TinyVector<Real, 3> e1, e2, e3;
			Real min_corner_jacobian = 2.0, min_corner_jacobian_ext = 2.0;
			for (int j = 0; j < 8; j++)
			{
				e1 = vertices[elements[i]->indices[corner_vert[j][1]]] - vertices[elements[i]->indices[corner_vert[j][0]]];
				e2 = vertices[elements[i]->indices[corner_vert[j][2]]] - vertices[elements[i]->indices[corner_vert[j][0]]];
				e3 = vertices[elements[i]->indices[corner_vert[j][3]]] - vertices[elements[i]->indices[corner_vert[j][0]]];
				e1.Normalize(); e2.Normalize(); e3.Normalize();
				Real tmp_corner_jacobian = element_scaledjacobian(j, e1, e2, e3);
				if (tmp_corner_jacobian < min_corner_jacobian)
				{
					min_corner_jacobian = tmp_corner_jacobian;
				}
			}
			min_corner_jacobian_ext = min_corner_jacobian;
			m_hex_jacobians[i] = min_corner_jacobian;
			avg_jac += m_hex_jacobians[i];
			if (min_corner_jacobian <= 0) neg_hexes++;

			for (int j = 8; j < 32; j++)
			{
				e1 = vertices[elements[i]->indices[corner_vert[j][1]]] - vertices[elements[i]->indices[corner_vert[j][0]]];
				e2 = vertices[elements[i]->indices[corner_vert[j][2]]] - vertices[elements[i]->indices[corner_vert[j][0]]];
				e3 = vertices[elements[i]->indices[corner_vert[j][3]]] - vertices[elements[i]->indices[corner_vert[j][0]]];
				e1.Normalize(); e2.Normalize(); e3.Normalize();
				Real tmp_corner_jacobian = element_scaledjacobian(j, e1, e2, e3);
				if (tmp_corner_jacobian < min_corner_jacobian_ext)
				{
					min_corner_jacobian_ext = tmp_corner_jacobian;
				}
			}
			if (min_corner_jacobian_ext <= 0) neg_hexes_ext++;
			m_hex_jacobians_ext[i] = min_corner_jacobian_ext;
		}
		max_jac = *std::max_element(m_hex_jacobians.begin(), m_hex_jacobians.end());
		min_jac = *std::min_element(m_hex_jacobians.begin(), m_hex_jacobians.end());
		min_jac_ext = *std::min_element(m_hex_jacobians_ext.begin(), m_hex_jacobians_ext.end());
		avg_jac /= elements.size();
		Real std_jac = 0.0;
#pragma omp parallel for reduction(+:std_jac)
		for (ptrdiff_t i = 0; i < (ptrdiff_t)elements.size(); i++)
		{
			Real dv = m_hex_jacobians[i] - avg_jac;
			std_jac += dv * dv;
		}
		std_jac = std::sqrt(std_jac / elements.size());

		Real ave_dis = 0;
		size_t num = 0;
		std::vector<Real> max_dis_array(omp_get_max_threads(), 0);
		int m_sample_rate = 4;
#pragma omp parallel for reduction(+:ave_dis,num)
		for (ptrdiff_t i = 0; i < (ptrdiff_t)boundaryfacets.size(); i++)
		{
			const size_t& id = boundaryfacets[i].cell_id;
			const int& face = boundaryfacets[i].face;
			const size_t* indices = elements[id]->indices;
			int thread = omp_get_thread_num();
			TinyVector<Real, 3> T, P;
			for (int xi = 0; xi <= m_sample_rate; xi++)
			{
				Real x = (double)xi / m_sample_rate;
				for (int yi = 0; yi <= m_sample_rate; yi++)
				{
					Real y = (double)yi / m_sample_rate;

					T = ((1 - x) * (1 - y)) * vertices[indices[face_id[face][0]]] + (x * (1 - y)) * vertices[indices[face_id[face][1]]] +
						((1 - x) * y) * vertices[indices[face_id[face][3]]] + (x * y) * vertices[indices[face_id[face][2]]];
					m_aabb->find_Nearest_Point(T, P);
					Real tmp_dis = (T - P).SquaredLength();
					max_dis_array[thread] = std::max(max_dis_array[thread], tmp_dis);
					ave_dis += tmp_dis;
					num++;
				}
			}
		}
		ave_dis = sqrt(ave_dis / num);
		Real max_dis = sqrt(*std::max_element(max_dis_array.begin(), max_dis_array.end()));

		std::cout << "------ quality assessment (scaled jacobian) -------" << std::endl;
		std::cout << "# avg jacobian: " << avg_jac << std::endl;
		std::cout << "# min jacobian: " << min_jac << std::endl;
		std::cout << "# max jacobian: " << max_jac << std::endl;
		std::cout << "# stdev jacobian: " << std_jac << std::endl;
		std::cout << "# stdev/avg jacobian: " << std_jac / avg_jac << std::endl;
		std::cout << "# flipped hexes: " << neg_hexes << std::endl;
		std::cout << "# min jacobian ext: " << min_jac_ext << std::endl;
		std::cout << "# flipped extend hexes: " << neg_hexes_ext << std::endl;
		std::cout << "# avg bnd distance: " << ave_dis << std::endl;
		std::cout << "# max bnd distance: " << max_dis << std::endl;
		std::cout << "---------------------------------------------------" << std::endl;
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::graph_coloring()
	{
		std::vector< TinyVector<Real, 3> >& vertices = m_hexmesh->get_vertices();
		std::vector< HexElement<Real>* >& elements = m_hexmesh->get_elements();

		std::vector<std::set<size_t>> edges(vertices.size());
		int numbound = use_extend_jac ? 32 : 8;
		for (size_t i = 0; i < elements.size(); i++)
		{
			for (int j = 0; j < numbound; j++)
			{
				const size_t& id = elements[i]->indices[corner_vert[j][0]];
				const size_t& v1 = elements[i]->indices[corner_vert[j][1]];
				const size_t& v2 = elements[i]->indices[corner_vert[j][2]];
				const size_t& v3 = elements[i]->indices[corner_vert[j][3]];
				edges[id].insert(v1); edges[id].insert(v2); edges[id].insert(v3);
				edges[v1].insert(id); edges[v1].insert(v2); edges[v1].insert(v3);
				edges[v2].insert(id); edges[v2].insert(v1); edges[v2].insert(v3);
				edges[v3].insert(id); edges[v3].insert(v1); edges[v3].insert(v2);
			}
		}

		greedy_graph_coloring(vertices.size(), edges, colored_vertices);

		std::vector<std::set<size_t>> facet_edges(m_quadmesh->get_num_of_faces());
		std::vector<size_t> neighbor;
		for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
		{
			HE_vert<Real>* hv = m_quadmesh->get_vertex(i);
			HE_edge<Real>* he = hv->edge;
			neighbor.resize(0);
			do
			{
				if (he->face)
					neighbor.push_back((size_t)he->face->id);
				he = he->pair->next;
			} while (he != hv->edge);
			for (size_t j = 0; j < neighbor.size(); j++)
			{
				for (size_t k = 0; k < neighbor.size(); k++)
				{
					if (j != k)
						facet_edges[neighbor[j]].insert(neighbor[k]);
				}
			}
		}
		greedy_graph_coloring(m_quadmesh->get_num_of_faces(), facet_edges, colored_facets);
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	Real HexaSmoothBase<Real>::element_scaledjacobian(const int corner_type, const TinyVector<Real, 3>& V0, const TinyVector<Real, 3>& V1, const TinyVector<Real, 3>& V2)
	{
		return V2.Dot(V0.Cross(V1));

		/*
		Real K = 1;
		if (corner_type >= 8 && corner_type < 20)
		{
			K = K1;
		}
		else if (corner_type >= 20 && corner_type < 32)
		{
			K = K2;
		}
		else
			K = 1;

		Real tmp_corner_jacobian = V2.Dot(V0.Cross(V1));
		if (tmp_corner_jacobian > K)
		{
			tmp_corner_jacobian = 1 + K - tmp_corner_jacobian;
		}
		else if (tmp_corner_jacobian < -K)
		{
			tmp_corner_jacobian = -(1 + K) - tmp_corner_jacobian;
		}
		else
		{
			tmp_corner_jacobian /= K;
		}
		return tmp_corner_jacobian;
		*/
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::reorganize(const std::vector < std::vector<size_t>>& lines, std::vector < std::vector<size_t>>& newlines)
	{
		newlines.resize(0);
		std::map<size_t, std::vector<size_t>> vertex2line;
		for (size_t i = 0; i < lines.size(); i++)
		{
			vertex2line[lines[i].front()].push_back(i);
			vertex2line[lines[i].back()].push_back(i);
		}
		std::vector<bool> linetag(lines.size(), false);

		for (size_t i = 0; i < lines.size(); i++)
		{
			if (linetag[i]) continue;

			std::vector<size_t> newline;
			if (vertex2line[lines[i].back()].size() != 2)
				newline.assign(lines[i].rbegin(), lines[i].rend());
			else if (vertex2line[lines[i].front()].size() != 2)
				newline.assign(lines[i].begin(), lines[i].end());
			else
				continue;

			linetag[i] = true;
			size_t prev_line = i;
			bool find = false;
			do
			{
				find = false;
				size_t linebegin = newline.back();
				const std::vector<size_t>& shortline = vertex2line[linebegin];
				size_t degree = shortline.size();
				if (degree == 2)
				{
					if (shortline.back() == prev_line)
					{
						if (!linetag[shortline.front()])
						{
							prev_line = shortline.front();
							find = true;
						}
					}
					else
					{
						if (!linetag[shortline.back()])
						{
							prev_line = shortline.back();
							find = true;
						}
					}
					if (find)
					{
						linetag[prev_line] = true;
						if (lines[prev_line].front() == linebegin)
							newline.insert(newline.end(), lines[prev_line].begin() + 1, lines[prev_line].end());
						else
							newline.insert(newline.end(), lines[prev_line].rbegin() + 1, lines[prev_line].rend());
					}
				}
			} while (find);
			newlines.push_back(newline);
		}
		for (size_t i = 0; i < lines.size(); i++)
		{
			if (linetag[i]) continue;

			std::vector<size_t> newline;
			newline.assign(lines[i].begin(), lines[i].end());

			linetag[i] = true;
			size_t prev_line = i;
			bool find = false;
			do
			{
				find = false;
				size_t linebegin = newline.back();
				const std::vector<size_t>& shortline = vertex2line[linebegin];
				size_t degree = shortline.size();
				if (degree == 2)
				{
					if (shortline.back() == prev_line)
					{
						if (!linetag[shortline.front()])
						{
							prev_line = shortline.front();
							find = true;
						}
					}
					else
					{
						if (!linetag[shortline.back()])
						{
							prev_line = shortline.back();
							find = true;
						}
					}
					if (find)
					{
						linetag[prev_line] = true;
						if (lines[prev_line].front() == linebegin)
							newline.insert(newline.end(), lines[prev_line].begin() + 1, lines[prev_line].end());
						else
							newline.insert(newline.end(), lines[prev_line].rbegin() + 1, lines[prev_line].rend());
					}
				}
			} while (find);

			newlines.push_back(newline);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::reorganize_all_features(
		const std::vector < std::vector<size_t>>& trilines,
		std::vector < std::vector<size_t>>& trinewlines,
		const std::vector < std::vector<size_t>>& hexlines,
		std::vector < std::vector<size_t>>& hexnewlines
	)
	{
		trinewlines.resize(0); hexnewlines.resize(0);

		std::map<size_t, std::vector<size_t>> trivertex2line;
		std::map<size_t, std::vector<size_t>> neighbor;
		for (size_t i = 0; i < trilines.size(); i++)
		{
			trivertex2line[trilines[i].front()].push_back(i);
			trivertex2line[trilines[i].back()].push_back(i);
			neighbor[trilines[i].front()].push_back(trilines[i][1]);
			neighbor[trilines[i].back()].push_back(trilines[i][trilines[i].size() - 2]);
		}

		for (auto iter = trivertex2line.begin(); iter != trivertex2line.end(); iter++)
		{
			if (iter->second.size() == 2)
			{
				HE_vert<Real>* hv = m_bndmesh->get_vertex(iter->first);
				HE_vert<Real>* hv0 = m_bndmesh->get_vertex(neighbor[iter->first][0]);
				HE_vert<Real>* hv1 = m_bndmesh->get_vertex(neighbor[iter->first][1]);
				TinyVector<Real, 3> V0 = hv0->pos - hv->pos, V1 = hv1->pos - hv->pos;
				V0.Normalize(); V1.Normalize();
				//if (std::fabs(V0.Dot(V1)) < 0.9) //cos(25 degree)
				if (std::fabs(V0.Dot(V1)) < 0.5) //cos(60 degree)
				{
					iter->second.push_back(iter->first);// change its size != 2
					//std::cout << "angle " << acos(V0.Dot(V1)) * 180 / 3.14 << std::endl;
					//std::cout << hv0->id << " " << hv1->id << std::endl;
				}
			}
		}

		std::vector<bool> trilinetag(trilines.size(), false);

		for (size_t i = 0; i < trilines.size(); i++)
		{
			if (trilinetag[i]) continue;

			std::vector<size_t> newline, hexline;
			if (trivertex2line[trilines[i].back()].size() != 2)
			{
				newline.assign(trilines[i].rbegin(), trilines[i].rend());
				hexline.assign(hexlines[i].rbegin(), hexlines[i].rend());
			}
			else if (trivertex2line[trilines[i].front()].size() != 2)
			{
				newline.assign(trilines[i].begin(), trilines[i].end());
				hexline.assign(hexlines[i].begin(), hexlines[i].end());
			}
			else
				continue;

			trilinetag[i] = true;
			size_t prev_line = i;
			bool find = false;
			do
			{
				find = false;
				size_t linebegin = newline.back();
				const std::vector<size_t>& shortline = trivertex2line[linebegin];
				size_t degree = shortline.size();
				if (degree == 2)
				{
					if (shortline.back() == prev_line)
					{
						if (!trilinetag[shortline.front()])
						{
							prev_line = shortline.front();
							find = true;
						}
					}
					else
					{
						if (!trilinetag[shortline.back()])
						{
							prev_line = shortline.back();
							find = true;
						}
					}
					if (find)
					{
						trilinetag[prev_line] = true;
						if (trilines[prev_line].front() == linebegin)
						{
							newline.insert(newline.end(), trilines[prev_line].begin() + 1, trilines[prev_line].end());
							hexline.insert(hexline.end(), hexlines[prev_line].begin() + 1, hexlines[prev_line].end());
						}
						else
						{
							newline.insert(newline.end(), trilines[prev_line].rbegin() + 1, trilines[prev_line].rend());
							hexline.insert(hexline.end(), hexlines[prev_line].rbegin() + 1, hexlines[prev_line].rend());
						}
					}
				}
			} while (find);

			trinewlines.push_back(newline);
			hexnewlines.push_back(hexline);
		}

		for (size_t i = 0; i < trilines.size(); i++)
		{
			if (trilinetag[i]) continue;

			std::vector<size_t> newline, hexline;
			newline.assign(trilines[i].begin(), trilines[i].end());
			hexline.assign(hexlines[i].begin(), hexlines[i].end());
			trilinetag[i] = true;
			size_t prev_line = i;
			bool find = false;
			do
			{
				find = false;
				size_t linebegin = newline.back();
				const std::vector<size_t>& shortline = trivertex2line[linebegin];
				size_t degree = shortline.size();
				if (degree == 2)
				{
					if (shortline.back() == prev_line)
					{
						if (!trilinetag[shortline.front()])
						{
							prev_line = shortline.front();
							find = true;
						}
					}
					else
					{
						if (!trilinetag[shortline.back()])
						{
							prev_line = shortline.back();
							find = true;
						}
					}
					if (find)
					{
						trilinetag[prev_line] = true;
						if (trilines[prev_line].front() == linebegin)
						{
							newline.insert(newline.end(), trilines[prev_line].begin() + 1, trilines[prev_line].end());
							hexline.insert(hexline.end(), hexlines[prev_line].begin() + 1, hexlines[prev_line].end());
						}
						else
						{
							newline.insert(newline.end(), trilines[prev_line].rbegin() + 1, trilines[prev_line].rend());
							hexline.insert(hexline.end(), hexlines[prev_line].rbegin() + 1, hexlines[prev_line].rend());
						}
					}
				}
			} while (find);

			trinewlines.push_back(newline);
			hexnewlines.push_back(hexline);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::laplacian_smooth_init()
	{
		std::vector< TinyVector<Real, 3> >& vertices = m_hexmesh->get_vertices();
		std::vector< HexElement<Real>* >& elements = m_hexmesh->get_elements();

		label_feature_edges();

		if (neighor_vertices.empty())
		{
			neighor_vertices.resize(vertices.size());
			for (size_t i = 0; i < elements.size(); i++)
			{
				for (int j = 0; j < 12; j++)
				{
					size_t id0 = elements[i]->indices[edgeid[j][0]];
					size_t id1 = elements[i]->indices[edgeid[j][1]];

					if (boundaryvertexlabel[id0] && boundaryvertexlabel[id1])
					{
						int tag0 = vertex_feature_tags[vert_map_ids[id0]];
						int tag1 = vertex_feature_tags[vert_map_ids[id1]];
						if (tag0 == -1 || (tag0 == 1 && tag1 > -1))
							neighor_vertices[id0].insert(id1);
						if (tag1 == -1 || (tag1 == 1 && tag0 > -1))
							neighor_vertices[id1].insert(id0);
					}
					else if (boundaryvertexlabel[id0] && boundaryvertexlabel[id1] == false)
					{
						neighor_vertices[id1].insert(id0);
					}
					else if (boundaryvertexlabel[id0] == false && boundaryvertexlabel[id1])
					{
						neighor_vertices[id0].insert(id1);
					}
					else
					{
						neighor_vertices[id0].insert(id1);
						neighor_vertices[id1].insert(id0);
					}
				}
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::fast_laplacian_smooth(unsigned int iter, bool only_boundary)
	{
		std::vector< TinyVector<Real, 3> >& vertices = m_hexmesh->get_vertices();
		std::vector< HexElement<Real>* >& elements = m_hexmesh->get_elements();

		Real TaubinWeight[2] = { 0.4507499669, -0.4720265626 };
		//if (feature_sensitive)
		TaubinWeight[0] = TaubinWeight[1] = 1;
		std::vector<TinyVector<Real, 3>> new_vertices;

		if (has_feature())
		{
#pragma omp parallel for
			for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
			{
				if (vertex_feature_tags[i] == 0)
					corner_projection(i, vertices[reverse_ids[i]], vertices[reverse_ids[i]]);
			}

			for (unsigned int miter = 0; miter < 2 * iter; miter++)
			{
				new_vertices.assign(vertices.size(), TinyVector<Real, 3>(0, 0, 0));
#pragma omp parallel for schedule(dynamic)
				for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
				{
					if (vertex_feature_tags[i] == 0) continue;
					HE_vert<Real>* hv = m_quadmesh->get_vertex(i);
					HE_edge<Real>* he = hv->edge;
					int d = 0;
					do
					{
						TinyVector<Real, 3> M = 0.5 * (vertices[reverse_ids[i]] + vertices[reverse_ids[he->vert->id]]);
						if (vertex_feature_tags[i] == 1)
							feature_projection(i, M, M);
						new_vertices[i] += M;
						d++;
						he = he->pair->next;
					} while (he != hv->edge);
					new_vertices[i] /= d;
				}

#pragma omp parallel for schedule(dynamic)
				for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
				{
					if (vertex_feature_tags[i] == 0) continue;

					new_vertices[i] = ((Real)1.0 - TaubinWeight[miter % 2]) * vertices[reverse_ids[i]] + TaubinWeight[miter % 2] * new_vertices[i];

					if (vertex_feature_tags[i] == 1)
						feature_projection(i, new_vertices[i], vertices[reverse_ids[i]]);
					else
					{
						if (miter < iter)
						{
							TinyVector<Real, 3> T;
							m_aabb->find_Nearest_Point(new_vertices[i], T);
							if ((vertices[reverse_ids[i]] - T).SquaredLength() < (vertices[reverse_ids[i]] - new_vertices[i]).SquaredLength())
								vertices[reverse_ids[i]] = T;
							else
								vertices[reverse_ids[i]] = new_vertices[i];
						}
						else
							m_aabb->find_Nearest_Point(new_vertices[i], vertices[reverse_ids[i]]);
					}
				}
			}

			//feature-aware
			for (unsigned int miter = 0; miter < 2 * iter; miter++)
			{
				new_vertices.assign(vertices.size(), TinyVector<Real, 3>(0, 0, 0));
#pragma omp parallel for schedule(dynamic)
				for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
				{
					if (vertex_feature_tags[i] != 1) continue;
					HE_vert<Real>* hv = m_quadmesh->get_vertex(i);
					HE_edge<Real>* he = hv->edge;
					int d = 0;
					do
					{
						if (featureedgelabels[he->id])
						{
							TinyVector<Real, 3> M = 0.5 * (vertices[reverse_ids[i]] + vertices[reverse_ids[he->vert->id]]);
							feature_projection(i, M, M);
							new_vertices[i] += M;
							d++;
						}
						he = he->pair->next;
					} while (he != hv->edge);
					new_vertices[i] /= d;
				}

#pragma omp parallel for schedule(dynamic)
				for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
				{
					if (vertex_feature_tags[i] != 1) continue;
					new_vertices[i] = ((Real)1.0 - TaubinWeight[miter % 2]) * vertices[reverse_ids[i]] + TaubinWeight[miter % 2] * new_vertices[i];
					feature_projection(i, new_vertices[i], vertices[reverse_ids[i]]);
				}
			}

			for (unsigned int miter = 0; miter < 2 * iter; miter++)
			{
				new_vertices.assign(vertices.size(), TinyVector<Real, 3>(0, 0, 0));
#pragma omp parallel for schedule(dynamic)
				for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
				{
					if (vertex_feature_tags[i] != -1) continue;
					HE_vert<Real>* hv = m_quadmesh->get_vertex(i);
					HE_edge<Real>* he = hv->edge;
					do
					{
						TinyVector<Real, 3> M = 0.5 * (vertices[reverse_ids[i]] + vertices[reverse_ids[he->vert->id]]);
						new_vertices[i] += M;
						he = he->pair->next;
					} while (he != hv->edge);
					new_vertices[i] /= hv->degree;
				}

#pragma omp parallel for schedule(dynamic)
				for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
				{
					if (vertex_feature_tags[i] != -1) continue;
					new_vertices[i] = ((Real)1.0 - TaubinWeight[miter % 2]) * vertices[reverse_ids[i]] + TaubinWeight[miter % 2] * new_vertices[i];
					m_aabb->find_Nearest_Point(new_vertices[i], vertices[reverse_ids[i]]);
				}
			}
		}
		else
		{
			for (unsigned int miter = 0; miter < 2 * iter; miter++)
			{
				new_vertices.assign(vertices.size(), TinyVector<Real, 3>(0, 0, 0));
#pragma omp parallel for schedule(dynamic)
				for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
				{
					HE_vert<Real>* hv = m_quadmesh->get_vertex(i);
					HE_edge<Real>* he = hv->edge;
					do
					{
						TinyVector<Real, 3> M = 0.5 * (vertices[reverse_ids[i]] + vertices[reverse_ids[he->vert->id]]);
						m_aabb->find_Nearest_Point(M, M);
						new_vertices[i] += M;
						he = he->pair->next;
					} while (he != hv->edge);
					new_vertices[i] /= hv->degree;
				}

#pragma omp parallel for schedule(dynamic)
				for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
				{
					if (vertex_feature_tags[i] != -1) continue;
					new_vertices[i] = ((Real)1.0 - TaubinWeight[miter % 2]) * vertices[reverse_ids[i]] + TaubinWeight[miter % 2] * new_vertices[i];
					m_aabb->find_Nearest_Point(new_vertices[i], vertices[reverse_ids[i]]);
				}
			}
		}

		if (!only_boundary)
		{
			//interior laplacian smooth
			TaubinWeight[0] = TaubinWeight[1] = 1;
			for (unsigned int miter = 0; miter < 2 * iter; miter++)
			{
				new_vertices.assign(vertices.size(), TinyVector<Real, 3>(0, 0, 0));
#pragma omp parallel for schedule(dynamic)
				for (ptrdiff_t i = 0; i < (ptrdiff_t)vertices.size(); i++)
				{
					if (boundaryvertexlabel[i]) continue;
					for (std::unordered_set<size_t>::iterator iter = neighor_vertices[i].begin(); iter != neighor_vertices[i].end(); iter++)
						new_vertices[i] += 0.5 * (vertices[i] + vertices[*iter]);
				}
#pragma omp parallel for schedule(dynamic)
				for (ptrdiff_t i = 0; i < (ptrdiff_t)vertices.size(); i++)
				{
					if (boundaryvertexlabel[i]) continue;
					vertices[i] = ((Real)1.0 - TaubinWeight[miter % 2]) * vertices[i] + TaubinWeight[miter % 2] / neighor_vertices[i].size() * new_vertices[i];
				}
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::laplacian_smooth(unsigned int iter, bool only_boundary, bool projection_only)
	{
		std::vector< TinyVector<Real, 3> >& vertices = m_hexmesh->get_vertices();
		std::vector< HexElement<Real>* >& elements = m_hexmesh->get_elements();

		//projection
#pragma omp parallel for
		for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
		{
			if (vertex_feature_tags[i] == 0)
				corner_projection(i, vertices[reverse_ids[i]], vertices[reverse_ids[i]]);
			//else if (vertex_feature_tags[i] == 1)
			//	feature_projection(i, vertices[reverse_ids[i]], vertices[reverse_ids[i]]);
			//else
			//	m_aabb->find_Nearest_Point(vertices[reverse_ids[i]], vertices[reverse_ids[i]]);
		}

		if (projection_only)
			return;

		////////////////////////////////////

		Real TaubinWeight[2] = { 0.4507499669, -0.4720265626 };
		if (feature_sensitive)
			TaubinWeight[0] = TaubinWeight[1] = 1;
		std::vector<TinyVector<Real, 3>> new_vertices;

		bool smooth_all = true, struct_smoothing = true, feature_smoothing = true, non_feature_smoothing = true;

		if (!has_feature())
			struct_smoothing = feature_smoothing = non_feature_smoothing = false;
		else
			;// smooth_all = struct_smoothing = false;

		//surface laplacian smooth
		if (!fix_boundary)
		{
			//smooth all vertices
			if (smooth_all)
			{
				for (unsigned int miter = 0; miter < 2 * iter; miter++)
				{
					new_vertices.assign(vertices.size(), TinyVector<Real, 3>(0, 0, 0));
#pragma omp parallel for schedule(dynamic)
					for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
					{
						if (vertex_feature_tags[i] == 0) continue;
						HE_vert<Real>* hv = m_quadmesh->get_vertex(i);
						HE_edge<Real>* he = hv->edge;
						int d = 0;
						do
						{
							TinyVector<Real, 3> M = 0.5 * (vertices[reverse_ids[i]] + vertices[reverse_ids[he->vert->id]]);
							new_vertices[i] += M;
							d++;
							he = he->pair->next;
						} while (he != hv->edge);
						new_vertices[i] /= d;
					}

#pragma omp parallel for schedule(dynamic)
					for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
					{
						if (vertex_feature_tags[i] == 0) continue;
						new_vertices[i] = ((Real)1.0 - TaubinWeight[miter % 2]) * vertices[reverse_ids[i]] + TaubinWeight[miter % 2] * new_vertices[i];

						if (vertex_feature_tags[i] == 1)
							feature_projection(i, new_vertices[i], vertices[reverse_ids[i]]);
						else
						{
							if (miter < iter)
							{
								TinyVector<Real, 3> T;
								m_aabb->find_Nearest_Point(new_vertices[i], T);
								if ((vertices[reverse_ids[i]] - T).SquaredLength() < (vertices[reverse_ids[i]] - new_vertices[i]).SquaredLength())
									vertices[reverse_ids[i]] = T;
								else
									vertices[reverse_ids[i]] = new_vertices[i];
							}
							else
								m_aabb->find_Nearest_Point(new_vertices[i], vertices[reverse_ids[i]]);
						}
					}
				}
			}

			//structline smoothing
			if (struct_smoothing)
			{
				std::vector<bool> irregularvertexlabel;
				extract_irrugular_lines(irregularvertexlabel, true);
#pragma omp parallel for
				for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
				{
					if (vertex_feature_tags[i] != -1)
						irregularvertexlabel[i] = true;
				}

				unsigned int struct_iter = 50;
				for (unsigned int miter = 0; miter < 2 * struct_iter; miter++)
				{
					new_vertices.assign(vertices.size(), TinyVector<Real, 3>(0, 0, 0));
#pragma omp parallel for schedule(dynamic)
					for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
					{
						if (!irregularvertexlabel[i] || vertex_feature_tags[i] == 0) continue;
						HE_vert<Real>* hv = m_quadmesh->get_vertex(i);
						HE_edge<Real>* he = hv->edge;
						int d = 0;
						do
						{
							if ((vertex_feature_tags[i] == -1 && he->tag) ||
								(vertex_feature_tags[i] == 1 && featureedgelabels[he->id]))
							{
								TinyVector<Real, 3> M = 0.5 * (vertices[reverse_ids[i]] + vertices[reverse_ids[he->vert->id]]);
								if (featureedgelabels[he->id])
									feature_projection(i, M, M);
								//else
								//	m_aabb->find_Nearest_Point(M, M);
								new_vertices[i] += M;

								//new_vertices[i] += vertices[reverse_ids[he->vert->id]];
								d++;
							}
							he = he->pair->next;
						} while (he != hv->edge);
						new_vertices[i] /= d;
					}

#pragma omp parallel for schedule(dynamic)
					for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
					{
						if (!irregularvertexlabel[i] || vertex_feature_tags[i] == 0) continue;
						new_vertices[i] = ((Real)1.0 - TaubinWeight[miter % 2]) * vertices[reverse_ids[i]] + TaubinWeight[miter % 2] * new_vertices[i];
						if (vertex_feature_tags[i] == 1)
							feature_projection(i, new_vertices[i], vertices[reverse_ids[i]]);
						else
							m_aabb->find_Nearest_Point(new_vertices[i], vertices[reverse_ids[i]]);
					}
				}
				for (unsigned int miter = 0; miter < 2 * struct_iter; miter++)
				{
					new_vertices.assign(vertices.size(), TinyVector<Real, 3>(0, 0, 0));
#pragma omp parallel for schedule(dynamic)
					for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
					{
						if (irregularvertexlabel[i]) continue;
						HE_vert<Real>* hv = m_quadmesh->get_vertex(i);
						HE_edge<Real>* he = hv->edge;
						do
						{
							TinyVector<Real, 3> M = 0.5 * (vertices[reverse_ids[i]] + vertices[reverse_ids[he->vert->id]]);
							//m_aabb->find_Nearest_Point(M, M);
							new_vertices[i] += M;

							//new_vertices[i] += vertices[reverse_ids[he->vert->id]];
							he = he->pair->next;
						} while (he != hv->edge);
						new_vertices[i] /= hv->degree;
					}

#pragma omp parallel for schedule(dynamic)
					for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
					{
						if (irregularvertexlabel[i]) continue;
						new_vertices[i] = ((Real)1.0 - TaubinWeight[miter % 2]) * vertices[reverse_ids[i]] + TaubinWeight[miter % 2] * new_vertices[i];
						m_aabb->find_Nearest_Point(new_vertices[i], vertices[reverse_ids[i]]);
					}
				}
			}

			//feature smoothing & non-feature smoothing
			if (feature_smoothing)
			{
				for (unsigned int miter = 0; miter < 2 * iter; miter++)
				{
					new_vertices.assign(vertices.size(), TinyVector<Real, 3>(0, 0, 0));
#pragma omp parallel for schedule(dynamic)
					for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
					{
						if (vertex_feature_tags[i] != 1) continue;
						HE_vert<Real>* hv = m_quadmesh->get_vertex(i);
						HE_edge<Real>* he = hv->edge;
						int d = 0;
						do
						{
							if (featureedgelabels[he->id])
							{
								TinyVector<Real, 3> M = 0.5 * (vertices[reverse_ids[i]] + vertices[reverse_ids[he->vert->id]]);
								feature_projection(i, M, M);
								new_vertices[i] += M;

								//new_vertices[i] += vertices[reverse_ids[he->vert->id]];
								d++;
							}
							he = he->pair->next;
						} while (he != hv->edge);
						new_vertices[i] /= d;
					}

#pragma omp parallel for schedule(dynamic)
					for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
					{
						if (vertex_feature_tags[i] != 1) continue;
						new_vertices[i] = ((Real)1.0 - TaubinWeight[miter % 2]) * vertices[reverse_ids[i]] + TaubinWeight[miter % 2] * new_vertices[i];
						feature_projection(i, new_vertices[i], vertices[reverse_ids[i]]);
					}
				}

				for (unsigned int miter = 0; miter < 2 * iter; miter++)
				{
					new_vertices.assign(vertices.size(), TinyVector<Real, 3>(0, 0, 0));
#pragma omp parallel for schedule(dynamic)
					for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
					{
						if (vertex_feature_tags[i] != -1) continue;
						HE_vert<Real>* hv = m_quadmesh->get_vertex(i);
						HE_edge<Real>* he = hv->edge;
						do
						{
							TinyVector<Real, 3> M = 0.5 * (vertices[reverse_ids[i]] + vertices[reverse_ids[he->vert->id]]);
							//m_aabb->find_Nearest_Point(M, M);
							new_vertices[i] += M;
							//new_vertices[i] += vertices[reverse_ids[he->vert->id]];
							he = he->pair->next;
						} while (he != hv->edge);
						new_vertices[i] /= hv->degree;
					}

#pragma omp parallel for schedule(dynamic)
					for (ptrdiff_t i = 0; i < m_quadmesh->get_num_of_vertices(); i++)
					{
						if (vertex_feature_tags[i] != -1) continue;
						new_vertices[i] = ((Real)1.0 - TaubinWeight[miter % 2]) * vertices[reverse_ids[i]] + TaubinWeight[miter % 2] * new_vertices[i];
						m_aabb->find_Nearest_Point(new_vertices[i], vertices[reverse_ids[i]]);
					}
				}
			}
		}

		if (!only_boundary)
		{
			//interior laplacian smooth
			TaubinWeight[0] = TaubinWeight[1] = 1;
			for (unsigned int miter = 0; miter < 2 * iter; miter++)
			{
				new_vertices.assign(vertices.size(), TinyVector<Real, 3>(0, 0, 0));
#pragma omp parallel for schedule(dynamic)
				for (ptrdiff_t i = 0; i < (ptrdiff_t)vertices.size(); i++)
				{
					if (boundaryvertexlabel[i]) continue;
					for (std::unordered_set<size_t>::iterator iter = neighor_vertices[i].begin(); iter != neighor_vertices[i].end(); iter++)
						new_vertices[i] += vertices[*iter];
				}
#pragma omp parallel for schedule(dynamic)
				for (ptrdiff_t i = 0; i < (ptrdiff_t)vertices.size(); i++)
				{
					if (boundaryvertexlabel[i]) continue;
					vertices[i] = ((Real)1.0 - TaubinWeight[miter % 2]) * vertices[i] + TaubinWeight[miter % 2] / neighor_vertices[i].size() * new_vertices[i];
				}
			}
		}
	}
	/////////////////////////////////////////////////////////////
	template <typename Real>
	bool HexaSmoothBase<Real>::has_feature()
	{
		std::cout << "please override this function" << std::endl;
		return false;
	}
	/////////////////////////////////////////////////////////////
	template <typename Real>
	void HexaSmoothBase<Real>::label_feature_edges()
	{
		featureedgelabels.assign(m_quadmesh->get_num_of_edges(), false);
		for (size_t i = 0; i < feature_edges.size(); i += 2)
		{
			HE_vert<Real>* hv = m_quadmesh->get_vertex(feature_edges[i]);
			HE_vert<Real>* hv2 = m_quadmesh->get_vertex(feature_edges[i + 1]);
			HE_edge<Real>* he = hv->edge;
			do
			{
				if (he->vert == hv2)
				{
					featureedgelabels[he->id] = featureedgelabels[he->pair->id] = true;
					break;
				}
				he = he->pair->next;
			} while (he != hv->edge);
		}
	}
	/////////////////////////////////////////////////////////////
	template class HexaSmoothBase<double>;
}