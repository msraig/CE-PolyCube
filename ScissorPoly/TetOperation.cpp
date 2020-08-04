#include "TetOperation.h"
#include <iostream>

namespace ig
{

	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	bool TetOperation<Real>::split_edge(Tetrahedron<Real>* tet, unsigned int vi, unsigned int vj, bool split_boundary)
	{
		TetVertex<Real>* nvi = tet->vertex[vi];
		TetVertex<Real>* nvj = tet->vertex[vj];

		bool is_boundary_edge = false;
		if (nvi->boundary && nvj->boundary)
		{
			if (split_boundary)
			{
				is_boundary_edge = true;
			}
		}

		std::vector<Tetrahedron<Real>*> adjtets;
		std::vector<TetVertex<Real>*> loop;
		tetstructure->find_adjacent_tets_by_edge(tet, vi, vj, adjtets, loop);
		TetVertex<Real>* nv = new TetVertex<Real>((nvi->pos + nvj->pos) / 2);
		nv->id = tetstructure->tetra_vertices.size();
		tetstructure->tetra_vertices.push_back(nv);
		if (is_boundary_edge)
		{
			nv->boundary = true;
		}
		//std::vector<Real> mvol;
		//for (size_t i = 0; i < adjtets.size(); i++)
		//{
		//	mvol.push_back(adjtets[i]->volume);
		//	assert(adjtets[i]->volume >= 0);
		//}

		TetVertex<Real>* dv[4];
		for (int i = 0; i < 4; i++)
			dv[i] = adjtets[0]->vertex[i];
		for (int i = 0; i < 4; i++)
		{
			if (adjtets[0]->vertex[i] == nvi && i != 0)
			{
				std::swap(dv[i], dv[0]);
				if (i == 1)
				{
					std::swap(dv[2], dv[3]);
				}
				else if (i == 2)
				{
					std::swap(dv[1], dv[3]);
				}
				else if (i == 3)
				{
					std::swap(dv[1], dv[2]);
				}
				break;
			}
		}
		if ((dv[1] == nvj && dv[2] == loop[0] && dv[3] == loop[1]) || (dv[3] == nvj && dv[1] == loop[0] && dv[2] == loop[1])
			|| (dv[2] == nvj && dv[1] == loop[1] && dv[3] == loop[0]))
		{
			;
		}
		else
		{
			std::swap(nvi, nvj);
		}

		//for (size_t i = 0; i < adjtets.size(); i++)
		//{
		//	DF_Tetrahedron<Real> vol_test;
		//	vol_test.vertex[0] = nvi;
		//	vol_test.vertex[1] = nvj;
		//	vol_test.vertex[2] = loop[i];
		//	vol_test.vertex[3] = loop[(i+1)];
		//	vol_test.compute_tet_volume();
		//	assert(vol_test.volume >= 0);
		//}

		//DF_Tetrahedron<Real> Vol_test(nvi, nvj, loop[0], loop[1]);

		//if (Vol_test.vertex[0] != nvi)
		//{
		//	std::swap(nvi, nvj);
		//}

		std::vector<Tetrahedron<Real>*> T(2 * adjtets.size());
		std::vector<int> matching(4 * adjtets.size());
		std::vector<unsigned int> localid(4 * adjtets.size()), opid(2 * adjtets.size());
		std::vector<Tetrahedron<Real>*> adjT(2 * adjtets.size());
		for (size_t i = 0; i < adjtets.size(); i++)
		{
			localid[4 * i] = adjtets[i]->get_local_vertex_id(nvi);
			localid[4 * i + 1] = adjtets[i]->get_local_vertex_id(nvj);
			localid[4 * i + 2] = adjtets[i]->get_local_vertex_id(loop[i]);
			localid[4 * i + 3] = adjtets[i]->get_local_vertex_id(loop[i + 1]);
			adjtets[i]->get_opposite_vertex(localid[4 * i], &opid[2 * i]);
			adjtets[i]->get_opposite_vertex(localid[4 * i + 1], &opid[2 * i + 1]);
			adjT[2 * i] = adjtets[i]->neighborTet[localid[4 * i]];
			adjT[2 * i + 1] = adjtets[i]->neighborTet[localid[4 * i + 1]];
			matching[4 * i] = cell_matching[4 * adjtets[i]->id + localid[4 * i]];
			matching[4 * i + 1] = cell_matching[4 * adjtets[i]->id + localid[4 * i + 1]];
			matching[4 * i + 2] = cell_matching[4 * adjtets[i]->id + localid[4 * i + 2]];
			matching[4 * i + 3] = cell_matching[4 * adjtets[i]->id + localid[4 * i + 3]];
		}

		for (size_t i = 0; i < adjtets.size(); i++)
		{
			T[2 * i] = new Tetrahedron<Real>;
			T[2 * i + 1] = adjtets[i];
			//T[2 * i]->layer = adjtets[i]->layer;
			T[2 * i]->vertex[0] = nvi; T[2 * i]->vertex[1] = nv; T[2 * i]->vertex[2] = loop[i]; T[2 * i]->vertex[3] = loop[i + 1];
			T[2 * i + 1]->vertex[0] = nv; T[2 * i + 1]->vertex[1] = nvj; T[2 * i + 1]->vertex[2] = loop[i]; T[2 * i + 1]->vertex[3] = loop[i + 1];
			T[2 * i]->id = tetstructure->tetras.size();	tetstructure->tetras.push_back(T[2 * i]);
		}

		size_t ns = adjtets.size();
		cell_matching.resize(cell_matching.size() + ns * 4);
		for (size_t i = 0; i < adjtets.size(); i++)
		{
			size_t b = (i - 1 + ns) % ns;
			size_t n = (i + 1) % ns;
			T[2 * i]->neighborTet[0] = T[2 * i + 1]; T[2 * i]->neighborTet[1] = adjT[2 * i + 1]; T[2 * i]->neighborTet[2] = is_boundary_edge && i == ns - 1 ? 0 : T[2 * n]; T[2 * i]->neighborTet[3] = is_boundary_edge && i == 0 ? 0 : T[2 * b];
			T[2 * i + 1]->neighborTet[0] = adjT[2 * i]; T[2 * i + 1]->neighborTet[1] = T[2 * i]; T[2 * i + 1]->neighborTet[2] = is_boundary_edge && i == ns - 1 ? 0 : T[2 * n + 1]; T[2 * i + 1]->neighborTet[3] = is_boundary_edge && i == 0 ? 0 : T[2 * b + 1];
			if (adjT[2 * i + 1])
			{
				adjT[2 * i + 1]->neighborTet[opid[2 * i + 1]] = T[2 * i];
			}
			if (adjT[2 * i])
			{
				adjT[2 * i]->neighborTet[opid[2 * i]] = T[2 * i + 1];
			}
			cell_matching[T[2 * i]->id * 4] = 0; cell_matching[T[2 * i]->id * 4 + 1] = matching[4 * i + 1]; cell_matching[T[2 * i]->id * 4 + 2] = matching[4 * i + 2]; cell_matching[T[2 * i]->id * 4 + 3] = matching[4 * i + 3];
			cell_matching[T[2 * i + 1]->id * 4] = matching[4 * i]; cell_matching[T[2 * i + 1]->id * 4 + 1] = 0; cell_matching[T[2 * i + 1]->id * 4 + 2] = matching[4 * i + 2]; cell_matching[T[2 * i + 1]->id * 4 + 3] = matching[4 * i + 3];
		}


		return true;
	}

	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	bool TetOperation<Real>::split_inner_edge_two_boundary_point(TetStructure<Real> &tet_mesh)
	{
		std::vector<ig::CVec<Real, 3>> points;
		std::vector<unsigned int> faces;

		int nv = (int)tet_mesh.tetra_vertices.size();
		int nc = (int)tet_mesh.tetras.size();
		std::vector<std::vector<unsigned int>> faces_array;
		std::vector<unsigned int> tmp_face;
		std::vector<Tetrahedron<Real>*> *tetras = &tet_mesh.get_tetras();



		//second step of processing
		std::set<std::pair<int, int>> boundary_edge_set; //first less than second, keep untouched

		//get_boundary_edge_set
		get_boundary_edge_set(tet_mesh, boundary_edge_set);


		std::set<std::pair<int, int>> error_edge_pairs;
		std::set<int> error_tet_set;

		for (int i = 0; i < nc; i++)
		{
			//if ((*tetras)[i]->is_onboundary()) continue;
			for (int j = 0; j < 3; j++)
			{
				for (size_t k = j + 1; k < 4; k++)
				{
					int vid0 = (int)(*tetras)[i]->vertex[j]->id;
					int vid1 = (int)(*tetras)[i]->vertex[k]->id;
					int minid = (int)std::min(vid0, vid1);
					int maxid = (int)std::max(vid0, vid1);
					if (tet_mesh.tetra_vertices[vid0]->boundary && tet_mesh.tetra_vertices[vid1]->boundary)
					{
						//both vert on boundary, not boundary edge
						std::pair<int, int> tmp_pair(minid, maxid);
						auto it = boundary_edge_set.find(tmp_pair);
						if (it == boundary_edge_set.end())
						{
							//not boundary edge
							error_edge_pairs.insert(tmp_pair);
							error_tet_set.insert(i);
						}
					}
				}
			}
		}
		std::cout << "Error Edge Pair: " << error_edge_pairs.size() << std::endl;

		std::vector<std::pair<int, int>> onecell_ee_lid;
		while (!error_tet_set.empty())
		{
			auto begin = error_tet_set.begin();
			//error edges
			onecell_ee_lid.clear();
			std::pair<int, int> first_delete_e(-1, -1);

			for (size_t i = 0; i < 3; i++)
			{
				for (size_t j = i + 1; j < 4; j++)
				{
					int vid0 = (int)(*tetras)[*begin]->vertex[i]->id;
					int vid1 = (int)(*tetras)[*begin]->vertex[j]->id;
					int minid = (int)std::min(vid0, vid1);
					int maxid = (int)std::max(vid0, vid1);
					if (tet_mesh.tetra_vertices[vid0]->boundary && tet_mesh.tetra_vertices[vid1]->boundary)
					{
						//both vert on boundary, not boundary edge
						std::pair<int, int> tmp_pair(minid, maxid);
						auto it = error_edge_pairs.find(tmp_pair);
						if (it != error_edge_pairs.end())
						{
							//error edge
							if (first_delete_e.first == -1)
								first_delete_e = *it;
							onecell_ee_lid.push_back(std::pair<int, int>(i, j));
						}
					}
				}
			}

			if (onecell_ee_lid.empty())
			{
				//empty, delete tet from set
				error_tet_set.erase(begin);
			}
			else
			{
				int cur_nc = (int)tetras->size();
				//delete the first edge
				error_edge_pairs.erase(first_delete_e);
				//split
				split_edge((*tetras)[*begin], onecell_ee_lid[0].first, onecell_ee_lid[0].second, false);

				if (onecell_ee_lid.size() == 1)
				{
					error_tet_set.erase(begin);
				}
				else
				{
					//more than one edge
					error_tet_set.insert(cur_nc);
				}

			}

		}

		//deal with remaining part of edges
		while (!error_edge_pairs.empty())
		{
			auto it = error_edge_pairs.begin();
			int tetid, local_id[2];
			get_edge_local_id(*tetras, *it, tetid, local_id);
			assert(tetid != -1);
			//split
			split_edge((*tetras)[tetid], local_id[0], local_id[1], false);

			error_edge_pairs.erase(it);
		}
		std::cout << "Split Edge finished" << std::endl;
		return true;
	}

	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	bool TetOperation<Real>::split_cell_four_boundary_point(TetStructure<Real> &tet_mesh)
	{
		std::vector<ig::CVec<Real, 3>> points;
		std::vector<unsigned int> faces;

		std::vector<std::vector<unsigned int>> faces_array;
		std::vector<unsigned int> tmp_face;
		std::vector<Tetrahedron<Real>*> *tetras = &tet_mesh.get_tetras();

		get_vert_cell(&tet_mesh, points, faces);
		unsigned int nc = (unsigned int)tet_mesh.tetras.size();
		unsigned int nv = (unsigned int)tet_mesh.tetra_vertices.size();
		std::set<int> del_face_id;
		for (int i = 0; i < (int)nc; i++)
		{
			tmp_face.clear();
			int boundary_count = 0;
			for (int j = 0; j < 4; j++)
			{
				tmp_face.push_back(faces[4 * i + j]);
				if ((*tetras)[i]->vertex[j]->boundary)
				{
					boundary_count++;
				}
			}
			faces_array.push_back(tmp_face);
			if (boundary_count == 4)
				del_face_id.insert(i);
		}

		std::cout << "deleled face size: " << del_face_id.size() << std::endl;

		std::vector<unsigned int> new_faces;
		for (int i = 0; i < (int)nc; i++)
		{
			auto it = del_face_id.find(i);
			if (it == del_face_id.end())
			{
				//not found
				for (size_t j = 0; j < 4; j++)
				{
					new_faces.push_back(faces[4 * i + j]);
				}
			}
		}

		static int s_tet_id[4][3] = { { 1, 2, 3}, {2, 0, 3}, {0, 1, 3}, {1, 0, 2} }; //to outer part

		for (auto it : del_face_id)
		{
			int new_pt_id = (int)points.size();
			CVec<Real, 3> new_pt(0.0, 0.0, 0.0);

			for (size_t i = 0; i < 4; i++)
			{
				new_pt += points[faces_array[it][i]];
			}
			new_pt = new_pt / 4.0;
			points.push_back(new_pt);

			for (size_t i = 0; i < 4; i++)
			{
				for (size_t j = 0; j < 3; j++)
				{
					new_faces.push_back(faces_array[it][s_tet_id[i][2 - j]]);
				}
				new_faces.push_back(new_pt_id);
			}

		}

		std::cout << "Split Tet Finished" << std::endl;

		tet_mesh.load_tet(points, new_faces);
		
		return true;
	}

	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void TetOperation<Real>::get_vert_cell(TetStructure<Real>* tet_mesh, std::vector<CVec<Real, 3>> &points, std::vector<unsigned int> &indices)
	{
		points.clear();
		indices.clear();

		const std::vector<Tetrahedron<Real>*> &tetras = tet_mesh->tetras;
		const std::vector<TetVertex<Real>*> &tetra_vertices = tet_mesh->tetra_vertices;

		int nv = (int)tetra_vertices.size();
		int nc = (int)tetras.size();

		for (size_t i = 0; i < nv; i++)
		{
			CVec<Real, 3> temp = tetra_vertices[i]->pos;
			points.push_back(temp);
		}

		for (size_t i = 0; i < nc; i++)
		{
			for (size_t j = 0; j < 4; j++)
			{
				unsigned id = (unsigned)tetras[i]->vertex[j]->id;
				indices.push_back(id);
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void TetOperation<Real>::get_boundary_vert_face(TetStructure<Real>* tet_mesh, std::vector<Real> &points, std::vector<int> &faces, std::vector<int> &s2v, bool ori_vert_flag, std::map<std::vector<int>, int> *face2cellidx)
	{
		//s2v: surface to volume vert map
		points.clear();
		faces.clear();
		s2v.clear();
		face2cellidx->clear();
		//VolumeSurfaceVertex.clear();
		std::vector<int> tmp_face;
		std::map<int, int> VolumeSurfaceVertex;
		std::map<int, int>::iterator map_it;
		int indexOnSurface = 0;

		int num_tet = (int)tet_mesh->tetras.size();
		std::vector<Tetrahedron<Real>*> tetras = tet_mesh->get_tetras();

		static int s_tet_id[4][3] = { { 1, 2, 3}, {2, 0, 3}, {0, 1, 3}, {1, 0, 2} }; //to outer part

		if (!ori_vert_flag)
		{
			for (size_t i = 0; i < num_tet; i++)
			{
				if (!tetras[i]->is_onboundary())
				{
					continue;
				}
				//BoundaryFaceFlag[i] = true;

				int face_id = -1;
				for (size_t j = 0; j < 4; j++)
				{
					TetVertex<Real>* tv_temp = tetras[i]->vertex[j];
					if (!(tv_temp->boundary))
					{
						face_id = (int)j;
						break;
					}
				}

				if (face_id == -1)
				{
					std::cout << "read vtk error" << std::endl;
					return;
				}

				//faceVertexHandle.clear();
				tmp_face.clear();
				for (size_t j = 0; j < 3; j++)
				{
					map_it = VolumeSurfaceVertex.find((int)tetras[i]->vertex[s_tet_id[face_id][j]]->id);
					if (map_it == VolumeSurfaceVertex.end())
					{
						VolumeSurfaceVertex.insert(std::pair<int, int>(tetras[i]->vertex[s_tet_id[face_id][j]]->id, indexOnSurface));
						s2v.push_back((int)tetras[i]->vertex[s_tet_id[face_id][j]]->id);
						Real px, py, pz;
						px = tetras[i]->vertex[s_tet_id[face_id][j]]->pos[0];
						py = tetras[i]->vertex[s_tet_id[face_id][j]]->pos[1];
						pz = tetras[i]->vertex[s_tet_id[face_id][j]]->pos[2];
						//SurfaceMesh::Point v2(px, py, pz);
						//SurfaceMesh::Point v2((px - cx)/max_diff, (py - cy)/max_diff, (pz - cz)/max_diff);
						points.push_back(px);
						points.push_back(py);
						points.push_back(pz);
						//vertexHandleVec.push_back(boundaryMesh.add_vertex(v2));
						//boundaryMesh.data(vertexHandleVec.back()).set_tet_vertex_id(tetras[i]->vertex[s_tet_id[face_id][j]]->id);
						//faceVertexHandle.push_back(vertexHandleVec.back());
						tmp_face.push_back(indexOnSurface);
						++indexOnSurface;
					}
					else
					{
						//faceVertexHandle.push_back(vertexHandleVec[map_it->second]);
						tmp_face.push_back(map_it->second);
					}
				}
				//OpenMesh::FaceHandle fh = boundaryMesh.add_face(faceVertexHandle);

				//boundaryMesh.data(fh).set_cell_id(i);

				faces.push_back(tmp_face[0]);
				faces.push_back(tmp_face[1]);
				faces.push_back(tmp_face[2]);


			}
		}
		else
		{
			//initialize points
			for (size_t i = 0; i < tet_mesh->tetra_vertices.size(); i++)
			{
				Real px, py, pz;
				px = tet_mesh->get_vertices()[i]->pos[0];
				py = tet_mesh->get_vertices()[i]->pos[1];
				pz = tet_mesh->get_vertices()[i]->pos[2];
				points.push_back(px);
				points.push_back(py);
				points.push_back(pz);
			}

			for (size_t i = 0; i < num_tet; i++)
			{
				if (!tetras[i]->is_onboundary())
				{
					continue;
				}
				//BoundaryFaceFlag[i] = true;

				int face_id = -1;
				for (size_t j = 0; j < 4; j++)
				{
					TetVertex<Real>* tv_temp = tetras[i]->vertex[j];
					if (!(tv_temp->boundary))
					{
						face_id = (int)j;
						break;
					}
				}

				if (face_id == -1)
				{
					std::cout << "read vtk error" << std::endl;
					return;
				}

				//faceVertexHandle.clear();
				tmp_face.clear();
				for (size_t j = 0; j < 3; j++)
				{
					tmp_face.push_back((int)tetras[i]->vertex[s_tet_id[face_id][j]]->id);
				}
				//OpenMesh::FaceHandle fh = boundaryMesh.add_face(faceVertexHandle);

				//boundaryMesh.data(fh).set_cell_id(i);

				faces.push_back(tmp_face[0]);
				faces.push_back(tmp_face[1]);
				faces.push_back(tmp_face[2]);

				if (face2cellidx != NULL)
				{
					//reorder tmp_face;
					sort(tmp_face.begin(), tmp_face.end());
					(*face2cellidx)[tmp_face] = (int)i;
				}
			}


		}

		//boundaryMesh.update_normals();


	}



	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void TetOperation<Real>::get_boundary_edge_set(const TetStructure<Real> &tet_mesh, std::set<std::pair<int, int>> &boundary_edge_set)
	{
		boundary_edge_set.clear();
		int nc = (int)tet_mesh.tetras.size();
		int nv = (int)tet_mesh.tetra_vertices.size();

		static int s_tet_id[4][3] = { { 1, 2, 3}, {2, 0, 3}, {0, 1, 3}, {1, 0, 2} };

		for (size_t i = 0; i < nc; i++)
		{
			if (tet_mesh.tetras[i]->is_onboundary())
			{
				//iterate over faces
				for (size_t j = 0; j < 4; j++)
				{
					if (tet_mesh.tetras[i]->neighborTet[j] == 0)
					{
						//boundary face
						for (size_t k1 = 0; k1 < 2; k1++)
						{
							for (size_t k2 = k1 + 1; k2 < 3; k2++)
							{
								int vid0 = (int)tet_mesh.tetras[i]->vertex[s_tet_id[j][k1]]->id;
								int vid1 = (int)tet_mesh.tetras[i]->vertex[s_tet_id[j][k2]]->id;
								int minid = (int)std::min(vid0, vid1);
								int maxid = (int)std::max(vid0, vid1);
								assert(minid != maxid);
								boundary_edge_set.insert(std::pair<int, int>(minid, maxid));
							}
						}

					}
				}
			}
		}


	}

	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void TetOperation<Real>::get_edge_local_id(const std::vector<Tetrahedron<Real>* > &tetras, std::pair<int, int> edge, int &tet_id, int localid[])
	{
		tet_id = -1;
		int nc = (int)tetras.size();
		for (size_t i = 0; i < nc; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				for (size_t k = j + 1; k < 4; k++)
				{
					int vid0 = (int)(tetras)[i]->vertex[j]->id;
					int vid1 = (int)(tetras)[i]->vertex[k]->id;
					int minid = std::min(vid0, vid1);
					int maxid = std::max(vid0, vid1);
					if (minid == edge.first && maxid == edge.second)
					{
						tet_id = (int)i;
						localid[0] = (int)j;
						localid[1] = (int)k;
						return;
					}
				}
			}
		}

	}


	/////////////////////////////////////
	template class TetOperation<double>;

}