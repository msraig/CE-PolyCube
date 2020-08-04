#include "VHexMesh.h"

#include <fstream>
#include <string>
#include <map>
#include <queue>
#include <set>
#include "MyTuple.h"
#include "ANN/ANN.h"
#include <omp.h>

namespace MeshLib
{
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	VHexMesh<Real>::~VHexMesh()
	{
		clear_mem();
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::clear_mem()
	{
		for (size_t i = 0; i < elements.size(); i++)
		{
			delete elements[i];
		}
		vertices.resize(0);
		elements.resize(0);
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::load_hex(const char filename[])
	{
		std::ifstream hexfile(filename);
		if (hexfile.is_open())
		{
			clear_mem();
			std::string str;
			try
			{
				size_t nv, ne;
				hexfile >> nv >> str >> ne >> str;
				vertices.resize(nv);
				elements.reserve(ne);
				for (size_t i = 0; i < nv; i++)
				{
					hexfile >> vertices[i][0] >> vertices[i][1] >> vertices[i][2];
				}
				size_t c_id;
				size_t ind[8];
				for (size_t i = 0; i < ne; i++)
				{
					hexfile >> c_id;
					if (c_id == 8)
					{
						for (int j = 0; j < 8; j++)
						{
							hexfile >> ind[j];
						}
						std::swap(ind[0], ind[1]);
						std::swap(ind[4], ind[5]);
						elements.push_back(new HexElement<Real>(ind));
					}
				}
			}
			catch (...)
			{
				vertices.resize(0);
				elements.resize(0);
				hexfile.close();
				return;
			}

			hexfile.close();
			remove_unused_vertices();
			update_boundingbox();
			build_connectivity();
			reorientation();
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::update_boundingbox()
	{
		if (vertices.size() > 0)
		{
			boundingbox[0][0] = boundingbox[0][1] = vertices[0][0];
			boundingbox[1][0] = boundingbox[1][1] = vertices[0][1];
			boundingbox[2][0] = boundingbox[2][1] = vertices[0][2];
			for (size_t i = 1; i < (size_t)vertices.size(); i++)
			{
				boundingbox[0][0] = std::min(boundingbox[0][0], vertices[i][0]);
				boundingbox[0][1] = std::max(boundingbox[0][1], vertices[i][0]);
				boundingbox[1][0] = std::min(boundingbox[1][0], vertices[i][1]);
				boundingbox[1][1] = std::max(boundingbox[1][1], vertices[i][1]);
				boundingbox[2][0] = std::min(boundingbox[2][0], vertices[i][2]);
				boundingbox[2][1] = std::max(boundingbox[2][1], vertices[i][2]);
			}
			Real width;
			for (int j = 0; j < 3; j++)
			{
				width = boundingbox[j][1] - boundingbox[j][0];
				boundingbox[j][1] += 0.5 * width;
				boundingbox[j][0] -= 0.5 * width;
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::load_vtk(const char filename[])
	{
		std::ifstream hexfile(filename);
		if (hexfile.is_open())
		{
			clear_mem();
			std::string str;
			try
			{
				do
				{
					hexfile >> str;
				} while (str != "DATASET");
				hexfile >> str;
				if (str != "UNSTRUCTURED_GRID")
				{
					return;
				}
				do
				{
					hexfile >> str;
				} while (str != "POINTS");

				size_t nv, ne, num_ne_type;

				hexfile >> nv >> str;

				vertices.resize(nv);
				for (size_t i = 0; i < nv; i++)
				{
					hexfile >> vertices[i][0] >> vertices[i][1] >> vertices[i][2];
				}

				hexfile >> str >> ne >> num_ne_type;

				elements.reserve(ne);
				size_t ind[8];
				int index[8] = { 4, 5, 1, 0, 7, 6, 2, 3 };
				for (size_t i = 0; i < ne; i++)
				{
					hexfile >> ind[0]; // 8
					for (int j = 0; j < 8; j++)
					{
						hexfile >> ind[index[j]];
					}
					elements.push_back(new HexElement<Real>(ind));
				}
			}
			catch (...)
			{
				vertices.resize(0);
				elements.resize(0);
				hexfile.close();
				return;
			}

			hexfile.close();
			remove_unused_vertices();
			update_boundingbox();
			build_connectivity();
			reorientation();
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::save_hex(const char filename[])
	{
		std::ofstream hexfile(filename);
		if (hexfile.is_open())
		{
			hexfile << vertices.size() << " vertices" << std::endl;
			hexfile << elements.size() << " elements" << std::endl;
			hexfile.precision(16);
			for (size_t i = 0; i < (size_t)vertices.size(); i++)
			{
				hexfile << std::scientific << vertices[i] << std::endl;
			}
			int index[8] = { 1, 0, 2, 3, 5, 4, 6, 7 };
			for (size_t i = 0; i < (size_t)elements.size(); i++)
			{
				hexfile << "8";
				for (int j = 0; j < 8; j++)
				{
					hexfile << " " << elements[i]->indices[index[j]];
				}
				hexfile << std::endl;
			}
			hexfile.close();
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::save_reversed_hex(const char filename[])
	{
		std::ofstream hexfile(filename);
		if (hexfile.is_open())
		{
			hexfile << vertices.size() << " vertices" << std::endl;
			hexfile << elements.size() << " elements" << std::endl;
			hexfile.precision(16);
			for (size_t i = 0; i < (size_t)vertices.size(); i++)
			{
				hexfile << std::scientific << vertices[i] << std::endl;
			}
			int index[8] = { 1, 0, 2, 3, 5, 4, 6, 7 };
			for (size_t i = 0; i < (size_t)elements.size(); i++)
			{
				hexfile << "8";
				for (int j = 4; j < 8; j++)
				{
					hexfile << " " << elements[i]->indices[index[j]];
				}
				for (int j = 0; j < 4; j++)
				{
					hexfile << " " << elements[i]->indices[index[j]];
				}
				hexfile << std::endl;
			}
			hexfile.close();
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::save_vtk(const char filename[])
	{
		std::ofstream hexfile(filename);
		if (hexfile.is_open())
		{
			hexfile << "# vtk DataFile Version 3.0\n" << "Volume mesh\n" << "ASCII\n" << "DATASET UNSTRUCTURED_GRID\n"
				<< "POINTS " << vertices.size() << " double\n";
			hexfile.precision(16);
			for (size_t i = 0; i < (size_t)vertices.size(); i++)
			{
				hexfile << std::scientific << vertices[i] << "\n";
			}
			hexfile << "CELLS " << elements.size() << " " << elements.size() * 9 << "\n";
			int index[8] = { 4, 5, 1, 0, 7, 6, 2, 3 };
			for (size_t i = 0; i < (size_t)elements.size(); i++)
			{
				hexfile << "8";
				for (int j = 0; j < 8; j++)
				{
					hexfile << " " << elements[i]->indices[index[j]];
				}
				hexfile << "\n";
			}
			hexfile << "CELL_TYPES " << elements.size() << "\n";
			for (size_t i = 0; i < (size_t)elements.size(); i++)
			{
				hexfile << "12 \n";
			}
			hexfile.close();
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::load_data(const std::vector<TinyVector<double, 3>>& hex_vertices, const std::vector<size_t>& hex_indices, bool vtk_order)
	{
		clear_mem();
		vertices.assign(hex_vertices.begin(), hex_vertices.end());
		elements.resize(hex_indices.size() / 8);

		if (vtk_order)
		{
#pragma omp parallel for
			for (ptrdiff_t i = 0; i < (ptrdiff_t)hex_indices.size(); i += 8)
			{
				elements[i / 8] = new HexElement<Real>(hex_indices[i + 4], hex_indices[i + 5], hex_indices[i + 1], hex_indices[i],
					hex_indices[i + 7], hex_indices[i + 6], hex_indices[i + 2], hex_indices[i + 3]);
			}
		}
		else
		{
#pragma omp parallel for
			for (ptrdiff_t i = 0; i < (ptrdiff_t)hex_indices.size(); i += 8)
			{
				elements[i / 8] = new HexElement<Real>(hex_indices[i], hex_indices[i + 1], hex_indices[i + 2], hex_indices[i + 3],
					hex_indices[i + 4], hex_indices[i + 5], hex_indices[i + 6], hex_indices[i + 7]);
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::build_connectivity()
	{
		static const unsigned int face_id[][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 },
		{ 7, 6, 2, 3 }, { 5, 4, 0, 1 }, { 5, 1, 2, 6 }, { 4, 7, 3, 0 }
		};

		QuadFace qf;
		std::map<QuadFace, std::vector< std::pair<size_t, int> > > m_quadmap;
		for (size_t i = 0; i < elements.size(); i++)
		{
			elements[i]->id = i;
			for (int k = 0; k < 6; k++)
			{
				qf.set_vertices(elements[i]->indices[face_id[k][0]], elements[i]->indices[face_id[k][1]], elements[i]->indices[face_id[k][2]], elements[i]->indices[face_id[k][3]]);
				m_quadmap[qf].push_back(std::pair<size_t, int>(i, k));
				elements[i]->neighborHex[k].first = 0;
				elements[i]->neighborHex[k].second = 0;
			}
		}

		num_quad_faces = (size_t)m_quadmap.size();
		typename std::map<QuadFace, std::vector< std::pair<size_t, int> > >::iterator miter = m_quadmap.begin();
		for (; miter != m_quadmap.end(); miter++)
		{
			std::vector< std::pair<size_t, int> >& plist = miter->second;
			if (plist.size() == 1)
			{
				elements[plist[0].first]->neighborHex[plist[0].second] = std::pair<HexElement<Real>*, int>((HexElement<Real>*)(0), 0);
			}
			else if (plist.size() == 2) // size = 2
			{
				elements[plist[0].first]->neighborHex[plist[0].second] = std::pair<HexElement<Real>*, int>(elements[plist[1].first], plist[1].second);
				elements[plist[1].first]->neighborHex[plist[1].second] = std::pair<HexElement<Real>*, int>(elements[plist[0].first], plist[0].second);
			}
			else
			{
				assert(false);
			}
		}
		remove_unused_vertices(false);
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::remove_unused_vertices(bool removed)
	{
		std::vector<bool> vertices_tag(vertices.size(), false);
		for (size_t i = 0; i < elements.size(); i++)
		{
			for (int j = 0; j < 8; j++)
			{
				vertices_tag[elements[i]->indices[j]] = true;
			}
		}
		if (removed == false)
		{
			return;
		}

		std::vector<size_t> new_vert_id(vertices.size());
		size_t counter = 0;
		for (size_t i = 0; i < vertices_tag.size(); i++)
		{
			if (vertices_tag[i])
			{
				new_vert_id[i] = counter;
				counter++;
			}
		}
		if (counter != vertices_tag.size())
		{
			for (ptrdiff_t i = (ptrdiff_t)vertices.size() - 1; i >= 0; i--)
			{
				if (vertices_tag[i] == false)
				{
					vertices.erase(vertices.begin() + i);
				}
			}
			for (size_t i = 0; i < elements.size(); i++)
			{
				for (int j = 0; j < 8; j++)
				{
					elements[i]->indices[j] = new_vert_id[elements[i]->indices[j]];
				}
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::reorientation()
	{
		static const int face_id[][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 },
		{ 7, 6, 2, 3 }, { 5, 4, 0, 1 }, { 5, 1, 2, 6 }, { 4, 7, 3, 0 }
		};

		std::vector<bool> cell_tag(elements.size(), false);
		std::vector<bool> reversed_tag(elements.size(), false);
		bool has_reversed = false;

		for (size_t i = 0; i < elements.size(); i++)
		{
			elements[i]->id = i;
		}

		for (size_t i = 0; i < cell_tag.size(); i++)
		{
			if (cell_tag[i] == false)
			{
				std::queue<HexElement<Real>*> hqueue;
				hqueue.push(elements[i]);
				cell_tag[i] = true;
				while (!hqueue.empty())
				{
					HexElement<Real>* hc = hqueue.front();
					hqueue.pop();
					for (int j = 0; j < 6; j++)
					{
						HexElement<Real>* hadj = hc->neighborHex[j].first;
						if (hadj && cell_tag[hadj->id] == false)
						{
							cell_tag[hadj->id] = true;
							hqueue.push(hadj);
							QuadFace_O qh(hc->indices[face_id[j][0]], hc->indices[face_id[j][1]], hc->indices[face_id[j][2]], hc->indices[face_id[j][3]]);
							int k = hc->neighborHex[j].second;
							QuadFace_O qadj(hadj->indices[face_id[k][0]], hadj->indices[face_id[k][1]], hadj->indices[face_id[k][2]], hadj->indices[face_id[k][3]]);
							reversed_tag[hadj->id] = qh.quad_vert[1] == qadj.quad_vert[1] ? !reversed_tag[hc->id] : reversed_tag[hc->id];
							if (reversed_tag[hadj->id])
							{
								has_reversed = true;
							}
						}
					}
				}
			}
		}

		if (has_reversed)
		{
			for (size_t i = 0; i < elements.size(); i++)
			{
				if (reversed_tag[i])
				{
					for (int j = 0; j < 4; j++)
					{
						std::swap(elements[i]->indices[j], elements[i]->indices[j + 4]);
					}
				}
			}
			build_connectivity();
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	Mesh3D<Real>* VHexMesh<Real>::extract_boundary_quad_mesh(std::map<size_t, size_t>& vertexmap)
	{
		static const unsigned int face_id[][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 },
		{ 7, 6, 2, 3 }, { 5, 4, 0, 1 }, { 5, 1, 2, 6 }, { 4, 7, 3, 0 }
		};

		Mesh3D<Real>* bmesh = new Mesh3D<Real>;
		size_t count = 0;
		for (size_t i = 0; i < (size_t)elements.size(); i++)
		{
			for (int j = 0; j < 6; j++)
			{
				if (elements[i]->neighborHex[j].first != 0) { continue; }
				for (int k = 0; k < 4; k++)
				{
					size_t vertID = elements[i]->indices[face_id[j][k]];
					if (vertexmap.count(vertID) == 0)
					{
						vertexmap[vertID] = count;
						bmesh->insert_vertex(vertices[vertID]);
						count++;
					}
				}
			}
		}

		typename Mesh3D<Real>::VERTEX_LIST vlist(4);
		for (size_t i = 0; i < elements.size(); i++)
		{
			for (unsigned int j = 0; j < 6; j++)
			{
				if (elements[i]->neighborHex[j].first == 0)
				{
					vlist[0] = bmesh->get_vertex((int)vertexmap[elements[i]->indices[face_id[j][0]]]);
					vlist[1] = bmesh->get_vertex((int)vertexmap[elements[i]->indices[face_id[j][1]]]);
					vlist[2] = bmesh->get_vertex((int)vertexmap[elements[i]->indices[face_id[j][2]]]);
					vlist[3] = bmesh->get_vertex((int)vertexmap[elements[i]->indices[face_id[j][3]]]);
					bmesh->insert_face(vlist);
				}
			}
		}

		bmesh->update_mesh();
		return bmesh;
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	bool VHexMesh<Real>::detect_orientation()
	{
		int face_id[][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 },
{ 7, 6, 2, 3 }, { 5, 4, 0, 1 }, { 5, 1, 2, 6 }, { 4, 7, 3, 0 } };

		int pos_counter = 0, neg_counter = 0;
		for (size_t i = 0; i < elements.size(); i++)
		{
			Vector3d C = (vertices[elements[i]->indices[0]] + vertices[elements[i]->indices[1]] + vertices[elements[i]->indices[2]] + vertices[elements[i]->indices[3]]
				+ vertices[elements[i]->indices[4]] + vertices[elements[i]->indices[5]] + vertices[elements[i]->indices[6]] + vertices[elements[i]->indices[7]]) / 8;

			for (int j = 0; j < 6; j++)
			{
				Vector3d& V0 = vertices[elements[i]->indices[face_id[j][0]]];
				Vector3d& V1 = vertices[elements[i]->indices[face_id[j][1]]];
				Vector3d& V2 = vertices[elements[i]->indices[face_id[j][2]]];
				Vector3d& V3 = vertices[elements[i]->indices[face_id[j][3]]];

				double vol0 = (C - V2).Dot((C - V0).Cross(C - V1));
				double vol1 = (C - V3).Dot((C - V0).Cross(C - V2));
				if (vol0 > 0) pos_counter++; else neg_counter++;
				if (vol1 > 0) pos_counter++; else neg_counter++;
			}
		}
		return pos_counter > neg_counter;
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::load_ovm(const char filename[])
	{
		std::ifstream ovmfile(filename);

		std::set<MySortedTuple<2>> edgeset;

		if (ovmfile.is_open())
		{
			clear_mem();
			std::string str;
			size_t vsize = 0, esize = 0, fsize = 0, csize = 0;
			std::vector<size_t> edgevec, facevec;
			try
			{
				while (!ovmfile.eof())
				{
					ovmfile >> str;
					if (str == "Vertices")
					{
						ovmfile >> vsize;
						vertices.resize(vsize);
						for (size_t i = 0; i < vsize; i++)
						{
							ovmfile >> vertices[i];
						}
					}
					else if (str == "Edges")
					{
						ovmfile >> esize;
						edgevec.resize(2 * esize);
						for (size_t i = 0; i < edgevec.size(); i++)
							ovmfile >> edgevec[i];
						for (size_t i = 0; i < 2 * esize; i += 2)
							edgeset.insert(MySortedTuple<2>(edgevec[i], edgevec[i + 1]));
					}
					else if (str == "Faces")
					{
						ovmfile >> fsize;
						facevec.resize(4 * fsize);
						for (size_t i = 0; i < fsize; i++)
						{
							int npoly = 0;
							ovmfile >> npoly;
							if (npoly != 4)
								throw 1;
							for (int j = 0; j < 4; j++)
								ovmfile >> facevec[4 * i + j];
						}
					}
					else if (str == "Polyhedra")
					{
						ovmfile >> csize;
						elements.reserve(csize);
						int faceid[6];
						size_t ind[8], right[4];
						for (size_t i = 0; i < csize; i++)
						{
							int npoly = 0;
							ovmfile >> npoly;
							if (npoly != 6)
								throw 1;
							ovmfile >> faceid[0] >> faceid[1] >> faceid[2] >> faceid[3] >> faceid[4] >> faceid[5];

							size_t e0 = facevec[4 * (faceid[0] / 2)];
							size_t e1 = facevec[4 * (faceid[0] / 2) + 1];
							size_t e2 = facevec[4 * (faceid[0] / 2) + 2];
							size_t e3 = facevec[4 * (faceid[0] / 2) + 3];

							if (faceid[0] % 2 == 0)
							{
								ind[0] = edgevec[e0]; ind[1] = edgevec[e1]; ind[2] = edgevec[e2]; ind[3] = edgevec[e3];
							}
							else
							{
								ind[3] = edgevec[e0]; ind[2] = edgevec[e1]; ind[1] = edgevec[e2]; ind[0] = edgevec[e3];
							}
							bool find = false;
							for (int j = 1; j < 6; j++)
							{
								size_t ne0 = facevec[4 * (faceid[j] / 2)];
								if (ne0 / 2 == e0 / 2 || ne0 / 2 == e1 / 2 || ne0 / 2 == e2 / 2 || ne0 / 2 == e3 / 2) continue;
								size_t ne1 = facevec[4 * (faceid[j] / 2) + 1];
								if (ne1 / 2 == e0 / 2 || ne1 / 2 == e1 / 2 || ne1 / 2 == e2 / 2 || ne1 / 2 == e3 / 2) continue;
								size_t ne2 = facevec[4 * (faceid[j] / 2) + 2];
								if (ne2 / 2 == e0 / 2 || ne2 / 2 == e1 / 2 || ne2 / 2 == e2 / 2 || ne2 / 2 == e3 / 2) continue;
								size_t ne3 = facevec[4 * (faceid[j] / 2) + 3];
								if (ne3 / 2 == e0 / 2 || ne3 / 2 == e1 / 2 || ne3 / 2 == e2 / 2 || ne3 / 2 == e3 / 2) continue;

								if (faceid[j] % 2 == 0)
								{
									right[0] = edgevec[ne0]; right[1] = edgevec[ne1]; right[2] = edgevec[ne2]; right[3] = edgevec[ne3];
								}
								else
								{
									right[3] = edgevec[ne0]; right[2] = edgevec[ne1]; right[1] = edgevec[ne2]; right[0] = edgevec[ne3];
								}

								for (int k = 0; k < 4; k++)
								{
									MySortedTuple<2> e(ind[0], right[k]);
									if (edgeset.find(e) != edgeset.end())
									{
										ind[4] = right[k]; ind[7] = right[(k + 1) % 4]; ind[6] = right[(k + 2) % 4]; ind[5] = right[(k + 3) % 4];

										elements.push_back(new HexElement<Real>(ind));
										find = true;
										break;
									}
								}
								if (find)
									break;
							}
						}

						break;
					}
				}
			}
			catch (...)
			{
				std::cout << "catch error!" << std::endl;
				vertices.resize(0);
				elements.resize(0);
				ovmfile.close();
				return;
			}

			ovmfile.close();
			remove_unused_vertices();
			update_boundingbox();
			build_connectivity();
			reorientation();
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::save_ovm(const char filename[])
	{
		std::ofstream ovmfile(filename);
		if (!ovmfile.is_open()) return;

		std::map<MySortedTuple<2>, size_t> edges;
		std::map<QuadFace, size_t> faces;

		static const unsigned int face_id[][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 },{ 7, 6, 2, 3 }, { 5, 4, 0, 1 }, { 5, 1, 2, 6 }, { 4, 7, 3, 0 } };
		static const unsigned int edge_id[][2] = { {0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7} };
		for (size_t i = 0; i < elements.size(); i++)
		{
			for (int j = 0; j < 12; j++)
			{
				edges[MySortedTuple<2>(elements[i]->indices[edge_id[j][0]], elements[i]->indices[edge_id[j][1]])] = 0;
			}
			for (int j = 0; j < 6; j++)
			{
				faces[QuadFace(elements[i]->indices[face_id[j][0]], elements[i]->indices[face_id[j][1]], elements[i]->indices[face_id[j][2]], elements[i]->indices[face_id[j][3]])] = 0;
			}
		}

		ovmfile << "OVM ASCII" << std::endl;
		ovmfile << "Vertices" << std::endl << vertices.size() << std::endl;
		for (size_t i = 0; i < vertices.size(); i++)
		{
			ovmfile << vertices[i] << std::endl;
		}
		ovmfile << "Edges" << std::endl << edges.size() << std::endl;
		size_t ecounter = 0;
		for (typename std::map<MySortedTuple<2>, size_t>::iterator eiter = edges.begin(); eiter != edges.end(); eiter++)
		{
			ovmfile << eiter->first.sorted_vert[0] << ' ' << eiter->first.sorted_vert[1] << std::endl;
			eiter->second = ecounter;
			ecounter += 2;
		}
		ovmfile << "Faces" << std::endl << faces.size() << std::endl;
		size_t fcounter = 0;
		for (typename std::map<QuadFace, size_t>::iterator qiter = faces.begin(); qiter != faces.end(); qiter++)
		{
			ovmfile << "4";
			for (int i = 0; i < 4; i++)
			{
				MyTuple<2> e(qiter->first.quad_vert[i], qiter->first.quad_vert[(i + 1) % 4]);
				MySortedTuple<2> se(qiter->first.quad_vert[i], qiter->first.quad_vert[(i + 1) % 4]);
				size_t eindex = edges[se];
				ovmfile << ' ' << (e.vert[0] == se.sorted_vert[0] ? eindex : eindex + 1);
			}
			ovmfile << std::endl;
			qiter->second = fcounter;
			fcounter += 2;
		}
		ovmfile << "Polyhedra" << std::endl << elements.size() << std::endl;
		for (size_t i = 0; i < elements.size(); i++)
		{
			ovmfile << "6";
			for (int j = 0; j < 6; j++)
			{
				QuadFace q(elements[i]->indices[face_id[j][0]], elements[i]->indices[face_id[j][1]], elements[i]->indices[face_id[j][2]], elements[i]->indices[face_id[j][3]]);
				QuadFace_O qo(elements[i]->indices[face_id[j][0]], elements[i]->indices[face_id[j][1]], elements[i]->indices[face_id[j][2]], elements[i]->indices[face_id[j][3]]);
				size_t eindex = faces[q];
				ovmfile << ' ' << (q.quad_vert[1] == qo.quad_vert[1] ? eindex + 1 : eindex);
			}
			ovmfile << std::endl;
		}

		ovmfile.close();
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::load_mesh(const char filename[])
	{
		std::ifstream meshfile(filename);
		if (meshfile.is_open())
		{
			clear_mem();
			std::string str;
			int tag;
			try
			{
				size_t vsize = 0, fsize = 0, csize = 0;
				while (!meshfile.eof())
				{
					meshfile >> str;
					if (str == "Vertices")
					{
						meshfile >> vsize;
						vertices.resize(vsize);
						for (size_t i = 0; i < vsize; i++)
						{
							meshfile >> vertices[i] >> tag;
						}
					}
					else if (str == "Quadrilaterals")
					{
						meshfile >> fsize;
						for (size_t i = 0; i < fsize; i++)
							meshfile >> tag >> tag >> tag >> tag >> tag;
					}
					else if (str == "Hexahedra")
					{
						meshfile >> csize;
						elements.reserve(csize);

						size_t ind[8];
						for (size_t i = 0; i < csize; i++)
						{
							meshfile >> ind[0] >> ind[1] >> ind[2] >> ind[3] >> ind[4] >> ind[5] >> ind[6] >> ind[7] >> tag;
							for (int j = 0; j < 8; j++) ind[j]--;
							elements.push_back(new HexElement<Real>(ind));
						}
						break;
					}
					else if (str == "End")
					{
						break;
					}
				}
			}
			catch (...)
			{
				std::cout << "catch error!" << std::endl;
				vertices.resize(0);
				elements.resize(0);
				meshfile.close();
				return;
			}

			meshfile.close();
			remove_unused_vertices();
			update_boundingbox();
			build_connectivity();
			reorientation();
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::save_mesh(const char filename[])
	{
		std::ofstream meshfile(filename);
		if (meshfile.is_open())
		{
			meshfile << "MeshVersionFormatted 1" << std::endl;
			//meshfile << "Dimension 3" << std::endl;
			meshfile << "Dimension" << std::endl << "3" << std::endl;

			//meshfile << "Vertices " << vertices.size() << std::endl;
			meshfile << "Vertices" << std::endl << vertices.size() << std::endl;
			for (size_t i = 0; i < vertices.size(); i++)
				meshfile << vertices[i] << " 0" << std::endl;

			static const unsigned int face_id[][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 },
{ 7, 6, 2, 3 }, { 5, 4, 0, 1 }, { 5, 1, 2, 6 }, { 4, 7, 3, 0 }
			};

			Mesh3D<Real>* bmesh = new Mesh3D<Real>;
			size_t count = 0;
			for (size_t i = 0; i < (size_t)elements.size(); i++)
			{
				for (int j = 0; j < 6; j++)
				{
					if (elements[i]->neighborHex[j].first == 0)
						count++;
				}
			}
			//meshfile << "Quadrilaterals " << count << std::endl;
			meshfile << "Quadrilaterals" << std::endl << count << std::endl;

			for (size_t i = 0; i < elements.size(); i++)
			{
				for (unsigned int j = 0; j < 6; j++)
				{
					if (elements[i]->neighborHex[j].first == 0)
					{
						meshfile
							<< elements[i]->indices[face_id[j][0]] + 1 << ' '
							<< elements[i]->indices[face_id[j][1]] + 1 << ' '
							<< elements[i]->indices[face_id[j][2]] + 1 << ' '
							<< elements[i]->indices[face_id[j][3]] + 1 << " 0" << std::endl;
					}
				}
			}

			//meshfile << "Hexahedra " << elements.size() << std::endl;
			meshfile << "Hexahedra" << std::endl << elements.size() << std::endl;
			for (size_t i = 0; i < elements.size(); i++)
			{
				for (int j = 0; j < 8; j++)
					meshfile << elements[i]->indices[j] + 1 << ' ';
				meshfile << "0" << std::endl;
			}
			meshfile << "End" << std::endl;
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::merge_vertices(double DIST_THRES)
	{
		if (vertices.size() < 8) return;

		ANNpointArray dataPts;
		int n_pt = (int)vertices.size();
		dataPts = annAllocPts(n_pt, 3);
		for (int k = 0; k < n_pt; k++)
		{
			dataPts[k][0] = vertices[k][0];
			dataPts[k][1] = vertices[k][1];
			dataPts[k][2] = vertices[k][2];
		}
		ANNkd_tree* data_kdTree = new ANNkd_tree(dataPts, n_pt, 3);
		std::vector<ptrdiff_t> mapid(vertices.size(), -1);

		std::vector<ptrdiff_t> merge2uniqueID_map(vertices.size(), -1), back2overlapID_map(vertices.size(), -1);

		std::vector<int> nnIdx; // allocate near neigh indices
		std::vector<double> dists; // allocate near neighbor dists
		size_t unique_pt_counter = 0;

		std::vector<TinyVector<Real, 3>> new_vertices;
		new_vertices.reserve(vertices.size());

		for (size_t i = 0; i < vertices.size(); i++)
		{
			if (merge2uniqueID_map[i] >= 0) { continue; }
			bool done = true;
			int K = 4;
			while (done)
			{
				nnIdx.resize(K);
				dists.resize(K);
				data_kdTree->annkSearch(			// search
					&vertices[i][0],				// query point
					K, 								// number of near neighbors
					&nnIdx[0], 						// nearest neighbors (returned)
					&dists[0], 						// distance (returned)
					0.0);
				if (dists[K - 1] < DIST_THRES)
				{
					//if all the found K points are within distance threshold, there might be more overlap points.
					K = 2 * K;
					continue;
				}
				while (dists[K - 1] > DIST_THRES)
				{
					K--;
				}
				//now there are K overlapped points
				for (int i = 0; i < K; i++)
				{
					merge2uniqueID_map[nnIdx[i]] = unique_pt_counter;
				}
				break;
			}
			new_vertices.push_back(vertices[i]);
			back2overlapID_map[unique_pt_counter] = i;
			unique_pt_counter++;
		}
		delete data_kdTree;
		annDeallocPts(dataPts);
		annClose();

		if (vertices.size() > new_vertices.size())
		{
			vertices.assign(new_vertices.begin(), new_vertices.end());

			for (size_t i = 0; i < elements.size(); i++)
			{
				for (int j = 0; j < 8; j++)
				{
					elements[i]->indices[j] = merge2uniqueID_map[elements[i]->indices[j]];
				}
			}
			build_connectivity();
			reorientation();
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::merge_hexes(VHexMesh<Real>* hexmesh, double DIST_THRES)
	{
		if (hexmesh == 0) return;

		ANNpointArray dataPts;
		int n_pt = (int)vertices.size();
		dataPts = annAllocPts(n_pt, 3);
		for (int k = 0; k < n_pt; k++)
		{
			dataPts[k][0] = vertices[k][0];
			dataPts[k][1] = vertices[k][1];
			dataPts[k][2] = vertices[k][2];
		}
		ANNkd_tree* data_kdTree = new ANNkd_tree(dataPts, n_pt, 3);
		std::vector<ptrdiff_t> mapid(hexmesh->vertices.size(), -1);
		for (size_t i = 0; i < hexmesh->vertices.size(); i++)
		{
			int nnIdx;
			double dists;

			data_kdTree->annkSearch(			// search
				&hexmesh->vertices[i][0],		// query point
				1, 								// number of near neighbors
				&nnIdx, 						// nearest neighbors (returned)
				&dists, 						// distance (returned)
				0.0);
			if (dists < DIST_THRES)
			{
				mapid[i] = nnIdx;
			}
			else
			{
				//mapid[i] = vertices.size();
				//vertices.push_back(hexmesh->vertices[i]);
			}
		}

		delete data_kdTree;
		annDeallocPts(dataPts);
		annClose();

		std::set<MySortedTuple<8>> cells;
		for (size_t i = 0; i < elements.size(); i++)
		{
			cells.insert(MySortedTuple<8>(elements[i]->indices[0], elements[i]->indices[1], elements[i]->indices[2], elements[i]->indices[3],
				elements[i]->indices[4], elements[i]->indices[5], elements[i]->indices[6], elements[i]->indices[7]));
		}

		for (size_t i = 0; i < hexmesh->elements.size(); i++)
		{
			HexElement<Real>* he = hexmesh->elements[i];

			bool invalid = false;
			for (int j = 0; j < 8; j++)
			{
				if (mapid[he->indices[j]] == -1)
				{
					invalid = true; break;
				}
			}
			if (invalid) continue;

			MySortedTuple<8> cell(mapid[he->indices[0]], mapid[he->indices[1]], mapid[he->indices[2]], mapid[he->indices[3]],
				mapid[he->indices[4]], mapid[he->indices[5]], mapid[he->indices[6]], mapid[he->indices[7]]);

			if (cells.find(cell) == cells.end())
				elements.push_back(new HexElement<Real>(mapid[he->indices[0]], mapid[he->indices[1]], mapid[he->indices[2]], mapid[he->indices[3]],
					mapid[he->indices[4]], mapid[he->indices[5]], mapid[he->indices[6]], mapid[he->indices[7]]));
		}
		update_boundingbox();
		build_connectivity();
		reorientation();
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::subdivide()
	{
		static const int face_id[][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 },
		{ 7, 6, 2, 3 }, { 5, 4, 0, 1 }, { 5, 1, 2, 6 }, { 4, 7, 3, 0 } };

		static const int edge_id[][2] = { {0,3},{3,2},{2,1},{1,0},{4,5},{5,6},{6,7},{7,4},
		{0,4},{3,7},{2,6},{1,5} };

		std::vector< HexElement<Real>* > sub_elements;
		vertices.reserve(vertices.size() * 3);
		sub_elements.reserve(elements.size() * 8);

		std::map< std::pair<size_t, size_t>, size_t> edge_vertices_map;
		std::map< QuadFace, size_t> face_vertices_map;
		size_t center_id, face_vid[6], edge_vid[12];
		for (size_t i = 0; i < elements.size(); i++)
		{
			TinyVector<Real, 3> Z(0, 0, 0);
			for (int j = 0; j < 8; j++)
			{
				Z += vertices[elements[i]->indices[j]];
			}
			Z /= 8;
			center_id = vertices.size();
			vertices.push_back(Z);

			//face vertices
			for (int j = 0; j < 6; j++)
			{
				QuadFace qf(elements[i]->indices[face_id[j][0]], elements[i]->indices[face_id[j][1]], elements[i]->indices[face_id[j][2]], elements[i]->indices[face_id[j][3]]);
				std::map< QuadFace, size_t>::iterator iter = face_vertices_map.find(qf);
				if (iter != face_vertices_map.end())
				{
					face_vid[j] = iter->second;
				}
				else
				{
					Z = (vertices[qf.quad_vert[0]] + vertices[qf.quad_vert[1]] + vertices[qf.quad_vert[2]] + vertices[qf.quad_vert[3]]) / 4;
					face_vid[j] = vertices.size();
					vertices.push_back(Z);
					face_vertices_map[qf] = face_vid[j];
				}
			}

			//edge vertices
			for (int j = 0; j < 12; j++)
			{
				size_t i1 = std::min(elements[i]->indices[edge_id[j][0]], elements[i]->indices[edge_id[j][1]]);
				size_t i2 = std::max(elements[i]->indices[edge_id[j][0]], elements[i]->indices[edge_id[j][1]]);
				std::pair<size_t, size_t> qf(i1, i2);
				std::map< std::pair<size_t, size_t>, size_t>::iterator iter = edge_vertices_map.find(qf);
				if (iter != edge_vertices_map.end())
				{
					edge_vid[j] = iter->second;
				}
				else
				{
					Z = (vertices[qf.first] + vertices[qf.second]) / 2;
					edge_vid[j] = vertices.size();
					vertices.push_back(Z);
					edge_vertices_map[qf] = edge_vid[j];
				}
			}

			HexElement<Real>* H0 = new HexElement<Real>(face_vid[3], edge_vid[11], face_vid[4], center_id, edge_vid[4], elements[i]->indices[5], edge_vid[5], face_vid[1]);
			HexElement<Real>* H1 = new HexElement<Real>(edge_vid[3], elements[i]->indices[1], edge_vid[2], face_vid[0], face_vid[3], edge_vid[11], face_vid[4], center_id);
			HexElement<Real>* H2 = new HexElement<Real>(edge_vid[8], face_vid[3], center_id, face_vid[5], elements[i]->indices[4], edge_vid[4], face_vid[1], edge_vid[7]);
			HexElement<Real>* H3 = new HexElement<Real>(elements[i]->indices[0], edge_vid[3], face_vid[0], edge_vid[0], edge_vid[8], face_vid[3], center_id, face_vid[5]);
			HexElement<Real>* H4 = new HexElement<Real>(center_id, face_vid[4], edge_vid[10], face_vid[2], face_vid[1], edge_vid[5], elements[i]->indices[6], edge_vid[6]);
			HexElement<Real>* H5 = new HexElement<Real>(face_vid[0], edge_vid[2], elements[i]->indices[2], edge_vid[1], center_id, face_vid[4], edge_vid[10], face_vid[2]);
			HexElement<Real>* H6 = new HexElement<Real>(face_vid[5], center_id, face_vid[2], edge_vid[9], edge_vid[7], face_vid[1], edge_vid[6], elements[i]->indices[7]);
			HexElement<Real>* H7 = new HexElement<Real>(edge_vid[0], face_vid[0], edge_vid[1], elements[i]->indices[3], face_vid[5], center_id, face_vid[2], edge_vid[9]);
			sub_elements.push_back(H0); sub_elements.push_back(H1); sub_elements.push_back(H2); sub_elements.push_back(H3);
			sub_elements.push_back(H4); sub_elements.push_back(H5); sub_elements.push_back(H6); sub_elements.push_back(H7);
		}

		for (size_t i = 0; i < elements.size(); i++)
		{
			delete elements[i];
		}
		vertices.shrink_to_fit();
		elements.resize(sub_elements.size());
		std::copy(sub_elements.begin(), sub_elements.end(), elements.begin());
		build_connectivity();
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void VHexMesh<Real>::padding()
	{
		std::map<size_t, size_t> vertexmap;
		Mesh3D<Real>* bndmesh = extract_boundary_quad_mesh(vertexmap);

		size_t counter = vertices.size(), orginal_size = vertices.size();
		std::vector<size_t> reverse_ids(vertexmap.size()), map_ids(vertexmap.size());

		vertices.resize(counter + vertexmap.size());
		for (typename std::map<size_t, size_t>::iterator iter = vertexmap.begin(); iter != vertexmap.end(); iter++)
		{
			reverse_ids[iter->second] = iter->first;
			map_ids[iter->second] = counter;
			HE_vert<Real>* hv = bndmesh->get_vertex(iter->second);
			HE_edge<Real>* he = hv->edge;
			Real ave_length = 0;
			do
			{
				ave_length += (he->vert->pos - he->pair->vert->pos).SquaredLength();
				he = he->pair->next;
			} while (he != hv->edge);
			ave_length = sqrt(ave_length / hv->degree);
			vertices[counter] = vertices[iter->first];
			vertices[iter->first] -= (Real)0.1 * ave_length * bndmesh->get_vertex(iter->second)->normal;
			counter++;
		}

		counter = elements.size();
		elements.resize(counter + bndmesh->get_num_of_faces());
		for (ptrdiff_t i = 0; i < bndmesh->get_num_of_faces(); i++)
		{
			HE_face<Real>* hf = bndmesh->get_face(i);
			HE_edge<Real>* he = hf->edge;
			elements[counter + i] = new HexElement<Real>(
				reverse_ids[he->next->next->vert->id], reverse_ids[he->next->vert->id], reverse_ids[he->vert->id], reverse_ids[he->pair->vert->id],
				map_ids[he->next->next->vert->id], map_ids[he->next->vert->id], map_ids[he->vert->id], map_ids[he->pair->vert->id]
				);
			//elements[counter + i] = new HexElement<Real>(
			//	reverse_ids[he->pair->vert->id], reverse_ids[he->vert->id], reverse_ids[he->next->vert->id], reverse_ids[he->next->next->vert->id],
			//	map_ids[he->pair->vert->id], map_ids[he->vert->id], map_ids[he->next->vert->id], map_ids[he->next->next->vert->id]
			//	);
		}
		build_connectivity();
		reorientation();
		delete bndmesh;
	}
	//////////////////////////////////////////////////////////////////////////
	template class VHexMesh<double>;
}