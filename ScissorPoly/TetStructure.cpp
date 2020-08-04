#include "TetStructure.h"
#include <fstream>
#include <map>
#include <string>

namespace ig
{
	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	TetStructure<Real>::TetStructure()
		:m_kdtree(nullptr)
	{
	}
	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	TetStructure<Real>::~TetStructure()
	{
		clear();
	}
	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	void TetStructure<Real>::clear()
	{
		for (typename std::vector<Tetrahedron<Real>*>::iterator liter =
			tetras.begin(); liter != tetras.end(); liter++)
		{
			delete *liter;
		}

		for (typename std::vector<TetVertex<Real>*>::iterator viter =
			tetra_vertices.begin(); viter != tetra_vertices.end(); viter++)
		{
			delete *viter;
		}

		tetras.clear();
		tetra_vertices.clear();
		if (m_kdtree)
		{
			delete m_kdtree;
			m_kdtree = 0;
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	void TetStructure<Real>::build_connectivity()
	{
		//collect adjacent information

		std::map<MyIntTriple<Real>, std::pair<Tetrahedron<Real>*, Tetrahedron<Real>*> > facemap;
		static unsigned int sortid[4][3] = { { 0, 2, 1 }, { 2, 3, 1 }, { 1, 3, 0 }, { 0, 3, 2 } };

		for (size_t i = 0; i < tetras.size(); i++)
		{
			for (unsigned int j = 0; j < 4; j++)
			{
				MyIntTriple<Real> mt(tetras[i]->vertex[sortid[j][0]], tetras[i]->vertex[sortid[j][1]], tetras[i]->vertex[sortid[j][2]]);

				typename std::map<MyIntTriple<Real>, std::pair<	Tetrahedron<Real>*, Tetrahedron<Real>*> >::iterator
					fiter = facemap.find(mt);

				if (fiter != facemap.end())
				{
					fiter->second.second = tetras[i];
				}
				else
				{
					facemap[mt] = std::pair<Tetrahedron<Real>*, Tetrahedron<
						Real>*>(tetras[i], (Tetrahedron<Real>*) (0));
				}
			}
		}

		for (typename std::map<MyIntTriple<Real>, std::pair<Tetrahedron<Real>*, Tetrahedron<Real>*> >::iterator miter = facemap.begin();
			miter != facemap.end(); miter++)
		{
			Tetrahedron<Real>* i1 = miter->second.first;
			Tetrahedron<Real>* i2 = miter->second.second;

			if (i1 && i2)
			{
				i1->set_neighbor_tet(i2);
				i2->set_neighbor_tet(i1);
			}
			else if (i1 == 0 || i2 == 0)
			{
				miter->first.id[0]->boundary = true;
				miter->first.id[1]->boundary = true;
				miter->first.id[2]->boundary = true;
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	void TetStructure<Real>::load_tet(const std::vector<CVec<Real, 3> >& points,
		const std::vector<unsigned int>& indices)
	{
		clear();
		size_t n_tetra_vertices = points.size();
		size_t n_tetras = indices.size() / 4;

		tetras.reserve(n_tetras);
		tetra_vertices.reserve(n_tetra_vertices);

		for (size_t i = 0; i < n_tetra_vertices; i++)
		{
			TetVertex<Real>* new_v = new TetVertex<Real>(points[i]);
			new_v->id = i;
			tetra_vertices.push_back(new_v);
		}

		for (size_t i = 0; i < n_tetras; i++)
		{
			Tetrahedron<Real>* new_tet = new Tetrahedron<Real>(
				tetra_vertices[indices[4 * i]], tetra_vertices[indices[4 * i
				+ 1]], tetra_vertices[indices[4 * i + 2]],
				tetra_vertices[indices[4 * i + 3]]);
			new_tet->id = i;
			tetras.push_back(new_tet);
		}

		build_connectivity();
	}

	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	void TetStructure<Real>::load_tet(const char filename[])
	{
		clear();

		std::ifstream mfile(filename);

		if (!mfile.is_open())
		{
			return;
		}

		std::string str1, str2;
		unsigned int n_tetra_vertices, n_tetras;
		mfile >> n_tetra_vertices >> str1 >> n_tetras >> str2;

		tetras.reserve(n_tetras);
		tetra_vertices.reserve(n_tetra_vertices);

		for (size_t i = 0; i < n_tetra_vertices; i++)
		{
			TetVertex<Real>* new_v = new TetVertex<Real>();
			mfile >> new_v->pos[0] >> new_v->pos[1] >> new_v->pos[2];
			new_v->id = i;
			tetra_vertices.push_back(new_v);
		}

		unsigned int index[4];

		unsigned int four;

		for (size_t i = 0; i < n_tetras; i++)
		{
			mfile >> four;

			for (unsigned int j = 0; j < four; j++)
			{
				mfile >> index[j];
			}

			if (four != 4)
			{
				continue;
			}

			Tetrahedron<Real>* new_tet = new Tetrahedron<Real>(
				tetra_vertices[index[0]], tetra_vertices[index[1]],
				tetra_vertices[index[2]], tetra_vertices[index[3]]);
			new_tet->id = i;
			tetras.push_back(new_tet);
		}

		mfile.close();
		build_connectivity();
	}
	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	void TetStructure<Real>::save_tet(const char filename[])
	{
		std::ofstream mfile(filename);

		if (mfile.is_open())
		{
			mfile << tetra_vertices.size() << " vertices\n" << tetras.size() << " cells\n";
			//mfile.precision(16);
			std::map<TetVertex<Real>*, size_t> mapvid;
			//size_t count = 0;

			for (size_t i = 0; i < tetra_vertices.size(); i++)
			{
				mfile << std::scientific << tetra_vertices[i]->pos << "\n";
			}

			for (size_t i = 0; i < tetras.size(); i++)
			{
				mfile << "4 " << tetras[i]->vertex[0]->id << ' '
					<< tetras[i]->vertex[1]->id << ' '
					<< tetras[i]->vertex[2]->id << ' '
					<< tetras[i]->vertex[3]->id << std::endl;
			}

			mfile.close();
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	void TetStructure<Real>::load_vtk(const char filename[])
	{
		clear();

		std::ifstream vtkfile(filename);

		if (!vtkfile.is_open())
		{
			return;
		}

		std::string str;
		try
		{
			do
			{
				vtkfile >> str;
			} while (str != "DATASET");
			vtkfile >> str;
			if (str != "UNSTRUCTURED_GRID")
			{
				throw 1;
			}
			do
			{
				vtkfile >> str;
			} while (str != "POINTS");

			unsigned int n_tetra_vertices, n_tetras;

			vtkfile >> n_tetra_vertices >> str;

			tetra_vertices.reserve(n_tetra_vertices);

			for (size_t i = 0; i < n_tetra_vertices; i++)
			{
				TetVertex<Real>* new_v = new TetVertex<Real>();
				vtkfile >> new_v->pos[0] >> new_v->pos[1] >> new_v->pos[2];
				new_v->id = i;
				tetra_vertices.push_back(new_v);
			}
			do
			{
				vtkfile >> str;
			} while (str != "CELLS");

			unsigned int num_ne_type;

			vtkfile >> n_tetras >> num_ne_type;

			tetras.reserve(n_tetras);

			unsigned int index[4];

			unsigned int four;

			for (size_t i = 0; i < n_tetras; i++)
			{
				vtkfile >> four;

				for (unsigned int j = 0; j < four; j++)
				{
					vtkfile >> index[j];
				}

				if (four != 4)
				{
					continue;
				}

				Tetrahedron<Real>* new_tet = new Tetrahedron<Real>(
					tetra_vertices[index[0]], tetra_vertices[index[1]],
					tetra_vertices[index[2]], tetra_vertices[index[3]]);
				new_tet->id = i;
				tetras.push_back(new_tet);
			}

			vtkfile.close();
			build_connectivity();
		}
		catch (...)
		{
			//...
			vtkfile.close();
		}
	}
	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	void TetStructure<Real>::save_vtk(const char filename[])
	{
		std::ofstream vtkfile(filename);
		if (!vtkfile.is_open()) return;
		vtkfile.precision(16);
		vtkfile << "# vtk DataFile Version 2.0" << std::endl
			<< "Tetrahedral Mesh" << std::endl
			<< "ASCII" << std::endl
			<< "DATASET UNSTRUCTURED_GRID\n"
			<< "POINTS " << tetra_vertices.size() << " double\n";
		for (size_t i = 0; i < tetra_vertices.size(); i++)
		{
			vtkfile << std::scientific << tetra_vertices[i]->pos << std::endl;
		}
		vtkfile << "CELLS " << tetras.size() << " " << tetras.size() * 5 << "\n";
		for (size_t i = 0; i < tetras.size(); i++)
		{
			vtkfile << "4 " << tetras[i]->vertex[0]->id << " "
				<< tetras[i]->vertex[1]->id << " "
				<< tetras[i]->vertex[2]->id << " "
				<< tetras[i]->vertex[3]->id << std::endl;
		}
		vtkfile << "CELL_TYPES " << tetras.size() << std::endl;
		for (size_t i = 0; i < tetras.size(); i++)
			vtkfile << 10 << std::endl;
        vtkfile.close();
	}
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void TetStructure<Real>::find_adjacent_tets_by_edge(Tetrahedron<Real> * tet, const unsigned int vi, const unsigned int vj, std::vector<Tetrahedron<Real>*> & adjacent_tets)
	{
		static unsigned int table[][4] =
		{
			{0, 1, 2, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {1, 2, 0, 3}, {1, 3, 2, 0}, {2, 3, 0, 1},
			{1, 0, 3, 2}, {2, 0, 1, 3}, {3, 0, 2, 1}, {2, 1, 3, 0}, {3, 1, 0, 2}, {3, 2, 1, 0}
		};
		//find other two vertices' ids
		unsigned int choice = 0;
		for (int i = 0; i < 12; i++)
		{
			if (vi == table[i][0] && vj == table[i][1])
			{
				choice = i;
				break;
			}
		}

		adjacent_tets.resize(0);
		adjacent_tets.push_back(tet);

		std::vector<TetVertex<Real>*> loop;
		loop.resize(0);
		loop.push_back(tet->vertex[table[choice][2]]);
		loop.push_back(tet->vertex[table[choice][3]]);

		Tetrahedron<Real>* prev_tet = tet;
		Tetrahedron<Real>* start_tet = tet;
		TetVertex<Real>* start_v = 0;
		for (size_t j = 0; j < (size_t)loop.size(); j++)
		{
			unsigned int local_id = start_tet->get_local_vertex_id(loop[j]);
			start_v = start_tet->get_opposite_vertex(local_id);
			start_tet = start_tet->neighborTet[local_id];
			if (start_tet == 0 || start_tet == tet)
			{
				break;
			}
			adjacent_tets.push_back(start_tet);
			loop.push_back(start_v);
			prev_tet = start_tet;
		}
		if (start_tet == 0) // the edge is a boundary edge, we need to turn around to find other tets
		{
			TetVertex<Real>* a = loop.back();
			TetVertex<Real>* b = loop[(unsigned int)loop.size() - 2];
			loop.resize(0);
			loop.push_back(a);
			loop.push_back(b);
			start_tet = prev_tet;
			adjacent_tets.resize(0);
			adjacent_tets.push_back(start_tet);
			for (size_t j = 0; j < (size_t)loop.size(); j++)
			{
				unsigned int local_id = start_tet->get_local_vertex_id(loop[j]);
				start_v = start_tet->get_opposite_vertex(local_id);
				start_tet = start_tet->neighborTet[local_id];
				if (start_tet == 0)
				{
					break;
				}
				adjacent_tets.push_back(start_tet);
				loop.push_back(start_v);
			}
		}
	}
	
	//////////////////////////////////////////////////////////////////////////
	template <typename Real>
	void TetStructure<Real>::find_adjacent_tets_by_edge(Tetrahedron<Real> * tet, const unsigned int vi, const unsigned int vj,
		std::vector<Tetrahedron<Real>*> & adjacent_tets, std::vector<TetVertex<Real>*> & loop)
	{
		static unsigned int table[][4] =
		{
			{0, 1, 2, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {1, 2, 0, 3}, {1, 3, 2, 0}, {2, 3, 0, 1},
			{1, 0, 3, 2}, {2, 0, 1, 3}, {3, 0, 2, 1}, {2, 1, 3, 0}, {3, 1, 0, 2}, {3, 2, 1, 0}
		};
		//find other two vertices' ids
		unsigned int choice = 0;
		for (int i = 0; i < 12; i++)
		{
			if (vi == table[i][0] && vj == table[i][1])
			{
				choice = i;
				break;
			}
		}

		adjacent_tets.resize(0);
		adjacent_tets.push_back(tet);

		loop.resize(0);
		loop.push_back(tet->vertex[table[choice][2]]);
		loop.push_back(tet->vertex[table[choice][3]]);

		Tetrahedron<Real>* prev_tet = tet;
		Tetrahedron<Real>* start_tet = tet;
		TetVertex<Real>* start_v = 0;
		for (size_t j = 0; j < (size_t)loop.size(); j++)
		{
			unsigned int local_id = start_tet->get_local_vertex_id(loop[j]);
			start_v = start_tet->get_opposite_vertex(local_id);
			if (start_v == 0 && start_tet->neighborTet[local_id])
			{
				assert(false);
			}
			start_tet = start_tet->neighborTet[local_id];

			if (start_tet == 0 || start_tet == tet)
			{
				break;
			}
			assert(start_tet!=0);
			adjacent_tets.push_back(start_tet);
			loop.push_back(start_v);
			prev_tet = start_tet;
		}
		if (start_tet == 0) // the edge is a boundary edge, we need to turn around to find other tets
		{
			TetVertex<Real>* a = loop.back();
			TetVertex<Real>* b = loop[(unsigned int)loop.size() - 2];
			loop.resize(0);
			loop.push_back(a);
			loop.push_back(b);
			start_tet = prev_tet;
			adjacent_tets.resize(0);
			adjacent_tets.push_back(start_tet);
			for (size_t j = 0; j < (size_t)loop.size(); j++)
			{
				unsigned int local_id = start_tet->get_local_vertex_id(loop[j]);
				start_v = start_tet->get_opposite_vertex(local_id);
				start_tet = start_tet->neighborTet[local_id];
				if (start_tet == 0)
				{
					break;
				}
				assert(start_tet!=0);
				adjacent_tets.push_back(start_tet);
				loop.push_back(start_v);
			}
		}
	}


	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	void TetStructure<Real>::build_kdtree()
	{
		if (m_kdtree)
		{
			delete m_kdtree;
			m_kdtree = 0;
		}
		if (tetras.empty()) return;

		std::vector <CVec<Real, 3>> tetcenters(tetras.size());
		std::vector< ptrdiff_t> tetids(tetras.size());
		for (size_t i = 0; i < tetras.size(); i++)
		{
			tetcenters[i] = (Real)0.25*(tetras[i]->vertex[0]->pos + tetras[i]->vertex[1]->pos +
				tetras[i]->vertex[2]->pos + tetras[i]->vertex[3]->pos);
			tetids[i] = i;
		}
		m_kdtree = new KD_Tree<Real, 3>(tetcenters, tetids);
	}
	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	Tetrahedron<Real>* TetStructure<Real>::nearest_tetra_center(const CVec<Real, 3>& point)
	{
		if (m_kdtree == nullptr) return nullptr;
		CVec<Real, 3> nearestpoint;
		ptrdiff_t tet_id = 0;
		m_kdtree->find_Nearest_Point(point, nearestpoint, &tet_id);
		return tetras[tet_id];
	}
	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	Tetrahedron<Real>* TetStructure<Real>::search_tetra(Tetrahedron<Real>* tet, const CVec<Real, 3>& point)
	{
		if (tet == nullptr) return nullptr;
		Real msignvolume[4];
		if (tet->is_inside(point, msignvolume))
		{
			return tet;
		}
		else
		{
			int fid = 0;
			Real vol = msignvolume[0];
			for (int i = 1; i < 4; i++)
			{
				if (msignvolume[i] < vol)
				{
					fid = i;
					vol = msignvolume[i];
				}
			}
			if (vol < -1.0e-10 && tet->neighborTet[fid])
			{
				Tetrahedron<Real>* mtet = search_tetra(tet->neighborTet[fid], point);
				if (mtet)
				{
					return mtet;
				}
			}
			else
				return tet;
		}
		return nullptr;
	}
	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	Tetrahedron<Real>* TetStructure<Real>::locate_point(const CVec<Real, 3>& point)
	{
		if (m_kdtree == nullptr)
		{
			build_kdtree();
		}
		if (m_kdtree == nullptr)
			return nullptr;

		//Tetrahedron<Real>* nt = nearest_tetra_center(point)£»
		//Tetrahedron<Real>* find_t = search_tetra(nt, point);
		//return find_t ? find_t : nt;

		return search_tetra(nearest_tetra_center(point), point);
	}

	template<typename Real>
	bool TetStructure<Real>::is_inside(const CVec<Real, 3>& point)
	{
		Tetrahedron<Real>* mt = locate_point(point);
		//Real v[4];
		Real msignvolume[4];


		CVec<Real, 3> dir[4];
		dir[0] = mt->vertex[0]->pos - point;
		dir[1] = mt->vertex[1]->pos - point;
		dir[2] = mt->vertex[2]->pos - point;
		dir[3] = mt->vertex[3]->pos - point;

		msignvolume[0] = dir[3].Dot(dir[1].Cross(dir[2]));
		msignvolume[1] = -dir[3].Dot(dir[0].Cross(dir[2]));
		msignvolume[2] = dir[3].Dot(dir[0].Cross(dir[1]));
		msignvolume[3] = -dir[2].Dot(dir[0].Cross(dir[1]));

		/*if (msignvolume[0] >= -1.0e-10 &&
			msignvolume[1] >= -1.0e-10 &&
			msignvolume[2] >= -1.0e-10 &&
			msignvolume[3] >= -1.0e-10)*/
		if (msignvolume[0] >= 0 &&
			msignvolume[1] >= 0 &&
			msignvolume[2] >= 0 &&
			msignvolume[3] >= 0)
		{
			return true;
		}
		return false;



		//return mt->is_inside(point, v);


	}


	//////////////////////////////////////////////////////////////////////////
	template class TetVertex<double>;
	template class Tetrahedron<double>;
	template class TetStructure<double>;
	template class TetVertex<float>;
	template class Tetrahedron<float>;
	template class TetStructure<float>;
}