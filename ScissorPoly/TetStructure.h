#ifndef TET_STRUCTURE_H
#define TET_STRUCTURE_H

#include <vector>
#include <algorithm>
#include <cassert>
#include "SmallVec.h"
#include "KD_Tree.h"

namespace ig
{
	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	class TetVertex
	{
	public:
		TetVertex() :
			boundary(false), id(-1)
		{
		}
		//////////////////////////////////////////////////////////////////////////
		TetVertex(const CVec<Real, 3>& mpos) :
			pos(mpos), boundary(false), id(-1)
		{
		}
	public:
		bool boundary;
		ptrdiff_t id;
		CVec<Real, 3> pos;
	};
	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	class Tetrahedron
	{
	public:
		//////////////////////////////////////////////////////////////////////////
		Tetrahedron()
		{
			vertex[0] = vertex[1] = vertex[2] = vertex[3] = 0;
			neighborTet[0] = neighborTet[1] = neighborTet[2] = neighborTet[3] = 0;
		}
		//////////////////////////////////////////////////////////////////////////
		Tetrahedron(TetVertex<Real>* v0, TetVertex<Real>* v1, TetVertex <
			Real > * v2, TetVertex<Real>* v3)
		{
			set_vertices(v0, v1, v2, v3);
		}
		//////////////////////////////////////////////////////////////////////////
		inline void set_vertices(TetVertex<Real>* v0, TetVertex<Real>* v1, TetVertex<Real>* v2, TetVertex<Real>* v3, bool swap = false)
		{
			vertex[0] = v0;
			vertex[1] = v1;
			vertex[2] = v2;
			vertex[3] = v3;
			neighborTet[0] = neighborTet[1] = neighborTet[2] = neighborTet[3] = 0;
			if (swap && compute_tet_volume() < 0)
			{
				std::swap(vertex[0], vertex[1]);
			}
		}
		//////////////////////////////////////////////////////////////////////////
		inline Real compute_tet_volume()
		{
			const CVec<Real, 3> v10 = vertex[1]->pos - vertex[0]->pos;
			const CVec<Real, 3> v20 = vertex[2]->pos - vertex[0]->pos;
			const CVec<Real, 3> v30 = vertex[3]->pos - vertex[0]->pos;
			return v30.Dot(v10.Cross(v20)) / (Real)6.0;
		}
		//////////////////////////////////////////////////////////////////////////
		inline void set_neighbor_tet(Tetrahedron<Real>* t)
		{
			assign_adjacent_id(t);
		}
		//////////////////////////////////////////////////////////////////////////
		inline void set_neighbor_tet(Tetrahedron<Real>* t, unsigned int id) //use this function carefully
		{
			neighborTet[id] = t;
		}
		//////////////////////////////////////////////////////////////////////////
		inline unsigned int get_local_vertex_id(TetVertex<Real>* v)
		{
			for (unsigned int i = 0; i < 4; i++)
			{
				if (v == vertex[i])
				{
					return i;
				}
			}

			assert(false);
			return 0;
		}
		//////////////////////////////////////////////////////////////////////////
		inline bool has_vertex(TetVertex<Real>* v)
		{
			for (unsigned int i = 0; i < 4; i++)
			{
				if (v == vertex[i])
				{
					return true;
				}
			}

			return false;
		}
		//////////////////////////////////////////////////////////////////////////
		inline TetVertex<Real>* get_opposite_vertex(const unsigned int id,
			unsigned int* opposite_id = 0)
		{
			if (id > 3 || neighborTet[id] == 0)
			{
				return 0;
			}

			for (unsigned int i = 0; i < 4; i++)
			{
				if (neighborTet[id]->vertex[i] != vertex[(id + 1) % 4]
					&& neighborTet[id]->vertex[i] != vertex[(id + 2) % 4]
					&& neighborTet[id]->vertex[i] != vertex[(id + 3) % 4])
				{
					if (opposite_id)
					{
						*opposite_id = i;
					}

					return neighborTet[id]->vertex[i];
				}
			}

			return 0;
		}
		//////////////////////////////////////////////////////////////////////////
		inline bool is_onboundary()
		{
			return !(neighborTet[0] && neighborTet[1] && neighborTet[2]
				&& neighborTet[3]);
		}
		//////////////////////////////////////////////////////////////////////////
		inline Tetrahedron<Real>* get_neighborTet_by_vertex(TetVertex<Real>* v)
		{
			return neighborTet[get_local_vertex_id(v)];
		}
		//////////////////////////////////////////////////////////////////////////
		bool is_inside(const CVec<Real, 3>& point, Real msignvolume[])
		{
			if (compute_tet_volume() <= 0)
				return true;

			CVec<Real, 3> dir[4];
			dir[0] = vertex[0]->pos - point;
			dir[1] = vertex[1]->pos - point;
			dir[2] = vertex[2]->pos - point;
			dir[3] = vertex[3]->pos - point;

			msignvolume[0] = dir[3].Dot(dir[1].Cross(dir[2]));
			msignvolume[1] = -dir[3].Dot(dir[0].Cross(dir[2]));
			msignvolume[2] = dir[3].Dot(dir[0].Cross(dir[1]));
			msignvolume[3] = -dir[2].Dot(dir[0].Cross(dir[1]));

			if (msignvolume[0] >= 0 &&
				msignvolume[1] >= 0 &&
				msignvolume[2] >= 0 &&
				msignvolume[3] >= 0)
			{
				return true;
			}
			return false;
		}
	private:
		inline bool assign_adjacent_id(Tetrahedron<Real>* t) //assume t is an adjacent tetra
		{
			if (t)
			{
				for (unsigned int i = 0; i < 4; i++)
				{
					bool find_shared_vertex = false;

					for (unsigned int j = 0; j < 4; j++)
					{
						if (vertex[i] == t->vertex[j])
						{
							find_shared_vertex = true;
							break;
						}
					}

					if (find_shared_vertex == false)
					{
						neighborTet[i] = t;
						return true;
					}
				}
			}

			return false;
		}
		//////////////////////////////////////////////////////////////////////////
	public:
		TetVertex<Real>* vertex[4];
		Tetrahedron<Real>* neighborTet[4];
		size_t id;
	};
	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	class MyIntTriple
	{
	public:
		MyIntTriple(TetVertex<Real>* i, TetVertex<Real>* j,
			TetVertex<Real>* k)
		{
			id[0] = i;
			id[1] = j;
			id[2] = k;
			std::sort(id, id + 3);
		}
		//! The compare function for sorting
		inline bool operator< (const MyIntTriple<Real>& m_r) const
		{
			for (unsigned int i = 0; i < 3; i++)
			{
				if (id[i] < m_r.id[i])
				{
					return true;
				}

				if (id[i] > m_r.id[i])
				{
					return false;
				}
			}

			return false;
		}

		inline bool operator == (const MyIntTriple<Real>& m_r) const
		{
			for (unsigned int i = 0; i < 3; i++)
			{
				if (id[i] != m_r.id[i])
				{
					return false;
				}
			}

			return true;
		}

		inline bool operator != (const MyIntTriple<Real>& m_r) const
		{
			return id[0] != m_r.id[0] || id[1] != m_r.id[1] || id[2] != m_r.id[2];
		}
	public:
		TetVertex<Real>* id[3];
	};
	//////////////////////////////////////////////////////////////////////////
	template<typename Real>
	class TetStructure
	{
	public:
		TetStructure();
		~TetStructure();

		void load_tet(const std::vector<CVec<Real, 3> >& points, const std::vector <unsigned int>& indices);
		void load_tet(const char filename[]);
		void save_tet(const char filename[]);
		void load_vtk(const char filename[]);
		void save_vtk(const char filename[]);

		void find_adjacent_tets_by_edge(Tetrahedron<Real> *tet, const unsigned int vi, const unsigned int vj, std::vector<Tetrahedron<Real>*> &adjacent_tets);
		void find_adjacent_tets_by_edge(Tetrahedron<Real> *tet, const unsigned int vi, const unsigned int vj, std::vector<Tetrahedron<Real>*> &adjacent_tets, std::vector<TetVertex<Real>*> &loop);

		const std::vector<Tetrahedron<Real>*>& get_tetras() const
		{
			return tetras;
		}
		std::vector<Tetrahedron<Real>*>& get_tetras()
		{
			return tetras;
		}
		const std::vector<TetVertex<Real>*>& get_vertices() const
		{
			return tetra_vertices;
		}
		std::vector<TetVertex<Real>*>& get_vertices()
		{
			return tetra_vertices;
		}

		void build_kdtree();
		Tetrahedron<Real>* locate_point(const CVec<Real, 3>& point);
		bool is_inside(const CVec<Real, 3>& point);

	private:
		void build_connectivity();
		void clear();

		Tetrahedron<Real>* nearest_tetra_center(const CVec<Real, 3>& point);
		Tetrahedron<Real>* search_tetra(Tetrahedron<Real>* tet, const CVec<Real, 3>& point);
	public:
		std::vector<Tetrahedron<Real>*> tetras;
		std::vector<TetVertex<Real>*> tetra_vertices;

	private:
		KD_Tree<Real, 3>* m_kdtree;
	};
}
#endif
