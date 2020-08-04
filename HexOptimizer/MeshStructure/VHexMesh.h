#ifndef VHEXMESH_H
#define VHEXMESH_H

#include "Mesh3D.h"
#include <set>

namespace MeshLib
{
	/*
				   v
			3----------2
			|\     ^   |\
			| \    |   | \
			|  \   |   |  \
			|   7------+---6
			|   |  +-- |-- | -> u
			0---+---\--1   |
			 \  |    \  \  |
			  \ |     \  \ |
			   \|      w  \|
				4----------5
	*/
	template <typename Real>
	class HexElement
	{
	public:
		HexElement()
		{
			memset(indices, 0, sizeof(size_t) * 8);
			clear_bnd();
		}
		HexElement(size_t ind[])
		{
			memcpy(indices, ind, sizeof(size_t) * 8);
			clear_bnd();
		}

		HexElement(size_t ind0, size_t ind1, size_t ind2, size_t ind3, size_t ind4, size_t ind5, size_t ind6, size_t ind7)
		{
			indices[0] = ind0;
			indices[1] = ind1;
			indices[2] = ind2;
			indices[3] = ind3;
			indices[4] = ind4;
			indices[5] = ind5;
			indices[6] = ind6;
			indices[7] = ind7;
			clear_bnd();
		}

		inline HexElement<Real> &operator =(const HexElement<Real> &tc)
		{
			memcpy(indices, tc.indices, sizeof(size_t) * 8);
			memcpy(neighborHex, tc.neighborHex, sizeof(std::pair<HexElement<Real>*, int>) * 6);
			return *this;
		}
		const int get_local_id(size_t id) const
		{
			for (int j = 0; j < 8; j++)
			{
				if (indices[j] == id)
				{
					return j;
				}
			}
			return -1;
		}
		bool is_on_boundary()
		{
			for (int j = 0; j < 6; j++)
			{
				if (neighborHex[j].first == 0)
				{
					return true;
				}
			}
			return false;
		}
		size_t indices[8];
		std::pair<HexElement<Real>*, int> neighborHex[6];
		size_t id; // only used in orientation
	private:
		void clear_bnd()
		{
			for (int i = 0; i < 6; i++)
			{
				neighborHex[i].first = 0;
				neighborHex[i].second = -1;
			}
		}
	};

	class QuadFace
	{
	public:
		QuadFace()
		{
			quad_vert[0] = quad_vert[1] = quad_vert[2] = quad_vert[3] = 0;
		}
		QuadFace(const size_t A, const size_t B, const size_t C,
			const size_t D)
		{
			set_vertices(A, B, C, D);
		}
		QuadFace(const QuadFace &qf)
		{
			memcpy(quad_vert, qf.quad_vert, sizeof(size_t) * 4);
		}
		inline void set_vertices(const size_t A, const size_t B,
			const size_t C, const size_t D)
		{
			size_t vert[4];
			vert[0] = A;
			vert[1] = B;
			vert[2] = C;
			vert[3] = D;
			unsigned int min_id = 0;
			size_t min_V = vert[0];
			for (int i = 1; i < 4; i++)
			{
				if (vert[i] < min_V)
				{
					min_V = vert[i];
					min_id = i;
				}
			}
			quad_vert[0] = vert[min_id];
			if (vert[(min_id + 1) % 4] < vert[(min_id + 3) % 4])
			{
				quad_vert[1] = vert[(min_id + 1) % 4];
				quad_vert[2] = vert[(min_id + 2) % 4];
				quad_vert[3] = vert[(min_id + 3) % 4];
			}
			else
			{
				quad_vert[1] = vert[(min_id + 3) % 4];
				quad_vert[2] = vert[(min_id + 2) % 4];
				quad_vert[3] = vert[(min_id + 1) % 4];
			}
		}
		inline bool operator<(const QuadFace &m_r) const
		{
			for (int i = 0; i < 4; i++)
			{
				if (quad_vert[i] < m_r.quad_vert[i])
				{
					return true;
				}
				if (quad_vert[i] > m_r.quad_vert[i])
				{
					return false;
				}
			}
			return false;
		}
		size_t quad_vert[4];
	};

	class QuadFace_O
	{
	public:
		QuadFace_O()
		{
			quad_vert[0] = quad_vert[1] = quad_vert[2] = quad_vert[3] = 0;
		}
		QuadFace_O(const size_t A, const size_t B, const size_t C,
			const size_t D)
		{
			set_vertices(A, B, C, D);
		}
		QuadFace_O(const QuadFace_O &qf)
		{
			memcpy(quad_vert, qf.quad_vert, sizeof(size_t) * 4);
		}
		inline void set_vertices(const size_t A, const size_t B,
			const size_t C, const size_t D)
		{
			size_t vert[4];
			vert[0] = A;
			vert[1] = B;
			vert[2] = C;
			vert[3] = D;
			unsigned int min_id = 0;
			size_t min_V = vert[0];
			for (int i = 1; i < 4; i++)
			{
				if (vert[i] < min_V)
				{
					min_V = vert[i];
					min_id = i;
				}
			}
			quad_vert[0] = vert[min_id];
			quad_vert[1] = vert[(min_id + 1) % 4];
			quad_vert[2] = vert[(min_id + 2) % 4];
			quad_vert[3] = vert[(min_id + 3) % 4];
		}
		inline bool operator<(const QuadFace_O &m_r) const
		{
			for (int i = 0; i < 4; i++)
			{
				if (quad_vert[i] < m_r.quad_vert[i])
				{
					return true;
				}
				if (quad_vert[i] > m_r.quad_vert[i])
				{
					return false;
				}
			}
			return false;
		}
		size_t quad_vert[4];
	};


	template <typename Real>
	class VHexMesh
	{
	public:
		VHexMesh()
		{
			num_quad_faces = 0;
			boundingbox[0][0] = boundingbox[2][0] = boundingbox[1][0] = 0;
			boundingbox[0][1] = boundingbox[2][1] = boundingbox[1][1] = 1;
		}
		~VHexMesh();
		void load_hex(const char filename[]);
		void save_hex(const char filename[]);
		void save_reversed_hex(const char filename[]);
		void load_vtk(const char filename[]);
		void save_vtk(const char filename[]);
		void load_ovm(const char filename[]);
		void save_ovm(const char filename[]);
		void load_mesh(const char filename[]);
		void save_mesh(const char filename[]);
		void load_data(const std::vector<TinyVector<double, 3>>& hex_vertices, const std::vector<size_t>& hex_indices, bool vtk_order = false);
		inline const std::vector< TinyVector<Real, 3> > &get_vertices() const { return vertices; }
		inline const std::vector< HexElement<Real>* > &get_elements() const { return elements; }
		inline std::vector< TinyVector<Real, 3> > &get_vertices() { return vertices; }
		inline std::vector< HexElement<Real>* > &get_elements() { return elements; }
		inline const size_t get_num_quad_faces() const { return num_quad_faces; }
		inline void get_boundingbox(Real &xmin, Real &xmax, Real &ymin, Real &ymax, Real &zmin, Real &zmax)
		{
			xmin = boundingbox[0][0];
			xmax = boundingbox[0][1];
			ymin = boundingbox[1][0];
			ymax = boundingbox[1][1];
			zmin = boundingbox[2][0];
			zmax = boundingbox[2][1];
		}
		void rebuild_mesh() { build_connectivity(); reorientation(); update_boundingbox(); }
		Mesh3D<Real> *extract_boundary_quad_mesh(std::map<size_t, size_t> &vertexmap);
		bool detect_orientation();
		void merge_hexes(VHexMesh<Real> *hexmesh, double DIST_THRES = 1.0e-10);
		void merge_vertices(double DIST_THRES = 1.0e-10);
		void subdivide();
		void padding();
		//////////////////////////////////////////////////////////////////////////
		void remove_unused_vertices(bool removed = true);
	private:
		void clear_mem();
		void build_connectivity();
		void update_boundingbox();
		void reorientation();
	private:
		std::vector< TinyVector<Real, 3> > vertices;
		std::vector< HexElement<Real>* > elements;
		size_t num_quad_faces;
		Real boundingbox[3][2];
	};

	typedef MeshLib::VHexMesh<float> VHexMesh3f;
	typedef MeshLib::VHexMesh<double> VHexMesh3d;
}

#endif
