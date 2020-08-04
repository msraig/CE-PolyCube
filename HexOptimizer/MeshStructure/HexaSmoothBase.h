#pragma once

#include "CellSmooth.h"
#include "VHexMesh.h"
#include "AABB_Tree.h"
#include "AABB_Segment_Tree.h"
#include <unordered_set>

namespace MeshLib
{
	template <typename Real>
	class OptData
	{
	public:
		OptData() { proj[0] = proj[1] = proj[2] = 0; edgelen = 1; weight = 1e8; }
	public:
		TinyVector<Real, 3> proj;
		Real edgelen, weight;
	};

	template <typename Real>
	class HexaSmoothBase
	{
	public:
		HexaSmoothBase(VHexMesh<Real>* hexmesh = 0, Mesh3D<Real>* bndmesh = 0, int max_iter = 20, bool fixboundary = false, bool surf_preserve = false,
			double threshold = 0, bool use_ejac = false, bool opt_illegal_only = false);
		virtual ~HexaSmoothBase();
		void smooth(bool laplaciansmooth, bool verbose = true);
		void save_quadmesh(const char quadname[]);
		void save_trimesh(const char triname[]);
		void save_trifeature(const char feaname[]);
		void save_quadfeature(const char feaname[]);
		virtual bool load_features(const char trifeaturename[], const char hexfeaturename[]);
	protected:
		void initialize();
		bool load_feature_edge(const char filename[], ig::KD_Tree<Real, 3>& kdtree, std::vector<std::vector<size_t>>& alllines);
		bool load_feature_edge(const char filename[], std::vector<std::vector<size_t>>& alllines);
		Real compute_distortion_and_gradient(const size_t id, const TinyVector<Real, 3>& vertex, const OptData<Real>& data,
			bool& has_flip, TinyVector<Real, 3>* gradient = 0);
		void optimize(int max_iter = 100, bool fixboundary = false, bool verbose = true);
		void print_scalejacobian();
		void graph_coloring();
		virtual void feature_projection(size_t surf_vertex_id, const TinyVector<Real, 3>& cur_point, TinyVector<Real, 3>& projection_point);
		virtual void corner_projection(size_t surf_vertex_id, const TinyVector<Real, 3>& cur_point, TinyVector<Real, 3>& projection_point);
		Real element_scaledjacobian(const int corner_type, const TinyVector<Real, 3>& V0, const TinyVector<Real, 3>& V1, const TinyVector<Real, 3>& V2);
		void reorganize(const std::vector < std::vector<size_t>>& lines, std::vector < std::vector<size_t>>& newlines);
		void reorganize_all_features(
			const std::vector < std::vector<size_t>>& trilines,
			std::vector < std::vector<size_t>>& trinewlines,
			const std::vector < std::vector<size_t>>& hexlines,
			std::vector < std::vector<size_t>>& hexnewlines
		);
		void laplacian_smooth(unsigned int iter, bool only_boundary, bool projection_only = false);
		void fast_laplacian_smooth(unsigned int iter, bool only_boundary);
		void laplacian_smooth_init();
		void extract_irrugular_lines(std::vector<bool>& vertexirregularlabel, bool detect_structure = false);
		virtual bool has_feature();
		void label_feature_edges();
	protected:
		const unsigned int corner_vert[32][4] = {
			{0,1,3,4 }, { 1,2,0,5 },{ 3,0,2,7}, { 2,3,1,6 }, { 4,5,0,7 }, { 5,1,4,6 },{7,6,4,3},{6,5,7,2},
		{0,1,2,6},{1,2,3,7},{2,3,0,4},{3,0,1,5},{4,7,6,2},{7,6,5,1},
		{6,5,4,0},{5,4,7,3},{1,0,4,7},{3,2,6,5},{0,3,7,6},{2,1,5,4},
		{0,1,2,4},{1,2,3,5},{2,3,0,6},{3,0,1,7},{4,7,6,0},{7,6,5,3},
		{6,5,4,2},{5,4,7,1},{4,0,3,5},{2,6,7,1},{3,7,4,2},{1,5,6,0}
		};
		const unsigned int face_id[6][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 7, 6, 2, 3 }, { 5, 4, 0, 1 }, { 5, 1, 2, 6 }, { 4, 7, 3, 0 } };
		const int edgeid[12][2] = { {0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7} };
		const Real K1 = sqrt(6.0) / 6;
		const Real K2 = sqrt(2.0) / 2;
		const Real threshold = 1.0e-12;
		Real bound = 0.2;
		int edge_sample = 5;
		Real ave_edgelength;
		bool fix_boundary, feature_sensitive, use_extend_jac, opt_illegal_vertex_only;
		int max_iters, sample_rate = 3;
		VHexMesh<Real>* m_hexmesh;
		Mesh3D<Real>* m_bndmesh, * m_quadmesh;
		ig::AABB_Tree<Real>* m_aabb;
		std::vector<std::vector<CellCorner>> vertex2cornerlist;
		std::vector<bool> boundaryvertexlabel;
		std::vector<CellFace> boundaryfacets;
		std::vector<std::vector<size_t>> colored_vertices, colored_facets;
		std::map<size_t, size_t> vert_map_ids;
		std::vector<size_t> reverse_ids;
		//vertex_feature_tags: -1(normal surf point), 0 (end point), 1 (edge point)
		std::vector<int> vertex_feature_tags;
		std::vector<size_t> feature_edges, trifeature_edges;
		std::vector<std::unordered_set<size_t>> neighor_vertices;
		std::vector<bool> featureedgelabels;
	};
}