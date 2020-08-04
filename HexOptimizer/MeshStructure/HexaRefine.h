#pragma once

#include "HexaSmoothBase.h"

namespace MeshLib
{
	enum FeatureLineType { UNDEFINED, CORNER, SEGMENT, LOOP };
	template <typename Real>
	class FeatureLine
	{
	public:
		FeatureLine()
			:segment_tree(0), type(UNDEFINED)
		{
		}
		~FeatureLine()
		{
			if (segment_tree)
				delete segment_tree;
		}
		bool build_segment_tree(const std::vector < TinyVector<Real, 3>>& vertices)
		{
			if (vertices.size() == 1)
			{
				type = CORNER;
				edges.resize(0);
				edges.push_back(vertices[0]);
			}
			else if (vertices.size() >= 2 && vertices.front() != vertices.back())
			{
				type = SEGMENT;
				edges.resize(0);
				edges.reserve(2 * vertices.size() - 2);
				for (size_t i = 0; i < vertices.size() - 1; i++)
				{
					edges.push_back(vertices[i]);
					edges.push_back(vertices[i + 1]);
				}
				TinyVector<Real, 3> E[2];
				E[0] = edges.front(), E[1] = edges.back();
				edges.front() = edges.front() + 0.01 * (edges[1] - edges.front());
				edges.back() = edges.back() + 0.01 * (edges[edges.size() - 2] - edges.back());
				segment_tree = new ig::AABB_Segment_Tree<Real>(edges);
				edges.front() = E[0], edges.back() = E[1];
			}
			else if (vertices.size() > 2 && vertices.front() == vertices.back())
			{
				type = LOOP;
				edges.resize(0);
				edges.reserve(2 * vertices.size());
				for (size_t i = 0; i < vertices.size() - 1; i++)
				{
					edges.push_back(vertices[i]);
					edges.push_back(vertices[i + 1]);
				}
				segment_tree = new ig::AABB_Segment_Tree<Real>(edges);
			}
			else
				return false;

			return true;
		}
		bool project(const TinyVector<Real, 3>& p, TinyVector<Real, 3>& vec, ptrdiff_t* seg_id = 0)
		{
			if (type == UNDEFINED) return false;
			else if (type == CORNER)
			{
				vec = edges.front(); if (seg_id) *seg_id = 0;
				return true;
			}
			else
			{
				if (!segment_tree) return false;
				segment_tree->find_Nearest_Point(p, vec, seg_id);
				return true;
			}
		}
		const std::vector<size_t>& get_vert_indices() const { return vert_indices; }
		std::vector<size_t>& get_vert_indices() { return vert_indices; }
		const std::vector< TinyVector<Real, 3>>& get_edges() const { return edges; }
		const FeatureLineType& get_type() const { return type; }
		FeatureLineType& get_type() { return type; }
	private:
		std::vector<size_t> vert_indices;
		std::vector< TinyVector<Real, 3>> edges;
		ig::AABB_Segment_Tree<Real>* segment_tree;
		FeatureLineType type;
	};

	template <typename Real>
	class HexaRefine : public HexaSmoothBase<Real>
	{
	public:
		HexaRefine(VHexMesh<Real>* hexmesh, Mesh3D<Real>* bndmesh = 0, int max_iter = 20, 
			bool fixboundary = false, bool surf_preserve = false, double threshold = 0, 
			bool use_ejac = false, bool opt_illegal_only = false);
		~HexaRefine();
		bool load_features(const char trifeaturename[], const char hexfeaturename[]);
		void dump(const char meshfile[], const char trifeaturename[], const char quadname[], const char hexfeaturename[]);
		void save_trifeature_polyline(const char trifeaturename[]);
		void save_hexfeature_polyline(const char hexfeaturename[]);
	protected:
		void clear_feature();
		void prepare_feature_edges_for_trimesh(const std::vector<std::vector<size_t>>& newlines);
		void prepare_feature_edges_for_hexmesh(const std::vector<std::vector<size_t>>& newlines);
		void reorganize_linesegment(const std::vector < std::vector<size_t>>& shortlines, std::vector < std::vector<size_t>>& lines);
		void feature_projection(size_t surf_vertex_id, const TinyVector<Real, 3>& cur_point, TinyVector<Real, 3>& projection_point);
		void corner_projection(size_t surf_vertex_id, const TinyVector<Real, 3>& cur_point, TinyVector<Real, 3>& projection_point);
		bool has_feature();
	private:
		std::vector<FeatureLine<Real>*> trimesh_featurelines, hexmesh_featurelines;
		std::vector<size_t> feature_edge_targetedgeID;
		std::vector<ptrdiff_t> feature_point_targetID;
	};
}