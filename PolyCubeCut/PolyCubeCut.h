#ifndef _POLYCUBE_CUT_H_
#define _POLYCUBE_CUT_H_

#include "PolyCubeCutHelper.h"
#include "algebra/SmallVec.h"
#include "PolyCubeCutInterface.h"

class PolyCubeCut
{
public:
	~PolyCubeCut();

    PolyCubeCut(
        const std::vector<ig::CVec<double, 3>> &v_input_point,
        const std::vector<unsigned int> &v_input_face,
        const std::vector<int> &v_input_chart,
        const std::vector<int> &v_input_label
    );

    void get_cuttable_edge(std::vector<unsigned int> &v_cuttable_edge);

	void get_associated_edge(
		std::map<
			std::pair<unsigned int, unsigned int>,
			std::vector<std::pair<unsigned int, unsigned int>>
		> &map_associated_edge
	);

	void get_convex_map(std::map<std::pair<unsigned int, unsigned int>, bool> &map_is_convex_edge);

    bool is_equivalent_cut(
        const std::vector<unsigned int> &v_cut_edge_a,
        const std::vector<double> &v_cut_depth_a,
        const std::vector<unsigned int> &v_cut_edge_b,
        const std::vector<double> &v_cut_depth_b
    );

	void cut_test(
		const std::vector<unsigned int> &v_cut_edge,
		const std::vector<double> &v_cut_depth,
		const double &max_volume,
		const double &tet_quality
	);

	void cut(
        const std::vector<unsigned int> &v_cut_edge,
        const std::vector<double> &v_cut_depth,
        const double &max_volume,
        const double &tet_quality,
        std::vector<ig::CVec<double, 3>> &v_output_point,
        std::vector<unsigned int> &v_output_tet,
        std::vector<int> &v_output_chart,
        std::vector<int> &v_output_label,
        std::vector<int> &v_cut_type,
        std::vector<int> &v_cut_section_chart,
        std::vector<int> &v_seam_edge_vertex_id,
        std::vector<std::vector<int>> &vv_congruent_face,
        std::vector<ThreeCut> &v_three_cut
    );

    void cut(
        const std::vector<Halfedge_handle> &v_cut_edge,
        const std::vector<double> &v_cut_depth,
        const double &max_volume,
        const double &tet_quality,
        std::vector<ig::CVec<double, 3>> &v_output_point,
        std::vector<unsigned int> &v_output_tet,
        std::vector<int> &v_output_chart,
        std::vector<int> &v_output_label,
        std::vector<int> &v_cut_type,
        std::vector<int> &v_cut_section_chart,
        std::vector<int> &v_seam_edge_vertex_id,
        std::vector<std::vector<int>> &vv_congruent_face,
        std::vector<ThreeCut> &v_three_cut
    );

    void get_cut_edge_map(std::map<std::pair<unsigned int, unsigned int>, std::vector<std::pair<unsigned int, unsigned int>>>& um_cut_edge_map);

private:
	Polyhedron &m_polymesh;
    double m_threshold;
	int m_initial_chart_size;
	std::vector<Vertex_handle> m_v_id2vertex;
    std::vector<Halfedge_handle> m_v_id2halfedge;
	static const std::array<std::map<int, int>, 4> sm_am_label2cutidx;
	static const std::array<Vector, 12> sm_a_cut_idx2dir;
	static const std::array<Vector, 12> sm_a_cut_idx2normal;
    static const std::array<Vector, 6> sm_a_label2normal;
	static const std::array<int, 24> sm_a_cutidx2label;
	std::vector<int> m_v_chart2label;
	std::vector<int> m_v_face2chart;
	std::vector<int> m_v_face2label;
    std::vector<Halfedge_handle> m_v_polycube_edge;
    std::vector<std::array<Halfedge_handle, 3>> m_va_vertex2_a_polycube_edge;
    std::map<Halfedge_handle, bool> m_map_is_convex_cut; // for a cuttable edge, true if it is convex, false otherwise
    std::map<Halfedge_handle, std::set<Halfedge_handle>> m_ms_associated_edge;
    std::vector<CutFacet> m_v_cut_facet;
    std::vector<std::array<int, 2>> m_v_duplicated_vertex_map;
    std::vector<std::vector<int>> m_vv_congruent_face;

    void index2halfedge(const std::vector<unsigned int> &v_vertex_index, std::vector<Halfedge_handle> *v_halfedge);

    void get_complete_cut_edge(
        const std::vector<Halfedge_handle> &v_cut_edge,
        const std::vector<double> &v_cut_depth,
        std::map<unsigned int, double> &map_cut2depth,
        std::vector<bool> &v_is_cut_edge,
        std::set<unsigned int> &s_complete_cut_edge
    );

    int vertex_cut_count(
        const std::vector<bool> &v_is_cut_edge,
        const std::map<unsigned int, double> &map_cut2depth,
        const Vertex_handle &v,
        Halfedge_handle &new_cut_edge,
        double &depth
    );

    void check_three_cut(
        const Vertex_handle &v,
        std::vector<bool> &v_is_cut_edge,
        std::set<unsigned int> &s_complete_cut_edge,
        std::map<unsigned int, double> &map_cut2depth
    );

    void adjust_cut_depth(
        const std::vector<bool> &v_is_cut_edge,
        std::map<unsigned int, double> &map_cut2depth
    );
    
    void build_triangle_mesh(
        const std::vector<ig::CVec<double, 3>> &v_input_point,
        const std::vector<unsigned int> &v_input_face
    );

    void find_polycube_edge();

    bool is_on_surface(const Halfedge_handle &he, const Vector &cut_dir);

    Halfedge_handle halfedge2edge(const Halfedge_handle &he);

    void construct_PLC(
        const std::set<unsigned int> &s_complete_cut_edge,
        const std::map<unsigned int, double> &map_cut2depth,
        const std::vector<int> &v_cut_type,
        const std::vector<bool> &v_is_cut_edge,
        std::vector<int> *v_cut_section_chart,
        std::map<int, ThreeCut> *map_vertex2three_cut,
        tetgenio *tetin
    );

    int opposite_label(int label);

    void add_cut(
        const Halfedge_handle &cut_he,
        const double &cut_depth,
        const int &cut_type,
        const std::vector<bool> &v_is_cut_edge,
        std::map<size_t, unsigned int> *map_vertex2cut_end_vertex_id,
        std::vector<bool> *v_is_interior,
        std::map<int, ThreeCut> *map_vertex2three_cut,
        TetgenHelper *tetgen_helper
    );

    unsigned int create_cut_end_vertex(
        const Vertex_handle &cut_start_vertex,
        const Halfedge_handle &cut_he,
        const double &cut_depth,
        const int &cut_type,
        const std::vector<bool> &v_is_cut_edge,
        std::map<int, ThreeCut> *map_vertex2three_cut,
        bool *interior,
        TetgenHelper *tetgen_helper
    );

	std::array<int, 2> find_cut_chart(const Halfedge_handle &cut_he);

    int three_cut_direction2index(const Vector &dir);

	bool is_facet_too_close(
		const CutFacet& fa,
		const CutFacet& fb
	);

	void check_cut_spacing(
		const CutFacet& new_cf
	);

	void check_intersection(
		const Halfedge_handle &cut_he,
		const std::array<Point, 2> &a_new_vertex
	);

	bool find_intersection_point_in_polygon(
		const Halfedge_handle &he_begin,
		const Halfedge_handle &he_end,
		const Segment &cut_segment,
		Halfedge_handle *he_intersected,
		Point *intersection_point
	);

	void do_cut_in_tetgen_helper(
		const std::array<Vertex_handle, 2> &a_cut_start_vertex_handle,
		const std::array<unsigned int, 2> &a_cut_end_vertex_id,
		const std::array<Point, 2> &a_cut_end_vertex,
		const Vector &cut_dir,
		const std::array<int, 2> &a_cut_chart,
		TetgenHelper *tetgen_helper
	);

	void cut_chart_in_tetgen_helper(
		const Vertex_handle &cut_start_vertex_handle,
		const unsigned int &cut_end_vertex_id,
		const int &cut_chart,
		const Vector &cut_dir,
		const Segment &cut_segment,
		TetgenHelper *tetgen_helper
	);

	void split_edge_in_tetgen_helper(const Halfedge_handle &he, const unsigned int &start, const unsigned int &end, TetgenHelper *tetgen_helper);

	void split_edge_in_polygon_helper(const Halfedge_handle &he, const unsigned int &start, const unsigned int &end, PolygonHelper *p);

	void insert_intersection_point(
		const int &n_intersection_point,
		const unsigned int &intersection_point_id,
		const unsigned int &start_or_end_vertex_id,
		const std::vector<Point> &pointlist,
		FacetHelper *f
	);

	void tetrahedron_subdivision(tetgenio *tetout);

	void update_boundary_tet_chart(int marker_max_index, tetgenio *tetout, std::vector<int> *v_tet_chart);

	void separate_cut_facet(
        const std::vector<int> &v_cut_type,
        const std::vector<int> &v_cut_section_chart,
        std::vector<int> *v_tet_chart,
        std::vector<std::vector<int>> *vv_congruent_face,
        std::map<int, ThreeCut> *map_vertex2three_cut,
        tetgenio *tetout
	);

	void duplicate_cut_facet_vertex(
        const std::set<int> &s_id_vertex_duplicate_1,
        const std::set<int> &s_id_vertex_duplicate_2,
        tetgenio *tetout
	);

	void update_tetrahedron(
		const std::vector<int> &v_cut_section_chart,
		std::vector<int> *v_tet_chart,
        std::vector<std::vector<int>> *vv_congruent_face,
		tetgenio *tetout
	);

    int region_on_tangent_plane(const unsigned int &vid, const Vector &normal, const Point &p);
};

#endif