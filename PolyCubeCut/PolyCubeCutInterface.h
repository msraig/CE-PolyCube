#ifndef _POLYCUBE_CUT_INTERFACE_H_
#define _POLYCUBE_CUT_INTERFACE_H_

#include <vector>
#include <set>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "algebra/SmallVec.h"

class PolyCubeCut;

class ThreeCut
{
public:
    int cut_type;
    unsigned int shared_vertex_id;
    std::array<int, 3> a_adjacent_one_cut_index;
    std::vector<std::array<unsigned int, 3>> va_corresponding_vertex;
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3> Polyhedron;
typedef Polyhedron::Halfedge_handle Halfedge_handle;

class PolyCubeCutInterface
{
public:
    PolyCubeCutInterface(
        const std::vector<ig::CVec<double, 3>> &v_input_point,
        const std::vector<unsigned int> &v_input_face,
        const std::vector<int> &v_input_chart,
        const std::vector<int> &v_input_label
    );
    ~PolyCubeCutInterface();

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
    PolyCubeCut *m_pcc;
};

#endif