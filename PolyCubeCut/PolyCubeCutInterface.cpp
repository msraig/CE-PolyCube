#include "PolyCubeCutInterface.h"
#include "PolyCubeCut.h"

PolyCubeCutInterface::PolyCubeCutInterface(
    const std::vector<ig::CVec<double, 3>> &v_input_point,
    const std::vector<unsigned int> &v_input_face,
    const std::vector<int> &v_input_chart,
    const std::vector<int> &v_input_label
)
{
    m_pcc = new PolyCubeCut(v_input_point, v_input_face, v_input_chart, v_input_label);
}

PolyCubeCutInterface::~PolyCubeCutInterface()
{
    if (m_pcc)
        delete m_pcc;
}

void PolyCubeCutInterface::get_cuttable_edge(std::vector<unsigned int> &v_cuttable_edge)
{
    m_pcc->get_cuttable_edge(v_cuttable_edge);
}

void PolyCubeCutInterface::get_associated_edge(
	std::map<
		std::pair<unsigned int, unsigned int>,
		std::vector<std::pair<unsigned int, unsigned int>>
	> &map_associated_edge
)
{
	m_pcc->get_associated_edge(map_associated_edge);
}

void PolyCubeCutInterface::get_convex_map(std::map<std::pair<unsigned int, unsigned int>, bool> &map_is_convex_edge)
{
	m_pcc->get_convex_map(map_is_convex_edge);
}

bool PolyCubeCutInterface::is_equivalent_cut(
    const std::vector<unsigned int> &v_cut_edge_a,
    const std::vector<double> &v_cut_depth_a,
    const std::vector<unsigned int> &v_cut_edge_b,
    const std::vector<double> &v_cut_depth_b
)
{
    return m_pcc->is_equivalent_cut(v_cut_edge_a, v_cut_depth_a, v_cut_edge_b, v_cut_depth_b);
}

void PolyCubeCutInterface::cut_test(
	const std::vector<unsigned int> &v_cut_edge,
	const std::vector<double> &v_cut_depth,
	const double &max_volume,
	const double &tet_quality
)
{
	m_pcc->cut_test(v_cut_edge, v_cut_depth, max_volume, tet_quality);
}

void PolyCubeCutInterface::cut(
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
)
{
    m_pcc->cut(
        v_cut_edge, v_cut_depth, max_volume, tet_quality,
        v_output_point, v_output_tet, v_output_chart, v_output_label, v_cut_type,
        v_cut_section_chart, v_seam_edge_vertex_id, vv_congruent_face, v_three_cut
    );
}

void PolyCubeCutInterface::cut(
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
)
{
    m_pcc->cut(
        v_cut_edge, v_cut_depth, max_volume, tet_quality,
        v_output_point, v_output_tet, v_output_chart, v_output_label, v_cut_type,
        v_cut_section_chart, v_seam_edge_vertex_id, vv_congruent_face, v_three_cut
    );
}

void PolyCubeCutInterface::get_cut_edge_map(std::map<std::pair<unsigned int, unsigned int>, std::vector<std::pair<unsigned int, unsigned int>>>& um_cut_edge_map)
{
    m_pcc->get_cut_edge_map(um_cut_edge_map);
}
