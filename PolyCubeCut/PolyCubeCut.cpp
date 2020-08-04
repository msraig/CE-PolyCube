#include <queue>
#include <algorithm>
#include <unordered_set>

#include "PolyCubeCut.h"
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/intersections.h>

//#define RECONSTRUCTION_DETAIL

//void output_vtk(
//	const std::vector<ig::CVec<double, 3>>& v_output_point,
//	const std::vector<unsigned int>& v_output_tet,
//	const std::string& filename
//)
//{
//	std::ofstream output(filename);
//	output << "# vtk DataFile Version 2.0\n" << "Mesh saved from PolyCubeCut\n";
//	output << "ASCII\n" << "DATASET UNSTRUCTURED_GRID\n";
//	output << "POINTS " << v_output_point.size() << " float\n";
//	for (int i = 0; i < v_output_point.size(); ++i)
//	{
//		for (int j = 0; j < 3; ++j)
//			output << v_output_point[i][j] << " ";
//		output << std::endl;
//	}
//	int n_tet = static_cast<int>(v_output_tet.size() >> 2);
//	output << "CELLS " << n_tet << " " << 5 * n_tet << std::endl;
//	for (int i = 0; i < n_tet; ++i)
//	{
//		output << "4";
//		for (int j = 0; j < 4; ++j)
//			output << " " << v_output_tet[4 * i + j];
//		output << std::endl;
//	}
//	output << "CELL_TYPES " << n_tet << std::endl;
//	for (int i = 0; i < n_tet; ++i)
//		output << "10\n";
//	output.close();
//}

PolyCubeCut::~PolyCubeCut()
{
    delete &m_polymesh;
}

const std::array<std::map<int, int>, 4> PolyCubeCut::sm_am_label2cutidx{
    {
        {
            {2, 0},
            {3, 1},
            {4, 4},
            {5, 5}
        },
        {
            {2, 2},
            {3, 3},
            {4, 6},
            {5, 7}
        },
        {
            {4, 8},
            {5, 9}
        },
        {
            {4, 10},
            {5, 11}
        }
    }
};

const std::array<Vector, 12> PolyCubeCut::sm_a_cut_idx2dir{
	Vector(-1, -1,  0),
	Vector(-1,  1,  0),
	Vector( 1, -1,  0),
	Vector( 1,  1,  0),
	Vector(-1,  0, -1),
	Vector(-1,  0,  1),
	Vector( 1,  0, -1),
	Vector( 1,  0,  1),
	Vector( 0, -1, -1),
	Vector( 0, -1,  1),
	Vector( 0,  1, -1),
	Vector( 0,  1,  1)
};

const std::array<Vector, 12> PolyCubeCut::sm_a_cut_idx2normal{
	Vector( 1, -1,  0),
	Vector( 1,  1,  0),
	Vector( 1,  1,  0),
	Vector( 1, -1,  0),
	Vector(-1,  0,  1),
	Vector( 1,  0,  1),
	Vector( 1,  0,  1),
	Vector(-1,  0,  1),
	Vector( 0, -1,  1),
	Vector( 0,  1,  1),
	Vector( 0,  1,  1),
	Vector( 0, -1,  1)
};

const std::array<Vector, 6> PolyCubeCut::sm_a_label2normal{
    Vector( 1,  0,  0),
    Vector(-1,  0,  0),
    Vector( 0,  1,  0),
    Vector( 0, -1,  0),
    Vector( 0,  0,  1),
    Vector( 0,  0, -1)
};

const std::array<int, 24> PolyCubeCut::sm_a_cutidx2label{
	2, 0,
	3, 0,
	1, 2,
	1, 3,
	0, 4,
	5, 0,
	1, 4,
	5, 1,
	2, 4,
	5, 2,
	3, 4,
	5, 3
};

PolyCubeCut::PolyCubeCut(
    const std::vector<ig::CVec<double, 3>> &v_input_point,
    const std::vector<unsigned int> &v_input_face,
    const std::vector<int> &v_input_chart,
    const std::vector<int> &v_input_label
) : m_polymesh(*new Polyhedron())
{
    try
    {
        build_triangle_mesh(v_input_point, v_input_face);
    }
    catch (const std::string &s)
    {
        throw s;
    }

    for (const auto &vit : vertices(m_polymesh))
        m_v_id2vertex.push_back(vit);
    for (const auto &hit : halfedges(m_polymesh))
        m_v_id2halfedge.push_back(hit);

	m_initial_chart_size = (*std::max_element(v_input_chart.begin(), v_input_chart.end())) + 1;
	m_v_chart2label.resize(m_initial_chart_size);
    m_v_face2label.resize(m_polymesh.size_of_facets());
    m_v_face2chart.resize(m_polymesh.size_of_facets());

    for (unsigned int i = 0; i < m_polymesh.size_of_facets(); ++i)
    {
		m_v_chart2label[v_input_chart[i]] = v_input_label[i];
        m_v_face2label[i] = v_input_label[i];
        m_v_face2chart[i] = v_input_chart[i];
    }

    find_polycube_edge();

	// compute min_edge_length and threshold
	double min_edge_length_2 = DBL_MAX;
	for (Edge_iterator eit = m_polymesh.edges_begin(); eit != m_polymesh.edges_end(); ++eit)
	{
		double temp = (eit->vertex()->point() - eit->opposite()->vertex()->point()).squared_length();
		if (temp < min_edge_length_2)
			min_edge_length_2 = temp;
	}
	m_threshold = std::max(1e-6, min_edge_length_2 / 16);
}

void PolyCubeCut::get_cuttable_edge(std::vector<unsigned int> &v_cuttable_edge)
{
    v_cuttable_edge.clear();

    for (const auto &hit : m_v_polycube_edge)
    {
        m_ms_associated_edge[hit].insert(hit);

        std::array<int, 2> a_adjacent_label{ m_v_face2label[hit->face()->id()], m_v_face2label[hit->opposite()->face()->id()] };
        std::array<Vector, 2> a_adjacent_normal{ sm_a_label2normal[a_adjacent_label[0]], sm_a_label2normal[a_adjacent_label[1]] };

        Vector cut_dir;
        if (a_adjacent_label[0] < a_adjacent_label[1])
            cut_dir = sm_a_cut_idx2dir[sm_am_label2cutidx[a_adjacent_label[0]].at(a_adjacent_label[1])];
        else
            cut_dir = sm_a_cut_idx2dir[sm_am_label2cutidx[a_adjacent_label[1]].at(a_adjacent_label[0])];

        Vector chart_dir = (hit->next()->vertex()->point() - hit->vertex()->point()) * a_adjacent_normal[1] > 0 ? a_adjacent_normal[1] : -a_adjacent_normal[1];

        if (chart_dir * cut_dir > 0)
            m_map_is_convex_cut[hit] = true;
        else
            m_map_is_convex_cut[hit] = false;

        std::array<bool, 2> end_point_on_surface{ is_on_surface(hit, cut_dir), is_on_surface(hit->opposite(), cut_dir) };
        if (end_point_on_surface[0] && end_point_on_surface[1])
        {
            v_cuttable_edge.push_back(static_cast<int>(hit->vertex()->id()));
            v_cuttable_edge.push_back(static_cast<int>(hit->opposite()->vertex()->id()));
        }
        else
        {
            if (!end_point_on_surface[0])
            {
                for (int i = 0; i < 3; ++i)
                    if (m_va_vertex2_a_polycube_edge[hit->vertex()->id()][i] != hit)
                        m_ms_associated_edge[hit].insert(halfedge2edge(m_va_vertex2_a_polycube_edge[hit->vertex()->id()][i]));
            }
            if (!end_point_on_surface[1])
            {
                for (int i = 0; i < 3; ++i)
                    if (m_va_vertex2_a_polycube_edge[hit->opposite()->vertex()->id()][i] != hit->opposite())
                        m_ms_associated_edge[hit].insert(halfedge2edge(m_va_vertex2_a_polycube_edge[hit->opposite()->vertex()->id()][i]));
            }
        }
    }

    bool edited;
    do
    {
        edited = false;
        for (auto &pair : m_ms_associated_edge)
        {
            size_t old_size = pair.second.size();
            for (const auto &hit : pair.second)
                pair.second.insert(m_ms_associated_edge[hit].begin(), m_ms_associated_edge[hit].end());
            if (old_size != pair.second.size())
                edited = true;
        }
    } while (edited);

    //std::cout << "Associated Edges:" << m_ms_associated_edge.size() << std::endl;
    //for (const auto &pair : m_ms_associated_edge)
    //{
    //    std::cout << "\nHalfedge #" << pair.first->id() << ": " << pair.first->vertex()->id() << " -> " << pair.first->opposite()->vertex()->id() << ":\n";
    //    for (const auto &hit: pair.second)
    //        std::cout << "#" << hit->id() << ": " << hit->vertex()->id() << " -> " << hit->opposite()->vertex()->id() << "\n";
    //    std::cout << std::endl;
    //}
}

void PolyCubeCut::get_associated_edge(
	std::map<
		std::pair<unsigned int, unsigned int>,
		std::vector<std::pair<unsigned int, unsigned int>>
	> &map_associated_edge // first id is less than the second
)
{
	for (const auto &pc_pair : m_ms_associated_edge)
	{
		std::pair<unsigned int, unsigned int> edge_key;
		if (pc_pair.first->vertex()->id() < pc_pair.first->opposite()->vertex()->id())
		{
			edge_key.first = static_cast<unsigned int>(pc_pair.first->vertex()->id());
			edge_key.second = static_cast<unsigned int>(pc_pair.first->opposite()->vertex()->id());
		}
		else
		{
			edge_key.first = static_cast<unsigned int>(pc_pair.first->opposite()->vertex()->id());
			edge_key.second = static_cast<unsigned int>(pc_pair.first->vertex()->id());
		}

		if (pc_pair.second.size() != 1)
			for (const auto &hit : m_ms_associated_edge.at(pc_pair.first))
			{
				if (hit == pc_pair.first)
					continue;
				if (hit->vertex()->id() < hit->opposite()->vertex()->id())
					map_associated_edge[edge_key].push_back(std::make_pair(hit->vertex()->id(), hit->opposite()->vertex()->id()));
				else
					map_associated_edge[edge_key].push_back(std::make_pair(hit->opposite()->vertex()->id(), hit->vertex()->id()));
			}
		else
			map_associated_edge[edge_key] = std::vector<std::pair<unsigned int, unsigned int>>();
	}
}

void PolyCubeCut::get_convex_map(std::map<std::pair<unsigned int, unsigned int>, bool> &map_is_convex_edge)
{

	if (m_map_is_convex_cut.empty())
		get_cuttable_edge(std::vector<unsigned int>());

	for (const auto &hit : m_v_polycube_edge)
	{
		std::pair<unsigned int, unsigned int> edge_key;
		if (hit->vertex()->id() < hit->opposite()->vertex()->id())
		{
			edge_key.first = static_cast<unsigned int>(hit->vertex()->id());
			edge_key.second = static_cast<unsigned int>(hit->opposite()->vertex()->id());
		}
		else
		{
			edge_key.first = static_cast<unsigned int>(hit->opposite()->vertex()->id());
			edge_key.second = static_cast<unsigned int>(hit->vertex()->id());
		}
		map_is_convex_edge[edge_key] = m_map_is_convex_cut[hit];
	}
}

// return true if cut a and b is equivalent
bool PolyCubeCut::is_equivalent_cut(
    const std::vector<unsigned int> &v_cut_edge_a,
    const std::vector<double> &v_cut_depth_a,
    const std::vector<unsigned int> &v_cut_edge_b,
    const std::vector<double> &v_cut_depth_b
)
{
    std::vector<Halfedge_handle> v_cut_halfedge_a, v_cut_halfedge_b;
    index2halfedge(v_cut_edge_a, &v_cut_halfedge_a);
    index2halfedge(v_cut_edge_b, &v_cut_halfedge_b);

    std::vector<bool> v_is_cut_edge(m_polymesh.size_of_halfedges(), false);
    std::set<unsigned int> s_complete_cut_edge_a, s_complete_cut_edge_b;
    std::map<unsigned int, double> map_cut2depth_a, map_cut2depth_b;

    get_complete_cut_edge(v_cut_halfedge_a, v_cut_depth_a, map_cut2depth_a, v_is_cut_edge, s_complete_cut_edge_a);
    
    v_is_cut_edge.assign(m_polymesh.size_of_halfedges(), false);
    get_complete_cut_edge(v_cut_halfedge_b, v_cut_depth_b, map_cut2depth_b, v_is_cut_edge, s_complete_cut_edge_b);

    if (s_complete_cut_edge_a == s_complete_cut_edge_b) // same cut edge
    {
        for (const auto &pair : map_cut2depth_a)
            if (fabs(pair.second - map_cut2depth_b[pair.first]) >= EPSILON)
                return false;
        return true;
    }
    else
        return false;
}

void PolyCubeCut::cut_test(
	const std::vector<unsigned int> &v_cut_edge,
	const std::vector<double> &v_cut_depth,
	const double &max_volume,
	const double &tet_quality
)
{
	m_v_chart2label.resize(m_initial_chart_size);
	std::vector<Halfedge_handle> v_cut_halfedge;
	std::vector<int> v_cut_type;
	std::vector<int> v_cut_section_chart;

	try
	{
		index2halfedge(v_cut_edge, &v_cut_halfedge);

		tetgenio tetin, tetout;

		m_v_cut_facet.clear();

		std::map<int, ThreeCut> map_vertex2three_cut;
		std::set<unsigned int> s_complete_cut_edge;
		std::map<unsigned int, double> map_cut2depth;
		std::vector<bool> v_is_cut_edge(m_polymesh.size_of_halfedges(), false);
		std::map<Halfedge_handle, int> map_cut_edge2index;

		get_complete_cut_edge(v_cut_halfedge, v_cut_depth, map_cut2depth, v_is_cut_edge, s_complete_cut_edge);

		// compute cut types
		int i = 0;
		for (const auto &hid : s_complete_cut_edge)
		{
			std::array<int, 2> a_adjacent_label{ m_v_face2label[m_v_id2halfedge[hid]->face()->id()], m_v_face2label[m_v_id2halfedge[hid]->opposite()->face()->id()] };
			if (a_adjacent_label[0] < a_adjacent_label[1])
				v_cut_type.push_back(sm_am_label2cutidx[a_adjacent_label[0]].at(a_adjacent_label[1]));
			else
				v_cut_type.push_back(sm_am_label2cutidx[a_adjacent_label[1]].at(a_adjacent_label[0]));
			map_cut_edge2index[m_v_id2halfedge[hid]] = i++;
		}

		try
		{
			construct_PLC(s_complete_cut_edge, map_cut2depth, v_cut_type, v_is_cut_edge, &v_cut_section_chart, &map_vertex2three_cut, &tetin);
		}
		catch (const std::string &s)
		{
			throw s;
		}

#ifdef _DEBUG
		for (int i = 0; i < m_v_cut_facet.size(); ++i)
		{
			std::cout << "[Info] Cut facet #" << i << ": ";
			for (int j = 0; j < 4; ++j)
				std::cout << m_v_cut_facet[i].a_vertex[j] << " ";
			std::cout << std::endl;
		}
		std::cout << "-------------------------------------------\n";
#endif

		std::cout << "[Info] Tetrahedralizing...\n-------------------------------------------\n";

		char switches[256] = "pfnnQT1e-14";
		char tmpchar[256];

		if (tet_quality != 0)
		{
			sprintf_s(tmpchar, "q%f", tet_quality);
			strcat_s(switches, tmpchar);
		}
		if (max_volume != 0)
		{
			sprintf_s(tmpchar, "a%e", max_volume);
			strcat_s(switches, tmpchar);
		}

#ifdef _DEBUG
		// for debug
		tetin.save_nodes("test");
		tetin.save_poly("test");

		//time_t now = time(0);
		//std::ofstream output(std::to_string(now) + ".off");
		//tetgenio::facet *f;
		//tetgenio::polygon *p;
		//int n_face = 0;

		//for (int i = 0; i < tetin.numberoffacets; ++i)
		//{
		//	f = &tetin.facetlist[i];
		//	for (int j = 0; j < f->numberofpolygons; ++j)
		//	{
		//		p = &f->polygonlist[j];
		//		//if (p->numberofvertices < 3)
		//			//continue;
		//		++n_face;
		//	}
		//}

		//output << "OFF\n" << tetin.numberofpoints << " " << n_face << " 0\n";
		//output << std::setiosflags(std::ios::fixed) << std::setprecision(14);
		//for (int i = 0; i < tetin.numberofpoints; ++i)
		//	output << tetin.pointlist[3 * i] << " " << tetin.pointlist[3 * i + 1] << " " << tetin.pointlist[3 * i + 2] << std::endl;
		//
		//for (int i = 0; i < tetin.numberoffacets; ++i)
		//{
		//	f = &tetin.facetlist[i];
		//	for (int j = 0; j < f->numberofpolygons; ++j)
		//	{
		//		p = &f->polygonlist[j];
		//		//if (p->numberofvertices < 3)
		//			//continue;
		//		output << p->numberofvertices;
		//		for (int k = 0; k < p->numberofvertices; ++k)
		//			output << " " << p->vertexlist[k];
		//		output << std::endl;
		//	}
		//}
		//output.close();
#endif

		try
		{
			tetrahedralize(switches, &tetin, &tetout);
		}
		catch (const int &x)
		{
			std::string err_msg = "[Tetgen error] ";
			switch (x)
			{
			case 1: // Out of memory.
				err_msg += "Error:  Out of memory.\n";
				break;
			case 2: // Encounter an internal error.
				err_msg += "Please report this bug to Hang.Si@wias-berlin.de. Include ";
				err_msg += "the message above, your input data set, and the exact ";
				err_msg += "command line you used to run this program, thank you.\n";
				break;
			case 3:
				err_msg += "A self-intersection was detected. Program stopped. ";
				err_msg += "Hint: use -d option to detect all self-intersections.\n";
				break;
			case 4:
				err_msg += "A very small input feature size was detected. Program stopped.\n";
				break;
			case 5:
				err_msg += "Two very close input facets were detected. Program stopped. ";
				err_msg += "Hint: use -Y option to avoid adding Steiner points in boundary.\n";
				break;
			case 10:
				err_msg += "An input error was detected. Program stopped.\n";
				break;
			default:
				err_msg += "Error code: " + std::to_string(x);
				break;
			}
			throw err_msg;
		}
		catch (const std::string &s)
		{
			throw s;
		}
	}
	catch (const std::string &s)
	{
		throw s;
	}
}

void PolyCubeCut::cut(
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
	m_v_chart2label.resize(m_initial_chart_size);

    std::vector<Halfedge_handle> v_cut_halfedge;

    try
    {
        index2halfedge(v_cut_edge, &v_cut_halfedge);

        cut(
            v_cut_halfedge, v_cut_depth, max_volume, tet_quality,
            v_output_point, v_output_tet, v_output_chart, v_output_label,
            v_cut_type, v_cut_section_chart, v_seam_edge_vertex_id, vv_congruent_face, v_three_cut
        );
    }
    catch (const std::string &s)
    {
        throw s;
    }

    m_vv_congruent_face = vv_congruent_face;

	//output_vtk(v_output_point, v_output_tet, "test_out.vtk");
}

void PolyCubeCut::index2halfedge(const std::vector<unsigned int> &v_vertex_index, std::vector<Halfedge_handle> *v_halfedge)
{
    assert((v_vertex_index.size() & 1) == 0);
    int n = static_cast<int>(v_vertex_index.size()) >> 1;
    for (int k = 0; k < n; ++k)
    {
        auto find_he = halfedge(
            m_v_id2vertex[v_vertex_index[k << 1]],
            m_v_id2vertex[v_vertex_index[(k << 1) + 1]],
            m_polymesh
        ); // find halfedge va -> vb
        
        Halfedge_handle he = find_he.first;
        if (he == NULL || !find_he.second)
        {
            std::string err_msg = "[Error] Cannot find halfedge: " +
                std::to_string(v_vertex_index[k << 1]) + " -> " +
                std::to_string(v_vertex_index[(k << 1) + 1]) + "\n";
            throw err_msg;
        }

        v_halfedge->push_back(he);
    }
}

void PolyCubeCut::get_complete_cut_edge(
    const std::vector<Halfedge_handle> &v_cut_edge,
    const std::vector<double> &v_cut_depth,
    std::map<unsigned int, double> &map_cut2depth,
    std::vector<bool> &v_is_cut_edge,
    std::set<unsigned int> &s_complete_cut_edge
)
{
    for (int k = 0; k < v_cut_edge.size(); ++k)
    {
        for (const auto &hit : m_ms_associated_edge[halfedge2edge(v_cut_edge[k])])
        {
            s_complete_cut_edge.insert(static_cast<unsigned int>(hit->id()));
            map_cut2depth[static_cast<unsigned int>(hit->id())] = v_cut_depth[k];
            v_is_cut_edge[hit->id()] = true;
            v_is_cut_edge[hit->opposite()->id()] = true;
        }
    }

    for (const auto &vit : vertices(m_polymesh))
        check_three_cut(vit, v_is_cut_edge, s_complete_cut_edge, map_cut2depth);

    adjust_cut_depth(v_is_cut_edge, map_cut2depth);
}

int PolyCubeCut::vertex_cut_count(
    const std::vector<bool> &v_is_cut_edge,
    const std::map<unsigned int, double> &map_cut2depth,
    const Vertex_handle &v,
    Halfedge_handle &new_cut_edge,
    double &depth
)
{
    int n_vertex_cut = 0;
    for (int i = 0; i < 3; ++i)
    {
        if (!v_is_cut_edge[m_va_vertex2_a_polycube_edge[v->id()][i]->id()])
        {
            new_cut_edge = m_va_vertex2_a_polycube_edge[v->id()][i];
        }
        else
        {
            depth = map_cut2depth.at(static_cast<unsigned int>(halfedge2edge(m_va_vertex2_a_polycube_edge[v->id()][i])->id()));
            ++n_vertex_cut;
        }
    }
    return n_vertex_cut;
}

void PolyCubeCut::check_three_cut(
    const Vertex_handle &v,
    std::vector<bool> &v_is_cut_edge,
    std::set<unsigned int> &s_complete_cut_edge,
    std::map<unsigned int, double> &map_cut2depth
)
{
    Halfedge_handle new_cut_edge;
    double depth;

    if (vertex_cut_count(v_is_cut_edge, map_cut2depth, v, new_cut_edge, depth) == 2)
    {
        for (const auto &hit : m_ms_associated_edge[halfedge2edge(new_cut_edge)])
        {
            if (!v_is_cut_edge[hit->id()])
            {
                s_complete_cut_edge.insert(static_cast<unsigned int>(hit->id()));
                map_cut2depth[static_cast<unsigned int>(hit->id())] = depth;
                v_is_cut_edge[hit->id()] = true;
                v_is_cut_edge[hit->opposite()->id()] = true;
                check_three_cut(hit->vertex(), v_is_cut_edge, s_complete_cut_edge, map_cut2depth);
                check_three_cut(hit->opposite()->vertex(), v_is_cut_edge, s_complete_cut_edge, map_cut2depth);
            }
        }
    }
}

void PolyCubeCut::adjust_cut_depth(
    const std::vector<bool> &v_is_cut_edge,
    std::map<unsigned int, double> &map_cut2depth
)
{
    std::vector<bool> v_visited(m_polymesh.size_of_vertices(), false);

    for (const auto &vit : vertices(m_polymesh))
    {
        std::queue<Vertex_handle> q_vertex;
        std::set<Halfedge_handle> s_joint_cut;
        double min_depth = DBL_MAX;

        if (!v_visited[vit->id()])
            q_vertex.push(vit);

        while (!q_vertex.empty())
        {
            Vertex_handle v = q_vertex.front();
            q_vertex.pop();
            v_visited[v->id()] = true;

            Halfedge_handle new_cut_edge;
            double depth;

            if (vertex_cut_count(v_is_cut_edge, map_cut2depth, v, new_cut_edge, depth) == 3) // current vertex is 3-cut
            {
                for (int i = 0; i < 3; ++i)
                {
                    Halfedge_handle he = m_va_vertex2_a_polycube_edge[v->id()][i];
                    s_joint_cut.insert(halfedge2edge(he));
                    if (map_cut2depth.at(static_cast<unsigned int>(halfedge2edge(he)->id())) < min_depth)
                        min_depth = map_cut2depth.at(static_cast<unsigned int>(halfedge2edge(he)->id()));
                    if (!v_visited[he->opposite()->vertex()->id()])
                        q_vertex.push(he->opposite()->vertex());
                }
            }
        }

        for (const auto &hit : s_joint_cut)
            map_cut2depth[static_cast<unsigned int>(hit->id())] = min_depth;
    }
}

bool vertex_compare(
    const std::pair<double, std::array<unsigned int, 3>> &a,
    const std::pair<double, std::array<unsigned int, 3>> &b
)
{
    return a.first > b.first;
}

//void call_tetgen(std::string str_switches, tetgenio &in, tetgenio &out)
//{
//	char switches[256];
//	memcpy(switches, str_switches.c_str(), 256);
//	tetrahedralize(switches, &in, &out);
//}

void PolyCubeCut::cut(
    const std::vector<Halfedge_handle> &v_cut_edge,
    const std::vector<double> &v_cut_depth,
    const double  &max_volume,
    const double  &tet_quality,
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
    v_output_point.clear();
    v_output_tet.clear();
    v_output_chart.clear();
    v_output_label.clear();
    v_cut_type.clear();
    v_cut_section_chart.clear();
    v_seam_edge_vertex_id.clear();
    vv_congruent_face.clear();
    v_three_cut.clear();

    tetgenio tetin, tetout;

    m_v_cut_facet.clear();

    std::map<int, ThreeCut> map_vertex2three_cut;
    std::set<unsigned int> s_complete_cut_edge;
    std::map<unsigned int, double> map_cut2depth;
    std::vector<bool> v_is_cut_edge(m_polymesh.size_of_halfedges(), false);
    std::map<Halfedge_handle, int> map_cut_edge2index;

    get_complete_cut_edge(v_cut_edge, v_cut_depth, map_cut2depth, v_is_cut_edge, s_complete_cut_edge);

    // compute cut types
    int i = 0;
    for (const auto &hid : s_complete_cut_edge)
    {
        std::array<int, 2> a_adjacent_label{ m_v_face2label[m_v_id2halfedge[hid]->face()->id()], m_v_face2label[m_v_id2halfedge[hid]->opposite()->face()->id()] };
        if (a_adjacent_label[0] < a_adjacent_label[1])
            v_cut_type.push_back(sm_am_label2cutidx[a_adjacent_label[0]].at(a_adjacent_label[1]));
        else
            v_cut_type.push_back(sm_am_label2cutidx[a_adjacent_label[1]].at(a_adjacent_label[0]));
        map_cut_edge2index[m_v_id2halfedge[hid]] = i++;
    }

    try
    {
        construct_PLC(s_complete_cut_edge, map_cut2depth, v_cut_type, v_is_cut_edge, &v_cut_section_chart, &map_vertex2three_cut, &tetin);
    }
    catch (const std::string &s)
    {
        throw s;
    }

#ifdef _DEBUG
    for (int i = 0; i < m_v_cut_facet.size(); ++i)
    {
        std::cout << "[Info] Cut facet #" << i << ": ";
        for (int j = 0; j < 4; ++j)
            std::cout << m_v_cut_facet[i].a_vertex[j] << " ";
        std::cout << std::endl;
    }
    std::cout << "-------------------------------------------\n";
#endif

    std::cout << "[Info] Tetrahedralizing...\n-------------------------------------------\n";

    char switches[256] = "pfnnQT1e-14";
    char tmpchar[256];

    if (tet_quality != 0)
    {
        sprintf_s(tmpchar, "q%f", tet_quality);
        strcat_s(switches, tmpchar);
    }
    if (max_volume != 0)
    {
        sprintf_s(tmpchar, "a%e", max_volume);
        strcat_s(switches, tmpchar);
    }

#ifdef _DEBUG
    // for debug
    tetin.save_nodes("test");
    tetin.save_poly("test");
#endif

	try
	{
		//std::string str_switches(switches);
		//std::thread tetgen_thread(call_tetgen, std::ref(str_switches), std::ref(tetin), std::ref(tetout));

		//std::mutex mutex_wait;
		//std::condition_variable cond_wait;
		//std::unique_lock<std::mutex> lock_wait(mutex_wait);
		//std::cv_status cv_wait = cond_wait.wait_for(lock_wait, std::chrono::seconds(120));
		//if (cv_wait == std::cv_status::timeout)
		//{
		//	tetgen_thread.join();
		//	throw std::string("[Error] A timeout occurred while tetrahedralizing.\n");
		//}
		//tetgen_thread.join();
		tetrahedralize(switches, &tetin, &tetout);
    }
    catch (const int &x)
    {
        std::string err_msg = "[Tetgen error] ";
        switch (x)
        {
        case 1: // Out of memory.
            err_msg += "Error: Out of memory.\n";
            break;
        case 2: // Encounter an internal error.
            err_msg += "Please report this bug to Hang.Si@wias-berlin.de. Include ";
            err_msg += "the message above, your input data set, and the exact ";
            err_msg += "command line you used to run this program, thank you.\n";
            break;
        case 3:
            err_msg += "A self-intersection was detected. Program stopped. ";
            err_msg += "Hint: use -d option to detect all self-intersections.\n";
            break;
        case 4:
            err_msg += "A very small input feature size was detected. Program stopped.\n";
            break;
        case 5:
            err_msg += "Two very close input facets were detected. Program stopped. ";
            err_msg += "Hint: use -Y option to avoid adding Steiner points in boundary.\n";
            break;
        case 10:
            err_msg += "An input error was detected. Program stopped.\n";
            break;
        default:
            err_msg += "Error code: " + std::to_string(x);
            break;
        }
        throw err_msg;
    }
	catch (const std::string &s)
	{
		throw s;
	}

    // for debug
    //PolyCubeCutHelper::tetgen2vtk(&tetout, "test-before.vtk");
    //std::cout << "\nPoint:\n";
    //for (int i = 0; i < tetout.numberofpoints; ++i)
    //    std::cout << i << "\t" << tetout.pointlist[3 * i] << "\t" << tetout.pointlist[3 * i + 1] << "\t" << tetout.pointlist[3 * i + 2] << std::endl;
    //std::cout << "\nFacet:\n";
    //for (int i = 0; i < tetout.numberoftrifaces; ++i)
    //{
    //    std::cout << "face #" << i << "\tmarker: " << tetout.trifacemarkerlist[i];
    //    std::cout << std::endl;
    //    for (int j = 0; j < 3; ++j)
    //    {
    //        std::cout << tetout.trifacelist[3 * i + j] << "\t";
    //    }
    //    std::cout << std::endl;
    //}
    //std::cout << "\nTetrahedron:\n";
    //for (int i = 0; i < tetout.numberoftetrahedra; ++i)
    //{
    //    std::cout << "Tet #" << i << std::endl;
    //    for (int j = 0; j < 4; ++j)
    //    {
    //        std::cout << tetout.tetrahedronlist[4 * i + j] << "\t";
    //    }
    //    std::cout << std::endl;
    //}
#ifdef _DEBUG
// for debug
	tetout.save_nodes("tetout");
#endif

    // Before subdivision, record the number of faces that have boundary markers.
    // Boundary markers of newly added faces are all 0 (not on the boundary)
    int marker_max_index = tetout.numberoftrifaces;

    std::cout << "[Info] Tetrahedron subdividing...\n-------------------------------------------\n";

    tetrahedron_subdivision(&tetout);

    v_output_chart.assign(tetout.numberoftetrahedra, -1);
    v_output_label.resize(tetout.numberoftetrahedra);

    std::cout << "[Info] Chart updating...\n-------------------------------------------\n";

    update_boundary_tet_chart(marker_max_index, &tetout, &v_output_chart);

    size_t n_cut = v_cut_type.size();
    vv_congruent_face.resize(n_cut);

    std::cout << "[Info] Cut facet separating...\n-------------------------------------------\n";

    separate_cut_facet(v_cut_type, v_cut_section_chart, &v_output_chart, &vv_congruent_face, &map_vertex2three_cut, &tetout);

    // for debug
	//PolyCubeCutHelper::tetgen2vtk(&tetout, "test-after.vtk");
    //std::cout << "\nPoint:\n";
    //for (int i = 0; i < tetout.numberofpoints; ++i)
    //    std::cout << i << "\t" << tetout.pointlist[3 * i] << "\t" << tetout.pointlist[3 * i + 1] << "\t" << tetout.pointlist[3 * i + 2] << std::endl;
    //std::cout << "\nTetrahedron:\n";
    //for (int i = 0; i < tetout.numberoftetrahedra; ++i)
    //{
    //    std::cout << "Tet #" << i << std::endl;
    //    for (int j = 0; j < 4; ++j)
    //    {
    //        std::cout << tetout.tetrahedronlist[4 * i + j] << "\t";
    //    }
    //    std::cout << std::endl;
    //}
    //std::cout << "\nFacet:\n";
    //for (int i = 0; i < tetout.numberoftrifaces; ++i)
    //{
    //    std::cout << "face #" << i;
    //    if (i < marker_max_index)
    //        std::cout << "\tmarker: " << tetout.trifacemarkerlist[i];
    //    std::cout << std::endl;
    //    for (int j = 0; j < 3; ++j)
    //    {
    //        std::cout << tetout.trifacelist[3 * i + j] << "\t";
    //    }
    //    std::cout << std::endl;
    //}

    for (int i = 0; i < tetout.numberofpoints; ++i)
        v_output_point.push_back(
            ig::CVec<double, 3>(
                tetout.pointlist[3 * i],
                tetout.pointlist[3 * i + 1],
                tetout.pointlist[3 * i + 2]
                )
        );

    for (int i = 0; i < tetout.numberoftetrahedra; ++i)
    {
        for (int j = 0; j < 4; ++j)
            v_output_tet.push_back(tetout.tetrahedronlist[4 * i + j]);

        if (v_output_chart[i] >= 0)
            v_output_label[i] = m_v_chart2label[v_output_chart[i]];
        else
            v_output_label[i] = -1;
    }

    v_seam_edge_vertex_id.reserve(n_cut << 1);
    for (int k = 0; k < n_cut; ++k)
    {
        v_seam_edge_vertex_id.push_back(m_v_cut_facet[k].a_vertex[2]);
        v_seam_edge_vertex_id.push_back(m_v_cut_facet[k].a_vertex[3]);
    }

    for (const auto &tc : map_vertex2three_cut)
    {
        v_three_cut.push_back(tc.second);

        // sort corresponding vertices
        ThreeCut &three_cut = v_three_cut.back();

        Point shared_point(
            tetout.pointlist[3 * three_cut.shared_vertex_id],
            tetout.pointlist[3 * three_cut.shared_vertex_id + 1],
            tetout.pointlist[3 * three_cut.shared_vertex_id + 2]
        );
        
        std::vector<std::pair<double, std::array<unsigned int, 3>>> vp_corresponding_vertex_with_dis;
        for (const auto &cor_vert : three_cut.va_corresponding_vertex)
        {
            Point point(
                tetout.pointlist[3 * cor_vert[2]],
                tetout.pointlist[3 * cor_vert[2] + 1],
                tetout.pointlist[3 * cor_vert[2] + 2]
            );
            vp_corresponding_vertex_with_dis.push_back(std::make_pair((point - shared_point).squared_length(), cor_vert));
        }

        std::sort(vp_corresponding_vertex_with_dis.begin(), vp_corresponding_vertex_with_dis.end(), vertex_compare);

        for (int i = 0; i < vp_corresponding_vertex_with_dis.size(); ++i)
            three_cut.va_corresponding_vertex[i] = vp_corresponding_vertex_with_dis[i].second;

        three_cut.a_adjacent_one_cut_index[0] = map_cut_edge2index.at(halfedge2edge(m_va_vertex2_a_polycube_edge[tc.second.va_corresponding_vertex[0][2]][1]));
        three_cut.a_adjacent_one_cut_index[1] = map_cut_edge2index.at(halfedge2edge(m_va_vertex2_a_polycube_edge[tc.second.va_corresponding_vertex[0][2]][0]));
        three_cut.a_adjacent_one_cut_index[2] = map_cut_edge2index.at(halfedge2edge(m_va_vertex2_a_polycube_edge[tc.second.va_corresponding_vertex[0][2]][2]));
    }
}

void PolyCubeCut::get_cut_edge_map(std::map<std::pair<unsigned int, unsigned int>, std::vector<std::pair<unsigned int, unsigned int>>>& um_cut_edge_map)
{
    um_cut_edge_map.clear();

    std::unordered_set<int> us_cut_corner; // corners of all cut edges

    for (const auto& cf : m_v_cut_facet)
    {
        us_cut_corner.insert(cf.a_vertex[0]);
        us_cut_corner.insert(cf.a_vertex[1]);
        us_cut_corner.insert(m_v_duplicated_vertex_map.at(cf.a_vertex[0])[0]);
        us_cut_corner.insert(m_v_duplicated_vertex_map.at(cf.a_vertex[1])[0]);
        if (m_v_duplicated_vertex_map.at(cf.a_vertex[0])[1] != -1)
            us_cut_corner.insert(m_v_duplicated_vertex_map.at(cf.a_vertex[0])[1]);
        if (m_v_duplicated_vertex_map.at(cf.a_vertex[1])[1] != -1)
            us_cut_corner.insert(m_v_duplicated_vertex_map.at(cf.a_vertex[1])[1]);
    }

    for (const auto& v_congruent_face : m_vv_congruent_face)
    {
        std::vector<std::pair<int, int>> vp_corner_pair;
        int first = -1;
        for (int i = 0; i < v_congruent_face.size(); ++i)
        {
            if (i != 0 && i % 3 == 0) i += 3;
            if (i >= v_congruent_face.size()) break;
            if (us_cut_corner.find(v_congruent_face[i]) != us_cut_corner.end() && (first == -1 || (first != -1 && v_congruent_face[i] != first)))
            {
                vp_corner_pair.push_back(std::make_pair(v_congruent_face[i], v_congruent_face[i + 3]));
                first = v_congruent_face[i];
            }
        }

        // find the original indices of two corners
        int o0 = vp_corner_pair[0].first > m_v_duplicated_vertex_map[vp_corner_pair[0].first][0] ? m_v_duplicated_vertex_map[vp_corner_pair[0].first][0] : vp_corner_pair[0].first;
        int o1 = vp_corner_pair[1].first > m_v_duplicated_vertex_map[vp_corner_pair[1].first][0] ? m_v_duplicated_vertex_map[vp_corner_pair[1].first][0] : vp_corner_pair[1].first;

        std::vector<std::pair<unsigned int, unsigned int>> vp_new_edge;
        if (o0 < o1)
        {
            vp_new_edge.push_back(std::make_pair(vp_corner_pair[0].first, vp_corner_pair[1].first));
            vp_new_edge.push_back(std::make_pair(vp_corner_pair[0].second, vp_corner_pair[1].second));
            um_cut_edge_map[std::make_pair(o0, o1)] = vp_new_edge;
        }
        else
        {
            vp_new_edge.push_back(std::make_pair(vp_corner_pair[1].first, vp_corner_pair[0].first));
            vp_new_edge.push_back(std::make_pair(vp_corner_pair[1].second, vp_corner_pair[0].second));
            um_cut_edge_map[std::make_pair(o1, o0)] = vp_new_edge;
        }
    }
}

void PolyCubeCut::build_triangle_mesh(const std::vector<ig::CVec<double, 3>>& v_input_point, const std::vector<unsigned int>& v_input_face)
{
    typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
    typedef typename HalfedgeDS::Vertex Vertex;

    CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> builder(m_polymesh.hds());
    size_t n_face = v_input_face.size() / 3;
    builder.begin_surface(v_input_point.size(), n_face);

    for (const auto &p : v_input_point)
        builder.add_vertex(Vertex::Point(p[0], p[1], p[2]));

    for (size_t i = 0; i < n_face; ++i)
    {
        builder.begin_facet();
        for (int j = 0; j < 3; ++j)
            builder.add_vertex_to_facet(v_input_face[3 * i + j]);
        builder.end_facet();
    }

    if (builder.error())
        throw std::string("[Error] Error occurred when building the surface mesh.\n");

    builder.end_surface();
    
    if (!m_polymesh.is_valid())
        throw std::string("[Error] Invalid input surface mesh.\n");
    if (!CGAL::is_triangle_mesh(m_polymesh))
        throw std::string("[Error] Input is not a triangle mesh.\n");

    PolyCubeCutHelper::assign_cgal_id(&m_polymesh);
}

void PolyCubeCut::find_polycube_edge()
{
    m_v_polycube_edge.clear();
    for (Edge_iterator eit = m_polymesh.edges_begin(); eit != m_polymesh.edges_end(); ++eit)
        if (m_v_face2chart[eit->face()->id()] != m_v_face2chart[eit->opposite()->face()->id()])
            m_v_polycube_edge.push_back(eit);

    m_va_vertex2_a_polycube_edge.resize(m_polymesh.size_of_vertices());
    for (auto const &vit : vertices(m_polymesh))
    {
        // collect all 3 polycube incident halfedges with respect to this vertex (clockwise)
        Halfedge_handle hit = vit->halfedge();
        for (int i = 0; i < 3;)
        {
            if (m_v_face2chart[hit->face()->id()] != m_v_face2chart[hit->opposite()->face()->id()])
                m_va_vertex2_a_polycube_edge[vit->id()][i++] = hit;
            hit = hit->next_on_vertex();
        }
    }
}

bool PolyCubeCut::is_on_surface(const Halfedge_handle &he, const Vector &cut_dir)
{
    Halfedge_handle in_he, out_he; // on orthogonal chart
    for (int i = 0; i < 3; ++i)
        if (m_va_vertex2_a_polycube_edge[he->vertex()->id()][i] == he)
        {
            in_he = m_va_vertex2_a_polycube_edge[he->vertex()->id()][(i + 1) % 3];
            out_he = m_va_vertex2_a_polycube_edge[he->vertex()->id()][(i + 2) % 3]->opposite();
        }
    int current_chart = m_v_face2chart[in_he->face()->id()];
    int current_label = m_v_chart2label[current_chart];
    Vector in_vec(in_he->vertex()->point() - in_he->opposite()->vertex()->point());
    Vector out_vec(out_he->vertex()->point() - out_he->opposite()->vertex()->point());
    if (CGAL::cross_product(in_vec, out_vec) * sm_a_label2normal[current_label] > 0) // convex point
    {
        if (-in_vec * cut_dir > 0 && out_vec * cut_dir > 0)
            return true;
        else
            return false;
    }
    else
    {
        if (-in_vec * cut_dir < 0 && out_vec * cut_dir < 0)
            return true;
        else
            return false;
    }
}

Halfedge_handle PolyCubeCut::halfedge2edge(const Halfedge_handle &he)
{
    return he->id() & 1 ? he->opposite() : he;
}

// construct PLC (Piecewise Linear Complexes) in tetin
void PolyCubeCut::construct_PLC(
    const std::set<unsigned int> &s_complete_cut_edge,
    const std::map<unsigned int, double> &map_cut2depth,
    const std::vector<int> &v_cut_type,
    const std::vector<bool> &v_is_cut_edge,
    std::vector<int> *v_cut_section_chart,
    std::map<int, ThreeCut> *map_vertex2three_cut,
    tetgenio *tetin
)
{
    TetgenHelper tetgen_helper;

    int n_cut = static_cast<int>(s_complete_cut_edge.size());

    v_cut_section_chart->reserve(s_complete_cut_edge.size());

    tetgen_helper.pointlist.reserve(3 * (m_polymesh.size_of_vertices() + (n_cut << 1)));
    tetgen_helper.facetlist.reserve(m_polymesh.size_of_facets() + n_cut);

    PolyCubeCutHelper::cgal2tetgen_helper(m_polymesh, m_v_face2chart, &tetgen_helper);

    std::map<size_t, unsigned int> map_vertex2cut_end_vertex_id;
    std::vector<bool> v_is_interior(m_polymesh.size_of_vertices(), false);
    
    auto hid_it = s_complete_cut_edge.begin();
    for (int k = 0; k < n_cut; ++k, ++hid_it)
    {
        int new_chart_idx = static_cast<int>(m_v_chart2label.size());
        v_cut_section_chart->push_back(-999); // placeholder
        v_cut_section_chart->push_back(new_chart_idx);
        v_cut_section_chart->push_back(new_chart_idx + 1);

        try
        {
			add_cut(
				m_v_id2halfedge[*hid_it], map_cut2depth.at(*hid_it), v_cut_type[k], v_is_cut_edge,
				&map_vertex2cut_end_vertex_id, &v_is_interior, map_vertex2three_cut, &tetgen_helper
			);
        }
        catch (const std::string &s)
        {
            throw s;
        }

        Vector next_he_dir(m_v_id2halfedge[*hid_it]->next()->vertex()->point() - m_v_id2halfedge[*hid_it]->vertex()->point());
        if (m_v_cut_facet.back().normal * next_he_dir > 0)
        {
            (*v_cut_section_chart)[v_cut_section_chart->size() - 3] = m_v_face2chart[m_v_id2halfedge[*hid_it]->face()->id()];
            v_cut_section_chart->push_back(m_v_face2chart[m_v_id2halfedge[*hid_it]->opposite()->face()->id()]);
        }
        else
        {
            (*v_cut_section_chart)[v_cut_section_chart->size() - 3] = m_v_face2chart[m_v_id2halfedge[*hid_it]->opposite()->face()->id()];
            v_cut_section_chart->push_back(m_v_face2chart[m_v_id2halfedge[*hid_it]->face()->id()]);
        }

        if (m_v_cut_facet.back().is_convex)
        {
			m_v_chart2label.push_back(sm_a_cutidx2label[v_cut_type[k] << 1]);
			m_v_chart2label.push_back(sm_a_cutidx2label[(v_cut_type[k] << 1) + 1]);
        }
        else
        {
			m_v_chart2label.push_back(opposite_label(sm_a_cutidx2label[(v_cut_type[k] << 1) + 1]));
			m_v_chart2label.push_back(opposite_label(sm_a_cutidx2label[v_cut_type[k] << 1]));
        }
    }

    PolyCubeCutHelper::tetgen_helper2tetgen(tetgen_helper, tetin);
}

// return the label which is opposite to the input
int PolyCubeCut::opposite_label(int label)
{
    return label & 1 ? label - 1 : label + 1;
}

void PolyCubeCut::add_cut(
    const Halfedge_handle &cut_he,
    const double &cut_depth,
    const int &cut_type,
    const std::vector<bool> &v_is_cut_edge,
    std::map<size_t, unsigned int> *map_vertex2cut_end_vertex_id,
    std::vector<bool> *v_is_interior,
    std::map<int, ThreeCut> *map_vertex2three_cut,
    TetgenHelper *tetgen_helper
)
{
	CutFacet cut_facet;
    cut_facet.is_convex = m_map_is_convex_cut[cut_he];
    std::array<int, 2> a_cut_chart = find_cut_chart(cut_he);
    std::array<Vertex_handle, 2> a_cut_start_vertex{
        cut_he->opposite()->vertex(),
        cut_he->vertex()
    };
    std::array<Point, 2> a_cut_end_vertex; // two end vertices of the cut
    std::array<unsigned int, 2> a_cut_end_vertex_id;
    for (int i = 0; i < 2; ++i)
    {
        if (map_vertex2cut_end_vertex_id->find(a_cut_start_vertex[i]->id()) == map_vertex2cut_end_vertex_id->end())
        {
            bool interior = false;
            (*map_vertex2cut_end_vertex_id)[a_cut_start_vertex[i]->id()] = create_cut_end_vertex(
                a_cut_start_vertex[i],
                cut_he,
                cut_depth,
                cut_type,
                v_is_cut_edge,
                map_vertex2three_cut,
                &interior,
                tetgen_helper
            );
            (*v_is_interior)[a_cut_start_vertex[i]->id()] = interior;

            const Point &p_new = tetgen_helper->pointlist[(*map_vertex2cut_end_vertex_id)[a_cut_start_vertex[i]->id()]];
#ifdef _DEBUG
            std::cout << "[Info] New vertex #" << (*map_vertex2cut_end_vertex_id)[a_cut_start_vertex[i]->id()] <<
                " created: " << p_new << "\n-------------------------------------------\n";
#endif
        }
        cut_facet.a_is_interior[i] = (*v_is_interior)[a_cut_start_vertex[i]->id()];
        if (cut_facet.a_is_interior[i] == true)
            a_cut_chart[i] = -1;
        a_cut_end_vertex[i] = tetgen_helper->pointlist[(*map_vertex2cut_end_vertex_id)[a_cut_start_vertex[i]->id()]];
        a_cut_end_vertex_id[i] = (*map_vertex2cut_end_vertex_id)[a_cut_start_vertex[i]->id()];
    }

	// find 4 vertex positions of the current facet
	// last 2 vertices of each facet form a seam edge
	cut_facet.a_vertex[0] = static_cast<unsigned int>(a_cut_start_vertex[0]->id());
	cut_facet.a_vertex[1] = static_cast<unsigned int>(a_cut_start_vertex[1]->id());
	cut_facet.a_vertex[2] = a_cut_end_vertex_id[1];
	cut_facet.a_vertex[3] = a_cut_end_vertex_id[0];
	for (int i = 0; i < 4; ++i)
		cut_facet.a_point[i] = tetgen_helper->pointlist[cut_facet.a_vertex[i]];
	for (int i = 0; i < 4; ++i)
	{
		cut_facet.a_edge_vector[i] = Vector(cut_facet.a_point[(i + 1) % 4] - cut_facet.a_point[i]);
		cut_facet.a_edge_length_2[i] = cut_facet.a_edge_vector[i].squared_length();
	}
	cut_facet.normal = sm_a_cut_idx2normal[cut_type];
	cut_facet.plane = Plane(cut_facet.a_point[0], cut_facet.normal);

    try
    {
        check_intersection(cut_he, a_cut_end_vertex);
		check_cut_spacing(cut_facet);
    }
    catch (const std::string &s)
    {
        throw s;
    }

    do_cut_in_tetgen_helper(a_cut_start_vertex, a_cut_end_vertex_id, a_cut_end_vertex, sm_a_cut_idx2dir[cut_type], a_cut_chart, tetgen_helper);

	m_v_cut_facet.push_back(cut_facet);
}

// Create corresponding end vertex of the start vertex with respect to cut_he in tetgen_helper.
// Return the index of the newly created vertex.
unsigned int PolyCubeCut::create_cut_end_vertex(
    const Vertex_handle &cut_start_vertex,
    const Halfedge_handle &cut_he,
    const double &cut_depth,
    const int &cut_type,
    const std::vector<bool> &v_is_cut_edge,
    std::map<int, ThreeCut> *map_vertex2three_cut,
    bool *interior,
    TetgenHelper *tetgen_helper
)
{
    Halfedge_handle incident_cut_he = cut_he->vertex() == cut_start_vertex ? cut_he : cut_he->opposite();
    Vector cut_dir(0, 0, 0);
    for (const auto &pc_he : m_va_vertex2_a_polycube_edge[cut_start_vertex->id()])
    {
        if (pc_he != incident_cut_he && v_is_cut_edge[pc_he->id()])
            *interior = true;
        cut_dir += sm_a_label2normal[m_v_face2label[pc_he->face()->id()]];
    }

    Point cut_end_vertex;

    if (*interior)
    {
        cut_dir = -cut_dir;
        (*map_vertex2three_cut)[static_cast<int>(cut_start_vertex->id())].cut_type = three_cut_direction2index(cut_dir);
        (*map_vertex2three_cut)[static_cast<int>(cut_start_vertex->id())].shared_vertex_id = static_cast<unsigned int>(tetgen_helper->pointlist.size());
    }
    else
        cut_dir = sm_a_cut_idx2dir[cut_type];

    cut_dir /= CGAL::sqrt(cut_dir.squared_length()); // normalize
    if (*interior)
        cut_end_vertex = cut_start_vertex->point() + cut_dir * cut_depth * 3 / sqrt(6);
    else
        cut_end_vertex = cut_start_vertex->point() + cut_dir * cut_depth;

    tetgen_helper->pointlist.push_back(cut_end_vertex);

    return static_cast<unsigned int>(tetgen_helper->pointlist.size()) - 1;
}

// find the 2 chart that the cut plane intersect with
std::array<int, 2> PolyCubeCut::find_cut_chart(const Halfedge_handle &cut_he)
{
	std::array<int, 2> a_cut_chart{ -1, -1 };
	int cut_he_chart = m_v_face2chart[cut_he->face()->id()];
	
	Halfedge_handle hit = cut_he->prev()->opposite();
	int current_chart = m_v_face2chart[hit->face()->id()];
	while (a_cut_chart[0] == -1)
	{
		if (current_chart != cut_he_chart)
			a_cut_chart[0] = current_chart;
		hit = hit->prev()->opposite();
		current_chart = m_v_face2chart[hit->face()->id()];
	}

	// another side
	hit = cut_he->next()->opposite();
	current_chart = m_v_face2chart[hit->face()->id()];
	while (a_cut_chart[1] == -1)
	{
		if (current_chart != cut_he_chart)
			a_cut_chart[1] = current_chart;
		hit = hit->next()->opposite();
		current_chart = m_v_face2chart[hit->face()->id()];
	}

	return a_cut_chart;
}

int PolyCubeCut::three_cut_direction2index(const Vector &dir)
{
    if (dir.x() > 0 && dir.y() > 0 && dir.z() > 0)
        return 0;
    else if (dir.x() > 0 && dir.y() > 0 && dir.z() < 0)
        return 1;
    else if (dir.x() > 0 && dir.y() < 0 && dir.z() > 0)
        return 2;
    else if (dir.x() > 0 && dir.y() < 0 && dir.z() < 0)
        return 3;
    else if (dir.x() < 0 && dir.y() > 0 && dir.z() > 0)
        return 4;
    else if (dir.x() < 0 && dir.y() > 0 && dir.z() < 0)
        return 5;
    else if (dir.x() < 0 && dir.y() < 0 && dir.z() > 0)
        return 6;
    else if (dir.x() < 0 && dir.y() < 0 && dir.z() < 0)
        return 7;
    assert(false);
    return -1;
}

bool PolyCubeCut::is_facet_too_close(
	const CutFacet& fa,
	const CutFacet& fb
)
{
	for (int i = 2; i < 3; ++i)
	{
		const Point& p = fa.a_point[i];
		Point proj_p = fb.plane.projection(p);
		if (PolyCubeCutHelper::point_inside_quadrilateral(proj_p, fb.a_point))
		{
			return (p - proj_p).squared_length() <= m_threshold;
		}
		else
		{
			for (int j = 0; j < 4; ++j)
			{
				// compute the length between proj_p(P) and edge AB(a_point[j] -> a_point[(j + 1) % 4])
				// https://blog.csdn.net/angelazy/article/details/38489293
				Vector ap = p - fb.a_point[j]; // AP
				double dot = ap * fb.a_edge_vector[j]; // AP * AB
				if (dot <= 0)
				{
					if (ap.squared_length() <= m_threshold) // |AP|
						return true;
				}
				else if (dot >= fb.a_edge_length_2[j])
				{
					if ((p - fb.a_point[(j + 1) % 4]).squared_length() <= m_threshold) // |BP|
						return true;
				}
				else
				{
					Vector ac = (dot / fb.a_edge_length_2[j]) * fb.a_edge_vector[j]; // AC
					if ((ap - ac).squared_length() <= m_threshold)
						return true;
				}
			}
			return false;
		}
	}
    return false;
}

// check whether the new cut facet is too close to previous facets
void PolyCubeCut::check_cut_spacing(
	const CutFacet& new_cf
)
{
	for (const auto& old_cf : m_v_cut_facet)
	{
		if (new_cf.a_vertex[0] == old_cf.a_vertex[0] || new_cf.a_vertex[0] == old_cf.a_vertex[1] ||
			new_cf.a_vertex[1] == old_cf.a_vertex[0] || new_cf.a_vertex[1] == old_cf.a_vertex[1])
			continue; // skip 3-cut
		if (is_facet_too_close(new_cf, old_cf) || is_facet_too_close(old_cf, new_cf))
		{
			std::string err_msg = "[Error] Cut facet (";
			for (int i = 0; i < 3; ++i)
				err_msg += "#" + std::to_string(new_cf.a_vertex[i]) + ", ";
			err_msg += "#" + std::to_string(new_cf.a_vertex[3]) + ") and (";
			for (int i = 0; i < 3; ++i)
				err_msg += "#" + std::to_string(old_cf.a_vertex[i]) + ", ";
			err_msg += "#" + std::to_string(old_cf.a_vertex[3]) + ") are too close.\n";
			throw err_msg;
		}
	}
}

// check whether the cut leads to intersection
void PolyCubeCut::check_intersection(
	const Halfedge_handle &cut_he,
	const std::array<Point, 2> &a_cut_end_vertex
)
{
	std::array<Triangle, 2> a_cut_triangle{
		Triangle(cut_he->opposite()->vertex()->point(), cut_he->vertex()->point(), a_cut_end_vertex[0]),
		Triangle(cut_he->vertex()->point(), a_cut_end_vertex[1], a_cut_end_vertex[0])
	};

	// check whether the cut intersect mesh faces
    std::set<int> s_related_chart;
    for (int i = 0; i < 3; ++i)
    {
        s_related_chart.insert(m_v_face2chart[m_va_vertex2_a_polycube_edge[cut_he->vertex()->id()][i]->face()->id()]);
        s_related_chart.insert(m_v_face2chart[m_va_vertex2_a_polycube_edge[cut_he->opposite()->vertex()->id()][i]->face()->id()]);
    }
    for (const auto &fit : faces(m_polymesh))
	{
		int current_chart = m_v_face2chart[fit->id()];
        if (s_related_chart.find(current_chart) == s_related_chart.end())
        {
            Halfedge_handle he = fit->halfedge();
            Triangle current_face(
                he->vertex()->point(),
                he->opposite()->vertex()->point(),
                he->next()->vertex()->point()
            );
            for (int i = 0; i < 2; ++i)
            {
                bool is_intersected = CGAL::do_intersect(current_face, a_cut_triangle[i]);
                if (is_intersected)
                {
                    std::string err_msg = "[Error] Face #" + std::to_string(fit->id()) + " intersects the triangle:\n";
                    for (int j = 0; j < 3; ++j)
                        err_msg += "#" + std::to_string(j) + "\t" +
                        std::to_string(a_cut_triangle[i][j].x()) + "\t" +
                        std::to_string(a_cut_triangle[i][j].y()) + "\t" +
                        std::to_string(a_cut_triangle[i][j].z()) + "\n";
                    throw err_msg;
                }
            }
        }
	}

	// check whether the cut intersect other cuts
	for (const auto cf : m_v_cut_facet)
	{
		bool is_adjacent = false;
		for (int i = 0; i < 2; ++i)
			if (cf.a_vertex[i] == cut_he->vertex()->id() || cf.a_vertex[i] == cut_he->opposite()->vertex()->id())
				is_adjacent = true;
		if (is_adjacent)
			continue;

		std::array<Triangle, 2> a_cf_triangle{
			Triangle(cf.a_point[0], cf.a_point[1], cf.a_point[2]),
			Triangle(cf.a_point[1], cf.a_point[2], cf.a_point[3])
		};
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j)
			{
				bool is_intersected = CGAL::do_intersect(a_cut_triangle[i], a_cf_triangle[j]);
				if (is_intersected)
				{
					std::string err_msg = "[Error] The cut facets of Edge(#" + std::to_string(cf.a_vertex[0]) + ", #" +
						std::to_string(cf.a_vertex[1]) + ") and Edge(#" + std::to_string(cut_he->vertex()->id()) + ", #" +
						std::to_string(cut_he->opposite()->vertex()->id()) + ") intersect each other.\n";
					throw err_msg;
				}
			}
	}
}

bool PolyCubeCut::find_intersection_point_in_polygon(
	const Halfedge_handle &he_begin,
	const Halfedge_handle &he_end,
	const Segment &cut_segment,
	Halfedge_handle *he_intersected,
	Point *intersection_point
)
{
	for (Halfedge_handle hit = he_begin; hit != he_end; hit = hit->next())
	{
		Segment current_edge(hit->vertex()->point(), hit->opposite()->vertex()->point());
		auto result = CGAL::intersection(cut_segment, current_edge);
		if (result)
		{
			*he_intersected = hit;
			*intersection_point = *boost::get<Point>(&*result);
			return true;
		}
	}
	return false;
}

void PolyCubeCut::do_cut_in_tetgen_helper(
	const std::array<Vertex_handle, 2> &a_cut_start_vertex_handle,
	const std::array<unsigned int, 2> &a_cut_end_vertex_id,
	const std::array<Point, 2> &a_cut_end_vertex,
	const Vector &cut_dir,
	const std::array<int, 2> &a_cut_chart,
	TetgenHelper *tetgen_helper
)
{
	std::array<unsigned int, 2> a_new_vertex_start_id;

	for (int i = 0; i < 2; ++i)
	{
        a_new_vertex_start_id[i] = static_cast<unsigned int>(tetgen_helper->pointlist.size());
        if (a_cut_chart[i] != -1)
        {
            Segment cut_segment(a_cut_start_vertex_handle[i]->point(), a_cut_end_vertex[i]);
            cut_chart_in_tetgen_helper(a_cut_start_vertex_handle[i], a_cut_end_vertex_id[i], a_cut_chart[i], cut_dir, cut_segment, tetgen_helper);
        }
	}

	// add the cut facet
	tetgen_helper->facetlist.push_back(FacetHelper());
	tetgen_helper->facetlist.back().polygonlist.push_back(PolygonHelper());
	tetgen_helper->facetlist.back().chart = -2;
	
	PolygonHelper &p = tetgen_helper->facetlist.back().polygonlist.back();
	p.vertexlist.push_back(static_cast<unsigned int>(a_cut_start_vertex_handle[0]->id()));
	for (unsigned int i = a_new_vertex_start_id[0]; i < a_new_vertex_start_id[1]; ++i)
		p.vertexlist.push_back(i);
	p.vertexlist.push_back(a_cut_end_vertex_id[0]);
	p.vertexlist.push_back(a_cut_end_vertex_id[1]);
	for (unsigned int i = static_cast<unsigned int>(tetgen_helper->pointlist.size()) - 1; i >= a_new_vertex_start_id[1]; --i)
		p.vertexlist.push_back(i);
	p.vertexlist.push_back(static_cast<unsigned int>(a_cut_start_vertex_handle[1]->id()));

#ifdef RECONSTRUCTION_DETAIL
	std::cout << "-------------------Cut facet-------------------\n";
	std::cout << "Polygon:";
	for (const auto &v : p.vertexlist)
		std::cout << " " << v;
	std::cout << std::endl;
#endif
}

void PolyCubeCut::cut_chart_in_tetgen_helper(
	const Vertex_handle &cut_start_vertex_handle,
	const unsigned int &cut_end_vertex_id,
	const int &cut_chart,
	const Vector &cut_dir,
	const Segment &cut_segment,
	TetgenHelper *tetgen_helper
)
{
	Halfedge_handle he_intersected, hit = cut_start_vertex_handle->halfedge();
	unsigned int last_vertex_id = static_cast<unsigned int>(cut_start_vertex_handle->id());
	bool cut_finished = false;

	do
	{
		int current_face = static_cast<int>(hit->face()->id());
		int current_chart = m_v_face2chart[current_face];
		if (current_chart == cut_chart)
		{
			std::array<Vector, 2> a_edge_vector{ // two direction vectors of the edges at the start point
				hit->next()->vertex()->point() - hit->vertex()->point(),
				hit->prev()->vertex()->point() - hit->vertex()->point()
			};

			std::array<Vector, 2> a_cp_edge_cut{
				CGAL::cross_product(a_edge_vector[0], cut_dir),
				CGAL::cross_product(a_edge_vector[1], cut_dir)
			};

			if (a_cp_edge_cut[0].squared_length() < EPSILON) // the cut is coincident with hit->next()
			{
				split_edge_in_tetgen_helper(hit->next(), static_cast<unsigned int>(cut_start_vertex_handle->id()), cut_end_vertex_id, tetgen_helper);
				return;
			}

			if (a_cp_edge_cut[1].squared_length() < EPSILON) // the cut is coincident with hit
			{
				split_edge_in_tetgen_helper(hit, static_cast<unsigned int>(cut_start_vertex_handle->id()), cut_end_vertex_id, tetgen_helper);
				return;
			}

			// find the first polygon that contains the cut
			if (a_cp_edge_cut[0] * a_cp_edge_cut[1] < 0) // the cut is between two edges
			{
				Point intersection_point;
				if (find_intersection_point_in_polygon(hit->next()->next(), hit, cut_segment, &he_intersected, &intersection_point))
				{
					
					if ((intersection_point - tetgen_helper->pointlist[cut_end_vertex_id]).squared_length() < EPSILON)
					{
						// The intersection point is too close to the end vertex. Discard the intersection point.
						last_vertex_id = cut_end_vertex_id;
					}
					else
					{
						// add new vertex to point list
						last_vertex_id = static_cast<unsigned int>(tetgen_helper->pointlist.size());
						tetgen_helper->pointlist.push_back(intersection_point);
#ifdef _DEBUG
						std::cout << "[Info] Intersection point #" << last_vertex_id << ": " << intersection_point << std::endl;
#endif
					}
					insert_intersection_point(
						1,
						last_vertex_id,
						static_cast<unsigned int>(cut_start_vertex_handle->id()),
						tetgen_helper->pointlist,
						&tetgen_helper->facetlist[current_face]
					);
				}
				else
					cut_finished = true;
				break;
			}
		}

		hit = hit->next_on_vertex();
	} while (hit != cut_start_vertex_handle->halfedge());

	while (!cut_finished)
	{
		hit = he_intersected->opposite();
		Halfedge_handle last_he_intersected = hit;
		Point intersection_point;
		if (find_intersection_point_in_polygon(hit->next(), hit, cut_segment, &he_intersected, &intersection_point))
		{
			last_vertex_id = static_cast<unsigned int>(tetgen_helper->pointlist.size());
			
			if ((intersection_point - tetgen_helper->pointlist[cut_end_vertex_id]).squared_length() < EPSILON)
			{
				// The intersection point is too close to the end vertex. Discard the intersection point.
				insert_intersection_point(
					2,
					cut_end_vertex_id,
					last_vertex_id - 1,
					tetgen_helper->pointlist,
					&tetgen_helper->facetlist[hit->face()->id()]
				);
				last_vertex_id = cut_end_vertex_id;
			}
			else
			{
				// add new vertex to point list
				tetgen_helper->pointlist.push_back(intersection_point);

#ifdef _DEBUG
				std::cout << "[Info] Intersection point #" << last_vertex_id << ": " << intersection_point << std::endl;
#endif
				insert_intersection_point(
					2,
					last_vertex_id,
					last_vertex_id - 1,
					tetgen_helper->pointlist,
					&tetgen_helper->facetlist[hit->face()->id()]
				);
			}
		}
		else
			cut_finished = true;
	}

	if (last_vertex_id == cut_start_vertex_handle->id()) // no intersection with edges of the current face
	{
		unsigned int current_face = static_cast<unsigned int>(hit->face()->id());
		tetgen_helper->facetlist[current_face].polygonlist.push_back(PolygonHelper());
		PolygonHelper &p_seg = tetgen_helper->facetlist[current_face].polygonlist.back();
		p_seg.vertexlist.push_back(last_vertex_id);
		p_seg.vertexlist.push_back(cut_end_vertex_id);
	}
	else
	{
		he_intersected = he_intersected->opposite();
		unsigned int current_face = static_cast<unsigned int>(he_intersected->face()->id());
		insert_intersection_point(
			1,
			last_vertex_id,
			cut_end_vertex_id,
			tetgen_helper->pointlist,
			&tetgen_helper->facetlist[current_face]
		);
	}
}

void PolyCubeCut::split_edge_in_tetgen_helper(const Halfedge_handle &he, const unsigned int &start, const unsigned int &end, TetgenHelper *tetgen_helper)
{
	split_edge_in_polygon_helper(he, start, end, &tetgen_helper->facetlist[he->face()->id()].polygonlist[0]);
	split_edge_in_polygon_helper(he->opposite(), start, end, &tetgen_helper->facetlist[he->opposite()->face()->id()].polygonlist[0]);
}

void PolyCubeCut::split_edge_in_polygon_helper(const Halfedge_handle &he, const unsigned int &start,  const unsigned int &end, PolygonHelper *p)
{
	if (he->vertex()->id() == start)
	{
		if (p->vertexlist[0] == start)
			p->vertexlist.push_back(end);
		else
		{
			auto it = p->vertexlist.begin() + 1;
			for (; it != p->vertexlist.end(); ++it)
				if (*it == start)
					break;
			p->vertexlist.insert(it, end);
		}
	}
	else
	{
		auto it = p->vertexlist.begin();
		for (; it != p->vertexlist.end(); ++it)
			if (*it == start)
				break;
		p->vertexlist.insert(it + 1, end);
	}
}

void PolyCubeCut::insert_intersection_point(
	const int &n_intersection_point,
	const unsigned int &intersection_point_id,
	const unsigned int &start_or_end_vertex_id,
	const std::vector<Point> &pointlist,
	FacetHelper *f
)
{
#ifdef RECONSTRUCTION_DETAIL
	std::cout << "-------------Before reconstruction-------------\n";
	for (const auto &p : f->polygonlist)
	{
		std::cout << "Polygon:";
		for (const auto &v : p.vertexlist)
			std::cout << " " << v;
		std::cout << std::endl;
	}
#endif

	// insert the intersection point(s) to this polygon
	std::vector<unsigned int> &vl = f->polygonlist[0].vertexlist;
	std::array<unsigned int, 2> a_index{ intersection_point_id, start_or_end_vertex_id };

	for (int i = 0; i < n_intersection_point; ++i)
	{
        unsigned int new_id = a_index[i];
		const Point &p_new = pointlist[new_id];

		for (int j = 0; j < vl.size(); ++j)
		{
			const Point &p = pointlist[vl[j]], &p_next = pointlist[vl[(j + 1) % vl.size()]];
			std::array<Vector, 2> a_dir_vector{
				p_new - p,
				p_new - p_next
			};

			if (CGAL::cross_product(a_dir_vector[0], a_dir_vector[1]).squared_length() < EPSILON && // collinear
				a_dir_vector[0] * a_dir_vector[1] < 0) 
				vl.insert(vl.begin() + j + 1, new_id);
		}
	}

	// add the cut segment
	if (intersection_point_id != start_or_end_vertex_id)
	{
		f->polygonlist.push_back(PolygonHelper());
		PolygonHelper &p_seg = f->polygonlist.back();
		p_seg.vertexlist.push_back(start_or_end_vertex_id);
		p_seg.vertexlist.push_back(intersection_point_id);
	}

#ifdef RECONSTRUCTION_DETAIL
	std::cout << "--------------After reconstruction-------------\n";
	for (const auto &p : f->polygonlist)
	{
		std::cout << "Polygon:";
		for (const auto &v : p.vertexlist)
			std::cout << " " << v;
		std::cout << std::endl;
	}
#endif
}

// ensure that there is at most only one boundary face in each tetrahedron by subdivision
void PolyCubeCut::tetrahedron_subdivision(tetgenio *tetout)
{
    // indices of the tetrahedrons which should be subdivided
    std::vector<int> v_subdivision_tet_id;

    std::vector<bool> v_is_boundary_vertex(tetout->numberofpoints, false);

    for (int i = 0; i < tetout->numberoftrifaces; ++i)
        if (tetout->trifacemarkerlist[i] != 0)
            for (int j = 0; j < 3; ++j)
                v_is_boundary_vertex[tetout->trifacelist[3 * i + j]] = true;

	for (int i = 0; i < tetout->numberoftetrahedra; ++i) // for each tet
	{
        // a tet should be subdivided if and only if all four vertices of this tet are on the boundary
        int j = 0;
		for (; j < 4; ++j) // for each vertex
            if (!v_is_boundary_vertex[tetout->tetrahedronlist[4 * i + j]])
                break;

		if (j == 4) // all vertices of the tet are on the boundary
			v_subdivision_tet_id.push_back(i);
	}

	int n_subdivision_tet = static_cast<int>(v_subdivision_tet_id.size());

	// each subdivision operation results in one new vertex in pointlist
	double *new_pointlist = new double[3 * (tetout->numberofpoints + n_subdivision_tet)];
	// copy old pointlist to the new one
	memcpy(new_pointlist, tetout->pointlist, 3 * tetout->numberofpoints * sizeof(double));

	// each subdivision operation results in six new faces in trifacelist
	int *new_trifacelist = new int[3 * (tetout->numberoftrifaces + 6 * n_subdivision_tet)];
	// copy old trifacelist to the new one
	memcpy(new_trifacelist, tetout->trifacelist, 3 * tetout->numberoftrifaces * sizeof(int));

	// each subdivision operation results in three new tetrahedrons in tetrahedronlist
	int *new_tetrahedronlist = new int[4 * (tetout->numberoftetrahedra + 3 * n_subdivision_tet)];
	// copy old tetrahedronlist to the new one
	memcpy(new_tetrahedronlist, tetout->tetrahedronlist, 4 * tetout->numberoftetrahedra * sizeof(int));

	// each subdivision operation results in three new tetrahedrons
	int *new_tet2facelist = new int[4 * (tetout->numberoftetrahedra + 3 * n_subdivision_tet)];
	// copy old tet2facelist to the new one
	memcpy(new_tet2facelist, tetout->tet2facelist, 4 * tetout->numberoftetrahedra * sizeof(int));

	delete[] tetout->pointlist;
	delete[] tetout->trifacelist;
	delete[] tetout->tetrahedronlist;
	delete[] tetout->tet2facelist;
	tetout->pointlist = new_pointlist;
	tetout->trifacelist = new_trifacelist;
	tetout->tetrahedronlist = new_tetrahedronlist;
	tetout->tet2facelist = new_tet2facelist;

	for (int tid : v_subdivision_tet_id)
	{
		int new_vert_id = tetout->numberofpoints;

		// compute the centroid of the tetrahedra
		double tet_centroid[3] = { 0.0, 0.0, 0.0 };
		for (int j = 0; j < 4; ++j) // for each vertex
		{
			int vidx = new_tetrahedronlist[4 * tid + j]; // the index of j-th vertex of tid-th tetrahedra
			for (int k = 0; k < 3; ++k)
				tet_centroid[k] += new_pointlist[3 * vidx + k];
		}
		for (int k = 0; k < 3; ++k)
			tet_centroid[k] /= 4;

		// add the centroid to new pointlist
		memcpy(new_pointlist + 3 * tetout->numberofpoints, &tet_centroid, 3 * sizeof(double));
		++tetout->numberofpoints;

		// each pair forms a new face with the centroid
		static int s_face_sort_id[6][2] = { {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3} };
		int first_new_face_id = tetout->numberoftrifaces;
		for (int j = 0; j < 6; ++j)
		{
			int new_face_id = tetout->numberoftrifaces;
			new_trifacelist[3 * new_face_id] = new_tetrahedronlist[4 * tid + s_face_sort_id[j][0]];
			new_trifacelist[3 * new_face_id + 1] = new_tetrahedronlist[4 * tid + s_face_sort_id[j][1]];
			new_trifacelist[3 * new_face_id + 2] = new_vert_id;
			++tetout->numberoftrifaces;
		}

		// each triplet forms a new tetrahedron with the centroid
		static int s_tet_sort_id[4][3] = { {2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2} }; // 4 faces of tet abcd: cbd acd bad abc
		static int s_tet2face_offset[4][3] = { {4, 5, 3}, {5, 2, 1}, {2, 4, 0}, {3, 1, 0} };
		for (int j = 0; j < 3; ++j) // add three new tetrahedrons
		{
			int new_tet_id = tetout->numberoftetrahedra;

			new_tet2facelist[4 * new_tet_id + 3] = new_tet2facelist[4 * tid + j];
			for (int k = 0; k < 3; ++k)
			{
				new_tetrahedronlist[4 * new_tet_id + k] = new_tetrahedronlist[4 * tid + s_tet_sort_id[j][k]];
				new_tet2facelist[4 * new_tet_id + k] = first_new_face_id + s_tet2face_offset[j][k];
			}
			new_tetrahedronlist[4 * new_tet_id + 3] = new_vert_id;

			++tetout->numberoftetrahedra;
		}
		// update the original tetrahedron as the last one
		new_tetrahedronlist[4 * tid + 3] = new_vert_id;
		new_tet2facelist[4 * tid + 3] = new_tet2facelist[4 * tid + 3];
		for (int k = 0; k < 3; ++k)
			new_tet2facelist[4 * tid + k] = first_new_face_id + s_tet2face_offset[3][k];
	}
}

void PolyCubeCut::update_boundary_tet_chart(int marker_max_index, tetgenio *tetout, std::vector<int> *v_tet_chart)
{
	for (int i = 0; i < tetout->numberoftetrahedra; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			int fidx = tetout->tet2facelist[4 * i + j];
			if (fidx < marker_max_index && tetout->trifacemarkerlist[fidx] != 0)
			{
				(*v_tet_chart)[i] = tetout->trifacemarkerlist[fidx] - 1;
				break;
			}
		}
	}
}

void PolyCubeCut::separate_cut_facet(
	const std::vector<int> &v_cut_type,
	const std::vector<int> &v_cut_section_chart,
	std::vector<int> *v_tet_chart,
	std::vector<std::vector<int>> *vv_congruent_face,
    std::map<int, ThreeCut> *map_vertex2three_cut,
	tetgenio *tetout
)
{
	// the number of cut facet
	int n_cut = static_cast<int>(v_cut_type.size());

	// indices of the faces on each cut facet
	std::vector<std::vector<int>> vv_cut_face_id(n_cut);

	// indices of the points on each cut facet (except the points on seam edges)
	std::vector<std::set<int>> vs_cut_facet_vertex_id(n_cut);

    // key: polycube vertex index
    // value: int - the number of concave cut, set - indices of the vertices on the interior edge
    std::map<int, std::pair<int, std::set<int>>> ms_vertex2interior_edge_info;
    m_v_duplicated_vertex_map.clear();
    std::set<int> s_id_vertex_duplicate_1; // create one duplicate
    std::set<int> s_id_vertex_duplicate_2; // create two duplicates

	for (int k = 0; k < n_cut; ++k)
	{
        CutFacet &facet = (m_v_cut_facet)[k];

        for (int i = 0; i < tetout->numberofpoints; ++i)
        {
            Point p(tetout->pointlist[3 * i], tetout->pointlist[3 * i + 1], tetout->pointlist[3 * i + 2]);

            PolyCubeCutHelper::LOCATION_RESULT location = PolyCubeCutHelper::locate_point_in_cut_facet(p, facet);

            if (location == PolyCubeCutHelper::ON_INTERIOR_EDGE_0 || location == PolyCubeCutHelper::ON_INTERIOR_EDGE_1)
            {
                s_id_vertex_duplicate_2.insert(i);
                if (ms_vertex2interior_edge_info.find(facet.a_vertex[location]) == ms_vertex2interior_edge_info.end())
                    ms_vertex2interior_edge_info[facet.a_vertex[location]].first = !facet.is_convex;
                else
                    ms_vertex2interior_edge_info[facet.a_vertex[location]].first += !facet.is_convex;
                ms_vertex2interior_edge_info[facet.a_vertex[location]].second.insert(i);
            }
            else if (location == PolyCubeCutHelper::BORDER_OR_INSIDE)
                s_id_vertex_duplicate_1.insert(i);
        }
	}

	duplicate_cut_facet_vertex(s_id_vertex_duplicate_1, s_id_vertex_duplicate_2, tetout);
	update_tetrahedron(v_cut_section_chart, v_tet_chart, vv_congruent_face, tetout);

    for (const auto &pair : ms_vertex2interior_edge_info)
    {
        const auto &info = pair.second;
        if (info.first & 1)
            for (const auto &id : info.second)
                (*map_vertex2three_cut)[pair.first].va_corresponding_vertex.push_back(
                    {
                        static_cast<unsigned int>(m_v_duplicated_vertex_map.at(id)[1]),
                        static_cast<unsigned int>(m_v_duplicated_vertex_map.at(id)[0]),
                        static_cast<unsigned int>(id)
                    }
                );
        else
            for (const auto &id : info.second)
                (*map_vertex2three_cut)[pair.first].va_corresponding_vertex.push_back(
                    {
                        static_cast<unsigned int>(m_v_duplicated_vertex_map.at(id)[0]),
                        static_cast<unsigned int>(m_v_duplicated_vertex_map.at(id)[1]),
                        static_cast<unsigned int>(id)
                    }
                );
    }
}

// copy old vertices and duplicate the vertices on the cut facets to a new pointlist
// for each old vertex, record its corresponding new vertex index in the map
void PolyCubeCut::duplicate_cut_facet_vertex(
    const std::set<int> &s_id_vertex_duplicate_1,
    const std::set<int> &s_id_vertex_duplicate_2,
	tetgenio *tetout
)
{
	// the number of new points
    int n_new_point = static_cast<int>(s_id_vertex_duplicate_1.size() + (s_id_vertex_duplicate_2.size() << 1));

	double *new_pointlist = new double[3 * (n_new_point + tetout->numberofpoints)];

	// copy old vertices
	memcpy(new_pointlist, tetout->pointlist, 3 * sizeof(double) * tetout->numberofpoints);

    m_v_duplicated_vertex_map.resize(n_new_point + tetout->numberofpoints);
    for (int i = 0; i < tetout->numberofpoints; ++i)
        m_v_duplicated_vertex_map[i] = { i, -1 }; // map the vertex to itself if it is not duplicated

	delete[] tetout->pointlist; // delete old pointlist to avoid memory leaks
	tetout->pointlist = new_pointlist;

	int i = tetout->numberofpoints;
    for (const auto &id : s_id_vertex_duplicate_1)
    {
        memcpy(new_pointlist + 3 * i, &new_pointlist[3 * id], 3 * sizeof(double));
        m_v_duplicated_vertex_map[id] = { i, -1 };
        m_v_duplicated_vertex_map[i] = { id, -1 };
        ++i;
    }

    for (const auto &id : s_id_vertex_duplicate_2)
    {
        memcpy(new_pointlist + 3 * i++, &new_pointlist[3 * id], 3 * sizeof(double));
        memcpy(new_pointlist + 3 * i, &new_pointlist[3 * id], 3 * sizeof(double));
        m_v_duplicated_vertex_map[id] = { i - 1, i };
        m_v_duplicated_vertex_map[i - 1] = { id, -1 };
        m_v_duplicated_vertex_map[i] = { id, -1 };
        ++i;
    }

	tetout->numberofpoints = i;
}

void PolyCubeCut::update_tetrahedron(
	const std::vector<int> &v_cut_section_chart,
	std::vector<int> *v_tet_chart,
    std::vector<std::vector<int>> *vv_congruent_face,
	tetgenio *tetout
)
{
    // indices of tets on the cut facet (at least one vertex on an interior edge)
    // array for two side of the facet, 0 - positive, 1 - negative
    std::vector<std::array<std::vector<int>, 2>> vav_tet_on_facet(m_v_cut_facet.size());
    // tet index -> centroid of the face on cut facet
    std::map<int, Point> map_tet2face_centroid;
    // tet index -> 3 vertices of the face on cut facet
    std::map<int, std::vector<int>> mv_tet2face;

    std::vector<Point> v_tet_point;
    v_tet_point.reserve(tetout->numberofpoints);
    for (int i = 0; i < tetout->numberofpoints; ++i)
        v_tet_point.push_back(
            Point(
                tetout->pointlist[3 * i],
                tetout->pointlist[3 * i + 1],
                tetout->pointlist[3 * i + 2]
            )
        );

	for (int k = 0; k < m_v_cut_facet.size(); ++k) // for each cut facet
	{
        const CutFacet &facet = m_v_cut_facet[k];
		for (int i = 0; i < tetout->numberoftetrahedra; ++i) // for each tetrahedra
		{
			std::vector<int> v_offset_vertex_on_cut_facet; // offset (0-3) list in tet_i of the vertices that are on the cut facet
            std::array<std::vector<int>, 2> av_offset_vertex_on_interior_edge;
			
            for (int j = 0; j < 4; ++j) // for each vertex of the tetrahedra
            {
                int vidx = tetout->tetrahedronlist[4 * i + j]; // the index of j-th vertex of i-th tetrahedra
                Point p(tetout->pointlist[3 * vidx], tetout->pointlist[3 * vidx + 1], tetout->pointlist[3 * vidx + 2]);

                PolyCubeCutHelper::LOCATION_RESULT location = PolyCubeCutHelper::locate_point_in_cut_facet(p, facet);
                if (location == PolyCubeCutHelper::BORDER_OR_INSIDE || location == PolyCubeCutHelper::ON_SEAM_EDGE)
                    v_offset_vertex_on_cut_facet.push_back(j);
                else if (location == PolyCubeCutHelper::ON_INTERIOR_EDGE_0)
                    av_offset_vertex_on_interior_edge[0].push_back(j);
                else if (location == PolyCubeCutHelper::ON_INTERIOR_EDGE_1)
                    av_offset_vertex_on_interior_edge[1].push_back(j);
			}

			if (!(v_offset_vertex_on_cut_facet.empty() &&
                av_offset_vertex_on_interior_edge[0].empty() &&
                av_offset_vertex_on_interior_edge[1].empty())) // at least one vertex on the cut facet
			{
                std::array<Point, 4> a_tet_point; // 4 vertex positions of the current tetrahedra
                for (int j = 0; j < 4; ++j) // for each vertex of the tetrahedra
                {
                    int vidx = tetout->tetrahedronlist[4 * i + j]; // the index of j-th vertex of i-th tetrahedra
                    a_tet_point[j] = Point(
                        tetout->pointlist[3 * vidx],
                        tetout->pointlist[3 * vidx + 1],
                        tetout->pointlist[3 * vidx + 2]
                    );
                }
                // compute the centroid of the tetrahedra
				Point tet_centroid(
					(a_tet_point[0].x() + a_tet_point[1].x() + a_tet_point[2].x() + a_tet_point[3].x()) / 4,
					(a_tet_point[0].y() + a_tet_point[1].y() + a_tet_point[2].y() + a_tet_point[3].y()) / 4,
					(a_tet_point[0].z() + a_tet_point[1].z() + a_tet_point[2].z() + a_tet_point[3].z()) / 4
				);

                bool has_congruent_face = false;
                std::vector<int> old_face, new_face;

                Point point_on_facet(
                    tetout->pointlist[3 * facet.a_vertex[0]],
                    tetout->pointlist[3 * facet.a_vertex[0] + 1],
                    tetout->pointlist[3 * facet.a_vertex[0] + 2]
                );
                // compute the relative position of the centroid
                if ((tet_centroid - point_on_facet) * facet.normal > 0)
                {
                    if ((v_offset_vertex_on_cut_facet.size() +
                        av_offset_vertex_on_interior_edge[0].size() +
                        av_offset_vertex_on_interior_edge[1].size()) == 3)
                    {
                        // record the congruent face
                        has_congruent_face = true;

                        if ((*v_tet_chart)[i] == -2) // cut facet
                            (*v_tet_chart)[i] = v_cut_section_chart[(k << 2) + 1]; // assign a new chart

                        if (!(av_offset_vertex_on_interior_edge[0].empty() && av_offset_vertex_on_interior_edge[1].empty()))
                            vav_tet_on_facet[k][0].push_back(i);
                    }

                    for (int j : v_offset_vertex_on_cut_facet)
                    {
                        int &vidx = tetout->tetrahedronlist[4 * i + j];
                        if (vidx > m_v_duplicated_vertex_map[vidx][0])
                            vidx = m_v_duplicated_vertex_map[vidx][0];
                        if (has_congruent_face)
                            old_face.push_back(vidx);
                        vidx = m_v_duplicated_vertex_map[vidx][0];
                        if (has_congruent_face)
                            new_face.push_back(vidx);
                    }
                }
                else
                {
                    if ((v_offset_vertex_on_cut_facet.size() +
                        av_offset_vertex_on_interior_edge[0].size() +
                        av_offset_vertex_on_interior_edge[1].size()) == 3 &&
                        (*v_tet_chart)[i] == -2) // cut facet
                    {
                        (*v_tet_chart)[i] = v_cut_section_chart[(k << 2) + 2]; // assign a new chart

                        if (!(av_offset_vertex_on_interior_edge[0].empty() && av_offset_vertex_on_interior_edge[1].empty()))
                            vav_tet_on_facet[k][1].push_back(i);
                    }
                }

                for (int j : av_offset_vertex_on_interior_edge[0])
                {
                    int &vidx = tetout->tetrahedronlist[4 * i + j];
                    if (vidx > m_v_duplicated_vertex_map[vidx][0])
                        vidx = m_v_duplicated_vertex_map[vidx][0];
                    Point va(tetout->pointlist[3 * facet.a_vertex[0]], tetout->pointlist[3 * facet.a_vertex[0] + 1], tetout->pointlist[3 * facet.a_vertex[0] + 2]);
                    Point vb(tetout->pointlist[3 * facet.a_vertex[3]], tetout->pointlist[3 * facet.a_vertex[3] + 1], tetout->pointlist[3 * facet.a_vertex[3] + 2]);
                    int region = region_on_tangent_plane(facet.a_vertex[0], va - vb, tet_centroid);
                    if (region < 2)
                        vidx = m_v_duplicated_vertex_map[vidx][region];
                }

                for (int j : av_offset_vertex_on_interior_edge[1])
                {
                    int &vidx = tetout->tetrahedronlist[4 * i + j];
                    if (vidx > m_v_duplicated_vertex_map[vidx][0])
                        vidx = m_v_duplicated_vertex_map[vidx][0];
                    Point va(tetout->pointlist[3 * facet.a_vertex[1]], tetout->pointlist[3 * facet.a_vertex[1] + 1], tetout->pointlist[3 * facet.a_vertex[1] + 2]);
                    Point vb(tetout->pointlist[3 * facet.a_vertex[2]], tetout->pointlist[3 * facet.a_vertex[2] + 1], tetout->pointlist[3 * facet.a_vertex[2] + 2]);
                    int region = region_on_tangent_plane(facet.a_vertex[1], va - vb, tet_centroid);
                    if (region < 2)
                        vidx = m_v_duplicated_vertex_map[vidx][region];
                }

                if (has_congruent_face && old_face.size() == 3)
                {
                    (*vv_congruent_face)[k].insert((*vv_congruent_face)[k].end(), old_face.begin(), old_face.end());
                    (*vv_congruent_face)[k].insert((*vv_congruent_face)[k].end(), new_face.begin(), new_face.end());
                }
			}
		}
	}

    // find congruent faces on interior edges
    for (int k = 0; k < m_v_cut_facet.size(); ++k)
    {
        const CutFacet &facet = m_v_cut_facet[k];
        for (int side = 0; side < 2; ++side)
        {
            for (const auto &i : vav_tet_on_facet[k][side])
            {
                double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
                for (int j = 0; j < 4; ++j) // for each vertex of the tetrahedra
                {
                    int vidx = tetout->tetrahedronlist[4 * i + j]; // the index of j-th vertex of i-th tetrahedra
                    Point p(tetout->pointlist[3 * vidx], tetout->pointlist[3 * vidx + 1], tetout->pointlist[3 * vidx + 2]);

                    PolyCubeCutHelper::LOCATION_RESULT location = PolyCubeCutHelper::locate_point_in_cut_facet(p, facet);
                    if (location != PolyCubeCutHelper::OUTSIDE)
                    {
                        int vidx = tetout->tetrahedronlist[4 * i + j]; // the index of j-th vertex of i-th tetrahedra
                        mv_tet2face[i].push_back(vidx);
                        sum_x += tetout->pointlist[3 * vidx];
                        sum_y += tetout->pointlist[3 * vidx + 1];
                        sum_z += tetout->pointlist[3 * vidx + 2];
                    }
                }
                assert(3 == mv_tet2face[i].size());
                map_tet2face_centroid[i] = Point(sum_x / 3.0, sum_y / 3.0, sum_z / 3.0);
            }
        }

        for (const auto &i : vav_tet_on_facet[k][0])
        {
            for (const auto &j : vav_tet_on_facet[k][1])
            {
                if ((map_tet2face_centroid[i] - map_tet2face_centroid[j]).squared_length() < EPSILON)
                {
                    std::array<int, 3> a_ordered_face;
                    for (int x = 0; x < 3; ++x)
                    {
                        (*vv_congruent_face)[k].push_back(mv_tet2face[j][x]);
                        int origin_id = mv_tet2face[j][x] > m_v_duplicated_vertex_map[mv_tet2face[j][x]][0] ? m_v_duplicated_vertex_map[mv_tet2face[j][x]][0] : mv_tet2face[j][x];
                        for (int y = 0; y < 3; ++y)
                            if (mv_tet2face[i][y] == m_v_duplicated_vertex_map[origin_id][0] ||
                                mv_tet2face[i][y] == m_v_duplicated_vertex_map[origin_id][1] ||
                                mv_tet2face[i][y] == origin_id)
                            {
                                a_ordered_face[x] = mv_tet2face[i][y];
                                break;
                            }
                    }
                    (*vv_congruent_face)[k].insert((*vv_congruent_face)[k].end(), a_ordered_face.begin(), a_ordered_face.end());
                }
            }
        }
    }
}

int PolyCubeCut::region_on_tangent_plane(const unsigned int &vid, const Vector &normal, const Point &tet_centroid)
{
    Point origin(m_v_id2vertex[vid]->point());
    Plane tan_plane(origin, normal);
    Point proj_tet_centroid = tan_plane.projection(tet_centroid);
    std::array<Point, 3> proj_polycube_edge_point;
    for (int i = 0; i < 3; ++i)
        proj_polycube_edge_point[i] = tan_plane.projection(m_va_vertex2_a_polycube_edge[vid][i]->opposite()->vertex()->point());
    int remaining_area;
    for (int i = 0; i < 3; ++i)
    {
        Vector origin2centroid(proj_tet_centroid - origin);
        origin2centroid /= CGAL::sqrt(origin2centroid.squared_length());
        std::array<Vector, 2> edge_dir{
            proj_polycube_edge_point[i] - origin,
            proj_polycube_edge_point[(i + 1) % 3] - origin
        };
        for (int j = 0; j < 2; ++j)
            edge_dir[j] /= CGAL::sqrt(edge_dir[j].squared_length());
        if (edge_dir[0] * edge_dir[1] < 0.5 - EPSILON)
            remaining_area = i;
        if (CGAL::cross_product(proj_polycube_edge_point[i] - origin, origin2centroid) * normal < 0 &&
            CGAL::cross_product(proj_polycube_edge_point[(i + 1) % 3] - origin, origin2centroid) * normal > 0)
            return i;
    }
    return remaining_area;
}
