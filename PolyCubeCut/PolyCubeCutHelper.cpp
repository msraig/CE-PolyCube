#include "PolyCubeCutHelper.h"

void PolyCubeCutHelper::off2cgal(std::string filename, Polyhedron *polymesh)
{
	std::ifstream input(filename);
	if (!input || !(input >> *polymesh) || polymesh->empty())
	{
		input.close();
		throw "[Error] Invalid input file: " + filename + ".\n";
	}
	input.close();

	if (!CGAL::is_triangle_mesh(*polymesh))
		throw "[Error] " + filename + " is not a triangle mesh.\n";

    assign_cgal_id(polymesh);
}

void PolyCubeCutHelper::assign_cgal_id(Polyhedron *polymesh)
{
    // associate indices to the vertices using the "id()" field of the vertex.
    int idx = 0;
    for (const auto &vit : vertices(*polymesh))
        vit->id() = idx++;

    // associate indices to the faces using the "id()" field of the face.
    idx = 0;
    for (const auto &fit : faces(*polymesh))
        fit->id() = idx++;

    // associate indices to the halfedges using the "id()" field of the halfedge.
    idx = 0;
    for (const auto &hit : halfedges(*polymesh))
        hit->id() = idx++;
}

void PolyCubeCutHelper::tetgen_helper2tetgen(const TetgenHelper &tetgen_helper, tetgenio *tetin)
{
	// copy vertices
	tetin->numberofpoints = (int)tetgen_helper.pointlist.size();
	tetin->pointlist = new double[tetgen_helper.pointlist.size() * 3];
    for (int i = 0; i < tetin->numberofpoints; ++i)
	    memcpy(tetin->pointlist + 3 * i, &tetgen_helper.pointlist[i], 3 * sizeof(double));

	// copy facets
	tetin->numberoffacets = (int)tetgen_helper.facetlist.size();
	tetin->facetlist = new tetgenio::facet[tetin->numberoffacets];
	tetin->facetmarkerlist = new int[tetin->numberoffacets];
	tetgenio::facet *f;
	tetgenio::polygon *p;
	for (int i = 0; i < tetin->numberoffacets; ++i)
	{
		const FacetHelper &my_f = tetgen_helper.facetlist[i];
		tetin->facetmarkerlist[i] = my_f.chart + 1;
		f = &tetin->facetlist[i];
		tetin->init(f);
		f->numberofpolygons = (int)my_f.polygonlist.size();
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		for (int j = 0; j < f->numberofpolygons; ++j)
		{
			const PolygonHelper &my_p = my_f.polygonlist[j];
			p = &f->polygonlist[j];
			tetin->init(p);
			p->numberofvertices = (int)my_p.vertexlist.size();
			p->vertexlist = new int[p->numberofvertices];
			for (int k = 0; k < p->numberofvertices; ++k)
				p->vertexlist[k] = (int)my_p.vertexlist[k];
		}
	}

	// for debug
	//std::cout << "\nPoint:\n";
	//for (int i = 0; i < tetin->numberofpoints; ++i)
	//	std::cout << i << "\t" << tetin->pointlist[3 * i] << "\t" << tetin->pointlist[3 * i + 1] << "\t" << tetin->pointlist[3 * i + 2] << std::endl;
	//std::cout << "\nFacet:\n";
	//for (int i = 0; i < tetin->numberoffacets; ++i)
	//{
	//	std::cout << "facet #" << i << "\tmarker: " << tetin->facetmarkerlist[i] << std::endl;
	//	for (int j = 0; j < tetin->facetlist[i].numberofpolygons; ++j)
	//	{
	//		std::cout << "polygon #" << j << ":";
	//		for (int k = 0; k < tetin->facetlist[i].polygonlist[j].numberofvertices; ++k)
	//			std::cout << " " << tetin->facetlist[i].polygonlist[j].vertexlist[k];
	//		std::cout << std::endl;
	//	}
	//}
}

void PolyCubeCutHelper::tetgen2vtk(tetgenio *tetout, std::string filename)
{
	std::ofstream output(filename);
	output << "# vtk DataFile Version 2.0\n" << "Mesh saved from PolyCubeCut\n";
	output << "ASCII\n" << "DATASET UNSTRUCTURED_GRID\n";
	output << "POINTS " << tetout->numberofpoints << " float\n";
	for (int i = 0; i < tetout->numberofpoints; ++i)
	{
		for (int j = 0; j < 3; ++j)
			output << tetout->pointlist[3 * i + j] << " ";
		output << std::endl;
	}
	output << "CELLS " << tetout->numberoftetrahedra << " " << 5 * tetout->numberoftetrahedra << std::endl;
	for (int i = 0; i < tetout->numberoftetrahedra; ++i)
	{
		output << "4";
		for (int j = 0; j < 4; ++j)
			output << " " << tetout->tetrahedronlist[4 * i + j];
		output << std::endl;
	}
	output << "CELL_TYPES " << tetout->numberoftetrahedra << std::endl;
	for (int i = 0; i < tetout->numberoftetrahedra; ++i)
		output << "10\n";
	output.close();
}

void PolyCubeCutHelper::cgal2tetgen_helper(const Polyhedron &polymesh, const std::vector<int> &v_face2chart, TetgenHelper *tetgen_helper)
{
	using namespace CGAL::Polygon_mesh_processing;

    for (const auto &vit : vertices(polymesh))
		tetgen_helper->pointlist.push_back(vit->point());

    for (const auto &fit : faces(polymesh))
	{
		tetgen_helper->facetlist.push_back(FacetHelper());
		tetgen_helper->facetlist.back().polygonlist.push_back(PolygonHelper());
		tetgen_helper->facetlist.back().chart = v_face2chart[fit->id()];

		PolygonHelper &p = tetgen_helper->facetlist.back().polygonlist.back();
		Halfedge_handle h = fit->halfedge();
		do
		{
			p.vertexlist.push_back((int)h->vertex()->id());
			h = h->next();
		} while (h != fit->halfedge());
	}
}

// test whether the point p lies inside the quadrilateral quad[0->1->2->3] (counter-clockwise)
// assume that the point and the rectangle are in the same plane
bool PolyCubeCutHelper::point_inside_quadrilateral(const Point &p, const std::array<Point, 4> &v_quad)
{
    if (CGAL::cross_product(v_quad[1] - v_quad[0], p - v_quad[0]) *
        CGAL::cross_product(v_quad[3] - v_quad[2], p - v_quad[2]) > 0 &&
        CGAL::cross_product(v_quad[2] - v_quad[1], p - v_quad[1]) *
        CGAL::cross_product(v_quad[0] - v_quad[3], p - v_quad[3]) > 0)
        return true;

    return false;
}

// test whether the point q is on the segment pa->pb
bool PolyCubeCutHelper::point_on_segment(const Point &q, const Point &pa, const Point &pb)
{
	if (CGAL::cross_product(pa - pb, q - pb).squared_length() < EPSILON &&
		std::min(pa.x(), pb.x()) - q.x() <= EPSILON && q.x() - std::max(pa.x(), pb.x()) <= EPSILON &&
		std::min(pa.y(), pb.y()) - q.y() <= EPSILON && q.y() - std::max(pa.y(), pb.y()) <= EPSILON &&
		std::min(pa.z(), pb.z()) - q.z() <= EPSILON && q.z() - std::max(pa.z(), pb.z()) <= EPSILON)
		return true;
    return false;
}

PolyCubeCutHelper::LOCATION_RESULT PolyCubeCutHelper::locate_point_in_cut_facet(const Point &p, const CutFacet &facet)
{
	if (point_on_segment(p, facet.a_point[2], facet.a_point[3]))
	{
		return ON_SEAM_EDGE;
	}
	else if (point_on_segment(p, facet.a_point[0], facet.a_point[3]))
	{
		if (facet.a_is_interior[0])
			return ON_INTERIOR_EDGE_0;
		else
			return BORDER_OR_INSIDE;
	}
	else if (point_on_segment(p, facet.a_point[1], facet.a_point[2]))
	{
		if (facet.a_is_interior[1])
			return ON_INTERIOR_EDGE_1;
		else
			return BORDER_OR_INSIDE;
	}
	else if (point_on_segment(p, facet.a_point[0], facet.a_point[1]))
	{
		return BORDER_OR_INSIDE;
	}
	else if (abs(facet.normal * (p - facet.a_point[0])) < EPSILON && PolyCubeCutHelper::point_inside_quadrilateral(p, facet.a_point))
	{
		return BORDER_OR_INSIDE;
	}
	return OUTSIDE;
}
