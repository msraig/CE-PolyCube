#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <map>
#include <set>
#include <queue>
#include "rapidxml_print.hpp"
#include "rapidxml.hpp"
#include "rapidxml_iterators.hpp"
#include "rapidxml_utils.hpp"

#include "TetStructure.h"


using namespace std;
using namespace rapidxml;
using ig::TetStructure;
using ig::CVec;
using ig::Tetrahedron;
using ig::TetVertex;

void get_boundary_vert_face(TetStructure<double>* tet_mesh, std::vector<double> &points, std::vector<int> &faces, std::vector<int> &s2v, bool ori_vert_flag = false, std::map<std::vector<int>, int> *face2cellidx = NULL)
{
	points.clear();
	faces.clear();
	s2v.clear();
	//VolumeSurfaceVertex.clear();
	std::vector<int> tmp_face;
	std::map<int, int> VolumeSurfaceVertex;
	std::map<int, int>::iterator map_it;
	int indexOnSurface = 0;

	int num_tet = (int)tet_mesh->tetras.size();
	std::vector<Tetrahedron<double>*> tetras = tet_mesh->get_tetras();

	static int s_tet_id[4][3] = { { 1, 2, 3}, {2, 0, 3}, {0, 1, 3}, {1, 0, 2} }; //to outer part

	if (!ori_vert_flag)
	{
		for (size_t i = 0; i < num_tet; i++)
		{
			if (!tetras[i]->is_onboundary())
			{
				continue;
			}
			//BoundaryFaceFlag[i] = true;

			int face_id = -1;
			for (size_t j = 0; j < 4; j++)
			{
				TetVertex<double>* tv_temp = tetras[i]->vertex[j];
				if (!(tv_temp->boundary))
				{
					face_id = (int)j;
					break;
				}
			}

			if (face_id == -1)
			{
				std::cout << "read vtk error" << std::endl;
				return;
			}

			//faceVertexHandle.clear();
			tmp_face.clear();
			for (size_t j = 0; j < 3; j++)
			{
				map_it = VolumeSurfaceVertex.find((int)tetras[i]->vertex[s_tet_id[face_id][j]]->id);
				if (map_it == VolumeSurfaceVertex.end())
				{
					VolumeSurfaceVertex.insert(std::pair<int, int>(tetras[i]->vertex[s_tet_id[face_id][j]]->id, indexOnSurface));
					s2v.push_back((int)tetras[i]->vertex[s_tet_id[face_id][j]]->id);
					double px, py, pz;
					px = tetras[i]->vertex[s_tet_id[face_id][j]]->pos[0];
					py = tetras[i]->vertex[s_tet_id[face_id][j]]->pos[1];
					pz = tetras[i]->vertex[s_tet_id[face_id][j]]->pos[2];
					//SurfaceMesh::Point v2(px, py, pz);
					//SurfaceMesh::Point v2((px - cx)/max_diff, (py - cy)/max_diff, (pz - cz)/max_diff);
					points.push_back(px);
					points.push_back(py);
					points.push_back(pz);
					//vertexHandleVec.push_back(boundaryMesh.add_vertex(v2));
					//boundaryMesh.data(vertexHandleVec.back()).set_tet_vertex_id(tetras[i]->vertex[s_tet_id[face_id][j]]->id);
					//faceVertexHandle.push_back(vertexHandleVec.back());
					tmp_face.push_back(indexOnSurface);
					++indexOnSurface;
				}
				else
				{
					//faceVertexHandle.push_back(vertexHandleVec[map_it->second]);
					tmp_face.push_back(map_it->second);
				}
			}
			//OpenMesh::FaceHandle fh = boundaryMesh.add_face(faceVertexHandle);

			//boundaryMesh.data(fh).set_cell_id(i);

			faces.push_back(tmp_face[0]);
			faces.push_back(tmp_face[1]);
			faces.push_back(tmp_face[2]);


		}
	}
	else
	{
		//initialize points
		for (size_t i = 0; i < tet_mesh->tetra_vertices.size(); i++)
		{
			double px, py, pz;
			px = tet_mesh->get_vertices()[i]->pos[0];
			py = tet_mesh->get_vertices()[i]->pos[1];
			pz = tet_mesh->get_vertices()[i]->pos[2];
			points.push_back(px);
			points.push_back(py);
			points.push_back(pz);
		}

		for (size_t i = 0; i < num_tet; i++)
		{
			if (!tetras[i]->is_onboundary())
			{
				continue;
			}
			//BoundaryFaceFlag[i] = true;

			int face_id = -1;
			for (size_t j = 0; j < 4; j++)
			{
				TetVertex<double>* tv_temp = tetras[i]->vertex[j];
				if (!(tv_temp->boundary))
				{
					face_id = (int)j;
					break;
				}
			}

			if (face_id == -1)
			{
				std::cout << "read vtk error" << std::endl;
				return;
			}

			//faceVertexHandle.clear();
			tmp_face.clear();
			for (size_t j = 0; j < 3; j++)
			{
				tmp_face.push_back((int)tetras[i]->vertex[s_tet_id[face_id][j]]->id);
			}
			//OpenMesh::FaceHandle fh = boundaryMesh.add_face(faceVertexHandle);

			//boundaryMesh.data(fh).set_cell_id(i);

			faces.push_back(tmp_face[0]);
			faces.push_back(tmp_face[1]);
			faces.push_back(tmp_face[2]);

			if (face2cellidx != NULL)
			{
				//reorder tmp_face;
				sort(tmp_face.begin(), tmp_face.end());
				(*face2cellidx)[tmp_face] = (int)i;
			}
		}


	}

	//boundaryMesh.update_normals();


}


void string2vec_int(const string &s, vector<int> &vec)
{
	vec.clear();
	int posp = 0, posn = 0;

	while (posn < s.length())
	{
		posp = posn;
		//decide posp
		while ((s[posp] == ' ' || s[posp] == '\n') && posp < s.length())
		{
			posp++;
		}

		//decide posn:
		posn = posp;

		while ((s[posn] != ' ' && s[posn] != '\n') && posn < s.length())
		{
			posn++;
		}
		if (posn < s.length() && posn != posp)
		{
			int tmp = stoi(s.substr(posp, posn - posp));
			vec.push_back(tmp);
		}

	}
}

void string2vec_double(const string &s, vector<double> &vec)
{
	vec.clear();
	int posp = 0, posn = 0;

	while (posn < s.length())
	{
		posp = posn;
		//decide posp
		while ((s[posp] == ' ' || s[posp] == '\n') && posp < s.length())
		{
			posp++;
		}

		//decide posn:
		posn = posp;

		while ((s[posn] != ' ' && s[posn] != '\n') && posn < s.length())
		{
			posn++;
		}
		if (posn < s.length() && posn != posp)
		{
			double tmp = stod(s.substr(posp, posn - posp));
			vec.push_back(tmp);
		}

	}
}


int main(int argc, char** argv)
{
	if (argc < 4)
	{
		cout << "Usage: " << argv[0] << " input.vtu output.vtk outputchartlabel.txt" << endl;
		return 0;
	}
	
	std::cout << argv[1] << std::endl;
	std::ifstream in(argv[1]);

	std::string text((std::istreambuf_iterator<char>(in)),
		std::istreambuf_iterator<char>());

	//std::cout << text << std::endl;

	//string text;
	char *cstr = new char[text.length() + 1];
	strcpy(cstr, text.c_str());

	xml_document<>  doc;
	doc.parse<0>(cstr);

	cout << "Name of my first node is: " << doc.first_node()->name() << "\n";
	xml_node<> *vtkfile = doc.first_node("VTKFile");
	/*cout << "Node foobar has value " << node->value() << "\n";
	for (xml_attribute<> *attr = node->first_attribute();
		attr; attr = attr->next_attribute())
	{
		cout << "Node foobar has attribute " << attr->name() << " ";
		cout << "with value " << attr->value() << "\n";
	}*/

	
	xml_node<> *celldata = vtkfile->first_node("UnstructuredGrid")->first_node("Piece")->first_node("CellData");

	map<string, xml_node<> *> name2node;
	
	for (auto it = celldata->first_node(); it; it = it->next_sibling())
	{
		cout << "Name: " << it->first_attribute("Name")->value() << std::endl;
		name2node[it->first_attribute("Name")->value()] = it;
	}
	
	cout << "map size: " << name2node.size() << std::endl;
	
	xml_node<> *tetsel = NULL, *color = NULL, *color_dir = NULL;
	
	std::map<string, xml_node<> *>::iterator mapit = name2node.find("tet_selector");
	assert(mapit != name2node.end());
	if (mapit != name2node.end())
		tetsel = mapit->second;
		//std::cout << "found" << std::endl;
	
	mapit = name2node.find("color");
	assert(mapit != name2node.end());
	if (mapit != name2node.end())
		color = mapit->second;

	mapit = name2node.find("color_dir");
	assert(mapit != name2node.end());
	if (mapit != name2node.end())
		color_dir = mapit->second;

	//convert string to vector
	vector<int> v_tetsel, v_color, v_color_dir;
	
	string str_tetsel(tetsel->value()), str_color(color->value()), str_color_dir(color_dir->value());
	//while (pos < str_trisel.length())
	string2vec_int(str_tetsel, v_tetsel);
	string2vec_int(str_color, v_color);
	string2vec_int(str_color_dir, v_color_dir);

	assert(v_tetsel.size() == v_color.size() && v_color.size() == v_color_dir.size());
	std::cout << "Cell Size: " << v_color.size() << std::endl;

	//PointData
	xml_node<> *pointdata = vtkfile->first_node("UnstructuredGrid")->first_node("Piece")->first_node("PointData");
	for (auto it = pointdata->first_node(); it; it = it->next_sibling())
	{
		cout << "Name: " << it->first_attribute("Name")->value() << std::endl;
		name2node[it->first_attribute("Name")->value()] = it;
	}
	
	

	//point part
	xml_node<> *points = vtkfile->first_node("UnstructuredGrid")->first_node("Piece")->first_node("Points")->first_node("DataArray");
	vector<double> v_points;
	string str_points(points->value());

	//str changed to para here
	//para
	xml_node<> *params = NULL;
	mapit = name2node.find("param");
	if (mapit != name2node.end())
	{
		params = mapit->second;
		str_points = params->value();
	}


	string2vec_double(str_points, v_points);
	std::cout << "Point Array Size: " << v_points.size() << std::endl;

	//cells part
	xml_node<> *cells = vtkfile->first_node("UnstructuredGrid")->first_node("Piece")->first_node("Cells");

	for (auto it = cells->first_node(); it; it = it->next_sibling())
	{
		cout << "Name: " << it->first_attribute("Name")->value() << std::endl;
		name2node[it->first_attribute("Name")->value()] = it;
	}

	xml_node<> *connectivity = NULL, *offsets = NULL, *types = NULL;
	mapit = name2node.find("connectivity");
	assert(mapit != name2node.end());
	if (mapit != name2node.end())
		connectivity = mapit->second;
	mapit = name2node.find("offsets");
	assert(mapit != name2node.end());
	if (mapit != name2node.end())
		offsets = mapit->second;
	mapit = name2node.find("types");
	assert(mapit != name2node.end());
	if (mapit != name2node.end())
		types = mapit->second;

	string str_connectivity(connectivity->value());
	vector<int> v_connectivity;
	string2vec_int(str_connectivity, v_connectivity);
	std::cout << "Connectivity Size: " << v_connectivity.size() << std::endl;

	//tet might in front or back
	int tet_size = 0;
	for (size_t i = 0; i < v_tetsel.size(); i++)
	{
		if (v_tetsel[i] == 1) tet_size++;
	}

	TetStructure<double> tet_mesh;
	
	vector<CVec<double, 3>> points_cvec;
	CVec<double, 3> onepoint;
	for (size_t i = 0; i < v_points.size() / 3; i++)
	{
		onepoint[0] = v_points[3 * i + 0];
		onepoint[1] = v_points[3 * i + 1];
		onepoint[2] = v_points[3 * i + 2];
		points_cvec.push_back(onepoint);
	}

	vector<unsigned int> v_connectivity_cut;

	/*for (size_t i = 0; i < 4 * tet_size; i++)
	{
		v_connectivity_cut.push_back((unsigned int)v_connectivity[i]);
	}*/
	
	int conid = 0;
	for (size_t i = 0; i < v_tetsel.size(); i++)
	{
		if (v_tetsel[i] == 1)
		{
			v_connectivity_cut.push_back(v_connectivity[conid++]);
			v_connectivity_cut.push_back(v_connectivity[conid++]);
			v_connectivity_cut.push_back(v_connectivity[conid++]);
			v_connectivity_cut.push_back(v_connectivity[conid++]);
		}
		else
		{
			conid = conid + 3;
		}
	}
	

	tet_mesh.load_tet(points_cvec, v_connectivity_cut);
	tet_mesh.save_vtk(argv[2]);

	std::map<std::vector<int>, int> face2cellidx;
	
	std::vector<double> bf_pts;
	std::vector<int> bf_idx;
	std::vector<int> vert_s2v;
	
	get_boundary_vert_face(&tet_mesh, bf_pts, bf_idx, vert_s2v, true, &face2cellidx);

	//get label info
	std::vector<int> cell_label(tet_size, -1);
	
	int bf_size = (int)v_color.size() - tet_size;
	vector<int> tmp_face;

	conid = 0;
	for (size_t i = 0; i < v_tetsel.size(); i++)
	{
		if (v_tetsel[i] == 1)
		{
			conid = conid + 4;
		}
		else
		{
			tmp_face.clear();
			tmp_face.push_back(v_connectivity[conid++]);
			tmp_face.push_back(v_connectivity[conid++]);
			tmp_face.push_back(v_connectivity[conid++]);

			sort(tmp_face.begin(), tmp_face.end());
			auto it = face2cellidx.find(tmp_face);
			assert(it != face2cellidx.end());
			int cellid = it->second;
			int label = 2 * v_color[i] + v_color_dir[i];
			cell_label[cellid] = label;
		}
	}

	//cluster to form chart
	std::vector<std::vector<int>> v2c(v_points.size() / 3);
	std::vector<std::vector<int>> bcc(tet_size);
	const std::vector<Tetrahedron<double>*> tetras = tet_mesh.get_tetras();
	assert(tetras.size() == tet_size);

	//construct v2c
	for (size_t i = 0; i < tet_size; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			v2c[tetras[i]->vertex[j]->id].push_back((int)i);
		}
	}


	for (int i = 0; i < tet_size; i++)
	{
		if (tetras[i]->is_onboundary())
		{
			for (size_t j = 0; j < 4; j++)
			{
				if (tetras[i]->vertex[j]->boundary)
				{
					int vid = (int)tetras[i]->vertex[j]->id;
					for (size_t k = 0; k < v2c[vid].size(); k++)
					{
						int cid = v2c[vid][k];
						if (cid != i && tetras[cid]->is_onboundary())
						{
							bcc[i].push_back(cid);
						}
					}
				}
			}
		}
	}
	
	//clusterpart
	set<int> allbc;
	for (int i = 0; i < bcc.size(); i++)
	{
		if (bcc[i].size() > 0)
		{
			assert(cell_label[i] != -1);
			allbc.insert(i);
		}
	}
	
	//vector<vector<int>> chart_array;
	vector<int> cell_chart(tet_size, -1);
	int chart_count = 0;
	while (!allbc.empty())
	{
		vector<int> onechart;
		//init onechart
		auto it = allbc.begin();
		queue<int> q;
		q.push(*it);
		vector<int> cell_color(bcc.size(), 0);
		cell_color[*it] = 1;
		int cur_label = cell_label[*it];
		
		while (!q.empty())
		{
			int front = q.front();
			onechart.push_back(front);
			q.pop();
			for (auto n : bcc[front])
			{
				if (cell_color[n] != 1 && cell_label[n] == cur_label)
				{
					//not colored
					cell_color[n] = 1;
					q.push(n);
				}
			}
		}
		
		//delete elements from onechart
		for (auto e : onechart)
		{
			cell_chart[e] = chart_count;
			allbc.erase(e);
		}
		chart_count++;
	}


	//of stream
	std::ofstream ofs(argv[3]);

	for (size_t i = 0; i < tet_size; i++)
	{
		ofs << i << " " << cell_chart[i] << " " << cell_label[i] << endl;
	}
	
	ofs.close();
	in.close();
	//ofs << doc;
	/*std::string s;
	print(std::back_inserter(s), doc, 0);*/


	//system("pause");

	return 1;
}