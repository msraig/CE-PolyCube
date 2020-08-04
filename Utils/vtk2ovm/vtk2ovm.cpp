#include <OpenVolumeMesh/FileManager/FileManager.hh>
//#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "TetStructure.h"
//#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

using namespace OpenVolumeMesh;
using OpenVolumeMesh::Geometry::Vec3d;
using namespace std;
using ig::TetStructure;
using ig::CVec;
using ig::Tetrahedron;
using ig::TetVertex;


int tag_check(const vector<int> &v1, const vector<int> &v2)
{
	//if v1 and v2 follow the same order, return 0
	//else return 1;
	int tag = 0;
	int id0 = 0;
	for (; id0 < 3; ++id0)
	{
		if (v2[id0] == v1[0])
			break;
	}

	int id1 = (id0 + 1) % 3;
	if (v2[id1] != v1[1])
		tag = 1;
	return tag;

}


int main(int argc, char** argv) {

	if (argc != 3)
	{
		std::cout << "Input Format: " << argv[0] << " input.vtk output.ovm" << std::endl;
		return 0;
	}

	std::string input = argv[1];
	std::string output = argv[2];

	//read vtk file
	TetStructure<double> tet_mesh;
	tet_mesh.load_vtk(argv[1]);
	int n_cells = (int)tet_mesh.tetras.size();
	vector<vector<int>> cells(n_cells);
	for (size_t i = 0; i < n_cells; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			cells[i].push_back((int)tet_mesh.tetras[i]->vertex[j]->id);
		}
	}
	
	int n_points = (int)tet_mesh.tetra_vertices.size();

	GeometricPolyhedralMeshV3d ovm;
	//add vertices
	for (size_t i = 0; i < n_points; i++)
	{
		CVec<double, 3> p = tet_mesh.tetra_vertices[i]->pos;
		ovm.add_vertex(Vec3d(p[0], p[1], p[2]));
	}

	//construct face;
	vector<OpenVolumeMesh::VertexHandle> vertices;
	vector<vector<int>> tet2faceid(n_cells, vector<int>(4, -1));
	vector<vector<int>> tet2facetag(n_cells, vector<int>(4, -1));
	//vector<int> face_tags; //0 for same direction and 1 for opposite, corresponds to tet2face

	vector<vector<int>> total_fs(n_cells * 4);
	vector<tuple<int, int, int, int, int, int, int>> tempF(n_cells * 4);
	vector<int> oneface(3, -1);
	
	//to interior
	//to interior, same order : tag 0, df order : tag 1.
	static int f2v[4][3] = { {2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2} };

	//acceleration

	for (size_t i = 0; i < n_cells; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			for (size_t k = 0; k < 3; k++)
			{
				oneface[k] = cells[i][f2v[j][k]];
			}
			int id = (int)(4 * i + j);
			total_fs[id] = oneface;
			std::sort(oneface.begin(), oneface.end());
			//tag
			int tag = 0;
			vector<int> tmp = total_fs[id];
			tag = tag_check(total_fs[id], oneface);
			tempF[id] = std::make_tuple(oneface[0], oneface[1], oneface[2], id, i, j, tag);
		}
	}
	
	std::sort(tempF.begin(), tempF.end());
	
	vector<vector<int>> face2tet;
	int F_num = 0;
	for (size_t i = 0; i < tempF.size(); i++)
	{
		if (i == 0 || (i != 0 &&
			(std::get<0>(tempF[i]) != std::get<0>(tempF[i - 1]) || std::get<1>(tempF[i]) != std::get<1>(tempF[i - 1]) ||
				std::get<2>(tempF[i]) != std::get<2>(tempF[i - 1])))) {
			//first face || (at least one vertex df from the former)
			//new face created
			F_num++;
			vector<int> tmpface = total_fs[std::get<3>(tempF[i])];
			sort(tmpface.begin(), tmpface.end());
			vertices.clear();
			for (size_t j = 0; j < 3; j++)
			{
				vertices.push_back(VertexHandle(tmpface[j]));
			}
			ovm.add_face(vertices);
			
			tet2faceid[std::get<4>(tempF[i])][std::get<5>(tempF[i])] = F_num - 1;
			tet2facetag[std::get<4>(tempF[i])][std::get<5>(tempF[i])] = std::get<6>(tempF[i]);

		}
		else if (i != 0 && (std::get<0>(tempF[i]) == std::get<0>(tempF[i - 1]) && std::get<1>(tempF[i]) == std::get<1>(tempF[i - 1]) &&
			std::get<2>(tempF[i]) == std::get<2>(tempF[i - 1])))
		{
			tet2faceid[std::get<4>(tempF[i])][std::get<5>(tempF[i])] = F_num - 1;
			tet2facetag[std::get<4>(tempF[i])][std::get<5>(tempF[i])] = std::get<6>(tempF[i]);
		}
	}

	
	//construct cells
	std::vector<OpenVolumeMesh::HalfFaceHandle> halffaces;
	for (size_t i = 0; i < n_cells; i++)
	{
		halffaces.clear();
		for (size_t j = 0; j < 4; j++)
		{
			assert(tet2faceid[i][j] != -1 && tet2facetag[i][j] != -1);
			halffaces.push_back(ovm.halfface_handle(FaceHandle(tet2faceid[i][j]), tet2facetag[i][j]));
		}
		ovm.add_cell(halffaces);
	}

	OpenVolumeMesh::IO::FileManager fileManager;
	std::string filename(output);

	fileManager.writeFile(filename, ovm);
	std::cout << "write file successfully" << std::endl;
	//system("pause");

	return 0;
}