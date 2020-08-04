#include <iostream>
#include <fstream>
#include <map>
#include "TetStructure.h"

using ig::TetStructure;
using ig::CVec;
using ig::Tetrahedron;
using ig::TetVertex;

void get_boundary_vert_face(TetStructure<double>* tet_mesh, std::vector<double> &points, std::vector<int> &faces, std::vector<int> &s2v, bool ori_vert_flag = false, std::map<std::vector<int>, int> *face2cellidx = NULL)
{
	//s2v: surface to volume vert map
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
				map_it = VolumeSurfaceVertex.find((const int)tetras[i]->vertex[s_tet_id[face_id][j]]->id);
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

int main(int argc, char** argv)
{
	if (argc < 3)
	{
		std::cout << "input format: " << "input.vtk output.obj ori_vert_flag [OPTIONAL]input_volume_chartlabel.txt [OPTIONAL]output_surface_label.txt" << std::endl;
		return 1;
	}

	bool ori_vert_flag = false;
	if (argc > 3)
	{
		//redefine ori_vert_flag
		ori_vert_flag = atoi(argv[3]);
	}
	
	TetStructure<double> tet_mesh;
	tet_mesh.load_vtk(argv[1]);
	std::cout << "input file: " << argv[1] << " output file: " << argv[2] << std::endl;
	std::cout << "n vert: " << tet_mesh.tetra_vertices.size() << std::endl;

	std::vector<double> points; 
	std::vector<int> faces; 
	std::vector<int> s2v; //surface vertex to volume vertex
	
	std::map<std::vector<int>, int> face2cell;
	
	get_boundary_vert_face(&tet_mesh, points, faces, s2v, ori_vert_flag);
	/*if (argc < 5)
		get_boundary_vert_face(&tet_mesh, points, faces, s2v, ori_vert_flag);
	else
		get_boundary_vert_face(&tet_mesh, points, faces, s2v, ori_vert_flag, &face2cell);*/


	std::ofstream outputfile(argv[2]);
	outputfile << "g object" << std::endl;
	int n_point = (int)points.size() / 3;
	int n_face = (int)faces.size() / 3;
	for (size_t i = 0; i < n_point; i++)
	{
		outputfile << "v " << points[3 * i] << " " << points[3 * i + 1] << " " << points[3 * i + 2] << std::endl;
	}
	
	for (size_t i = 0; i < n_face; i++)
	{
		outputfile << "f " << faces[3 * i] + 1 << " " << faces[3 * i + 1] + 1 << " " << faces[3 * i + 2] + 1 << std::endl;
	}


	outputfile.close();

	if (argc == 6)
	{
		if (ori_vert_flag == 0) return 1; //ori vert must be preserved
		std::ifstream input_chartlabel(argv[4]);
		//load chart label
		std::vector<int> face_label;
		int id, chart, label;
		std::vector<int> face2volume;
		for (size_t i = 0; i < tet_mesh.tetras.size(); i++)
		{
			input_chartlabel >> id >> chart >> label;
			if (!tet_mesh.tetras[i]->is_onboundary()) continue;
			assert(label >= 0);
			face_label.push_back(label);
			face2volume.push_back((int)i);
		}
		
		input_chartlabel.close();
		assert(face_label.size() * 3 == faces.size());

		std::ofstream output_label(argv[5]);
		
		output_label << face_label.size() << std::endl;
		for (size_t i = 0; i < face_label.size(); i++)
		{
			output_label << i << " " << face_label[i] << std::endl;
		}
		output_label << tet_mesh.tetras.size() << std::endl;
		for (size_t i = 0; i < face2volume.size(); i++)
		{
			output_label << i << " " << face2volume[i] << std::endl;
		}
		

		output_label.close();

	}


	return 1;
}