#include <iostream>
#include <fstream>
#include <string>

// Include the file manager header
#include <OpenVolumeMesh/FileManager/FileManager.hh>
// Include the polyhedral mesh header
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

//do not consider cell order, might need to be changed
int main(int argc, char** argv)
{
	//input must be two file name
	if (argc != 3)
	{
		std::cout << "input format: " << "input.ovm output.vtk" << std::endl;
	}
	OpenVolumeMesh::GeometricPolyhedralMeshV3d myMesh;
	// Create file manager object
	OpenVolumeMesh::IO::FileManager fileManager;
	// Read mesh from file "myMesh.ovm" in the current directory
	fileManager.readFile(argv[1], myMesh);

	/*std::ifstream inputfile(argv[1]);
	std::ofstream outputfile(argv[2]);*/
	std::vector<double> x, y, z;
	std::vector<std::vector<int>> indices;
	
	
	for (OpenVolumeMesh::VertexIter v_it = myMesh.vertices_begin(); v_it != myMesh.vertices_end(); ++v_it)
	{
		int v_id = v_it->idx();
		OpenVolumeMesh::Geometry::Vec3d p = myMesh.vertex(*v_it);
		//dpx[v_id] = p[0]; dpy[v_id] = p[1]; dpz[v_id] = p[2];
		x.push_back(p[0]);
		y.push_back(p[1]);
		z.push_back(p[2]);
	}

	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_vec;

	for (OpenVolumeMesh::CellIter c_it = myMesh.cells_begin(); c_it != myMesh.cells_end(); ++c_it)
	{
		double cv_count = 0.0; int c_id = c_it->idx();
		std::vector<int> tmp_indices;
		const std::vector<OpenVolumeMesh::VertexHandle>& vertices_ = myMesh.cell(*c_it).vertices();
		hff_vec = myMesh.cell(*c_it).halffaces();
		assert(hff_vec.size() == 4);

		OpenVolumeMesh::HalfFaceVertexIter hfv_it = myMesh.hfv_iter(hff_vec[0]);
		//unsigned int i = 0;
		std::set<int> face_set;
		for (hfv_it; hfv_it; ++hfv_it)
		{
			tmp_indices.push_back(hfv_it->idx());
			face_set.insert(hfv_it->idx());
		}

		for (unsigned i = 0; i < vertices_.size(); ++i)
		{
			int v_id = vertices_[i];

			auto it = face_set.find(v_id);
			if (it == face_set.end())
			{
				tmp_indices.push_back(v_id);
			}
		}
		assert(tmp_indices.size() == 4);
		indices.push_back(tmp_indices);
	}

	std::ofstream outputfile(argv[2]);
	outputfile << "# vtk DataFile Version 3.0\n"
		<< "Volume mesh\n"
		<< "ASCII\n"
		<< "DATASET UNSTRUCTURED_GRID\n";

	int n_point = x.size();
	int n_cell = indices.size();
	outputfile << "POINTS " << n_point << " double" << std::endl;

	for (size_t i = 0; i < n_point; i++)
	{
		outputfile << x[i] << " " << y[i] << " " << z[i] << std::endl;
	}
	outputfile << "CELLS " << n_cell << " " << 5 * n_cell << std::endl;
	for (size_t i = 0; i < n_cell; i++)
	{
		outputfile << "4 ";
		for (size_t j = 0; j < 4; j++)
		{
			outputfile << indices[i][j] << " ";
		}
		outputfile << std::endl;
	}

	outputfile << "CELL_TYPES " << n_cell << std::endl;
	for (size_t i = 0; i < n_cell; i++)
	{
		outputfile << "10" << std::endl;
	}

	outputfile.close();
	
}