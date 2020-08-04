#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <queue>
#include "TetStructure.h"
#include "TetOperation.h"
#include "TetStructure.cpp"

using ig::TetStructure;
using ig::CVec;
using ig::Tetrahedron;
using ig::TetVertex;

using ig::TetOperation;



int main(int argc, char** argv)
{
	if (argc != 3)
	{
		std::cout << "input format: " << "SplitTet input.vtk output.vtk" << std::endl;
	}

	TetStructure<double> tet_mesh;
	tet_mesh.load_vtk(argv[1]);
	std::cout << "input file: " << argv[1] << " output file: " << argv[2] << std::endl;
	std::cout << "n vert: " << tet_mesh.tetra_vertices.size() << std::endl;
	std::cout << "n cell: " << tet_mesh.tetras.size() << std::endl;

	TetOperation<double> tet_op;
	tet_op.set_tet_structure(&tet_mesh);
	tet_op.split_inner_edge_two_boundary_point(tet_mesh);
	tet_op.split_cell_four_boundary_point(tet_mesh);
	tet_mesh.save_vtk(argv[2]);


	return 1;
}