#include <iostream>
#include <fstream>
#include <string>
#include <vector>

//hex to vtk

int main(int argc, char** argv)
{
	//input must be two file name'
	//tet version
	if (argc != 3)
	{
		std::cout << "input format: " << "input.vtk output.mesh" << std::endl;
	}

	std::ifstream inputfile(argv[1]);
	std::ofstream outputfile(argv[2]);
	std::string str;
	unsigned int n_point, n_tet;
	std::vector<double> coord_x, coord_y, coord_z;

	do
	{
		inputfile >> str;
	} while (str != "POINTS");


	inputfile >> n_point >> str;
	//hexfile >> n_hex >> str;
	//std::cout << "n point: " << n_point << std::endl;

	coord_x.resize(n_point);
	coord_y.resize(n_point);
	coord_z.resize(n_point);

	double x, y, z;
	for (size_t i = 0; i < n_point; i++)
	{
		inputfile >> x >> y >> z;
		//pts.push_back(ig::CVec<double, 3>(x, y, z));
		coord_x[i] = x;
		coord_y[i] = y;
		coord_z[i] = z;
	}

	do
	{
		inputfile >> str;
	} while (str != "CELLS");

	int hex_type;
	inputfile >> n_tet >> hex_type;

	outputfile << "MeshVersionFormatted 1" << std::endl;
	outputfile << "Dimension 3" << std::endl;
	outputfile << "Vertices " << std::endl;
	outputfile << n_point << std::endl;

	for (size_t i = 0; i < n_point; i++)
	{
		outputfile << coord_x[i] << " " << coord_y[i] << " " << coord_z[i] << " 0" << std::endl;
	}

	outputfile << "Tetrahedra " << std::endl;
	outputfile << n_tet << std::endl;

	//int vtk_id[] = { 2, 3, 6, 7, 1, 0, 5, 4 };
	int vtk_id[] = { 0,1,2,3,4,5,6,7 };

	for (size_t i = 0; i < n_tet; i++)
	{
		int tmp;
		std::vector<int> id(4, 0);
		inputfile >> tmp;
		for (size_t j = 0; j < 4; j++)
		{
			inputfile >> id[j];
		}
		//outputfile << "8";
		for (size_t j = 0; j < 4; j++)
		{
			outputfile << id[vtk_id[j]] + 1 << " ";
		}
		outputfile << "0" << std::endl;
	}
	outputfile << "End" << std::endl;




	inputfile.close();
	outputfile.close();

	return 1;
}