#include "VHexMesh.h"
#include "cxxopts.hpp"

std::string GetFileExtension(const std::string& FileName)
{
	if (FileName.find_last_of(".") != std::string::npos)
		return FileName.substr(FileName.find_last_of(".") + 1);
	return "";
}

int main(int argc, char** argv)
{
	using namespace MeshLib;

	try
	{
		cxxopts::Options options("HexaConvert", "Hexahedral mesh format converter (author: Yang Liu, Email: yangliu@microsoft.com)");
		options
			.positional_help("[optional args]")
			.show_positional_help()
			.allow_unrecognised_options()
			.add_options()
			("i,input", "input hexahedral mesh (hex/vtk/ovm/mesh format)", cxxopts::value<std::string>())
			("o,output", "output hexahedral mesh (hex/vtk/ovm/mesh format)", cxxopts::value<std::string>())
			("m,merge", "merge closest vertices", cxxopts::value<bool>())
			("p,closeness", "closeness threshold (default 1.0e-8)", cxxopts::value<double>())
			("s,subdiv", "subdivision", cxxopts::value<int>())
			("g,padding", "global padding", cxxopts::value<bool>())
			("h,help", "print help");

		bool is_merge = false, padding = false;
		double closeness = 1.0e-8;

		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			std::cout << options.help({ "", "Group" }) << std::endl;
			exit(0);
		}

		if (result.count("i") && result.count("o"))
		{
			auto& inputfile = result["i"].as<std::string>();
			auto& outputfile = result["o"].as<std::string>();

			std::string inputext = GetFileExtension(inputfile);
			std::string outputext = GetFileExtension(outputfile);

			VHexMesh3d m_hexmesh;

			if (inputext == "hex")
			{
				m_hexmesh.load_hex(inputfile.c_str());
			}
			else if (inputext == "ovm")
			{
				m_hexmesh.load_ovm(inputfile.c_str());
			}
			else if (inputext == "vtk")
			{
				m_hexmesh.load_vtk(inputfile.c_str());
			}
			else if (inputext == "mesh")
			{
				m_hexmesh.load_mesh(inputfile.c_str());
			}
			else
			{
				throw(cxxopts::OptionException("the input file does not exist!"));
			}

			if (result.count("m"))
				is_merge = result["m"].as<bool>();
			if (is_merge)
			{
				if (result.count("p"))
				{
					closeness = std::abs(result["p"].as<double>());
				}
				size_t nv = m_hexmesh.get_vertices().size();
				m_hexmesh.merge_vertices(closeness);
				size_t dnv = nv - m_hexmesh.get_vertices().size();
				std::cout << dnv << " vertices have been merged!" << std::endl;
			}

			if (result.count("s"))
			{
				int num = result["s"].as<int>();
				if (num > 0) std::cout << "subdivide " << num << " times" << std::endl;
				for (int i = 0; i < num; i++)
					m_hexmesh.subdivide();
			}

			if (result.count("g"))
			{
				padding = result["g"].as<bool>();
				if (padding)
					m_hexmesh.padding();
			}

			if (outputext == "hex")
			{
				m_hexmesh.save_hex(outputfile.c_str());
			}
			else if (outputext == "ovm")
			{
				m_hexmesh.save_ovm(outputfile.c_str());
			}
			else if (outputext == "vtk")
			{
				m_hexmesh.save_vtk(outputfile.c_str());
			}
			else if (outputext == "mesh")
			{
				m_hexmesh.save_mesh(outputfile.c_str());
			}
			else
			{
				throw(cxxopts::OptionException("the output file extension is not valid!"));
			}
		}
		else
		{
			std::cout << "Please check your input!" << std::endl;
		}
	}
	catch (const cxxopts::OptionException & e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}

	return 0;
}