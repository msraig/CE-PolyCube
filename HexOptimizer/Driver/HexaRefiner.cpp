#include "HexaRefine.h"
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
	int max_iter = 20;
	double threshold = 0;
	bool fixbnd = false, surf_preserve = false, use_extend_jac = false, opt_illegal_only = false;
	try
	{
		cxxopts::Options options("HexaRefiner", "Hexahedral mesh smoothing (author: Yang Liu, Email: yangliu@microsoft.com)");
		options
			.positional_help("[optional args]")
			.show_positional_help()
			.allow_unrecognised_options()
			.add_options()
			("i,input", "input hexahedral mesh (hex/vtk/ovm/mesh format)", cxxopts::value<std::string>())
			("o,output", "output hexahedral mesh (hex/vtk/ovm/mesh format)", cxxopts::value<std::string>())
			("b,bnd", "boundary mesh (obj/off format)", cxxopts::value<std::string>())
			("n,iter", "max num of iteration (default: 20)", cxxopts::value<int>())
			("t,threshold", "threshold value (0~0.99, default: 0)", cxxopts::value<double>())
			("f,fixbnd", "fix boundary vertices (default: false)", cxxopts::value<bool>())
			("p,surf", "surface preserving (default: true)", cxxopts::value<bool>())
			("k,bf", "boundary feature (tfe file)", cxxopts::value<std::string>())
			("l,hf", "hex feature (hfe file)", cxxopts::value<std::string>())
			("d,dump", "dump features for HexaSmooth (default: false)", cxxopts::value<bool>())
			("e,ejac", "extended jacobian (default: false)", cxxopts::value<bool>())
			("r,restrict", "restrict optimizaton on illege cells only (default: false)", cxxopts::value<bool>())
			("s,nosmooth", "disable Laplacian smoothing (default:false)", cxxopts::value<bool>())
			("h,help", "print help");

		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			std::cout << options.help({ "", "Group" }) << std::endl;
			exit(0);
		}

		if (result.count("i"))
		{
			auto& inputfile = result["i"].as<std::string>();
			std::string inputext = GetFileExtension(inputfile);

			std::string outputfile = inputfile.substr(0, inputfile.find_last_of('.')) + "_opt." + inputext;
			if (result.count("o"))
			{
				outputfile = result["o"].as<std::string>();
			}

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

			if (result.count("n"))
			{
				max_iter = std::max(1, result["n"].as<int>());
			}

			Mesh3d bndmesh;
			bool fileloadOK = false;
			if (result.count("b"))
			{
				auto& bndmeshfile = result["b"].as<std::string>();
				std::string meshext = GetFileExtension(bndmeshfile);
				if (meshext == "obj")
				{
					fileloadOK = bndmesh.load_obj(bndmeshfile.c_str());
				}
				else if (meshext == "off")
				{
					fileloadOK = bndmesh.load_off(bndmeshfile.c_str());
				}
				else
				{
					std::cout << "Unrecogonized mesh format: " << bndmeshfile << " !" << std::endl;
				}
			}

			Mesh3d* bmesh = 0;
			if (fileloadOK)
				bmesh = &bndmesh;

			if (result.count("t"))
				threshold = std::min(std::max(0.0, result["t"].as<double>()), 0.99);

			if (result.count("f"))
				fixbnd = result["f"].as<bool>();

			if (result.count("e"))
			{
				use_extend_jac = result["e"].as<bool>();
			}

			if (result.count("p"))
			{
				surf_preserve = true;
			}
			if (result.count("r"))
				opt_illegal_only = true;

			MeshLib::HexaRefine<double> m_smooth(&m_hexmesh, bmesh, max_iter, fixbnd, surf_preserve, threshold, use_extend_jac, opt_illegal_only);

			if (result.count("k") && result.count("l"))
			{
				auto& trifeature = result["k"].as<std::string>();
				auto& hexfeature = result["l"].as<std::string>();
				if (!m_smooth.load_features(trifeature.c_str(), hexfeature.c_str()))
				{
					std::cout << "the feature numbers do not match!" << std::endl;
					return 0;
				}
			}

			if (result.count("d"))
			{
				auto& bndmeshfile = result["b"].as<std::string>();
				std::string filename = bndmeshfile.substr(0, bndmeshfile.find_last_of('.'));
				std::string objname = filename + "_tri.obj";
				std::string trifeaname = filename + "_tri.fea";
				std::string quadname = filename + "_quad.obj";
				std::string hexfeaname = filename + "_quad.fea";
				std::string trifeaname2 = filename + "_tri.polyline";
				std::string hexfeaname2 = filename + "_hex.polyline";

				m_smooth.save_quadmesh(quadname.c_str());
				m_smooth.save_trimesh(objname.c_str());
				m_smooth.save_trifeature(trifeaname.c_str());
				m_smooth.save_quadfeature(hexfeaname.c_str());
				m_smooth.save_trifeature_polyline(trifeaname2.c_str());
				m_smooth.save_hexfeature_polyline(hexfeaname2.c_str());
				//m_smooth.dump(objname.c_str(), trifeaname.c_str(), quadname.c_str(), hexfeaname.c_str());
			}

			bool laplacian_smoothing = result.count("s")? false : true;
			
			m_smooth.smooth(laplacian_smoothing, true);

			if (result.count("d"))
			{
				auto& bndmeshfile = result["b"].as<std::string>();
				std::string filename = bndmeshfile.substr(0, bndmeshfile.find_last_of('.'));
				std::string quadname = filename + "_quad_final.obj";
				m_smooth.save_quadmesh(quadname.c_str());
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