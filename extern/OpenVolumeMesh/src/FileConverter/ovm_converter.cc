#include <fstream>
#include <sstream>

#include <boost/spirit/include/qi.hpp>

#include <OpenVolumeMesh/FileManager/FileManager.hh>

#include "MeshGenerator.hpp"

#include "Grammars/Tetmesh.hpp"
#include "Grammars/Netgen.hpp"

void printUsage() {
    std::cerr << "You need to specify a format and a source file to convert!" << std::endl << std::endl;
    std::clog << "Usage: file_converter <format> <filename> [output_filename.ovm]" << std::endl << std::endl;
    std::clog << "Available file formats:" << std::endl;
    std::clog << "  -t\tTetmesh" << std::endl;
    std::clog << "  -n\tNetgen" << std::endl;
    std::clog << std::endl;
}

int main(int _argc, char* _argv[]) {

    if(_argc < 3 || _argc > 4 ||
            (_argc > 1 && (std::strcmp(_argv[1], "--help") == 0 || std::strcmp(_argv[1], "-h") == 0))) {
        printUsage();
        return -1;
    }

    std::ifstream iff(_argv[2], std::ios::in);

    if(!iff.good()) {
        std::cerr << "Could not open file " << _argv[1] << " for reading!" << std::endl;
        return -1;
    }

    MeshGenerator::PolyhedralMesh mesh;
    MeshGenerator generator(mesh);

    std::ostringstream oss;
    oss << iff.rdbuf();

    std::string fileContent = oss.str();

    // Instantiate grammar objects
    tetmesh_grammar<std::string::iterator> tetmesh_grammar(generator);
    netgen_grammar<std::string::iterator> netgen_grammar(generator);

    bool r = false;

    if(std::strcmp(_argv[1], "-t") == 0) {

        // Tetmesh format
        r = boost::spirit::qi::phrase_parse(fileContent.begin(), fileContent.end(), tetmesh_grammar, qi::space);

    } else if(std::strcmp(_argv[1], "-n") == 0) {

        // Netgen format
        r = boost::spirit::qi::phrase_parse(fileContent.begin(), fileContent.end(), netgen_grammar, qi::space);

    } else {
        printUsage();
        return -1;
    }

    if(r) {
        std::cout << "Successfully read file data!" << std::endl;
    } else {
        std::cout << "Parsing failed!" << std::endl;
    }

    std::cerr << "Converted " << mesh.n_vertices() << " vertices," << std::endl;
    std::cerr << "\t  " << mesh.n_edges() << " edges," << std::endl;
    std::cerr << "\t  " << mesh.n_faces() << " faces," << std::endl;
    std::cerr << "\t  " << mesh.n_cells() << " cells!" << std::endl;

    OpenVolumeMesh::IO::FileManager fileManager;

    std::string filename;

    if(_argc == 3) {
        filename = _argv[2];
        std::string::size_type idx = filename.rfind('.');

        filename = filename.substr(0, idx);
        filename.append(".ovm");
    } else {
        filename = _argv[3];
    }

    // Write mesh to file
    fileManager.writeFile(filename.c_str(), mesh);

    return 0;
}
