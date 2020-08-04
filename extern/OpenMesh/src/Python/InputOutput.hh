#ifndef OPENMESH_PYTHON_INPUTOUTPUT_HH
#define OPENMESH_PYTHON_INPUTOUTPUT_HH

#include "Python/Bindings.hh"

namespace OpenMesh {
namespace Python {

const IO::Options::Flag FLAG_DEFAULT        = IO::Options::Default;
const IO::Options::Flag FLAG_BINARY         = IO::Options::Binary;
const IO::Options::Flag FLAG_MSB            = IO::Options::MSB;
const IO::Options::Flag FLAG_LSB            = IO::Options::LSB;
const IO::Options::Flag FLAG_SWAP           = IO::Options::Swap;
const IO::Options::Flag FLAG_VERTEXNORMAL   = IO::Options::VertexNormal;
const IO::Options::Flag FLAG_VERTEXCOLOR    = IO::Options::VertexColor;
const IO::Options::Flag FLAG_VERTEXTEXCOORD = IO::Options::VertexTexCoord;
const IO::Options::Flag FLAG_EDGECOLOR      = IO::Options::EdgeColor;
const IO::Options::Flag FLAG_FACENORMAL     = IO::Options::FaceNormal;
const IO::Options::Flag FLAG_FACECOLOR      = IO::Options::FaceColor;
const IO::Options::Flag FLAG_FACETEXCOORD   = IO::Options::FaceTexCoord;
const IO::Options::Flag FLAG_COLORALPHA     = IO::Options::ColorAlpha;
const IO::Options::Flag FLAG_COLORFLOAT     = IO::Options::ColorFloat;

BOOST_PYTHON_FUNCTION_OVERLOADS(read_mesh_overloads, IO::read_mesh, 3, 4)
BOOST_PYTHON_FUNCTION_OVERLOADS(write_mesh_overloads, IO::write_mesh, 2, 4)

/**
 * Expose the input/output functions and options to Python.
 */
void expose_io() {

	//======================================================================
	//  Functions
	//======================================================================

	bool (*read_mesh_poly        )(PolyMesh&, const std::string&                    ) = &IO::read_mesh;
	bool (*read_mesh_poly_options)(PolyMesh&, const std::string&, IO::Options&, bool) = &IO::read_mesh;
	bool (*read_mesh_tri         )(TriMesh&,  const std::string&                    ) = &IO::read_mesh;
	bool (*read_mesh_tri_options )(TriMesh&,  const std::string&, IO::Options&, bool) = &IO::read_mesh;

	bool (*write_mesh_poly)(const PolyMesh&, const std::string&, IO::Options, std::streamsize) = &IO::write_mesh;
	bool (*write_mesh_tri )(const TriMesh&,  const std::string&, IO::Options, std::streamsize) = &IO::write_mesh;

	def("read_mesh", read_mesh_poly);
	def("read_mesh", read_mesh_poly_options, read_mesh_overloads());
	def("read_mesh", read_mesh_tri);
	def("read_mesh", read_mesh_tri_options, read_mesh_overloads());

	def("write_mesh", write_mesh_poly, write_mesh_overloads());
	def("write_mesh", write_mesh_tri, write_mesh_overloads());

	//======================================================================
	//  Options
	//======================================================================

	scope scope_options = class_<IO::Options>("Options")
		.def(init<IO::Options::Flag>())
		.def("cleanup", &IO::Options::cleanup)
		.def("clear", &IO::Options::clear)
		.def("is_empty", &IO::Options::is_empty)
		.def("check", &IO::Options::check)
		.def("is_binary", &IO::Options::is_binary)
		.def("vertex_has_normal", &IO::Options::vertex_has_normal)
		.def("vertex_has_color", &IO::Options::vertex_has_color)
		.def("vertex_has_texcoord", &IO::Options::vertex_has_texcoord)
		.def("edge_has_color", &IO::Options::edge_has_color)
		.def("face_has_normal", &IO::Options::face_has_normal)
		.def("face_has_color", &IO::Options::face_has_color)
		.def("face_has_texcoord", &IO::Options::face_has_texcoord)
		.def("color_has_alpha", &IO::Options::color_has_alpha)
		.def("color_is_float", &IO::Options::color_is_float)

		.def(self == self)
		.def(self != self)
		.def(self -= IO::Options::Flag())
		.def(self += IO::Options::Flag())

		.def_readonly("Default", &FLAG_DEFAULT)
		.def_readonly("Binary", &FLAG_BINARY)
		.def_readonly("MSB", &FLAG_MSB)
		.def_readonly("LSB", &FLAG_LSB)
		.def_readonly("Swap", &FLAG_SWAP)
		.def_readonly("VertexNormal", &FLAG_VERTEXNORMAL)
		.def_readonly("VertexColor", &FLAG_VERTEXCOLOR)
		.def_readonly("VertexTexCoord", &FLAG_VERTEXTEXCOORD)
		.def_readonly("EdgeColor", &FLAG_EDGECOLOR)
		.def_readonly("FaceNormal", &FLAG_FACENORMAL)
		.def_readonly("FaceColor", &FLAG_FACECOLOR)
		.def_readonly("FaceTexCoord", &FLAG_FACETEXCOORD)
		.def_readonly("ColorAlpha", &FLAG_COLORALPHA)
		.def_readonly("ColorFloat", &FLAG_COLORFLOAT)
		;

	enum_<IO::Options::Flag>("Flag")
		.value("Default", IO::Options::Default)
		.value("Binary", IO::Options::Binary)
		.value("MSB", IO::Options::MSB)
		.value("LSB", IO::Options::LSB)
		.value("Swap", IO::Options::Swap)
		.value("VertexNormal", IO::Options::VertexNormal)
		.value("VertexColor", IO::Options::VertexColor)
		.value("VertexTexCoord", IO::Options::VertexTexCoord)
		.value("EdgeColor", IO::Options::EdgeColor)
		.value("FaceNormal", IO::Options::FaceNormal)
		.value("FaceColor", IO::Options::FaceColor)
		.value("FaceTexCoord", IO::Options::FaceTexCoord)
		.value("ColorAlpha", IO::Options::ColorAlpha)
		.value("ColorFloat", IO::Options::ColorFloat)
		;
}

} // namespace OpenMesh
} // namespace Python

#endif
