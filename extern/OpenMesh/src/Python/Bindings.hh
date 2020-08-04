/** @file */

#ifndef OPENMESH_PYTHON_BINDINGS_HH
#define OPENMESH_PYTHON_BINDINGS_HH

#include <boost/python.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/reference_existing_object.hpp>
#include <boost/python/copy_const_reference.hpp>

#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include "OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh"

using namespace boost::python;

namespace OpenMesh {

/**
 * This namespace contains classes and functions that are used to expose
 * %OpenMesh to %Python.
 */
namespace Python {

/**
 * Return value policy for functions that return references to objects that are
 * managed by %OpenMesh.
 */
#define OPENMESH_PYTHON_DEFAULT_POLICY return_value_policy<copy_const_reference>()

struct MeshTraits : public OpenMesh::DefaultTraits {
	/** Use double precision points */
	typedef OpenMesh::Vec3d Point;

	/** Use double precision normals */
	typedef OpenMesh::Vec3d Normal;

	/** Use RGBA colors */
	typedef OpenMesh::Vec4f Color;
};

typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> TriMesh;
typedef OpenMesh::PolyMesh_ArrayKernelT<MeshTraits> PolyMesh;

} // namespace OpenMesh
} // namespace Python

#endif
