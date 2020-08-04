#include "Python/Bindings.hh"
#include "Python/Vector.hh"
#include "Python/Mesh.hh"
#include "Python/Iterator.hh"
#include "Python/Circulator.hh"
#include "Python/PropertyManager.hh"
#include "Python/InputOutput.hh"

namespace OpenMesh {
namespace Python {

/**
 * Expose mesh items to %Python.
 */
void expose_items() {
	class_<ArrayItems::Vertex>("Vertex");
	class_<ArrayItems::Halfedge>("Halfedge");
	class_<ArrayItems::Edge>("Edge");
	class_<ArrayItems::Face>("Face");
}

/**
 * Expose item and property handles to %Python.
 */
void expose_handles() {
	class_<BaseHandle>("BaseHandle", init<optional<int> >())
		.def("idx", &BaseHandle::idx)
		.def("is_valid", &BaseHandle::is_valid)
		.def("reset", &BaseHandle::reset)
		.def("invalidate", &BaseHandle::invalidate)
		.def(self == self)
		.def(self != self)
		.def(self < self)
		;

	class_<VertexHandle, bases<BaseHandle> >("VertexHandle", init<optional<int> >());
	class_<HalfedgeHandle, bases<BaseHandle> >("HalfedgeHandle", init<optional<int> >());
	class_<EdgeHandle, bases<BaseHandle> >("EdgeHandle", init<optional<int> >());
	class_<FaceHandle, bases<BaseHandle> >("FaceHandle", init<optional<int> >());

	class_<BasePropHandleT<object>, bases<BaseHandle> >("BasePropHandle", init<optional<int> >());

	class_<VPropHandleT<object>, bases<BasePropHandleT<object> > >("VPropHandle", init<optional<int> >())
		.def(init<const BasePropHandleT<object>&>());
	class_<HPropHandleT<object>, bases<BasePropHandleT<object> > >("HPropHandle", init<optional<int> >())
		.def(init<const BasePropHandleT<object>&>());
	class_<EPropHandleT<object>, bases<BasePropHandleT<object> > >("EPropHandle", init<optional<int> >())
		.def(init<const BasePropHandleT<object>&>());
	class_<FPropHandleT<object>, bases<BasePropHandleT<object> > >("FPropHandle", init<optional<int> >())
		.def(init<const BasePropHandleT<object>&>());
	class_<MPropHandleT<object>, bases<BasePropHandleT<object> > >("MPropHandle", init<optional<int> >())
		.def(init<const BasePropHandleT<object>&>());
}


/**
 * Expose the StatusBits enum and StatusInfo class to %Python.
 */
void expose_status_bits_and_info() {
	using OpenMesh::Attributes::StatusBits;
	using OpenMesh::Attributes::StatusInfo;

	enum_<StatusBits>("StatusBits")
		.value("DELETED", OpenMesh::Attributes::DELETED)
		.value("LOCKED", OpenMesh::Attributes::LOCKED)
		.value("SELECTED", OpenMesh::Attributes::SELECTED)
		.value("HIDDEN", OpenMesh::Attributes::HIDDEN)
		.value("FEATURE", OpenMesh::Attributes::FEATURE)
		.value("TAGGED", OpenMesh::Attributes::TAGGED)
		.value("TAGGED2", OpenMesh::Attributes::TAGGED2)
		.value("FIXEDNONMANIFOLD", OpenMesh::Attributes::FIXEDNONMANIFOLD)
		.value("UNUSED", OpenMesh::Attributes::UNUSED)
		;

	class_<StatusInfo>("StatusInfo")
		.def("deleted", &StatusInfo::deleted)
		.def("set_deleted", &StatusInfo::set_deleted)
		.def("locked", &StatusInfo::locked)
		.def("set_locked", &StatusInfo::set_locked)
		.def("selected", &StatusInfo::selected)
		.def("set_selected", &StatusInfo::set_selected)
		.def("hidden", &StatusInfo::hidden)
		.def("set_hidden", &StatusInfo::set_hidden)
		.def("feature", &StatusInfo::feature)
		.def("set_feature", &StatusInfo::set_feature)
		.def("tagged", &StatusInfo::tagged)
		.def("set_tagged", &StatusInfo::set_tagged)
		.def("tagged2", &StatusInfo::tagged2)
		.def("set_tagged2", &StatusInfo::set_tagged2)
		.def("fixed_nonmanifold", &StatusInfo::fixed_nonmanifold)
		.def("set_fixed_nonmanifold", &StatusInfo::set_fixed_nonmanifold)
		.def("bits", &StatusInfo::bits)
		.def("set_bits", &StatusInfo::set_bits)
		.def("is_bit_set", &StatusInfo::is_bit_set)
		.def("set_bit", &StatusInfo::set_bit)
		.def("unset_bit", &StatusInfo::unset_bit)
		.def("change_bit", &StatusInfo::change_bit)
		;
}

BOOST_PYTHON_MODULE(openmesh) {
	expose_items();
	expose_handles();
	expose_status_bits_and_info();

	expose_vec<float,  2>("Vec2f");
	expose_vec<float,  3>("Vec3f");
	expose_vec<float,  4>("Vec4f");
	expose_vec<double, 2>("Vec2d");
	expose_vec<double, 3>("Vec3d");
	expose_vec<double, 4>("Vec4d");

	expose_mesh<PolyMesh>("PolyMesh");
	expose_mesh<TriMesh>("TriMesh");

	expose_iterator<OpenMesh::PolyConnectivity::VertexIter, &OpenMesh::ArrayKernel::n_vertices>("VertexIter");
	expose_iterator<OpenMesh::PolyConnectivity::HalfedgeIter, &OpenMesh::ArrayKernel::n_halfedges>("HalfedgeIter");
	expose_iterator<OpenMesh::PolyConnectivity::EdgeIter, &OpenMesh::ArrayKernel::n_edges>("EdgeIter");
	expose_iterator<OpenMesh::PolyConnectivity::FaceIter, &OpenMesh::ArrayKernel::n_faces>("FaceIter");

	expose_circulator<OpenMesh::PolyConnectivity::VertexVertexIter, VertexHandle>("VertexVertexIter");
	expose_circulator<OpenMesh::PolyConnectivity::VertexIHalfedgeIter, VertexHandle>("VertexIHalfedgeIter");
	expose_circulator<OpenMesh::PolyConnectivity::VertexOHalfedgeIter, VertexHandle>("VertexOHalfedgeIter");
	expose_circulator<OpenMesh::PolyConnectivity::VertexEdgeIter, VertexHandle>("VertexEdgeIter");
	expose_circulator<OpenMesh::PolyConnectivity::VertexFaceIter, VertexHandle>("VertexFaceIter");

	expose_circulator<OpenMesh::PolyConnectivity::FaceVertexIter, FaceHandle>("FaceVertexIter");
	expose_circulator<OpenMesh::PolyConnectivity::FaceHalfedgeIter, FaceHandle>("FaceHalfedgeIter");
	expose_circulator<OpenMesh::PolyConnectivity::FaceEdgeIter, FaceHandle>("FaceEdgeIter");
	expose_circulator<OpenMesh::PolyConnectivity::FaceFaceIter, FaceHandle>("FaceFaceIter");

	expose_circulator<OpenMesh::PolyConnectivity::HalfedgeLoopIter, HalfedgeHandle>("HalfedgeLoopIter");

	typedef IteratorWrapperT<PolyConnectivity::VertexIter, &ArrayKernel::n_vertices> VertexIterWrapper;
	typedef IteratorWrapperT<PolyConnectivity::HalfedgeIter, &ArrayKernel::n_halfedges> HalfedgeIterWrapper;
	typedef IteratorWrapperT<PolyConnectivity::EdgeIter, &ArrayKernel::n_edges> EdgeIterWrapper;
	typedef IteratorWrapperT<PolyConnectivity::FaceIter, &ArrayKernel::n_faces> FaceIterWrapper;

	expose_property_manager<VPropHandleT<object>, VertexHandle, VertexIterWrapper>("VPropertyManager");
	expose_property_manager<HPropHandleT<object>, HalfedgeHandle, HalfedgeIterWrapper>("HPropertyManager");
	expose_property_manager<EPropHandleT<object>, EdgeHandle, EdgeIterWrapper>("EPropertyManager");
	expose_property_manager<FPropHandleT<object>, FaceHandle, FaceIterWrapper>("FPropertyManager");

	expose_io();
}

} // namespace Python
} // namespace OpenMesh
