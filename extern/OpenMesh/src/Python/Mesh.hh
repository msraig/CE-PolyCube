#ifndef OPENMESH_PYTHON_MESH_HH
#define OPENMESH_PYTHON_MESH_HH

#include "Python/Bindings.hh"
#include "Python/Iterator.hh"
#include "Python/Circulator.hh"

#include <boost/python/stl_iterator.hpp>

namespace OpenMesh {
namespace Python {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(garbage_collection_overloads, garbage_collection, 0, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(add_property_overloads, add_property, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(copy_all_properties_overloads, copy_all_properties, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(delete_vertex_overloads, delete_vertex, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(delete_edge_overloads, delete_edge, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(delete_face_overloads, delete_face, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(is_boundary_overloads, is_boundary, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(find_feature_edges_overloads, find_feature_edges, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(update_normal_overloads, update_normal, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(update_halfedge_normals_overloads, update_halfedge_normals, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(calc_halfedge_normal_overloads, calc_halfedge_normal, 1, 2)

/**
 * Set the status of an item.
 *
 * @tparam Mesh A mesh type.
 * @tparam PropHandle A handle type.
 *
 * @param _self The mesh instance that is to be used.
 * @param _h The handle of the item whose status is to be set.
 * @param _info The status to be set.
 *
 * Depending on @ref OPENMESH_PYTHON_DEFAULT_POLICY, Mesh::status may
 * return by value instead of reference. This function ensures that the
 * status of an item can be changed nonetheless.
 */
template <class Mesh, class IndexHandle>
void set_status(Mesh& _self, IndexHandle _h, const OpenMesh::Attributes::StatusInfo& _info) {
	_self.status(_h) = _info;
}

/**
 * Set the value of a property of an item.
 *
 * @tparam Mesh A mesh type.
 * @tparam PropHandle A property handle type.
 * @tparam IndexHandle The appropriate handle type.
 *
 * @param _self The mesh instance that is to be used.
 * @param _ph The property that is to be set.
 * @param _h The handle of the item whose property is to be set.
 * @param _value The value to be set.
 *
 * Depending on @ref OPENMESH_PYTHON_DEFAULT_POLICY, Mesh::property may
 * return by value instead of reference. This function ensures that the
 * property value of an item can be changed nonetheless.
 */
template <class Mesh, class PropHandle, class IndexHandle>
void set_property(Mesh& _self, PropHandle _ph, IndexHandle _h, const object& _value) {
	_self.property(_ph, _h) = _value;
}

/**
 * Set the value of a mesh property.
 *
 * @tparam Mesh A mesh type.
 * @tparam PropHandle A property handle type.
 *
 * @param _self The mesh instance that is to be used.
 * @param _ph The property that is to be set.
 * @param _value The value to be set.
 *
 * Depending on @ref OPENMESH_PYTHON_DEFAULT_POLICY, Mesh::property may
 * return by value instead of reference. This function ensures that the
 * property value of an item can be changed nonetheless.
 */
template <class Mesh, class PropHandle>
void set_property(Mesh& _self, PropHandle _ph, const object& _value) {
	_self.property(_ph) = _value;
}

/**
 * Thin wrapper for assign_connectivity.
 *
 * @tparam Mesh A mesh type.
 * @tparam OtherMesh A mesh type.
 *
 * @param _self The mesh instance that is to be used.
 * @param _other The mesh from which the connectivity is to be copied.
 */
template <class Mesh, class OtherMesh>
void assign_connectivity(Mesh& _self, const OtherMesh& _other) {
	_self.assign_connectivity(_other);
}

/**
 * Get an iterator.
 */
template <class Mesh, class Iterator, size_t (ArrayKernel::*n_items)() const>
IteratorWrapperT<Iterator, n_items> get_iterator(Mesh& _self) {
	return IteratorWrapperT<Iterator, n_items>(_self, typename Iterator::value_type(0));
}

/**
 * Get a skipping iterator.
 */
template <class Mesh, class Iterator, size_t (ArrayKernel::*n_items)() const>
IteratorWrapperT<Iterator, n_items> get_skipping_iterator(Mesh& _self) {
	return IteratorWrapperT<Iterator, n_items>(_self, typename Iterator::value_type(0), true);
}

/**
 * Get a circulator.
 *
 * @tparam Mesh A Mesh type.
 * @tparam Circulator A circulator type.
 * @tparam CenterEntityHandle The appropriate handle type.
 *
 * @param _self The mesh instance that is to be used.
 * @param _handle The handle of the item to circulate around.
 */
template <class Mesh, class Circulator, class CenterEntityHandle>
CirculatorWrapperT<Circulator, CenterEntityHandle> get_circulator(Mesh& _self, CenterEntityHandle _handle) {
	return CirculatorWrapperT<Circulator, CenterEntityHandle>(_self, _handle);
}

/**
 * Garbage collection using lists instead of vectors to keep track of a set of
 * handles.
 *
 * @tparam Mesh A Mesh type.
 *
 * @param _self The mesh instance that is to be used.
 * @param _vh_to_update The list of vertex handles to be updated.
 * @param _hh_to_update The list of halfedge handles to be updated.
 * @param _fh_to_update The list of face handles to be updated.
 * @param _v Remove deleted vertices?
 * @param _e Remove deleted edges?
 * @param _f Remove deleted faces?
 */
template <class Mesh>
void garbage_collection(Mesh& _self, list& _vh_to_update, list& _hh_to_update, list& _fh_to_update, bool _v = true, bool _e = true, bool _f = true) {
	// Convert list of handles to vector of pointers
	stl_input_iterator<VertexHandle*> vh_begin(_vh_to_update);
	stl_input_iterator<VertexHandle*> vh_end;
	std::vector<VertexHandle*> vh_vector;
	vh_vector.insert(vh_vector.end(), vh_begin, vh_end);

	// Convert list of handles to vector of pointers
	stl_input_iterator<HalfedgeHandle*> hh_begin(_hh_to_update);
	stl_input_iterator<HalfedgeHandle*> hh_end;
	std::vector<HalfedgeHandle*> hh_vector;
	hh_vector.insert(hh_vector.end(), hh_begin, hh_end);

	// Convert list of handles to vector of pointers
	stl_input_iterator<FaceHandle*> fh_begin(_fh_to_update);
	stl_input_iterator<FaceHandle*> fh_end;
	std::vector<FaceHandle*> fh_vector;
	fh_vector.insert(fh_vector.end(), fh_begin, fh_end);

	// Call garbage collection
	_self.garbage_collection(vh_vector, hh_vector, fh_vector, _v, _e, _f);
}

/**
 * Add a new face from a %Python list of vertex handles.
 *
 * @tparam Mesh A Mesh type.
 *
 * @param _self The mesh instance that is to be used.
 * @param _vhandles The list of vertex handles.
 */
template<class Mesh>
FaceHandle add_face(Mesh& _self, const list& _vhandles) {
	stl_input_iterator<VertexHandle> begin(_vhandles);
	stl_input_iterator<VertexHandle> end;

	std::vector<VertexHandle> vector;
	vector.insert(vector.end(), begin, end);

	return _self.add_face(vector);
}

/**
 * This function template is used to expose mesh member functions that are only
 * available for a specific type of mesh (i.e. they are available for polygon
 * meshes or triangle meshes, but not both).
 *
 * @tparam Class A boost::python::class type.
 *
 * @param _class The boost::python::class instance for which the member
 * functions are to be defined.
 */
template <class Class>
void expose_type_specific_functions(Class& _class) {
	// See the template specializations below
}

/**
 * Function template specialization for polygon meshes.
 */
template <>
void expose_type_specific_functions(class_<PolyMesh>& _class) {
	typedef PolyMesh::Scalar Scalar;
	typedef PolyMesh::Point  Point;
	typedef PolyMesh::Normal Normal;
	typedef PolyMesh::Color  Color;

	FaceHandle (PolyMesh::*add_face_3_vh)(VertexHandle, VertexHandle, VertexHandle              ) = &PolyMesh::add_face;
	FaceHandle (PolyMesh::*add_face_4_vh)(VertexHandle, VertexHandle, VertexHandle, VertexHandle) = &PolyMesh::add_face;
	FaceHandle (*add_face_list)(PolyMesh&, const list&) = &add_face;

	void (PolyMesh::*split_eh_pt)(EdgeHandle, const Point&) = &PolyMesh::split;
	void (PolyMesh::*split_eh_vh)(EdgeHandle, VertexHandle) = &PolyMesh::split;
	void (PolyMesh::*split_fh_pt)(FaceHandle, const Point&) = &PolyMesh::split;
	void (PolyMesh::*split_fh_vh)(FaceHandle, VertexHandle) = &PolyMesh::split;

	Normal (PolyMesh::*calc_face_normal_pt)(const Point&, const Point&, const Point&) const = &PolyMesh::calc_face_normal;

	_class
		.def("add_face", add_face_3_vh)
		.def("add_face", add_face_4_vh)
		.def("add_face", add_face_list)

		.def("split", split_eh_pt)
		.def("split", split_eh_vh)
		.def("split", split_fh_pt)
		.def("split", split_fh_vh)

		.def("split_copy", &PolyMesh::split_copy)

		.def("calc_face_normal", calc_face_normal_pt)

		.def("insert_edge", &PolyMesh::insert_edge)
		;
}

/**
 * Function template specialization for triangle meshes.
 */
template <>
void expose_type_specific_functions(class_<TriMesh>& _class) {
	typedef TriMesh::Scalar Scalar;
	typedef TriMesh::Point  Point;
	typedef TriMesh::Normal Normal;
	typedef TriMesh::Color  Color;

	FaceHandle (TriMesh::*add_face_3_vh)(VertexHandle, VertexHandle, VertexHandle) = &TriMesh::add_face;
	FaceHandle (*add_face_list)(TriMesh&, const list&) = &add_face;

	VertexHandle (TriMesh::*split_eh_pt)(EdgeHandle, const Point&) = &TriMesh::split;
	void         (TriMesh::*split_eh_vh)(EdgeHandle, VertexHandle) = &TriMesh::split;
	VertexHandle (TriMesh::*split_fh_pt)(FaceHandle, const Point&) = &TriMesh::split;
	void         (TriMesh::*split_fh_vh)(FaceHandle, VertexHandle) = &TriMesh::split;

	VertexHandle (TriMesh::*split_copy_eh_pt)(EdgeHandle, const Point&) = &TriMesh::split_copy;
	void         (TriMesh::*split_copy_eh_vh)(EdgeHandle, VertexHandle) = &TriMesh::split_copy;
	VertexHandle (TriMesh::*split_copy_fh_pt)(FaceHandle, const Point&) = &TriMesh::split_copy;
	void         (TriMesh::*split_copy_fh_vh)(FaceHandle, VertexHandle) = &TriMesh::split_copy;

	HalfedgeHandle (TriMesh::*vertex_split_pt)(Point,        VertexHandle, VertexHandle, VertexHandle) = &TriMesh::vertex_split;
	HalfedgeHandle (TriMesh::*vertex_split_vh)(VertexHandle, VertexHandle, VertexHandle, VertexHandle) = &TriMesh::vertex_split;

	_class
		.def("add_face", add_face_3_vh)
		.def("add_face", add_face_list)

		.def("split", split_eh_pt)
		.def("split", split_eh_vh)
		.def("split", split_fh_pt)
		.def("split", split_fh_vh)

		.def("split_copy", split_copy_eh_pt)
		.def("split_copy", split_copy_eh_vh)
		.def("split_copy", split_copy_fh_pt)
		.def("split_copy", split_copy_fh_vh)

		.def("opposite_vh", &TriMesh::opposite_vh)
		.def("opposite_he_opposite_vh", &TriMesh::opposite_he_opposite_vh)

		.def("vertex_split", vertex_split_pt)
		.def("vertex_split", vertex_split_vh)

		.def("is_flip_ok", &TriMesh::is_flip_ok)
		.def("flip", &TriMesh::flip)
		;
}


/**
 * Expose a mesh type to %Python.
 *
 * @tparam Mesh A mesh type.
 *
 * @param _name The name of the mesh type to be exposed.
 */
template <class Mesh>
void expose_mesh(const char *_name) {
	using OpenMesh::Attributes::StatusInfo;

	typedef typename Mesh::Scalar Scalar;
	typedef typename Mesh::Point  Point;
	typedef typename Mesh::Normal Normal;
	typedef typename Mesh::Color  Color;

	//======================================================================
	//  KernelT Function Pointers
	//======================================================================

	// Get the i'th item
	VertexHandle   (Mesh::*vertex_handle_uint  )(unsigned int) const = &Mesh::vertex_handle;
	HalfedgeHandle (Mesh::*halfedge_handle_uint)(unsigned int) const = &Mesh::halfedge_handle;
	EdgeHandle     (Mesh::*edge_handle_uint    )(unsigned int) const = &Mesh::edge_handle;
	FaceHandle     (Mesh::*face_handle_uint    )(unsigned int) const = &Mesh::face_handle;

	// Delete items
	void (Mesh::*garbage_collection_bools)(bool, bool, bool) = &Mesh::garbage_collection;
	void (*garbage_collection_lists_bools)(Mesh&, list&, list&, list&, bool, bool, bool) = &garbage_collection;

	// Vertex connectivity
	HalfedgeHandle (Mesh::*halfedge_handle_vh)(VertexHandle) const = &Mesh::halfedge_handle;
	HalfedgeHandle (Mesh::*halfedge_handle_fh)(FaceHandle  ) const = &Mesh::halfedge_handle;

	// Halfedge connectivity
	FaceHandle     (Mesh::*face_handle_hh         )(HalfedgeHandle) const = &Mesh::face_handle;
	HalfedgeHandle (Mesh::*prev_halfedge_handle_hh)(HalfedgeHandle) const = &Mesh::prev_halfedge_handle;
	EdgeHandle     (Mesh::*edge_handle_hh         )(HalfedgeHandle) const = &Mesh::edge_handle;

	// Edge connectivity
	HalfedgeHandle (Mesh::*halfedge_handle_eh_uint)(EdgeHandle, unsigned int) const = &Mesh::halfedge_handle;

	// Set halfedge
	void (Mesh::*set_halfedge_handle_vh_hh)(VertexHandle, HalfedgeHandle) = &Mesh::set_halfedge_handle;
	void (Mesh::*set_halfedge_handle_fh_hh)(FaceHandle, HalfedgeHandle  ) = &Mesh::set_halfedge_handle;

	// Handle -> Item
	const typename Mesh::Vertex&   (Mesh::*vertex  )(VertexHandle  ) const = &Mesh::vertex;
	const typename Mesh::Halfedge& (Mesh::*halfedge)(HalfedgeHandle) const = &Mesh::halfedge;
	const typename Mesh::Edge&     (Mesh::*edge    )(EdgeHandle    ) const = &Mesh::edge;
	const typename Mesh::Face&     (Mesh::*face    )(FaceHandle    ) const = &Mesh::face;

	// Item -> Handle
	VertexHandle   (Mesh::*handle_v)(const typename Mesh::Vertex&  ) const = &Mesh::handle;
	HalfedgeHandle (Mesh::*handle_h)(const typename Mesh::Halfedge&) const = &Mesh::handle;
	EdgeHandle     (Mesh::*handle_e)(const typename Mesh::Edge&    ) const = &Mesh::handle;
	FaceHandle     (Mesh::*handle_f)(const typename Mesh::Face&    ) const = &Mesh::handle;

	// Get value of a standard property (point, normal, color)
	const typename Mesh::Point&  (Mesh::*point_vh )(VertexHandle  ) const = &Mesh::point;
	const typename Mesh::Normal& (Mesh::*normal_vh)(VertexHandle  ) const = &Mesh::normal;
	const typename Mesh::Normal& (Mesh::*normal_hh)(HalfedgeHandle) const = &Mesh::normal;
	const typename Mesh::Normal& (Mesh::*normal_fh)(FaceHandle    ) const = &Mesh::normal;
	const typename Mesh::Color&  (Mesh::*color_vh )(VertexHandle  ) const = &Mesh::color;
	const typename Mesh::Color&  (Mesh::*color_hh )(HalfedgeHandle) const = &Mesh::color;
	const typename Mesh::Color&  (Mesh::*color_eh )(EdgeHandle    ) const = &Mesh::color;
	const typename Mesh::Color&  (Mesh::*color_fh )(FaceHandle    ) const = &Mesh::color;

	// Get value of a standard property (texture coordinate)
	const typename Mesh::TexCoord1D& (Mesh::*texcoord1D_vh)(VertexHandle  ) const = &Mesh::texcoord1D;
	const typename Mesh::TexCoord1D& (Mesh::*texcoord1D_hh)(HalfedgeHandle) const = &Mesh::texcoord1D;
	const typename Mesh::TexCoord2D& (Mesh::*texcoord2D_vh)(VertexHandle  ) const = &Mesh::texcoord2D;
	const typename Mesh::TexCoord2D& (Mesh::*texcoord2D_hh)(HalfedgeHandle) const = &Mesh::texcoord2D;
	const typename Mesh::TexCoord3D& (Mesh::*texcoord3D_vh)(VertexHandle  ) const = &Mesh::texcoord3D;
	const typename Mesh::TexCoord3D& (Mesh::*texcoord3D_hh)(HalfedgeHandle) const = &Mesh::texcoord3D;

	// Get value of a standard property (status)
	const StatusInfo& (Mesh::*status_vh)(VertexHandle  ) const = &Mesh::status;
	const StatusInfo& (Mesh::*status_hh)(HalfedgeHandle) const = &Mesh::status;
	const StatusInfo& (Mesh::*status_eh)(EdgeHandle    ) const = &Mesh::status;
	const StatusInfo& (Mesh::*status_fh)(FaceHandle    ) const = &Mesh::status;

	// Set value of a standard property (point, normal, color)
	void (Mesh::*set_normal_vh)(VertexHandle,   const typename Mesh::Normal&) = &Mesh::set_normal;
	void (Mesh::*set_normal_hh)(HalfedgeHandle, const typename Mesh::Normal&) = &Mesh::set_normal;
	void (Mesh::*set_normal_fh)(FaceHandle,     const typename Mesh::Normal&) = &Mesh::set_normal;
	void (Mesh::*set_color_vh )(VertexHandle,   const typename Mesh::Color& ) = &Mesh::set_color;
	void (Mesh::*set_color_hh )(HalfedgeHandle, const typename Mesh::Color& ) = &Mesh::set_color;
	void (Mesh::*set_color_eh )(EdgeHandle,     const typename Mesh::Color& ) = &Mesh::set_color;
	void (Mesh::*set_color_fh )(FaceHandle,     const typename Mesh::Color& ) = &Mesh::set_color;

	// Set value of a standard property (texture coordinate)
	void (Mesh::*set_texcoord1D_vh)(VertexHandle,   const typename Mesh::TexCoord1D&) = &Mesh::set_texcoord1D;
	void (Mesh::*set_texcoord1D_hh)(HalfedgeHandle, const typename Mesh::TexCoord1D&) = &Mesh::set_texcoord1D;
	void (Mesh::*set_texcoord2D_vh)(VertexHandle,   const typename Mesh::TexCoord2D&) = &Mesh::set_texcoord2D;
	void (Mesh::*set_texcoord2D_hh)(HalfedgeHandle, const typename Mesh::TexCoord2D&) = &Mesh::set_texcoord2D;
	void (Mesh::*set_texcoord3D_vh)(VertexHandle,   const typename Mesh::TexCoord3D&) = &Mesh::set_texcoord3D;
	void (Mesh::*set_texcoord3D_hh)(HalfedgeHandle, const typename Mesh::TexCoord3D&) = &Mesh::set_texcoord3D;

	// Set value of a standard property (status)
	void (*set_status_vh)(Mesh&, VertexHandle,   const StatusInfo&) = &set_status;
	void (*set_status_hh)(Mesh&, HalfedgeHandle, const StatusInfo&) = &set_status;
	void (*set_status_eh)(Mesh&, EdgeHandle,     const StatusInfo&) = &set_status;
	void (*set_status_fh)(Mesh&, FaceHandle,     const StatusInfo&) = &set_status;

	// Property management - add property
	void (Mesh::*add_property_vph)(VPropHandleT<object>&, const std::string&) = &Mesh::add_property;
	void (Mesh::*add_property_eph)(EPropHandleT<object>&, const std::string&) = &Mesh::add_property;
	void (Mesh::*add_property_hph)(HPropHandleT<object>&, const std::string&) = &Mesh::add_property;
	void (Mesh::*add_property_fph)(FPropHandleT<object>&, const std::string&) = &Mesh::add_property;
	void (Mesh::*add_property_mph)(MPropHandleT<object>&, const std::string&) = &Mesh::add_property;

	// Property management - remove property
	void (Mesh::*remove_property_vph)(VPropHandleT<object>&) = &Mesh::remove_property;
	void (Mesh::*remove_property_eph)(EPropHandleT<object>&) = &Mesh::remove_property;
	void (Mesh::*remove_property_hph)(HPropHandleT<object>&) = &Mesh::remove_property;
	void (Mesh::*remove_property_fph)(FPropHandleT<object>&) = &Mesh::remove_property;
	void (Mesh::*remove_property_mph)(MPropHandleT<object>&) = &Mesh::remove_property;

	// Property management - get property by name
	bool (Mesh::*get_property_handle_vph)(VPropHandleT<object>&, const std::string&) const = &Mesh::get_property_handle;
	bool (Mesh::*get_property_handle_eph)(EPropHandleT<object>&, const std::string&) const = &Mesh::get_property_handle;
	bool (Mesh::*get_property_handle_hph)(HPropHandleT<object>&, const std::string&) const = &Mesh::get_property_handle;
	bool (Mesh::*get_property_handle_fph)(FPropHandleT<object>&, const std::string&) const = &Mesh::get_property_handle;
	bool (Mesh::*get_property_handle_mph)(MPropHandleT<object>&, const std::string&) const = &Mesh::get_property_handle;

	// Property management - get property value for an item
	const object& (Mesh::*property_vertex  )(VPropHandleT<object>, VertexHandle  ) const = &Mesh::property;
	const object& (Mesh::*property_edge    )(EPropHandleT<object>, EdgeHandle    ) const = &Mesh::property;
	const object& (Mesh::*property_halfedge)(HPropHandleT<object>, HalfedgeHandle) const = &Mesh::property;
	const object& (Mesh::*property_face    )(FPropHandleT<object>, FaceHandle    ) const = &Mesh::property;
	const object& (Mesh::*property_mesh    )(MPropHandleT<object>                ) const = &Mesh::property;

	// Property management - set property value for an item
	void (*set_property_vertex  )(Mesh&, VPropHandleT<object>, VertexHandle,   const object&) = &set_property;
	void (*set_property_edge    )(Mesh&, EPropHandleT<object>, EdgeHandle,     const object&) = &set_property;
	void (*set_property_halfedge)(Mesh&, HPropHandleT<object>, HalfedgeHandle, const object&) = &set_property;
	void (*set_property_face    )(Mesh&, FPropHandleT<object>, FaceHandle,     const object&) = &set_property;
	void (*set_property_mesh    )(Mesh&, MPropHandleT<object>,                 const object&) = &set_property;

	// Low-level adding new items
	VertexHandle (Mesh::*new_vertex_void )(void                        ) = &Mesh::new_vertex;
	VertexHandle (Mesh::*new_vertex_point)(const typename Mesh::Point& ) = &Mesh::new_vertex;
	FaceHandle   (Mesh::*new_face_void   )(void                        ) = &Mesh::new_face;
	FaceHandle   (Mesh::*new_face_face   )(const typename Mesh::Face&  ) = &Mesh::new_face;

	// Kernel item iterators
	IteratorWrapperT<typename Mesh::VertexIter,   &Mesh::n_vertices > (*vertices )(Mesh&) = &get_iterator;
	IteratorWrapperT<typename Mesh::HalfedgeIter, &Mesh::n_halfedges> (*halfedges)(Mesh&) = &get_iterator;
	IteratorWrapperT<typename Mesh::EdgeIter,     &Mesh::n_edges    > (*edges    )(Mesh&) = &get_iterator;
	IteratorWrapperT<typename Mesh::FaceIter,     &Mesh::n_faces    > (*faces    )(Mesh&) = &get_iterator;

	IteratorWrapperT<typename Mesh::VertexIter,   &Mesh::n_vertices > (*svertices )(Mesh&) = &get_skipping_iterator;
	IteratorWrapperT<typename Mesh::HalfedgeIter, &Mesh::n_halfedges> (*shalfedges)(Mesh&) = &get_skipping_iterator;
	IteratorWrapperT<typename Mesh::EdgeIter,     &Mesh::n_edges    > (*sedges    )(Mesh&) = &get_skipping_iterator;
	IteratorWrapperT<typename Mesh::FaceIter,     &Mesh::n_faces    > (*sfaces    )(Mesh&) = &get_skipping_iterator;

	//======================================================================
	//  BaseKernel Function Pointers
	//======================================================================

	// Copy property
	void (Mesh::*copy_property_vprop)(VPropHandleT<object>&, VertexHandle,   VertexHandle  ) = &Mesh::copy_property;
	void (Mesh::*copy_property_hprop)(HPropHandleT<object>,  HalfedgeHandle, HalfedgeHandle) = &Mesh::copy_property;
	void (Mesh::*copy_property_eprop)(EPropHandleT<object>,  EdgeHandle,     EdgeHandle    ) = &Mesh::copy_property;
	void (Mesh::*copy_property_fprop)(FPropHandleT<object>,  FaceHandle,     FaceHandle    ) = &Mesh::copy_property;

	// Copy all properties
	void (Mesh::*copy_all_properties_vh_vh_bool)(VertexHandle,   VertexHandle,   bool) = &Mesh::copy_all_properties;
	void (Mesh::*copy_all_properties_hh_hh_bool)(HalfedgeHandle, HalfedgeHandle, bool) = &Mesh::copy_all_properties;
	void (Mesh::*copy_all_properties_eh_eh_bool)(EdgeHandle,     EdgeHandle,     bool) = &Mesh::copy_all_properties;
	void (Mesh::*copy_all_properties_fh_fh_bool)(FaceHandle,     FaceHandle,     bool) = &Mesh::copy_all_properties;

	//======================================================================
	//  PolyConnectivity Function Pointers
	//======================================================================

	// Assign connectivity
	void (*assign_connectivity_poly)(Mesh&, const PolyMesh&) = &assign_connectivity;
	void (*assign_connectivity_tri )(Mesh&, const TriMesh& ) = &assign_connectivity;

	// Vertex and face valence
	unsigned int (Mesh::*valence_vh)(VertexHandle) const = &Mesh::valence;
	unsigned int (Mesh::*valence_fh)(FaceHandle  ) const = &Mesh::valence;

	// Triangulate face or mesh
	void (Mesh::*triangulate_fh  )(FaceHandle) = &Mesh::triangulate;
	void (Mesh::*triangulate_void)(          ) = &Mesh::triangulate;

	// Deleting mesh items and other connectivity/topology modifications
	void (Mesh::*delete_vertex)(VertexHandle, bool) = &Mesh::delete_vertex;
	void (Mesh::*delete_edge  )(EdgeHandle,   bool) = &Mesh::delete_edge;
	void (Mesh::*delete_face  )(FaceHandle,   bool) = &Mesh::delete_face;

	// Vertex and Face circulators
	CirculatorWrapperT<typename Mesh::VertexVertexIter,    VertexHandle  > (*vv )(Mesh&, VertexHandle  ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::VertexIHalfedgeIter, VertexHandle  > (*vih)(Mesh&, VertexHandle  ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::VertexOHalfedgeIter, VertexHandle  > (*voh)(Mesh&, VertexHandle  ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::VertexEdgeIter,      VertexHandle  > (*ve )(Mesh&, VertexHandle  ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::VertexFaceIter,      VertexHandle  > (*vf )(Mesh&, VertexHandle  ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::FaceVertexIter,      FaceHandle    > (*fv )(Mesh&, FaceHandle    ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::FaceHalfedgeIter,    FaceHandle    > (*fh )(Mesh&, FaceHandle    ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::FaceEdgeIter,        FaceHandle    > (*fe )(Mesh&, FaceHandle    ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::FaceFaceIter,        FaceHandle    > (*ff )(Mesh&, FaceHandle    ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::HalfedgeLoopIter,    HalfedgeHandle> (*hl )(Mesh&, HalfedgeHandle) = &get_circulator;

	// Boundary and manifold tests
	bool (Mesh::*is_boundary_hh)(HalfedgeHandle  ) const = &Mesh::is_boundary;
	bool (Mesh::*is_boundary_eh)(EdgeHandle      ) const = &Mesh::is_boundary;
	bool (Mesh::*is_boundary_vh)(VertexHandle    ) const = &Mesh::is_boundary;
	bool (Mesh::*is_boundary_fh)(FaceHandle, bool) const = &Mesh::is_boundary;

	// Generic handle derefertiation
	const typename Mesh::Vertex&   (Mesh::*deref_vh)(VertexHandle  ) const = &Mesh::deref;
	const typename Mesh::Halfedge& (Mesh::*deref_hh)(HalfedgeHandle) const = &Mesh::deref;
	const typename Mesh::Edge&     (Mesh::*deref_eh)(EdgeHandle    ) const = &Mesh::deref;
	const typename Mesh::Face&     (Mesh::*deref_fh)(FaceHandle    ) const = &Mesh::deref;

	//======================================================================
	//  PolyMeshT Function Pointers
	//======================================================================

	void (Mesh::*calc_edge_vector_eh_normal)(EdgeHandle,     Normal&) const = &Mesh::calc_edge_vector;
	void (Mesh::*calc_edge_vector_hh_normal)(HalfedgeHandle, Normal&) const = &Mesh::calc_edge_vector;

	Normal (Mesh::*calc_edge_vector_eh)(EdgeHandle    ) const = &Mesh::calc_edge_vector;
	Normal (Mesh::*calc_edge_vector_hh)(HalfedgeHandle) const = &Mesh::calc_edge_vector;

	Scalar (Mesh::*calc_edge_length_eh)(EdgeHandle    ) const = &Mesh::calc_edge_length;
	Scalar (Mesh::*calc_edge_length_hh)(HalfedgeHandle) const = &Mesh::calc_edge_length;

	Scalar (Mesh::*calc_edge_sqr_length_eh)(EdgeHandle    ) const = &Mesh::calc_edge_sqr_length;
	Scalar (Mesh::*calc_edge_sqr_length_hh)(HalfedgeHandle) const = &Mesh::calc_edge_sqr_length;

	Scalar (Mesh::*calc_dihedral_angle_fast_hh)(HalfedgeHandle) const = &Mesh::calc_dihedral_angle_fast;
	Scalar (Mesh::*calc_dihedral_angle_fast_eh)(EdgeHandle    ) const = &Mesh::calc_dihedral_angle_fast;

	Scalar (Mesh::*calc_dihedral_angle_hh)(HalfedgeHandle) const = &Mesh::calc_dihedral_angle;
	Scalar (Mesh::*calc_dihedral_angle_eh)(EdgeHandle    ) const = &Mesh::calc_dihedral_angle;

	unsigned int (Mesh::*find_feature_edges)(Scalar) = &Mesh::find_feature_edges;

	void (Mesh::*split_fh_vh)(FaceHandle, VertexHandle) = &Mesh::split;
	void (Mesh::*split_eh_vh)(EdgeHandle, VertexHandle) = &Mesh::split;

	void (Mesh::*update_normal_fh)(FaceHandle            ) = &Mesh::update_normal;
	void (Mesh::*update_normal_hh)(HalfedgeHandle, double) = &Mesh::update_normal;
	void (Mesh::*update_normal_vh)(VertexHandle          ) = &Mesh::update_normal;

	void (Mesh::*update_halfedge_normals)(double) = &Mesh::update_halfedge_normals;

	Normal (Mesh::*calc_face_normal    )(FaceHandle            ) const = &Mesh::calc_face_normal;
	Normal (Mesh::*calc_halfedge_normal)(HalfedgeHandle, double) const = &Mesh::calc_halfedge_normal;

	void  (Mesh::*calc_face_centroid_fh_point)(FaceHandle, Point&) const = &Mesh::calc_face_centroid;
	Point (Mesh::*calc_face_centroid_fh      )(FaceHandle        ) const = &Mesh::calc_face_centroid;

	//======================================================================
	//  Mesh Type
	//======================================================================

	class_<Mesh> class_mesh(_name);

	class_mesh

		//======================================================================
		//  KernelT
		//======================================================================

		.def("reserve", &Mesh::reserve)

		.def("vertex", vertex, return_value_policy<reference_existing_object>())
		.def("halfedge", halfedge, return_value_policy<reference_existing_object>())
		.def("edge", edge, return_value_policy<reference_existing_object>())
		.def("face", face, return_value_policy<reference_existing_object>())

		.def("handle", handle_v)
		.def("handle", handle_h)
		.def("handle", handle_e)
		.def("handle", handle_f)

		.def("vertex_handle", vertex_handle_uint)
		.def("halfedge_handle", halfedge_handle_uint)
		.def("edge_handle", edge_handle_uint)
		.def("face_handle", face_handle_uint)

		.def("clear", &Mesh::clear)
		.def("clean", &Mesh::clean)
		.def("garbage_collection", garbage_collection_bools, garbage_collection_overloads())
		.def("garbage_collection", garbage_collection_lists_bools)

		.def("n_vertices", &Mesh::n_vertices)
		.def("n_halfedges", &Mesh::n_halfedges)
		.def("n_edges", &Mesh::n_edges)
		.def("n_faces", &Mesh::n_faces)
		.def("vertices_empty", &Mesh::vertices_empty)
		.def("halfedges_empty", &Mesh::halfedges_empty)
		.def("edges_empty", &Mesh::edges_empty)
		.def("faces_empty", &Mesh::faces_empty)

		.def("halfedge_handle", halfedge_handle_vh)
		.def("set_halfedge_handle", set_halfedge_handle_vh_hh)

		.def("to_vertex_handle", &Mesh::to_vertex_handle)
		.def("from_vertex_handle", &Mesh::from_vertex_handle)
		.def("set_vertex_handle", &Mesh::set_vertex_handle)
		.def("face_handle", face_handle_hh)
		.def("set_face_handle", &Mesh::set_face_handle)
		.def("next_halfedge_handle", &Mesh::next_halfedge_handle)
		.def("set_next_halfedge_handle", &Mesh::set_next_halfedge_handle)
		.def("prev_halfedge_handle", prev_halfedge_handle_hh)
		.def("opposite_halfedge_handle", &Mesh::opposite_halfedge_handle)
		.def("ccw_rotated_halfedge_handle", &Mesh::ccw_rotated_halfedge_handle)
		.def("cw_rotated_halfedge_handle", &Mesh::cw_rotated_halfedge_handle)
		.def("edge_handle", edge_handle_hh)

		.def("halfedge_handle", halfedge_handle_eh_uint)

		.def("halfedge_handle", halfedge_handle_fh)
		.def("set_halfedge_handle", set_halfedge_handle_fh_hh)

		.def("point", point_vh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_point", &Mesh::set_point)
		.def("normal", normal_vh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_normal", set_normal_vh)
		.def("normal", normal_hh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_normal", set_normal_hh)
		.def("color", color_vh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_color", set_color_vh)
		.def("texcoord1D", texcoord1D_vh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_texcoord1D", set_texcoord1D_vh)
		.def("texcoord2D", texcoord2D_vh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_texcoord2D", set_texcoord2D_vh)
		.def("texcoord3D", texcoord3D_vh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_texcoord3D", set_texcoord3D_vh)
		.def("texcoord1D", texcoord1D_hh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_texcoord1D", set_texcoord1D_hh)
		.def("texcoord2D", texcoord2D_hh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_texcoord2D", set_texcoord2D_hh)
		.def("texcoord3D", texcoord3D_hh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_texcoord3D", set_texcoord3D_hh)
		.def("status", status_vh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_status", set_status_vh)
		.def("status", status_hh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_status", set_status_hh)
		.def("color", color_hh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_color", set_color_hh)
		.def("color", color_eh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_color", set_color_eh)
		.def("status", status_eh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_status", set_status_eh)
		.def("normal", normal_fh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_normal", set_normal_fh)
		.def("color", color_fh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_color", set_color_fh)
		.def("status", status_fh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_status", set_status_fh)

		.def("request_vertex_normals", &Mesh::request_vertex_normals)
		.def("request_vertex_colors", &Mesh::request_vertex_colors)
		.def("request_vertex_texcoords1D", &Mesh::request_vertex_texcoords1D)
		.def("request_vertex_texcoords2D", &Mesh::request_vertex_texcoords2D)
		.def("request_vertex_texcoords3D", &Mesh::request_vertex_texcoords3D)
		.def("request_vertex_status", &Mesh::request_vertex_status)
		.def("request_halfedge_status", &Mesh::request_halfedge_status)
		.def("request_halfedge_normals", &Mesh::request_halfedge_normals)
		.def("request_halfedge_colors", &Mesh::request_halfedge_colors)
		.def("request_halfedge_texcoords1D", &Mesh::request_halfedge_texcoords1D)
		.def("request_halfedge_texcoords2D", &Mesh::request_halfedge_texcoords2D)
		.def("request_halfedge_texcoords3D", &Mesh::request_halfedge_texcoords3D)
		.def("request_edge_status", &Mesh::request_edge_status)
		.def("request_edge_colors", &Mesh::request_edge_colors)
		.def("request_face_normals", &Mesh::request_face_normals)
		.def("request_face_colors", &Mesh::request_face_colors)
		.def("request_face_status", &Mesh::request_face_status)
		.def("request_face_texture_index", &Mesh::request_face_texture_index)

		.def("release_vertex_normals", &Mesh::release_vertex_normals)
		.def("release_vertex_colors", &Mesh::release_vertex_colors)
		.def("release_vertex_texcoords1D", &Mesh::release_vertex_texcoords1D)
		.def("release_vertex_texcoords2D", &Mesh::release_vertex_texcoords2D)
		.def("release_vertex_texcoords3D", &Mesh::release_vertex_texcoords3D)
		.def("release_vertex_status", &Mesh::release_vertex_status)
		.def("release_halfedge_status", &Mesh::release_halfedge_status)
		.def("release_halfedge_normals", &Mesh::release_halfedge_normals)
		.def("release_halfedge_colors", &Mesh::release_halfedge_colors)
		.def("release_halfedge_texcoords1D", &Mesh::release_halfedge_texcoords1D)
		.def("release_halfedge_texcoords2D", &Mesh::release_halfedge_texcoords2D)
		.def("release_halfedge_texcoords3D", &Mesh::release_halfedge_texcoords3D)
		.def("release_edge_status", &Mesh::release_edge_status)
		.def("release_edge_colors", &Mesh::release_edge_colors)
		.def("release_face_normals", &Mesh::release_face_normals)
		.def("release_face_colors", &Mesh::release_face_colors)
		.def("release_face_status", &Mesh::release_face_status)
		.def("release_face_texture_index", &Mesh::release_face_texture_index)

		.def("has_vertex_normals", &Mesh::has_vertex_normals)
		.def("has_vertex_colors", &Mesh::has_vertex_colors)
		.def("has_vertex_texcoords1D", &Mesh::has_vertex_texcoords1D)
		.def("has_vertex_texcoords2D", &Mesh::has_vertex_texcoords2D)
		.def("has_vertex_texcoords3D", &Mesh::has_vertex_texcoords3D)
		.def("has_vertex_status", &Mesh::has_vertex_status)
		.def("has_halfedge_status", &Mesh::has_halfedge_status)
		.def("has_halfedge_normals", &Mesh::has_halfedge_normals)
		.def("has_halfedge_colors", &Mesh::has_halfedge_colors)
		.def("has_halfedge_texcoords1D", &Mesh::has_halfedge_texcoords1D)
		.def("has_halfedge_texcoords2D", &Mesh::has_halfedge_texcoords2D)
		.def("has_halfedge_texcoords3D", &Mesh::has_halfedge_texcoords3D)
		.def("has_edge_status", &Mesh::has_edge_status)
		.def("has_edge_colors", &Mesh::has_edge_colors)
		.def("has_face_normals", &Mesh::has_face_normals)
		.def("has_face_colors", &Mesh::has_face_colors)
		.def("has_face_status", &Mesh::has_face_status)
		.def("has_face_texture_index", &Mesh::has_face_texture_index)

		.def("add_property", add_property_vph, add_property_overloads())
		.def("add_property", add_property_eph, add_property_overloads())
		.def("add_property", add_property_hph, add_property_overloads())
		.def("add_property", add_property_fph, add_property_overloads())
		.def("add_property", add_property_mph, add_property_overloads())

		.def("remove_property", remove_property_vph)
		.def("remove_property", remove_property_eph)
		.def("remove_property", remove_property_hph)
		.def("remove_property", remove_property_fph)
		.def("remove_property", remove_property_mph)

		.def("get_property_handle", get_property_handle_vph)
		.def("get_property_handle", get_property_handle_eph)
		.def("get_property_handle", get_property_handle_hph)
		.def("get_property_handle", get_property_handle_fph)
		.def("get_property_handle", get_property_handle_mph)

		.def("property", property_vertex, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("property", property_edge, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("property", property_halfedge, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("property", property_face, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("property", property_mesh, OPENMESH_PYTHON_DEFAULT_POLICY)

		.def("set_property", set_property_vertex)
		.def("set_property", set_property_edge)
		.def("set_property", set_property_halfedge)
		.def("set_property", set_property_face)
		.def("set_property", set_property_mesh)

		.def("new_vertex", new_vertex_void)
		.def("new_vertex", new_vertex_point)
		.def("new_edge", &Mesh::new_edge)
		.def("new_face", new_face_void)
		.def("new_face", new_face_face)

		.def("vertices", vertices)
		.def("halfedges", halfedges)
		.def("edges", edges)
		.def("faces", faces)

		.def("svertices", svertices)
		.def("shalfedges", shalfedges)
		.def("sedges", sedges)
		.def("sfaces", sfaces)

		//======================================================================
		//  BaseKernel
		//======================================================================

		.def("copy_property", copy_property_vprop)
		.def("copy_property", copy_property_hprop)
		.def("copy_property", copy_property_eprop)
		.def("copy_property", copy_property_fprop)

		.def("copy_all_properties", copy_all_properties_vh_vh_bool, copy_all_properties_overloads())
		.def("copy_all_properties", copy_all_properties_hh_hh_bool, copy_all_properties_overloads())
		.def("copy_all_properties", copy_all_properties_eh_eh_bool, copy_all_properties_overloads())
		.def("copy_all_properties", copy_all_properties_fh_fh_bool, copy_all_properties_overloads())

		//======================================================================
		//  PolyConnectivity
		//======================================================================

		.def("assign_connectivity", assign_connectivity_poly)
		.def("assign_connectivity", assign_connectivity_tri)

		.def("opposite_face_handle", &Mesh::opposite_face_handle)
		.def("adjust_outgoing_halfedge", &Mesh::adjust_outgoing_halfedge)
		.def("find_halfedge", &Mesh::find_halfedge)
		.def("valence", valence_vh)
		.def("valence", valence_fh)
		.def("collapse", &Mesh::collapse)
		.def("is_simple_link", &Mesh::is_simple_link)
		.def("is_simply_connected", &Mesh::is_simply_connected)
		.def("remove_edge", &Mesh::remove_edge)
		.def("reinsert_edge", &Mesh::reinsert_edge)
		.def("triangulate", triangulate_fh)
		.def("triangulate", triangulate_void)
		.def("split_edge", &Mesh::split_edge)
		.def("split_edge_copy", &Mesh::split_edge_copy)

		.def("add_vertex", &Mesh::add_vertex)

		.def("is_collapse_ok",  &Mesh::is_collapse_ok)
		.def("delete_vertex", delete_vertex, delete_vertex_overloads())
		.def("delete_edge", delete_edge, delete_edge_overloads())
		.def("delete_face", delete_face, delete_face_overloads())

		.def("vv", vv)
		.def("vih", vih)
		.def("voh", voh)
		.def("ve", ve)
		.def("vf", vf)

		.def("fv", fv)
		.def("fh", fh)
		.def("fe", fe)
		.def("ff", ff)

		.def("hl", hl)

		.def("is_boundary", is_boundary_hh)
		.def("is_boundary", is_boundary_eh)
		.def("is_boundary", is_boundary_vh)
		.def("is_boundary", is_boundary_fh, is_boundary_overloads())
		.def("is_manifold", &Mesh::is_manifold)

		.def("deref", deref_vh, return_value_policy<reference_existing_object>())
		.def("deref", deref_hh, return_value_policy<reference_existing_object>())
		.def("deref", deref_eh, return_value_policy<reference_existing_object>())
		.def("deref", deref_fh, return_value_policy<reference_existing_object>())

		.def("is_triangles", &Mesh::is_triangles)
		.staticmethod("is_triangles")

		.def_readonly("InvalidVertexHandle", &Mesh::InvalidVertexHandle)
		.def_readonly("InvalidHalfedgeHandle", &Mesh::InvalidHalfedgeHandle)
		.def_readonly("InvalidEdgeHandle", &Mesh::InvalidEdgeHandle)
		.def_readonly("InvalidFaceHandle", &Mesh::InvalidFaceHandle)

		//======================================================================
		//  PolyMeshT
		//======================================================================

		.def("add_vertex", &Mesh::add_vertex)

		.def("calc_edge_vector", calc_edge_vector_eh_normal)
		.def("calc_edge_vector", calc_edge_vector_eh)
		.def("calc_edge_vector", calc_edge_vector_hh_normal)
		.def("calc_edge_vector", calc_edge_vector_hh)

		.def("calc_edge_length", calc_edge_length_eh)
		.def("calc_edge_length", calc_edge_length_hh)
		.def("calc_edge_sqr_length", calc_edge_sqr_length_eh)
		.def("calc_edge_sqr_length", calc_edge_sqr_length_hh)

		.def("calc_sector_vectors", &Mesh::calc_sector_vectors)
		.def("calc_sector_angle", &Mesh::calc_sector_angle)
		.def("calc_sector_normal", &Mesh::calc_sector_normal)
		.def("calc_sector_area", &Mesh::calc_sector_area)

		.def("calc_dihedral_angle_fast", calc_dihedral_angle_fast_hh)
		.def("calc_dihedral_angle_fast", calc_dihedral_angle_fast_eh)
		.def("calc_dihedral_angle", calc_dihedral_angle_hh)
		.def("calc_dihedral_angle", calc_dihedral_angle_eh)

		.def("find_feature_edges", find_feature_edges, find_feature_edges_overloads())

		.def("split", split_fh_vh)
		.def("split", split_eh_vh)

		.def("update_normals", &Mesh::update_normals)
		.def("update_normal", update_normal_fh)
		.def("update_face_normals", &Mesh::update_face_normals)

		.def("calc_face_normal", calc_face_normal)

		.def("calc_face_centroid", calc_face_centroid_fh_point)
		.def("calc_face_centroid", calc_face_centroid_fh)

		.def("update_normal", update_normal_hh, update_normal_overloads())
		.def("update_halfedge_normals", update_halfedge_normals, update_halfedge_normals_overloads())

		.def("calc_halfedge_normal", calc_halfedge_normal, calc_halfedge_normal_overloads())

		.def("is_estimated_feature_edge", &Mesh::is_estimated_feature_edge)

		.def("update_normal", update_normal_vh)
		.def("update_vertex_normals", &Mesh::update_vertex_normals)

		.def("calc_vertex_normal", &Mesh::calc_vertex_normal)
		.def("calc_vertex_normal_fast", &Mesh::calc_vertex_normal_fast)
		.def("calc_vertex_normal_correct", &Mesh::calc_vertex_normal_correct)
		.def("calc_vertex_normal_loop", &Mesh::calc_vertex_normal_loop)

		.def("is_polymesh", &Mesh::is_polymesh)
		.staticmethod("is_polymesh")

		.def("is_trimesh", &Mesh::is_trimesh)
		.staticmethod("is_trimesh")
		;

	expose_type_specific_functions(class_mesh);

	//======================================================================
	//  Nested Types
	//======================================================================

	// Enter mesh scope
	scope scope_mesh = class_mesh;

	// Point
	const boost::python::type_info point_info = type_id<typename Mesh::Point>();
	const converter::registration * point_registration = converter::registry::query(point_info);
	scope_mesh.attr("Point") = handle<>(point_registration->m_class_object);

	// Normal
	const boost::python::type_info normal_info = type_id<typename Mesh::Normal>();
	const converter::registration * normal_registration = converter::registry::query(normal_info);
	scope_mesh.attr("Normal") = handle<>(normal_registration->m_class_object);

	// Color
	const boost::python::type_info color_info = type_id<typename Mesh::Color>();
	const converter::registration * color_registration = converter::registry::query(color_info);
	scope_mesh.attr("Color") = handle<>(color_registration->m_class_object);

	// TexCoord2D
	const boost::python::type_info texcoord2d_info = type_id<typename Mesh::TexCoord2D>();
	const converter::registration * texcoord2d_registration = converter::registry::query(texcoord2d_info);
	scope_mesh.attr("TexCoord2D") = handle<>(texcoord2d_registration->m_class_object);

	// TexCoord3D
	const boost::python::type_info texcoord3d_info = type_id<typename Mesh::TexCoord3D>();
	const converter::registration * texcoord3d_registration = converter::registry::query(texcoord3d_info);
	scope_mesh.attr("TexCoord3D") = handle<>(texcoord3d_registration->m_class_object);
}

} // namespace OpenMesh
} // namespace Python

#endif
