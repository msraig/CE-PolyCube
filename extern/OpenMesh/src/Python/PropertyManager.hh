#ifndef OPENMESH_PYTHON_PROPERTYMANAGER_HH
#define OPENMESH_PYTHON_PROPERTYMANAGER_HH

#include "Python/Bindings.hh"
#include "OpenMesh/Core/Utils/PropertyManager.hh"

namespace OpenMesh {
namespace Python {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(retain_overloads, retain, 0, 1)

/**
 * Implementation of %Python's \_\_getitem\_\_ magic method.
 *
 * @tparam PropertyManager A property manager type.
 * @tparam IndexHandle The appropriate handle type.
 *
 * @param _self The property manager instance that is to be used.
 * @param _handle The index of the property value to be returned.
 *
 * @return The requested property value.
 */
template <class PropertyManager, class IndexHandle>
object propman_get_item(PropertyManager& _self, IndexHandle _handle) {
	return _self[_handle];
}

/**
 * Implementation of %Python's \_\_setitem\_\_ magic method.
 *
 * @tparam PropertyManager A property manager type.
 * @tparam IndexHandle The appropriate handle type.
 *
 * @param _self The property manager instance that is to be used.
 * @param _handle The index of the property value to be set.
 * @param _value The property value to be set.
 */
template <class PropertyManager, class IndexHandle>
void propman_set_item(PropertyManager& _self, IndexHandle _handle, object _value) {
	_self[_handle] = _value;
}

/**
 * Conveniently set the property value for an entire range of mesh items
 * using a %Python iterator.
 *
 * @tparam PropertyManager A property manager type.
 * @tparam Iterator A %Python iterator type.
 *
 * @param _self The property manager instance that is to be used.
 * @param _it An iterator that iterates over the items in the range.
 * @param _value The value the range will be set to.
 */
template <class PropertyManager, class Iterator>
void propman_set_range(PropertyManager& _self, Iterator _it, object _value) {
	try {
		while (true) {
			_self[_it.next()] = _value;
		}
	}
	catch (error_already_set exception) {
		// This is expected behavior
		PyErr_Clear();
	}
}

/**
 * Thin wrapper for propertyExists.
 *
 * @tparam PropertyManager A property manager type.
 * @tparam Mesh A mesh type.
 *
 * @param _mesh The mesh that is used to check if the property exists.
 * @param _propname The name of the property.
 */
template <class PropertyManager, class Mesh>
bool property_exists(Mesh& _mesh, const char *_propname) {
	return PropertyManager::propertyExists(_mesh, _propname);
}

/**
 * Expose a property manager type to %Python.
 *
 * This function template is used to expose property managers to %Python. The
 * template parameters are used to instantiate the appropriate property manager
 * type.
 *
 * @tparam PropHandle A property handle type (e.g. %VPropHandle\<object\>).
 * @tparam IndexHandle The appropriate handle type (e.g. %VertexHandle for
 * %VPropHandle\<object\>).
 * @tparam Iterator A %Python iterator type. This type is used to instantiate
 * the propman_set_range function.
 *
 * @param _name The name of the property manager type to be exposed.
 */
template <class PropHandle, class IndexHandle, class Iterator>
void expose_property_manager(const char *_name) {
	// Convenience typedef
	typedef PropertyManager<PropHandle, PolyConnectivity> PropertyManager;

	// Function pointers
	void (PropertyManager::*retain)(bool) = &PropertyManager::retain;

	object (*getitem)(PropertyManager&, IndexHandle        ) = &propman_get_item;
	void   (*setitem)(PropertyManager&, IndexHandle, object) = &propman_set_item;

	void (*set_range)(PropertyManager&, Iterator, object) = &propman_set_range;

	bool (*property_exists_poly)(PolyMesh&, const char *) = &property_exists<PropertyManager, PolyMesh>;
	bool (*property_exists_tri )(TriMesh&,  const char *) = &property_exists<PropertyManager, TriMesh >;

	// Expose property manager
	class_<PropertyManager, boost::noncopyable>(_name)
		.def(init<PolyMesh&, const char *, optional<bool> >())
		.def(init<TriMesh&,  const char *, optional<bool> >())

		.def("swap", &PropertyManager::swap)
		.def("is_valid", &PropertyManager::isValid)

		.def("__bool__", &PropertyManager::operator bool)
		.def("__nonzero__", &PropertyManager::operator bool)

		.def("get_raw_property", &PropertyManager::getRawProperty, return_value_policy<copy_const_reference>())
		.def("get_name", &PropertyManager::getName, return_value_policy<copy_const_reference>())
		.def("get_mesh", &PropertyManager::getMesh, return_value_policy<reference_existing_object>())

		.def("retain", retain, retain_overloads())

		.def("__getitem__", getitem)
		.def("__setitem__", setitem)

		.def("set_range", set_range)

		.def("property_exists", property_exists_poly)
		.def("property_exists", property_exists_tri)
		.staticmethod("property_exists")
		;
}

} // namespace OpenMesh
} // namespace Python

#endif
