#ifndef OPENMESH_PYTHON_CIRCULATOR_HH
#define OPENMESH_PYTHON_CIRCULATOR_HH

#include "Python/Bindings.hh"

namespace OpenMesh {
namespace Python {

/**
 * Wrapper for circulators.
 *
 * This class template is used to wrap circulators for %Python. It implements
 * %Python's iterator protocol (the magic methods \_\_iter\_\_ and
 * \_\_next\_\_).
 *
 * @tparam Circulator A circulator type.
 */
template<class Circulator, class CenterEntityHandle>
class CirculatorWrapperT {
	public:

		/**
		 * Constructor
		 *
		 * @param _mesh The mesh that contains the items to iterate over.
		 * @param _center The handle to the center item.
		 */
		CirculatorWrapperT(PolyMesh& _mesh, CenterEntityHandle _center) :
			circulator_(_mesh, _center) {
		}

		/**
		 * Constructor
		 *
		 * @param _mesh The mesh that contains the items to iterate over.
		 * @param _center The handle to the center item.
		 */
		CirculatorWrapperT(TriMesh& _mesh, CenterEntityHandle _center) :
			circulator_(_mesh, _center) {
		}

		/**
		 * Implementation of %Python's \_\_iter\_\_ magic method.
		 *
		 * @return This circulator.
		 */
		CirculatorWrapperT iter() const {
			return *this;
		}

		/**
		 * Implementation of %Python's \_\_next\_\_ magic method.
		 *
		 * @return The next item. Raises a %Python StopIteration exception if
		 * there are no more items.
		 */
		typename Circulator::value_type next() {
			if (circulator_.is_valid()) {
				typename Circulator::value_type res = *circulator_;
				++circulator_;
				return res;
			}
			else {
				PyErr_SetString(PyExc_StopIteration, "No more data.");
				boost::python::throw_error_already_set();
			}
			return typename Circulator::value_type();
		}

	private:
		Circulator circulator_;
};

/**
 * Expose a circulator type to %Python.
 *
 * @tparam Circulator A circulator type.
 *
 * @param _name The name of the circulator type to be exposed.
 *
 * @note Circulators are wrapped by CirculatorWrapperT before they are exposed
 * to %Python, i.e. they are not exposed directly. This means that circulators
 * that are passed from %Python to C++ are instances of CirculatorWrapperT.
 */
template<class Circulator, class CenterEntityHandle>
void expose_circulator(const char *_name) {
	class_<CirculatorWrapperT<Circulator, CenterEntityHandle> >(_name, init<TriMesh&, CenterEntityHandle>())
		.def(init<PolyMesh&, CenterEntityHandle>())
		.def("__iter__", &CirculatorWrapperT<Circulator, CenterEntityHandle>::iter)
		.def("__next__", &CirculatorWrapperT<Circulator, CenterEntityHandle>::next)
		.def("next", &CirculatorWrapperT<Circulator, CenterEntityHandle>::next)
		;
}

} // namespace OpenMesh
} // namespace Python

#endif
