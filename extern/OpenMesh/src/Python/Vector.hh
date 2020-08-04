#ifndef OPENMESH_PYTHON_VECTOR_HH
#define OPENMESH_PYTHON_VECTOR_HH

#include "Python/Bindings.hh"

namespace OpenMesh {
namespace Python {

template <class Vector, class Scalar>
void set_item(Vector& _vec, int _index, Scalar _value) {
	if (_index < 0) {
		_index += _vec.size();
	}

	if ((size_t)_index < _vec.size()) {
		_vec[_index] = _value;
	}
	else {
		PyErr_SetString(PyExc_IndexError, "Index out of range.");
		throw_error_already_set();
	}
}

template <class Vector, class Scalar>
Scalar get_item(Vector& _vec, int _index) {
	if (_index < 0) {
		_index += _vec.size();
	}

	if ((size_t)_index < _vec.size()) {
		return _vec[_index];
	}
	else {
		PyErr_SetString(PyExc_IndexError, "Index out of range.");
		throw_error_already_set();
	}

	return 0.0;
}

/**
 * Expose a vector type to %Python.
 *
 * This function template is used to expose vectors to %Python. The template
 * parameters are used to instantiate the appropriate vector type.
 *
 * @tparam Scalar A scalar type.
 * @tparam N The dimension of the vector.
 *
 * @param _name The name of the vector type to be exposed.
 *
 * @note N must be either 2, 3 or 4.
 */
template<class Scalar, int N>
void expose_vec(const char *_name) {
	typedef OpenMesh::VectorT<Scalar, N> Vector;

	Scalar (Vector::*min_void)() const = &Vector::min;
	Scalar (Vector::*max_void)() const = &Vector::max;

	Vector (Vector::*max_vector)(const Vector&) const = &Vector::max;
	Vector (Vector::*min_vector)(const Vector&) const = &Vector::min;

	class_<Vector> classVector = class_<Vector>(_name);

	classVector
		.def("__setitem__", &set_item<Vector, Scalar>)
		.def("__getitem__", &get_item<Vector, Scalar>)
		.def(self == self)
		.def(self != self)
		.def(self *= Scalar())
		.def(self /= Scalar())
		.def(self * Scalar())
		.def(Scalar() * self)
		.def(self / Scalar())
		.def(self *= self)
		.def(self /= self)
		.def(self -= self)
		.def(self += self)
		.def(self * self)
		.def(self / self)
		.def(self + self)
		.def(self - self)
		.def(-self)
		.def(self | self)
		.def("vectorize", &Vector::vectorize, return_internal_reference<>())
		.def(self < self)

		.def("norm", &Vector::norm)
		.def("length", &Vector::length)
		.def("sqrnorm", &Vector::sqrnorm)
		.def("normalize", &Vector::normalize, return_internal_reference<>())
		.def("normalized", &Vector::normalized)
		.def("normalize_cond", &Vector::normalize_cond, return_internal_reference<>())

		.def("l1_norm", &Vector::l1_norm)
		.def("l8_norm", &Vector::l8_norm)

		.def("max", max_void)
		.def("max_abs", &Vector::max_abs)
		.def("min", min_void)
		.def("min_abs", &Vector::min_abs)
		.def("mean", &Vector::mean)
		.def("mean_abs", &Vector::mean_abs)
		.def("minimize", &Vector::minimize, return_internal_reference<>())
		.def("minimized", &Vector::minimized)
		.def("maximize", &Vector::maximize, return_internal_reference<>())
		.def("maximized", &Vector::maximized)
		.def("min", min_vector)
		.def("max", max_vector)

		.def("size", &Vector::size)
		.staticmethod("size")
		.def("vectorized", &Vector::vectorized)
		.staticmethod("vectorized")
		;

	typedef OpenMesh::VectorT<Scalar, 2> Vector2;
	typedef OpenMesh::VectorT<Scalar, 3> Vector3;
	typedef OpenMesh::VectorT<Scalar, 4> Vector4;

	struct Factory {
		static Vector2 *vec2_default() {
			return new Vector2(Scalar(), Scalar());
		}
		static Vector2 *vec2_user_defined(const Scalar& _v0, const Scalar& _v1) {
			return new Vector2(_v0, _v1);
		}
		static Vector3 *vec3_default() {
			return new Vector3(Scalar(), Scalar(), Scalar());
		}
		static Vector3 *vec3_user_defined(const Scalar& _v0, const Scalar& _v1, const Scalar& _v2) {
			return new Vector3(_v0, _v1, _v2);
		}
		static Vector4 *vec4_default() {
			return new Vector4(Scalar(), Scalar(), Scalar(), Scalar());
		}
		static Vector4 *vec4_user_defined(const Scalar& _v0, const Scalar& _v1, const Scalar& _v2, const Scalar& _v3) {
			return new Vector4(_v0, _v1, _v2, _v3);
		}
	};

	if (N == 2) {
		classVector
			.def("__init__", make_constructor(&Factory::vec2_default))
			.def("__init__", make_constructor(&Factory::vec2_user_defined))
			;
	}

	if (N == 3) {
		classVector
			.def("__init__", make_constructor(&Factory::vec3_default))
			.def("__init__", make_constructor(&Factory::vec3_user_defined))
			.def("__mod__", &Vector3::operator%)
			;

		def("cross", &Vector3::operator%);
	}

	if (N == 4) {
		classVector
			.def("__init__", make_constructor(&Factory::vec4_default))
			.def("__init__", make_constructor(&Factory::vec4_user_defined))
			;
	}

	def("dot", &Vector::operator|);
}

} // namespace OpenMesh
} // namespace Python

#endif
