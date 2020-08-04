/*===========================================================================*\
 *                                                                           *
 *                            OpenVolumeMesh                                 *
 *        Copyright (C) 2011 by Computer Graphics Group, RWTH Aachen         *
 *                        www.openvolumemesh.org                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of OpenVolumeMesh.                                     *
 *                                                                           *
 *  OpenVolumeMesh is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenVolumeMesh is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenVolumeMesh.  If not,                              *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *
 *   $Revision$                                                         *
 *   $Date$                    *
 *   $LastChangedBy$                                                *
 *                                                                           *
\*===========================================================================*/

#ifndef OPENVOLUMEMESHPROPERTY_HH
#define OPENVOLUMEMESHPROPERTY_HH

//== INCLUDES =================================================================

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cassert>

#include "OpenVolumeMeshBaseProperty.hh"
#include "OpenVolumeMeshHandle.hh"

namespace OpenVolumeMesh {

//== CLASS DEFINITION =========================================================

/** \class OpenVolumeMeshPropertyT
 *
 *  \brief Default property class for any type T.
 *
 *  The default property class for any type T.
 */

template<class T>
class OpenVolumeMeshPropertyT: public OpenVolumeMeshBaseProperty {
public:

    template <class PropT, class HandleT> friend class PropertyPtr;

	typedef T 										Value;
	typedef std::vector<T> 				            vector_type;
	typedef T 										value_type;
	typedef typename vector_type::reference 		reference;
	typedef typename vector_type::const_reference 	const_reference;

public:

	/// Default constructor
	OpenVolumeMeshPropertyT(const std::string& _name = "<unknown>", const T _def = T()) :
		OpenVolumeMeshBaseProperty(_name),
		def_(_def) {
	}

	/// Copy constructor
	OpenVolumeMeshPropertyT(const OpenVolumeMeshPropertyT& _rhs) :
		OpenVolumeMeshBaseProperty(_rhs),
		data_(_rhs.data_),
		def_(_rhs.def_) {
	}

public:
	// inherited from OpenVolumeMeshBaseProperty
	virtual void reserve(size_t _n) {
		data_.reserve(_n);
	}
	virtual void resize(size_t _n) {
		data_.resize(_n, def_);
	}
	virtual void clear() {
		data_.clear();
		vector_type().swap(data_);
	}
	virtual void push_back() {
		data_.push_back(def_);
	}
	virtual void swap(size_t _i0, size_t _i1) {
		std::swap(data_[_i0], data_[_i1]);
	}

	void delete_element(size_t _idx) {
	    data_.erase(data_.begin() + _idx);
	}

public:

	virtual size_t n_elements() const {
		return data_.size();
	}
	virtual size_t element_size() const {
		return sizeof(T);
	}

#ifndef DOXY_IGNORE_THIS
	struct plus {
		size_t operator ()(size_t _b, const T& /*_v*/) {
			return _b + sizeof(T);
		}
	};
#endif

	virtual size_t size_of() const {
		if (element_size() != OpenVolumeMeshBaseProperty::UnknownSize)
			return this->OpenVolumeMeshBaseProperty::size_of(n_elements());
		return std::accumulate(data_.begin(), data_.end(), 0, plus());
	}

	virtual size_t size_of(size_t _n_elem) const {
		return this->OpenVolumeMeshBaseProperty::size_of(_n_elem);
	}

	// Function to serialize a property
    virtual void serialize(std::ostream& _ostr) const {
        OpenVolumeMeshBaseProperty::serialize(_ostr);
        for(typename vector_type::const_iterator it = data_.begin();
                it != data_.end(); ++it) {
            _ostr << *it << std::endl;
        }
    }

    // Function to deserialize a property
    virtual void deserialize(std::istream& _istr) {
        for(unsigned int i = 0; i < n_elements(); ++i) {
            value_type val;
            _istr >> val;
            data_[i] = val;
        }
    }

public:
	// data access interface

	/// Get pointer to array (does not work for T==bool)
	const T* data() const {

		if (data_.empty())
			return 0;

		return &data_[0];
	}

	/// Get reference to property vector (be careful, improper usage, e.g. resizing, may crash)
	vector_type& data_vector() {

		return data_;
	}

	/// Access the i'th element. No range check is performed!
	reference operator[](int _idx) {
		assert(size_t(_idx) < data_.size());
		return data_[_idx];
	}

	/// Const access to the i'th element. No range check is performed!
	const_reference operator[](int _idx) const {
		assert(size_t(_idx) < data_.size());
		return data_[_idx];
	}

	/// Make a copy of self.
	OpenVolumeMeshPropertyT<T>* clone() const {
		OpenVolumeMeshPropertyT<T>* p = new OpenVolumeMeshPropertyT<T>(*this);
		return p;
	}

	typename vector_type::const_iterator begin() const { return data_.begin(); }

	typename vector_type::iterator begin() { return data_.begin(); }

	typename vector_type::const_iterator end() const { return data_.end(); }

    typename vector_type::iterator end() { return data_.end(); }

protected:

    /// Delete multiple entries in list
    virtual void delete_multiple_entries(const std::vector<bool>& _tags) {

        assert(_tags.size() == data_.size());
        vector_type new_data;
        typename vector_type::iterator d_it = data_.begin();
        std::vector<bool>::const_iterator t_it = _tags.begin();
        std::vector<bool>::const_iterator t_end = _tags.end();
        for(; t_it != t_end; ++t_it, ++d_it) {
            if(!*t_it) {
                new_data.push_back(*d_it);
            }
        }
        data_.swap(new_data);
    }

private:

	vector_type data_;

	const T def_;
};

//-----------------------------------------------------------------------------


/**
 * Property specialization for bool type.
 */
template<>
class OpenVolumeMeshPropertyT<bool> : public OpenVolumeMeshBaseProperty {
public:

    template <class PropT, class HandleT> friend class PropertyPtr;

	typedef std::vector<bool> 				vector_type;
	typedef bool 							value_type;
	typedef vector_type::reference 			reference;
	typedef vector_type::const_reference 	const_reference;

public:

	OpenVolumeMeshPropertyT(const std::string& _name = "<unknown>", const bool _def = bool()) :
		OpenVolumeMeshBaseProperty(_name),
		def_(_def) {
	}

public:
	// inherited from OpenVolumeMeshBaseProperty

	virtual void reserve(size_t _n) {
		data_.reserve(_n);
	}
	virtual void resize(size_t _n) {
		data_.resize(_n, def_);
	}
	virtual void clear() {
		data_.clear();
		vector_type().swap(data_);
	}
	virtual void push_back() {
		data_.push_back(def_);
	}
	virtual void swap(size_t _i0, size_t _i1) {
		bool t(data_[_i0]);
		data_[_i0] = data_[_i1];
		data_[_i1] = t;
	}

	void delete_element(size_t _idx) {
        data_.erase(data_.begin() + _idx);
    }

public:

	virtual size_t n_elements() const {
		return data_.size();
	}
	virtual size_t element_size() const {
		return OpenVolumeMeshBaseProperty::UnknownSize;
	}
	virtual size_t size_of() const {
		return size_of(n_elements());
	}
	virtual size_t size_of(size_t _n_elem) const {
		return _n_elem / 8 + ((_n_elem % 8) != 0);
	}

public:

	/// Access the i'th element. No range check is performed!
	reference operator[](int _idx) {
		assert(size_t(_idx) < data_.size());
		return data_[_idx];
	}

	/// Const access to the i'th element. No range check is performed!
	const_reference operator[](int _idx) const {
		assert(size_t(_idx) < data_.size());
		return data_[_idx];
	}

	/// Make a copy of self.
	OpenVolumeMeshPropertyT<bool>* clone() const {
		OpenVolumeMeshPropertyT<bool>* p = new OpenVolumeMeshPropertyT<bool> (
				*this);
		return p;
	}

	vector_type::const_iterator begin() const { return data_.begin(); }

    vector_type::iterator begin() { return data_.begin(); }

    vector_type::const_iterator end() const { return data_.end(); }

    vector_type::iterator end() { return data_.end(); }

protected:

    /// Delete multiple entries in list
    virtual void delete_multiple_entries(const std::vector<bool>& _tags) {

        assert(_tags.size() == data_.size());
        vector_type new_data;
        vector_type::iterator d_it = data_.begin();
        std::vector<bool>::const_iterator t_it = _tags.begin();
        std::vector<bool>::const_iterator t_end = _tags.end();
        for(; t_it != t_end; ++t_it, ++d_it) {
            if(!*t_it) {
                new_data.push_back(*d_it);
            }
        }
        data_.swap(new_data);
    }

private:

	vector_type data_;

	const bool def_;
};

//-----------------------------------------------------------------------------


/**
 * Property specialization for std::string type.
 */
template<>
class OpenVolumeMeshPropertyT<std::string> : public OpenVolumeMeshBaseProperty {
public:

    template <class PropT, class HandleT> friend class PropertyPtr;

	typedef std::string 					Value;
	typedef std::vector<std::string> 		vector_type;
	typedef std::string 					value_type;
	typedef vector_type::reference 			reference;
	typedef vector_type::const_reference 	const_reference;

public:

	OpenVolumeMeshPropertyT(const std::string& _name = "<unknown>", const std::string _def = std::string()) :
		OpenVolumeMeshBaseProperty(_name),
		def_(_def) {
	}

public:
	// inherited from OpenVolumeMeshBaseProperty

	virtual void reserve(size_t _n) {
		data_.reserve(_n);
	}
	virtual void resize(size_t _n) {
		data_.resize(_n, def_);
	}
	virtual void clear() {
		data_.clear();
		vector_type().swap(data_);
	}
	virtual void push_back() {
		data_.push_back(def_);
	}
	virtual void swap(size_t _i0, size_t _i1) {
		std::swap(data_[_i0], data_[_i1]);
	}

	virtual void delete_element(size_t _idx) {
        data_.erase(data_.begin() + _idx);
    }

public:

	virtual size_t n_elements() const {
		return data_.size();
	}
	virtual size_t element_size() const {
		return OpenVolumeMeshBaseProperty::UnknownSize;
	}
	virtual size_t size_of() const {
		return sizeof(data_);
	}

	virtual size_t size_of(size_t /* _n_elem */) const {
		return OpenVolumeMeshBaseProperty::UnknownSize;
	}

	virtual void stats(std::ostream& _ostr) const {
		for(vector_type::const_iterator it = data_.begin();
			it != data_.end(); ++it) {
				_ostr << *it << " ";
		}
	}

public:

	const value_type* data() const {
		if (data_.empty())
			return 0;

		return (value_type*) &data_[0];
	}

	/// Access the i'th element. No range check is performed!
	reference operator[](int _idx) {
		assert(size_t(_idx) < data_.size());
		return ((value_type*) &data_[0])[_idx];
	}

	/// Const access the i'th element. No range check is performed!
	const_reference operator[](int _idx) const {
		assert(size_t(_idx) < data_.size());
		return ((value_type*) &data_[0])[_idx];
	}

	OpenVolumeMeshPropertyT<value_type>* clone() const {
		OpenVolumeMeshPropertyT<value_type>* p = new OpenVolumeMeshPropertyT<
				value_type> (*this);
		return p;
	}

	vector_type::const_iterator begin() const { return data_.begin(); }

    vector_type::iterator begin() { return data_.begin(); }

    vector_type::const_iterator end() const { return data_.end(); }

    vector_type::iterator end() { return data_.end(); }

protected:

    /// Delete multiple entries in list
    virtual void delete_multiple_entries(const std::vector<bool>& _tags) {

        assert(_tags.size() == data_.size());
        vector_type new_data;
        vector_type::iterator d_it = data_.begin();
        std::vector<bool>::const_iterator t_it = _tags.begin();
        std::vector<bool>::const_iterator t_end = _tags.end();
        for(; t_it != t_end; ++t_it, ++d_it) {
            if(!*t_it) {
                new_data.push_back(*d_it);
            }
        }
        data_.swap(new_data);
    }

private:

	vector_type data_;

	const std::string def_;
};

} // Namespace OpenVolumeMesh

//=============================================================================
#endif // OPENVOLUMEMESHPROPERTY_HH defined
//=============================================================================
