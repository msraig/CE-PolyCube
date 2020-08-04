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


#ifndef OPENVOLUMEMESHBASEPROPERTY_HH
#define OPENVOLUMEMESHBASEPROPERTY_HH

#include <limits>
#include <string>
#include <iostream>
#include <vector>

#include "OpenVolumeMeshHandle.hh"

//== CLASS DEFINITION =========================================================

/** \class OpenVolumeMeshBaseProperty

 Abstract class defining the basic interface of a dynamic property.

 **/

namespace OpenVolumeMesh {

class OpenVolumeMeshBaseProperty {
public:

    friend class ResourceManager;
    template <class PropT, class HandleT> friend class PropertyPtr;

	/// Indicates an error when a size is returned by a member.
	static const size_t UnknownSize;

public:

	OpenVolumeMeshBaseProperty(const std::string& _name = "<unknown>") :
		name_(_name), persistent_(false), handle_(-1) {
	}

	OpenVolumeMeshBaseProperty(const OpenVolumeMeshBaseProperty& _rhs) :
		name_(_rhs.name_), persistent_(_rhs.persistent_), handle_(_rhs.handle_.idx()) {
	}

	virtual ~OpenVolumeMeshBaseProperty() {}

public:

	/// Reserve memory for n elements.
	virtual void reserve(size_t _n) = 0;

	/// Resize storage to hold n elements.
	virtual void resize(size_t _n) = 0;

	/// Clear all elements and free memory.
	virtual void clear() = 0;

	/// Extend the number of elements by one.
	virtual void push_back() = 0;

	/// Let two elements swap their storage place.
	virtual void swap(size_t _i0, size_t _i1) = 0;

	/// Erase an element of the vector
	virtual void delete_element(size_t _idx) = 0;

	/// Return a deep copy of self.
	virtual OpenVolumeMeshBaseProperty* clone() const = 0;

	/// Return the name of the property
	const std::string& name() const {
		return name_;
	}

	// Function to serialize a property
	virtual void serialize(std::ostream& _ostr) const {
	    _ostr << "\"" << name_ << "\"" << std::endl;
	}

	// Function to deserialize a property
    virtual void deserialize(std::istream& /*_istr*/) {}
	// I/O support

	void set_persistent(bool _persistent) { persistent_ = _persistent; }

	bool persistent() const { return persistent_; }

	/// Number of elements in property
	virtual size_t n_elements() const = 0;

	/// Size of one element in bytes or UnknownSize if not known.
	virtual size_t element_size() const = 0;

	/// Return size of property in bytes
	virtual size_t size_of() const {
		return size_of(n_elements());
	}

	/// Estimated size of property if it has _n_elem elements.
	/// The member returns UnknownSize if the size cannot be estimated.
	virtual size_t size_of(size_t _n_elem) const {
		return (element_size() != UnknownSize) ? (_n_elem * element_size())
				: UnknownSize;
	}

	const OpenVolumeMeshHandle& handle() const { return handle_; }

	void set_handle(const OpenVolumeMeshHandle& _handle) { handle_.idx(_handle.idx()); }

protected:

	/// Delete multiple entries in list
    virtual void delete_multiple_entries(const std::vector<bool>&) = 0;

private:

	std::string name_;

	bool persistent_;

	OpenVolumeMeshHandle handle_;
};

} // Namespace OpenVolumeMesh

#endif //OPENVOLUMEMESHBASEPROPERTY_HH

