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

#ifndef PROPERTYPTR_HH_
#define PROPERTYPTR_HH_

#include <iostream>
#include <string>
#include <cassert>

#include "BaseProperty.hh"
#include "PropertyHandles.hh"

#include "../System/MemoryInclude.hh"

namespace OpenVolumeMesh {

class ResourceManager;

/**
 * \class PropertyPtr
 *
 * A smart-pointer-like class that counts the encapsulated
 * object's references and automatically deletes the memory
 * as soon as the object is not used anymore.
 */

template <class PropT, class HandleT>
class PropertyPtr : public ptr::shared_ptr<PropT>, public BaseProperty {
public:

    friend class ResourceManager;

    typedef typename PropT::value_type                  value_type;
    typedef typename PropT::vector_type::const_iterator const_iterator;
    typedef typename PropT::vector_type::iterator       iterator;
    typedef typename PropT::reference                   reference;
    typedef typename PropT::const_reference             const_reference;

    /// Constructor
    PropertyPtr(PropT* _ptr, ResourceManager& _resMan, HandleT _handle);

    /// Destructor
    virtual ~PropertyPtr();

    virtual const std::string& name() const;

    virtual void delete_element(size_t _idx);

    const_iterator begin() const { return ptr::shared_ptr<PropT>::get()->begin(); }
    iterator begin() { return ptr::shared_ptr<PropT>::get()->begin(); }

    const_iterator end() const { return ptr::shared_ptr<PropT>::get()->end(); }
    iterator end() { return ptr::shared_ptr<PropT>::get()->end(); }

    reference operator[](size_t _idx) { return (*ptr::shared_ptr<PropT>::get())[_idx]; }
    const_reference operator[](size_t _idx) const { return (*ptr::shared_ptr<PropT>::get())[_idx]; }

    reference operator[](const OpenVolumeMeshHandle& _h) { return (*ptr::shared_ptr<PropT>::get())[_h.idx()]; }
    const_reference operator[](const OpenVolumeMeshHandle& _h) const { return (*ptr::shared_ptr<PropT>::get())[_h.idx()]; }

    virtual OpenVolumeMeshHandle handle() const;

    virtual bool persistent() const { return ptr::shared_ptr<PropT>::get()->persistent(); }

    virtual bool anonymous() const { return ptr::shared_ptr<PropT>::get()->name().empty(); }

protected:

    virtual void delete_multiple_entries(const std::vector<bool>& _tags);

    virtual void resize(unsigned int _size);

    virtual void set_handle(const OpenVolumeMeshHandle& _handle);
};

} // Namespace OpenVolumeMesh

#if defined(INCLUDE_TEMPLATES) && !defined(PROPERTYPTRT_CC)
#include "PropertyPtrT.cc"
#endif

#endif /* PROPERTYPTR_HH_ */
