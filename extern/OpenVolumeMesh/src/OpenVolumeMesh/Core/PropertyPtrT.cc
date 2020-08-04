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

#define PROPERTYPTRT_CC

#include "PropertyPtr.hh"
#include "ResourceManager.hh"
#include "PropertyDefines.hh"

namespace OpenVolumeMesh {

template <class PropT, class HandleT>
PropertyPtr<PropT,HandleT>::PropertyPtr(PropT* _ptr, ResourceManager& _resMan, HandleT _handle) :
    ptr::shared_ptr<PropT>(_ptr), BaseProperty(_resMan) {
    ptr::shared_ptr<PropT>::get()->set_handle(_handle);
}

template <class PropT, class HandleT>
PropertyPtr<PropT,HandleT>::~PropertyPtr() {

    /*
     * If use count is 2 and prop is not set persistent,
     * remove it, since the resource manager is the
     * only one who stores the property.
     */
    if(!locked() && !persistent() && ptr::shared_ptr<PropT>::use_count() == 2) {
        resMan_.release_property(HandleT(handle().idx()));
        unlock();
    }
}

template <class PropT, class HandleT>
void PropertyPtr<PropT,HandleT>::resize(unsigned int _size) {
    ptr::shared_ptr<PropT>::get()->resize(_size);
}

template <class PropT, class HandleT>
const std::string& PropertyPtr<PropT,HandleT>::name() const {
    return ptr::shared_ptr<PropT>::get()->name();
}

template <class PropT, class HandleT>
void PropertyPtr<PropT,HandleT>::delete_element(size_t _idx) {
    ptr::shared_ptr<PropT>::get()->delete_element(_idx);
}

template <class PropT, class HandleT>
void PropertyPtr<PropT,HandleT>::set_handle(const OpenVolumeMeshHandle& _handle) {
    return ptr::shared_ptr<PropT>::get()->set_handle(_handle);
}

template <class PropT, class HandleT>
OpenVolumeMeshHandle PropertyPtr<PropT,HandleT>::handle() const {
    return ptr::shared_ptr<PropT>::get()->handle();
}

template <class PropT, class HandleT>
void PropertyPtr<PropT,HandleT>::delete_multiple_entries(const std::vector<bool>& _tags) {
    ptr::shared_ptr<PropT>::get()->delete_multiple_entries(_tags);
}

} // Namespace OpenVolumeMesh
