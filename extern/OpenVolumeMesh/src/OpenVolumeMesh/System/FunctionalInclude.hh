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
 *   $Revision: 36 $                                                         *
 *   $Date: 2012-01-10 18:00:06 +0100 (Di, 10 Jan 2012) $                    *
 *   $LastChangedBy: kremer $                                                *
 *                                                                           *
\*===========================================================================*/

#ifndef FUNCTIONALINCLUDE_HH_
#define FUNCTIONALINCLUDE_HH_

#if (__cplusplus >= 201103L)
   // C++11:
   #include <functional>
   namespace fun = std;
#elif defined(__GXX_EXPERIMENTAL_CXX0X__)
   // C++11 via -std=c++0x on gcc:
   #include <functional>
   namespace fun = std;
#else
   // C++98 and TR1:
   #if (_MSC_VER >= 1600)
     // VStudio 2010 supports some C++11 features
     #include <functional>
     namespace fun = std;
   #elif (_MSC_VER >= 1500)
     // hope for TR1 equivalents
     #include <functional>
	 namespace fun = std::tr1;
   #else
     // hope for TR1 equivalents
     #include <tr1/functional>
     namespace fun = std::tr1;
   #endif
#endif

#endif /* FUNCTIONALINCLUDE_HH_ */
