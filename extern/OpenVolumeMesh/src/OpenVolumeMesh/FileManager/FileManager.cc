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
 *   $Date$                   *
 *   $LastChangedBy$                                                *
 *                                                                           *
\*===========================================================================*/

#define FILEMANAGERT_CC

#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <typeinfo>

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

#include "FileManager.hh"

namespace OpenVolumeMesh {

namespace IO {

//==================================================

FileManager::FileManager() {

}

//==================================================

FileManager::~FileManager() {

}

//==================================================

void FileManager::trimString(std::string& _string) const {

    // Trim Both leading and trailing spaces
    size_t start = _string.find_first_not_of(" \t\r\n");
    size_t end = _string.find_last_not_of(" \t\r\n");

    if((std::string::npos == start) || (std::string::npos == end)) {
        _string = "";
    } else {
        _string = _string.substr(start, end - start + 1);
    }
}

//==================================================

void FileManager::extractQuotedText(std::string& _string) const {

    // Trim Both leading and trailing quote marks
    size_t start = _string.find_first_of("\""); ++start;
    size_t end = _string.find_last_not_of("\"");

    if((std::string::npos == start) || (std::string::npos == end)) {
        _string = "";
    } else {
        _string = _string.substr(start, end - start + 1);
    }
}

//==================================================

bool FileManager::getCleanLine(std::istream& _ifs, std::string& _string, bool _skipEmptyLines) const {

    // While we are not at the end of the file
    while(true) {

        // Get the current line:
        std::getline(_ifs, _string);

        // Remove whitespace at beginning and end
        trimString(_string);

        // Check if string is not empty ( otherwise we continue
        if(_string.size() != 0) {

            // Check if string is a comment ( starting with # )
            if(_string[0] != '#') {
                return true;
            }

        } else {
            if(!_skipEmptyLines)
                return true;
        }

        if(_ifs.eof()) {
            //std::cerr << "End of file reached while searching for input!" << std::endl;
            return false;
        }
    }

    return false;
}

//==================================================

bool FileManager::isHexahedralMesh(const std::string& _filename) const {

  std::ifstream iff(_filename.c_str(), std::ios::in);

  if(!iff.good()) {
    std::cerr << "Could not open file " << _filename << " for reading!" << std::endl;
    iff.close();
    return false;
  }

  std::string s;
  unsigned int n = 0u;

  // Skip until we find polyhedra section
  while (true || !iff.eof()) {
    iff >> s;
    if (s == "Polyhedra") {
      break;
    }
  }

  if (iff.eof()) {
    // Polyhedra section not found in file. Defaulting to polyhedral type.
    iff.close();
    return false;
  }

  // Read in number of cells
  iff >> n;
  if(n == 0) return false;
  unsigned int v = 0;
  char tmp[256];
  for (unsigned int i = 0; i < n; ++i) {
    iff >> v;
    iff.getline(tmp, 256);
    if (v != 6u) {
      iff.close();
      return false;
    }
  }
  iff.close();
  return true;
}

//==================================================

} // Namespace IO
} // Namespace OpenVolumeMesh
