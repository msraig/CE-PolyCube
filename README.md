# CE-PolyCube: Cut-enhanced PolyCube-Maps for Feature-aware All-Hex Meshing
By Haoxiang Guo, Xiaohan Liu, Dong-Ming Yan, Yang Liu

The program works in 64-bit windows.
### Prerequisite
+ Boost (1.66): https://sourceforge.net/projects/boost/files/boost-binaries/
+ CGAL (4.13): https://www.cgal.org/releases.html
+ MOSEK 8: https://www.mosek.com/downloads/
+ VTK (8.2): https://vtk.org/download/
+ Eigen3: http://eigen.tuxfamily.org/index.php?title=Main_Page#Download
+ PolyCut: http://www.cs.ubc.ca/labs/imager/tr/2018/HexDemo/, download and extract it to a folder, make sure the one-month license has not expired.
+ HexEx: https://www.graphics.rwth-aachen.de/media/resource_files/HexEx_Windows_1_01.zip, download and extract it to a folder.

The libraries and softwares above need to be installed **manually**, please set BOOST_ROOT CGAL_DIR EIGEN3_INCLUDE_DIR in your system path. Besides, our program has the following build-in dependencies:
+ OpenMesh: https://www.graphics.rwth-aachen.de/software/openmesh/
+ OpenVolumeMesh: https://www.graphics.rwth-aachen.de/software/openvolumemesh/, we modify the source code slightly.
+ gflags: https://github.com/gflags/gflags.git
+ rapidxml: http://rapidxml.sourceforge.net/
+ SUS: https://www.dca.iusiani.ulpgc.es/SUScode/, we modify the source code slightly.
+ ANN: http://www.cs.umd.edu/~mount/ANN/
+ TetGen: http://wias-berlin.de/software/index.jsp?id=TetGen&lang=1, we modify the source code slightly.
+ SuiteSparse: http://faculty.cse.tamu.edu/davis/suitesparse.html

### Data Preparation
Our program takes a tetrahedral mesh (\*.vtk) and its feature edges (\*.fea) as input.  The feature file starts with the number of feature segments **n**, followed by **n** lines, where each line contains the two vertex indices of a feature edge segment. You can generate tetrahedral meshes from surface meshes using [TetGen](http://wias-berlin.de/software/index.jsp?id=TetGen&lang=1) or download our preprocessed data from: https://drive.google.com/open?id=1KUysv4rmUIX1F5nnHKHaek5UboeCclE8

### Installation & Usage
First, clone this repository with its submodules:
```
git clone https://github.com/guohaoxiang/CE-PolyCube.git --recurse-submodules
cd CE-PolyCube
```
Apply patch changes:
```
sh ./apply_patch.sh
```
Run CMake:
```
mkdir Build
cd Build
cmake ..
```
Build the solution in .\Build folder, then you can find the generated executable files in .\Bin folder.
We also provide a script for automating our pipeline. To use it, please first set DATA_ROOT_PATH (path contains model folders), HEXEX_PATH (path containing hexex.exe) and POLYCUT_PATH (path containing polycut.exe) in .\Script\gen_hex.bat with no trailing slash. Then you can process a specific model, e.g sculpt by running:
```
cd ../Scripts
gen_hex.bat sculpt
```
There will be three new files in DATA_ROOT_PATH/sculpt folder:
+ sculpt_deform_polycube.vtk: the initial PolyCube mesh.
+ sculpt_cut_flattening.vtk: the CE-PolyCube mesh.
+ sculpt_hex_opt.vtk: the final optimized all-hex mesh.

### Cite(to do)
+ Our paper
+ SUS
+ PolyCut
+ HexEx

### License(to do)

