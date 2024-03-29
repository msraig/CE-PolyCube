# CE-PolyCube: Cut-enhanced PolyCube-Maps for Feature-aware All-Hex Meshing
Haoxiang Guo, Xiaohan Liu, Dong-Ming Yan, Yang Liu.
ACM Transactions on Graphics (SIGGRAPH2020).

If you are interested in the algorithm details, please refer to our **[paper](https://app.box.com/s/e6nb0ert440zbul6i84gl28u060b027t)**. This program works on 64-bit windows.

### Prerequisite
+ **Boost** (1.66 or later): https://sourceforge.net/projects/boost/files/boost-binaries/
+ **CGAL** (4.13 or later): https://www.cgal.org/releases.html
+ **MOSEK 8**: https://www.mosek.com/downloads/, MOSEK license is required.
+ **PolyCut**: http://www.cs.ubc.ca/labs/imager/tr/2018/HexDemo/, download and extract it to a folder, be aware that it has one-month license.
+ **HexEx**: https://www.graphics.rwth-aachen.de/media/resource_files/HexEx_Windows_1_01.zip, download and extract it to a folder.

The libraries and softwares above need to be installed **manually**, please set BOOST_ROOT CGAL_DIR in your system path if Boost and CGAL installers do not set them sucessfully. 

Our source code has the following build-in dependencies:
+ **Eigen3**: http://eigen.tuxfamily.org/index.php?title=Main_Page#Download
+ **OpenMesh**: https://www.graphics.rwth-aachen.de/software/openmesh/
+ **OpenVolumeMesh**: https://www.graphics.rwth-aachen.de/software/openvolumemesh/, we modified the source code.
+ **gflags**: https://github.com/gflags/gflags.git
+ **rapidxml**: http://rapidxml.sourceforge.net/
+ **SUS**: https://www.dca.iusiani.ulpgc.es/SUScode/, we modified the source code.
+ **ANN**: http://www.cs.umd.edu/~mount/ANN/
+ **TetGen**: http://wias-berlin.de/software/index.jsp?id=TetGen&lang=1, we modified the source code slightly.
+ **SuiteSparse**: http://faculty.cse.tamu.edu/davis/suitesparse.html
+ **VTK (8.2)**: https://vtk.org/download/,  One of utility programs -- vtudecode requires VTK, we have provided the binary file and it is not necessary to recompile it.
+ **CXXOPTS**: https://github.com/jarro2783/cxxopts
+ **LAPACK**: https://www.netlib.org/lapack/

### Data Preparation
Our program takes a tetrahedral mesh (**\*.vtk**) and its feature edges (**\*.fea**) as input.  The feature file starts with the number of feature segments **n**, followed by **n** lines, where each line contains the two vertex indices of a feature edge segment. You can generate tetrahedral meshes from surface meshes using **[TetGen](http://wias-berlin.de/software/index.jsp?id=TetGen&lang=1)** or download our preprocessed data from **[here](https://drive.google.com/open?id=1KUysv4rmUIX1F5nnHKHaek5UboeCclE8)**

### Installation & Usage
First, clone this repository:
```
git clone https://github.com/msraig/CE-PolyCube.git
cd CE-PolyCube
```
We provide both the source code and the compiled executable files(in **Bin** folder). To use the executable files directly, please first set DATA_ROOT_PATH (path contains model folders), HEXEX_PATH (path containing **hexex.exe**) and POLYCUT_PATH (path containing **polycut.exe**) in **.\Script\gen_hex.bat** with no trailing slash. Then you can process a specific model, e.g sculpt by running:
```
cd Scripts
gen_hex.bat sculpt
```
If everything goes well, there will be three new files in DATA_ROOT_PATH/sculpt folder:
+ **sculpt_deform_polycube.vtk**: the initial PolyCube mesh.
+ **sculpt_cut_flattening.vtk**: the CE-PolyCube mesh.
+ **sculpt_hex_opt.vtk**: the final optimized all-hex mesh.

Or if you want to recompile the project, please go to the root directory, then run:
```
mkdir Build
cd Build
cmake ..
```
Build the solution in .\Build folder, then you can find the generated executable files in .\Bin folder. Then modify and run **.\Script\gen_hex.bat** as previously mentioned.

We also provide the processed data **[here](https://drive.google.com/file/d/1g1RwWSkPRhl4HpcstE5hJALF3Zpq61pQ/view)**.

### Update
**2022.5.11** If you want to apply our algorithm on your own model, please firstly normalize the model so that it is within a bounding box of [0, 10]^3, and prepare the feature file as described above. If our algorithm fails to produce a valid PolyCube or hex mesh, here are some tips on parameter tuning:
+ If the initial PolyCube '*_deform_polycube.vtk' contains flipped/flattened tets, you can try to generate the polycube using PolyCut directly, with  '**.\Script\gen_hex_nodeform.bat**'.
+ For paramters in Line 32 of 'gen_hex.bat': 'bs' means how many edges are cut simultaneously, increasing 'bs' can sometimes improve hex quality; 'mr' means the mininum cut depth ratio, for models like 'fandisk', 'mr' needs to be set as 0.3 to ensure all non-feature edges are cut. However, if 'mr' is too small, the TetGen module might break down. 'sr' controls the size of hex element, larger value corresponds to larger hex. Decreasing 'sr' can sometimes help to generate a valid cut-enhanced PolyCube.

### Citation
```
@article{Guo2020Cut,
  title={Cut-enhanced PolyCube-Maps for Feature-aware All-Hex Meshing},
  author={Guo, Hao-Xiang and Liu, Xiaohan and Yan, Dong-Ming and Yang, Liu},
  journal={ACM Transactions on Graphics (TOG)},
  volume={39},
  number={4},
  pages={106:1--106:14},
  year={2020},
  publisher={ACM New York, NY, USA}
}
```
Please contact us (Haoxiang Guo ghx17@mails.tsinghua.edu.cn, Xiaohan Liu xh.liu.tech@gmail.com, Yang Liu yangliu@microsoft.com) or file an issue 
if you meet problems in using our code. 

