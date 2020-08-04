#pragma once
//Created by Haoxiang Guo 2019.9.20 version 1
//ghx17@mails.tsinghua.edu.cn
//The only exterior dependency are Eigen
//remove dependency on Eigen, use CVec instead
//-----------------------------------------------
//             7_________6
//            /|        /|
//           4_|_______5 |
//           | |       | |
//           | 3_______|_2
//           |/        |/
//           0_________1


//#include <Eigen/Dense>
#include <vector>
#include "GlobalType.h"
#include "SmallVec.h"
//using namespace Eigen;
using std::vector;

namespace HexStructure
{
	//template <class T, class Real>
	struct Vertex
	{
		Int id, mvid, svid;
		ig::CVec<Real, 3> v;
		
		vector<Int> n_vs; // neighboring vs
		vector<Int> n_es;
		vector<Int> n_fs;
		vector<Int> n_cs;
		
		bool boundary;
		
	};

	//template <class T, class Real>
	struct Edge
	{
		Int id;
		vector<Int> vs;
		vector<Int> n_fs;
		vector<Int> n_cs;

		bool boundary;

	};

	//template <class T, class Real>
	struct Face
	{
		Int id;

		vector<Int> vs;
		vector<Int> es;
		vector<Int> n_cs;
		
		bool boundary;

	};

	//template <class T, class Real>
	struct Cell
	{
		Int id;
		vector<Int> vs;
		vector<Int> es;
		vector<Int> fs;
		
		Real quality; // hex quality of a specific type
		
		bool boundary;
	};

	struct Mesh
	{
		short type; ////temporatily only support hex mesh
		//MatrixXd V; //coordinates 3 * n_vert
		vector<Vertex> Vs; //Order consistent with vtk file format
		vector<Edge> Es;
		vector<Face> Fs;
		vector<Cell> Cs;
	};


	//metamesh and singularity remains to be defined


	struct MetaMesh_V
	{
		Int id, hid, svid; //metamesh id, hex mesh id, singularity id
		//Int what_type;
		vector<Int>  n_fvs;
		vector<Int>  n_fes;
		vector<Int>  n_ffs;
		vector<Int>  n_fhs;
		bool boundary;
	};
	struct MetaMesh_E
	{
		Int id;
		bool singular = false;
		std::vector<Int> vs;
		bool boundary;
		vector<Int>  vs_onmesh;//vertex of mesh
		vector<Int>  es_onmesh;//edge of mesh
		vector<Int>  n_fes;
		vector<Int>  n_ffs;
		vector<Int>  n_fhs;
	};
	struct MetaMesh_F
	{
		Int id;
		bool boundary;
		//Int F_location;

		vector<Int> vs; //idx in the meta mesh space
		vector<Int> es;
		vector<Int>  fvs_grid;
		vector<Int>  ffs_grid;
		vector<Int>  n_ffs;
		vector<Int>  n_fhs;
		Int Color_ID;
	};
	struct MetaMesh_H
	{
		Int id;

		std::vector<Int> vs;
		std::vector<Int> es;
		std::vector<Int> fs;
		vector<vector<vector<Int> >> vs_grid;
		vector<Int>  fs_grid;
		vector<Int>  hs_grid;
		vector<Int>  n_fhs;//neighboring cube	
		Int Color_ID;
	};

	struct MetaMesh
	{
		vector<MetaMesh_V> FVs;
		vector<MetaMesh_E> FEs;
		vector<MetaMesh_F> FFs;
		vector<MetaMesh_H> FHs;
	};

	struct Singularity_V
	{
		Int id, oid; // oid means id in the hex mesh
		
		vector<Int> n_svs;
		vector<Int> n_ses;
		
		bool fake;

		bool boundary;

	};

	struct Singularity_E
	{
		Int id;
		vector<Int> vs;
		vector<Int> es_onmesh;
		vector<Int> vs_onmesh;
		
		bool boundary;
		vector<Int> n_ses;
		bool circle;
	};


	struct Singularity
	{
		vector<Singularity_V> SVs;
		vector<Singularity_E> SEs;
	};

	struct HexFrame
	{
		//contains the Hex mesh with its related singularity structure and meta mesh
		Mesh mesh;
		MetaMesh meta_mesh;
		Singularity singularity;
		
	};
}


