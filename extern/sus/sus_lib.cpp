/* $Id: sus.cpp 35 2010-09-24 10:49:40Z erodriguez $ */

/* -------------------------------------------------------------
SUS code
Simultaneous Untangling and Smoothing
Programa de desenredado y suavizado de mallas 3D

sus.cpp

Copyright (C) 2009, 2010 ULPGC
Universidad de Las Palmas de Gran Canaria
Institute of Intelligent Systems and
Numerical Applications in Engineering - SIANI
Jose Maria Escobar Sanchez
Eduardo Rodriguez Barrera
Rafael Montenegro Armas
Gustavo Montero Garcia
Jose Maria Gonzalez Yuste

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------*/

#pragma once
#include "sus_lib.h"
#include "Uncmin.h"
#include "vec.h"	// Simple C++ Numeric Toolkit (SCPPNT) vector class
#include "cmat.h"	// Simple C++ Numeric Toolkit (SCPPNT) matrix class
#include "sus.h"

#include <iostream>     // for cout
#include <cstdio>	// for printf
#include <string.h>     // para strerror
#include <assert.h>     
#include <cfloat>
#include <errno.h>
#include "getopt.h"
#include <vector>
#include <Eigen/Dense>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/Status.hh>

//#include <list>
//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/AABB_tree.h>
//#include <CGAL/AABB_traits.h>
//#include <CGAL/AABB_triangle_primitive.h>
//#include <CGAL/Polyhedron_incremental_builder_3.h>
//#include <CGAL/Polyhedron_3.h>
//#include <CGAL/Polyhedron_items_with_id_3.h>
//#include <CGAL/AABB_face_graph_triangle_primitive.h>
//
//typedef CGAL::Simple_cartesian<double> Ker;
//typedef Ker::FT FT;
//typedef Ker::Point_3 Point3;
//typedef Ker::Point_2 Point2;
//typedef std::vector<FT> Scalar_vector;
////typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Ker> Triangle_coordinates;
//
//typedef Ker::Segment_3 Segment;
//typedef CGAL::Polyhedron_3<Ker, CGAL::Polyhedron_items_with_id_3> Polyhedron;
//typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
//typedef CGAL::AABB_traits<Ker, Primitive> Traits;
//typedef CGAL::AABB_tree<Traits> Tree;
//typedef Tree::Point_and_primitive_id Point_and_primitive_id;
//typedef Polyhedron::HalfedgeDS             HalfedgeDS;

template <class HDS>
class Build_triangle : public CGAL::Modifier_base<HDS> {
public:
	Build_triangle() {}
	void ConstructList(std::vector<Point3> PointList, std::vector<std::vector<int>> FaceList)
	{
		point_list = PointList;
		face_list = FaceList;
	}
	void operator()(HDS& hds) {
		// Postcondition: hds is a valid polyhedral surface.
		CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
		B.begin_surface(point_list.size(), face_list.size());
		typedef typename HDS::Vertex   Vertex;
		typedef typename Vertex::Point Point;
		for (int i = 0; i < point_list.size(); i++)
		{
			B.add_vertex(point_list[i]);
		}

		for (int i = 0; i < face_list.size(); i++)
		{
			B.begin_facet();
			for (int j = 0; j < face_list[i].size(); j++)
			{
				B.add_vertex_to_facet(face_list[i][j]);
			}
			B.end_facet();
		}
		B.end_surface();
	}
	std::vector<Point3> point_list;
	std::vector<std::vector<int>> face_list;
};

using Eigen::Matrix3d;
using Eigen::Vector3d;

using namespace std;

#define MASCARA_X (1)
#define MASCARA_Y (2)
#define MASCARA_Z (4)
#define MAX_CARAS_VISTAS 100000
#define BUFFER_SIZE 2048




char const *version = "1.02"; // Program version

double W0_inversa[3][3] = { {1.0,  -1.0 / rd3 ,  -1.0 / rd6} ,
				{0.0,   2.0 / rd3 ,  -1.0 / rd6} ,
				{0.0,   0.0     ,  rd3 / rd2} };

double DosTercios = 2.0 / 3.0;
double Norma = 1.0;
int nodo;         // nodo actual en la minimizacion

T_Malla LaMalla;   // Declara la malla
T_Malla *m;        // 指向LaMalla的指针
struct mis_parametros params;
double peor_sigma = DBL_MAX;
double deltaglobal;

int MaxNumIter = 100;          // Max. number of iterations allowed
int NumIterSuavizado = 20;     // Number of smoothing iterations 

// -------- Some important global variables (see README) ---------
double asimptotic_epsilon = 1.0e-10; // Asymptotic value of 
									 // h(sigma) when sigma tends to -inf
double EpsSafetyFactor = 1.0e6;      // "Safety" factor for epsilon
double meps = 2.220446e-16;         // Machine epsilon
double EpsilonEfectivo = meps * EpsSafetyFactor; // Effective epsilon

std::vector<double> s, sigma;
std::vector<double> boundary_neighbor_avg_coord;

//Tree *search_tree = NULL;
//Polyhedron polyhedron;

//std::vector<Vector3d> proj_normals;


Tree* project_points_prepare(const std::vector<double> &input_pts, const std::vector<int> &input_faces,  std::vector<Vector3d>& normals)
{
	//proj_pts_on faces;
	std::cout << "project points on faces" << std::endl;

	Polyhedron polyhedron;
	std::vector<Point3> PointList;
	std::vector<std::vector<int>> FaceList;
	std::vector<int> face;

	int input_pts_size = (int)input_pts.size() / 3;
	int input_face_size = (int)input_faces.size() / 3;
	for (size_t i = 0; i < input_pts_size; i++)
	{
		//PointList.push_back(Point3(source_points[i][0], source_points[i][1], source_points[i][2]));
		PointList.push_back(Point3(input_pts[3 * i], input_pts[3 * i + 1], input_pts[3 * i + 2]));
	}
	//compute normals
	normals.clear();
	for (size_t i = 0; i < input_face_size; i++)
	{
		face.clear();
		/*face.push_back(source_faces[i][0]);
		face.push_back(source_faces[i][1]);
		face.push_back(source_faces[i][2]);*/
		face.push_back(input_faces[3 * i + 0]);
		face.push_back(input_faces[3 * i + 1]);
		face.push_back(input_faces[3 * i + 2]);

		Vector3d v1, v2, tmp_normal;
		v1[0] = input_pts[3 * face[1] + 0] - input_pts[3 * face[0] + 0];
		v1[1] = input_pts[3 * face[1] + 1] - input_pts[3 * face[0] + 1];
		v1[2] = input_pts[3 * face[1] + 2] - input_pts[3 * face[0] + 2];
					
		v2[0] = input_pts[3 * face[2] + 0] - input_pts[3 * face[0] + 0];
		v2[1] = input_pts[3 * face[2] + 1] - input_pts[3 * face[0] + 1];
		v2[2] = input_pts[3 * face[2] + 2] - input_pts[3 * face[0] + 2];

		tmp_normal = v1.cross(v2);

		tmp_normal = tmp_normal / tmp_normal.norm();
		/*normals.push_back(tmp_normal[0]);
		normals.push_back(tmp_normal[1]);
		normals.push_back(tmp_normal[2]);*/
		normals.push_back(tmp_normal);


		

		FaceList.push_back(face);
	}

	/*face.push_back(0);
	face.push_back(1);
	face.push_back(2);
	FaceList.push_back(face);
	face.clear();
	face.push_back(1);
	face.push_back(3);
	face.push_back(2);
	FaceList.push_back(face);*/

	Build_triangle<HalfedgeDS> triangles;
	triangles.ConstructList(PointList, FaceList);
	polyhedron.delegate(triangles);
	int count = 0;
	for (Polyhedron::Face_iterator i = polyhedron.facets_begin(); i != polyhedron.facets_end(); ++i)
		i->id() = count++;
	count = 0;
	for (Polyhedron::Vertex_iterator viter = polyhedron.vertices_begin(); viter != polyhedron.vertices_end(); ++viter)
		viter->id() = count++;

	//CGAL_assertion(polyhedron.is_triangle(polyhedron.halfedges_begin()));
	// constructs AABB tree and computes internal KD-tree 
	// data structure to accelerate distance queries
	Tree *tree = new Tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	//Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree->accelerate_distance_queries();


	//Tree temp_(faces(polyhedron).first, faces(polyhedron).second, polyhedron);


	Point3 query(input_pts[0], input_pts[1], input_pts[2]);
	Point_and_primitive_id pp = tree->closest_point_and_primitive(query);
	int face_id = (int)pp.second->id();

	Point3 closest = pp.first;
	
	return tree;
	
	
	
	
	// query point
	/*face_idx_list.clear();
	barycentric_coord.clear();*/
	//int proj_pts_size = proj_pts.size() / 3;
	//for (size_t i = 0; i < proj_pts_size; i++)
	//{

	//	//Point3 query(target_points[i][0], target_points[i][1], target_points[i][2]);
	//	Point3 query(proj_pts[3 * i], proj_pts[3 * i + 1], proj_pts[3 * i + 2]);
	//	// computes squared distance from query
	//	/*FT sqd = tree.squared_distance(query);
	//	std::cout << "squared distance: " << sqd << std::endl;*/
	//	// computes closest point
	//	Point3 closest = tree.closest_point(query);
	//	proj_pts[3 * i] = closest[0];
	//	proj_pts[3 * i + 1] = closest[1];
	//	proj_pts[3 * i + 2] = closest[2];

	//}

}




int
ExtraeNumero(char **linea)
{
	int i;
	int b;
	char c;

	i = (int)strcspn(*linea, "0123456789");
	(*linea) += i;
	i = (int)strspn(*linea, "0123456789");
	if (i == 0) return(-1);
	c = (*linea)[i];
	(*linea)[i] = '\0';
	b = atoi(*linea);
	(*linea)[i] = c;
	(*linea) += i;
	return(b);
}/* ExtraeNumero */

int
LeeLinea(FILE *f, char *buffer)
{
	/* Se salta las lineas vacias o con comentarios */
	/* Los comentarios tienen el caracter # en la columna 0 */
	if (feof(f)) return(0);
	do
	{
		fgets(buffer, BUFFER_SIZE, f);
		buffer[strlen(buffer) - 1] = '\0';
	} while (((buffer[0] == '#') || (buffer[0] == '\0')) && (!feof(f)));
	if ((buffer[0] == '\n') && feof(f)) return(0);
	if (feof(f))
		return(0);
	else
		return(1);
}/* LeeLinea */


/* Calcula el numero de caras vistas desde cada nodo */
void
CuentaCarasVistas(T_Malla *Malla)
{
	int i, j;

	for (i = 0; i < Malla->Num_nodos; i++) Malla->nodo[i].NCarasVistas = 0;

	for (i = 0; i < Malla->Num_tetra; i++)
		for (j = 0; j < 4; j++)
		{
			int k = Malla->tetra[i].nodo[j];
			Malla->nodo[k].NCarasVistas++;
		}
} /* CuentaCarasVistas */

/* Calcula las caras vistas desde cada nodo */
void
CalculaCarasVistas(T_Malla *Malla)
{
	int i, j;
	int *cont;
	int orden[4][3] = { {1, 2, 3}, {0, 3, 2},
				{0, 1, 3}, {0, 2, 1} };

	cont = (int *)malloc(Malla->Num_nodos * sizeof(int));
	if (cont == NULL)
	{
		fprintf(stderr, "CalculaCarasVistas: No se puede reservar memoria\n");
		abort();
	}

	for (i = 0; i < Malla->Num_nodos; i++) cont[i] = 0;

	for (i = 0; i < Malla->Num_tetra; i++)
	{
		for (j = 0; j < 4; j++)
		{
			int n;
			int k = Malla->tetra[i].nodo[j];
			int idx = cont[k];

			for (n = 0; n < 3; n++)
			{
				int nodo = Malla->tetra[i].nodo[orden[j][n]];
				Malla->nodo[k].CarasVistas[n][idx] = nodo;
			}
			cont[k]++;
		}
	}

	free(cont);
} /* CalculaCarasVistas */

void
CalculaLasConexiones(T_Malla *Malla)
{
	int i;
	int *m;

	/* Calcula el numero de caras vista desde cada nodo (las conexiones) */
	CuentaCarasVistas(Malla);
	/* Reserva memoria para las caras vistas */
	for (i = 0; i < Malla->Num_nodos; i++)
	{
		if ((m = (int *)malloc(3 * sizeof(int) * Malla->nodo[i].NCarasVistas))
			== NULL) {
			fprintf(stderr,
				"LeerMalla: no puedo reservar memoria para las conexiones: %s\n",
				strerror(errno));
			abort();
		}
		Malla->nodo[i].CarasVistas[0] = m;
		Malla->nodo[i].CarasVistas[1] = m + Malla->nodo[i].NCarasVistas;
		Malla->nodo[i].CarasVistas[2] = m + 2 * Malla->nodo[i].NCarasVistas;
	}
	/* Calcula las caras vistas desde cada nodo */
	CalculaCarasVistas(Malla);

} /* CalculaLasConexiones */

/* ----------------------------------------------------------------- */
/* 12/01/2007                                                        */
/*                                                                   */
/* funcion LeerMalla                                                 */
/* Lee los ficheros de coordenadas y elementos de una malla y        */
/* calcula las conexiones, para evitar tener que tener un fichero    */
/* de conexiones.                                                    */
/*                                                                   */
/* Entradas:                                                         */
/*   f_coor: nombre del fichero de coordenadas de los nodos          */
/*   f_elem: nombre del fichero de nodos que componen cada tetraedro */
/* Salidas:                                                          */
/*   Malla: Rellena la estructura con los datos leidos               */
/*   Devuelve el numero de nodos leidos                              */
/* ----------------------------------------------------------------- */

int InitMalla(const std::vector<double> &nodes_coord, const std::vector<int> &ref_nodes, const std::vector<int> &elems, T_Malla *Malla)
{
	if (nodes_coord.size() != 3 * ref_nodes.size())
	{
		std::cout << "Initialize Malla error" << std::endl;
		return 0;
	}

	//input node information
	int N_nodes = (int)ref_nodes.size();
	Malla->Num_nodos = N_nodes;
	if ((Malla->nodo = (T_Nodo *)malloc(sizeof(T_Nodo) * Malla->Num_nodos))
		== NULL) {
		fprintf(stderr,
			"LeerMalla: no puedo reservar memoria para los nodos: %s\n",
			strerror(errno));
		exit(1);
	}
	for (size_t i = 0; i < N_nodes; i++)
	{
		Malla->nodo[i].x = nodes_coord[3 * i];
		Malla->nodo[i].y = nodes_coord[3 * i + 1];
		Malla->nodo[i].z = nodes_coord[3 * i + 2];
		Malla->nodo[i].nr = ref_nodes[i];

	}


	//input tet information
	Malla->Num_tetra = (int)elems.size() / 4;
	if ((Malla->tetra = (T_Tetra *)malloc(sizeof(T_Tetra) * Malla->Num_tetra))
		== NULL) {
		fprintf(stderr,
			"LeerMalla: no puedo reservar memoria para los tetraedros: %s\n",
			strerror(errno));
		exit(1);
	}

	for (size_t i = 0; i < Malla->Num_tetra; i++)
	{
		//index start from 0
		Malla->tetra[i].nodo[0] = elems[4 * i];
		Malla->tetra[i].nodo[1] = elems[4 * i + 1];
		Malla->tetra[i].nodo[2] = elems[4 * i + 2];
		Malla->tetra[i].nodo[3] = elems[4 * i + 3];
	}

	CalculaLasConexiones(Malla);
	
	
	return (Malla->Num_nodos);

	
}



int InitMalla(const std::vector<double> &nodes_deform, const std::vector<int> &ref_nodes, const std::vector<unsigned int> &elems, const std::vector<double> &nodes_ori, const std::vector<int> &normal_tag,  T_Malla *Malla)
{
	if (nodes_deform.size() != 3 * ref_nodes.size())
	{
		std::cout << "Initialize Malla error" << std::endl;
		return 0;
	}

	//input node information
	int N_nodes = (int)ref_nodes.size();
	Malla->Num_nodos = N_nodes;
	if ((Malla->nodo = (T_Nodo *)malloc(sizeof(T_Nodo) * Malla->Num_nodos))
		== NULL) {
		fprintf(stderr,
			"LeerMalla: no puedo reservar memoria para los nodos: %s\n",
			strerror(errno));
		exit(1);
	}
	for (size_t i = 0; i < N_nodes; i++)
	{
		Malla->nodo[i].x = nodes_deform[3 * i];
		Malla->nodo[i].y = nodes_deform[3 * i + 1];
		Malla->nodo[i].z = nodes_deform[3 * i + 2];

		Malla->nodo[i].o_x = nodes_ori[3 * i];
		Malla->nodo[i].o_y = nodes_ori[3 * i + 1];
		Malla->nodo[i].o_z = nodes_ori[3 * i + 2];

		Malla->nodo[i].normal_x = normal_tag[3 * i];
		Malla->nodo[i].normal_y = normal_tag[3 * i + 1];
		Malla->nodo[i].normal_z = normal_tag[3 * i + 2];

		Malla->nodo[i].nr = ref_nodes[i];

	}


	//input tet information
	Malla->Num_tetra = (int)elems.size() / 4;
	if ((Malla->tetra = (T_Tetra *)malloc(sizeof(T_Tetra) * Malla->Num_tetra))
		== NULL) {
		fprintf(stderr,
			"LeerMalla: no puedo reservar memoria para los tetraedros: %s\n",
			strerror(errno));
		exit(1);
	}

	for (size_t i = 0; i < Malla->Num_tetra; i++)
	{
		//index start from 0
		Malla->tetra[i].nodo[0] = elems[4 * i];
		Malla->tetra[i].nodo[1] = elems[4 * i + 1];
		Malla->tetra[i].nodo[2] = elems[4 * i + 2];
		Malla->tetra[i].nodo[3] = elems[4 * i + 3];
		T_Tetra *t = &Malla->tetra[i];
		//std::cout << "calculating: " << i << std::endl;
		Matrix3d s_o;
		Matrix3d s_o_I;
		for (int j = 1; j < 4; j++) {
			s_o(0, j - 1) = Malla->nodo[t->nodo[j]].o_x - Malla->nodo[t->nodo[0]].o_x;
			s_o(1, j - 1) = Malla->nodo[t->nodo[j]].o_y - Malla->nodo[t->nodo[0]].o_y;
			s_o(2, j - 1) = Malla->nodo[t->nodo[j]].o_z - Malla->nodo[t->nodo[0]].o_z;
		}
		//std::cout << "i: " << i << " det: " << s_o.determinant() << std::endl;
		s_o_I = s_o.inverse();
		for (size_t j = 0; j < 3; j++)
		{
			for (size_t k = 0; k < 3; k++)
			{
				t->IS[j][k] = s_o_I(j, k);
				//t->IS[j][k] = s_o(j, k);
			}
		}


	}

	


	
	

	


	CalculaLasConexiones(Malla);





	return (Malla->Num_nodos);


}




int
LeerMalla(char *f_coor, char *f_elem, T_Malla *Malla)
{
	FILE *f1, *f2;
	char *buffer;
	int i, foo;

	/* Tengo que usar memoria dinamica para el buffer porque si no */
	/* no funciona bien ExtreNumero */

	if ((buffer = (char *)malloc(BUFFER_SIZE * sizeof(char))) == NULL) {
		fprintf(stderr,
			"LeerMalla: no puedo reservar memoria para el buffer: %s\n",
			strerror(errno));
		exit(1);
	}

	if ((f1 = fopen(f_coor, "r")) == NULL) {
		fprintf(stderr, "LeerMalla: error abriendo fichero %s: %s\n",
			f_coor, strerror(errno));
		exit(1);
	}
	if ((f2 = fopen(f_elem, "r")) == NULL) {
		fprintf(stderr, "LeerMalla: error abriendo fichero %s: %s\n",
			f_elem, strerror(errno));
		exit(1);
	}

	/* Reserva memoria para las coordenadas de los nodos y las */
	/* condiciones de contorno */
	LeeLinea(f1, buffer);
	foo = sscanf(buffer, "%d", &Malla->Num_nodos);
	assert(foo);
	assert(foo != EOF);
	if ((Malla->nodo = (T_Nodo *)malloc(sizeof(T_Nodo) * Malla->Num_nodos))
		== NULL) {
		fprintf(stderr,
			"LeerMalla: no puedo reservar memoria para los nodos: %s\n",
			strerror(errno));
		exit(1);
	}
	/* Lectura de las coordenadas de los nodos y sus */
	/* numeros de referencia */
	i = 0;
	while ((LeeLinea(f1, buffer) != 0) && (i < Malla->Num_nodos)) {
		sscanf(buffer, "%lf %lf %lf %d", &Malla->nodo[i].x,
			&Malla->nodo[i].y, &Malla->nodo[i].z,
			&Malla->nodo[i].nr
		);
		i++;
	}
	fclose(f1);

	if (i < Malla->Num_nodos) {
		fprintf(stderr, "LeerMalla: Problema al leer el fichero de nodos:\n");
		fprintf(stderr, "debo leer %d nodos y solo he podido leer %d\n",
			Malla->Num_nodos, i - 1);
		exit(EXIT_FAILURE);
	}

	/* Reserva memoria para los tetraedros */
	LeeLinea(f2, buffer);
	foo = sscanf(buffer, "%d", &Malla->Num_tetra);
	assert(foo);
	assert(foo != EOF);
	if ((Malla->tetra = (T_Tetra *)malloc(sizeof(T_Tetra) * Malla->Num_tetra))
		== NULL) {
		fprintf(stderr,
			"LeerMalla: no puedo reservar memoria para los tetraedros: %s\n",
			strerror(errno));
		exit(1);
	}

	/* Lectura de los nodos que componen cada tetraedro */
	i = 0;
	while ((LeeLinea(f2, buffer) != 0) && (i < Malla->Num_tetra)) {
		sscanf(buffer, "%d %d %d %d %d", &foo, &Malla->tetra[i].nodo[0],
			&Malla->tetra[i].nodo[1], &Malla->tetra[i].nodo[2],
			&Malla->tetra[i].nodo[3]);
		/* OJO : resto 1 porque la malla de escobar empieza en el nodo 1 (FORTRAN)
		   mientras que la mia empieza en 0 (C) */
		Malla->tetra[i].nodo[0]--; Malla->tetra[i].nodo[1]--;
		Malla->tetra[i].nodo[2]--; Malla->tetra[i].nodo[3]--;
		i++;
	}
	fclose(f2);
	if (i < Malla->Num_tetra) {
		fprintf(stderr, "LeerMalla: Problema al leer el fichero %s:\n", f_elem);
		fprintf(stderr, "debo leer %d tetraedros y solo he podido leer %d\n",
			Malla->Num_tetra, i - 1);
		exit(EXIT_FAILURE);
	}

	CalculaLasConexiones(Malla);

	free(buffer);
	return(Malla->Num_nodos);
}/* LeerMalla */



/*
  Conviene que la malla este centrada en el origen y escalada, lo que
  viene siendo una normalizacion, vamos
*/
void
NormalizarLaMalla(T_Malla *Malla)
{
	int i;
	double cx, cy, cz;     // Centro de gravedad de la malla
	double max, max_radio;

	cx = cy = cz = 0.0;
	// Calculo del centro de gravedad
	for (i = 0; i < Malla->Num_nodos; i++)
	{
		cx += Malla->nodo[i].x;
		cy += Malla->nodo[i].y;
		cz += Malla->nodo[i].z;
	}
	cx = cx / Malla->Num_nodos;
	cy = cy / Malla->Num_nodos;
	cz = cz / Malla->Num_nodos;

	// Traslado la malla, restando a cada coordenada la componente
	// correspondiente del centro de gravedad
	for (i = 0; i < Malla->Num_nodos; i++)
	{
		Malla->nodo[i].x = Malla->nodo[i].x - cx;
		Malla->nodo[i].y = Malla->nodo[i].y - cy;
		Malla->nodo[i].z = Malla->nodo[i].z - cz;
	}

	max_radio = 0.0;
	// Calculo del factor de escala (esfera envolvente)
	for (i = 0; i < Malla->Num_nodos; i++)
	{
		max = Malla->nodo[i].x * Malla->nodo[i].x;
		max += Malla->nodo[i].y * Malla->nodo[i].y;
		max += Malla->nodo[i].z * Malla->nodo[i].z;
		if (max > max_radio)
			max_radio = max;
	}

	// Escalado de la malla
	max_radio = sqrt(max_radio);
	max_radio = max_radio / 10.0;
	for (i = 0; i < Malla->Num_nodos; i++)
	{
		Malla->nodo[i].x = Malla->nodo[i].x / max_radio;
		Malla->nodo[i].y = Malla->nodo[i].y / max_radio;
		Malla->nodo[i].z = Malla->nodo[i].z / max_radio;
	}

	// Guarda la informacion en la malla para poder des-normalizar tras
	// el suavizado
	Malla->cx = cx;    Malla->cy = cy;    Malla->cz = cz;
	Malla->max_radio = max_radio;
}


/*
  Devuelve una malla normalizada a su forma y lugar original
*/
void
DesnormalizarLaMalla(T_Malla *Malla)
{
	int i;
	double cx, cy, cz;     // Centro de gravedad de la malla
	double max, max_radio;

	// Recupera la informacion necesaria para des-normalizar
	cx = Malla->cx;    cy = Malla->cy;    cz = Malla->cz;
	max_radio = Malla->max_radio;

	// des-escalado de la malla
	for (i = 0; i < Malla->Num_nodos; i++)
	{
		Malla->nodo[i].x = Malla->nodo[i].x * max_radio;
		Malla->nodo[i].y = Malla->nodo[i].y * max_radio;
		Malla->nodo[i].z = Malla->nodo[i].z * max_radio;
	}

	// des-traslado la malla, sumando a cada coordenada la componente
	// correspondiente del centro de gravedad
	for (i = 0; i < Malla->Num_nodos; i++)
	{
		Malla->nodo[i].x = Malla->nodo[i].x + cx;
		Malla->nodo[i].y = Malla->nodo[i].y + cy;
		Malla->nodo[i].z = Malla->nodo[i].z + cz;
	}
}




/* ----------------------------------------------------------------- */
/* funcion LiberarMalla                                              */
/*                                                                   */
/* Libera la memoria usada por los elementos que componen una malla  */
/*                                                                   */
/* Entrada:                                                          */
/*  Malla: puntero a la malla que hay que liberar                    */
/* ----------------------------------------------------------------- */
void
LiberarMalla(T_Malla *Malla)
{
	free(Malla->nodo);
	free(Malla->tetra);
}/* LiberarMalla */



// -------------------------------------------------------------
// -------------------------------------------------------------
//           DEFINICION DE LAS FUNCION DE KNUPP
// -------------------------------------------------------------
// -------------------------------------------------------------

/* ---------------------------------------------------------------------------
Det
Calcula el determinante de una matriz de 3x3 con pivotacion total
--------------------------------------------------------------------------- */
double
//Det(double m[3][3])
Det(double *m)
{
	int ncont = 0;
	int ir, ifil, icol, i, j, k;
	double amax, aux;
	double A[3][3];

	memcpy(A, m, 9 * sizeof(double));
	/* Determinacion del pivote maximo */
	for (ir = 0; ir < 2; ir++) {
		amax = fabs(A[ir][ir]);
		ifil = icol = ir;
		for (i = ir; i < 3; i++)
			for (j = ir; j < 3; j++)
				if (fabs(A[i][j]) > amax) {
					ifil = i;
					icol = j;
					amax = fabs(A[i][j]);
				}
		/* Intercambio de la fila ir por la fila ifil */
		if (ifil != ir) {
			ncont++;
			for (k = ir; k < 3; k++) {
				aux = A[ir][k];
				A[ir][k] = A[ifil][k];
				A[ifil][k] = aux;
			}
		}
		/* Intercambio de columnas */
		if (icol != ir) {
			ncont++;
			for (k = 0; k < 3; k++) {
				aux = A[k][ir];
				A[k][ir] = A[k][icol];
				A[k][icol] = aux;
			}
		}

		for (i = ir + 1; i < 3; i++)
			for (j = ir + 1; j < 3; j++)
				A[i][j] = A[i][j] - A[i][ir] * A[ir][j] / A[ir][ir];
	}

	aux = A[0][0] * A[1][1] * A[2][2];
	if (ncont == 1 || ncont == 3)
		aux = -aux;
	return(aux);
}/* Det */


// 21/11
double
h_sigma(const double d, double DeltaCuadrado, double epsilon_sigma)
{
	double h;
	h = 0.5 * (d + sqrt(((d - 2.0*epsilon_sigma) * (d - 2.0*epsilon_sigma))
		+ (4.0 * DeltaCuadrado)));

	return(h);
}


/*
  Funcion objetivo modificada
*/
double
//K_Knupp_Modif_directa(const T_Malla *m, double s[3][3], 
K_Knupp_Modif_directa(const T_Malla *m, double *s,
	const double det, const double Delta, double *nfr,
	double epsilon_sigma) {
	double nf = 0.0;
	double  d = det;
	double DeltaCuadrado = Delta * Delta;

	d = h_sigma(det, DeltaCuadrado, epsilon_sigma);
	d = cbrt(d * d);

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		{
			nf = nf + (s[i * 3 + j] * s[i * 3 + j]);
		}
	*nfr = sqrt(nf);

	double aux = nf / (3.0 * d);
	return(aux);
} /* K_Knupp_Modif_directa */

/*
Construye la matriz jacobiana S
Entradas:
   m : malla
   t : tetraedro
Salidas:
   s: matriz jacobiana
*/
void

//建立S
//ConstruyeS(const T_Malla *m, T_Tetra t, double s[3][3])
ConstruyeS(const T_Malla *m, T_Tetra &t, double *s)
{
	double aux[3][3];

	/* 构建雅可比矩阵 */
#if 0
	cout << "tetraedro = "
		<< t.nodo[0] << " "
		<< t.nodo[1] << " " << t.nodo[2] << " "
		<< t.nodo[3] << " " << endl;

	for (int i = 0; i < 4; i++) {
		cout << "nodo " << i << "="
			<< m->nodo[t.nodo[i]].x << " "
			<< m->nodo[t.nodo[i]].y << " "
			<< m->nodo[t.nodo[i]].z << " "
			<< endl;
	}
#endif

	for (int i = 1; i < 4; i++) {
		s[0 * 3 + i - 1] = m->nodo[t.nodo[i]].x - m->nodo[t.nodo[0]].x;
		s[1 * 3 + i - 1] = m->nodo[t.nodo[i]].y - m->nodo[t.nodo[0]].y;
		s[2 * 3 + i - 1] = m->nodo[t.nodo[i]].z - m->nodo[t.nodo[0]].z;

	}

#ifdef DEBUG
	/*cout << "s[]=" << s[0][0] << " " << s[0][1] << " " << s[0][2] << " " << endl;
	cout << "s[]=" << s[1][0] << " " << s[1][1] << " " << s[1][2] << " " << endl;
	cout << "s[]=" << s[2][0] << " " << s[2][1] << " " << s[2][2] << " " << endl;*/
#endif


#if STANDARD_TET
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			double p = 0.0;
			for (int k = 0; k < 3; k++)
			{
				p += s[i * 3 + k] * W0_inversa[k][j];
				//p += s[i * 3 + k] * t.IS[k][j];
				//t.IS;
				//t 非引用，无法使用IS
				//p += s[i * 3 + k] * s_o_I(k, j);
			}
			aux[i][j] = p;
		}
	}
#else
	Matrix3d s_o;
	Matrix3d s_o_I;

	for (int i = 1; i < 4; i++) {
		s_o(0, i - 1) = m->nodo[t.nodo[i]].o_x - m->nodo[t.nodo[0]].o_x;
		s_o(1, i - 1) = m->nodo[t.nodo[i]].o_y - m->nodo[t.nodo[0]].o_y;
		s_o(2, i - 1) = m->nodo[t.nodo[i]].o_z - m->nodo[t.nodo[0]].o_z;
	}

	s_o_I = s_o.inverse();


	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			double p = 0.0;
			for (int k = 0; k < 3; k++)
			{
				//p += s[i * 3 + k] * W0_inversa[k][j];
				//p += s[i * 3 + k] * t.IS[k][j];
				//t.IS;
				//t 非引用，无法使用IS
				p += s[i * 3 + k] * s_o_I(k,j);
			}
			aux[i][j] = p;
		}
	}
#endif


	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			s[i * 3 + j] = aux[i][j];
		}
	}

	//  Producto(s, W0_inversa);
} /* ConstruyeS */


double
Calcula_Delta(const void *params, double EpsilonEfectivo)
{
	struct mis_parametros *parametros = (struct mis_parametros *)params;
	double Delta = 1.0e-3;  // delta por defecto
	int j, k, nodo;
	double sigma_min = DBL_MAX;
	//double sigma[1000];
	//std::vector<double> sigma(1000);
	//double s[1000][3][3];
   // std::vector<double> s(1000 * 9);
	T_Malla *m;
	T_Tetra t;

	double dmin;

	m = parametros->malla;
	nodo = parametros->nodo;

	double sigma_medio = 0.0;
	for (j = 0; j < m->nodo[nodo].NCarasVistas; j++) {
		/* El nodo 0 del tetraedro es el nodo i considerado */
		t.nodo[0] = nodo;
		for (k = 1; k < 4; k++) t.nodo[k] = m->nodo[nodo].CarasVistas[k - 1][j];
		ConstruyeS(m, t, &s[j * 9]);
		sigma[j] = Det(&s[j * 9]);
		if (sigma[j] < sigma_min) sigma_min = sigma[j];
		sigma_medio += fabs(sigma[j]);
	}

	sigma_medio = sigma_medio / m->nodo[nodo].NCarasVistas;
#ifdef DEBUG
	cout << "SIGMA_MIN=" << sigma_min << endl;
#endif

	double radicando = EpsilonEfectivo * (EpsilonEfectivo - sigma_min);
	if (radicando <= 0.0)
		dmin = 0.0;
	else
		dmin = sqrt(radicando);

	double umbral = 1.0e-4 * sigma_medio;

	Delta = (dmin >= umbral) ? dmin : umbral;

	return(Delta);
}
/* Calcula_Delta */



/* 计算Knupp的修正F，其梯度和Hessian */
/* 使用牛顿法直接最小化 */
void
new_Fknupp_con_grad_y_hess(const double v[3], void *params,
	double *f, double gradient[3],
	double H[][3], double Delta)
{
	struct mis_parametros *parametros = (struct mis_parametros *)params;
	T_Malla *m;
	int nodo;
	T_Tetra t;
	double grad[3];
	//double s[MAX_CARAS_VISTAS][3][3];
	//std::vector<double> s(MAX_CARAS_VISTAS * 3 * 3);
	int j, k;
	double S_S[3];
	double sig_[3];
	double K;
	double K_[3];
	double K__[3][3];
	//double sigma[MAX_CARAS_VISTAS];
	//std::vector<double> sigma(MAX_CARAS_VISTAS);
	double nf = -1.0;
	double sigmacuadrado;
	double aux1, aux2, aux3;
	double acum = 0.0;
	double nf2, LaRaiz, ElFactor, ElOtroFactor;
	int i = 0;
	double sigma_min = 9.9e199;
	double hess[3][3];
	const char *nfun = "new_Fknupp_con_grad_y_hess: ";


	m = parametros->malla;
	nodo = parametros->nodo;
	parametros->delta_calculada = Delta;
	double DeltaCuadrado = Delta * Delta;

#ifdef DEBUG
	cout << nfun << endl << "\tDelta=" << Delta << endl;
#endif

	for (i = 0; i < 3; i++)
		grad[i] = hess[0][i] = hess[1][i] = hess[2][i] = 0.0;

	i = 0;
	if (parametros->mascara & MASCARA_X) {
		m->nodo[nodo].x = v[i]; i++;
	}
	else
		m->nodo[nodo].x = parametros->var[0];

	if (parametros->mascara & MASCARA_Y) {
		m->nodo[nodo].y = v[i]; i++;
	}
	else
		m->nodo[nodo].y = parametros->var[1];

	if (parametros->mascara & MASCARA_Z) {
		m->nodo[nodo].z = v[i]; i++;
	}
	else
		m->nodo[nodo].z = parametros->var[2];
#ifdef DEBUG
	cout <<
		"\tx = " << m->nodo[nodo].x <<
		"  y = " << m->nodo[nodo].y <<
		"  z = " << m->nodo[nodo].z << endl;
#endif
	acum = 0.0;

	for (j = 0; j < m->nodo[nodo].NCarasVistas; j++) {
		/* El nodo 0 del tetraedro es el nodo i considerado */
		t.nodo[0] = nodo;
		for (k = 1; k < 4; k++) t.nodo[k] = m->nodo[nodo].CarasVistas[k - 1][j];
		//ConstruyeS(m, t, s[j]);
		ConstruyeS(m, t, &s[j * 9]);
		//sigma[j] = Det(s[j]);
		sigma[j] = Det(&s[j * 9]);
		if (sigma[j] < sigma_min) sigma_min = sigma[j];
	}

	for (j = 0; j < m->nodo[nodo].NCarasVistas; j++) {
		//K = K_Knupp_Modif_directa(m, s[j], sigma[j], Delta, &nf, asimptotic_epsilon);
		K = K_Knupp_Modif_directa(m, &s[j * 9], sigma[j], Delta, &nf, asimptotic_epsilon);
		acum += K;
		sigmacuadrado = sigma[j] * sigma[j];
		nf *= nf;
		nf2 = nf; /* nf2: norma frob. al cuadrado */

		/* El nodo 0 del tetraedro es el nodo i considerado */
		t.nodo[0] = nodo;
		for (k = 1; k < 4; k++) t.nodo[k] = m->nodo[nodo].CarasVistas[k - 1][j];

		i = 0;
		if (parametros->mascara & MASCARA_X) {
			sig_[i] = sigma_x(m->nodo[t.nodo[1]].y, m->nodo[t.nodo[2]].y,
				m->nodo[t.nodo[3]].y, m->nodo[t.nodo[1]].z,
				m->nodo[t.nodo[2]].z, m->nodo[t.nodo[3]].z);
			S_S[i] = PES_gS(m->nodo[t.nodo[0]].x, m->nodo[t.nodo[1]].x,
				m->nodo[t.nodo[2]].x, m->nodo[t.nodo[3]].x);
			i++;
		}
		if (parametros->mascara & MASCARA_Y) {
			sig_[i] = sigma_y(m->nodo[t.nodo[1]].x, m->nodo[t.nodo[2]].x,
				m->nodo[t.nodo[3]].x, m->nodo[t.nodo[1]].z,
				m->nodo[t.nodo[2]].z, m->nodo[t.nodo[3]].z);
			S_S[i] = PES_gS(m->nodo[t.nodo[0]].y, m->nodo[t.nodo[1]].y,
				m->nodo[t.nodo[2]].y, m->nodo[t.nodo[3]].y);
			i++;
		}
		if (parametros->mascara & MASCARA_Z) {
			sig_[i] = sigma_z(m->nodo[t.nodo[1]].x, m->nodo[t.nodo[2]].x,
				m->nodo[t.nodo[3]].x, m->nodo[t.nodo[1]].y,
				m->nodo[t.nodo[2]].y, m->nodo[t.nodo[3]].y);
			S_S[i] = PES_gS(m->nodo[t.nodo[0]].z, m->nodo[t.nodo[1]].z,
				m->nodo[t.nodo[2]].z, m->nodo[t.nodo[3]].z);
			i++;
		}

		double aux = sigma[j] - 2.0 * asimptotic_epsilon;
		aux = aux * aux;
		LaRaiz = sqrt(aux + 4.0 * DeltaCuadrado);
		double h_sg = h_sigma(sigma[j], DeltaCuadrado, asimptotic_epsilon);
		ElFactor = (DosTercios * (1 - asimptotic_epsilon / h_sg)) / LaRaiz;
		ElOtroFactor = 2.0 / nf;

		for (k = 0; k < parametros->ndim; k++)
			K_[k] = K * ((ElOtroFactor * S_S[k]) - (ElFactor * sig_[k]));

		/* Calculo del gradiente */
		for (i = 0; i < parametros->ndim; i++) {
			aux1 = K_[i];
			if (Norma > 1.0)
				aux1 = Norma * pow(K, Norma - 1.0) * aux1;
			grad[i] += aux1;
		}

		{
			/* Hessian的计算 */
			//double aux = sigma[j] - 2*asimptotic_epsilon;
			//aux = aux*aux;
			nf *= nf; /* norma de Frobenius a la cuarta potencia */
			ElFactor = 4.0 / nf;
			//double h_sg = h_sigma(sigma[j], DeltaCuadrado, asimptotic_epsilon);
			aux1 = aux + 4.0 * DeltaCuadrado;
			aux1 = aux1 * aux1 * aux1;
			ElOtroFactor = DosTercios / sqrt(aux1);
			double termino = (1 - asimptotic_epsilon / h_sg);
			termino = termino * ((1 + asimptotic_epsilon / h_sg) * sigma[j] - 4.0 * asimptotic_epsilon);
			ElOtroFactor = ElOtroFactor * termino;

			for (i = 0; i < parametros->ndim; i++) {
				for (k = i; k < parametros->ndim; k++) { /* la matriz K__ es simetrica */
					K__[i][k] = K_[i] * K_[k] / K; /* se calcula la triangular superior */
					if (i != k)
						aux1 = 0.0;
					else
						aux1 = 3.0 / nf2;

					aux2 = ElFactor * S_S[i] * S_S[k];
					aux3 = ElOtroFactor * sig_[i] * sig_[k];
					K__[i][k] = K__[i][k] + K * (aux1 - aux2 + aux3);
				}
			}
			/* Calculo del Hessiano */
			for (i = 0; i < parametros->ndim; i++) {
				for (k = i; k < parametros->ndim; k++) {
					aux1 = K__[i][k];
					if (Norma > 1.0)
						aux1 = pow(K, Norma - 1.0) * aux1;
					aux2 = 0.0;
					if (Norma >= 2.0)
						aux2 = (Norma - 1.0) * pow(K, Norma - 2.0) * K_[i] * K_[k];
					hess[i][k] = hess[i][k] + Norma * (aux2 + aux1);
				}
			}
		}
	}

	hess[1][0] = hess[0][1];
	hess[2][0] = hess[0][2];
	hess[2][1] = hess[1][2];

	for (i = 0; i < parametros->ndim; i++)
	{
		gradient[i] = grad[i];
		for (j = 0; j < parametros->ndim; j++)
			H[i][j] = hess[i][j];
	}
	*f = acum;


	{
#ifdef DEBUG
		// Solo para DEBUG
		cout << "\tFUNCION DE KNUPP VALE " << acum << endl;

		cout << "\tGRADIENTE VALE "
			<< grad[0] << " "
			<< grad[1] << " "
			<< grad[2] << endl;

		cout << "\tHESSIANO VALE " << endl
			<< "\t" << hess[0][0] << "  " << hess[0][1] << "  " << hess[0][2] << endl
			<< "\t" << hess[1][0] << "  " << hess[1][1] << "  " << hess[1][2] << endl
			<< "\t" << hess[2][0] << "  " << hess[2][1] << "  " << hess[2][2] << endl;
#endif
	}
}/* new_Fknupp_con_grad_y_hess */


/* 计算边界Knupp的修正F，其梯度和Hessian */
/* 使用牛顿法直接最小化 */
/* 增加项： 顶点位置减去邻居点位置平均值*/

void
new_Fknupp_con_grad_y_hess_boundary(const double v[3], void *params,
	double *f, double gradient[3],
	double H[][3], double Delta)
{
	struct mis_parametros *parametros = (struct mis_parametros *)params;
	T_Malla *m;
	int nodo;
	T_Tetra t;
	double grad[3];
	//double s[MAX_CARAS_VISTAS][3][3];
	//std::vector<double> s(MAX_CARAS_VISTAS * 3 * 3);
	int j, k;
	double S_S[3];
	double sig_[3];
	double K;
	double K_[3];
	double K__[3][3];
	//double sigma[MAX_CARAS_VISTAS];
	//std::vector<double> sigma(MAX_CARAS_VISTAS);
	double nf = -1.0;
	double sigmacuadrado;
	double aux1, aux2, aux3;
	double acum = 0.0;
	double nf2, LaRaiz, ElFactor, ElOtroFactor;
	int i = 0;
	double sigma_min = 9.9e199;
	double hess[3][3];
	const char *nfun = "new_Fknupp_con_grad_y_hess: ";


	m = parametros->malla;
	nodo = parametros->nodo;
	parametros->delta_calculada = Delta;
	double DeltaCuadrado = Delta * Delta;

#ifdef DEBUG
	cout << nfun << endl << "\tDelta=" << Delta << endl;
#endif

	for (i = 0; i < 3; i++)
		grad[i] = hess[0][i] = hess[1][i] = hess[2][i] = 0.0;

	i = 0;
	if (parametros->mascara & MASCARA_X) {
		m->nodo[nodo].x = v[i]; i++;
	}
	else
		m->nodo[nodo].x = parametros->var[0];

	if (parametros->mascara & MASCARA_Y) {
		m->nodo[nodo].y = v[i]; i++;
	}
	else
		m->nodo[nodo].y = parametros->var[1];

	if (parametros->mascara & MASCARA_Z) {
		m->nodo[nodo].z = v[i]; i++;
	}
	else
		m->nodo[nodo].z = parametros->var[2];
#ifdef DEBUG
	cout <<
		"\tx = " << m->nodo[nodo].x <<
		"  y = " << m->nodo[nodo].y <<
		"  z = " << m->nodo[nodo].z << endl;
#endif
	acum = 0.0;

	for (j = 0; j < m->nodo[nodo].NCarasVistas; j++) {
		/* El nodo 0 del tetraedro es el nodo i considerado */
		t.nodo[0] = nodo;
		for (k = 1; k < 4; k++) t.nodo[k] = m->nodo[nodo].CarasVistas[k - 1][j];
		//ConstruyeS(m, t, s[j]);
		ConstruyeS(m, t, &s[j * 9]);
		//sigma[j] = Det(s[j]);
		sigma[j] = Det(&s[j * 9]);
		if (sigma[j] < sigma_min) sigma_min = sigma[j];
	}

	for (j = 0; j < m->nodo[nodo].NCarasVistas; j++) {
		//K = K_Knupp_Modif_directa(m, s[j], sigma[j], Delta, &nf, asimptotic_epsilon);
		K = K_Knupp_Modif_directa(m, &s[j * 9], sigma[j], Delta, &nf, asimptotic_epsilon);
		acum += K;
		sigmacuadrado = sigma[j] * sigma[j];
		nf *= nf;
		nf2 = nf; /* nf2: norma frob. al cuadrado */

		/* El nodo 0 del tetraedro es el nodo i considerado */
		t.nodo[0] = nodo;
		for (k = 1; k < 4; k++) t.nodo[k] = m->nodo[nodo].CarasVistas[k - 1][j];

		i = 0;
		if (parametros->mascara & MASCARA_X) {
			sig_[i] = sigma_x(m->nodo[t.nodo[1]].y, m->nodo[t.nodo[2]].y,
				m->nodo[t.nodo[3]].y, m->nodo[t.nodo[1]].z,
				m->nodo[t.nodo[2]].z, m->nodo[t.nodo[3]].z);
			S_S[i] = PES_gS(m->nodo[t.nodo[0]].x, m->nodo[t.nodo[1]].x,
				m->nodo[t.nodo[2]].x, m->nodo[t.nodo[3]].x);
			i++;
		}
		if (parametros->mascara & MASCARA_Y) {
			sig_[i] = sigma_y(m->nodo[t.nodo[1]].x, m->nodo[t.nodo[2]].x,
				m->nodo[t.nodo[3]].x, m->nodo[t.nodo[1]].z,
				m->nodo[t.nodo[2]].z, m->nodo[t.nodo[3]].z);
			S_S[i] = PES_gS(m->nodo[t.nodo[0]].y, m->nodo[t.nodo[1]].y,
				m->nodo[t.nodo[2]].y, m->nodo[t.nodo[3]].y);
			i++;
		}
		if (parametros->mascara & MASCARA_Z) {
			sig_[i] = sigma_z(m->nodo[t.nodo[1]].x, m->nodo[t.nodo[2]].x,
				m->nodo[t.nodo[3]].x, m->nodo[t.nodo[1]].y,
				m->nodo[t.nodo[2]].y, m->nodo[t.nodo[3]].y);
			S_S[i] = PES_gS(m->nodo[t.nodo[0]].z, m->nodo[t.nodo[1]].z,
				m->nodo[t.nodo[2]].z, m->nodo[t.nodo[3]].z);
			i++;
		}

		double aux = sigma[j] - 2.0 * asimptotic_epsilon;
		aux = aux * aux;
		LaRaiz = sqrt(aux + 4.0 * DeltaCuadrado);
		double h_sg = h_sigma(sigma[j], DeltaCuadrado, asimptotic_epsilon);
		ElFactor = (DosTercios * (1 - asimptotic_epsilon / h_sg)) / LaRaiz;
		ElOtroFactor = 2.0 / nf;

		for (k = 0; k < parametros->ndim; k++)
			K_[k] = K * ((ElOtroFactor * S_S[k]) - (ElFactor * sig_[k]));

		/* Calculo del gradiente */
		for (i = 0; i < parametros->ndim; i++) {
			aux1 = K_[i];
			if (Norma > 1.0)
				aux1 = Norma * pow(K, Norma - 1.0) * aux1;
			grad[i] += aux1;
		}

		{
			/* Hessian的计算 */
			//double aux = sigma[j] - 2*asimptotic_epsilon;
			//aux = aux*aux;
			nf *= nf; /* norma de Frobenius a la cuarta potencia */
			ElFactor = 4.0 / nf;
			//double h_sg = h_sigma(sigma[j], DeltaCuadrado, asimptotic_epsilon);
			aux1 = aux + 4.0 * DeltaCuadrado;
			aux1 = aux1 * aux1 * aux1;
			ElOtroFactor = DosTercios / sqrt(aux1);
			double termino = (1 - asimptotic_epsilon / h_sg);
			termino = termino * ((1 + asimptotic_epsilon / h_sg) * sigma[j] - 4.0 * asimptotic_epsilon);
			ElOtroFactor = ElOtroFactor * termino;

			for (i = 0; i < parametros->ndim; i++) {
				for (k = i; k < parametros->ndim; k++) { /* la matriz K__ es simetrica */
					K__[i][k] = K_[i] * K_[k] / K; /* se calcula la triangular superior */
					if (i != k)
						aux1 = 0.0;
					else
						aux1 = 3.0 / nf2;

					aux2 = ElFactor * S_S[i] * S_S[k];
					aux3 = ElOtroFactor * sig_[i] * sig_[k];
					K__[i][k] = K__[i][k] + K * (aux1 - aux2 + aux3);
				}
			}
			/* Calculo del Hessiano */
			for (i = 0; i < parametros->ndim; i++) {
				for (k = i; k < parametros->ndim; k++) {
					aux1 = K__[i][k];
					if (Norma > 1.0)
						aux1 = pow(K, Norma - 1.0) * aux1;
					aux2 = 0.0;
					if (Norma >= 2.0)
						aux2 = (Norma - 1.0) * pow(K, Norma - 2.0) * K_[i] * K_[k];
					hess[i][k] = hess[i][k] + Norma * (aux2 + aux1);
				}
			}
		}
	}

	hess[1][0] = hess[0][1];
	hess[2][0] = hess[0][2];
	hess[2][1] = hess[1][2];

	//test part
	/*acum = 0.0;
	for (size_t i = 0; i < 3; ++i)
	{
		grad[i] = 0.0;
		hess[i][0] = 0.0;
		hess[i][1] = 0.0;
		hess[i][2] = 0.0;
	}*/
	
	double weight = 0.0;

	weight += grad[0] * grad[0];
	weight += grad[1] * grad[1];
	weight += grad[2] * grad[2];
	weight = sqrt(weight);
	//weight = 100.0;

	

	for (i = 0; i < parametros->ndim; i++)
	{
		//modify here
		acum += weight * (v[i] - boundary_neighbor_avg_coord[3 * nodo + i]) * (v[i] - boundary_neighbor_avg_coord[3 * nodo + i]);


		gradient[i] = grad[i] + weight * 2.0 * (v[i] - boundary_neighbor_avg_coord[3 * nodo + i]);
		for (j = 0; j < parametros->ndim; j++)
			H[i][j] = hess[i][j];

		H[i][i] += weight * 2.0;
	}

	

	*f = acum;


	{
#ifdef DEBUG
		// Solo para DEBUG
		cout << "\tFUNCION DE KNUPP VALE " << acum << endl;

		cout << "\tGRADIENTE VALE "
			<< grad[0] << " "
			<< grad[1] << " "
			<< grad[2] << endl;

		cout << "\tHESSIANO VALE " << endl
			<< "\t" << hess[0][0] << "  " << hess[0][1] << "  " << hess[0][2] << endl
			<< "\t" << hess[1][0] << "  " << hess[1][1] << "  " << hess[1][2] << endl
			<< "\t" << hess[2][0] << "  " << hess[2][1] << "  " << hess[2][2] << endl;
#endif
	}
}/* new_Fknupp_con_grad_y_hess */




// -------------------------------------------------------------
// -------------------------------------------------------------
//        FIN DE DEFINICION DE LAS FUNCION DE KNUPP
// -------------------------------------------------------------
// -------------------------------------------------------------

/* ------------------------------------------------------------------- */
/* Calcula la norma de Frobenius de una matriz                         */
/* ------------------------------------------------------------------- */
double
NFrobenius(double m[3][3])
{
	int i, j;
	double nf = 0.0;

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			nf = nf + (m[i][j] * m[i][j]);
		}
	return(sqrt(nf));
}
/* NFrobenius */

double
Q_Knupp(const T_Malla *m, T_Tetra t)
{
	double s[3][3];
	//double s_o[3][3];
	
	double nf, d, sg;
	int i;

	/* Construye la matriz jacobiana */
	for (i = 1; i < 4; i++) {
		s[0][i - 1] = m->nodo[t.nodo[i]].x - m->nodo[t.nodo[0]].x;
		s[1][i - 1] = m->nodo[t.nodo[i]].y - m->nodo[t.nodo[0]].y;
		s[2][i - 1] = m->nodo[t.nodo[i]].z - m->nodo[t.nodo[0]].z;
	}

	/*Matrix3d s_o;
	Matrix3d s_o_I;

	for (i = 1; i < 4; i++) {
		s_o(0, i - 1) = m->nodo[t.nodo[i]].o_x - m->nodo[t.nodo[0]].o_x;
		s_o(1, i - 1) = m->nodo[t.nodo[i]].o_y - m->nodo[t.nodo[0]].o_y;
		s_o(2, i - 1) = m->nodo[t.nodo[i]].o_z - m->nodo[t.nodo[0]].o_z;
	}
	
	s_o_I = s_o.inverse();*/


	{
		double aux[3][3];

#if STANDARD_TET
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				double p = 0.0;
				for (int k = 0; k < 3; k++)
				{
					p += s[i][k] * W0_inversa[k][j];
					//p += s[i][k] * t.IS[k][j];
					//p += s[i][k] * s_o_I(k,j);
				}
				aux[i][j] = p;
			}
		}
#else
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				double p = 0.0;
				for (int k = 0; k < 3; k++)
				{
					//p += s[i][k] * W0_inversa[k][j];
					p += s[i][k] * t.IS[k][j];
					//p += s[i][k] * s_o_I(k,j);
				}
				aux[i][j] = p;
			}
		}
#endif

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				s[i][j] = aux[i][j];
			}
		}
	}

	//  Producto(s, W0_inversa);
	d = sg = Det(&s[0][0]);
	if (sg < peor_sigma) peor_sigma = sg;
	d = d * d;
	d = cbrt(d);
	nf = NFrobenius(s);
	nf = nf * nf;
	if (sg < 0.0) d = -d;
	return((3.0 * d) / nf);
} /* Q_Knupp */

/* 计算网格的所有四面体的质量 */
/* 如果list不为null，则将质量存储在列表中 */
/* 返回质量低于阈值的四面体数 */
int
CalculaCalidades(const T_Malla *m, double *lista, double umbral)
{
	int i, flag = 0;
	double q;

	for (i = 0; i < m->Num_tetra; i++) {
		q = Q_Knupp(m, m->tetra[i]);
		if (q < umbral) flag++;
		if (lista != NULL) lista[i] = q;
		LaMalla.tetra[i].q = q;
	}

	return(flag);
}/* CalculaCalidades */


/* Funcion para comparar doubles, usada por qsort */
int
compara_doubles(const void *a, const void *b)
{
	const double *da = (const double *)a;
	const double *db = (const double *)b;

	return (*da > *db) - (*da < *db);
}


void
OrdenaCalidades(const T_Malla *malla, double *buffer)
{
	int i;

	for (i = 0; i < malla->Num_tetra; i++) buffer[i] = malla->tetra[i].q;
	qsort((void *)buffer, malla->Num_tetra, sizeof(double), compara_doubles);
}




// -------------------------------------------------------------
/* Graba un fichero con los nodos de la malla suavizada */
void
GrabaFicheroSalida(const char *nombre, const T_Malla *malla)
{
	int i;
	FILE *FiCa;

	if ((FiCa = fopen(nombre, "w")) == NULL) {
		fprintf(stderr, "GrabaFicheroSalida: error en el fichero %s\n%s",
			nombre, strerror(errno));
		abort();
	}
	//fprintf(FiCa, "# Numero de nodos de la malla\n%d\n", malla->Num_nodos);
	fprintf(FiCa, "%d\n", malla->Num_nodos);
	for (i = 0; i < malla->Num_nodos; i++)
		fprintf(FiCa, "%.16f  %.16f  %.16f  %d\n",
			malla->nodo[i].x, malla->nodo[i].y, malla->nodo[i].z,
			malla->nodo[i].nr);
	fclose(FiCa);
}/* GrabaFicheroSalida */


/* Graba un fichero con las calidades 记录具有质量的文件 */
void
GrabaFicheroCalidades(double *lista, int numtetra, const char *prefijo)
{
	static int numiter = 0;
	int i;
	char nf_calidad[1024]; /* prefijo del nombre de los ficheros de calidades */
	FILE *FiCa;

	OrdenaCalidades(&LaMalla, lista);
	sprintf(nf_calidad, "%s%03d", prefijo, numiter++);
	if ((FiCa = fopen(nf_calidad, "w")) == NULL) {
		fprintf(stderr, "GrabaFicheroCalidades: error en el fichero %s\n%s\n",
			nf_calidad, strerror(errno));
		abort();
	}
	for (i = 0; i < numtetra; i++)
		fprintf(FiCa, "%.04g\n", (lista[i] < 0.0) ? 0.0 : lista[i]);
	fclose(FiCa);
} /* GrabaFicheroCalidades */


void
GrabaEstadisticas(FILE *f_estadis, double *lista, int numtetra)
{
	int i;
	double media = 0;
	double max = -2;
	double min = 2;
	double aux;
	int deformes = 0;

	OrdenaCalidades(&LaMalla, lista);
	for (i = 0; i < numtetra; i++) {
		aux = lista[i];
		if (aux < 0.0) {
			deformes++;
			aux = 0.0;
		}
		media += aux;
		if (aux < min) min = aux;
		if (aux > max) max = aux;
	}
	media = media / numtetra;

	fprintf(f_estadis, "%.6e %.6e %.6e %d\n", min, media, max, deformes);
	fflush(f_estadis);

} /* GrabaEstadisticas */







#ifdef BOOST_NO_STDC_NAMESPACE
// Needed for compilers like Visual C++ 6 in which the standard
// C library functions are not in the std namespace.
namespace std { using ::printf; }
#endif

// Vector and Matrix classes from SCPPNT
typedef SCPPNT::Vector<double> dvec;
typedef SCPPNT::Matrix<double> dmat;




// Class to use as function template argument to Uncmin for numeric minimization
class BivarFunctionOnly
{
public:

	BivarFunctionOnly();
	~BivarFunctionOnly();

	// Function to minimize
	double f_to_minimize(dvec &p);

	double f_to_minimize(dvec &p, int normal_x, int normal_y, int normal_z);

	// Gradient of function to minimize
	void gradient(dvec &p, dvec &g);

	// Hessian of function
	void hessian(dvec &p, dmat &h);

	// Indicates analytic gradient is used (0 = do not use)
	int HasAnalyticGradient(); //{return 0;}

	// Indicates analytic hessian is not used (0 = do not use)
	int HasAnalyticHessian();// {return 0;}

	// Any real vector will contain valid parameter values
	int ValidParameters(dvec &x) { return 1; }

	// Dimension of problem (2 parameters in quadratic function)
	int dim() { return 3; }

	// 最后计算的梯度
	double Gradiente[3];

	// 最后的Hessian
	double Hessiano[3][3];
};

// Constructor 
BivarFunctionOnly::BivarFunctionOnly()
{
}

// Destructor
BivarFunctionOnly::~BivarFunctionOnly()
{
}

// Function to minimize



double BivarFunctionOnly::f_to_minimize(dvec &p)
{
	double f;
	/*
	gsl_vector *v;
	gsl_vector *grad;
	gsl_matrix *H;
	*/
	double v[3];
	double grad[3];
	double H[3][3];

	/*
	v = gsl_vector_alloc(3);
	gradiente = gsl_vector_alloc(3);
	H = gsl_matrix_alloc(3,3);

	gsl_vector_set(v, 0, p[0]);
	gsl_vector_set(v, 1, p[1]);
	gsl_vector_set(v, 2, p[2]);
	*/

	v[0] = p[0];
	v[1] = p[1];
	v[2] = p[2];


	new_Fknupp_con_grad_y_hess(v, (void *)&params, &f, grad, H, deltaglobal);

	// 保持梯度的结果
	for (int i = 0; i < 3; i++)
		Gradiente[i] = grad[i];
	//Gradiente[i] = gsl_vector_get(gradiente, i);

	// 保持Hessian矩阵（只有下三角）的结果
	for (int i = 0; i < 3; i++)
		for (int j = 0; j <= i; j++)
			Hessiano[i][j] = H[i][j];
	//      Hessiano[i][j] = gsl_matrix_get(H, i, j);

	/*
	gsl_vector_free(v);
	gsl_vector_free(gradiente);
	gsl_matrix_free(H);
	*/
	return  f;
}

#if USING_PROJ

double BivarFunctionOnly::f_to_minimize(dvec &p, int normal_x, int normal_y, int normal_z)
{
	double f;
	/*
	gsl_vector *v;
	gsl_vector *grad;
	gsl_matrix *H;
	*/
	double v[3];
	double grad[3];
	double H[3][3];

	/*
	v = gsl_vector_alloc(3);
	gradiente = gsl_vector_alloc(3);
	H = gsl_matrix_alloc(3,3);

	gsl_vector_set(v, 0, p[0]);
	gsl_vector_set(v, 1, p[1]);
	gsl_vector_set(v, 2, p[2]);
	*/

	//v might need to be changed here when considering boundary mesh
	int face_id = -1;

	if (normal_x > 0)
	{
		//proj
		//std::cout << "boundary pts: " << std::endl;
		Point3 query(p[0], p[1], p[2]);
		Point_and_primitive_id pp = search_tree->closest_point_and_primitive(query);
		face_id = (int)pp.second->id();
		//closet point
		v[0] = pp.first[0];
		v[1] = pp.first[1];
		v[2] = pp.first[2];

		new_Fknupp_con_grad_y_hess_boundary(v, (void *)&params, &f, grad, H, deltaglobal);

		/*std::cout << "po: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
		std::cout << "pn: " << v[0] << " " << v[1] << " " << v[2] << std::endl;*/
		
		/*v[0] = p[0];
		v[1] = p[1];
		v[2] = p[2];*/
	}
	else
	{
		v[0] = p[0];
		v[1] = p[1];
		v[2] = p[2];
		new_Fknupp_con_grad_y_hess(v, (void *)&params, &f, grad, H, deltaglobal);
	}

	


	//new_Fknupp_con_grad_y_hess(v, (void *)&params, &f, grad, H, deltaglobal);

	// 保持梯度的结果
	

	//if (normal_x > 0)
	//{
	//	assert(face_id >= 0);
	//	
	//	Vector3d tmp_g(grad[0], grad[1], grad[2]);
	//	
	///*	if (face_id > proj_normals.size())
	//		std::cout << "proj normals: " << proj_normals[face_id][0] << " " << proj_normals[face_id][1] << " " << proj_normals[face_id][2] << std::endl;*/
	//	double l = proj_normals[face_id].dot(tmp_g);
	//	Vector3d final_g = tmp_g - l * proj_normals[face_id];
	//	for (int i = 0; i < 3; i++)
	//		Gradiente[i] = final_g[i];
	//	
	//	/*for (int i = 0; i < 3; i++)
	//		Gradiente[i] = grad[i];*/
	//}
	//else
	{
		for (int i = 0; i < 3; i++)
			Gradiente[i] = grad[i];
	}

	//	for (int i = 0; i < 3; i++)
	//		Gradiente[i] = 0;

	//if (normal_x == 0)
	//	Gradiente[0] = grad[0];
	//if (normal_y == 0)
	//	Gradiente[1] = grad[1];
	//if (normal_z == 0)
	//	Gradiente[2] = grad[2];

	//Gradiente[i] = gsl_vector_get(gradiente, i);

	// 保持Hessian矩阵（只有下三角）的结果
	for (int i = 0; i < 3; i++)
		for (int j = 0; j <= i; j++)
			Hessiano[i][j] = H[i][j];
			//Hessiano[i][j] = 0;

	//if (normal_x == 1)
	//{
	//	Hessiano[0][0] = 0;
	//	Hessiano[1][0] = 0;
	//	Hessiano[2][0] = 0;
	//}

	//if (normal_y == 1)
	//{
	//	Hessiano[1][1] = 0;
	//	Hessiano[1][0] = 0;
	//	Hessiano[2][1] = 0;
	//}

	//if (normal_z == 1)
	//{
	//	Hessiano[2][0] = 0;
	//	Hessiano[2][1] = 0;
	//	Hessiano[2][2] = 0;
	//}


	/*for (int i = 0; i < 3; i++)
		for (int j = 0; j <= i; j++)*/
	//      Hessiano[i][j] = gsl_matrix_get(H, i, j);

	/*
	gsl_vector_free(v);
	gsl_vector_free(gradiente);
	gsl_matrix_free(H);
	*/
	return  f;
}

#else
double BivarFunctionOnly::f_to_minimize(dvec &p, int normal_x, int normal_y, int normal_z)
{
	double f;
	/*
	gsl_vector *v;
	gsl_vector *grad;
	gsl_matrix *H;
	*/
	double v[3];
	double grad[3];
	double H[3][3];

	/*
	v = gsl_vector_alloc(3);
	gradiente = gsl_vector_alloc(3);
	H = gsl_matrix_alloc(3,3);

	gsl_vector_set(v, 0, p[0]);
	gsl_vector_set(v, 1, p[1]);
	gsl_vector_set(v, 2, p[2]);
	*/

	//v might need to be changed here when considering boundary mesh

	v[0] = p[0];
	v[1] = p[1];
	v[2] = p[2];


	new_Fknupp_con_grad_y_hess(v, (void *)&params, &f, grad, H, deltaglobal);

	// 保持梯度的结果
	/*for (int i = 0; i < 3; i++)
		Gradiente[i] = grad[i];*/

	for (int i = 0; i < 3; i++)
		Gradiente[i] = 0;

	if (normal_x == 0)
		Gradiente[0] = grad[0];
	if (normal_y == 0)
		Gradiente[1] = grad[1];
	if (normal_z == 0)
		Gradiente[2] = grad[2];

	//Gradiente[i] = gsl_vector_get(gradiente, i);

	// 保持Hessian矩阵（只有下三角）的结果
	for (int i = 0; i < 3; i++)
		for (int j = 0; j <= i; j++)
			Hessiano[i][j] = H[i][j];
	//Hessiano[i][j] = 0;

//if (normal_x == 1)
//{
//	Hessiano[0][0] = 0;
//	Hessiano[1][0] = 0;
//	Hessiano[2][0] = 0;
//}

//if (normal_y == 1)
//{
//	Hessiano[1][1] = 0;
//	Hessiano[1][0] = 0;
//	Hessiano[2][1] = 0;
//}

//if (normal_z == 1)
//{
//	Hessiano[2][0] = 0;
//	Hessiano[2][1] = 0;
//	Hessiano[2][2] = 0;
//}


/*for (int i = 0; i < 3; i++)
	for (int j = 0; j <= i; j++)*/
	//      Hessiano[i][j] = gsl_matrix_get(H, i, j);

	/*
	gsl_vector_free(v);
	gsl_vector_free(gradiente);
	gsl_matrix_free(H);
	*/
	return  f;
}
#endif

// Gradient of function to minimize
void BivarFunctionOnly::gradient(dvec &p, dvec &g)
{
	// Not provided for this function class
}

// Hessian of function to minimize
void BivarFunctionOnly::hessian(dvec &p, dmat &h)
{
	// Not provided for this function class
}

// Indicate whether analytic gradient is used (0 = do not use)
int BivarFunctionOnly::HasAnalyticGradient() { return 0; }

// Indicate whether analytic hessian is used (0 = do not use)
int BivarFunctionOnly::HasAnalyticHessian() { return 0; }


// Function with analytic gradient

// Class to use as function template argument to Uncmin for minimization with analytic first derivative

class BivarFunctionAndGradient : public BivarFunctionOnly
{
public:

	BivarFunctionAndGradient();
	~BivarFunctionAndGradient();

	// Gradient of function to minimize
	void gradient(dvec &p, dvec &g);

	// Indicates analytic gradient is used (0 = do not use)
	int HasAnalyticGradient();

	//private: // No private data structures in this function!
};

// Constructor
BivarFunctionAndGradient::BivarFunctionAndGradient()
{
}

// Destructor
BivarFunctionAndGradient::~BivarFunctionAndGradient()
{
}

// Gradient of function to minimize
void BivarFunctionAndGradient::gradient(dvec &p, dvec &g)
{
}

// Indicates whether analytic gradient is used (0 = do not use)
int BivarFunctionAndGradient::HasAnalyticGradient() { return 1; }

// Function with both analytic Gradient and Hessian

// Class to use as function template argument to Uncmin for minimization with analytic first adn second derivatives
class BivarFunctionGradientAndHessian : public BivarFunctionOnly
{
public:
	BivarFunctionGradientAndHessian();
	~BivarFunctionGradientAndHessian();

	// Gradient of function to minimize
	void gradient(dvec &p, dvec &g);

	// Indicates analytic gradient is used (0 = do not use)
	int HasAnalyticGradient();

	// Hessian of function
	void hessian(dvec &p, dmat &h);

	// Indicates analytic hessian is not used (0 = do not use)
	int HasAnalyticHessian();

	//private: // No private data structures in this function!
};

// Constructor
BivarFunctionGradientAndHessian::BivarFunctionGradientAndHessian()
{
	// Pone a 0 la triangular superior del Hessiano
	Hessiano[0][1] = Hessiano[0][2] = Hessiano[1][2] = 0.0;
}

// Destructor
BivarFunctionGradientAndHessian::~BivarFunctionGradientAndHessian()
{
}

// Definicion del gradiente
void BivarFunctionGradientAndHessian::gradient(dvec &p, dvec &g)
{
	// Devuelvo el ultimo gradiente calculado
	g[0] = Gradiente[0];
	g[1] = Gradiente[1];
	g[2] = Gradiente[2];

#ifdef DEBUG
	cout << "Gradient is: " <<
		g[0] << " " << g[1] << " " << g[2] << endl;
#endif
}

// Hessian of function to minimize
void BivarFunctionGradientAndHessian::hessian(dvec &p, dmat &h)
{
	// Devuelvo el ultimo Hessiano calculado
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			h[i][j] = Hessiano[i][j];

#ifdef DEBUG
	cout << "Hessian is: " <<
		h[0][0] << " " << h[0][1] << " " << h[0][2] << endl <<
		h[1][0] << " " << h[1][1] << " " << h[1][2] << endl <<
		h[2][0] << " " << h[2][1] << " " << h[2][2] << endl;
#endif
}

// Indicates whether analytic hessian is  used (0 = do not use)
int BivarFunctionGradientAndHessian::HasAnalyticHessian() { return 1; }
int BivarFunctionGradientAndHessian::HasAnalyticGradient() { return 1; }



void uso(void)
{
	cout << "sus " << version << endl;
	cout << "Usage:" << endl;
	cout << "sus -n fnodes -e felems [-o fsmooth] [-s fstat] [-d] [-q qprefix] "
		<< "[-i max_iter] [-m smooth_iter]" << endl;
	cout << "Mandatory arguments:" << endl;
	cout << "  fnodes  : file of nodes" << endl;
	cout << "  felems  : file of elements" << endl;
	cout << "Optional arguments" << endl;
	cout << "    fsmooth : file of smoothed nodes" << endl;
	cout << "     fstats : file of statistics" << endl;
	cout << "    qprefix : file name prefix for quality files" << endl;
	cout << "   max_iter : maximum number of iterations (default, "
		<< MaxNumIter << ")" << endl;
	cout << "smooth_iter : number of smoothing iterations (default, "
		<< NumIterSuavizado << ")" << endl;
	cout << "         -d : do not normalize input mesh" << endl;
}


//bool test()
//{
//	return true;
//}


bool SUS(const std::vector<double> &nodes_deform, const std::vector<int> &ref_nodes, const std::vector<unsigned int> &elems, std::vector<int> &node_normal_tag, const std::vector<double> &node_ori, std::vector<double> &output_coord, int max_smooth_iter, int max_iter)
{
	//std::cout << "test SUS" << std::endl;
	bool Normalized = true;
	bool SaveQuality = false;
	char *fnodes = NULL; // Nodes file (input)
	char *felems = NULL; // Element file (input)
	char *fstat = NULL; // Smooth statistics file (output)
	char *fsmoot = NULL; // Smoothed mesh file (output). Only nodes are changed during smoothing
	char *fquality = NULL; // File(s) quality prefix
	char *ftetfile = NULL; // tet file (input)
	char *fotetfile = NULL; // tet file (output)
	int optchar;


	/*Eigen::Matrix2d temp_m;
	Eigen::Matrix2d m_I;
	temp_m(0, 0) = 0;
	temp_m(1, 0) = 0;
	temp_m(0, 1) = 0;
	temp_m(1, 1) = temp_m(1, 0) + temp_m(0, 1);
	m_I = temp_m.inverse();*/




	s.resize(MAX_CARAS_VISTAS * 9); sigma.resize(MAX_CARAS_VISTAS);

	/*if ((fnodes == NULL) || (felems == NULL))
	{
		cout << "Missing file(s) name(s)" << endl;
		uso();
		exit(EXIT_FAILURE);
	}*/

	// Dimensionality of the problem 
	const int dim = 3;

	// starting values
	double start_values[3] = { 0, 0, 0 };

	// Find filehandle of standard output stream (used for intermediate minimization output)
	FILE* fh = (FILE*)stdout;

	// xpls will contain solution, gpls will contain
	// gradient at solution after call to min_<f>.Minimize
	dvec xpls(dim), gpls(dim);

	// hpls provides the space for the hessian
	dmat hpls(dim, dim);
	// there are also these derived matrices
	dmat cholesky(dim, dim), Hessian(dim, dim);

	// fpls contains the value of the function at
	// the solution given by xpls.
	double fpls;

	// Read input mesh
	//int n = LeerMalla(fnodes, felems, &LaMalla);
	//int n = InitMalla(nodes_deform, ref_nodes, elems, &LaMalla);
	int n = InitMalla(nodes_deform, ref_nodes, elems, node_ori, node_normal_tag, &LaMalla);

	if (n < 1)
	{
		cerr << "Error reading the mesh" << endl;
		exit(EXIT_FAILURE);
	}

	cout << "Mesh data:" << endl;
	cout << "   nodes   : " << LaMalla.Num_nodos << endl;
	cout << "   elements: " << LaMalla.Num_tetra << endl;

	// Mesh normalization
	if (Normalized)
		NormalizarLaMalla(&LaMalla);

	m = &LaMalla;

	// Set typical magnitude of parameter values (off by one in this example)
	double typical_size[3] = { 0., 0., 0. }; //
	{
		// Toma el tama�o tipico de cada variable como la media
		// de los valores absolutos de las coordenadas de los nodos.
		// Esto sirve para el escalado de las variables de la funcion
		// objetivo
		//
		// Typical magnitude is obtained as the mean absolute value
		for (int i = 0; i < m->Num_nodos; i++)
		{
			typical_size[0] += fabs(m->nodo[i].x);
			typical_size[1] += fabs(m->nodo[i].y);
			typical_size[2] += fabs(m->nodo[i].z);
		}
		typical_size[0] /= m->Num_nodos;
		typical_size[1] /= m->Num_nodos;
		typical_size[2] /= m->Num_nodos;
	}
	dvec typ(typical_size, typical_size + dim);

	// Part 2: Load and minimimize the function. Function, first and
	// second-order derivatives supplied
	//
	// Create function object
	BivarFunctionGradientAndHessian FunctionGradientAndHessian;
	// create Uncmin object
	Uncmin<dvec, dmat, BivarFunctionGradientAndHessian> min_fgh(&FunctionGradientAndHessian);
	// Set typical magnitude of function values at the minimum
	min_fgh.SetScaleFunc(10);

	// Set typical magnitude of argument values at function minimum
	if (min_fgh.SetScaleArg(typ))
	{
		cout << endl << "Error in setting typical value" << endl;
		exit(EXIT_FAILURE);
	}

	// Diagnostic printout piped to stdout (fh, 1, 1)
	// Only results piped to stdout (fh, 1, 0)
	// No info piped to stdout (0, *, *)
	min_fgh.SetPrint(0, 1, 0);

	min_fgh.SetTolerances(-1, -1);

	T_Malla *m = &LaMalla;

	params.mascara = 7;   // minimiza en las 3 dimensiones
	params.ndim = 3;
	params.malla = m;
	params.MargenEpsilon = EpsSafetyFactor; // "Safety" factor for epsilon

	/* Memoria para guardar las calidades de los tetraedros */
	double *lista = (double *)malloc(LaMalla.Num_tetra * sizeof(double));

	// Crea un fichero para grabar las estadisticas
	FILE *f_estadisticas = NULL;
	if (fstat != NULL)
	{
		if ((f_estadisticas = fopen(fstat, "w")) == NULL)
		{
			cout << "Can't create statistics file " << fstat << endl;
		}
		else
		{
			// Statistics file header
			fprintf(f_estadisticas, "   Min.         Median       Max.      # Invalid \n");
			fprintf(f_estadisticas, "   Quality      Quality      Quality   Elements  \n");
			fprintf(f_estadisticas, "--------------------------------------------------\n");
		}
	}

	// Print statistics file header
	fprintf(stdout, "\n             Mesh quality statistics\n");
	fprintf(stdout, "   Min.         Median       Max.      # Invalid \n");
	fprintf(stdout, "   Quality      Quality      Quality   Elements  \n");
	fprintf(stdout, "--------------------------------------------------\n");

	// 计算四面体的质量
	n = CalculaCalidades(&LaMalla, lista, 0.0);

	// 记录质量文件
	if (fquality != NULL)
		GrabaFicheroCalidades(lista, LaMalla.Num_tetra, fquality);

	// 记录网格的统计信息
	if (f_estadisticas != NULL)
		GrabaEstadisticas(f_estadisticas, lista, LaMalla.Num_tetra);

	// 打印统计信息
	GrabaEstadisticas(stdout, lista, LaMalla.Num_tetra);

	int iter_suav = -1;            // 平滑迭代计数器
	// 在网格上迭代循环

	MaxNumIter = max_iter;
	NumIterSuavizado = max_smooth_iter;

	for (int iteraciones = 1; iteraciones <= MaxNumIter; iteraciones++)
	{
		// 在节点中循环
		for (nodo = 0; nodo < m->Num_nodos; nodo++)
		{
			double coor_antes[3];
			double coor_despues[3];
			params.nodo = nodo;
			start_values[0] = params.var[0] = coor_antes[0] = m->nodo[nodo].x;
			start_values[1] = params.var[1] = coor_antes[1] = m->nodo[nodo].y;
			start_values[2] = params.var[2] = coor_antes[2] = m->nodo[nodo].z;

			int normal_x, normal_y, normal_z;
			normal_x = m->nodo[nodo].normal_x;
			normal_y = m->nodo[nodo].normal_y;
			normal_z = m->nodo[nodo].normal_z;

			// 仅最小化参考编号为0的节点
			if (m->nodo[nodo].nr != 0) continue;
			deltaglobal = Calcula_Delta(&params, EpsilonEfectivo);
#ifdef DEBUG
			cout << "Minimizando el nodo " << nodo
				<< " (nref " << LaMalla.nodo[nodo].nr << ")" << endl;
#endif

			//Minimize the function
			int iMethod = 1;  // Line Search 
			// Set Minimization method 
			min_fgh.SetMethod(iMethod);

			dvec start(start_values, start_values + dim);
			// Saves the initial value of the objective function in f_ini
			double f_ini;
			{
				BivarFunctionOnly v;
				//f_ini = v.f_to_minimize(start);
				f_ini = v.f_to_minimize(start, normal_x, normal_y, normal_z);
			}
#ifdef DEBUG
			cout << "Inicialmente la funcion vale " << f_ini << endl;
#endif
			// Minimize function
			//min_fgh.Minimize(start, xpls, fpls, gpls, hpls);			
			min_fgh.Minimize(start, xpls, fpls, gpls, hpls, normal_x, normal_y, normal_z);

			// Find status of function minimization.
			// For a description of the values returned by
			// GetMessage see the documentation at the beginning
			// of Uncmin.h.
			int msg = min_fgh.GetMessage();
#ifdef DEBUG  
			cout << endl << endl << "Message returned from Uncmin: " << msg << endl;
			cout << endl << "Function minimum: " << endl << fpls << endl;
			cout << endl << "Parameters: " << xpls << endl;
			cout << "Gradient: " << gpls << endl;
#endif
			// Get lower triangle extract of hpls, with upper off-diagonal values set to zero.
			cholesky = hpls;
			for (int icolumn = 1; icolumn < cholesky.num_columns(); icolumn++)
			{
				for (int irow = 0; irow < icolumn; irow++)
					cholesky[irow][icolumn] = 0.0;
			}
			Hessian = cholesky * transpose(cholesky);
#ifdef DEBUG
			cout << "Hessian:  " << Hessian << endl << endl;


			cout << "Coords Antes: " << coor_antes[0] << "  "
				<< coor_antes[1] << "  "
				<< coor_antes[2] << endl;

			cout << "Coords Despues: "
				<< m->nodo[nodo].x << "  "
				<< m->nodo[nodo].y << "  "
				<< m->nodo[nodo].z << endl;

			cout << "Coords xpls: " << xpls[0] << "  "
				<< xpls[1] << "  "
				<< xpls[2] << endl;
#endif
			// 如果找到了解决方案，请移动节点
			if ((msg >= 1) && (msg <= 3))
			{
				//std::cout << "solution found" << std::endl;
#ifdef DEBUG
				cout << "main" << endl
					<< "\tEncontrado el optimo (msg=" << msg << ")" << endl
					<< "\tMovido el nodo " << nodo << " de "
					<< start_values[0] << "  "
					<< start_values[1] << "  "
					<< start_values[2] << "  a  "
					<< m->nodo[nodo].x << "  "
					<< m->nodo[nodo].y << "  "
					<< m->nodo[nodo].z << endl;
#endif
			}
			else
			{
				//std::cout << "not found" << std::endl;
				// 优化器尚未找到最佳解决方案：
				// entonces, solo muevo el nodo si la solucion dada por
				// el optimizador es valida y la nueva calidad es mejor
				// que la de partida.
				//
				// Optimum not found. Actual node is moved to new
				// (non-optimum) position iif:
				// a) new position is valid and
				// b) local mesh quality is better than original
#ifdef DEBUG
				cout << "\tNo encontrado el optimo (msg=" << msg
					<< ")" << endl;
#endif	      
				if ((!(std::isfinite(m->nodo[nodo].x) && std::isfinite(m->nodo[nodo].y)
					&& std::isfinite(m->nodo[nodo].z))) || (fpls > f_ini))
				{
#ifdef DEBUG
					cout << "main" << endl
						<< "\tNo se mueve el nodo " << nodo
						<< " de su posicion inicial  "
						<< coor_antes[0] << "  "
						<< coor_antes[1] << "  "
						<< coor_antes[2] << endl;
#endif
					std::cout << "no perfect solution found :" << std::endl;
					m->nodo[nodo].x = coor_antes[0];
					m->nodo[nodo].y = coor_antes[1];
					m->nodo[nodo].z = coor_antes[2];
				}
				else
				{
#ifdef DEBUG
					cout << "\tAun asi, se mueve el nodo " << nodo << " de "
						<< start_values[0] << "  "
						<< start_values[1] << "  "
						<< start_values[2] << "  a  "
						<< m->nodo[nodo].x << "  "
						<< m->nodo[nodo].y << "  "
						<< m->nodo[nodo].z << endl;
#endif
				}
			}

		}
		// Quality of all mesh elements
		n = CalculaCalidades(&LaMalla, lista, 0.0);

		// If all elements are valid (not tangled) count this iteration
		// as an smooth iteration
		if (n == 0)
			iter_suav++;

		// Save file of element qualities
		if (fquality != NULL)
			GrabaFicheroCalidades(lista, LaMalla.Num_tetra, fquality);

		// Save statistics file
		if (f_estadisticas != NULL)
			GrabaEstadisticas(f_estadisticas, lista, LaMalla.Num_tetra);

		// Print statistics
		// Here median means average
		GrabaEstadisticas(stdout, lista, LaMalla.Num_tetra);

		// Quit the for loop if all smooth iterations have been done
		if (iter_suav == NumIterSuavizado)
			break;

	} // main loop

  // Mesh denormalization
	if (Normalized)
		DesnormalizarLaMalla(m);

	// Save untangled mesh to a file
	if (fsmoot != NULL)
		GrabaFicheroSalida(fsmoot, m);
	output_coord.clear();
	for (size_t i = 0; i < m->Num_nodos; i++)
	{
		output_coord.push_back(m->nodo[i].x);
		output_coord.push_back(m->nodo[i].y);
		output_coord.push_back(m->nodo[i].z);
	}

	if (f_estadisticas != NULL)
		fclose(f_estadisticas);
	LiberarMalla(m);
	free(lista);
	return true;
}

bool SUS_proj(const std::vector<double> &nodes_deform, const std::vector<int> &ref_nodes, const std::vector<unsigned int> &elems, std::vector<int> &node_normal_tag, const std::vector<double> &node_ori, std::vector<double> &proj_pts, const std::vector<int> &proj_faces, std::vector<double> &output_coord, int max_smooth_iter, int max_iter)
{
	//search_tree = project_points_prepare(proj_pts, proj_faces, proj_normals);

	
	//std::cout << "test SUS" << std::endl;
	bool Normalized = true;
	bool SaveQuality = false;
	char *fnodes = NULL; // Nodes file (input)
	char *felems = NULL; // Element file (input)
	char *fstat = NULL; // Smooth statistics file (output)
	char *fsmoot = NULL; // Smoothed mesh file (output). Only nodes are changed during smoothing
	char *fquality = NULL; // File(s) quality prefix
	char *ftetfile = NULL; // tet file (input)
	char *fotetfile = NULL; // tet file (output)
	int optchar;


	/*Eigen::Matrix2d temp_m;
	Eigen::Matrix2d m_I;
	temp_m(0, 0) = 0;
	temp_m(1, 0) = 0;
	temp_m(0, 1) = 0;
	temp_m(1, 1) = temp_m(1, 0) + temp_m(0, 1);
	m_I = temp_m.inverse();*/




	s.resize(MAX_CARAS_VISTAS * 9); sigma.resize(MAX_CARAS_VISTAS);

	/*if ((fnodes == NULL) || (felems == NULL))
	{
		cout << "Missing file(s) name(s)" << endl;
		uso();
		exit(EXIT_FAILURE);
	}*/

	// Dimensionality of the problem 
	const int dim = 3;

	// starting values
	double start_values[3] = { 0, 0, 0 };

	// Find filehandle of standard output stream (used for intermediate minimization output)
	FILE* fh = (FILE*)stdout;

	// xpls will contain solution, gpls will contain
	// gradient at solution after call to min_<f>.Minimize
	dvec xpls(dim), gpls(dim);

	// hpls provides the space for the hessian
	dmat hpls(dim, dim);
	// there are also these derived matrices
	dmat cholesky(dim, dim), Hessian(dim, dim);

	// fpls contains the value of the function at
	// the solution given by xpls.
	double fpls;

	// Read input mesh
	//int n = LeerMalla(fnodes, felems, &LaMalla);
	//int n = InitMalla(nodes_deform, ref_nodes, elems, &LaMalla);
	int n = InitMalla(nodes_deform, ref_nodes, elems, node_ori, node_normal_tag, &LaMalla);

	if (n < 1)
	{
		cerr << "Error reading the mesh" << endl;
		exit(EXIT_FAILURE);
	}

	cout << "Mesh data:" << endl;
	cout << "   nodes   : " << LaMalla.Num_nodos << endl;
	cout << "   elements: " << LaMalla.Num_tetra << endl;

	// Mesh normalization
	if (Normalized)
	{
		NormalizarLaMalla(&LaMalla);
		double cx(LaMalla.cx), cy(LaMalla.cy), cz(LaMalla.cz), max_ratio(LaMalla.max_radio);
	//normalize ori
#if 1
		int N_nodes = LaMalla.Num_nodos;
		for (size_t i = 0; i < N_nodes; i++)
		{

			LaMalla.nodo[i].o_x = (LaMalla.nodo[i].o_x - cx) / max_ratio;
			LaMalla.nodo[i].o_y = (LaMalla.nodo[i].o_y - cy) / max_ratio;
			LaMalla.nodo[i].o_z = (LaMalla.nodo[i].o_z - cz) / max_ratio;
		}

		for (size_t i = 0; i < LaMalla.Num_tetra; i++)
		{
			//index start from 0
			T_Tetra *t = &LaMalla.tetra[i];
			//std::cout << "calculating: " << i << std::endl;
			Matrix3d s_o;
			Matrix3d s_o_I;
			for (int j = 1; j < 4; j++) {
				s_o(0, j - 1) = LaMalla.nodo[t->nodo[j]].o_x - LaMalla.nodo[t->nodo[0]].o_x;
				s_o(1, j - 1) = LaMalla.nodo[t->nodo[j]].o_y - LaMalla.nodo[t->nodo[0]].o_y;
				s_o(2, j - 1) = LaMalla.nodo[t->nodo[j]].o_z - LaMalla.nodo[t->nodo[0]].o_z;
			}
			//std::cout << "i: " << i << " det: " << s_o.determinant() << std::endl;
			s_o_I = s_o.inverse();
			for (size_t j = 0; j < 3; j++)
			{
				for (size_t k = 0; k < 3; k++)
				{
					t->IS[j][k] = s_o_I(j, k);
					//t->IS[j][k] = s_o(j, k);
				}
			}
		}

#endif 


		//normalize proj pts
		for (size_t i = 0; i < proj_pts.size() / 3; i++)
		{
			proj_pts[3 * i] = (proj_pts[3 * i] - cx) / max_ratio;
			proj_pts[3 * i + 1] = (proj_pts[3 * i + 1] - cy) / max_ratio;
			proj_pts[3 * i + 2] = (proj_pts[3 * i + 2] - cz) / max_ratio;
		}
	}
		
	//for (size_t i = 0; i < 100; i++)
	//{

	//	/*LaMalla.nodo[i].o_x = (LaMalla.nodo[i].o_x - cx) / max_ratio;
	//	LaMalla.nodo[i].o_y = (LaMalla.nodo[i].o_y - cy) / max_ratio;
	//	LaMalla.nodo[i].o_z = (LaMalla.nodo[i].o_z - cz) / max_ratio;*/

	//	std::cout << "ori pos: " << LaMalla.nodo[i].o_x << " " << LaMalla.nodo[i].o_y << " " << LaMalla.nodo[i].o_z << std::endl;
	//}
	
	


	if (proj_pts.size() != 0)
	{
		//Polyhedron polyhedron;
		polyhedron.clear();
		proj_normals.clear();
		std::vector<Point3> PointList;
		std::vector<std::vector<int>> FaceList;
		std::vector<int> face;

		int input_pts_size = (int)proj_pts.size() / 3;
		int input_face_size = (int)proj_faces.size() / 3;
		for (size_t i = 0; i < input_pts_size; i++)
		{
			//PointList.push_back(Point3(source_points[i][0], source_points[i][1], source_points[i][2]));
			PointList.push_back(Point3(proj_pts[3 * i], proj_pts[3 * i + 1], proj_pts[3 * i + 2]));
		}
		//compute normals
		
		for (size_t i = 0; i < input_face_size; i++)
		{
			face.clear();
			/*face.push_back(source_faces[i][0]);
			face.push_back(source_faces[i][1]);
			face.push_back(source_faces[i][2]);*/
			face.push_back(proj_faces[3 * i + 0]);
			face.push_back(proj_faces[3 * i + 1]);
			face.push_back(proj_faces[3 * i + 2]);

			Vector3d v1, v2, tmp_normal;
			v1[0] = proj_pts[3 * face[1] + 0] - proj_pts[3 * face[0] + 0];
			v1[1] = proj_pts[3 * face[1] + 1] - proj_pts[3 * face[0] + 1];
			v1[2] = proj_pts[3 * face[1] + 2] - proj_pts[3 * face[0] + 2];

			v2[0] = proj_pts[3 * face[2] + 0] - proj_pts[3 * face[0] + 0];
			v2[1] = proj_pts[3 * face[2] + 1] - proj_pts[3 * face[0] + 1];
			v2[2] = proj_pts[3 * face[2] + 2] - proj_pts[3 * face[0] + 2];

			tmp_normal = v1.cross(v2);

			tmp_normal = tmp_normal / tmp_normal.norm();
			/*normals.push_back(tmp_normal[0]);
			normals.push_back(tmp_normal[1]);
			normals.push_back(tmp_normal[2]);*/
			proj_normals.push_back(tmp_normal);




			FaceList.push_back(face);
		}

		/*face.push_back(0);
		face.push_back(1);
		face.push_back(2);
		FaceList.push_back(face);
		face.clear();
		face.push_back(1);
		face.push_back(3);
		face.push_back(2);
		FaceList.push_back(face);*/

		Build_triangle<HalfedgeDS> triangles;
		triangles.ConstructList(PointList, FaceList);
		polyhedron.delegate(triangles);
		int count = 0;
		for (Polyhedron::Face_iterator i = polyhedron.facets_begin(); i != polyhedron.facets_end(); ++i)
			i->id() = count++;
		count = 0;
		for (Polyhedron::Vertex_iterator viter = polyhedron.vertices_begin(); viter != polyhedron.vertices_end(); ++viter)
			viter->id() = count++;

		//CGAL_assertion(polyhedron.is_triangle(polyhedron.halfedges_begin()));
		// constructs AABB tree and computes internal KD-tree 
		// data structure to accelerate distance queries
		search_tree = new Tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		//Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		search_tree->accelerate_distance_queries();

		Point3 query(3.0149025908252440, -1.4642328057617426, 3.2269443344044211);

		Point_and_primitive_id pp = search_tree->closest_point_and_primitive(query);
		int face_id = (int)pp.second->id();

		Point3 closest = pp.first;

	}




	m = &LaMalla;

	// Set typical magnitude of parameter values (off by one in this example)
	double typical_size[3] = { 0., 0., 0. }; //
	{
		// Toma el tama�o tipico de cada variable como la media
		// de los valores absolutos de las coordenadas de los nodos.
		// Esto sirve para el escalado de las variables de la funcion
		// objetivo
		//
		// Typical magnitude is obtained as the mean absolute value
		for (int i = 0; i < m->Num_nodos; i++)
		{
			typical_size[0] += fabs(m->nodo[i].x);
			typical_size[1] += fabs(m->nodo[i].y);
			typical_size[2] += fabs(m->nodo[i].z);
		}
		typical_size[0] /= m->Num_nodos;
		typical_size[1] /= m->Num_nodos;
		typical_size[2] /= m->Num_nodos;
	}
	dvec typ(typical_size, typical_size + dim);

	// Part 2: Load and minimimize the function. Function, first and
	// second-order derivatives supplied
	//
	// Create function object
	BivarFunctionGradientAndHessian FunctionGradientAndHessian;
	// create Uncmin object
	Uncmin<dvec, dmat, BivarFunctionGradientAndHessian> min_fgh(&FunctionGradientAndHessian);
	// Set typical magnitude of function values at the minimum
	min_fgh.SetScaleFunc(10);

	// Set typical magnitude of argument values at function minimum
	if (min_fgh.SetScaleArg(typ))
	{
		cout << endl << "Error in setting typical value" << endl;
		exit(EXIT_FAILURE);
	}

	// Diagnostic printout piped to stdout (fh, 1, 1)
	// Only results piped to stdout (fh, 1, 0)
	// No info piped to stdout (0, *, *)
	min_fgh.SetPrint(0, 1, 0);

	min_fgh.SetTolerances(-1, -1);

	T_Malla *m = &LaMalla;

	params.mascara = 7;   // minimiza en las 3 dimensiones 最小化3个维度
	params.ndim = 3;
	params.malla = m;
	params.MargenEpsilon = EpsSafetyFactor; // "Safety" factor for epsilon

	/* Memoria para guardar las calidades de los tetraedros */
	double *lista = (double *)malloc(LaMalla.Num_tetra * sizeof(double));

	// Crea un fichero para grabar las estadisticas
	FILE *f_estadisticas = NULL;
	if (fstat != NULL)
	{
		if ((f_estadisticas = fopen(fstat, "w")) == NULL)
		{
			cout << "Can't create statistics file " << fstat << endl;
		}
		else
		{
			// Statistics file header
			fprintf(f_estadisticas, "   Min.         Median       Max.      # Invalid \n");
			fprintf(f_estadisticas, "   Quality      Quality      Quality   Elements  \n");
			fprintf(f_estadisticas, "--------------------------------------------------\n");
		}
	}

	// Print statistics file header
	fprintf(stdout, "\n             Mesh quality statistics\n");
	fprintf(stdout, "   Min.         Median       Max.      # Invalid \n");
	fprintf(stdout, "   Quality      Quality      Quality   Elements  \n");
	fprintf(stdout, "--------------------------------------------------\n");

	// 计算四面体的质量
	n = CalculaCalidades(&LaMalla, lista, 0.0);

	// 记录质量文件
	if (fquality != NULL)
		GrabaFicheroCalidades(lista, LaMalla.Num_tetra, fquality);

	// 记录网格的统计信息
	if (f_estadisticas != NULL)
		GrabaEstadisticas(f_estadisticas, lista, LaMalla.Num_tetra);

	// 打印统计信息
	GrabaEstadisticas(stdout, lista, LaMalla.Num_tetra);

	int iter_suav = -1;            // 平滑迭代计数器
	// 在网格上迭代循环

	MaxNumIter = max_iter;
	NumIterSuavizado = max_smooth_iter;

	for (int iteraciones = 1; iteraciones <= MaxNumIter; iteraciones++)
	{
		// 在节点中循环
		//first boundary
		for (nodo = 0; nodo < m->Num_nodos; nodo++)
		{
			double coor_antes[3];
			double coor_despues[3];
			params.nodo = nodo;
			start_values[0] = params.var[0] = coor_antes[0] = m->nodo[nodo].x;
			start_values[1] = params.var[1] = coor_antes[1] = m->nodo[nodo].y;
			start_values[2] = params.var[2] = coor_antes[2] = m->nodo[nodo].z;

			int normal_x, normal_y, normal_z;
			normal_x = m->nodo[nodo].normal_x;
			normal_y = m->nodo[nodo].normal_y;
			normal_z = m->nodo[nodo].normal_z;

			// 仅最小化参考编号为0的节点
			if (m->nodo[nodo].nr != 0) continue;

			//boundary

			/*if (iteraciones % 2 == 0)
			{
				if (normal_x > 0) continue;
			}
			else*/
			{
				if (normal_x == 0) continue;
			}
			

			deltaglobal = Calcula_Delta(&params, EpsilonEfectivo);
#ifdef DEBUG
			cout << "Minimizando el nodo " << nodo
				<< " (nref " << LaMalla.nodo[nodo].nr << ")" << endl;
#endif

			//Minimize the function
			int iMethod = 1;  // Line Search 
			// Set Minimization method 
			min_fgh.SetMethod(iMethod);

			dvec start(start_values, start_values + dim);
			// Saves the initial value of the objective function in f_ini
			double f_ini;
			{
				BivarFunctionOnly v;
				//f_ini = v.f_to_minimize(start);
				f_ini = v.f_to_minimize(start, normal_x, normal_y, normal_z);
				//normal y can be used to determine thread
			}
#ifdef DEBUG
			cout << "Inicialmente la funcion vale " << f_ini << endl;
#endif
			// Minimize function
			//min_fgh.Minimize(start, xpls, fpls, gpls, hpls);			
			min_fgh.Minimize(start, xpls, fpls, gpls, hpls, normal_x, normal_y, normal_z);


			//project solution here
//#if USING_PROJ
//			
//#endif

			//continue;
			// Find status of function minimization.
			// For a description of the values returned by
			// GetMessage see the documentation at the beginning
			// of Uncmin.h.
			int msg = min_fgh.GetMessage();
#ifdef DEBUG  
			cout << endl << endl << "Message returned from Uncmin: " << msg << endl;
			cout << endl << "Function minimum: " << endl << fpls << endl;
			cout << endl << "Parameters: " << xpls << endl;
			cout << "Gradient: " << gpls << endl;
#endif
			// Get lower triangle extract of hpls, with upper off-diagonal values set to zero.
			cholesky = hpls;
			for (int icolumn = 1; icolumn < cholesky.num_columns(); icolumn++)
			{
				for (int irow = 0; irow < icolumn; irow++)
					cholesky[irow][icolumn] = 0.0;
			}
			Hessian = cholesky * transpose(cholesky);
#ifdef DEBUG
			cout << "Hessian:  " << Hessian << endl << endl;


			cout << "Coords Antes: " << coor_antes[0] << "  "
				<< coor_antes[1] << "  "
				<< coor_antes[2] << endl;

			cout << "Coords Despues: "
				<< m->nodo[nodo].x << "  "
				<< m->nodo[nodo].y << "  "
				<< m->nodo[nodo].z << endl;

			cout << "Coords xpls: " << xpls[0] << "  "
				<< xpls[1] << "  "
				<< xpls[2] << endl;
#endif
			// 如果找到了解决方案，请移动节点
			if ((msg >= 1) && (msg <= 3))
			{
				//std::cout << "solution found" << std::endl;
#ifdef DEBUG
				cout << "main" << endl
					<< "\tEncontrado el optimo (msg=" << msg << ")" << endl
					<< "\tMovido el nodo " << nodo << " de "
					<< start_values[0] << "  "
					<< start_values[1] << "  "
					<< start_values[2] << "  a  "
					<< m->nodo[nodo].x << "  "
					<< m->nodo[nodo].y << "  "
					<< m->nodo[nodo].z << endl;
#endif
			}
			else
			{
				//std::cout << "not found" << std::endl;
				// 优化器尚未找到最佳解决方案：
				// entonces, solo muevo el nodo si la solucion dada por
				// el optimizador es valida y la nueva calidad es mejor
				// que la de partida.
				//
				// Optimum not found. Actual node is moved to new
				// (non-optimum) position iif:
				// a) new position is valid and
				// b) local mesh quality is better than original
#ifdef DEBUG
				cout << "\tNo encontrado el optimo (msg=" << msg
					<< ")" << endl;
#endif	      
				if ((!(std::isfinite(m->nodo[nodo].x) && std::isfinite(m->nodo[nodo].y)
					&& std::isfinite(m->nodo[nodo].z))) || (fpls > f_ini))
				{
#ifdef DEBUG
					cout << "main" << endl
						<< "\tNo se mueve el nodo " << nodo
						<< " de su posicion inicial  "
						<< coor_antes[0] << "  "
						<< coor_antes[1] << "  "
						<< coor_antes[2] << endl;
#endif
					std::cout << "no perfect solution found :" << std::endl;
					m->nodo[nodo].x = coor_antes[0];
					m->nodo[nodo].y = coor_antes[1];
					m->nodo[nodo].z = coor_antes[2];
				}
				else
				{
#ifdef DEBUG
					cout << "\tAun asi, se mueve el nodo " << nodo << " de "
						<< start_values[0] << "  "
						<< start_values[1] << "  "
						<< start_values[2] << "  a  "
						<< m->nodo[nodo].x << "  "
						<< m->nodo[nodo].y << "  "
						<< m->nodo[nodo].z << endl;
#endif
				}
			}

		}
		// Quality of all mesh elements


		n = CalculaCalidades(&LaMalla, lista, 0.0);

		// If all elements are valid (not tangled) count this iteration
		// as an smooth iteration
		if (n == 0)
			iter_suav++;

		// Save file of element qualities
		if (fquality != NULL)
			GrabaFicheroCalidades(lista, LaMalla.Num_tetra, fquality);

		// Save statistics file
		if (f_estadisticas != NULL)
			GrabaEstadisticas(f_estadisticas, lista, LaMalla.Num_tetra);

		// Print statistics
		// Here median means average
		GrabaEstadisticas(stdout, lista, LaMalla.Num_tetra);

		// Quit the for loop if all smooth iterations have been done
		if (iter_suav == NumIterSuavizado)
			break;

		
		//then interior points


#if 1
		for (nodo = 0; nodo < m->Num_nodos; nodo++)
		{
			double coor_antes[3];
			double coor_despues[3];
			params.nodo = nodo;
			start_values[0] = params.var[0] = coor_antes[0] = m->nodo[nodo].x;
			start_values[1] = params.var[1] = coor_antes[1] = m->nodo[nodo].y;
			start_values[2] = params.var[2] = coor_antes[2] = m->nodo[nodo].z;

			int normal_x, normal_y, normal_z;
			normal_x = m->nodo[nodo].normal_x;
			normal_y = m->nodo[nodo].normal_y;
			normal_z = m->nodo[nodo].normal_z;

			// 仅最小化参考编号为0的节点
			if (m->nodo[nodo].nr != 0) continue;

			//interior points
			
			if (normal_x > 0) continue;

			deltaglobal = Calcula_Delta(&params, EpsilonEfectivo);
#ifdef DEBUG
			cout << "Minimizando el nodo " << nodo
				<< " (nref " << LaMalla.nodo[nodo].nr << ")" << endl;
#endif

			//Minimize the function
			int iMethod = 1;  // Line Search 
			// Set Minimization method 
			min_fgh.SetMethod(iMethod);

			dvec start(start_values, start_values + dim);
			// Saves the initial value of the objective function in f_ini
			double f_ini;
			{
				BivarFunctionOnly v;
				//f_ini = v.f_to_minimize(start);
				f_ini = v.f_to_minimize(start, normal_x, normal_y, normal_z);
				//normal y can be used to determine thread
			}
#ifdef DEBUG
			cout << "Inicialmente la funcion vale " << f_ini << endl;
#endif
			// Minimize function
			//min_fgh.Minimize(start, xpls, fpls, gpls, hpls);			
			min_fgh.Minimize(start, xpls, fpls, gpls, hpls, normal_x, normal_y, normal_z);

			// Find status of function minimization.
			// For a description of the values returned by
			// GetMessage see the documentation at the beginning
			// of Uncmin.h.
			int msg = min_fgh.GetMessage();
#ifdef DEBUG  
			cout << endl << endl << "Message returned from Uncmin: " << msg << endl;
			cout << endl << "Function minimum: " << endl << fpls << endl;
			cout << endl << "Parameters: " << xpls << endl;
			cout << "Gradient: " << gpls << endl;
#endif
			// Get lower triangle extract of hpls, with upper off-diagonal values set to zero.
			cholesky = hpls;
			for (int icolumn = 1; icolumn < cholesky.num_columns(); icolumn++)
			{
				for (int irow = 0; irow < icolumn; irow++)
					cholesky[irow][icolumn] = 0.0;
			}
			Hessian = cholesky * transpose(cholesky);
#ifdef DEBUG
			cout << "Hessian:  " << Hessian << endl << endl;


			cout << "Coords Antes: " << coor_antes[0] << "  "
				<< coor_antes[1] << "  "
				<< coor_antes[2] << endl;

			cout << "Coords Despues: "
				<< m->nodo[nodo].x << "  "
				<< m->nodo[nodo].y << "  "
				<< m->nodo[nodo].z << endl;

			cout << "Coords xpls: " << xpls[0] << "  "
				<< xpls[1] << "  "
				<< xpls[2] << endl;
#endif
			// 如果找到了解决方案，请移动节点
			if ((msg >= 1) && (msg <= 3))
			{
				//std::cout << "solution found" << std::endl;
#ifdef DEBUG
				cout << "main" << endl
					<< "\tEncontrado el optimo (msg=" << msg << ")" << endl
					<< "\tMovido el nodo " << nodo << " de "
					<< start_values[0] << "  "
					<< start_values[1] << "  "
					<< start_values[2] << "  a  "
					<< m->nodo[nodo].x << "  "
					<< m->nodo[nodo].y << "  "
					<< m->nodo[nodo].z << endl;
#endif
			}
			else
			{
				//std::cout << "not found" << std::endl;
				// 优化器尚未找到最佳解决方案：
				// entonces, solo muevo el nodo si la solucion dada por
				// el optimizador es valida y la nueva calidad es mejor
				// que la de partida.
				//
				// Optimum not found. Actual node is moved to new
				// (non-optimum) position iif:
				// a) new position is valid and
				// b) local mesh quality is better than original
#ifdef DEBUG
				cout << "\tNo encontrado el optimo (msg=" << msg
					<< ")" << endl;
#endif	      
				if ((!(std::isfinite(m->nodo[nodo].x) && std::isfinite(m->nodo[nodo].y)
					&& std::isfinite(m->nodo[nodo].z))) || (fpls > f_ini))
				{
#ifdef DEBUG
					cout << "main" << endl
						<< "\tNo se mueve el nodo " << nodo
						<< " de su posicion inicial  "
						<< coor_antes[0] << "  "
						<< coor_antes[1] << "  "
						<< coor_antes[2] << endl;
#endif
					std::cout << "no perfect solution found :" << std::endl;
					m->nodo[nodo].x = coor_antes[0];
					m->nodo[nodo].y = coor_antes[1];
					m->nodo[nodo].z = coor_antes[2];
				}
				else
				{
#ifdef DEBUG
					cout << "\tAun asi, se mueve el nodo " << nodo << " de "
						<< start_values[0] << "  "
						<< start_values[1] << "  "
						<< start_values[2] << "  a  "
						<< m->nodo[nodo].x << "  "
						<< m->nodo[nodo].y << "  "
						<< m->nodo[nodo].z << endl;
#endif
				}
			}

		}

		//data changed here
		n = CalculaCalidades(&LaMalla, lista, 0.0);

		// If all elements are valid (not tangled) count this iteration
		// as an smooth iteration
		if (n == 0)
			iter_suav++;

		// Save file of element qualities
		if (fquality != NULL)
			GrabaFicheroCalidades(lista, LaMalla.Num_tetra, fquality);

		// Save statistics file
		if (f_estadisticas != NULL)
			GrabaEstadisticas(f_estadisticas, lista, LaMalla.Num_tetra);

		// Print statistics
		// Here median means average
		GrabaEstadisticas(stdout, lista, LaMalla.Num_tetra);

		// Quit the for loop if all smooth iterations have been done
		if (iter_suav == NumIterSuavizado)
			break;

#endif


	} // main loop

  // Mesh denormalization
	if (Normalized)
		DesnormalizarLaMalla(m);

	// Save untangled mesh to a file
	if (fsmoot != NULL)
		GrabaFicheroSalida(fsmoot, m);
	output_coord.clear();
	for (size_t i = 0; i < m->Num_nodos; i++)
	{
		output_coord.push_back(m->nodo[i].x);
		output_coord.push_back(m->nodo[i].y);
		output_coord.push_back(m->nodo[i].z);
	}


	//if used multiple times, search_tree can be preserved
	if (!search_tree)
		delete search_tree;
	//if ()

	if (f_estadisticas != NULL)
		fclose(f_estadisticas);

	//free memory
	LiberarMalla(m);
	free(lista);
	return true;
}

typedef OpenMesh::PolyMesh_ArrayKernelT< OpenMesh::DefaultTraits>  Mesh;

void construct(Mesh *mesh, const std::vector<double> &pts, const std::vector<int> &faces)
{
	mesh->clear();
	int pts_size = pts.size() / 3;
	int face_size = faces.size() / 3;
	std::vector<Mesh::VertexHandle> vhandle;
	//std::vector<Mesh::FaceHandle> fhandle;
	vhandle.resize(pts_size);
	//fhandle.resize(face_size);

	for (int i = 0; i < pts_size; i++)
	{
		vhandle[i] = mesh->add_vertex(Mesh::Point(pts[3 * i], pts[3 * i + 1], pts[3 * i + 2]));
	}

	std::vector<Mesh::VertexHandle>  tmp_face_vhandles;
	for (size_t i = 0; i < face_size; i++)
	{
		tmp_face_vhandles.clear();
		tmp_face_vhandles.push_back(vhandle[faces[3 * i]]);
		tmp_face_vhandles.push_back(vhandle[faces[3 * i + 1]]);
		tmp_face_vhandles.push_back(vhandle[faces[3 * i + 2]]);
		mesh->add_face(tmp_face_vhandles);
	}

	mesh->update_normals();

}

void calc_weights(Mesh *mesh, OpenMesh::EPropHandleT<Mesh::Scalar>& eweight_)
{
	Mesh::VertexIter        v_it, v_end(mesh->vertices_end());
	Mesh::EdgeIter          e_it, e_end(mesh->edges_end());
	Mesh::VertexFaceIter    vf_it;
	Mesh::FaceVertexIter    fv_it;
	Mesh::HalfedgeHandle    h0, h1, h2;
	Mesh::VertexHandle      v0, v1;
	Mesh::Point             p0, p1, p2, d0, d1;
	Mesh::Scalar            w, area, b(0.99);



	for (e_it = mesh->edges_begin(); e_it != e_end; ++e_it)
	{
		w = 0.0;

		h0 = mesh->halfedge_handle(e_it.handle(), 0);
		v0 = mesh->to_vertex_handle(h0);
		p0 = mesh->point(v0);

		h1 = mesh->halfedge_handle(e_it.handle(), 1);
		v1 = mesh->to_vertex_handle(h1);
		p1 = mesh->point(v1);

		//if (!mesh->is_boundary(h0))
		{
			h2 = mesh->next_halfedge_handle(h0);
			p2 = mesh->point(mesh->to_vertex_handle(h2));
			d0 = (p0 - p2).normalize();
			d1 = (p1 - p2).normalize();
			w += 1.0 / tan(acos(std::max(-b, std::min(b, (d0 | d1)))));
		}

		//if (!mesh->is_boundary(h1))
		{
			h2 = mesh->next_halfedge_handle(h1);
			p2 = mesh->point(mesh->to_vertex_handle(h2));
			d0 = (p0 - p2).normalize();
			d1 = (p1 - p2).normalize();
			w += 1.0 / tan(acos(std::max(-b, std::min(b, (d0 | d1)))));
		}

		// force weights to be non-negative for higher robustness
		w = std::max(w, 0.0f);

		//weight(e_it) = w;
		mesh->property(eweight_, e_it) = w;
	}
}

void Laplacian_smooth(Mesh &mesh_, OpenMesh::EPropHandleT<Mesh::Scalar>& eweight_, int max_iter, double damp_ratio)
{
	mesh_.request_face_normals();
	Mesh::VertexIter        v_it, v_end(mesh_.vertices_end());
	Mesh::HalfedgeHandle    h;
	Mesh::EdgeHandle        e;
	Mesh::VertexVertexIter  vv_it;
	Mesh::Point             laplace(0.0, 0.0, 0.0);
	Mesh::Scalar            w, ww;
	//double            w, ww;

	OpenMesh::VPropHandleT<Mesh::Point>   vpos_;

	mesh_.add_property(vpos_);


	for (unsigned int iter = 0; iter < max_iter; ++iter)
	{

		// compute new vertex positions by Laplacian smoothing
		for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
		{
			laplace = Mesh::Point(0, 0, 0);
			ww = 0.0;

			if (!mesh_.is_boundary(v_it))
			{
				for (vv_it = mesh_.vv_iter(v_it); vv_it; ++vv_it)
				{
					h = vv_it.current_halfedge_handle();
					e = mesh_.edge_handle(h);
					w = mesh_.property(eweight_, e);
					//w = 1;
					ww += w;

					laplace += w * (mesh_.point(vv_it) - mesh_.point(v_it));
				}

				laplace /= ww;   // normalize by sum of weights
				laplace *= damp_ratio;  // damping
			}

			mesh_.property(vpos_, v_it) = mesh_.point(v_it) + laplace;

			/*if (iter == 4)
			{
				std::cout << "idx: " << v_it.handle().idx() << " : " << mesh_.point(v_it)[0] << " " << mesh_.point(v_it)[1] << " " << mesh_.point(v_it)[2] << std::endl;
				std::cout << "laplace: " << laplace[0] << " " << laplace[1] << " " << laplace[2] << std::endl;
				std::cout << "new pos: " << mesh_.property(vpos_, v_it)[0] << " " << mesh_.property(vpos_, v_it)[1] << " " << mesh_.property(vpos_, v_it)[2] << std::endl;

			}*/


		}



		// update vertex positions
		for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
			mesh_.set_point(v_it, mesh_.property(vpos_, v_it));
	}


	// update face and vertex normals
	mesh_.update_normals();
}


bool SUS_proj_laplacian(const std::vector<double> &nodes_deform, const std::vector<int> &ref_nodes, const std::vector<unsigned int> &elems, std::vector<int> &node_normal_tag, const std::vector<double> &node_ori, std::vector<double> &proj_pts, const std::vector<int> &proj_faces, std::vector<double> &output_coord, const std::vector<int> &boundary_idx, const std::vector<int> &boundary_faces, int max_smooth_iter, int max_iter)
{
	//search_tree = project_points_prepare(proj_pts, proj_faces, proj_normals);


	//std::cout << "test SUS" << std::endl;
	bool Normalized = true;
	bool SaveQuality = false;
	char *fnodes = NULL; // Nodes file (input)
	char *felems = NULL; // Element file (input)
	char *fstat = NULL; // Smooth statistics file (output)
	char *fsmoot = NULL; // Smoothed mesh file (output). Only nodes are changed during smoothing
	char *fquality = NULL; // File(s) quality prefix
	char *ftetfile = NULL; // tet file (input)
	char *fotetfile = NULL; // tet file (output)
	int optchar;


	/*Eigen::Matrix2d temp_m;
	Eigen::Matrix2d m_I;
	temp_m(0, 0) = 0;
	temp_m(1, 0) = 0;
	temp_m(0, 1) = 0;
	temp_m(1, 1) = temp_m(1, 0) + temp_m(0, 1);
	m_I = temp_m.inverse();*/




	s.resize(MAX_CARAS_VISTAS * 9); sigma.resize(MAX_CARAS_VISTAS);

	/*if ((fnodes == NULL) || (felems == NULL))
	{
		cout << "Missing file(s) name(s)" << endl;
		uso();
		exit(EXIT_FAILURE);
	}*/

	// Dimensionality of the problem 
	const int dim = 3;

	// starting values
	double start_values[3] = { 0, 0, 0 };

	// Find filehandle of standard output stream (used for intermediate minimization output)
	FILE* fh = (FILE*)stdout;

	// xpls will contain solution, gpls will contain
	// gradient at solution after call to min_<f>.Minimize
	dvec xpls(dim), gpls(dim);

	// hpls provides the space for the hessian
	dmat hpls(dim, dim);
	// there are also these derived matrices
	dmat cholesky(dim, dim), Hessian(dim, dim);

	// fpls contains the value of the function at
	// the solution given by xpls.
	double fpls;


	//construct surfacemesh
	Mesh mesh_;

	mesh_.request_face_status();
	mesh_.request_edge_status();
	mesh_.request_vertex_status();

	


	// Read input mesh
	//int n = LeerMalla(fnodes, felems, &LaMalla);
	//int n = InitMalla(nodes_deform, ref_nodes, elems, &LaMalla);
	int n = InitMalla(nodes_deform, ref_nodes, elems, node_ori, node_normal_tag, &LaMalla);

	if (n < 1)
	{
		cerr << "Error reading the mesh" << endl;
		exit(EXIT_FAILURE);
	}

	cout << "Mesh data:" << endl;
	cout << "   nodes   : " << LaMalla.Num_nodos << endl;
	cout << "   elements: " << LaMalla.Num_tetra << endl;

	// Mesh normalization
	if (Normalized)
	{
		NormalizarLaMalla(&LaMalla);
		double cx(LaMalla.cx), cy(LaMalla.cy), cz(LaMalla.cz), max_ratio(LaMalla.max_radio);
		//normalize ori
#if 1
		int N_nodes = LaMalla.Num_nodos;
		for (size_t i = 0; i < N_nodes; i++)
		{

			LaMalla.nodo[i].o_x = (LaMalla.nodo[i].o_x - cx) / max_ratio;
			LaMalla.nodo[i].o_y = (LaMalla.nodo[i].o_y - cy) / max_ratio;
			LaMalla.nodo[i].o_z = (LaMalla.nodo[i].o_z - cz) / max_ratio;
		}

		for (size_t i = 0; i < LaMalla.Num_tetra; i++)
		{
			//index start from 0
			T_Tetra *t = &LaMalla.tetra[i];
			//std::cout << "calculating: " << i << std::endl;
			Matrix3d s_o;
			Matrix3d s_o_I;
			for (int j = 1; j < 4; j++) {
				s_o(0, j - 1) = LaMalla.nodo[t->nodo[j]].o_x - LaMalla.nodo[t->nodo[0]].o_x;
				s_o(1, j - 1) = LaMalla.nodo[t->nodo[j]].o_y - LaMalla.nodo[t->nodo[0]].o_y;
				s_o(2, j - 1) = LaMalla.nodo[t->nodo[j]].o_z - LaMalla.nodo[t->nodo[0]].o_z;
			}
			//std::cout << "i: " << i << " det: " << s_o.determinant() << std::endl;
			s_o_I = s_o.inverse();
			for (size_t j = 0; j < 3; j++)
			{
				for (size_t k = 0; k < 3; k++)
				{
					t->IS[j][k] = s_o_I(j, k);
					//t->IS[j][k] = s_o(j, k);
				}
			}
		}

#endif 


		//normalize proj pts
		for (size_t i = 0; i < proj_pts.size() / 3; i++)
		{
			proj_pts[3 * i] = (proj_pts[3 * i] - cx) / max_ratio;
			proj_pts[3 * i + 1] = (proj_pts[3 * i + 1] - cy) / max_ratio;
			proj_pts[3 * i + 2] = (proj_pts[3 * i + 2] - cz) / max_ratio;
		}
	}


	std::vector<double> boundary_pts;

	for (size_t i = 0; i < boundary_idx.size(); i++)
	{
		//normalized version

		int idx = boundary_idx[i];
		boundary_pts.push_back(LaMalla.nodo[idx].x);
		boundary_pts.push_back(LaMalla.nodo[idx].y);
		boundary_pts.push_back(LaMalla.nodo[idx].z);


		/*boundary_pts.push_back(nodes_deform[3 * boundary_idx[i]]);
		boundary_pts.push_back(nodes_deform[3 * boundary_idx[i] + 1]);
		boundary_pts.push_back(nodes_deform[3 * boundary_idx[i] + 2]);*/
	}


	construct(&mesh_, boundary_pts, boundary_faces);

	//OpenMesh::IO::write_mesh(mesh_, "laplacian5_output.off");

	OpenMesh::EPropHandleT<Mesh::Scalar>  eweight_;

	mesh_.add_property(eweight_);


	//for (size_t i = 0; i < 100; i++)
	//{

	//	/*LaMalla.nodo[i].o_x = (LaMalla.nodo[i].o_x - cx) / max_ratio;
	//	LaMalla.nodo[i].o_y = (LaMalla.nodo[i].o_y - cy) / max_ratio;
	//	LaMalla.nodo[i].o_z = (LaMalla.nodo[i].o_z - cz) / max_ratio;*/

	//	std::cout << "ori pos: " << LaMalla.nodo[i].o_x << " " << LaMalla.nodo[i].o_y << " " << LaMalla.nodo[i].o_z << std::endl;
	//}




	if (proj_pts.size() != 0)
	{
		//Polyhedron polyhedron;
		polyhedron.clear();
		proj_normals.clear();
		std::vector<Point3> PointList;
		std::vector<std::vector<int>> FaceList;
		std::vector<int> face;

		int input_pts_size = proj_pts.size() / 3;
		int input_face_size = proj_faces.size() / 3;
		for (size_t i = 0; i < input_pts_size; i++)
		{
			//PointList.push_back(Point3(source_points[i][0], source_points[i][1], source_points[i][2]));
			PointList.push_back(Point3(proj_pts[3 * i], proj_pts[3 * i + 1], proj_pts[3 * i + 2]));
		}
		//compute normals

		for (size_t i = 0; i < input_face_size; i++)
		{
			face.clear();
			/*face.push_back(source_faces[i][0]);
			face.push_back(source_faces[i][1]);
			face.push_back(source_faces[i][2]);*/
			face.push_back(proj_faces[3 * i + 0]);
			face.push_back(proj_faces[3 * i + 1]);
			face.push_back(proj_faces[3 * i + 2]);

			Vector3d v1, v2, tmp_normal;
			v1[0] = proj_pts[3 * face[1] + 0] - proj_pts[3 * face[0] + 0];
			v1[1] = proj_pts[3 * face[1] + 1] - proj_pts[3 * face[0] + 1];
			v1[2] = proj_pts[3 * face[1] + 2] - proj_pts[3 * face[0] + 2];

			v2[0] = proj_pts[3 * face[2] + 0] - proj_pts[3 * face[0] + 0];
			v2[1] = proj_pts[3 * face[2] + 1] - proj_pts[3 * face[0] + 1];
			v2[2] = proj_pts[3 * face[2] + 2] - proj_pts[3 * face[0] + 2];

			tmp_normal = v1.cross(v2);

			tmp_normal = tmp_normal / tmp_normal.norm();
			/*normals.push_back(tmp_normal[0]);
			normals.push_back(tmp_normal[1]);
			normals.push_back(tmp_normal[2]);*/
			proj_normals.push_back(tmp_normal);




			FaceList.push_back(face);
		}

		/*face.push_back(0);
		face.push_back(1);
		face.push_back(2);
		FaceList.push_back(face);
		face.clear();
		face.push_back(1);
		face.push_back(3);
		face.push_back(2);
		FaceList.push_back(face);*/

		Build_triangle<HalfedgeDS> triangles;
		triangles.ConstructList(PointList, FaceList);
		polyhedron.delegate(triangles);
		int count = 0;
		for (Polyhedron::Face_iterator i = polyhedron.facets_begin(); i != polyhedron.facets_end(); ++i)
			i->id() = count++;
		count = 0;
		for (Polyhedron::Vertex_iterator viter = polyhedron.vertices_begin(); viter != polyhedron.vertices_end(); ++viter)
			viter->id() = count++;

		//CGAL_assertion(polyhedron.is_triangle(polyhedron.halfedges_begin()));
		// constructs AABB tree and computes internal KD-tree 
		// data structure to accelerate distance queries
		search_tree = new Tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		//Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		search_tree->accelerate_distance_queries();

		Point3 query(3.0149025908252440, -1.4642328057617426, 3.2269443344044211);

		Point_and_primitive_id pp = search_tree->closest_point_and_primitive(query);
		int face_id = pp.second->id();

		Point3 closest = pp.first;

	}




	m = &LaMalla;

	// Set typical magnitude of parameter values (off by one in this example)
	double typical_size[3] = { 0., 0., 0. }; //
	{
		// Toma el tama�o tipico de cada variable como la media
		// de los valores absolutos de las coordenadas de los nodos.
		// Esto sirve para el escalado de las variables de la funcion
		// objetivo
		//
		// Typical magnitude is obtained as the mean absolute value
		for (int i = 0; i < m->Num_nodos; i++)
		{
			typical_size[0] += fabs(m->nodo[i].x);
			typical_size[1] += fabs(m->nodo[i].y);
			typical_size[2] += fabs(m->nodo[i].z);
		}
		typical_size[0] /= m->Num_nodos;
		typical_size[1] /= m->Num_nodos;
		typical_size[2] /= m->Num_nodos;
	}
	dvec typ(typical_size, typical_size + dim);

	// Part 2: Load and minimimize the function. Function, first and
	// second-order derivatives supplied
	//
	// Create function object
	BivarFunctionGradientAndHessian FunctionGradientAndHessian;
	// create Uncmin object
	Uncmin<dvec, dmat, BivarFunctionGradientAndHessian> min_fgh(&FunctionGradientAndHessian);
	// Set typical magnitude of function values at the minimum
	min_fgh.SetScaleFunc(10);

	// Set typical magnitude of argument values at function minimum
	if (min_fgh.SetScaleArg(typ))
	{
		cout << endl << "Error in setting typical value" << endl;
		exit(EXIT_FAILURE);
	}

	// Diagnostic printout piped to stdout (fh, 1, 1)
	// Only results piped to stdout (fh, 1, 0)
	// No info piped to stdout (0, *, *)
	min_fgh.SetPrint(0, 1, 0);

	min_fgh.SetTolerances(-1, -1);

	T_Malla *m = &LaMalla;

	params.mascara = 7;   // minimiza en las 3 dimensiones 最小化3个维度
	params.ndim = 3;
	params.malla = m;
	params.MargenEpsilon = EpsSafetyFactor; // "Safety" factor for epsilon

	/* Memoria para guardar las calidades de los tetraedros */
	double *lista = (double *)malloc(LaMalla.Num_tetra * sizeof(double));

	// Crea un fichero para grabar las estadisticas
	FILE *f_estadisticas = NULL;
	if (fstat != NULL)
	{
		if ((f_estadisticas = fopen(fstat, "w")) == NULL)
		{
			cout << "Can't create statistics file " << fstat << endl;
		}
		else
		{
			// Statistics file header
			fprintf(f_estadisticas, "   Min.         Median       Max.      # Invalid \n");
			fprintf(f_estadisticas, "   Quality      Quality      Quality   Elements  \n");
			fprintf(f_estadisticas, "--------------------------------------------------\n");
		}
	}

	// Print statistics file header
	fprintf(stdout, "\n             Mesh quality statistics\n");
	fprintf(stdout, "   Min.         Median       Max.      # Invalid \n");
	fprintf(stdout, "   Quality      Quality      Quality   Elements  \n");
	fprintf(stdout, "--------------------------------------------------\n");

	// 计算四面体的质量
	n = CalculaCalidades(&LaMalla, lista, 0.0);

	// 记录质量文件
	if (fquality != NULL)
		GrabaFicheroCalidades(lista, LaMalla.Num_tetra, fquality);

	// 记录网格的统计信息
	if (f_estadisticas != NULL)
		GrabaEstadisticas(f_estadisticas, lista, LaMalla.Num_tetra);

	// 打印统计信息
	GrabaEstadisticas(stdout, lista, LaMalla.Num_tetra);

	int iter_suav = -1;            // 平滑迭代计数器
	// 在网格上迭代循环

	MaxNumIter = max_iter;
	NumIterSuavizado = max_smooth_iter;

	std::vector<double> avg_neighbor_edge_length;
	avg_neighbor_edge_length.resize(m->Num_nodos, -1.0);
	boundary_neighbor_avg_coord.clear();
	boundary_neighbor_avg_coord.resize(3 * m->Num_nodos, -1.0);



	for (int iteraciones = 1; iteraciones <= MaxNumIter; iteraciones++)
	{
		// 在节点中循环


		//first interior points
#if 1
		for (nodo = 0; nodo < m->Num_nodos; nodo++)
		{
			double coor_antes[3];
			double coor_despues[3];
			params.nodo = nodo;
			start_values[0] = params.var[0] = coor_antes[0] = m->nodo[nodo].x;
			start_values[1] = params.var[1] = coor_antes[1] = m->nodo[nodo].y;
			start_values[2] = params.var[2] = coor_antes[2] = m->nodo[nodo].z;

			int normal_x, normal_y, normal_z;
			normal_x = m->nodo[nodo].normal_x;
			normal_y = m->nodo[nodo].normal_y;
			normal_z = m->nodo[nodo].normal_z;

			// 仅最小化参考编号为0的节点
			if (m->nodo[nodo].nr != 0) continue;

			//interior points

			if (normal_x > 0) continue;

			deltaglobal = Calcula_Delta(&params, EpsilonEfectivo);
#ifdef DEBUG
			cout << "Minimizando el nodo " << nodo
				<< " (nref " << LaMalla.nodo[nodo].nr << ")" << endl;
#endif

			//Minimize the function
			int iMethod = 1;  // Line Search 
			// Set Minimization method 
			min_fgh.SetMethod(iMethod);

			dvec start(start_values, start_values + dim);
			// Saves the initial value of the objective function in f_ini
			double f_ini;
			{
				BivarFunctionOnly v;
				//f_ini = v.f_to_minimize(start);
				f_ini = v.f_to_minimize(start, normal_x, normal_y, normal_z);
				//normal y can be used to determine thread
			}
#ifdef DEBUG
			cout << "Inicialmente la funcion vale " << f_ini << endl;
#endif
			// Minimize function
			//min_fgh.Minimize(start, xpls, fpls, gpls, hpls);			
			min_fgh.Minimize(start, xpls, fpls, gpls, hpls, normal_x, normal_y, normal_z);

			// Find status of function minimization.
			// For a description of the values returned by
			// GetMessage see the documentation at the beginning
			// of Uncmin.h.
			int msg = min_fgh.GetMessage();
#ifdef DEBUG  
			cout << endl << endl << "Message returned from Uncmin: " << msg << endl;
			cout << endl << "Function minimum: " << endl << fpls << endl;
			cout << endl << "Parameters: " << xpls << endl;
			cout << "Gradient: " << gpls << endl;
#endif
			// Get lower triangle extract of hpls, with upper off-diagonal values set to zero.
			cholesky = hpls;
			for (int icolumn = 1; icolumn < cholesky.num_columns(); icolumn++)
			{
				for (int irow = 0; irow < icolumn; irow++)
					cholesky[irow][icolumn] = 0.0;
			}
			Hessian = cholesky * transpose(cholesky);
#ifdef DEBUG
			cout << "Hessian:  " << Hessian << endl << endl;


			cout << "Coords Antes: " << coor_antes[0] << "  "
				<< coor_antes[1] << "  "
				<< coor_antes[2] << endl;

			cout << "Coords Despues: "
				<< m->nodo[nodo].x << "  "
				<< m->nodo[nodo].y << "  "
				<< m->nodo[nodo].z << endl;

			cout << "Coords xpls: " << xpls[0] << "  "
				<< xpls[1] << "  "
				<< xpls[2] << endl;
#endif
			// 如果找到了解决方案，请移动节点
			if ((msg >= 1) && (msg <= 3))
			{
				//std::cout << "solution found" << std::endl;
#ifdef DEBUG
				cout << "main" << endl
					<< "\tEncontrado el optimo (msg=" << msg << ")" << endl
					<< "\tMovido el nodo " << nodo << " de "
					<< start_values[0] << "  "
					<< start_values[1] << "  "
					<< start_values[2] << "  a  "
					<< m->nodo[nodo].x << "  "
					<< m->nodo[nodo].y << "  "
					<< m->nodo[nodo].z << endl;
#endif
			}
			else
			{
				//std::cout << "not found" << std::endl;
				// 优化器尚未找到最佳解决方案：
				// entonces, solo muevo el nodo si la solucion dada por
				// el optimizador es valida y la nueva calidad es mejor
				// que la de partida.
				//
				// Optimum not found. Actual node is moved to new
				// (non-optimum) position iif:
				// a) new position is valid and
				// b) local mesh quality is better than original
#ifdef DEBUG
				cout << "\tNo encontrado el optimo (msg=" << msg
					<< ")" << endl;
#endif	      
				if ((!(std::isfinite(m->nodo[nodo].x) && std::isfinite(m->nodo[nodo].y)
					&& std::isfinite(m->nodo[nodo].z))) || (fpls > f_ini))
				{
#ifdef DEBUG
					cout << "main" << endl
						<< "\tNo se mueve el nodo " << nodo
						<< " de su posicion inicial  "
						<< coor_antes[0] << "  "
						<< coor_antes[1] << "  "
						<< coor_antes[2] << endl;
#endif
					std::cout << "no perfect solution found :" << std::endl;
					m->nodo[nodo].x = coor_antes[0];
					m->nodo[nodo].y = coor_antes[1];
					m->nodo[nodo].z = coor_antes[2];
				}
				else
				{
#ifdef DEBUG
					cout << "\tAun asi, se mueve el nodo " << nodo << " de "
						<< start_values[0] << "  "
						<< start_values[1] << "  "
						<< start_values[2] << "  a  "
						<< m->nodo[nodo].x << "  "
						<< m->nodo[nodo].y << "  "
						<< m->nodo[nodo].z << endl;
#endif
				}
			}

		}

		//data changed here
		n = CalculaCalidades(&LaMalla, lista, 0.0);

		// If all elements are valid (not tangled) count this iteration
		// as an smooth iteration
		if (n == 0)
			iter_suav++;

		// Save file of element qualities
		if (fquality != NULL)
			GrabaFicheroCalidades(lista, LaMalla.Num_tetra, fquality);

		// Save statistics file
		if (f_estadisticas != NULL)
			GrabaEstadisticas(f_estadisticas, lista, LaMalla.Num_tetra);

		// Print statistics
		// Here median means average
		GrabaEstadisticas(stdout, lista, LaMalla.Num_tetra);

		// Quit the for loop if all smooth iterations have been done
		if (iter_suav == NumIterSuavizado)
			break;

#endif




		//then boundary

		//get average length of neighbors of boundary pts
		Mesh::VertexIter v_iter = mesh_.vertices_begin();
		Mesh::VertexIter v_end = mesh_.vertices_end();
		Mesh::VertexVertexIter  vv_it;

		v_iter = mesh_.vertices_begin();
		//Mesh::VertexIter v_end = mesh_.vertices_end();

		for (; v_iter != v_end; ++v_iter)
		{
			int idx = boundary_idx[v_iter.handle().idx()];

			/*Point3 query(p[0], p[1], p[2]);
			Point_and_primitive_id pp = search_tree->closest_point_and_primitive(query);*/

			//mesh_.set_point(v_iter, Mesh::Point(m->nodo[idx].x, m->nodo[idx].y, m->nodo[idx].z));
			double edge_length = 0.0;
			double ax = 0.0, ay = 0.0, az = 0.0;
			double count = 0.0;

			Mesh::Point center = mesh_.point(v_iter);
			
			for (vv_it = mesh_.vv_iter(v_iter); vv_it; ++vv_it)
			{
				Mesh::Point pt = mesh_.point(vv_it);
				edge_length += (pt - center).norm();
				ax += pt[0];
				ay += pt[1];
				az += pt[2];
				count = count + 1.0;
			}

			avg_neighbor_edge_length[idx] = edge_length / count;
			boundary_neighbor_avg_coord[3 * idx] = ax / count;
			boundary_neighbor_avg_coord[3 * idx + 1] = ay / count;
			boundary_neighbor_avg_coord[3 * idx + 2] = az / count;


			

		}



		for (nodo = 0; nodo < m->Num_nodos; nodo++)
		{
			double coor_antes[3];
			double coor_despues[3];
			params.nodo = nodo;
			start_values[0] = params.var[0] = coor_antes[0] = m->nodo[nodo].x;
			start_values[1] = params.var[1] = coor_antes[1] = m->nodo[nodo].y;
			start_values[2] = params.var[2] = coor_antes[2] = m->nodo[nodo].z;

			int normal_x, normal_y, normal_z;
			normal_x = m->nodo[nodo].normal_x;
			normal_y = m->nodo[nodo].normal_y;
			normal_z = m->nodo[nodo].normal_z;

			// 仅最小化参考编号为0的节点
			if (m->nodo[nodo].nr != 0) continue;

			//boundary

			/*if (iteraciones % 2 == 0)
			{
				if (normal_x > 0) continue;
			}
			else*/
			{
				if (normal_x == 0) continue;
			}


			deltaglobal = Calcula_Delta(&params, EpsilonEfectivo);
#ifdef DEBUG
			cout << "Minimizando el nodo " << nodo
				<< " (nref " << LaMalla.nodo[nodo].nr << ")" << endl;
#endif

			//Minimize the function
			int iMethod = 1;  // Line Search 
			// Set Minimization method 
			min_fgh.SetMethod(iMethod);

			dvec start(start_values, start_values + dim);
			// Saves the initial value of the objective function in f_ini
			double f_ini;
			{
				BivarFunctionOnly v;
				//f_ini = v.f_to_minimize(start);
				f_ini = v.f_to_minimize(start, normal_x, normal_y, normal_z);
				//normal y can be used to determine thread
			}
#ifdef DEBUG
			cout << "Inicialmente la funcion vale " << f_ini << endl;
#endif
			// Minimize function
			//min_fgh.Minimize(start, xpls, fpls, gpls, hpls);			
			min_fgh.Minimize(start, xpls, fpls, gpls, hpls, normal_x, normal_y, normal_z);

			//if the stepsize is too large, then do not move
#if 1
			if (normal_x > 0)
			{
				double move_distance = 0.0;
				move_distance += (m->nodo[nodo].x - coor_antes[0]) * (m->nodo[nodo].x - coor_antes[0]);
				move_distance += (m->nodo[nodo].y - coor_antes[1]) * (m->nodo[nodo].y - coor_antes[1]);
				move_distance += (m->nodo[nodo].z - coor_antes[2]) * (m->nodo[nodo].z - coor_antes[2]);

				move_distance = sqrt(move_distance);
			
				assert(avg_neighbor_edge_length[nodo] > 0);

				if (move_distance > avg_neighbor_edge_length[nodo])
				{
					std::cout << "stepsize too big" << std::endl;
					m->nodo[nodo].x = coor_antes[0];
					m->nodo[nodo].y = coor_antes[1];
					m->nodo[nodo].z = coor_antes[2];
				}
				else
				{
					//std::cout << "normal stepsize" << std::endl;
				}
			}


			

			//else
			//{
			//	//stepsize not too big also move:
			//	double dx = m->nodo[nodo].x - coor_antes[0];
			//	double dy = m->nodo[nodo].y - coor_antes[1];
			//	double dz = m->nodo[nodo].z - coor_antes[2];

			//	m->nodo[nodo].x = coor_antes[0] + dx / 2.0;
			//	m->nodo[nodo].y = coor_antes[1] + dy / 2.0;
			//	m->nodo[nodo].z = coor_antes[2] + dz / 2.0;
			//	
			//}
			

#endif 


			
//#if USING_PROJ
//			
//#endif

			//continue;
			// Find status of function minimization.
			// For a description of the values returned by
			// GetMessage see the documentation at the beginning
			// of Uncmin.h.
			int msg = min_fgh.GetMessage();
#ifdef DEBUG  
			cout << endl << endl << "Message returned from Uncmin: " << msg << endl;
			cout << endl << "Function minimum: " << endl << fpls << endl;
			cout << endl << "Parameters: " << xpls << endl;
			cout << "Gradient: " << gpls << endl;
#endif
			// Get lower triangle extract of hpls, with upper off-diagonal values set to zero.
			cholesky = hpls;
			for (int icolumn = 1; icolumn < cholesky.num_columns(); icolumn++)
			{
				for (int irow = 0; irow < icolumn; irow++)
					cholesky[irow][icolumn] = 0.0;
			}
			Hessian = cholesky * transpose(cholesky);
#ifdef DEBUG
			cout << "Hessian:  " << Hessian << endl << endl;


			cout << "Coords Antes: " << coor_antes[0] << "  "
				<< coor_antes[1] << "  "
				<< coor_antes[2] << endl;

			cout << "Coords Despues: "
				<< m->nodo[nodo].x << "  "
				<< m->nodo[nodo].y << "  "
				<< m->nodo[nodo].z << endl;

			cout << "Coords xpls: " << xpls[0] << "  "
				<< xpls[1] << "  "
				<< xpls[2] << endl;
#endif
			// 如果找到了解决方案，请移动节点
			if ((msg >= 1) && (msg <= 3))
			{
				//std::cout << "solution found" << std::endl;
#ifdef DEBUG
				cout << "main" << endl
					<< "\tEncontrado el optimo (msg=" << msg << ")" << endl
					<< "\tMovido el nodo " << nodo << " de "
					<< start_values[0] << "  "
					<< start_values[1] << "  "
					<< start_values[2] << "  a  "
					<< m->nodo[nodo].x << "  "
					<< m->nodo[nodo].y << "  "
					<< m->nodo[nodo].z << endl;
#endif
			}
			else
			{
				//std::cout << "not found" << std::endl;
				// 优化器尚未找到最佳解决方案：
				// entonces, solo muevo el nodo si la solucion dada por
				// el optimizador es valida y la nueva calidad es mejor
				// que la de partida.
				//
				// Optimum not found. Actual node is moved to new
				// (non-optimum) position iif:
				// a) new position is valid and
				// b) local mesh quality is better than original
#ifdef DEBUG
				cout << "\tNo encontrado el optimo (msg=" << msg
					<< ")" << endl;
#endif	      
				if ((!(std::isfinite(m->nodo[nodo].x) && std::isfinite(m->nodo[nodo].y)
					&& std::isfinite(m->nodo[nodo].z))) || (fpls > f_ini))
				{
#ifdef DEBUG
					cout << "main" << endl
						<< "\tNo se mueve el nodo " << nodo
						<< " de su posicion inicial  "
						<< coor_antes[0] << "  "
						<< coor_antes[1] << "  "
						<< coor_antes[2] << endl;
#endif
					std::cout << "no perfect solution found :" << std::endl;
					m->nodo[nodo].x = coor_antes[0];
					m->nodo[nodo].y = coor_antes[1];
					m->nodo[nodo].z = coor_antes[2];
				}
				else
				{
#ifdef DEBUG
					cout << "\tAun asi, se mueve el nodo " << nodo << " de "
						<< start_values[0] << "  "
						<< start_values[1] << "  "
						<< start_values[2] << "  a  "
						<< m->nodo[nodo].x << "  "
						<< m->nodo[nodo].y << "  "
						<< m->nodo[nodo].z << endl;
#endif
				}
			}

		}
		// Quality of all mesh elements


		v_iter = mesh_.vertices_begin();
		//Mesh::VertexIter v_end = mesh_.vertices_end();

		for (; v_iter != v_end; ++v_iter)
		{
			int idx = boundary_idx[v_iter.handle().idx()];

			/*Point3 query(p[0], p[1], p[2]);
			Point_and_primitive_id pp = search_tree->closest_point_and_primitive(query);*/

			mesh_.set_point(v_iter, Mesh::Point(m->nodo[idx].x, m->nodo[idx].y, m->nodo[idx].z));

		}

		//Laplacian
#if 0
		

		calc_weights(&mesh_, eweight_);

		Laplacian_smooth(mesh_, eweight_, 1, 0.1);
		
		v_iter = mesh_.vertices_begin();
		
		for (; v_iter != v_end; ++v_iter)
		{
			int idx = boundary_idx[v_iter.handle().idx()];

			//projection here

			Point3 query(mesh_.point(v_iter)[0], mesh_.point(v_iter)[1], mesh_.point(v_iter)[2]);
			Point3 pp = search_tree->closest_point(query);

			//mesh_.set_point(v_iter, Mesh::Point(m->nodo[idx].x, m->nodo[idx].y, m->nodo[idx].z));
			/*m->nodo[idx].x = mesh_.point(v_iter)[0];
			m->nodo[idx].y = mesh_.point(v_iter)[1];
			m->nodo[idx].z = mesh_.point(v_iter)[2];*/
			m->nodo[idx].x = pp[0];
			m->nodo[idx].y = pp[1];
			m->nodo[idx].z = pp[2];


		}
#endif
		

		n = CalculaCalidades(&LaMalla, lista, 0.0);

		// If all elements are valid (not tangled) count this iteration
		// as an smooth iteration
		if (n == 0)
			iter_suav++;

		// Save file of element qualities
		if (fquality != NULL)
			GrabaFicheroCalidades(lista, LaMalla.Num_tetra, fquality);

		// Save statistics file
		if (f_estadisticas != NULL)
			GrabaEstadisticas(f_estadisticas, lista, LaMalla.Num_tetra);

		// Print statistics
		// Here median means average
		GrabaEstadisticas(stdout, lista, LaMalla.Num_tetra);

		// Quit the for loop if all smooth iterations have been done
		if (iter_suav == NumIterSuavizado)
			break;


		
		



	} // main loop

  // Mesh denormalization
	if (Normalized)
		DesnormalizarLaMalla(m);

	// Save untangled mesh to a file
	if (fsmoot != NULL)
		GrabaFicheroSalida(fsmoot, m);
	output_coord.clear();
	for (size_t i = 0; i < m->Num_nodos; i++)
	{
		output_coord.push_back(m->nodo[i].x);
		output_coord.push_back(m->nodo[i].y);
		output_coord.push_back(m->nodo[i].z);
	}


	//if used multiple times, search_tree can be preserved
	if (!search_tree)
		delete search_tree;
	//if ()

	if (f_estadisticas != NULL)
		fclose(f_estadisticas);

	//free memory
	LiberarMalla(m);
	free(lista);
	return true;
}



//bool SUS_standard(const std::vector<double> &nodes_coord, const std::vector<int> &ref_nodes, const std::vector<unsigned int> &elems, std::vector<double> &output_coord, int max_smooth_iter = 20, int max_iter = 40)
//{
//	//std::cout << "test SUS" << std::endl;
//	bool Normalized = true;
//	bool SaveQuality = false;
//	char *fnodes = NULL; // Nodes file (input)
//	char *felems = NULL; // Element file (input)
//	char *fstat = NULL; // Smooth statistics file (output)
//	char *fsmoot = NULL; // Smoothed mesh file (output). Only nodes are changed during smoothing
//	char *fquality = NULL; // File(s) quality prefix
//	char *ftetfile = NULL; // tet file (input)
//	char *fotetfile = NULL; // tet file (output)
//	int optchar;
//
//
//	/*Eigen::Matrix2d temp_m;
//	Eigen::Matrix2d m_I;
//	temp_m(0, 0) = 0;
//	temp_m(1, 0) = 0;
//	temp_m(0, 1) = 0;
//	temp_m(1, 1) = temp_m(1, 0) + temp_m(0, 1);
//	m_I = temp_m.inverse();*/
//
//
//
//
//	s.resize(MAX_CARAS_VISTAS * 9); sigma.resize(MAX_CARAS_VISTAS);
//
//	/*if ((fnodes == NULL) || (felems == NULL))
//	{
//		cout << "Missing file(s) name(s)" << endl;
//		uso();
//		exit(EXIT_FAILURE);
//	}*/
//
//	// Dimensionality of the problem 
//	const int dim = 3;
//
//	// starting values
//	double start_values[3] = { 0, 0, 0 };
//
//	// Find filehandle of standard output stream (used for intermediate minimization output)
//	FILE* fh = (FILE*)stdout;
//
//	// xpls will contain solution, gpls will contain
//	// gradient at solution after call to min_<f>.Minimize
//	dvec xpls(dim), gpls(dim);
//
//	// hpls provides the space for the hessian
//	dmat hpls(dim, dim);
//	// there are also these derived matrices
//	dmat cholesky(dim, dim), Hessian(dim, dim);
//
//	// fpls contains the value of the function at
//	// the solution given by xpls.
//	double fpls;
//
//	// Read input mesh
//	//int n = LeerMalla(fnodes, felems, &LaMalla);
//
//	std::vector<int> elems_int;
//	for (size_t i = 0; i < elems.size(); i++)
//	{
//		elems_int.push_back(elems[i]);
//	}
//	
//
//	int n = InitMalla(nodes_coord, ref_nodes, elems_int, &LaMalla);
//	//int n = InitMalla(nodes_deform, ref_nodes, elems, node_ori, node_normal_tag, &LaMalla);
//
//	if (n < 1)
//	{
//		cerr << "Error reading the mesh" << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	cout << "Mesh data:" << endl;
//	cout << "   nodes   : " << LaMalla.Num_nodos << endl;
//	cout << "   elements: " << LaMalla.Num_tetra << endl;
//
//	// Mesh normalization
//	if (Normalized)
//		NormalizarLaMalla(&LaMalla);
//
//	m = &LaMalla;
//
//	// Set typical magnitude of parameter values (off by one in this example)
//	double typical_size[3] = { 0., 0., 0. }; //
//	{
//		// Toma el tama�o tipico de cada variable como la media
//		// de los valores absolutos de las coordenadas de los nodos.
//		// Esto sirve para el escalado de las variables de la funcion
//		// objetivo
//		//
//		// Typical magnitude is obtained as the mean absolute value
//		for (int i = 0; i < m->Num_nodos; i++)
//		{
//			typical_size[0] += fabs(m->nodo[i].x);
//			typical_size[1] += fabs(m->nodo[i].y);
//			typical_size[2] += fabs(m->nodo[i].z);
//		}
//		typical_size[0] /= m->Num_nodos;
//		typical_size[1] /= m->Num_nodos;
//		typical_size[2] /= m->Num_nodos;
//	}
//	dvec typ(typical_size, typical_size + dim);
//
//	// Part 2: Load and minimimize the function. Function, first and
//	// second-order derivatives supplied
//	//
//	// Create function object
//	BivarFunctionGradientAndHessian FunctionGradientAndHessian;
//	// create Uncmin object
//	Uncmin<dvec, dmat, BivarFunctionGradientAndHessian> min_fgh(&FunctionGradientAndHessian);
//	// Set typical magnitude of function values at the minimum
//	min_fgh.SetScaleFunc(10);
//
//	// Set typical magnitude of argument values at function minimum
//	if (min_fgh.SetScaleArg(typ))
//	{
//		cout << endl << "Error in setting typical value" << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	// Diagnostic printout piped to stdout (fh, 1, 1)
//	// Only results piped to stdout (fh, 1, 0)
//	// No info piped to stdout (0, *, *)
//	min_fgh.SetPrint(0, 1, 0);
//
//	min_fgh.SetTolerances(-1, -1);
//
//	T_Malla *m = &LaMalla;
//
//	params.mascara = 7;   // minimiza en las 3 dimensiones
//	params.ndim = 3;
//	params.malla = m;
//	params.MargenEpsilon = EpsSafetyFactor; // "Safety" factor for epsilon
//
//	/* Memoria para guardar las calidades de los tetraedros */
//	double *lista = (double *)malloc(LaMalla.Num_tetra * sizeof(double));
//
//	// Crea un fichero para grabar las estadisticas
//	FILE *f_estadisticas = NULL;
//	if (fstat != NULL)
//	{
//		if ((f_estadisticas = fopen(fstat, "w")) == NULL)
//		{
//			cout << "Can't create statistics file " << fstat << endl;
//		}
//		else
//		{
//			// Statistics file header
//			fprintf(f_estadisticas, "   Min.         Median       Max.      # Invalid \n");
//			fprintf(f_estadisticas, "   Quality      Quality      Quality   Elements  \n");
//			fprintf(f_estadisticas, "--------------------------------------------------\n");
//		}
//	}
//
//	// Print statistics file header
//	fprintf(stdout, "\n             Mesh quality statistics\n");
//	fprintf(stdout, "   Min.         Median       Max.      # Invalid \n");
//	fprintf(stdout, "   Quality      Quality      Quality   Elements  \n");
//	fprintf(stdout, "--------------------------------------------------\n");
//
//	// 计算四面体的质量
//	n = CalculaCalidades(&LaMalla, lista, 0.0);
//
//	// 记录质量文件
//	if (fquality != NULL)
//		GrabaFicheroCalidades(lista, LaMalla.Num_tetra, fquality);
//
//	// 记录网格的统计信息
//	if (f_estadisticas != NULL)
//		GrabaEstadisticas(f_estadisticas, lista, LaMalla.Num_tetra);
//
//	// 打印统计信息
//	GrabaEstadisticas(stdout, lista, LaMalla.Num_tetra);
//
//	int iter_suav = -1;            // 平滑迭代计数器
//	// 在网格上迭代循环
//
//	MaxNumIter = max_iter;
//	NumIterSuavizado = max_smooth_iter;
//
//	for (int iteraciones = 1; iteraciones <= MaxNumIter; iteraciones++)
//	{
//		// 在节点中循环
//		for (nodo = 0; nodo < m->Num_nodos; nodo++)
//		{
//			double coor_antes[3];
//			double coor_despues[3];
//			params.nodo = nodo;
//			start_values[0] = params.var[0] = coor_antes[0] = m->nodo[nodo].x;
//			start_values[1] = params.var[1] = coor_antes[1] = m->nodo[nodo].y;
//			start_values[2] = params.var[2] = coor_antes[2] = m->nodo[nodo].z;
//
//			/*int normal_x, normal_y, normal_z;
//			normal_x = m->nodo[nodo].normal_x;
//			normal_y = m->nodo[nodo].normal_y;
//			normal_z = m->nodo[nodo].normal_z;*/
//
//			// 仅最小化参考编号为0的节点
//			if (m->nodo[nodo].nr != 0) continue;
//			deltaglobal = Calcula_Delta(&params, EpsilonEfectivo);
//#ifdef DEBUG
//			cout << "Minimizando el nodo " << nodo
//				<< " (nref " << LaMalla.nodo[nodo].nr << ")" << endl;
//#endif
//
//			//Minimize the function
//			int iMethod = 1;  // Line Search 
//			// Set Minimization method 
//			min_fgh.SetMethod(iMethod);
//
//			dvec start(start_values, start_values + dim);
//			// Saves the initial value of the objective function in f_ini
//			double f_ini;
//			{
//				BivarFunctionOnly v;
//				f_ini = v.f_to_minimize(start);
//				//f_ini = v.f_to_minimize(start, normal_x, normal_y, normal_z);
//			}
//#ifdef DEBUG
//			cout << "Inicialmente la funcion vale " << f_ini << endl;
//#endif
//			// Minimize function
//			min_fgh.Minimize(start, xpls, fpls, gpls, hpls);			
//			//min_fgh.Minimize(start, xpls, fpls, gpls, hpls, normal_x, normal_y, normal_z);
//
//			// Find status of function minimization.
//			// For a description of the values returned by
//			// GetMessage see the documentation at the beginning
//			// of Uncmin.h.
//			int msg = min_fgh.GetMessage();
//#ifdef DEBUG  
//			cout << endl << endl << "Message returned from Uncmin: " << msg << endl;
//			cout << endl << "Function minimum: " << endl << fpls << endl;
//			cout << endl << "Parameters: " << xpls << endl;
//			cout << "Gradient: " << gpls << endl;
//#endif
//			// Get lower triangle extract of hpls, with upper off-diagonal values set to zero.
//			cholesky = hpls;
//			for (int icolumn = 1; icolumn < cholesky.num_columns(); icolumn++)
//			{
//				for (int irow = 0; irow < icolumn; irow++)
//					cholesky[irow][icolumn] = 0.0;
//			}
//			Hessian = cholesky * transpose(cholesky);
//#ifdef DEBUG
//			cout << "Hessian:  " << Hessian << endl << endl;
//
//
//			cout << "Coords Antes: " << coor_antes[0] << "  "
//				<< coor_antes[1] << "  "
//				<< coor_antes[2] << endl;
//
//			cout << "Coords Despues: "
//				<< m->nodo[nodo].x << "  "
//				<< m->nodo[nodo].y << "  "
//				<< m->nodo[nodo].z << endl;
//
//			cout << "Coords xpls: " << xpls[0] << "  "
//				<< xpls[1] << "  "
//				<< xpls[2] << endl;
//#endif
//			// 如果找到了解决方案，请移动节点
//			if ((msg >= 1) && (msg <= 3))
//			{
//				//std::cout << "solution found" << std::endl;
//#ifdef DEBUG
//				cout << "main" << endl
//					<< "\tEncontrado el optimo (msg=" << msg << ")" << endl
//					<< "\tMovido el nodo " << nodo << " de "
//					<< start_values[0] << "  "
//					<< start_values[1] << "  "
//					<< start_values[2] << "  a  "
//					<< m->nodo[nodo].x << "  "
//					<< m->nodo[nodo].y << "  "
//					<< m->nodo[nodo].z << endl;
//#endif
//			}
//			else
//			{
//				//std::cout << "not found" << std::endl;
//				// 优化器尚未找到最佳解决方案：
//				// entonces, solo muevo el nodo si la solucion dada por
//				// el optimizador es valida y la nueva calidad es mejor
//				// que la de partida.
//				//
//				// Optimum not found. Actual node is moved to new
//				// (non-optimum) position iif:
//				// a) new position is valid and
//				// b) local mesh quality is better than original
//#ifdef DEBUG
//				cout << "\tNo encontrado el optimo (msg=" << msg
//					<< ")" << endl;
//#endif	      
//				if ((!(std::isfinite(m->nodo[nodo].x) && std::isfinite(m->nodo[nodo].y)
//					&& std::isfinite(m->nodo[nodo].z))) || (fpls > f_ini))
//				{
//#ifdef DEBUG
//					cout << "main" << endl
//						<< "\tNo se mueve el nodo " << nodo
//						<< " de su posicion inicial  "
//						<< coor_antes[0] << "  "
//						<< coor_antes[1] << "  "
//						<< coor_antes[2] << endl;
//#endif
//					std::cout << "no perfect solution found :" << std::endl;
//					m->nodo[nodo].x = coor_antes[0];
//					m->nodo[nodo].y = coor_antes[1];
//					m->nodo[nodo].z = coor_antes[2];
//				}
//				else
//				{
//#ifdef DEBUG
//					cout << "\tAun asi, se mueve el nodo " << nodo << " de "
//						<< start_values[0] << "  "
//						<< start_values[1] << "  "
//						<< start_values[2] << "  a  "
//						<< m->nodo[nodo].x << "  "
//						<< m->nodo[nodo].y << "  "
//						<< m->nodo[nodo].z << endl;
//#endif
//				}
//			}
//
//		}
//		// Quality of all mesh elements
//		n = CalculaCalidades(&LaMalla, lista, 0.0);
//
//		// If all elements are valid (not tangled) count this iteration
//		// as an smooth iteration
//		if (n == 0)
//			iter_suav++;
//
//		// Save file of element qualities
//		if (fquality != NULL)
//			GrabaFicheroCalidades(lista, LaMalla.Num_tetra, fquality);
//
//		// Save statistics file
//		if (f_estadisticas != NULL)
//			GrabaEstadisticas(f_estadisticas, lista, LaMalla.Num_tetra);
//
//		// Print statistics
//		// Here median means average
//		GrabaEstadisticas(stdout, lista, LaMalla.Num_tetra);
//
//		// Quit the for loop if all smooth iterations have been done
//		if (iter_suav == NumIterSuavizado)
//			break;
//
//	} // main loop
//
//  // Mesh denormalization
//	if (Normalized)
//		DesnormalizarLaMalla(m);
//
//	// Save untangled mesh to a file
//	if (fsmoot != NULL)
//		GrabaFicheroSalida(fsmoot, m);
//	output_coord.clear();
//	for (size_t i = 0; i < m->Num_nodos; i++)
//	{
//		output_coord.push_back(m->nodo[i].x);
//		output_coord.push_back(m->nodo[i].y);
//		output_coord.push_back(m->nodo[i].z);
//	}
//
//	if (f_estadisticas != NULL)
//		fclose(f_estadisticas);
//	LiberarMalla(m);
//	free(lista);
//	return true;
//}
//





