/* $Id: sus.h 21 2010-02-24 15:23:06Z erodriguez $ */
#ifndef __SUS_H
#define __SUS_H

/* -------------------------------------------------------------
SUS code
Simultaneous Untangling and Smoothing
Programa de desenredado y suavizado de mallas 3D

sus.h

Copyright (C) 2010 
University of Las Palmas de Gran Canaria - ULPGC
Institute of Intelligent Systems and 
Numerical Applications in Engineering - SIANI
Jose Maria Escobar Sanchez
Eduardo Rodriguez Barrera
Rafael Montenegro Armas
Gustavo Montero Garcia
Jose Maria Gonzalez Yuste
-------------------------------------------------------------*/


/* raiz cuadrada de 2 */
#ifndef rd2
#define rd2 (1.414213562373095633972752693807)
#endif
/* raiz cuadrada de 3 */
#ifndef rd3
#define rd3 (1.732050807568877637265813973499)
#endif
/* raiz cuadrada de 6 */
#ifndef rd6
#define rd6 (2.449489742783178098197284074706)
#endif

#ifndef sigma_x
#define sigma_x(y1, y2, y3, z1, z2, z3) (rd2 * (y3*(z2-z1) + y2*(z1-z3) + y1*(z3-z2)))
#endif
#ifndef sigma_y
#define sigma_y(x1, x2, x3, z1, z2, z3) (rd2 * (x3*(z1-z2) + x1*(z2-z3) + x2*(z3-z1)))
#endif
#ifndef sigma_z
#define sigma_z(x1, x2, x3, y1, y2, y3) (rd2 * (x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2)))
#endif
#ifndef PES_gS
#define PES_gS(x, x1, x2, x3) (0.5*(3*x -x1 -x2 -x3))
#endif


struct
S_Nodo
{
  double x,y,z;      /* Coordenadas del nodo */
  double o_x, o_y, o_z; //original coordinates
  int normal_x, normal_y, normal_z; //normal tag of three coordinate
  
  int nr;            /* numero de referencia 参考编号 */
  int NCarasVistas;  /* numero de caras de las que el nodo forma parte 节点所属的面数 */
  int *CarasVistas[3]; /* vector de nodos que forman cada cara构成每张面的节点矢量 (3xNCarasVistas)*/
};
typedef struct S_Nodo T_Nodo;

struct
S_Tetra
{
  int nodo[4];       /* nodos que forman el tetraedro */
  double q;          /* 测量四面体的质量 */
  double IS[3][3];
}; 
typedef struct S_Tetra T_Tetra;

struct               /* Estructura para guadar datos de los nodos */
S_Malla
{
  int Num_nodos;     /* 网格节点数 */
  int Num_tetra;     /* 网格的四面体数量 */
  T_Nodo *nodo;      /* 形成网格的节点的矢量 (Num_nodos)*/
  T_Tetra *tetra;    /* 矢量四面体形成网格 (Num_tetra) */
  /// Informacion para el normalizado
  double cx, cy, cz;    // 网格的重心
  double max_radio;     // 环绕球半径
};
typedef struct S_Malla T_Malla;

struct 
S_ParProg 
{
  char fich_coor[512]; /* fichero de coordenadas de los nodos de la malla */
  char fich_elem[512]; /* fichero de */
  char fich_cone[512];
  char fich_calidad[512];
  double norma;
  int verbose;
  char calidad[32]; /* medida de calidad a usar (Freitag o Knupp) */
  int NitDes; /* Numero maximo de iteraciones del desenrredado */
  int NitSuav; /* Numero maximo de iteraciones del suavizado */
  char fich_salida[512]; /* Fichero donde grabar la malla suavizada*/
  double margen_delta;  /* exponente del margen sobre el delta calculado */
  double tol_suav; /* tolerancia del algoritmo de optimizacion en */
				   /* la fase de suavizado */
  double tol_dese; /* tolerancia del algoritmo de optimizacion en */
				   /* la fase de desenrredado */
  double delta;    /* delta de la funcion modificada */
  double N_fich;   /* Numero de ficheros de calidades a escribir */
  char fich_estadis[512]; /* fichero de estadisticas */
  double paso_ini; /* Tamaño del paso inicial para el metodo de minimizacion */
				   /* de las librerias gsl (steepest descent, bfgs, ...) */
  int N_pasos_mini; /* Numero de pasos para el algoritmo de minimizacion */
  char fich_refina[512]; /* Fichero de salida con los tetraedros a refinar */
  int N_refina; /* Numero de tetraedros a refinar (los N_refina peores) */
  int margen_eps; /* Exponente de seguridad para el calculo de epsilon a */
  /* partir del epsilon de la maquina. Interviene en el calculo de la delta */
  double MargenEpsilon; /* 10^margen_eps */
};

typedef struct S_ParProg T_ParProg;

typedef double (*T_Func_Knupp)(const T_Malla *m, double s[3][3], 
				   const double det, const double Delta,
				   double *nfr);

struct mis_parametros {
  T_Malla *malla;
  int nodo;
  int ndim; /* 要优化的函数的维数 */
  double var[3]; /* x，y，z节点的原始 */
  int mascara; /* mascara de bits indicadora de que dimensiones 位掩码指示什么尺寸 */
  /* se permite mover en la minimizacion它被允许移动到最小化 (ver mi_OptimizaNodo) */
  T_ParProg *ParProg;
  double Delta; /* La delta para funcion de Knupp Knupp函数的delta*/
  double delta_calculada; /* La delta calculada */
  int salir;   /* 1 -> indica que hay que salir de las iteraciones de Newton 表明你必须摆脱牛顿的迭代*/
  double MargenEpsilon;
};


#endif
