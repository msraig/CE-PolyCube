#include "UtilityDefinitions.h"

const int Utility_Definition::colorDisplay[5][3] = { {229,16,59}, {33,75,203}, {33,203,75},{127,0,255},{255,70,0} };

const char Utility_Definition::Color_String[5][512] = {"color: rgb(229,16,59)",
	"color: rgb(33,75,203)",
	"color: rgb(33,203,75)",
	"color: rgb(127,0,255)",
	"color: rgb(255,70,0)"};
const int Utility_Definition::PointColor[2][3] = { {128,0,42}, {128,170,0} };
const double Utility_Definition::AxisColor[3][3] = { {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0} };
const char Utility_Definition::PointColor_String[2][512] = {"color: rgb(128,0,42)",
	"color: rgb(128,170,0)"};

const double Utility_Definition::Permutation_Matrix[24][9] = {
	{1,0,0,0,1,0,0,0,1},{1,0,0,0,0,1,0,-1,0},{1,0,0,0,0,-1,0,1,0},{1,0,0,0,-1,0,0,0,-1},
	{-1,0,0,0,1,0,0,0,-1},{-1,0,0,0,0,-1,0,-1,0},{-1,0,0,0,0,1,0,1,0},{-1,0,0,0,-1,0,0,0,1},
	{0,1,0,-1,0,0,0,0,1},{0,1,0,0,0,1,1,0,0},{0,1,0,0,0,-1,-1,0,0},{0,1,0,1,0,0,0,0,-1},
	{0,-1,0,1,0,0,0,0,1},{0,-1,0,0,0,-1,1,0,0},{0,-1,0,0,0,1,-1,0,0},{0,-1,0,-1,0,0,0,0,-1},
	{0,0,1,1,0,0,0,1,0},{0,0,1,0,1,0,-1,0,0},{0,0,1,0,-1,0,1,0,0},{0,0,1,-1,0,0,0,-1,0},
	{0,0,-1,-1,0,0,0,1,0},{0,0,-1,0,1,0,1,0,0},{0,0,-1,0,-1,0,-1,0,0},{0,0,-1,1,0,0,0,-1,0}};

const double Utility_Definition::max_length_slider = 100;
const double Utility_Definition::max_angle_slider  = 100;
const double Utility_Definition::max_cuboid_sheet_slider = 1000.0;

const char Utility_Definition::tensor_name[11][128] = {"Identity", "[10 1 1]", "X^2+y^2+z^4", "Exp(x^2+y^2+z^2)", "Spherical Shock", "Cylinder Shock", "Planar Shock","Interpolation" ,"Sine", "Sink", "Test"};

const char Utility_Definition::energy_name[7][128] = {"LCOT_SRC","LCOT_SRC2", "LCOT_SRC_OBO","LCOT_SRC_OBO2", "LCOR_SRC_LP", "LCOT_SRC_G", "LCOT_SRC_G_OBO" };

const char Utility_Definition::opt_tet_name[5][128] = { "Opt F Condition", "Opt F Condition OMP", "Opt F Condition CD", "Opt Trace Det", "Opt OGT" };

const double Utility_Definition::negative_volume_check_th = 1.0e-15;


