#ifndef UTILITY_DEFINITIONS_H
#define UTILITY_DEFINITIONS_H

#include <cmath>

namespace Utility_Definition
{	
	enum COLOR_Definition{RED=0,BLUE,GREEN,PURPLE,ORANGE,N_COLOR};
	extern const int colorDisplay[5][3];
	extern const char Color_String[5][512];
	extern const int PointColor[2][3];
	extern const double AxisColor[3][3];
	extern const char PointColor_String[2][512];
	
	extern const double Permutation_Matrix[24][9];

	struct Quaternion
	{
		Quaternion()
		{
			m_x = m_y = m_z = 0.0; m_w = 1.0;
		}
		void CreateMatrix(double *pMatrix)
		{
			// Make sure the matrix has allocated memory to store the rotation data
			if(!pMatrix) return;
			// First row
			pMatrix[ 0] = 1.0 - 2.0 * ( m_y * m_y + m_z * m_z ); 
			pMatrix[ 1] = 2.0 * (m_x * m_y + m_z * m_w);
			pMatrix[ 2] = 2.0 * (m_x * m_z - m_y * m_w);
			pMatrix[ 3] = 0.0;  
			// Second row
			pMatrix[ 4] = 2.0 * ( m_x * m_y - m_z * m_w );  
			pMatrix[ 5] = 1.0 - 2.0 * ( m_x * m_x + m_z * m_z ); 
			pMatrix[ 6] = 2.0 * (m_z * m_y + m_x * m_w );  
			pMatrix[ 7] = 0.0;  
			// Third row
			pMatrix[ 8] = 2.0 * ( m_x * m_z + m_y * m_w );
			pMatrix[ 9] = 2.0 * ( m_y * m_z - m_x * m_w );
			pMatrix[10] = 1.0 - 2.0 * ( m_x * m_x + m_y * m_y );  
			pMatrix[11] = 0.0;  
			// Fourth row
			pMatrix[12] = 0;
			pMatrix[13] = 0;
			pMatrix[14] = 0;
			pMatrix[15] = 1.0;
			// Now pMatrix[] is a 4x4 homogeneous matrix that can be applied to an OpenGL Matrix
		}
		void CreateFromAxisAngle(double x, double y, double z, double angle)
		{
			// Here we calculate the sin( theta / 2) once for optimization
			double result = std::sin( angle / 2.0 );
			// Calculate the w value by cos( theta / 2 )
			m_w = std::cos( angle / 2.0 );
			// Calculate the x, y and z of the quaternion
			m_x = x * result;
			m_y = y * result;
			m_z = z * result;
		}
		void CreateFromXYZW(double x, double y, double z, double w)
		{
			m_x = x; m_y = y; m_z = z; m_w = w;
		}
		Quaternion operator *(Quaternion q)
		{
			Quaternion r;
			r.m_w = m_w*q.m_w - m_x*q.m_x - m_y*q.m_y - m_z*q.m_z;
			r.m_x = m_w*q.m_x + m_x*q.m_w + m_y*q.m_z - m_z*q.m_y;
			r.m_y = m_w*q.m_y + m_y*q.m_w + m_z*q.m_x - m_x*q.m_z;
			r.m_z = m_w*q.m_z + m_z*q.m_w + m_x*q.m_y - m_y*q.m_x;

			return(r);
		}
		void conj()
		{
			m_z = -m_z;
			m_y = -m_y;
			m_x = -m_x;
		}
		double X(){return m_x;}
		double Y(){return m_y;}
		double Z(){return m_z;}
		double W(){return m_w;}

	private:
		double m_w;
		double m_z;
		double m_y;
		double m_x;
	};

	extern const double max_length_slider;
	extern const double max_angle_slider;
	extern const double max_cuboid_sheet_slider;

	enum TENSOR_NAME
	{ T_IDENTITY = 0, T_10_1_1, T_X2_Y2_Z4, T_EXP, T_SPHERICAL_SHOCK, T_CYLINDER_SHOCK, T_PLANAR_SHOCK, T_INTERPOLATION, T_SINE, T_SINK, T_TEST };
	extern const char tensor_name[11][128];
	enum ENERGY_NAME
	{
		ENERGY_LCOT = 0, ENERGY_LCOT2, ENERGY_LCOT_OBO, ENERGY_LCOT_OBO2, ENERGY_LCOT_LP, ENERGY_LCOT_G, ENERGY_LCOT_G_OBO, OPT_F_CONDITION, OPT_F_CONDITION_OMP, OPT_F_CONDITION_CD, OPT_TRACE_DET, OPT_OGT,
	};
	extern const char energy_name[7][128];

	extern const char opt_tet_name[5][128];

	extern const double negative_volume_check_th;
}


#endif