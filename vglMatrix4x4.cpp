
// vgl include
#include "stdafx.h"
#include "vglMatrix4x4.h"

// stl include
#include <math.h>

vglMatrix4x4::vglMatrix4x4()
{
	this ->Identity();
}

vglMatrix4x4::~vglMatrix4x4()
{

}

void vglMatrix4x4::Zero()
{
	int i,j;
	for (i = 0; i < 4; i++){
		for (j = 0; j < 4; j++){
			m_dElement[i][j] = 0.0;
      }
    }
}

void vglMatrix4x4::Identity()
{
	m_dElement[0][0] = m_dElement[1][1] = m_dElement[2][2] = m_dElement[3][3] = 1.0;
	m_dElement[0][1] = m_dElement[0][2] = m_dElement[0][3] = 
	m_dElement[1][0] = m_dElement[1][2] = m_dElement[1][3] =
	m_dElement[2][0] = m_dElement[2][1] = m_dElement[2][3] =
	m_dElement[3][0] = m_dElement[3][1] = m_dElement[3][2] =0.0;
}	

void vglMatrix4x4::Transpose(vglMatrix4x4 *in, vglMatrix4x4 *out)
{
	int i, j;

	for (i=0; i<4; i++){
		for(j=i; j<4; j++){
			out->m_dElement[j][i] = in->m_dElement[i][j];
      }
    }
}

void vglMatrix4x4::Invert(vglMatrix4x4 *in,vglMatrix4x4 *out)
{
	int i, j;
	float det;

	// calculate the 4x4 determinent
	// if the determinent is zero, 
	// then the inverse matrix is not unique.

	det = vglMatrix4x4::Determinant(in);
  
	if ( det == 0.0 ) 
    {
		out->Zero();
		return;
    }

	// calculate the adjoint matrix
	vglMatrix4x4::Adjoint(in,out);

	// scale the adjoint matrix to get the inverse
	for (i=0; i<4; i++)
    {
		for(j=0; j<4; j++)
		{
			out->m_dElement[i][j] /= det;
		}
    }

}

void vglMatrix4x4::Adjoint(vglMatrix4x4 *in,vglMatrix4x4 *out)
{
	float a1, a2, a3, a4, b1, b2, b3, b4;
	float c1, c2, c3, c4, d1, d2, d3, d4;
	
	// assign to individual variable names to aid
	// selecting correct values

	a1 = in->m_dElement[0][0]; b1 = in->m_dElement[0][1]; 
	c1 = in->m_dElement[0][2]; d1 = in->m_dElement[0][3];

	a2 = in->m_dElement[1][0]; b2 = in->m_dElement[1][1]; 
	c2 = in->m_dElement[1][2]; d2 = in->m_dElement[1][3];

	a3 = in->m_dElement[2][0]; b3 = in->m_dElement[2][1];
	c3 = in->m_dElement[2][2]; d3 = in->m_dElement[2][3];

	a4 = in->m_dElement[3][0]; b4 = in->m_dElement[3][1]; 
	c4 = in->m_dElement[3][2]; d4 = in->m_dElement[3][3];

	// row column labeling reversed since we transpose rows & columns

	out->m_dElement[0][0]  =  Determinant3x3( b2, b3, b4, c2, c3, c4, d2, d3, d4);
	out->m_dElement[1][0]  = -Determinant3x3( a2, a3, a4, c2, c3, c4, d2, d3, d4);
	out->m_dElement[2][0]  =  Determinant3x3( a2, a3, a4, b2, b3, b4, d2, d3, d4);
	out->m_dElement[3][0]  = -Determinant3x3( a2, a3, a4, b2, b3, b4, c2, c3, c4);

	out->m_dElement[0][1]  = -Determinant3x3( b1, b3, b4, c1, c3, c4, d1, d3, d4);
	out->m_dElement[1][1]  =  Determinant3x3( a1, a3, a4, c1, c3, c4, d1, d3, d4);
	out->m_dElement[2][1]  = -Determinant3x3( a1, a3, a4, b1, b3, b4, d1, d3, d4);
	out->m_dElement[3][1]  =  Determinant3x3( a1, a3, a4, b1, b3, b4, c1, c3, c4);
        
	out->m_dElement[0][2]  =  Determinant3x3( b1, b2, b4, c1, c2, c4, d1, d2, d4);
	out->m_dElement[1][2]  = -Determinant3x3( a1, a2, a4, c1, c2, c4, d1, d2, d4);
	out->m_dElement[2][2]  =  Determinant3x3( a1, a2, a4, b1, b2, b4, d1, d2, d4);
	out->m_dElement[3][2]  = -Determinant3x3( a1, a2, a4, b1, b2, b4, c1, c2, c4);
        
	out->m_dElement[0][3]  = -Determinant3x3( b1, b2, b3, c1, c2, c3, d1, d2, d3);
	out->m_dElement[1][3]  =  Determinant3x3( a1, a2, a3, c1, c2, c3, d1, d2, d3);
	out->m_dElement[2][3]  = -Determinant3x3( a1, a2, a3, b1, b2, b3, d1, d2, d3);
	out->m_dElement[3][3]  =  Determinant3x3( a1, a2, a3, b1, b2, b3, c1, c2, c3);

}

float vglMatrix4x4::Determinant(vglMatrix4x4 *mat)
{
	float a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4;

	// assign to individual variable names to aid selecting
	//  correct elements
	
	a1 = mat->m_dElement[0][0]; b1 = mat->m_dElement[0][1]; 
	c1 = mat->m_dElement[0][2]; d1 = mat->m_dElement[0][3];

	a2 = mat->m_dElement[1][0]; b2 = mat->m_dElement[1][1]; 
	c2 = mat->m_dElement[1][2]; d2 = mat->m_dElement[1][3];

	a3 = mat->m_dElement[2][0]; b3 = mat->m_dElement[2][1]; 
	c3 = mat->m_dElement[2][2]; d3 = mat->m_dElement[2][3];

	a4 = mat->m_dElement[3][0]; b4 = mat->m_dElement[3][1]; 
	c4 = mat->m_dElement[3][2]; d4 = mat->m_dElement[3][3];

  return a1 * Determinant3x3( b2, b3, b4, c2, c3, c4, d2, d3, d4)
       - b1 * Determinant3x3( a2, a3, a4, c2, c3, c4, d2, d3, d4)
       + c1 * Determinant3x3( a2, a3, a4, b2, b3, b4, d2, d3, d4)
       - d1 * Determinant3x3( a2, a3, a4, b2, b3, b4, c2, c3, c4);
}

float vglMatrix4x4::Determinant3x3(float a1,float a2,float a3, 
								   float b1,float b2,float b3, 
								   float c1,float c2,float c3)
{
	return a1*(b2*c3 - b3*c2) - a2*(b1*c3-b3*c1) + a3*(b1*c2 - b2*c1);
}

void vglMatrix4x4::MultiplyVector(float in[4], float out[4])
{
	
	out[0] = in[0] * m_dElement[0][0] + in[1] * m_dElement[0][1] + in[2] * m_dElement[0][2] + in[3] * m_dElement[0][3];
	out[1] = in[0] * m_dElement[1][0] + in[1] * m_dElement[1][1] + in[2] * m_dElement[1][2] + in[3] * m_dElement[1][3];
	out[2] = in[0] * m_dElement[2][0] + in[1] * m_dElement[2][1] + in[2] * m_dElement[2][2] + in[3] * m_dElement[2][3];
	out[3] = in[0] * m_dElement[3][0] + in[1] * m_dElement[3][1] + in[2] * m_dElement[3][2] + in[3] * m_dElement[3][3];
}
/*void vglMatrix4x4::MultiplyPoint(vglPoint3f in, vglPoint3f &out)
{
	float tmp[4];

	tmp[0] = in.x() * m_dElement[0][0] + in.y() * m_dElement[0][1] + in.z() * m_dElement[0][2] + m_dElement[0][3];
	tmp[1] = in.x() * m_dElement[1][0] + in.y() * m_dElement[1][1] + in.z() * m_dElement[1][2] + m_dElement[1][3];
	tmp[2] = in.x() * m_dElement[2][0] + in.y() * m_dElement[2][1] + in.z() * m_dElement[2][2] + m_dElement[2][3];
	tmp[3] = in.x() * m_dElement[3][0] + in.y() * m_dElement[3][1] + in.z() * m_dElement[3][2] + m_dElement[3][3];

	if(tmp[3] == 0) 
	{
		out.x(0);out.y(0);out.z(0);
		return;
	}

	if(tmp[3]!=1){
		tmp[0] /= tmp[3];
		tmp[1] /= tmp[3];
		tmp[2] /= tmp[3];
	}

	//copy to out
	out.x(tmp[0]);out.y(tmp[1]);out.z(tmp[2]);
}
*/

void vglMatrix4x4::Multiply(vglMatrix4x4 *a, vglMatrix4x4 *b, vglMatrix4x4 *c)
{
	int i, k;
	float Accum[4][4];

	for (i = 0; i < 4; i++) 
    {
		for (k = 0; k < 4; k++) 
		{
			Accum[i][k] = a->m_dElement[i][0] * b->m_dElement[0][k] +
                    a->m_dElement[i][1] * b->m_dElement[1][k] +
                    a->m_dElement[i][2] * b->m_dElement[2][k] +
                    a->m_dElement[i][3] * b->m_dElement[3][k];
		}
    }

	// Copy to final dest
	for (i = 0; i < 4; i++)
    {
		c->m_dElement[i][0] = Accum[i][0];
		c->m_dElement[i][1] = Accum[i][1];
		c->m_dElement[i][2] = Accum[i][2];
		c->m_dElement[i][3] = Accum[i][3];
    }

}

vglMatrix4x4& vglMatrix4x4::operator =(const vglMatrix4x4 &src)
{
	int i,j;
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
			m_dElement[i][j]=src.m_dElement[i][j];
	}
	return *this;
}

void vglMatrix4x4::Translated(float x, float y, float z)
{
	if (x == 0.0 && y == 0.0 && z == 0.0) 
    {
		return;
    }

	Identity();
	
	m_dElement[0][3] = x;
	m_dElement[1][3] = y;
	m_dElement[2][3] = z;
}

void vglMatrix4x4::Rotated(float angle, float x, float y, float z)
{
	if (angle == 0.0 || (x == 0.0 && y == 0.0 && z == 0.0)) 
    {
		return;
    }

	// convert to radians
	float angleToRadians = 3.1415926/180;
	angle = angle*angleToRadians;

	// make a normalized quaternion
	float w = cos(0.5*angle);
	float f = (sin(0.5*angle)/sqrt(x*x+y*y+z*z));
	x *= f;
	y *= f;
	z *= f;


	// convert the quaternion to a matrix
	Identity();

	float ww = w*w;
	float wx = w*x;
	float wy = w*y;
	float wz = w*z;

	float xx = x*x;
	float yy = y*y;
	float zz = z*z;

	float xy = x*y;
	float xz = x*z;
	float yz = y*z;

	float s = ww - xx - yy - zz;

	m_dElement[0][0] = xx*2 + s;
	m_dElement[1][0] = (xy + wz)*2;
	m_dElement[2][0] = (xz - wy)*2;

	m_dElement[0][1] = (xy - wz)*2;
	m_dElement[1][1] = yy*2 + s;
	m_dElement[2][1] = (yz + wx)*2;

	m_dElement[0][2] = (xz + wy)*2;
	m_dElement[1][2] = (yz - wx)*2;
	m_dElement[2][2] = zz*2 + s;
}


void vglMatrix4x4::Scaled(float x, float y, float z)
{
	if (x == 1.0 && y == 1.0 && z == 1.0) 
    {
		return;
    }

	Identity();
	
	m_dElement[0][0] = x;
	m_dElement[1][1] = y;
	m_dElement[2][2] = z;
}
