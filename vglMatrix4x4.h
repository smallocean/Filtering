
#ifndef _vglMatrix4x4_h
#define _vglMatrix4x4_h


class  vglMatrix4x4  
{

public:
	//! Default constructor
	vglMatrix4x4();

	//! Default destructor
	virtual ~vglMatrix4x4();

	// set all element to zero
	void Zero();	

	// set to identity
	void Identity();

	// Transpose the matrix and put it into out. 
	static void Transpose(vglMatrix4x4 *in, vglMatrix4x4 *out); 
	void Transpose() { vglMatrix4x4::Transpose(this,this); };
  
	// Matrix Inversion (adapted from Richard Carling in "Graphics Gems," 
	// Academic Press, 1990).	
	static void Invert(vglMatrix4x4 *in, vglMatrix4x4 *out); 
	void Invert() { vglMatrix4x4::Invert(this,this); }; 
   
	// Compute adjoint of the matrix and put it into out.
	static void Adjoint(vglMatrix4x4 *in, vglMatrix4x4 *out); 

	// Compute the determinant of the matrix and return it.
	static float Determinant(vglMatrix4x4 *mat);
	float Determinant(){ return vglMatrix4x4::Determinant(this); };

	void MultiplyVector(float in[4], float out[4]);
	
	// Multiply a homogeneous coordinate by this matrix, i.e. out = A*in.
	// The in[4] and out[4] can be the same array.


//	void MultiplyPoint(vglPoint3f in, vglPoint3f &out);

	// Multiplies matrices a and b and stores the result in c.
	static void Multiply(vglMatrix4x4 *a, vglMatrix4x4 *b, vglMatrix4x4 *c);

	// Create a translation matrix 
	void Translated(float x, float y, float z);

	// Create a rotation matrix 
	void Rotated(float angle, float x, float y, float z);
	
	// Create a scale matrix (i.e. set the diagonal elements to x, y, z)
	// and concatenate it with the current transformation according to
	// PreMultiply semantics.
	void Scaled(float x, float y, float z);

	//operator "="
	vglMatrix4x4& operator = (const vglMatrix4x4 &src);
	
private:
	// Calculate the determinant of a 3x3 matrix in the form:
	//     | a1,  b1,  c1 |
	//     | a2,  b2,  c2 |
	//     | a3,  b3,  c3 |
	static float Determinant3x3(float a1,float a2,float a3,
						  float b1,float b2,float b3,
						  float c1,float c2,float c3);

public: 
	float m_dElement[4][4];
};

#endif //_vglMatrix4x4_h