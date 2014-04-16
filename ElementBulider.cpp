#include "stdafx.h"
#include ".\elementbulider.h"
#define EPSILON 0.000001
#define CROSS(dest,v1,v2)\
	dest[0]=v1[1]*v2[2]-v1[2]*v2[1];\
	dest[1]=v1[2]*v2[0]-v1[0]*v2[2];\
	dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

ElementBulider::ElementBulider(void)
{
}

ElementBulider::~ElementBulider(void)
{
    
}
ElementBulider::GenerateElement(MeshTriangle MeshTriangleOne,MeshTriangle MeshTriangleTwo);
ElementBulider::RayIntersectTriangle(float orig[3],float direction[3],float vertice1[3],float vertice1[3],float vertice2[3],float* t,float * u,float* v);
int* ElementBulider::GetNumOfElement()
{
	return m_piElementPentahedronIndices;
}
int ElementBulider::GetNumOfVertice()
{
	return m_NumOfVertice;
}
int ElementBulider::GetNumOfElement()
{
	return m_NumOfElement;
}
float* ElementBulider::GetVerticeOfElement()
{
	return m_ppt3dVertices;
}
bool ElementBulider::RayIntersectTriangle(float orig[3],float direction[3],float vertice1[3],float vertice1[3],float vertice2[3],float* t,float * u,float* v)
{

}