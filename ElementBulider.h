#pragma once
#include <vector>
#include "ElementBulider.h"


struct MeshTriangle{
	int m_nPoint;
	int m_nTriangles;
	int* m_Triangles;
	float* m_Points;
	float* m_Normals;
};
struct Pentahedron{
	float m_Node[6];
};
typedef std::vector<Pentahedron> PentahedronVector;

class ElementBulider
{
public:
	ElementBulider(void);
	~ElementBulider(void);
	void GenerateElement(MeshTriangle MeshTriangleOne,MeshTriangle MeshTriangleTwo);
	int GetNumOfElement();
	int GetNumOfVertice();
	int* GetElement();
	float* GetVerticeOfElement();

protected:
	bool RayIntersectTriangle(float orig[3],float direction[3],float vertice1[3],float vertice1[3],float vertice2[3],float* t,float * u,float* v);

public:
	int m_NumOfElement;
	int m_NumOfVertice;
	int* m_piElementPentahedronIndices;
	float* m_ppt3dVertices;

protected:
	MeshTriangle* m_MeshTriangleOne;
	MeshTriangle* m_MeshTriangleTwo;
};
