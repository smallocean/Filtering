//////////////////////////////////////////////////////////////////////////
///生成一个四面体网格 试验数据  具有线性颜色变化纹理
#include <math.h>
#include "NormalStruct.h"

//struct KnearestField{
//	unsigned int m_numOfNearest;//邻域的数量，必须小于k
//	bool m_IfBoundary;//是否是边界 如果是为 true，否则为false
//	std::vector<int> m_nearest;//存放邻域点
//};// 保存点的邻域的结构
//
//struct Triangle{
//	int vertexIndex[3];
//};
//
//struct Color {
//	float pointColor[3];
//};
//
//struct POINTVECTOR3D {
//	float pointVector[3];
//};
//struct Triangle{
//	int vertexIndex[3];
//};
//
//struct Color {
//	float pointColor[3];
//};

//struct Point{
//	float vertexPosition[3];
//};

//Point m_originalVertexs[4];
//m_originalVertexs[0].vertexPosition[0]=-30;
//m_originalVertexs[0].vertexPosition[1]=-20;
//m_originalVertexs[0].vertexPosition[2]=0;
//m_originalVertexs[1].vertexPosition[0]=-30;
//m_originalVertexs[1].vertexPosition[1]=20;
//m_originalVertexs[1].vertexPosition[2]=0;
//m_originalVertexs[2].vertexPosition[0]=0;
//m_originalVertexs[2].vertexPosition[1]=40;
//m_originalVertexs[2].vertexPosition[2]=0;
//m_originalVertexs[3].vertexPosition[0]=0;
//m_originalVertexs[3].vertexPosition[1]=0;
//m_originalVertexs[3].vertexPosition[2]=40*sqrt(2);
//
// Triangle m_originalTriangles[4];
//m_originalTriangles[0].vertexIndex[0]=0;
//m_originalTriangles[0].vertexIndex[1]=1;
//m_originalTriangles[0].vertexIndex[2]=2;
//m_originalTriangles[1].vertexIndex[0]=1;
//m_originalTriangles[1].vertexIndex[1]=3;
//m_originalTriangles[1].vertexIndex[2]=2;
//m_originalTriangles[2].vertexIndex[0]=3;
//m_originalTriangles[2].vertexIndex[1]=1;
//m_originalTriangles[2].vertexIndex[2]=0;
//m_originalTriangles[3].vertexIndex[0]=2;
//m_originalTriangles[3].vertexIndex[1]=3;
//m_originalTriangles[3].vertexIndex[2]=0;
//////////////////////////////////////////////////////////////////////////
// 由一个三角形生成4个三角形

void subdivisionSingleTriangle(Point originalPoint[3],Triangle resultTriangle[4],Point resultPoint[6]);


//////////////////////////////////////////////////////////////////////////
// 循环产生网格

void GenerateMesh(int numOfVertex,int numOfMeshes,float* vertex,int* triangleMesh);

void CalculateBoundaryColorOfVertexs(int numOfPoints,float* pointSets,float* normals);
