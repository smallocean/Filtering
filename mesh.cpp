
#include "stdafx.h"
#include "mesh.h"
#include "Matrix.h"

void subdivisionSingleTriangle(Point originalPoint[3],Triangle resultTriangle[4],Point resultPoint[6])
{

	resultTriangle[0].vertexIndex[0]=0;
	resultTriangle[0].vertexIndex[1]=3;
	resultTriangle[0].vertexIndex[2]=5;
	resultTriangle[1].vertexIndex[0]=3;
	resultTriangle[1].vertexIndex[1]=1;
	resultTriangle[1].vertexIndex[2]=4;
	resultTriangle[2].vertexIndex[0]=3;
	resultTriangle[2].vertexIndex[1]=4;
	resultTriangle[2].vertexIndex[2]=5;
	resultTriangle[3].vertexIndex[0]=5;
	resultTriangle[3].vertexIndex[1]=4;
	resultTriangle[3].vertexIndex[2]=2;

	resultPoint[0].vertexPosition[0]=originalPoint[0].vertexPosition[0];
	resultPoint[0].vertexPosition[1]=originalPoint[0].vertexPosition[1];
	resultPoint[0].vertexPosition[2]=originalPoint[0].vertexPosition[2];
	resultPoint[1].vertexPosition[0]=originalPoint[1].vertexPosition[0];
	resultPoint[1].vertexPosition[1]=originalPoint[1].vertexPosition[1];
	resultPoint[1].vertexPosition[2]=originalPoint[1].vertexPosition[2];
	resultPoint[2].vertexPosition[0]=originalPoint[2].vertexPosition[0];
	resultPoint[2].vertexPosition[1]=originalPoint[2].vertexPosition[1];
	resultPoint[2].vertexPosition[2]=originalPoint[2].vertexPosition[2];
	resultPoint[3].vertexPosition[0]=(originalPoint[0].vertexPosition[0]+originalPoint[1].vertexPosition[0])/2;
	resultPoint[3].vertexPosition[1]=(originalPoint[0].vertexPosition[1]+originalPoint[1].vertexPosition[1])/2;
	resultPoint[3].vertexPosition[2]=(originalPoint[0].vertexPosition[2]+originalPoint[1].vertexPosition[2])/2;
	resultPoint[4].vertexPosition[0]=(originalPoint[2].vertexPosition[0]+originalPoint[1].vertexPosition[0])/2;
	resultPoint[4].vertexPosition[1]=(originalPoint[2].vertexPosition[1]+originalPoint[1].vertexPosition[1])/2;
	resultPoint[4].vertexPosition[2]=(originalPoint[2].vertexPosition[2]+originalPoint[1].vertexPosition[2])/2;
	resultPoint[5].vertexPosition[0]=(originalPoint[0].vertexPosition[0]+originalPoint[2].vertexPosition[0])/2;
	resultPoint[5].vertexPosition[1]=(originalPoint[0].vertexPosition[1]+originalPoint[2].vertexPosition[1])/2;
	resultPoint[5].vertexPosition[2]=(originalPoint[0].vertexPosition[2]+originalPoint[2].vertexPosition[2])/2;
	return;
}

//////////////////////////////////////////////////////////////////////////
// 循环产生网格
/*
 *	输入参数：
 *    int numOfMeshes  
 *    int numOfVertex
 *    float* vertex;
 *    int*   triangleMesh;
 */

void GenerateMesh(int numOfVertex,int numOfMeshes,float* vertex,int* triangleMesh)
{
	int numOfIteration;//交互次数numOfIteration＝7;
	int m_numOfOriginalMeshs;//原始三角片的数量（循环前）
	int m_numOfResultMeshs;//结果三角片的数量（循环后）
	int m_numOfOriginalPoints;//结果点的数量（循环前）
	int m_numOfResultPoints;//结果点的数量（循环后）
	std::map<int,Point> m_mapOfPoint;//点的map
	std::map<int,Triangle> m_mapOfTriangleOne;//原始三角形map
	std::map<int,Triangle> m_mapOfTriangleTwo;//结果三角形map
	std::map<int, Point>::iterator m_pointIterator;
	std::map<int,Triangle>::iterator m_triangleOneIterator;
	std::map<int,Triangle>::iterator m_triangleTwoIterator;
/*	Point m_originalVertexs[4];
	Triangle m_originalTriangles[4];
	numOfIteration=8;

	m_originalVertexs[0].vertexPosition[0]=-15*sqrt(3.0);
	m_originalVertexs[0].vertexPosition[1]=-15;
	m_originalVertexs[0].vertexPosition[2]=0;
	m_originalVertexs[1].vertexPosition[0]=15*sqrt(3.0);
	m_originalVertexs[1].vertexPosition[1]=-15;
	m_originalVertexs[1].vertexPosition[2]=0;
	m_originalVertexs[2].vertexPosition[0]=0;
	m_originalVertexs[2].vertexPosition[1]=30;
	m_originalVertexs[2].vertexPosition[2]=0;
	m_originalVertexs[3].vertexPosition[0]=0;
	m_originalVertexs[3].vertexPosition[1]=0;
	m_originalVertexs[3].vertexPosition[2]=45;


	m_originalTriangles[0].vertexIndex[0]=0;
	m_originalTriangles[0].vertexIndex[1]=1;
	m_originalTriangles[0].vertexIndex[2]=2;
	m_originalTriangles[1].vertexIndex[0]=1;
	m_originalTriangles[1].vertexIndex[1]=3;
	m_originalTriangles[1].vertexIndex[2]=2;
	m_originalTriangles[2].vertexIndex[0]=3;
	m_originalTriangles[2].vertexIndex[1]=1;
	m_originalTriangles[2].vertexIndex[2]=0;
	m_originalTriangles[3].vertexIndex[0]=2;
	m_originalTriangles[3].vertexIndex[1]=3;
	m_originalTriangles[3].vertexIndex[2]=0;
	*/

	numOfIteration=2;
	Point* m_originalVertexs;
	Triangle* m_originalTriangles;
	m_originalVertexs=new Point[numOfVertex];
	m_originalTriangles=new Triangle[numOfMeshes];

	for(int i=0;i<numOfVertex;i++){
		for(int j=0;j<3;j++)
		m_originalVertexs[i].vertexPosition[j]=vertex[i*3+j];
	}
	for(int i=0;i<numOfMeshes;i++){
		for(int j=0;j<3;j++){
			m_originalTriangles[i].vertexIndex[j]=triangleMesh[i*3+j];
		}
	}
	m_numOfOriginalPoints=numOfVertex;
	m_numOfOriginalMeshs=numOfMeshes;
	m_numOfResultPoints=numOfVertex;
	m_numOfResultMeshs=numOfMeshes;






	//////////////////////////////////////////////////////////////////////////
	// 初始化点map 和三角形map
	for(int i=0;i<m_numOfOriginalPoints;i++){
		m_mapOfPoint.insert(std::map<int, Point>::value_type(i,m_originalVertexs[i]));
	}
	for(int i=0;i<m_numOfOriginalMeshs;i++){
		m_mapOfTriangleOne.insert(std::map<int,Triangle>::value_type(i,m_originalTriangles[i]));
	}
	//m_numOfOriginalPoints=4;
	//m_numOfOriginalMeshs=4;
	//m_numOfResultPoints=4;
	//m_numOfResultMeshs=4;

	for(int i=0;i<numOfIteration;i++){
		m_numOfResultMeshs=0;
		m_numOfOriginalPoints=m_mapOfPoint.size();

		m_numOfOriginalMeshs=m_mapOfTriangleOne.size();
		for(int j=0;j<m_numOfOriginalMeshs;j++){
			Point tempPoint[3];
			Triangle resultTriangles[4];
			Point resultPoints[6];
			int resultPointsIndex[6];
			int tempNumOfPoints;
			m_triangleOneIterator=m_mapOfTriangleOne.find(j);


			for(int k=0;k<3;k++){
				m_pointIterator=m_mapOfPoint.find((*m_triangleOneIterator).second.vertexIndex[k]);
				resultPointsIndex[k]=(*m_triangleOneIterator).second.vertexIndex[k];
				for(int m=0;m<3;m++){
					tempPoint[k].vertexPosition[m]=(*m_pointIterator).second.vertexPosition[m];
				}
			}
			::subdivisionSingleTriangle(tempPoint,resultTriangles,resultPoints);
			//////////////////////////////////////////////////////////////////////////
			//插入点
			tempNumOfPoints=m_numOfResultPoints;
			for(int k=0;k<3;k++){
				float distance;
				bool flag=true;
				
				
				for(int m=m_numOfOriginalPoints;m<tempNumOfPoints;m++){
					m_pointIterator=m_mapOfPoint.find(m);
					distance=sqrt((resultPoints[k+3].vertexPosition[0]-(*m_pointIterator).second.vertexPosition[0])*(resultPoints[k+3].vertexPosition[0]-(*m_pointIterator).second.vertexPosition[0])
						+(resultPoints[k+3].vertexPosition[1]-(*m_pointIterator).second.vertexPosition[1])*(resultPoints[k+3].vertexPosition[1]-(*m_pointIterator).second.vertexPosition[1])
						+(resultPoints[k+3].vertexPosition[2]-(*m_pointIterator).second.vertexPosition[2])*(resultPoints[k+3].vertexPosition[2]-(*m_pointIterator).second.vertexPosition[2]));
					if(distance<0.1){
						resultPointsIndex[k+3]=m;
						flag=false;
						break;
					}
					
				}
				if(flag){
					m_mapOfPoint.insert(std::map<int,Point>::value_type(m_numOfResultPoints,resultPoints[k+3]));
					resultPointsIndex[k+3]=m_numOfResultPoints;
					m_numOfResultPoints++;
				}
			}

			//////////////////////////////////////////////////////////////////////////
			// 插入三角形
			for(int k=0;k<4;k++){
				Triangle tempTriangle;
				tempTriangle.vertexIndex[0]=resultPointsIndex[resultTriangles[k].vertexIndex[0]];
				tempTriangle.vertexIndex[1]=resultPointsIndex[resultTriangles[k].vertexIndex[1]];
				tempTriangle.vertexIndex[2]=resultPointsIndex[resultTriangles[k].vertexIndex[2]];
				m_mapOfTriangleTwo.insert(std::map<int,Triangle>::value_type(m_numOfResultMeshs,tempTriangle));
				m_numOfResultMeshs++;
			}
		}
		//////////////////////////////////////////////////////////////////////////
		// 重新给定原始网格，并删除结果网格
		std::map<int,Triangle>::iterator m_mapOfTriangleTwoIterator=m_mapOfTriangleTwo.begin();
		m_mapOfTriangleOne.swap(m_mapOfTriangleTwo);	
		m_mapOfTriangleTwo.clear();
	}

	//把网格数据写入一个文件
	CString filename_pw = "D:\\qin\\mesh data\\fandiskMesh.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		m_triangleOneIterator=m_mapOfTriangleOne.begin();
		for(;m_triangleOneIterator!=m_mapOfTriangleOne.end();m_triangleOneIterator++){
			int w[3];
			for(int i=0;i<3;i++){
				w[i]=(*m_triangleOneIterator).second.vertexIndex[i];
			}
			fprintf(fpout, "%d %d %d\n", w[0], w[1], w[2]);
		}
	}
	fclose(fpout);

	filename_pw = "D:\\qin\\mesh data\\Fandiskpoints.txt";
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		m_pointIterator=m_mapOfPoint.begin();
		for(;m_pointIterator!=m_mapOfPoint.end();m_pointIterator++){
			float w[3];
			for(int i=0;i<3;i++){
				w[i]=(*m_pointIterator).second.vertexPosition[i];
			}
			fprintf(fpout, "%f %f %f\n", w[0], w[1], w[2]);
		}
	}
	fclose(fpout);
}


double ComputeAreaOfTriangles(double position[3][3])
{
	double area;
	double edgeVectorOne[3],edgeVectorTwo[3];
	for(int i=0;i<3;i++){
		edgeVectorOne[i]=position[1][i]-position[0][i];
		edgeVectorTwo[i]=position[2][i]-position[0][i];
	}
	double areaVector[3];
	Cross3(edgeVectorOne,edgeVectorTwo,areaVector);
	area=Vector3Vector(areaVector,areaVector);
	area=sqrt(area)/2;
	return area;
}

void CalculateCentroidOrdinates(double trianglePosition[3][3],double anyPoint[3],double centroidOrdinate[3])
{
	double edgeVectorOne[3],edgeVectorTwo[3],edgeVectorThree[3],edgeVectorFour[3],areaVectorOne[3],areaVectorTwo[3];
	double areaOfTriangle=ComputeAreaOfTriangles(trianglePosition);
	double flag;
	double area;//与任意点构成的三角形的面积
	double centroid[3];
	for(int i=0;i<3;i++){
		centroid[i]=(trianglePosition[0][i]+trianglePosition[1][i]+trianglePosition[2][i])/3;
	}
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			edgeVectorOne[j]=trianglePosition[(int)(i+1)%3][j]-centroid[j];
			edgeVectorTwo[j]=trianglePosition[(int)(i+2)%3][j]-centroid[j];
			edgeVectorThree[j]=trianglePosition[(int)(i+1)%3][j]-anyPoint[j];
			edgeVectorFour[j]=trianglePosition[(int)(i+2)%3][j]-anyPoint[j];
		}
		Cross3(edgeVectorOne,edgeVectorTwo,areaVectorOne);
		Cross3(edgeVectorThree,edgeVectorFour,areaVectorTwo);
		flag=Vector3Vector(areaVectorOne,areaVectorTwo);
		area=Vector3Vector(areaVectorTwo,areaVectorTwo);
		area=sqrt(area)/2;
		if(flag>0||flag==0){
			centroidOrdinate[i]=area/areaOfTriangle;
		}
		else
			centroidOrdinate[i]=-area/areaOfTriangle;
	}
}

//////////////////////////////////////////////////////////////////////////
/// 对每个顶点分配颜色的函数
// 所对应的颜色为 四个顶点的颜色
void CalculateColorOfVertexs(int numOfPoints,float* pointSets,float* normals)
{
    int m_numOfPoints;
	float* colors;
	int m_numOfTriangles;
	float* m_pointSets;
	int* m_triangles;
	float m_pointColors[4][3];
	float* m_normals;
	Triangle triangleIndex[4];
	float triangleNormals[4][3];
	m_numOfPoints=numOfPoints;
	colors=new float[m_numOfPoints*3];
	m_pointSets=pointSets;
	m_normals=normals;
    // 初始话triangleIndex;
	triangleIndex[0].vertexIndex[0]=0;
	triangleIndex[0].vertexIndex[1]=1;
	triangleIndex[0].vertexIndex[2]=2;
	triangleIndex[1].vertexIndex[0]=1;
	triangleIndex[1].vertexIndex[1]=3;
	triangleIndex[1].vertexIndex[2]=2;
	triangleIndex[2].vertexIndex[0]=3;
	triangleIndex[2].vertexIndex[1]=1;
	triangleIndex[2].vertexIndex[2]=0;
	triangleIndex[3].vertexIndex[0]=2;
	triangleIndex[3].vertexIndex[1]=3;
	triangleIndex[3].vertexIndex[2]=0;
	// 初始化 triangleNormals



	// Calculate normals.
	for (int i = 0; i < 4; i++) {
		float vec1[3], vec2[3], normal[3];
		vec1[0]=m_pointSets[triangleIndex[i].vertexIndex[1]*3]-m_pointSets[triangleIndex[i].vertexIndex[0]*3];
		vec1[1]=m_pointSets[triangleIndex[i].vertexIndex[1]*3+1]-m_pointSets[triangleIndex[i].vertexIndex[0]*3+1];
		vec1[2]=m_pointSets[triangleIndex[i].vertexIndex[1]*3+2]-m_pointSets[triangleIndex[i].vertexIndex[0]*3+2];
		vec2[0]=m_pointSets[triangleIndex[i].vertexIndex[2]*3]-m_pointSets[triangleIndex[i].vertexIndex[0]*3];
		vec2[1]=m_pointSets[triangleIndex[i].vertexIndex[2]*3+1]-m_pointSets[triangleIndex[i].vertexIndex[0]*3+1];
		vec2[2]=m_pointSets[triangleIndex[i].vertexIndex[2]*3+2]-m_pointSets[triangleIndex[i].vertexIndex[0]*3+2];

		normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
		normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
		normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
		triangleNormals[i][0]=normal[0];
		triangleNormals[i][1]=normal[1];
		triangleNormals[i][2]=normal[2];

	}

	// Normalize normals.
	for (int i = 0; i < 4; i++) {
		float length = sqrt(triangleNormals[i][0]*triangleNormals[i][0]+triangleNormals[i][1]*triangleNormals[i][1]+triangleNormals[i][2]*triangleNormals[i][2]);
		triangleNormals[i][0]/=length;
		triangleNormals[i][1]/=length;
		triangleNormals[i][2]/=length;
	}

	//////////////////////////////////////////////////////////////////////////
	// 初始化 四个顶点的颜色
	m_pointColors[0][0]=0.8;
	m_pointColors[0][1]=0;
	m_pointColors[0][2]=0;

	m_pointColors[1][0]=0;
	m_pointColors[1][1]=0.8;
	m_pointColors[1][2]=0;

	m_pointColors[2][0]=0;
	m_pointColors[2][1]=0;
	m_pointColors[2][2]=0.8;

	m_pointColors[2][0]=0.5;
	m_pointColors[2][1]=0.5;
	m_pointColors[2][2]=0.5;

	//初始化colors 开始四个点的颜色
	for(int i=0;i<4;i++){
		colors[i*3]=m_pointColors[i][0];
		colors[i*3+1]=m_pointColors[i][1];
		colors[i*3+2]=m_pointColors[i][2];
	}

	for(int i=4;i<m_numOfPoints;i++){
		//////////////////////////////////////////////////////////////////////////
		//首先计算与四个三角形法向匹配情况
		int indexTriangle;// 目前点所在三角形的索引
		float distance=0;
		double tempTrianglePosition[3][3];
		double currentPoint[3];
		double currentCentroidCoordinate[3];
		for(int j=0;j<4;j++){
			if(m_normals[i*3]*triangleNormals[j][0]+m_normals[i*3+1]*triangleNormals[j][1]+m_normals[i*3+2]*triangleNormals[j][2]>distance){
				distance=m_normals[i*3]*triangleNormals[j][0]+m_normals[i*3+1]*triangleNormals[j][1]+m_normals[i*3+2]*triangleNormals[j][2];
				indexTriangle=j;
			}
		}

		for(int j=0;j<3;j++){
			tempTrianglePosition[j][0]=m_pointSets[triangleIndex[indexTriangle].vertexIndex[j]*3];
			tempTrianglePosition[j][1]=m_pointSets[triangleIndex[indexTriangle].vertexIndex[j]*3+1];
			tempTrianglePosition[j][2]=m_pointSets[triangleIndex[indexTriangle].vertexIndex[j]*3+2];
		}

		for(int j=0;j<3;j++){
			currentPoint[j]=m_pointSets[i*3+j];
		}

		CalculateCentroidOrdinates(tempTrianglePosition,currentPoint,currentCentroidCoordinate);
		for(int j=0;j<3;j++){
			colors[i*3+j]=currentCentroidCoordinate[0]*m_pointColors[triangleIndex[indexTriangle].vertexIndex[0]][j]
					+currentCentroidCoordinate[1]*m_pointColors[triangleIndex[indexTriangle].vertexIndex[1]][j]
					+currentCentroidCoordinate[2]*m_pointColors[triangleIndex[indexTriangle].vertexIndex[2]][j];
					
		}

	}
	CString filename_pw = "D:\\qin\\mesh data\\colors.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		
		for(int i=0;i<m_numOfPoints;i++){
			float w[3];
			for(int j=0;j<3;j++){
				w[i]=colors[i*3+j];
			}
			fprintf(fpout, "%f %f %f\n", w[0], w[1], w[2]);
		}
		fclose(fpout);
	}
	

}


void CalculateBoundaryColorOfVertexs(int numOfPoints,float* pointSets,float* normals)
{
	int m_numOfPoints;
	float* colors;
	int m_numOfTriangles;
	float* m_pointSets;
	int* m_triangles;
	float m_pointColors[4][3][3];//四个三角形三个顶点三种颜色
	float* m_normals;
	Triangle triangleIndex[4];
	float triangleNormals[4][3];
	m_numOfPoints=numOfPoints;
	colors=new float[m_numOfPoints*3];
	m_pointSets=pointSets;
	m_normals=normals;
	// 初始话triangleIndex;
	triangleIndex[0].vertexIndex[0]=0;
	triangleIndex[0].vertexIndex[1]=1;
	triangleIndex[0].vertexIndex[2]=2;
	triangleIndex[1].vertexIndex[0]=1;
	triangleIndex[1].vertexIndex[1]=3;
	triangleIndex[1].vertexIndex[2]=2;
	triangleIndex[2].vertexIndex[0]=3;
	triangleIndex[2].vertexIndex[1]=1;
	triangleIndex[2].vertexIndex[2]=0;
	triangleIndex[3].vertexIndex[0]=2;
	triangleIndex[3].vertexIndex[1]=3;
	triangleIndex[3].vertexIndex[2]=0;
	// 初始化 triangleNormals



	// Calculate normals.
	for (int i = 0; i < 4; i++) {
		float vec1[3], vec2[3], normal[3];
		vec1[0]=m_pointSets[triangleIndex[i].vertexIndex[1]*3]-m_pointSets[triangleIndex[i].vertexIndex[0]*3];
		vec1[1]=m_pointSets[triangleIndex[i].vertexIndex[1]*3+1]-m_pointSets[triangleIndex[i].vertexIndex[0]*3+1];
		vec1[2]=m_pointSets[triangleIndex[i].vertexIndex[1]*3+2]-m_pointSets[triangleIndex[i].vertexIndex[0]*3+2];
		vec2[0]=m_pointSets[triangleIndex[i].vertexIndex[2]*3]-m_pointSets[triangleIndex[i].vertexIndex[0]*3];
		vec2[1]=m_pointSets[triangleIndex[i].vertexIndex[2]*3+1]-m_pointSets[triangleIndex[i].vertexIndex[0]*3+1];
		vec2[2]=m_pointSets[triangleIndex[i].vertexIndex[2]*3+2]-m_pointSets[triangleIndex[i].vertexIndex[0]*3+2];

		normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
		normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
		normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
		triangleNormals[i][0]=normal[0];
		triangleNormals[i][1]=normal[1];
		triangleNormals[i][2]=normal[2];

	}

	// Normalize normals.
	for (int i = 0; i < 4; i++) {
		float length = sqrt(triangleNormals[i][0]*triangleNormals[i][0]+triangleNormals[i][1]*triangleNormals[i][1]+triangleNormals[i][2]*triangleNormals[i][2]);
		triangleNormals[i][0]/=length;
		triangleNormals[i][1]/=length;
		triangleNormals[i][2]/=length;
	}

	//////////////////////////////////////////////////////////////////////////
	// 初始化 四个顶点的颜色
	m_pointColors[0][0][0]=0.8;
	m_pointColors[0][0][1]=0;
	m_pointColors[0][0][2]=0;

	m_pointColors[0][1][0]=0.3;
	m_pointColors[0][1][1]=0;
	m_pointColors[0][1][2]=0;

	m_pointColors[0][2][0]=0.5;
	m_pointColors[0][2][1]=0;
	m_pointColors[0][2][2]=0;

	m_pointColors[1][0][0]=0;
	m_pointColors[1][0][1]=0.8;
	m_pointColors[1][0][2]=0;

	m_pointColors[1][1][0]=0;
	m_pointColors[1][1][1]=0.5;
	m_pointColors[1][1][2]=0;

	m_pointColors[1][2][0]=0;
	m_pointColors[1][2][1]=0.2;
	m_pointColors[1][2][2]=0;

	m_pointColors[2][0][0]=0;
	m_pointColors[2][0][1]=0;
	m_pointColors[2][0][2]=0.8;

	m_pointColors[2][1][0]=0;
	m_pointColors[2][1][1]=0;
	m_pointColors[2][1][2]=0.5;

	m_pointColors[2][2][0]=0;
	m_pointColors[2][2][1]=0;
	m_pointColors[2][2][2]=0.3;

	m_pointColors[3][0][0]=0.8;
	m_pointColors[3][0][1]=0.8;
	m_pointColors[3][0][2]=0;

	m_pointColors[3][1][0]=0.5;
	m_pointColors[3][1][1]=0.5;
	m_pointColors[3][1][2]=0;

	m_pointColors[3][2][0]=0.3;
	m_pointColors[3][2][1]=0.3;
	m_pointColors[3][2][2]=0;


	//初始化colors 开始四个点的颜色
	for(int i=0;i<4;i++){
		colors[i*3]=0.5;
		colors[i*3+1]=0.5;
		colors[i*3+2]=0.5;
	}

	for(int i=4;i<m_numOfPoints;i++){
		//////////////////////////////////////////////////////////////////////////
		//首先计算与四个三角形法向匹配情况
		int indexTriangle;// 目前点所在三角形的索引
		float distance=0;
		double tempTrianglePosition[3][3];
		double currentPoint[3];
		double currentCentroidCoordinate[3];

		if(i==131074){
			int d;
			d=0;
		}
		for(int j=0;j<4;j++){
			if(m_normals[i*3]*triangleNormals[j][0]+m_normals[i*3+1]*triangleNormals[j][1]+m_normals[i*3+2]*triangleNormals[j][2]>distance){
				distance=m_normals[i*3]*triangleNormals[j][0]+m_normals[i*3+1]*triangleNormals[j][1]+m_normals[i*3+2]*triangleNormals[j][2];
				indexTriangle=j;
			}
		}

		for(int j=0;j<3;j++){
			tempTrianglePosition[j][0]=m_pointSets[triangleIndex[indexTriangle].vertexIndex[j]*3];
			tempTrianglePosition[j][1]=m_pointSets[triangleIndex[indexTriangle].vertexIndex[j]*3+1];
			tempTrianglePosition[j][2]=m_pointSets[triangleIndex[indexTriangle].vertexIndex[j]*3+2];
		}

		for(int j=0;j<3;j++){
			currentPoint[j]=m_pointSets[i*3+j];
		}

		CalculateCentroidOrdinates(tempTrianglePosition,currentPoint,currentCentroidCoordinate);
		for(int j=0;j<3;j++){
			colors[i*3+j]=currentCentroidCoordinate[0]*m_pointColors[indexTriangle][0][j]
					+currentCentroidCoordinate[1]*m_pointColors[indexTriangle][1][j]
					+currentCentroidCoordinate[2]*m_pointColors[indexTriangle][2][j];

		}

	}
	CString filename_pw = "D:\\qin\\mesh data\\colors.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{

		for(int i=0;i<m_numOfPoints;i++){
			float w[3];
			for(int j=0;j<3;j++){
				w[j]=colors[i*3+j];
			}
			fprintf(fpout, "%f %f %f\n", w[0], w[1], w[2]);
		}
		fclose(fpout);
	}
	return;
}

