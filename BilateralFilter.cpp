
#include "stdafx.h"
#include <stdio.h>
#include ".\bilateralfilter.h"
#include "Delaunay.h"
#include "Matrix.h"
#include <math.h>
#include <cmath>
#include "GenerateNoise.h"
#include ".\mathlib\mathlib.h"
using namespace MATHLIB;

BilateralFilter::BilateralFilter(void)
{
	//输入数据
	m_originalPoints=NULL;//原始点
	m_originalNormals=NULL;// 原始点的法向
	//中间数据
	m_knearest=0;// k邻域的大小
	m_radius=0;// k 邻域的半径大小
	//m_mapKnearest=NULL;	
	//滤波时候用到的数据
	

	//m_LocalTriangles=NULL;//局部三角化


	//输出数据
    m_numOfPoints=0;//点的数量
	m_numOfTriangles=0;//三角片的数量
	m_triangles=NULL;// 三角片集合
	m_pointSets=NULL;//滤波后的点
	m_normals=NULL;//滤波后点的法向
	 m_meanShiftHposition=2;
	m_meanShiftHnormal=0.3;
	for(int i=0;i<3;i++)
		m_meanShiftHcolorGradient[i]=1;
	m_meanShiftStopNormal=0.01;
	m_meanShiftHcolor[0]=0.1;
	m_meanShiftHcolor[1]=0.1;
	m_meanShiftHcolor[2]=0.1;
	m_meanShiftStopColorGradient[0]=0.05;
	m_meanShiftStopColorGradient[1]=0.05;
	m_meanShiftStopColorGradient[2]=0.05;
}

BilateralFilter::~BilateralFilter(void)
{
	DeleteBilateralFilter();
}

void BilateralFilter::GetBilateralFilter(int numOfPointSets,float* pointSets,int numOfTriangle,int* triangles,float* vertexNormals,float* vertexColors)
{
	//////////////////////////////////////////////////////////////////////////
	//测试代码
	//long  idum;
	//idum=1;
	//for(int i=0;i<9;i++){
	//	float k,g;
	//	k=gasdev( &idum);
	//	g=ran1( &idum);
	//	g++;
	//}
	CString filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\localpoint.txt";
	FILE *fpout;

	//////////////////////////////////////////////////////////////////////////
	// 初始化原始数据

	m_numOfPoints=numOfPointSets;
	m_numOfTriangles=numOfTriangle;
	m_pointSets=new float[m_numOfPoints*3];
	m_originalPoints=new float[m_numOfPoints*3];
	m_originalNormals=new float[m_numOfPoints*3];
	m_originalFilterNormals=new float[m_numOfPoints*3];
	m_triangles=new int[m_numOfTriangles*3];
	m_originalColor=new float[m_numOfPoints*3];
	m_colors=new float[m_numOfPoints*3];
	//m_filterColor=new float[m_numOfPoints*3];
	m_colorGradient[0]=new float[m_numOfPoints*4];
	m_colorGradient[1]=new float[m_numOfPoints*4];
	m_colorGradient[2]=new float[m_numOfPoints*4];
	

	
	for(int i=0;i<m_numOfPoints;i++){
		for(int j=0;j<3;j++){
			m_originalPoints[i*3+j]=pointSets[i*3+j];
			m_originalNormals[i*3+j]=vertexNormals[i*3+j];
			m_originalColor[i*3+j]=vertexColors[i*3+j];
			for(int k=0;k<4;k++){
				m_colorGradient[j][i*4+k]=0;
			}
		}
		
	}
	for(int i=0;i<m_numOfTriangles*3;i++){
	       m_triangles[i]=triangles[i];
	}
	
  //  m_originalFilterNormals=m_originalNormals;


	//////////////////////////////////////////////////////////////////////////
	// 测试代码

	//////////////////////////////////////////////////////////////////////////
	

	//////////////////////////////////////////////////////////////////////////
	// 求邻域 （k＝20）

	m_knearest=10;
	m_radius=2;

	//////////////////////////////////////////////////////////////////////////
	// 加噪声
	AddPositionGaussianNoise();  
	/*for(int i=0;i<m_numOfPoints*3;i++){
		m_pointSets[i]=m_originalPoints[i];
	}*/
	AddColorGaussianNorse();
	//////////////////////////////////////////////////////////////////////////
	// 计算近邻点
	ComputeMapKnearest();
  
 
	//////////////////////////////////////////////////////////////////////////
	// 法向滤波
	for(int i=0;i<m_numOfPoints;i++){
		ComputeNormalFilter(i);
	}
	for(int i=0;i<m_numOfPoints;i++){
		for(int j=0;j<3;j++){
			m_originalNormals[i*3+j]=m_originalFilterNormals[i*3+j];
		}
	}
	delete[] m_originalFilterNormals;
	MeanShiftFilterNormal();

	//delete[] m_originalColor;
	
	CalculateFilterColorGradient();
	//////////////////////////////////////////////////////////////////////////
	// 测试代码
	filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\colorGradientOne.txt";

	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			float w[4];
			for(int k=0;k<4;k++){
				w[k]=m_filterColorGradient[0][i*4+k];
			}
			fprintf(fpout, "%f %f %f %f\n", w[0], w[1], w[2],w[3]);
		}
		fclose(fpout);
	}
	delete[] m_filterColorGradient[0];
	delete[] m_filterColorGradient[1];
	delete[] m_filterColorGradient[2];
	//////////////////////////////////////////////////////////////////////////
	
	MeanShiftFilterColorGradient();
	//////////////////////////////////////////////////////////////////////////
	// 测试代码
	filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\colorGradientTwo.txt";

	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			float w[4];
			for(int k=0;k<4;k++){
				w[k]=m_filterColorGradient[0][i*4+k];
			}
			fprintf(fpout, "%f %f %f %f\n", w[0], w[1], w[2],w[3]);
		}
		fclose(fpout);
	}

	//////////////////////////////////////////////////////////////////////////

//	ComputeFilterColor();
	
	//////////////////////////////////////////////////////////////////////////
	// 双边滤波
	for(int i=0;i<m_numOfPoints;i++){
		//求局部点
		if(!m_LocalPoints.empty())
			m_LocalPoints.clear();
		std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(i);		
		std::vector<int>::iterator mapNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
		m_LocalPoints.push_back(i);//把当前点作为第一个点放入
		for(;mapNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();mapNearestIterator++){
			m_LocalPoints.push_back((*mapNearestIterator));
		}
		//////////////////////////////////////////////////////////////////////////
		// 测试代码
		//std::vector<int>::iterator tempLocalPointIterator=m_LocalPoints.begin();
    	//int w;
		//for(;tempLocalPointIterator!=m_LocalPoints.end();tempLocalPointIterator++){
		//	w=(*tempLocalPointIterator);
		//	w++;
		//}
		//if(i==70){
		//	int d=0;
		//	d++;
		//}
		//
		//////////////////////////////////////////////////////////////////////////
		
		if(i==11){
			int dkj;
			dkj=0;
		}


		//对每个点求其局部三角片
		ComputeLocalTriangles();

		//求位置和颜色的一次预测
		ComputeOneOrderEstimation();

		
		//求滤波时候的系数
		ComputeVariation(i);
		//求双边滤波结果；
		ComputeBilateralFilter(i);
        //////////////////////////////////////////////////////////////////////////
        // 测试代码

		std::vector<Color>::iterator tempColorEstimation=m_LocalColorEstimation.begin();
		filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\localColor.txt";

		if((fpout = fopen(filename_pw, "w")) == NULL)
		{
			int dkjkd;
			//MessageBox("can't open the file!");
		}
		else
		{
			for(;tempColorEstimation!=m_LocalColorEstimation.end();tempColorEstimation++){
				float w[3];
				for(int i=0;i<3;i++){
					w[i]=(*tempColorEstimation).pointColor[i];
				}
				fprintf(fpout, "%f %f %f\n", w[0], w[1], w[2]);
			}
		}
		fclose(fpout);


	
				
		std::vector<int>::iterator  temptempLocalpoint=m_LocalPoints.begin();
		
		std::vector<POINTVECTOR3D>::iterator temptempLocalTriangleCentroid=m_LocalTriangleCentroid.begin();
		std::vector<POINTVECTOR3D>::iterator temptempLocalPositionEstimator=m_LocalPositionEstimation.begin();
		 filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\localpoint.txt";
		
		if((fpout = fopen(filename_pw, "w")) == NULL)
		{
			int dkjkd;
			//MessageBox("can't open the file!");
		}
		else
		{
			for(;temptempLocalpoint!=m_LocalPoints.end();temptempLocalpoint++){
				float w[3];
				for(int i=0;i<3;i++){
					w[i]=m_originalPoints[(*temptempLocalpoint)*3+i];
				}
				fprintf(fpout, "%f %f %f\n", w[0], w[1], w[2]);
			}
		}
		fclose(fpout);

		filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\trianglecentroid.txt";
		if((fpout = fopen(filename_pw, "w")) == NULL)
		{
			//MessageBox("can't open the file!");
			int lkd;
		}
		else
		{
			for(;temptempLocalTriangleCentroid!=m_LocalTriangleCentroid.end();temptempLocalTriangleCentroid++){
				float h[3];
				for(int i=0;i<3;i++){
					h[i]=(*temptempLocalTriangleCentroid).pointVector[i];
				}
				fprintf(fpout, "%f %f %f\n", h[0], h[1], h[2]);
			}
		}
		fclose(fpout);

		filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\pointprojectiontrianglecentroid.txt";
		if((fpout = fopen(filename_pw, "w")) == NULL)
		{
			//MessageBox("can't open the file!");
			int kdjk;
		}
		else
		{
			for(;temptempLocalPositionEstimator!=m_LocalPositionEstimation.end();temptempLocalPositionEstimator++){
				float wh[3];
				for(int i=0;i<3;i++){
					wh[i]=(*temptempLocalPositionEstimator).pointVector[i];
				}
				fprintf(fpout, "%f %f %f\n", wh[0], wh[1], wh[2]);
			}
		}
		fclose(fpout);
		std::vector<float>::iterator temptempLocalPointDistance=m_localPointDistance.begin();
		std::vector<float>::iterator temptempLocalProjectionDistance=m_localProjectionDistance.begin();
		filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\pointdistance.txt";
		if((fpout = fopen(filename_pw, "w")) == NULL)
		{
			//MessageBox("can't open the file!");
			int kdjk;
		}
		else
		{
			for(;temptempLocalPointDistance!=m_localPointDistance.end();temptempLocalPointDistance++){
				float whw;
				
					whw=(*temptempLocalPointDistance);
				
				fprintf(fpout, "%f \n", whw);
			}
		}
		fclose(fpout);

		filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\projectiondistance.txt";
		if((fpout = fopen(filename_pw, "w")) == NULL)
		{
			//MessageBox("can't open the file!");
			int kdjk;
		}
		else
		{
			for(;temptempLocalProjectionDistance!=m_localProjectionDistance.end();temptempLocalProjectionDistance++){
				float whwh;

				whwh=(*temptempLocalProjectionDistance);

				fprintf(fpout, "%f \n", whwh);
			}
		}
		fclose(fpout);

		std::vector<Triangle>::iterator temptempLocalTriangle=m_LocalTriangles.begin();
		filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\localtriangle.txt";
		if((fpout = fopen(filename_pw, "w")) == NULL)
		{
			//MessageBox("can't open the file!");
			int kdjk;
		}
		else
		{
			for(;temptempLocalTriangle!=m_LocalTriangles.end();temptempLocalTriangle++){
				int b[3];
				for(int i=0;i<3;i++){
					b[i]=(*temptempLocalTriangle).vertexIndex[i];
				}

				fprintf(fpout, "%d %d %d \n", b[0],b[1],b[2]);
			}
		}
		fclose(fpout);
		if(i==10){
			int d=0;
			d++;
		}

		

		if(i==70){
			int d=0;
			d++;
		}

		filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\colorestimate.txt";
		if((fpout = fopen(filename_pw, "w")) == NULL)
		{
			//MessageBox("can't open the file!");
			int kdjk;
		}
		else
		{
			std::vector<Color>::iterator tempLocalColorEstimateIterator=m_LocalColorEstimation.begin();
			for(;tempLocalColorEstimateIterator!=m_LocalColorEstimation.end();tempLocalColorEstimateIterator++){
				float whwh[3];

				for(int kd=0;kd<3;kd++){
					whwh[kd]=(*tempLocalColorEstimateIterator).pointColor[kd];
				}
				

				fprintf(fpout, "%f %f %f\n", whwh[0],whwh[1],whwh[2]);
			}
		}
		fclose(fpout);

		//////////////////////////////////////////////////////////////////////////
	


		m_LocalPoints.clear();//局部点集合，其中第一个点为要滤波的点
		m_LocalTriangleCentroid.clear();//局部三角形的重心
		m_originalCentroidNormals.clear();////局部三角形的重心法向，由三个顶点法向加权	
		m_LocalPositionEstimation.clear();//局部位置的一阶估计
		m_localPointDistance.clear();//点与点之间的距离，就是双边滤波中的第一式子
		m_localProjectionDistance.clear();//滤波点到投影点之间的距离
		//m_localColorDistance.clear();
		m_localColorDistance[0].clear();//颜色距离
		m_localColorDistance[1].clear();
		m_localColorDistance[2].clear();		
		m_LocalColorEstimation.clear();//局部颜色的一次估计
		m_LocalTriangles.clear();//局部三角化 三个点三个点存放
	}


	filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\originalpoint.txt";
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		//MessageBox("can't open the file!");
		int kdjk;
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			float b[3];
			for(int j=0;j<3;j++){
				b[j]=m_originalPoints[i*3+j];
			}

			fprintf(fpout, "%f %f %f \n", b[0],b[1],b[2]);
		}
	}
	fclose(fpout);

	filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\resultpoint.txt";
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		//MessageBox("can't open the file!");
		int kdjk;
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			float b[3];
			for(int j=0;j<3;j++){
				b[j]=m_pointSets[i*3+j];
			}

			fprintf(fpout, "%f %f %f \n", b[0],b[1],b[2]);
		}
	}
	fclose(fpout);

	filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\resultColor.txt";
	float sumOfColor[3];
	for(int i=0;i<3;i++){
		sumOfColor[i]=0;
	}
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		//MessageBox("can't open the file!");
		int kdjk;
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			float b[3];
			for(int j=0;j<3;j++){
				b[j]=m_colors[i*3+j];
				sumOfColor[j]+=b[j];
			}

			fprintf(fpout, "%f %f %f \n", b[0],b[1],b[2]);
		}
		float averageColor[3];
		for(int i=0;i<3;i++){
            averageColor[i]=sumOfColor[i]/m_numOfPoints;
		}
		fprintf(fpout, "%f %f %f \n", averageColor[0],averageColor[1],averageColor[2]);
		float variationColor[3];
		for(int i=0;i<3;i++){
			variationColor[i]=0;
		}
		for(int i=0;i<m_numOfPoints;i++){
			for(int j=0;j<3;j++)
			{
                variationColor[j]+=(averageColor[j]-m_colors[i*3+j])*(averageColor[j]-m_colors[i*3+j]);
			}
		}
		for(int i=0;i<3;i++){
			variationColor[i]/=m_numOfPoints;
		}
		fprintf(fpout, "%f %f %f \n", variationColor[0],variationColor[1],variationColor[2]);

	}
	
	fclose(fpout);

	delete[] m_originalPoints;//原始点
	delete[] m_originalColor;
	//delete[] m_originalColor;
	//delete[] m_originalNormals;// 原始点的法向*/
	//delete[] m_originalFilterNormals;//把原始点的法向滤波后的法向
	delete[] m_originalFilterNormals;
	ComputeNormals();
	return; 
}

void BilateralFilter::DeleteBilateralFilter()
{
	if(m_pointSets!=NULL)
		delete m_pointSets;
	m_pointSets=NULL;
	if(m_triangles!=NULL)
		delete m_triangles;
	m_triangles=NULL;
	if(m_normals!=NULL)
		delete m_normals;
	m_normals=NULL;
}
//求点的局部邻域，对三角形来说，就是求它的2－ring邻域
void BilateralFilter::ComputeMapKnearest()
{
	float radius=m_radius*m_radius;//存放近邻点到所求点的距离
	std::multimap<float ,int> mapKnearest;//对每一个点建立邻域vector时建立临时map
    int numOfKnearest;//统计一个点的临时邻域大小

	//////////////////////////////////////////////////////////////////////////
	//遍历所有点，并且按照距离由小到大存放存放到map中
	for(int i=0;i<m_numOfPoints;i++){
		numOfKnearest=0;
		for(int j=0;j<m_numOfPoints;j++){
			if(j==i)
				continue;
			float distancePointToPoint;
			distancePointToPoint=(m_originalPoints[i*3]-m_originalPoints[j*3])*(m_originalPoints[i*3]-m_originalPoints[j*3])
				+(m_originalPoints[i*3+1]-m_originalPoints[j*3+1])*(m_originalPoints[i*3+1]-m_originalPoints[j*3+1])
				+(m_originalPoints[i*3+2]-m_originalPoints[j*3+2])*(m_originalPoints[i*3+2]-m_originalPoints[j*3+2]);
			if(distancePointToPoint>radius)
				continue;
			mapKnearest.insert(std::multimap<float,int>::value_type(distancePointToPoint,j));
			numOfKnearest+=1;
			if(numOfKnearest>m_knearest){
				std::multimap<float,int>::iterator mapIterator = mapKnearest.end();
				mapIterator--;
				float w;
				w=(*mapIterator).first;
				mapKnearest.erase(mapIterator);
				numOfKnearest-=1;
			}
		}
		KnearestField fieldKnearest;//临时结构
		fieldKnearest.m_numOfNearest=numOfKnearest;
		fieldKnearest.m_IfBoundary=FALSE;
		std::multimap<float,int>::iterator mapIterator = mapKnearest.begin();
		for(;mapIterator!=mapKnearest.end();mapIterator++){
			fieldKnearest.m_nearest.push_back((*mapIterator).second);
		}
		m_mapKnearest.insert(std::map<int, KnearestField>::value_type(i,fieldKnearest));
		//////////////////////////////////////////////////////////////////////////
		// 测试代码
        mapIterator=mapKnearest.begin();
		for(;mapIterator!=mapKnearest.end();mapIterator++){
			float k;
			int g;
			k=(*mapIterator).first;
			g=(*mapIterator).second;
			k++;
			g++;
		}
		//////////////////////////////////////////////////////////////////////////
		mapKnearest.clear();
	}	
    //////////////////////////////////////////////////////////////////////////
    // 测试代码
	std::map<int,KnearestField>::iterator mapTempIterator=m_mapKnearest.begin();
	for(;mapTempIterator!=m_mapKnearest.end();mapTempIterator++){
		std::vector<int>::iterator tempTempIterator=(*mapTempIterator).second.m_nearest.begin();
        
		int k=(*mapTempIterator).first;
		int g=(*tempTempIterator);
		g++;
		k++;
	}	
	//////////////////////////////////////////////////////////////////////////
	
}

//求点的局部三角化
void BilateralFilter::ComputeLocalTriangles()
{
    //////////////////////////////////////////////////////////////////////////
    //初始化m_numOFLocalPoints
	std::vector<int>::iterator localPointIterator=m_LocalPoints.begin();
	//////////////////////////////////////////////////////////////////////////
	// 测试代码
	//CString filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\localpointIndex.txt";
	//FILE *fpout;
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(;localPointIterator!=m_LocalPoints.end();localPointIterator++){
	//		int a;
	//		a=(*localPointIterator);
	//		
	//		fprintf(fpout, "%d\n",a);
	//	}fclose(fpout);
	//}
	//
	//localPointIterator=m_LocalPoints.begin();
	//////////////////////////////////////////////////////////////////////////
	//测试代码
	
//	std::map<int,KnearestField>::iterator mapKnearestFieldIterator=m_mapKnearest.find((*localPointIterator));
//	m_numOFLocalPoints=(*mapKnearestFieldIterator).second.m_numOfNearest+1;
	m_numOFLocalPoints=m_LocalPoints.size();	
	//////////////////////////////////////////////////////////////////////////
	
	float localXaxis[3],localYaxis[3];//局部坐标系下的坐标轴
	float* dkjk= new float[50];
	float* temp2DPosition=new float[m_numOFLocalPoints*2];
	float localNormal[3];
	//Delaunay类中的一些变量
	vertexSet vSetFront;
	Delaunay delaunay;
	triangleSet tSetFront;
	tIterator tIt;

	for(int i=0;i<3;i++){
		//localNormal[i]=m_originalNormals[(*localPointIterator)*3+i];
		localNormal[i]=m_originalFilterNormals[(*localPointIterator)*3+i];
	}
	//清空局部三角形的vector
	if(!m_LocalTriangles.empty())
		m_LocalTriangles.clear();
	
	localPointIterator=m_LocalPoints.begin();
	//////////////////////////////////////////////////////////////////////////
	//测试代码
	if((*localPointIterator)==10){
		int dwh;
		dwh=0;
	}
	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	//首先把这些点投影到一个局部平面上
	//建立局部平面
	Perpendicular(localNormal,localXaxis);
	Cross3(localNormal,localXaxis,localYaxis);
	//投影并把每个点放到一个set中
	for(int i=0;i<m_numOFLocalPoints, localPointIterator!=m_LocalPoints.end();i++,localPointIterator++){
		float tempX,tempY;
		float tempPosition[3];
		int indexOfPoint;
		indexOfPoint=(*localPointIterator);
		
		for(int j=0;j<3;j++){
			tempPosition[j]=m_originalPoints[indexOfPoint*3+j];
		}
		tempX=Vector3Vector(tempPosition,localXaxis);
		tempY=Vector3Vector(tempPosition,localYaxis);
		temp2DPosition[i*2]=tempX;
		temp2DPosition[i*2+1]=tempY;
		vSetFront.insert(vertex(tempX,tempY));
	}
//	delaunay.Triangulate(vSetFront,tSetFront);//三角化
	m_numOfLocalTriangles=0;
    int i;
	
	for(i=0,tIt = tSetFront.begin();tIt!= tSetFront.end(); tIt++, i++){
		Triangle tempTriangle;

		for(int j=0;j<3;j++){
		    int k;
			for(k=0,localPointIterator=m_LocalPoints.begin();localPointIterator!=m_LocalPoints.end();localPointIterator++,k++){
				if(abs(tIt->m_Vertices[j]->GetX()-temp2DPosition[k*2])<0.0000001&&abs(tIt->m_Vertices[j]->GetY()-temp2DPosition[k*2+1])<0.0000001){
					tempTriangle.vertexIndex[j]=(*localPointIterator);					
					break;
				}
			}
		}
		
		m_LocalTriangles.push_back(tempTriangle);
		m_numOfLocalTriangles+=1;
	}

	vSetFront.clear();
	tSetFront.clear();
	delete[] temp2DPosition;
	return;
}

//求位置和颜色的一次预测
void BilateralFilter::ComputeOneOrderEstimation()
{
	std::vector<Triangle>::iterator triangleIterator=m_LocalTriangles.begin();
	//////////////////////////////////////////////////////////////////////////
	//测试代码
	CString filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\templocaltriangle.txt";
	FILE *fpout;

	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		//MessageBox("can't open the file!");
		int kdjk;
	}
	else
	{
		for(;triangleIterator!=m_LocalTriangles.end();triangleIterator++){
			int b[3];
			for(int i=0;i<3;i++){
				b[i]=(*triangleIterator).vertexIndex[i];
			}

			fprintf(fpout, "%d %d %d \n", b[0],b[1],b[2]);
		}
	}
	fclose(fpout);
	triangleIterator=m_LocalTriangles.begin();
	//////////////////////////////////////////////////////////////////////////
	 
	std::vector<int>::iterator pointIterator=m_LocalPoints.begin();
	for(;triangleIterator!=m_LocalTriangles.end();triangleIterator++){
	//	if((*triangleIterator).vertexIndex[0]!=(*pointIterator)&&(*triangleIterator).vertexIndex[1]!=(*pointIterator)&&(*triangleIterator).vertexIndex[0]!=(*pointIterator)){//判断是否是包含被滤波的点的三角形
			//计算重心和重心法向
			POINTVECTOR3D tempCentroid, tempProjection;
			float tempCentroidNormal[3];			
			Color tempColor;
			
			//初始化
			for(int i=0;i<3;i++){
				tempCentroid.pointVector[i]=0;tempProjection.pointVector[i]=0;tempCentroidNormal[i]=0;tempColor.pointColor[i]=0;
			}
			//计算重心重心法向
			float tempNuniformNomal[3];
			for(int i=0;i<3;i++){
				tempCentroid.pointVector[i]=(m_originalPoints[(*triangleIterator).vertexIndex[0]*3+i]
				              +m_originalPoints[(*triangleIterator).vertexIndex[1]*3+i]
							  +m_originalPoints[(*triangleIterator).vertexIndex[2]*3+i])/3;
				tempNuniformNomal[i]=(m_originalFilterNormals[(*triangleIterator).vertexIndex[0]*3+i]
								+m_originalFilterNormals[(*triangleIterator).vertexIndex[1]*3+i]
								+m_originalFilterNormals[(*triangleIterator).vertexIndex[2]*3+i])/3;

			}
			UnitVector(tempNuniformNomal,tempCentroidNormal);
			//求滤波点在当前三角片上的投影
			float distancePointToPlane;//滤波点到以过重心且垂直于重心法向的平面的距离
			float vectorPointToPoint[3];//滤波点到重心点的向量
			for(int i=0;i<3;i++){
				vectorPointToPoint[i]=m_originalPoints[(*pointIterator)*3+i]-tempCentroid.pointVector[i];
			}
            distancePointToPlane=(vectorPointToPoint[0]*tempCentroidNormal[0]
								+vectorPointToPoint[1]*tempCentroidNormal[1]
								+vectorPointToPoint[2]*tempCentroidNormal[2]);
			
			for(int i=0;i<3;i++){
				tempProjection.pointVector[i]=m_originalPoints[(*pointIterator)*3+i]-distancePointToPlane*tempCentroidNormal[i];
			}

			for(int i=0;i<3;i++){
				tempColor.pointColor[i]=(m_filterColorGradient[i][(*triangleIterator).vertexIndex[0]*4]+m_filterColorGradient[i][(*triangleIterator).vertexIndex[1]*4]+m_filterColorGradient[i][(*triangleIterator).vertexIndex[2]*4])/3*tempProjection.pointVector[0]
				+(m_filterColorGradient[i][(*triangleIterator).vertexIndex[0]*4+1]+m_filterColorGradient[i][(*triangleIterator).vertexIndex[1]*4+1]+m_filterColorGradient[i][(*triangleIterator).vertexIndex[2]*4+1])/3*tempProjection.pointVector[1]
				+(m_filterColorGradient[i][(*triangleIterator).vertexIndex[0]*4+2]+m_filterColorGradient[i][(*triangleIterator).vertexIndex[1]*4+2]+m_filterColorGradient[i][(*triangleIterator).vertexIndex[2]*4+2])/3*tempProjection.pointVector[2]
				+(m_filterColorGradient[i][(*triangleIterator).vertexIndex[0]*4+3]+m_filterColorGradient[i][(*triangleIterator).vertexIndex[1]*4+3]+m_filterColorGradient[i][(*triangleIterator).vertexIndex[2]*4+3])/3;
			}
			m_LocalPositionEstimation.push_back(tempProjection);
			m_LocalTriangleCentroid.push_back(tempCentroid);
			//m_LocalColorEstimation.push_back(tempColor);
	//	}
			//////////////////////////////////////////////////////////////////////////
			//估计颜色
			//常数估计
			/*for(int k=0;k<3;k++){
				tempColor.pointColor[k]=(m_originalColor[(*triangleIterator).vertexIndex[0]*3+k]
									+m_originalColor[(*triangleIterator).vertexIndex[1]*3+k]
									+m_originalColor[(*triangleIterator).vertexIndex[2]*3+k])/3;
			}*/
							
	//	   ComputeProjectionPointColor(*triangleIterator,tempCentroidNormal,tempCentroid,tempProjection,tempColor);

			
		
			m_LocalColorEstimation.push_back(tempColor);
	}
    return;
}

//求位置和颜色的1.5次预测
void BilateralFilter::ComputeOneHalfEstimation(int indexPoint)
{
}

//求双边滤波结果；
void BilateralFilter::ComputeBilateralFilter( int indexPoint)
{
	
	float sumOfWeightPosition, sumOfWeightColor[3];//加权平均值
	sumOfWeightPosition=0;
	//sumOfWeightColor=0;
	for(int i=0;i<3;i++)
        sumOfWeightColor[i]=0;
    float weightSpatial,weightInfluenceFunction,weightColor[3];//双边滤波里面的
	for(int i=0;i<3;i++){
		m_pointSets[indexPoint*3+i]=0;
		m_colors[indexPoint*3+i]=0;
	}
	
	std::vector<POINTVECTOR3D>::iterator localPositionEstimationIterator=m_LocalPositionEstimation.begin();
	std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();
	std::vector<float>::iterator localProjectionDistanceIterator=m_localProjectionDistance.begin();
	std::vector<float>::iterator localColorDistanceIterator[3];
	//localColorDistanceIterator=m_localColorDistance.begin();
	for(int k=0;k<3;k++){
		localColorDistanceIterator[k]=m_localColorDistance[k].begin();
	}
	std::vector<Color>::iterator localColorEstimation=m_LocalColorEstimation.begin();


	for(;localPositionEstimationIterator!=m_LocalPositionEstimation.end();localPositionEstimationIterator++,localPointDistanceIterator++,localProjectionDistanceIterator++,localColorDistanceIterator[0]++,localColorDistanceIterator[1]++,localColorDistanceIterator[2]++,localColorEstimation++){
		weightSpatial=exp(-0.5*(*localPointDistanceIterator)/m_pointVariation);
		weightInfluenceFunction=exp(-0.5*(*localProjectionDistanceIterator)/m_estimationVariation);
		//////////////////////////////////////////////////////////////////////////
		//假定颜色为常数时候的滤波


		//////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////
		//假定颜色为线性函数
		//weightColor=exp(-0.5*(*localColorDistanceIterator))/m_clolrVariation;
		
		weightColor[0]=exp(-0.5*(*localColorDistanceIterator[0])/ m_clolrVariation[0]);
		weightColor[1]=exp(-0.5*(*localColorDistanceIterator[1])/ m_clolrVariation[1]);
		weightColor[2]=exp(-0.5*(*localColorDistanceIterator[2])/ m_clolrVariation[2]);
		//////////////////////////////////////////////////////////////////////////
		

		float tempWeightPosition=weightSpatial*weightInfluenceFunction;
		float tempWeightColor[3];
		//tempWeightColor=tempWeightPosition*weightColor;
		for(int k=0;k<3;k++)
			tempWeightColor[k]=tempWeightPosition*weightColor[k];
		for(int k=0;k<3;k++){
			m_pointSets[indexPoint*3+k]+=(*localPositionEstimationIterator).pointVector[k]*tempWeightPosition;
			//////////////////////////////////////////////////////////////////////////
			//颜色为常数
			
			//m_colors[indexPoint*3+k]+=m_filterColor[indexPoint*3+k]*tempWeightPosition;
			//////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////
			//颜色为线性假定时
			m_colors[indexPoint*3+k]+=(*localColorEstimation).pointColor[k]*tempWeightPosition;
			//m_colors[indexPoint*3+k]+=(*localColorEstimation).pointColor[k]*tempWeightColor[k];
			//////////////////////////////////////////////////////////////////////////
			
		}
		sumOfWeightPosition=sumOfWeightPosition+tempWeightPosition;
		//sumOfWeightColor=sumOfWeightColor+tempWeightColor;
		for(int k=0;k<3;k++)
			sumOfWeightColor[k]=sumOfWeightColor[k]+tempWeightPosition;
            //sumOfWeightColor[k]=sumOfWeightColor[k]+tempWeightColor[k];
	}

	for(int i=0;i<3;i++){
		m_pointSets[indexPoint*3+i]+=m_originalPoints[indexPoint*3+i]*exp(0.0)*exp(0.0);
		//////////////////////////////////////////////////////////////////////////
		// 颜色为常数
	//	m_colors[indexPoint*3+i]+=m_filterColor[indexPoint*3+i]*exp(0.0)*exp(0.0);
		//////////////////////////////////////////////////////////////////////////
		//颜色为线性函数
		m_colors[indexPoint*3+i]+=m_originalColor[indexPoint*3+i]*exp(0.0)*exp(0.0);

		//m_colors[indexPoint*3+i]+=m_originalColor[indexPoint*3+i]*exp(0.0)*exp(0.0)*exp(0.0);
	}
	sumOfWeightPosition=sumOfWeightPosition+exp(0.0)*exp(0.0);
	//sumOfWeightColor+=exp(0.0)*exp(0.0)*exp(0.0);
	for(int k=0;k<3;k++)
		sumOfWeightColor[k]+=exp(0.0)*exp(0.0);
        //sumOfWeightColor[k]+=exp(0.0)*exp(0.0)*exp(0.0);
	for(int i=0;i<3;i++){
		m_pointSets[indexPoint*3+i]/=sumOfWeightPosition;
		// 颜色为常数
		//m_colors[indexPoint*3+i]/=sumOfWeightPosition;
		//颜色为线性函数
		m_colors[indexPoint*3+i]/=sumOfWeightColor[i];
	}
   return;
}

//求法向的双边滤波（molification)
void BilateralFilter::ComputeNormalFilter( int indexPoint)
{
	float sumweight,sumdistance,distanceVariance;
	sumweight=0;
	sumdistance=0;
	float* distance;
	
	std::map<int,KnearestField>::iterator mapIterator=m_mapKnearest.find(indexPoint);
	std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	distance=new float[(*mapIterator).second.m_nearest.size()];
	for(int i=0;i<3;i++){
		m_originalFilterNormals[indexPoint*3+i]=0;
	//	m_filterColor[indexPoint*3+i]=0;
	}
    
	for(int i=0;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,i++){
		
		distance[i]=(m_originalPoints[indexPoint*3]-m_originalPoints[(*vectorIterator)*3])*(m_originalPoints[indexPoint*3]-m_originalPoints[(*vectorIterator)*3])
			+(m_originalPoints[indexPoint*3+1]-m_originalPoints[(*vectorIterator)*3+1])*(m_originalPoints[indexPoint*3+1]-m_originalPoints[(*vectorIterator)*3+1])
			+(m_originalPoints[indexPoint*3+2]-m_originalPoints[(*vectorIterator)*3+2])*(m_originalPoints[indexPoint*3+2]-m_originalPoints[(*vectorIterator)*3+2]);
		sumdistance+=distance[i];
	}
	distanceVariance=sumdistance/((*mapIterator).second.m_nearest.size());
	vectorIterator=(*mapIterator).second.m_nearest.begin();
	for(int i=0;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,i++){
		float weight;
		weight=exp(-0.5*distance[i]/distanceVariance);
		for(int j=0;j<3;j++){
			m_originalFilterNormals[indexPoint*3+j]+=weight*m_originalNormals[(*vectorIterator)*3+j];
		//	m_filterColor[indexPoint*3+j]+=weight*m_originalColor[(*vectorIterator)*3+j];
		}
		sumweight+=weight;

	}
	for(int i=0;i<3;i++){
		m_originalFilterNormals[indexPoint*3+i]+=exp(0.0)*m_originalNormals[indexPoint*3+i];
	//	m_filterColor[indexPoint*3+i]+=exp(0.0)*m_originalColor[indexPoint*3+i];

	}
	sumweight+=exp(0.0);

	for(int i=0;i<3;i++){
		m_originalFilterNormals[indexPoint*3+i]/=sumweight;
	//	m_filterColor[indexPoint*3+i]/=sumweight;
	}
	float sum=0;

	for(int i=0;i<3;i++){
		sum+=m_originalFilterNormals[indexPoint*3+i]*m_originalFilterNormals[indexPoint*3+i];
	}
	sum=sqrt(sum);

	for(int i=0;i<3;i++){
		m_originalFilterNormals[indexPoint*3+i]/=sum;
	}

	delete[] distance;

	return;
}

//求滤波的函数的方差值 以及滤波中用到的距离的平方
void BilateralFilter::ComputeVariation( int indexPoint)
{
    //////////////////////////////////////////////////////////////////////////
    // 计算
	float localPointDistance,localProjectionDistance,localColorDistance[3];
	std::vector<POINTVECTOR3D>::iterator localTriangleCentroidIterator=m_LocalTriangleCentroid.begin();
	std::vector<POINTVECTOR3D>::iterator localPositionEstimationIterator=m_LocalPositionEstimation.begin();
	std::vector<Color>::iterator localColorEstimationIterator=m_LocalColorEstimation.begin();

	for(;localTriangleCentroidIterator!=m_LocalTriangleCentroid.end();localTriangleCentroidIterator++,localPositionEstimationIterator++,localColorEstimationIterator++){
		localPointDistance=((m_originalPoints[indexPoint*3]-(*localTriangleCentroidIterator).pointVector[0])*(m_originalPoints[indexPoint*3]-(*localTriangleCentroidIterator).pointVector[0])
			+(m_originalPoints[indexPoint*3+1]-(*localTriangleCentroidIterator).pointVector[1])*(m_originalPoints[indexPoint*3+1]-(*localTriangleCentroidIterator).pointVector[1])
			+(m_originalPoints[indexPoint*3+2]-(*localTriangleCentroidIterator).pointVector[2])*(m_originalPoints[indexPoint*3+2]-(*localTriangleCentroidIterator).pointVector[2]));
		localProjectionDistance=((m_originalPoints[indexPoint*3]-(*localPositionEstimationIterator).pointVector[0])*(m_originalPoints[indexPoint*3]-(*localPositionEstimationIterator).pointVector[0])
			+(m_originalPoints[indexPoint*3+1]-(*localPositionEstimationIterator).pointVector[1])*(m_originalPoints[indexPoint*3+1]-(*localPositionEstimationIterator).pointVector[1])
			+(m_originalPoints[indexPoint*3+2]-(*localPositionEstimationIterator).pointVector[2])*(m_originalPoints[indexPoint*3+2]-(*localPositionEstimationIterator).pointVector[2]));
	/*	localColorDistance=(m_originalColor[indexPoint*3]-(*localColorEstimationIterator).pointColor[0])*(m_originalColor[indexPoint*3]-(*localColorEstimationIterator).pointColor[0])
					+(m_originalColor[indexPoint*3+1]-(*localColorEstimationIterator).pointColor[1])*(m_originalColor[indexPoint*3+1]-(*localColorEstimationIterator).pointColor[1])
					+(m_originalColor[indexPoint*3+2]-(*localColorEstimationIterator).pointColor[2])*(m_originalColor[indexPoint*3+2]-(*localColorEstimationIterator).pointColor[2]);*/
		localColorDistance[0]=(m_originalColor[indexPoint*3]-(*localColorEstimationIterator).pointColor[0])*(m_originalColor[indexPoint*3]-(*localColorEstimationIterator).pointColor[0]);
		localColorDistance[1]=(m_originalColor[indexPoint*3+1]-(*localColorEstimationIterator).pointColor[1])*(m_originalColor[indexPoint*3+1]-(*localColorEstimationIterator).pointColor[1]);
		localColorDistance[2]=(m_originalColor[indexPoint*3+2]-(*localColorEstimationIterator).pointColor[2])*(m_originalColor[indexPoint*3+2]-(*localColorEstimationIterator).pointColor[2]);
       
		m_localPointDistance.push_back(localPointDistance);
		m_localProjectionDistance.push_back(localProjectionDistance);	
		//m_localColorDistance.push_back(localColorDistance);
		m_localColorDistance[0].push_back(localColorDistance[0]);
		m_localColorDistance[1].push_back(localColorDistance[1]);
		m_localColorDistance[2].push_back(localColorDistance[2]); 

	}

	std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();
	std::vector<float>::iterator localProjectionDistanceIterator=m_localProjectionDistance.begin();
	std::vector<float>::iterator localColorDistanceIterator[3];
	//localColorDistanceIterator=m_localColorDistance.begin();
	for(int i=0;i<3;i++){
		localColorDistanceIterator[i]=m_localColorDistance[i].begin();
	}
	m_pointVariation=0; 
	m_estimationVariation=0;
	//m_clolrVariation=0;
	for(int i=0;i<3;i++){
        m_clolrVariation[i]=0;
	}
	for(;localPointDistanceIterator!=m_localPointDistance.end();localPointDistanceIterator++,localProjectionDistanceIterator++,localColorDistanceIterator[0]++,localColorDistanceIterator[1]++,localColorDistanceIterator[2]++){
		m_pointVariation+=(*localPointDistanceIterator);
		m_estimationVariation+=(*localProjectionDistanceIterator);
		//m_clolrVariation+=(*localColorDistanceIterator);
		m_clolrVariation[0]+=(*localColorDistanceIterator[0]);
		m_clolrVariation[1]+=(*localColorDistanceIterator[1]);
		m_clolrVariation[2]+=(*localColorDistanceIterator[2]);
	}

	m_pointVariation=m_pointVariation/(int) m_localPointDistance.size();
	if(m_pointVariation==0)
		m_pointVariation=2;
	m_estimationVariation=m_estimationVariation/(int) m_localPointDistance.size();
	if(m_estimationVariation==0)
		m_estimationVariation=2;
	//m_clolrVariation=0.2;
	//m_clolrVariation=m_clolrVariation/(int)m_localPointDistance.size();
	//if(m_clolrVariation<0.000001)
	//	m_clolrVariation=0.05;
	m_clolrVariation[0]=m_clolrVariation[0]/(int)m_localPointDistance.size();
	m_clolrVariation[1]=m_clolrVariation[1]/(int)m_localPointDistance.size();
	m_clolrVariation[2]=m_clolrVariation[2]/(int)m_localPointDistance.size();
	for(int k=0;k<3;k++){
		if(m_clolrVariation[k]<0.000001)
			m_clolrVariation[k]=0.2;
	}


	return;
}



//加入位置高斯白噪声
void BilateralFilter::AddPositionGaussianNoise()
{
	//求邻近点的最小值

	//加入高斯噪声
	//////////////////////////////////////////////////////////////////////////
	// 首先加入 均值为1的 
	//float noise[3];
	//long idum[3];
	//idum[0]=1;
	//idum[1]=1;
	//idum[2]=1;

	//for(int i=0;i<m_numOfPoints;i++){
	//	for(int j=0;j<3;j++){
	//		noise[j]=gasdev(&idum[j]);
	//		m_pointSets[i*3+j]=m_originalPoints[i*3+j]+noise[j]*0.05;

	//	//	m_originalPoints[i*3+j]=m_originalPoints[i*3+j];//+noise[j]*0.01;
	//	}        
	//}
	//return;
}


//加入颜色高斯白噪声
void BilateralFilter::AddColorGaussianNorse()
{
	//把颜色转化为0-1之间的浮点数
	/*for(int i=0;i<m_numOfPoints*3;i++){
		m_originalFloatColor[i]=m_originalColor[i]/255;
	}
	delete [] m_originalColor;*/
	//加入高斯白噪声
	
	//float noise[3];
	//long idum[3];
	//idum[0]=2;
	//idum[1]=2;
	//idum[2]=2;

	//for(int i=0;i<m_numOfPoints;i++){
	//	for(int j=0;j<3;j++){
	//		noise[j]=gasdev(&idum[0]);
	//		m_originalColor[i*3+j]=m_originalColor[i*3+j]+noise[j]*0.05;
	//	
	//		/*if(m_originalColor[i*3+j]<0){
	//			m_originalColor[i*3+j]=0;
	//		}
	//		else
	//			if(m_originalColor[i*3+j]>1){
	//				m_originalColor[i*3+j]=1;

	//			}*/
	//				m_colors[i*3+j]=m_originalColor[i*3+j];
	//	}        
	//}

	
	FILE *fpout; 
	float sumOfColor[3];
	for(int i=0;i<3;i++){
		sumOfColor[i]=0;
	}
	CString filename_pw = "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\noiseColor.txt";
	if((fpout = fopen(filename_pw, "w")) == NULL)
	 {
		 //MessageBox("can't open the file!");
		 int kdjk;
	 }
	 else
	 {
		 for(int i=0;i<m_numOfPoints;i++){
			 float b[3];
			 for(int j=0;j<3;j++){
				 b[j]=m_originalColor[i*3+j];
				 sumOfColor[j]+=b[j];
			 }

			 fprintf(fpout, "%f %f %f \n", b[0],b[1],b[2]);
		 }
		 float averageColor[3];
		 for(int i=0;i<3;i++){
			 averageColor[i]=sumOfColor[i]/m_numOfPoints;
		 }
		 fprintf(fpout, "%f %f %f \n", averageColor[0],averageColor[1],averageColor[2]);
		 float variationColor[3];
		 for(int i=0;i<3;i++){
			 variationColor[i]=0;
		 }
		 for(int i=0;i<m_numOfPoints;i++){
			 for(int j=0;j<3;j++)
			 {
				 variationColor[j]+=(averageColor[j]-m_originalColor[i*3+j])*(averageColor[j]-m_originalColor[i*3+j]);
			 }
		 }
		 for(int i=0;i<3;i++){
			 variationColor[i]/=m_numOfPoints;
		 }
		 fprintf(fpout, "%f %f %f \n", variationColor[0],variationColor[1],variationColor[2]);
	 }
	 fclose(fpout);


	return;
}

//////////////////////////////////////////////////////////////////////////
//计算投影点的颜色，及计算颜色的一次估计
//已知条件为 局部三角形一个 重心法向一个 局部三角形三个顶点的颜色 重心坐标一个
//首先求得局部三角形三个顶点在重心法向平面的投影，构成在该平面上的三角形
//其次求得在该平面上投影点的颜色
void BilateralFilter::ComputeProjectionPointColor(Triangle localTriangle,float normals[3],POINTVECTOR3D centroid,POINTVECTOR3D estimatePoint,Color & colorEstimate)
{
    //首先计算三个点到法向所所确定的平面上的投影
	double originalVertex[3][3];
	double vertexProjection[3][3];
	double centroidOfTriangle[3];
	double projectionPoint[3];
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			originalVertex[i][j]=m_originalPoints[localTriangle.vertexIndex[i]*3+j];
		}
		centroidOfTriangle[i]=centroid.pointVector[i];
		projectionPoint[i]=estimatePoint.pointVector[i];
	}
	//求点的投影
	for(int i=0;i<3;i++){
		//求顶点重心矢量
		float vertexCentroid[3];//原始顶点到重心的矢量
		for(int j=0;j<3;j++){
			vertexCentroid[j]=originalVertex[i][j]-centroid.pointVector[j];
		}
		//求顶点到平面的距离
		float distance;
		distance=Vector3Vector(vertexCentroid,normals);
		//求顶点位置
		for(int j=0;j<3;j++){
			vertexProjection[i][j]=originalVertex[i][j]-distance*normals[j];
		}
	}

	//////////////////////////////////////////////////////////////////////////
	//求投影点的重心坐标
	
	double centroidOrdinate[3];
    CalculateCentroidOrdinate(vertexProjection,centroidOfTriangle,projectionPoint,centroidOrdinate);

	//for(int i=0;i<3;i++){
	//	float edgeVector[3];//表示三角形边的向量
	//	float centroidToVertex[3];//重心到顶点的向量
	//	float estimatePointToVertex[3];//投影点到顶点的向量
	//	for(int j=0;j<3;j++){
	//		edgeVector[j]=vertexProjection[(int)(i+2)%3][j]-vertexProjection[(int)(i+1)%3][j];
	//		centroidToVertex[j]=centroid.pointVector[j] -vertexProjection[(int)(i+1)%3][j];
	//		estimatePointToVertex[j]=estimatePoint.pointVector[j] -vertexProjection[(int)(i+1)%3][j];
	//	}
	//	float tempVectorOne[3],tempVectorTwo[3];
	//	Cross3(edgeVector,centroidToVertex,tempVectorOne);
	//	Cross3(tempVectorOne,edgeVector,tempVectorTwo);
	//	float uniformLocalNormals[3];
	//	UnitVector(tempVectorTwo,uniformLocalNormals);
	//	centroidOrdinate[i]=Vector3Vector(estimatePointToVertex,uniformLocalNormals);		
	//}
 //   //对重心坐标进行归一化
	//float sum=0;
	//for(int i=0;i<3;i++){
	//	sum+=centroidOrdinate[i];
	//}
	////if(sum<0)
	////	MessageBox("重心坐标出错了");
	//for(int i=0;i<3;i++){
	//	centroidOrdinate[i]=centroidOrdinate[i]/sum;
	//}

	//////////////////////////////////////////////////////////////////////////
	// 求颜色
	for(int i=0;i<3;i++){
		colorEstimate.pointColor[i]=m_filterColor[localTriangle.vertexIndex[0] *3+i]*centroidOrdinate[0]
							+m_filterColor[localTriangle.vertexIndex[1]*3+i]*centroidOrdinate[1]
							+m_filterColor[localTriangle.vertexIndex[2]*3+i]*centroidOrdinate[2];
		if(colorEstimate.pointColor[i]<0)
			colorEstimate.pointColor[i]=0;
		else
			if(colorEstimate.pointColor[i]>1)
				colorEstimate.pointColor[i]=1;
	}
	return;

}

//计算滤波后的法向
void BilateralFilter::ComputeNormals()
{
	int m_nNormals = m_numOfPoints;
	m_normals=new float[m_nNormals*3];

	// Set all normals to 0.
	for (int i = 0; i < m_nNormals*3; i++) {
		m_normals[i]=0;
	}

	// Calculate normals.
	for (int i = 0; i < m_numOfTriangles; i++) {
		VECTOR3D vec1, vec2, normal;
		int id0, id1, id2;
		id0 = m_triangles[i*3];
		id1 = m_triangles[i*3+1];
		id2 = m_triangles[i*3+2];
		vec1[0] = m_pointSets[id1*3]- m_pointSets[id0*3];
		vec1[1] = m_pointSets[id1*3+1] - m_pointSets[id0*3+1];
		vec1[2] = m_pointSets[id1*3+2] - m_pointSets[id0*3+2];
		vec2[0] = m_pointSets[id2*3] - m_pointSets[id0*3];
		vec2[1] = m_pointSets[id2*3+1] - m_pointSets[id0*3+1];
		vec2[2] = m_pointSets[id2*3+2] - m_pointSets[id0*3+2];
		normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
		normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
		normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
		m_normals[id0*3+0] += normal[0];
		m_normals[id0*3+1] += normal[1];
		m_normals[id0*3+2] += normal[2];
		m_normals[id1*3+0] += normal[0];
		m_normals[id1*3+1] += normal[1];
		m_normals[id1*3+2] += normal[2];
		m_normals[id2*3+0] += normal[0];
		m_normals[id2*3+1] += normal[1];
		m_normals[id2*3+2] += normal[2];
	}

	// Normalize normals.
	for (int i = 0; i < m_nNormals; i++) {
		float length = sqrt(m_normals[i*3+0]*m_normals[i*3+0] + m_normals[i*3+1]*m_normals[i*3+1] + m_normals[i*3+2]*m_normals[i*3+2]);
		m_normals[i*3] /= length;
		m_normals[i*3+1] /= length;
		m_normals[i*3+2] /= length;
	}
	return;
}

double BilateralFilter::ComputeAreaOfTriangle(double position[3][3])
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

void BilateralFilter::CalculateCentroidOrdinate(double trianglePosition[3][3],double centroid[3],double anyPoint[3],double centroidOrdinate[3])
{
	double edgeVectorOne[3],edgeVectorTwo[3],edgeVectorThree[3],edgeVectorFour[3],areaVectorOne[3],areaVectorTwo[3];
	double areaOfTriangle=ComputeAreaOfTriangle(trianglePosition);
	double flag;
	double area;//与任意点构成的三角形的面积
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
// 首先进行中值滤波
void BilateralFilter::ComputeFilterColor()
{
	//定义一个multimap，取其系列的中间值为滤波后的值
	int medianNearest=6;
	
	std::multimap<float,int> multimapColorFilter[3];//分别针对每一种颜色滤波
    std::map<int, KnearestField>::iterator mapIterator;
	//for(int i=0;i<m_numOfPoints*3;i++){
	//	m_filterColor[i]=m_originalColor[i];
	//}
	//////////////////////////////////////////////////////////////////////////
	//一般的高斯滤波
	for(int i=0;i<m_numOfPoints;i++){
        
	}

	//////////////////////////////////////////////////////////////////////////
	//中值滤波
	for(int i=0;i<m_numOfPoints;i++){
        mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		for(int k=0;vectorIterator!=(*mapIterator).second.m_nearest.end()&& k<medianNearest;vectorIterator++,k++){
			for(int j=0;j<3;j++){
				multimapColorFilter[j].insert(std::multimap<float,int>::value_type(m_originalColor[(*vectorIterator)*3+j],(*vectorIterator)));
			}
		}
		for(int j=0;j<3;j++){
			multimapColorFilter[j].insert(std::multimap<float,int>::value_type(m_originalColor[i*3+j],i));
		}
		
		for(int j=0;j<3;j++){
			std::multimap<float,int>::iterator multimapIterator=multimapColorFilter[j].begin();
			for(int k=0;k<medianNearest/2;k++){
				multimapIterator++;
			}
			m_filterColor[i*3+j]=m_originalColor[(*multimapIterator).second*3+j];		

		}
		for(int j=0;j<3;j++){
			multimapColorFilter[j].clear();
			
		}
	}
	delete[] m_originalColor;
	//////////////////////////////////////////////////////////////////////////
	
	for(int i=0;i<m_numOfPoints*3;i++){
		m_colors[i]=m_filterColor[i];
	}

	CString filename_pw= "F:\\marchingcube\\marchingcube\\marchingcube\\ceshi\\resultColor.txt";
	FILE *fpout; 
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		//MessageBox("can't open the file!");
		int kdjk;
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			float b[3];
			for(int j=0;j<3;j++){
				b[j]=m_filterColor[i*3+j];
			}

			fprintf(fpout, "%f %f %f \n", b[0],b[1],b[2]);
		}
	}
	fclose(fpout);
	return;
}

//估计颜色变化的梯度及原点的颜色值
void BilateralFilter::EstimateColorGradient(Triangle triangleIndex ,POINTVECTOR3D centroidNormals,float gradientColor[4],float attributeValue[3])
{
	mat_f8 m_matrix(4,4);
	vec_f8 m_vector(4);
	vec_f8 m_gradient(4);
	// 计算三角形顶点到法平面上的投影点
	float vertexProjection[3];
	float tempVertex[3];
	float tempNormal[3];
	for(int i=0;i<3;i++){
		tempNormal[i]=centroidNormals.pointVector[i];
	}

	for(int i=0;i<3;i++){
		float vertexProjection[3];
		float tempVertex[3];
		float projection_distance;
		for(int j=0;j<3;j++){
			tempVertex[j]=m_originalPoints[triangleIndex.vertexIndex[i]*3+j];
		}
		projection_distance=Vector3Vector(tempVertex,tempNormal);
		
		
		for(int j=0;j<3;j++){
			//加入投影点；
			
			m_matrix(i, j)=tempVertex[j]-projection_distance*tempNormal[j];//m_originalPoints[triangleIndex.vertexIndex[i]*3+j];
		}
	}
	for(int i =0;i<3;i++){
		m_matrix(3,i)=centroidNormals.pointVector[i];
	}
	m_matrix(0,3)=1;m_matrix(1,3)=1;
	m_matrix(2,3)=1;m_matrix(3,3)=0;

	for(int i=0;i<3;i++){
		m_vector(i)=attributeValue[i];
	}
	m_vector(3)=0;

	if(matrix_inverse(m_matrix)){
		m_gradient=m_matrix*m_vector;
		for(int i=0;i<4;i++){
			gradientColor[i]=m_gradient(i);

		}
	}
	else{
		for(int i=0;i<4;i++){
			gradientColor[i]=0;
		}
	}

	return;    
}

// 滤波颜色变化的梯度及原点的颜色值
void BilateralFilter::CalculateFilterColorGradient()
{
	float variation=0.6;
	for(int i=0;i<m_numOfPoints;i++){	

		float tempSumWeight=0;
		
		//求局部点
		if(!m_LocalPoints.empty())
			m_LocalPoints.clear();
		std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(i);		
		std::vector<int>::iterator mapNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
		m_LocalPoints.push_back(i);//把当前点作为第一个点放入
		for(;mapNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();mapNearestIterator++){
			m_LocalPoints.push_back((*mapNearestIterator));
		}
	//局部三角化
		ComputeLocalTriangles();
	//估计每个三角片的颜色梯度
		std::vector<Triangle>::iterator triangleIterator=m_LocalTriangles.begin();
		for(;triangleIterator!=m_LocalTriangles.end();triangleIterator++){
			float localTriangleCentroid[3];
			float area_triangle;
			float tempWeight;
			float distance_centroid_point;
			double vertexPosition[3][3];
			for(int j=0;j<3;j++){
				for(int k=0;k<3;k++){
					vertexPosition[j][k]=m_originalPoints[(*triangleIterator).vertexIndex[j]*3+k];
				}
				
			}
			for(int j=0;j<3;j++){
                localTriangleCentroid[j]=(vertexPosition[0][j]+vertexPosition[1][j]+vertexPosition[2][j])/3;
			}
			area_triangle=ComputeAreaOfTriangle(vertexPosition);
			distance_centroid_point=((localTriangleCentroid[0]-m_originalPoints[i*3])*(localTriangleCentroid[0]-m_originalPoints[i*3])
				             +(localTriangleCentroid[1]-m_originalPoints[i*3+1])*(localTriangleCentroid[1]-m_originalPoints[i*3+1])
							 +(localTriangleCentroid[2]-m_originalPoints[i*3+2])*(localTriangleCentroid[2]-m_originalPoints[i*3+2]));
			tempWeight=area_triangle*exp(-0.5*distance_centroid_point/variation);
			//
			POINTVECTOR3D tempCentroidNormal;
			float tempGradient[3][4];//三种颜色的梯度值和原点颜色值
			float tempAttribute[3][3];//三个颜色值三个顶点

			for(int j=0;j<3;j++){
				tempCentroidNormal.pointVector[j]=m_originalFilterNormals[i*3+j];
			}
			for(int j=0;j<3;j++){// j 代表颜色种类
				for(int k=0;k<3;k++){ //k 代表顶点
					tempAttribute[j][k]=m_originalColor[(*triangleIterator).vertexIndex[k]*3+j];
				}
				EstimateColorGradient((*triangleIterator),tempCentroidNormal,tempGradient[j],tempAttribute[j]);
				for(int k=0;k<4;k++){	// k 代表三个坐标轴和一个原点				
					m_colorGradient[j][i*4+k]+=tempWeight*tempGradient[j][k];
				}
			
			}	
			tempSumWeight+=tempWeight;
			
		
		}
		//估计每个点的颜色梯度和原点颜色
		for(int k=0;k<3;k++){
			for(int j=0;j<4;j++){
				m_colorGradient[k][i*4+j]/=tempSumWeight;
			}

		}
		
		m_LocalPoints.clear();
		m_LocalTriangles.clear();
	}
	//滤波每一点的颜色梯度
	for(int i=0;i<3;i++){
		m_filterColorGradient[i]=new float[m_numOfPoints*4];
	}
	for(int i=0;i<3;i++){
		for(int j=0;j<m_numOfPoints;j++){
			for(int k=0;k<4;k++){
				m_filterColorGradient[i][j*4+k]=0;
				//m_filterColorGradient[i][j*4+k]=m_colorGradient[i][j*4+k];
			}
		}
	}
	
	for(int i=0;i<m_numOfPoints;i++){
		float sumweight,sumdistance,distanceVariance;
		sumweight=0;
		sumdistance=0;
		float* distance;

		std::map<int,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		distance=new float[(*mapIterator).second.m_nearest.size()];

		for(int j=0;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,j++){

			distance[j]=(m_originalPoints[i*3]-m_originalPoints[(*vectorIterator)*3])*(m_originalPoints[i*3]-m_originalPoints[(*vectorIterator)*3])
				+(m_originalPoints[i*3+1]-m_originalPoints[(*vectorIterator)*3+1])*(m_originalPoints[i*3+1]-m_originalPoints[(*vectorIterator)*3+1])
				+(m_originalPoints[i*3+2]-m_originalPoints[(*vectorIterator)*3+2])*(m_originalPoints[i*3+2]-m_originalPoints[(*vectorIterator)*3+2]);
			sumdistance+=distance[j];
		}
		distanceVariance=sumdistance/((*mapIterator).second.m_nearest.size());
		vectorIterator=(*mapIterator).second.m_nearest.begin();
		for(int j=0;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,j++){
			float weight;
			weight=exp(-0.5*distance[j]/distanceVariance);
			for(int k=0;k<4;k++){
				m_filterColorGradient[0][i*4+k]+=weight*m_colorGradient[0][(*vectorIterator)*4+k];
				m_filterColorGradient[1][i*4+k]+=weight*m_colorGradient[1][(*vectorIterator)*4+k];
				m_filterColorGradient[2][i*4+k]+=weight*m_colorGradient[2][(*vectorIterator)*4+k];
			}
			sumweight+=weight;

		}
		for(int j=0;j<4;j++){
			m_filterColorGradient[0][i*4+j]+=exp(0.0)*m_colorGradient[0][i*4+j];
			m_filterColorGradient[1][i*4+j]+=exp(0.0)*m_colorGradient[1][i*4+j];
			m_filterColorGradient[2][i*4+j]+=exp(0.0)*m_colorGradient[2][i*4+j];


		}
		sumweight+=exp(0.0);

		for(int j=0;j<4;j++){
			m_filterColorGradient[0][i*4+j]/=sumweight;
			m_filterColorGradient[1][i*4+j]/=sumweight;
			m_filterColorGradient[2][i*4+j]/=sumweight;
		}
		delete[] distance;
	}

	for(int i=0;i<3;i++){
		for(int j=0;j<m_numOfPoints;j++){
			for(int k=0;k<4;k++){
				m_colorGradient[i][j*4+k]=m_filterColorGradient[i][j*4+k];
			}
		}
	}
	
	//delete[] m_filterColorGradient[0];
	//delete[] m_filterColorGradient[1];
	//delete[] m_filterColorGradient[2];
	return;
	
}

// mean shift 滤波法向
void BilateralFilter:: MeanShiftFilterNormal()
{
    float distance_position;
	float distance_normal;
	float position_kernel,range_kernel;
	float sum_kernel;
	float length_normal_vector;
	int num_stop=0;//已经停止的数
	bool* flag_stop=new bool[m_numOfPoints];
	for(int i=0;i<m_numOfPoints;i++){
		flag_stop[i]=false;
	}

    m_originalFilterNormals=new float[m_numOfPoints*3];
	while(num_stop<m_numOfPoints){
		for(int i=0;i<m_numOfPoints;i++){
			if(!flag_stop[i]){
				sum_kernel=0;
				for(int j=0;j<3;j++){
                    m_originalFilterNormals[i*3+j]=0;
				}
			
				std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(i);
				//int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
				std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
				for(;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++){
					distance_position=(m_pointSets[i*3]-m_pointSets[(*vectorNearestIterator)*3])*(m_pointSets[i*3]-m_pointSets[(*vectorNearestIterator)*3])
						+(m_pointSets[i*3+1]-m_pointSets[(*vectorNearestIterator)*3+1])*(m_pointSets[i*3+1]-m_pointSets[(*vectorNearestIterator)*3+1])
						+(m_pointSets[i*3+2]-m_pointSets[(*vectorNearestIterator)*3+2])*(m_pointSets[i*3+2]-m_pointSets[(*vectorNearestIterator)*3+2]);
					distance_normal=(m_originalNormals[i*3]-m_originalNormals[(*vectorNearestIterator)*3])*(m_originalNormals[i*3]-m_originalNormals[(*vectorNearestIterator)*3])
						+(m_originalNormals[i*3+1]-m_originalNormals[(*vectorNearestIterator)*3+1])*(m_originalNormals[i*3+1]-m_originalNormals[(*vectorNearestIterator)*3+1])
						+(m_originalNormals[i*3+2]-m_originalNormals[(*vectorNearestIterator)*3+2])*(m_originalNormals[i*3+2]-m_originalNormals[(*vectorNearestIterator)*3+2]);
					position_kernel=exp(-0.5*distance_position/m_meanShiftHposition);
					range_kernel=exp(-0.5*distance_normal/m_meanShiftHnormal);
					sum_kernel+=position_kernel*range_kernel;
					
					for(int j=0;j<3;j++){
						m_originalFilterNormals[i*3+j]+=m_originalNormals[(*vectorNearestIterator)*3+j]*position_kernel*range_kernel;
					}
				}
				sum_kernel+=exp(0.0)*exp(0.0);
				for(int j=0;j<3;j++){
					m_originalFilterNormals[i*3+j]+=m_originalNormals[i*3+j]*exp(0.0)*exp(0.0);
				}
				for(int j=0;j<3;j++){
					m_originalFilterNormals[i*3+j]/=sum_kernel;
				}
				//归一化 m_originalFilterNormals
				float tempNormal[3];
				for(int j=0;j<3;j++){
					tempNormal[j]=m_originalFilterNormals[i*3+j];
				}
				length_normal_vector=sqrt(Vector3Vector(tempNormal,tempNormal));
				for(int j=0;j<3;j++){
					m_originalFilterNormals[i*3+j]/=length_normal_vector;
				}

				length_normal_vector=((m_originalFilterNormals[i*3]-m_originalNormals[i*3])*(m_originalFilterNormals[i*3]-m_originalNormals[i*3])
					+(m_originalFilterNormals[i*3+1]-m_originalNormals[i*3+1])*(m_originalFilterNormals[i*3+1]-m_originalNormals[i*3+1])
					+(m_originalFilterNormals[i*3+2]-m_originalNormals[i*3+2])*(m_originalFilterNormals[i*3+2]-m_originalNormals[i*3+2]));
				
				if(((m_originalFilterNormals[i*3]-m_originalNormals[i*3])*(m_originalFilterNormals[i*3]-m_originalNormals[i*3])
					+(m_originalFilterNormals[i*3+1]-m_originalNormals[i*3+1])*(m_originalFilterNormals[i*3+1]-m_originalNormals[i*3+1])
					+(m_originalFilterNormals[i*3+2]-m_originalNormals[i*3+2])*(m_originalFilterNormals[i*3+2]-m_originalNormals[i*3+2]))<m_meanShiftStopNormal){
						flag_stop[i]=true;
						num_stop+=1;
					}
				else{
					m_originalNormals[i*3]=m_originalFilterNormals[i*3];
					m_originalNormals[i*3+1]=m_originalFilterNormals[i*3+1];
					m_originalNormals[i*3+2]=m_originalFilterNormals[i*3+2];
				}
				
			}
		}

	}
	delete[] m_originalNormals;

}
// mean shift 滤波颜色梯度
void BilateralFilter::MeanShiftFilterColorGradient()
{
	float distance_position;
	float distance_normal;
	float distance_color[3];
	float distance_color_gradient[3];
	float position_kernel,normal_kernel,color_kernel[3],color_gradient_kernel[3];
	float sum_kernel[3];
	//float length_normal_vector;
	int num_stop=0;//已经停止的数
	bool* flag_stop[3];
	flag_stop[0]=new bool[m_numOfPoints];
	flag_stop[1]=new bool[m_numOfPoints];
	flag_stop[2]=new bool[m_numOfPoints];
	for(int i=0;i<m_numOfPoints;i++){
		flag_stop[0][i]=false;
		flag_stop[1][i]=false;
		flag_stop[2][i]=false;
	}
    
	m_filterColorGradient[0]=new float[m_numOfPoints*4];
	m_filterColorGradient[1]=new float[m_numOfPoints*4];
	m_filterColorGradient[2]=new float[m_numOfPoints*4];
	
	while(num_stop<m_numOfPoints*3){
		for(int i=0;i<m_numOfPoints;i++){
			if(!flag_stop[0][i]||!flag_stop[1][i]||!flag_stop[2][i]){
				sum_kernel[0]=0;sum_kernel[1]=0;sum_kernel[2]=0;
				for(int k=0;k<3;k++){
					if(!flag_stop[k][i]){
						for(int j=0;j<4;j++){
							m_filterColorGradient[k][i*4+j]=0;
						}
					}
				}
				std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(i);
				//int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
				std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
				for(;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++){
					distance_position=(m_pointSets[i*3]-m_pointSets[(*vectorNearestIterator)*3])*(m_pointSets[i*3]-m_pointSets[(*vectorNearestIterator)*3])
						+(m_pointSets[i*3+1]-m_pointSets[(*vectorNearestIterator)*3+1])*(m_pointSets[i*3+1]-m_pointSets[(*vectorNearestIterator)*3+1])
						+(m_pointSets[i*3+2]-m_pointSets[(*vectorNearestIterator)*3+2])*(m_pointSets[i*3+2]-m_pointSets[(*vectorNearestIterator)*3+2]);
					distance_normal=(m_originalFilterNormals[i*3]-m_originalFilterNormals[(*vectorNearestIterator)*3])*(m_originalFilterNormals[i*3]-m_originalFilterNormals[(*vectorNearestIterator)*3])
						+(m_originalFilterNormals[i*3+1]-m_originalFilterNormals[(*vectorNearestIterator)*3+1])*(m_originalFilterNormals[i*3+1]-m_originalFilterNormals[(*vectorNearestIterator)*3+1])
						+(m_originalFilterNormals[i*3+2]-m_originalFilterNormals[(*vectorNearestIterator)*3+2])*(m_originalFilterNormals[i*3+2]-m_originalFilterNormals[(*vectorNearestIterator)*3+2]);
					position_kernel=exp(-0.5*distance_position/m_meanShiftHposition);
					normal_kernel=exp(-0.5*distance_normal/m_meanShiftHnormal);
					for(int j=0;j<3;j++){
						if(!flag_stop[j][i]){
							float temp_kernel;
							distance_color[j]=(m_originalColor[i*3+j]-m_originalColor[(*vectorNearestIterator)*3+j])*(m_originalColor[i*3+j]-m_originalColor[(*vectorNearestIterator)*3+j]);
							distance_color_gradient[j]=(m_colorGradient[j][i*4]-m_colorGradient[j][(*vectorNearestIterator)*4])*(m_colorGradient[j][i*4]-m_colorGradient[j][(*vectorNearestIterator)*4])
								+(m_colorGradient[j][i*4+1]-m_colorGradient[j][(*vectorNearestIterator)*4+1])*(m_colorGradient[j][i*4+1]-m_colorGradient[j][(*vectorNearestIterator)*4+1])
								+(m_colorGradient[j][i*4+2]-m_colorGradient[j][(*vectorNearestIterator)*4+2])*(m_colorGradient[j][i*4+2]-m_colorGradient[j][(*vectorNearestIterator)*4+2])
								+(m_colorGradient[j][i*4+3]-m_colorGradient[j][(*vectorNearestIterator)*4+3])*(m_colorGradient[j][i*4+3]-m_colorGradient[j][(*vectorNearestIterator)*4+3]);
							color_kernel[j]=exp(-0.5*distance_color[j]/m_meanShiftHcolor[j]);
							color_gradient_kernel[j]=exp(-0.5*distance_color_gradient[j]/m_meanShiftHcolorGradient[j]);
							temp_kernel=position_kernel*normal_kernel*color_kernel[j]*color_gradient_kernel[j];
							sum_kernel[j]+=temp_kernel;
							for(int k=0;k<4;k++){
								m_filterColorGradient[j][i*4+k]+=m_colorGradient[j][(*vectorNearestIterator)*4+k]*temp_kernel;
							}
						}
					}
					
				}
   
				for(int j=0;j<3;j++){
					if(!flag_stop[j][i]){
						float temp_kernel=exp(0.0)*exp(0.0)*exp(0.0)*exp(0.0);
						sum_kernel[j]+=temp_kernel;
						for(int k=0;k<4;k++){
							m_filterColorGradient[j][i*4+k]+=m_colorGradient[j][i*4+k]*temp_kernel;
							m_filterColorGradient[j][i*4+k]/=sum_kernel[j];
						}
						if(((m_filterColorGradient[j][i*4]-m_colorGradient[j][i*4])*(m_filterColorGradient[j][i*4]-m_colorGradient[j][i*4])
							+(m_filterColorGradient[j][i*4+1]-m_colorGradient[j][i*4+1])*(m_filterColorGradient[j][i*4+1]-m_colorGradient[j][i*4+1])
							+(m_filterColorGradient[j][i*4+2]-m_colorGradient[j][i*4+2])*(m_filterColorGradient[j][i*4+2]-m_colorGradient[j][i*4+2])
							+(m_filterColorGradient[j][i*4+3]-m_colorGradient[j][i*4+3])*(m_filterColorGradient[j][i*4+3]-m_colorGradient[j][i*4+3]))<m_meanShiftStopColorGradient[j]){
								flag_stop[j][i]=true;
								num_stop+=1;
							}
						else{
							for(int k=0;k<4;k++){
								m_colorGradient[j][i*4+k]=m_filterColorGradient[j][i*4+k];
							}
						}
					}

				}				
			}
		}

	}
	delete[] m_colorGradient[0];
	delete[] m_colorGradient[1];
	delete[] m_colorGradient[2];
}