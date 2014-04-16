/*********************************************************/
/*filename:    BilateralFilter.h                           */

/*description: This class is designed for point and mesh filtering */

/*author:      redstar    data: 2006/07/01                */
/*********************************************************/
#pragma once
#include <map>
#include <vector>
#include "Vectors.h"
//#include "mesh.h"
#include "NormalStruct.h"

class BilateralFilter
{
public:
	BilateralFilter(void);
	~BilateralFilter(void);
protected:
	//输入数据
	float* m_originalPoints;//原始点
	float* m_originalNormals;// 原始点的法向
	float* m_originalFilterNormals;//把原始点的法向滤波后的法向
	float* m_filterColorGradient[3];
//	unsigned int* m_originalColor;//把原来点的颜色转换为0---1之间
	float* m_originalColor;//原始点的颜色 
	float* m_filterColor;//一次滤波后的颜色（预处理）
	float* m_colorGradient[3];//顶点的颜色梯度和原点的估计值
	float  m_meanShiftHposition,m_meanShiftHnormal,m_meanShiftHcolorGradient[3];//用meanshift 滤波时候的核宽度；
	float m_meanShiftHcolor[3];
	float m_meanShiftStopNormal,m_meanShiftStopColorGradient[3];//mean Shift 停止标准


 //   float* m_originalNoiseColor;//把颜色加噪声//

protected:
	//中间数据
	int m_knearest;// k邻域的大小
    float m_radius;// k 邻域的半径大小
	std::map<int, KnearestField> m_mapKnearest;	
	//滤波时候用到的数据
	 int m_numOFLocalPoints;//局部点集的数量
	std::vector<int> m_LocalPoints;//局部点集合，其中第一个点为要滤波的点
	 int m_numOfLocalTriangles;//局部三角化后的三角片数
	std::vector<POINTVECTOR3D> m_LocalTriangleCentroid;//局部三角形的重心
	std::vector<POINTVECTOR3D> m_originalCentroidNormals;////局部三角形的重心法向，由三个顶点法向加权	
	std::vector<POINTVECTOR3D> m_LocalPositionEstimation;//局部位置的一阶估计
	std::vector<float> m_localPointDistance;//点与点之间的距离，就是双边滤波中的第一式子
	std::vector<float> m_localProjectionDistance;//滤波点到投影点之间的距离
	std::vector<float> m_localColorDistance[3];//颜色距离
	std::vector<Color> m_LocalColorEstimation;//局部颜色的一次估计
	std::vector<Triangle> m_LocalTriangles;//局部三角化 三个点三个点存放
	float m_pointVariation;//原始点的方差
	float m_estimationVariation;//投影距离的方差
	float m_clolrVariation[3];//颜色的方差
	
public:
	//输出数据
	int m_numOfPoints;//点的数量
	int m_numOfTriangles;//三角片的数量
	 int* m_triangles;// 三角片集合
	float* m_pointSets;//滤波后的点
	float* m_normals;//滤波后点的法向
	float* m_colors;//最后的颜色
public:
	void GetBilateralFilter(int numOfPointSets,float* pointSets, int numOfTriangle, int* triangles,float* vertexNormals,float* vertexColors);
    void DeleteBilateralFilter();
protected:
    //求点的局部邻域，对三角形来说，就是求它的2－ring邻域
	void ComputeMapKnearest();

	//求点的局部三角化
	void ComputeLocalTriangles();

	//求位置和颜色的一次预测
	void ComputeOneOrderEstimation();

	//求位置和颜色的1.5次预测
	void ComputeOneHalfEstimation( int indexPoint);

	//求双边滤波结果；
	void ComputeBilateralFilter( int indexPoint);

	//求法向的双边滤波（molification)
	void ComputeNormalFilter( int indexPoint);
	
   //求滤波的函数的方差值 以及滤波中用到的距离的平方
	void ComputeVariation(int indexPoint);

	//计算三角形的面积
	double ComputeAreaOfTriangle( double position[3][3]);

	//加入位置高斯噪声
	void AddPositionGaussianNoise();

	//加入颜色高斯噪声
	void AddColorGaussianNorse();
    
	//由三点颜色计算投影点颜色
	void ComputeProjectionPointColor(Triangle localTriangle,float normals[3],POINTVECTOR3D centroid,POINTVECTOR3D estimatePoint,Color & colorEstimate);

	//计算最滤波后点的法向
	void ComputeNormals();

	//计算滤波后的颜色（预处理）
	void ComputeFilterColor();

	//计算一点的重心坐标
	void CalculateCentroidOrdinate(double trianglePosition[3][3],double centroid[3],double anyPoint[3],double centroidOrdinate[3]);

	//估计颜色变化的剃度及原点的颜色值
	void EstimateColorGradient(Triangle triangleIndex,POINTVECTOR3D centroidNormals,float gradientColor[4],float attributeValue[3]);

	// 滤波颜色变化的剃度及原点的颜色值
	void CalculateFilterColorGradient();


	// mean shift 滤波法向
	void MeanShiftFilterNormal();
    // mean shift 滤波颜色梯度
	void MeanShiftFilterColorGradient();
};
