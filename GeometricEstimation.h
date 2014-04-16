#pragma once
//#include "GeometricFlow.h"
#include "Vectors.h"
//#include "MZsplating.h"

struct PointGeometry{
	//UINT8 m_ID;
	//POINT3D m_pointPosition;
	float m_gaussCurvature;
	float m_averageCurvature;
	VECTOR3D m_normal;
};//关于点的几何特征的结构

class GeometricEstimation
{
public:
	//当有原始法向时候求几何特征的构造函数
	GeometricEstimation(int numOfPointSet,float* pointSet,uint8 numOfKnearest,int* allPointsNearestNeighbor,float* allPointsNormal);
	//求邻域时的构造函数
	GeometricEstimation(int numOfPointSet,float* pointSet);
	//当没有原始法向时候求几何特征时的构造函数
	GeometricEstimation(int numOfPointSet,float* pointSet,UINT8 numOfKnearest,int* allPointsNearestNeighbor);
	~GeometricEstimation(void);
public:
	// 集合概念下的点集 邻域 法向 几何特征
	float* m_allPointsNormal;//所有点的法向
	int* m_allPointsNeighbor;//所有点的邻域
    float* m_pointSet;// 点集
	float* m_allPointsPricipalCurvatureOne;
	float* m_allPointsPricipalCurvatureTwo;
	float* m_allPointsPricipalDirectionOne;
	float* m_allPointsPricipalDirectionTwo;
	float* m_allPointsGaussCurvature;
	float* m_allPointsAverageCurvature;
	PointGeometry* m_allPointGeometry;//所估计点的几何特征
protected:    
	//单个点的邻域 法向 几何特征
    POINT3D m_position;//所估计点的位置
    int m_ID;//所估计点的ID
    int* m_pointNeighbor;//所估计点的邻域
    float m_averageCurvature;//所估计点的平均曲率
    PointGeometry m_pointGeometry;//所估计点的几何特征
	float m_gaussCurvature;//所估计点的高斯曲率
	float m_principalCurvature[2];//所估计点的主曲率
	VECTOR3D m_principalDirection[2];//所估计的两个主方向
    VECTOR3D m_normal;//所估计点的法向
    
	//共同参数
	float m_radius;//邻域半径
	int m_numOfPoint;//点的数目
	UINT8 m_numOfNeighbor;//每一个点的邻域数目
	
	//uint8 m_knearest;//k近邻的k
	//float m_radius;//搜索k近邻的球半径
protected:
	//法向估计
	void NormalEstimation();
	//曲率估计
	void CurvatureEstimation();
public:
	//初始化
	void SetInput(UINT8 ID);
	//释放空间
	void Delete();
	//邻域搜索
	void SearchNeighbor(UINT8 knearest,float radius,int* pointNeighbor);
	//确定近邻搜索的半径
	float EstimateSearchingRadius();
	//返回法向
	void GetNormal(VECTOR3D normal);
	//返回平均曲率
	float GetAverageCurvature();
	//返回高斯曲率
	float GetGaussCurvature();
	//返回几何特征（包括上面三项）
	void GetPointGeometricCharacter(PointGeometry pointGeometry);
	//返回用于rendering的椭圆
	MZSplatEllipse GetSplatEllipse();
	//////////////////////////////////////////////////////////////////////////
	//法向量全局化 
	void EstimateNormalDirection(int numOfPoints,float* pointSet,uint8 numOfEachPointNeighbor,int* allPointsNeighbor,float* allPointsNormal);
	//得到所有点的法向(没有初始法向)
	void GetAllPointsNormalFirst();
	//得到所有点的法向(有初始法向)
	void GetAllPointsNormalSecond();
	//得到所有点的几何特征
	void GetAllPointsGeometricChracter();
	void GetAllPointsNeighbor(uint8 knearest,float radius);
	void GetAllPointsGaussAverageCurvature();

};
