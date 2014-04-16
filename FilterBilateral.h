// FilterBilateral.h: interface for the FilterBilateral class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FILTERBILATERAL_H__0A901A27_1C99_49CA_ADE4_9D10C4B19C69__INCLUDED_)
#define AFX_FILTERBILATERAL_H__0A901A27_1C99_49CA_ADE4_9D10C4B19C69__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "NormalStruct.h"



//struct KnearestField{
//	unsigned int m_numOfNearest;//邻域的数量，必须小于k
//	bool m_IfBoundary;//是否是边界 如果是为 true，否则为false
//	std::vector<int> m_nearest;//存放邻域点
//};// 保存点的邻域的结构
//
//struct Triangle{
//	 int vertexIndex[3];
//};
//
//struct Color {
//	float pointColor[3];
//};
//
//struct POINTVECTOR3D {
//	float pointVector[3];
//};


class FilterBilateral  
{
public:
	FilterBilateral();
	virtual ~FilterBilateral();
protected:
	//输入数据
	float* m_originalPointSet;
	float* m_originalColors;
	float* m_originalNormals;
public:
	//输出数据
	int m_numOfPoints;
	float* m_resultPointSet;
	float* m_resultColors;
	float* m_resultNormals;
	long nearestTime;
	long filterTime;
protected:
	//过程数据
	std::map<int, KnearestField> m_mapKnearest;	
	float m_radius;
	int m_knearest;
	float  m_meanShiftHposition,m_meanShiftHnormal;
	float m_meanShiftStopNormal;
	float* m_filterNormals;
	float* m_originalColorGradients[3];
	float* m_resultColorGradients[3];
	float m_thresholdColor[3];
protected:
	//在计算双边滤波时候用到的数据
	std::vector<POINTVECTOR3D> m_localPositionEstimation;//局部位置的一阶估计
	std::vector<float> m_localPointDistance;//点与点之间的距离，就是双边滤波中的第一式子
	std::vector<float> m_localProjectionDistance;//滤波点到投影点之间的距离
	std::vector<float> m_localColorDistance[3];//颜色距离
	std::vector<float> m_localColorProjectionDistance[3];//投影颜色的距离
	std::vector<Color> m_localColorEstimation;//局部颜色的一次估计
	
	float m_pointVariation;//原始点的方差
	float m_estimationVariation;//投影距离的方差
	float m_colorVariation[3];//颜色的方差 
	float m_colorEstimationVariation[3];//颜色投影的方差
public:
	void GetFilterBilateral(int numOfPoints, float* pointSet, float* colors);
	void DeleteFilterBilateral();
protected:
	void CalculateVariation();
	void ComputeMapKnearest();
	void CalculateOriginalNormals();
	void MeanShiftNormals();
	void CalculateColorGradient();
	void CalculateBilateralFilterOne();
	void CalculateBilateralFilterTwo();
	void CalculateBilateralFilterThree();
	void CalculateBilateralFilterFour();
	void CalculateBilateralFilterFive();
	void CalculateBilateralFilterSix();//pauly 没有第二项 （投影距离加权项）  只有几何滤波
	void CalculateBilateralFilterSeven();// 没有meanshift 法向项   只有几何滤波
	void CalculateBilateralFilterEight();//有meanshift 项  只有几何滤波

	void ComputeOneOrderPositionEstimation(int indexPoint);
	void ComputeOneOrderColorEstimation(int indexPoint);
	//加入位置高斯噪声
	void AddPositionGaussianNoise();
	//加入颜色高斯噪声
	void AddColorGaussianNorse();
	void CalculateColorGradientFair();
			
};

#endif // !defined(AFX_FILTERBILATERAL_H__0A901A27_1C99_49CA_ADE4_9D10C4B19C69__INCLUDED_)
