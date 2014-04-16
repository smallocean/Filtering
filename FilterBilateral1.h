// FilterBilateral.h: interface for the FilterBilateral class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FILTERBILATERAL_H__0A901A27_1C99_49CA_ADE4_9D10C4B19C69__INCLUDED_)
#define AFX_FILTERBILATERAL_H__0A901A27_1C99_49CA_ADE4_9D10C4B19C69__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "NormalStruct.h"
#include <functional>
#include <algorithm>
#include <cmath>
//using namespace std;



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


class FilterBilateral1  
{
public:
	FilterBilateral1();
	virtual ~FilterBilateral1();
protected:
	//输入数据
	float* m_originalPointSet;
	float* m_originalColors;
	
public:
	//输出数据
	int m_numOfPoints;
	float* m_resultPointSet;
	float* m_resultColors;
	float* m_resultNormals;
	long nearestTime;
	long filterTime;
	float* m_originalNormals;
	std::map<int, KnearestField> m_mapKnearest;
protected:
	//过程数据
	//std::map<int, KnearestField> m_mapKnearest;	
	float m_radius;
	int m_knearest;
	float  m_meanShiftHposition,m_meanShiftHnormal;
	float m_meanShiftStopNormal;
	float* m_filterNormals;
	float* m_originalColorGradients[3];
	float* m_resultColorGradients[3];
	float m_thresholdColor[3];
	float m_functionWeight;
	float m_thresholdDistanceNormal;
	float m_thresholdDistanceGradientFunction;
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
	void GetFilterBilateral(int numOfPoints, float* pointSet,float * colors,float meanLength,float variationF,float variationG1,float variationG2,float variationH1,float variationH2,float variationH3,int knearest,int meanShift,int interactiveTimes,int indexFunction,float variationNstop,float functionVariation,float thresholdOfColor,float functionGradientWide,float functionGradientThreshold,float thresholdDistanceNormal,float thresholdDistanceGradientFunction);
	void DeleteFilterBilateral();
protected:
	void CalculateVariation();
	void ComputeMapKnearest();
	void CalculateOriginalNormals();
	void MeanShiftNormals();
	void CalculateColorGradient();
	void CalculateBilateralFilterOne(int interactiveTimes);
	void CalculateBilateralFilterTwo(int interactiveTimes);
	void CalculateBilateralFilterThree(int interactiveTimes);
	void CalculateBilateralFilterFour(int interactiveTimes);
	void CalculateBilateralFilterFive(int interactiveTimes);
	void CalculateBilateralFilterSix(int interactiveTimes);//pauly 没有第二项 （投影距离加权项）  只有几何滤波
	void CalculateBilateralFilterSeven(int interactiveTimes);// 没有meanshift 法向项   只有几何滤波
	void CalculateBilateralFilterEight(int interactiveTimes);//有meanshift 项  只有几何滤波
	void CalculateBilateralFilterTen(int interactiveTimes);//重心估计；
	void CalculateBilateralFilterEleven(int interactiveTimes);//重心估计；
	void CalculateBilateralFilterTwelve(int interactiveTimes);

	void ComputeOneOrderPositionEstimation(int indexPoint);
	void ComputeOneOrderColorEstimation(int indexPoint);
	//加入位置高斯噪声
	void AddPositionGaussianNoise();
	//加入颜色高斯噪声
	void AddColorGaussianNorse();
	void CalculateColorGradientFair();
	void MidianFilter();
	void MeanShiftColorGradient();		 
	void EstimateNormalDirection();
	void ComputerAverageRadius();
			
private:
	float m_gradientKernelWide;
	float m_thresholdColorGradient;
	float m_meanRadius;
};


struct lessFloat: public std::binary_function <float, float, bool> 
{
	bool operator()(const float& _Left, const float& _Right) const
	{
		if (std::fabs(_Left-_Right) < 1e-9)
			return false;
		else
			return _Left < _Right;
	}
};
#endif // !defined(AFX_FILTERBILATERAL_H__0A901A27_1C99_49CA_ADE4_9D10C4B19C69__INCLUDED_)
