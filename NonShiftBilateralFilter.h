#pragma once


#include "NormalStruct.h"
#include <functional>
#include <algorithm>
#include <cmath>
#include "FilterBilateral1.h"

class NonShiftBilateralFilter
{
public:
	NonShiftBilateralFilter(void);
	~NonShiftBilateralFilter(void);
protected:
	double* m_originalPointSet;
	double* m_originalColors;
	double m_minPosition,m_maxPosition;
	double* m_normOriginalPointSetAndColor;
	double* m_meanRadius;
	double* m_meanAreas;
	int m_numOfPoints;
	int m_kNearest;
	int m_IfFilterNormal;
	int m_ifVolumePreserve;
	int m_indexFunction;
	int m_tangentOrManifold;
	int m_ifNormalWeigtht;
	int m_ifAreaWeight;
	double m_meanShiftStopNormal;
	int m_ifVariationNormal;
	
	double* m_originalNormalColor;	
	double* m_filterNormals;
	double* m_resultNormalColor;	
	double* m_resultPointSetAndColor;
	double* m_volumeChange;
	double* m_distanceMove;
	double* m_meanCurvature;
	double m_radius;
	int m_interactiveTime;
	double m_meanShiftHnormal;
	double m_meanShiftHposition;
	double m_variationTangent;
	double m_variationNormal;
	double m_variationColor;
	double m_thresholdFunction;
	double m_thresholdGradient;
	float m_thresholdDistanceNormal;
	float m_thresholdDistanceGradientFunction;
	double m_timeStep;
	
	int m_numDimension;
	int m_numDimensionColor;
public:
	double* m_resultPointSet;
	double* m_resultColor;
	double* m_originalNormals;
	long nearestTime;
	long filterTime;
	std::map<int, KnearestField> m_mapKnearest;
public:
	void GetNonShiftBilateralFilter(int numOfPointSet,float* pointSet,int dimensionColor,float* colorSet,int kNearest, int ifFilterNormal, int ifPreserveVolume, int indexFunction,double thresholdNormal,double variationNormal,double variationTangent,int interactiveTime,double meanShiftHposition,double meanShiftHnormal,float variationColor,float thresholdFunction,float thresholdGradient,float thresholdDistanceNormal,float radius,float timeStep,int tangentOrManifold,int ifNormalWeight,int ifAreaWeight,int ifVariationNormal);
	void DeleteNonShiftBilateralFilter();
protected:
	void ComputeMapKnearest();
	void ComputeNormal();
	void ComputeNormalColor();
	void FilterNormal();
	void FilterNormalColor();
	void NormPositionAndColor();
	void ComputeFilterOne(int interractiveTime);// 几何滤波  一致分布 
	void ComputeFilterTwo(int interractiveTime);//几何滤波   非一致分布
	void ComputeFilterThree(int interractiveTime);// 高维滤波
	void ComputeFilterFour(int interractiveTime);//高维滤波 保护volume
	void ComputeFilterFive(int interractiveTime);// 类似于函数1但保护volume
	void ComputeFilterSix(int interractiveTime);//类似于函数2但保护volume
	void ComputeFilterSeven(int interractiveTime);//类似于1 但加了一个核函数
	void ComputeFilterTen(int interractiveTime);//Pauly 方法， isotropic 滤波
	void ComputeFilterElevn(int interractiveTime);// Marx 方法， MLS
	void ComputeFilterTwentyone(int interractiveTime);//定义在流形上的双边滤波
	void ComputeFilterTwentytwo(int interractiveTime);//定义在切空间上的双边滤波
	void ComputeFilterTwentyThree(int interractiveTime);//加上了法向距离
	void ComputeFilterTwentyFour(int interractiveTime);//定义在切空间上的双边滤波  保护体积
	void ComputeFilterTwentyFive(int interractiveTime);//加上了法向投影项的面积不一致性质  保护体积
	void ComputeFilterTwentySix(int interractiveTime);//加上了法向距离 和 面积
	void ComputeFilterThirty(int interractiveTime);//加上了法向的距离,方差取为实际方差 同23
	
	void RestoreVolume();
	void RestorPositinData();
	void CalculateMeanRadius();
	void EstimateNormalDirection();
	




};

//struct lessFloat: public std::binary_function <float, float, bool> 
//{
//	bool operator()(const float& _Left, const float& _Right) const
//	{
//		if (std::fabs(_Left-_Right) < 1e-9)
//			return false;
//		else
//			return _Left < _Right;
//	}
//};
