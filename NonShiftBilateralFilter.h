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
	void ComputeFilterOne(int interractiveTime);// �����˲�  һ�·ֲ� 
	void ComputeFilterTwo(int interractiveTime);//�����˲�   ��һ�·ֲ�
	void ComputeFilterThree(int interractiveTime);// ��ά�˲�
	void ComputeFilterFour(int interractiveTime);//��ά�˲� ����volume
	void ComputeFilterFive(int interractiveTime);// �����ں���1������volume
	void ComputeFilterSix(int interractiveTime);//�����ں���2������volume
	void ComputeFilterSeven(int interractiveTime);//������1 ������һ���˺���
	void ComputeFilterTen(int interractiveTime);//Pauly ������ isotropic �˲�
	void ComputeFilterElevn(int interractiveTime);// Marx ������ MLS
	void ComputeFilterTwentyone(int interractiveTime);//�����������ϵ�˫���˲�
	void ComputeFilterTwentytwo(int interractiveTime);//�������пռ��ϵ�˫���˲�
	void ComputeFilterTwentyThree(int interractiveTime);//�����˷������
	void ComputeFilterTwentyFour(int interractiveTime);//�������пռ��ϵ�˫���˲�  �������
	void ComputeFilterTwentyFive(int interractiveTime);//�����˷���ͶӰ��������һ������  �������
	void ComputeFilterTwentySix(int interractiveTime);//�����˷������ �� ���
	void ComputeFilterThirty(int interractiveTime);//�����˷���ľ���,����ȡΪʵ�ʷ��� ͬ23
	
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
