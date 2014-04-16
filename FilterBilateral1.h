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
//	unsigned int m_numOfNearest;//���������������С��k
//	bool m_IfBoundary;//�Ƿ��Ǳ߽� �����Ϊ true������Ϊfalse
//	std::vector<int> m_nearest;//��������
//};// ����������Ľṹ
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
	//��������
	float* m_originalPointSet;
	float* m_originalColors;
	
public:
	//�������
	int m_numOfPoints;
	float* m_resultPointSet;
	float* m_resultColors;
	float* m_resultNormals;
	long nearestTime;
	long filterTime;
	float* m_originalNormals;
	std::map<int, KnearestField> m_mapKnearest;
protected:
	//��������
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
	//�ڼ���˫���˲�ʱ���õ�������
	std::vector<POINTVECTOR3D> m_localPositionEstimation;//�ֲ�λ�õ�һ�׹���
	std::vector<float> m_localPointDistance;//�����֮��ľ��룬����˫���˲��еĵ�һʽ��
	std::vector<float> m_localProjectionDistance;//�˲��㵽ͶӰ��֮��ľ���
	std::vector<float> m_localColorDistance[3];//��ɫ����
	std::vector<float> m_localColorProjectionDistance[3];//ͶӰ��ɫ�ľ���
	std::vector<Color> m_localColorEstimation;//�ֲ���ɫ��һ�ι���
	
	float m_pointVariation;//ԭʼ��ķ���
	float m_estimationVariation;//ͶӰ����ķ���
	float m_colorVariation[3];//��ɫ�ķ��� 
	float m_colorEstimationVariation[3];//��ɫͶӰ�ķ���
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
	void CalculateBilateralFilterSix(int interactiveTimes);//pauly û�еڶ��� ��ͶӰ�����Ȩ�  ֻ�м����˲�
	void CalculateBilateralFilterSeven(int interactiveTimes);// û��meanshift ������   ֻ�м����˲�
	void CalculateBilateralFilterEight(int interactiveTimes);//��meanshift ��  ֻ�м����˲�
	void CalculateBilateralFilterTen(int interactiveTimes);//���Ĺ��ƣ�
	void CalculateBilateralFilterEleven(int interactiveTimes);//���Ĺ��ƣ�
	void CalculateBilateralFilterTwelve(int interactiveTimes);

	void ComputeOneOrderPositionEstimation(int indexPoint);
	void ComputeOneOrderColorEstimation(int indexPoint);
	//����λ�ø�˹����
	void AddPositionGaussianNoise();
	//������ɫ��˹����
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
