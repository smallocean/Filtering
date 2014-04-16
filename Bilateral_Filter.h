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
struct POINTVECTOR6D {
	float pointVector[6];
};


class Bilateral_Filter  
{
public:
	Bilateral_Filter();
	virtual ~Bilateral_Filter();
protected:
	//��������
	float* m_originalPointSet;
	float* m_originalColors;
	float* m_originalNormals;
public:
	//�������
	int m_numOfPoints;
	float* m_resultPointSet;
	float* m_resultColors;
	float* m_resultNormals;
	long nearestTime;
	long filterTime;
protected:
	//��������
	std::map<int, KnearestField> m_mapKnearest;	
	float m_radius;
	int m_knearest;
	float m_meanShiftHposition, m_meanShiftHnormal;
	float m_meanShiftStopNormal;
	float* m_filterNormals;
	float* m_originalColorGradients[3];
	float* m_resultColorGradients[3];
	float m_maxPosition,m_minPosition;
protected:
	//�ڼ���˫���˲�ʱ���õ�������
	std::vector<POINTVECTOR6D> m_localPositionEstimation;//�ֲ�λ�õ�һ�׹���
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
	void GetFilterBilateral(int numOfPoints, float* pointSet, float* colors);
	void DeleteFilterBilateral();
protected:
	void NormalPosition();//��һ��λ������
	void CalculateVariation();
	void ComputeMapKnearest();
	void CalculateOriginalNormals();
	void MeanShiftNormals();
	void CalculateColorGradient();
	void CalculateBilateralFilterOne();
	void CalculateBilateralFilterTwo();
	void CalculateBilateralFilterThree();
	void ComputeOneOrderPositionEstimation(int indexPoint);
	void ComputeOneOrderColorEstimation(int indexPoint);
	void RecoverPosition();//�ָ�λ������
			
};

#endif // !defined(AFX_FILTERBILATERAL_H__0A901A27_1C99_49CA_ADE4_9D10C4B19C69__INCLUDED_)
