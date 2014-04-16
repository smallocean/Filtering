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
};//���ڵ�ļ��������Ľṹ

class GeometricEstimation
{
public:
	//����ԭʼ����ʱ���󼸺������Ĺ��캯��
	GeometricEstimation(int numOfPointSet,float* pointSet,uint8 numOfKnearest,int* allPointsNearestNeighbor,float* allPointsNormal);
	//������ʱ�Ĺ��캯��
	GeometricEstimation(int numOfPointSet,float* pointSet);
	//��û��ԭʼ����ʱ���󼸺�����ʱ�Ĺ��캯��
	GeometricEstimation(int numOfPointSet,float* pointSet,UINT8 numOfKnearest,int* allPointsNearestNeighbor);
	~GeometricEstimation(void);
public:
	// ���ϸ����µĵ㼯 ���� ���� ��������
	float* m_allPointsNormal;//���е�ķ���
	int* m_allPointsNeighbor;//���е������
    float* m_pointSet;// �㼯
	float* m_allPointsPricipalCurvatureOne;
	float* m_allPointsPricipalCurvatureTwo;
	float* m_allPointsPricipalDirectionOne;
	float* m_allPointsPricipalDirectionTwo;
	float* m_allPointsGaussCurvature;
	float* m_allPointsAverageCurvature;
	PointGeometry* m_allPointGeometry;//�����Ƶ�ļ�������
protected:    
	//����������� ���� ��������
    POINT3D m_position;//�����Ƶ��λ��
    int m_ID;//�����Ƶ��ID
    int* m_pointNeighbor;//�����Ƶ������
    float m_averageCurvature;//�����Ƶ��ƽ������
    PointGeometry m_pointGeometry;//�����Ƶ�ļ�������
	float m_gaussCurvature;//�����Ƶ�ĸ�˹����
	float m_principalCurvature[2];//�����Ƶ��������
	VECTOR3D m_principalDirection[2];//�����Ƶ�����������
    VECTOR3D m_normal;//�����Ƶ�ķ���
    
	//��ͬ����
	float m_radius;//����뾶
	int m_numOfPoint;//�����Ŀ
	UINT8 m_numOfNeighbor;//ÿһ�����������Ŀ
	
	//uint8 m_knearest;//k���ڵ�k
	//float m_radius;//����k���ڵ���뾶
protected:
	//�������
	void NormalEstimation();
	//���ʹ���
	void CurvatureEstimation();
public:
	//��ʼ��
	void SetInput(UINT8 ID);
	//�ͷſռ�
	void Delete();
	//��������
	void SearchNeighbor(UINT8 knearest,float radius,int* pointNeighbor);
	//ȷ�����������İ뾶
	float EstimateSearchingRadius();
	//���ط���
	void GetNormal(VECTOR3D normal);
	//����ƽ������
	float GetAverageCurvature();
	//���ظ�˹����
	float GetGaussCurvature();
	//���ؼ��������������������
	void GetPointGeometricCharacter(PointGeometry pointGeometry);
	//��������rendering����Բ
	MZSplatEllipse GetSplatEllipse();
	//////////////////////////////////////////////////////////////////////////
	//������ȫ�ֻ� 
	void EstimateNormalDirection(int numOfPoints,float* pointSet,uint8 numOfEachPointNeighbor,int* allPointsNeighbor,float* allPointsNormal);
	//�õ����е�ķ���(û�г�ʼ����)
	void GetAllPointsNormalFirst();
	//�õ����е�ķ���(�г�ʼ����)
	void GetAllPointsNormalSecond();
	//�õ����е�ļ�������
	void GetAllPointsGeometricChracter();
	void GetAllPointsNeighbor(uint8 knearest,float radius);
	void GetAllPointsGaussAverageCurvature();

};
