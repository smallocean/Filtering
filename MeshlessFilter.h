#pragma once
#include "NormalStruct.h"

#if !defined (stiff_Matrix)
#define stiff_Matrix
	struct StiffMatrix{
		std::vector<float> m_elementStiffMarix;//��ŸնȾ����Ԫ��
	};// �����նȾ���Ľṹ
#endif

class MeshlessFilter
{
public:
	MeshlessFilter(void);
	~MeshlessFilter(void);
public:
	static const int m_rankOfInterpolate;
	static const int m_numOfInterpolatePoint;
	static const float m_localInterpolatePoint[9][2];
	static const float m_TestFunction[9];
	static const float m_TestFunctionDerivativeOne[9];
	static const float m_TestFunctionDerivativeTwo[9];
	static const float m_gaussItegral[9];

protected:
	float* m_originalPointSet;
	float* m_originalColors;
	int m_numOfPoints;
	int m_kNearest;
	int m_IfFilterNormal;
	float m_meanShiftStopNormal;
	float* m_originalNormals; //ԭʼ����
	float* m_filterNormals; //�˲���ķ���
	float* m_originalMajorDirection;//ԭʼ�����򣨶�Ӧ����������ֵ��
	float* m_originalMinorDirection;//ԭʼ�η��򣨶�Ӧ�ڵڶ�������ֵ��
	float* m_meanCurvature;//ƽ������
	float* m_gaussCurvature; //��˹����
	float m_radius;//���� m_mapKnearestʱ��������뾶
	int m_interactiveTime; // �˲��Ĵ���
	float* m_loadVector;//��������
	StiffMatrix* m_stiffMatrix;// �նȾ���	
	float m_namda;// ʱ�䲽��
	float m_loadConstant;//���س���
	float m_ebusainu;//ѭ��ֹͣ��������
protected://���ֵ�����йص�
	

protected:  //�������йصĲ���
	int m_numOfNeighbors;          //���������������ĵ�
	float m_localIntergralRadius;//      �ֲ����ְ뾶
	float m_localIntergralRadiusSquare;//�ֲ����ְ뾶��ƽ��
	float m_localWeightRadius;//         �ֲ���Ȩ�뾶
	float m_localWeightRadiusSquare;//   �ֲ���Ȩ�뾶��ƽ��
	float* m_localOrdinate;//             �������ľֲ�����
	float* m_localQmatrix;//              ��ֵ��������
	float* m_localQTmatrix;//             ��ֵ������ת�þ���
protected:  //������йصĲ���
	float* m_kernelShapeFunction;
    float* m_kernelDerivativeOneShapeFunctionOne;
	float* m_kernelDerivativeTwoShapeFunctionOne;
	float* m_kernelDerivativeOneShapeFunctionTwo;
	float* m_kernelDerivativeTwoShapeFunctionTwo;
protected:  //�κ������䵼��
	float* m_shapeFunction;//                ���������κ����ڸ��������ֵ
	float* m_shapeFunctionDerivativeOne;//   ���������κ����Ա���һ�ĵ����ڸ��������ֵ
	float* m_shapeFunctionDerivativeTwo;//   ���������κ����Ա������ĵ����ڸ��������ֵ
protected:
	float* m_determinentMetric;//   �ڲ����㴦���α�׼�Ĵ�������ʽ
	float* m_metricMatrix;//�ڲ����㴦���α�׼�ľ���
	float* m_metricMatrixInve;////�ڲ����㴦���α�׼�ľ������

public:
	float* m_resultPointSet;
	float* m_resultColor;
	long nearestTime;
	long filterTime;
	std::map<int, KnearestField> m_mapKnearest;
public:
	void GetMeshlessFilter(int numOfPointSet,float* pointSet,int kNearest, float radius, int ifFilterNormal,int interactiveTime,float meanShiftHposition,float meanShiftHnormal,	float namda,	float loadConstant,	float ebusainu);
	/*
	�������ܣ�  �ӿں���������㼯�˲�
	����˵����
	int numOfPointSet                        �������
	float* pointSet							ԭʼ�㼯
	int kNearest                             k���������
	float radius                             k���������İ뾶
	int ifFilterNormal                       �Ƿ�Է�������˲�
	int interactiveTime                      �˲��Ľ�������
	double meanShiftHposition                �����˲�ʱ��λ�÷���
	double meanShiftHnormal                  �����˲�ʱ�ķ��򷽲�
	*/
	
	void DeleteMeshlessFilter();
	/*
	*	�������ܣ� �ͷű����ڴ�
	*/
protected:
	
	void ComputeMapKnearest();
	/*
	*	�������ܣ�����ÿһ����� k����
	*/
	
	void ComputeNormal();
	/*
	*	�������ܣ�����ÿһ����ķ���
	*/
	
	void ComputeNearestParameter(int indexPointSet);
	/*
	*	�������ܣ� ������һ�㴦�Ļ��ְ뾶���������������ľֲ����꣬ ����ֲ�Q������ת��QT
	*  ����˵����
	*  int indexPointSet        ������
	*  float m_localIntergralRadius       �ֲ����ְ뾶
	*  float* m_localOrdinate             �������ľֲ�����
	*  float* m_localQmatrix              ��ֵ��������
	*  float* m_localQTmatrix             ��ֵ������ת�þ���
	*/
	void ComputeQmatrix(int indexPointSet);
	/*
	 *	�������ܣ� �������Q����ת��
	 */
	
	void ComputeMatrixRelatedSampling(int indexPointSet);
	/*
	*	�������ܣ� ������������йصľ���
	*  ����˵���� 
	*  int indexPointSet        ������
	*  float* m_kernelShapeFunction       
	*  float* m_kernelDerivativeOneShapeFunctionOne
	*  float* m_kernelDerivativeTwoShapeFunctionOne
	*  float* m_kernelDerivativeOneShapeFunctionTwo
	*  float* m_kernelDerivativeTwoShapeFunctionTwo
	*/
	
	void ComputeShapeFunction(int indexPointSet);
	/*
	*	�������ܣ������κ������䵼��
	*  ����˵����
	*  int indexPointSet        ������
	*  float* m_shapeFunction                ���������κ����ڸ��������ֵ
	*  float* m_shapeFunctionDerivativeOne   ���������κ����Ա���һ�ĵ����ڸ��������ֵ
	*  float* m_shapeFunctionDerivativeTwo   ���������κ����Ա������ĵ����ڸ��������ֵ
	*/
	
	void ComputeManifoldMetric(int indexPointSet);
	/*
	*	�������ܣ� �������α�׼�Ĵ�������ʽ
	*  ����˵����
	*  int indexPointSet        ������
	*  float* m_determinentMetric   �ڲ����㴦���α�׼�Ĵ�������ʽ
	*/
	
	void ComputeStiffMatrix(int indexPointSet);
	/*
	*	�������ܣ� ����նȾ���ͺ�������
	*  ����˵���� 
	*  int indexPointSet                    ������
	*  float* m_loadVector                  ��������
	*  StiffMatrix* m_stiffMatrix           �նȾ���
	*/

	void MeshlessFilter:: EstimateNormalDirection();
	/*
	 *	�������ܣ� ���Ʒ����׼ȷ����
	 */
	void MeshlessFilter::CalculateLinearSystem();
	/*
	 *	�������ܣ� �����Է���
	 */


};
