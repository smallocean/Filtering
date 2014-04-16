#pragma once
#include "NormalStruct.h"

#if !defined (stiff_Matrix)
#define stiff_Matrix
	struct StiffMatrix{
	std::vector<double> m_elementStiffMarix;//��ŸնȾ����Ԫ��
};// �����նȾ���Ľṹ
#endif

class OneFreedomMeshless
{
public:
	OneFreedomMeshless(void);
	~OneFreedomMeshless(void);
public:
	static const int m_rankOfInterpolate;
	static const int m_numOfIntegralPoint;
	//static const double m_localInterpolatePoint[9][2];
	//static const double m_TestFunction[9];
	//static const double m_TestFunctionDerivativeOne[9];
	//static const double m_TestFunctionDerivativeTwo[9];
	//static const double m_gaussItegral[9];
	static const double m_ordinateAndWeightCoff[120];

protected:
	double* m_originalPointSet;
	double* m_originalColors;
	int m_numOfPoints;
	int m_kNearest;
	int m_IfFilterNormal;
	double m_meanShiftStopNormal;
	double* m_originalNormals; //ԭʼ����
	double* m_filterNormals; //�˲���ķ���
	double* m_originalMajorDirection;//ԭʼ�����򣨶�Ӧ����������ֵ��
	double* m_originalMinorDirection;//ԭʼ�η��򣨶�Ӧ�ڵڶ�������ֵ��
	double* m_meanCurvature;//ƽ������
	double* m_gaussCurvature; //��˹����
	double m_radius;//���� m_mapKnearestʱ��������뾶
	int m_interactiveTime; // �˲��Ĵ���
	double* m_loadVector;//��������
	StiffMatrix* m_stiffMatrix;// �նȾ���
	StiffMatrix* m_massMatrix;// ��������
	double m_namda;// ʱ�䲽��
	double m_loadConstant;//���س���
	double m_ebusainu;//ѭ��ֹͣ��������
protected://���ֵ�����йص�


protected:  //�������йصĲ���
	int m_numOfNeighbors;          //���������������ĵ�
	double m_localIntergralRadius;//      �ֲ����ְ뾶
	double m_localIntergralRadiusSquare;//�ֲ����ְ뾶��ƽ��
	double m_localWeightRadius;//         �ֲ���Ȩ�뾶
	double m_localWeightRadiusC;//�ֲ���Ȩ�뾶��4������
	double m_localWeightRadiusCSquare;//�ֲ���Ȩ�뾶��4������ƽ��
	double m_localIntergralRadiusC;//�ֲ����ְ뾶��4������
	double m_localIntergralRadiusCSquare;//�ֲ����ְ뾶��4������ƽ��
	double m_localWeightRadiusSquare;//   �ֲ���Ȩ�뾶��ƽ��
	double* m_localOrdinate;//             �������ľֲ�����
	double* m_localIntegralPoints;  //���ֵ������
	double* m_localQmatrix;//              ��ֵ��������
	double* m_localQTmatrix;//             ��ֵ������ת�þ���
protected:  //������йصĲ���
	double* m_kernelShapeFunction;
	double* m_kernelDerivativeOneShapeFunctionOne;
	double* m_kernelDerivativeTwoShapeFunctionOne;
	double* m_kernelDerivativeOneShapeFunctionTwo;
	double* m_kernelDerivativeTwoShapeFunctionTwo;
	double* m_testFunction;
	double* m_testFunctionDerivativeOne;
	double* m_testFunctionDerivativeTwo;
protected:  //�κ������䵼��
	double* m_shapeFunction;//                ���������κ����ڸ��������ֵ
	double* m_shapeFunctionDerivativeOne;//   ���������κ����Ա���һ�ĵ����ڸ��������ֵ
	double* m_shapeFunctionDerivativeTwo;//   ���������κ����Ա������ĵ����ڸ��������ֵ
	double* m_localOrdinateVectorOne;//�ռ�����Ծֲ�����ĵ�����
	double* m_gradientShapeFunction;  // �κ������ݶ�
	double* m_gradientTestFunction; //    ���Ժ������ݶ�
	double* m_localOrdinateVectorTwo;
protected:
	double* m_determinentMetric;//   �ڲ����㴦���α�׼�Ĵ�������ʽ
	double* m_metricMatrix;//�ڲ����㴦���α�׼�ľ���
	double* m_metricMatrixInve;////�ڲ����㴦���α�׼�ľ������

public:
	double* m_resultPointSet;
	double* m_resultColor;
	long nearestTime;
	long filterTime;
	std::map<int, KnearestField> m_mapKnearest;
public:
	void GetMeshlessFilter(int numOfPointSet,float* pointSet,int kNearest, double radius, int ifFilterNormal,int interactiveTime,double meanShiftHposition,double meanShiftHnormal,	double namda,	double loadConstant,	double ebusainu);
	/*
	�������ܣ�  �ӿں���������㼯�˲�
	����˵����
	int numOfPointSet                        �������
	double* pointSet							ԭʼ�㼯
	int kNearest                             k���������
	double radius                             k���������İ뾶
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
	*  double m_localIntergralRadius       �ֲ����ְ뾶
	*  double* m_localOrdinate             �������ľֲ�����
	*  double* m_localQmatrix              ��ֵ��������
	*  double* m_localQTmatrix             ��ֵ������ת�þ���
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
	*  double* m_kernelShapeFunction       
	*  double* m_kernelDerivativeOneShapeFunctionOne
	*  double* m_kernelDerivativeTwoShapeFunctionOne
	*  double* m_kernelDerivativeOneShapeFunctionTwo
	*  double* m_kernelDerivativeTwoShapeFunctionTwo
	*/

	void ComputeShapeFunction(int indexPointSet);
	/*
	*	�������ܣ������κ������䵼��
	*  ����˵����
	*  int indexPointSet        ������
	*  double* m_shapeFunction                ���������κ����ڸ��������ֵ
	*  double* m_shapeFunctionDerivativeOne   ���������κ����Ա���һ�ĵ����ڸ��������ֵ
	*  double* m_shapeFunctionDerivativeTwo   ���������κ����Ա������ĵ����ڸ��������ֵ
	*/

	void ComputeManifoldMetric(int indexPointSet);
	/*
	*	�������ܣ� �������α�׼�Ĵ�������ʽ
	*  ����˵����
	*  int indexPointSet        ������
	*  double* m_determinentMetric   �ڲ����㴦���α�׼�Ĵ�������ʽ
	*/

	void ComputeStiffMatrix(int indexPointSet);
	/*
	*	�������ܣ� ����նȾ���ͺ�������
	*  ����˵���� 
	*  int indexPointSet                    ������
	*  double* m_loadVector                  ��������
	*  StiffMatrix* m_stiffMatrix           �նȾ���
	*/

	void EstimateNormalDirection();
	/*
	*	�������ܣ� ���Ʒ����׼ȷ����
	*/
	void CalculateLinearSystem();
	/*
	*	�������ܣ� �����Է���
	*/

	void  ComputeGradientShapeFunction();
	/*
	 *	�������ܣ� ����ÿһ�����ֵ��κ������ݶȺͲ��Ժ������ݶ�
	 *  ����˵����
	 *  double* m_gradientShapeFunction   �κ������ݶ�
	 *  double* m_gradientTestFunction    ���Ժ������ݶ�
	 */

	void CalculateLinearSystemGauss();
	/*
	 *	�������ܣ� ��˹��Ԫ���ⷽ�̣�
	 */
	void ComputeTestFunction(int indexPointSet);
	/*
	 *	������Ժ����ڻ��ֵ�ĺ���ֵ���䵼��ֵ
	 */
	void AssembleStiffMatrix();
	//////////////////////////////////////////////////////////////////////////	
	/*
	*	�������ܣ� ��װ�նȾ���
	*/
	//////////////////////////////////////////////////////////////////////////

};

