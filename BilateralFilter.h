/*********************************************************/
/*filename:    BilateralFilter.h                           */

/*description: This class is designed for point and mesh filtering */

/*author:      redstar    data: 2006/07/01                */
/*********************************************************/
#pragma once
#include <map>
#include <vector>
#include "Vectors.h"
//#include "mesh.h"
#include "NormalStruct.h"

class BilateralFilter
{
public:
	BilateralFilter(void);
	~BilateralFilter(void);
protected:
	//��������
	float* m_originalPoints;//ԭʼ��
	float* m_originalNormals;// ԭʼ��ķ���
	float* m_originalFilterNormals;//��ԭʼ��ķ����˲���ķ���
	float* m_filterColorGradient[3];
//	unsigned int* m_originalColor;//��ԭ�������ɫת��Ϊ0---1֮��
	float* m_originalColor;//ԭʼ�����ɫ 
	float* m_filterColor;//һ���˲������ɫ��Ԥ����
	float* m_colorGradient[3];//�������ɫ�ݶȺ�ԭ��Ĺ���ֵ
	float  m_meanShiftHposition,m_meanShiftHnormal,m_meanShiftHcolorGradient[3];//��meanshift �˲�ʱ��ĺ˿�ȣ�
	float m_meanShiftHcolor[3];
	float m_meanShiftStopNormal,m_meanShiftStopColorGradient[3];//mean Shift ֹͣ��׼


 //   float* m_originalNoiseColor;//����ɫ������//

protected:
	//�м�����
	int m_knearest;// k����Ĵ�С
    float m_radius;// k ����İ뾶��С
	std::map<int, KnearestField> m_mapKnearest;	
	//�˲�ʱ���õ�������
	 int m_numOFLocalPoints;//�ֲ��㼯������
	std::vector<int> m_LocalPoints;//�ֲ��㼯�ϣ����е�һ����ΪҪ�˲��ĵ�
	 int m_numOfLocalTriangles;//�ֲ����ǻ��������Ƭ��
	std::vector<POINTVECTOR3D> m_LocalTriangleCentroid;//�ֲ������ε�����
	std::vector<POINTVECTOR3D> m_originalCentroidNormals;////�ֲ������ε����ķ������������㷨���Ȩ	
	std::vector<POINTVECTOR3D> m_LocalPositionEstimation;//�ֲ�λ�õ�һ�׹���
	std::vector<float> m_localPointDistance;//�����֮��ľ��룬����˫���˲��еĵ�һʽ��
	std::vector<float> m_localProjectionDistance;//�˲��㵽ͶӰ��֮��ľ���
	std::vector<float> m_localColorDistance[3];//��ɫ����
	std::vector<Color> m_LocalColorEstimation;//�ֲ���ɫ��һ�ι���
	std::vector<Triangle> m_LocalTriangles;//�ֲ����ǻ� ��������������
	float m_pointVariation;//ԭʼ��ķ���
	float m_estimationVariation;//ͶӰ����ķ���
	float m_clolrVariation[3];//��ɫ�ķ���
	
public:
	//�������
	int m_numOfPoints;//�������
	int m_numOfTriangles;//����Ƭ������
	 int* m_triangles;// ����Ƭ����
	float* m_pointSets;//�˲���ĵ�
	float* m_normals;//�˲����ķ���
	float* m_colors;//������ɫ
public:
	void GetBilateralFilter(int numOfPointSets,float* pointSets, int numOfTriangle, int* triangles,float* vertexNormals,float* vertexColors);
    void DeleteBilateralFilter();
protected:
    //���ľֲ����򣬶���������˵������������2��ring����
	void ComputeMapKnearest();

	//���ľֲ����ǻ�
	void ComputeLocalTriangles();

	//��λ�ú���ɫ��һ��Ԥ��
	void ComputeOneOrderEstimation();

	//��λ�ú���ɫ��1.5��Ԥ��
	void ComputeOneHalfEstimation( int indexPoint);

	//��˫���˲������
	void ComputeBilateralFilter( int indexPoint);

	//�����˫���˲���molification)
	void ComputeNormalFilter( int indexPoint);
	
   //���˲��ĺ����ķ���ֵ �Լ��˲����õ��ľ����ƽ��
	void ComputeVariation(int indexPoint);

	//���������ε����
	double ComputeAreaOfTriangle( double position[3][3]);

	//����λ�ø�˹����
	void AddPositionGaussianNoise();

	//������ɫ��˹����
	void AddColorGaussianNorse();
    
	//��������ɫ����ͶӰ����ɫ
	void ComputeProjectionPointColor(Triangle localTriangle,float normals[3],POINTVECTOR3D centroid,POINTVECTOR3D estimatePoint,Color & colorEstimate);

	//�������˲����ķ���
	void ComputeNormals();

	//�����˲������ɫ��Ԥ����
	void ComputeFilterColor();

	//����һ�����������
	void CalculateCentroidOrdinate(double trianglePosition[3][3],double centroid[3],double anyPoint[3],double centroidOrdinate[3]);

	//������ɫ�仯����ȼ�ԭ�����ɫֵ
	void EstimateColorGradient(Triangle triangleIndex,POINTVECTOR3D centroidNormals,float gradientColor[4],float attributeValue[3]);

	// �˲���ɫ�仯����ȼ�ԭ�����ɫֵ
	void CalculateFilterColorGradient();


	// mean shift �˲�����
	void MeanShiftFilterNormal();
    // mean shift �˲���ɫ�ݶ�
	void MeanShiftFilterColorGradient();
};
