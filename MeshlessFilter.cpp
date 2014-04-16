#include "stdafx.h"
#include ".\meshlessfilter.h"
#include ".\mathlib\mathlib.h"
#include "Mcube.h"
#include "Matrix.h"
using namespace MATHLIB;

const int MeshlessFilter::m_numOfInterpolatePoint=9;
const int MeshlessFilter::m_rankOfInterpolate=6;

const float MeshlessFilter::m_localInterpolatePoint[9][2]={
	{-0.7745966692,-0.7745966692},
	{0,-0.7745966692},
	{0.7745966692,-0.7745966692},
	{-0.7745966692,0},
	{0,0},
	{0.7745966692,0},
	{-0.7745966692,0.7745966692},
	{0,0.7745966692},
	{0.7745966692,0.7745966692}
};

const float MeshlessFilter::m_TestFunction[9]={
	0.00646952576012,
	0.09626480766995,
	0.00646952576012,
	0.09626480766995,
	1.43239451226107,
	0.09626480766995,
	0.00646952576012,
	0.09626480766995,
	0.00646952576012
};
const float MeshlessFilter::m_TestFunctionDerivativeOne[9]={
	0.04510145794586,
	0,
	-0.04510145794586,
	 0.67109759444089,
	 0,
	 -0.67109759444089,
	 0.04510145794586,
	 0,
	 -0.04510145794586
};
const float MeshlessFilter::m_TestFunctionDerivativeTwo[9]={
	0.04510145794586, 
	0.67109759444089,
	0.04510145794586,
	0,
	0,
	0,
	-0.04510145794586,
	-0.67109759444089,
	-0.04510145794586
};
const float MeshlessFilter::m_gaussItegral[9]={
    0.30864197535802469136,
	0.49382716053950617284,
	0.30864197535802469136,
	0.49382716053950617284,
	0.79012345680987654321,
	0.49382716053950617284,
	0.30864197535802469136,
	0.49382716053950617284,
	0.30864197535802469136
};

MeshlessFilter::MeshlessFilter(void)
{
    m_originalPointSet=NULL;
	m_originalColors=NULL;
	m_numOfPoints=0;
	m_kNearest=0;
	m_IfFilterNormal=0;
	m_meanShiftStopNormal=100;
	m_originalNormals=NULL; //原始法向
	m_filterNormals=NULL; //滤波后的法向
	m_originalMajorDirection=NULL;//原始主方向（对应于最大的特征值）
	m_originalMinorDirection=NULL;//原始次方向（对应于第二大特征值）
	m_gaussCurvature=NULL; //高斯曲率
	m_radius=0;//计算 m_mapKnearest时候的搜索半径
	m_interactiveTime=0; // 滤波的次数
	m_loadVector=NULL;//荷载向量
	m_stiffMatrix=NULL;// 刚度矩阵	
	m_localIntergralRadius=0;//      局部积分半径
	m_localOrdinate=NULL;//             各邻域点的局部坐标
	m_localQmatrix=NULL;//              插值函数矩阵
	m_localQTmatrix=NULL;//             插值函数的转置矩阵
	m_kernelShapeFunction=NULL;
	m_kernelDerivativeOneShapeFunctionOne=NULL;
	m_kernelDerivativeTwoShapeFunctionOne=NULL;
	m_kernelDerivativeOneShapeFunctionTwo=NULL;
	m_kernelDerivativeTwoShapeFunctionTwo=NULL;
	m_shapeFunction=NULL;//                各邻域点的形函数在各采样点的值
	m_shapeFunctionDerivativeOne=NULL;//   各邻域点的形函数对变量一的倒数在各采样点的值
	m_shapeFunctionDerivativeTwo=NULL;//   各邻域点的形函数对变量二的倒数在各采样点的值
	m_determinentMetric=NULL;//   在采样点处流形标准的代数行列式
	m_resultPointSet=NULL;
	m_resultColor=NULL;
}

MeshlessFilter::~MeshlessFilter(void)
{
	DeleteMeshlessFilter();
}

void  MeshlessFilter::GetMeshlessFilter(int numOfPointSet,float* pointSet,int kNearest, float radius, int ifFilterNormal,int interactiveTime,float meanShiftHposition,float meanShiftHnormal,float namda,	float loadConstant,	float ebusainu)
{
	/*
	函数功能：  接口函数，计算点集滤波
	变量说明：
	int numOfPointSet                        点的数量
	float* pointSet							 原始点集
	int kNearest                             k邻域的数量
	float radius                             k邻域搜索的半径
	int ifFilterNormal                       是否对法向进行滤波
	int interactiveTime                      滤波的交互次数
	double meanShiftHposition                法向滤波时的位置方差
	double meanShiftHnormal                  法向滤波时的法向方差
	float namda;                             时间步长
	float loadConstant;                      荷载常数
	float ebusainu;                          循环停止的条件；
	*/
	//////////////////////////////////////////////////////////////////////////
	// 变量初始化
	m_numOfPoints=numOfPointSet;
	m_kNearest=kNearest;
	m_radius=radius;
	//m_IfFilterNormal=ifFilterNormal;
	m_interactiveTime=interactiveTime;
	m_namda=namda;
	m_loadConstant=loadConstant;
	m_ebusainu=ebusainu;
	m_originalPointSet=new float[m_numOfPoints*3];
	m_resultPointSet=new float[m_numOfPoints*3];
	m_stiffMatrix=new StiffMatrix[m_numOfPoints];
	m_loadVector=new float[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints*3;i++){
		m_originalPointSet[i]=pointSet[i];
		m_resultPointSet[i]=pointSet[i];
	}
	
	//////////////////////////////////////////////////////////////////////////
	// 求k 近邻
	ComputeMapKnearest();
	//////////////////////////////////////////////////////////////////////////
	// 求法向 和两个主方向
	ComputeNormal();
	//////////////////////////////////////////////////////////////////////////
	//循环求解
	for(int i=0;i<m_interactiveTime;i++){
		//////////////////////////////////////////////////////////////////////////
		// 建立刚度矩阵
		for(int j=0;j<m_numOfPoints;j++){
			if(j==6335){
				int m;
				m=0;
			}
			ComputeNearestParameter(j);
			ComputeQmatrix(j);
			ComputeMatrixRelatedSampling(j);
			ComputeShapeFunction(j);
			ComputeManifoldMetric(j);
			ComputeStiffMatrix(j);
		}
		//////////////////////////////////////////////////////////////////////////
		// 解线性系统
		CalculateLinearSystem();
		//////////////////////////////////////////////////////////////////////////
		// 释放内存
		if(m_localOrdinate!=NULL)
            delete[] m_localOrdinate;
		m_localOrdinate=NULL;
		if(m_localQmatrix!=NULL)
			delete[] m_localQmatrix;
		m_localQmatrix=NULL;
		if(m_localQTmatrix!=NULL)
			delete[] m_localQTmatrix;
		m_localQTmatrix=NULL;
		if(m_kernelShapeFunction!=NULL)
			delete[] m_kernelShapeFunction;
		m_kernelShapeFunction=NULL;
		if(m_kernelDerivativeOneShapeFunctionOne!=NULL)
			delete[] m_kernelDerivativeOneShapeFunctionOne;
		m_kernelDerivativeOneShapeFunctionOne=NULL;
		if(m_kernelDerivativeTwoShapeFunctionOne!=NULL)
			delete[] m_kernelDerivativeTwoShapeFunctionOne;
		m_kernelDerivativeTwoShapeFunctionOne=NULL;
		if(m_kernelDerivativeOneShapeFunctionTwo!=NULL)
			delete[] m_kernelDerivativeOneShapeFunctionTwo;
		m_kernelDerivativeOneShapeFunctionTwo=NULL;
		if(m_kernelDerivativeTwoShapeFunctionTwo!=NULL)
			delete[] m_kernelDerivativeTwoShapeFunctionTwo;
		m_kernelDerivativeTwoShapeFunctionTwo=NULL;
		if(m_shapeFunction!=NULL)
			delete[] m_shapeFunction;
		m_shapeFunction=NULL; //一行为一个邻域点的形函数在不同的积分点的值
		if(m_shapeFunctionDerivativeOne!=NULL)
			delete[] m_shapeFunctionDerivativeOne;
		m_shapeFunctionDerivativeOne=NULL;
		if(m_shapeFunctionDerivativeTwo!=NULL)
			delete[] m_shapeFunctionDerivativeTwo;
		m_shapeFunctionDerivativeTwo=NULL;
		if(m_determinentMetric!=NULL)
			delete[] m_determinentMetric;
		m_determinentMetric=NULL;
		for(int j=0;j<m_numOfPoints;j++){
			m_stiffMatrix[j].m_elementStiffMarix.clear();
		}

	}
	if(m_originalPointSet!=NULL)
		delete[] m_originalPointSet;
	m_originalPointSet=NULL;
	if(m_stiffMatrix!=NULL)
		delete[] m_stiffMatrix;
	m_stiffMatrix=NULL;
	if(m_loadVector!=NULL)
		delete[] m_loadVector;
	m_loadVector=NULL;

	return;
}
void MeshlessFilter:: DeleteMeshlessFilter()
{
	if(m_originalPointSet!=NULL)
		delete[] m_originalPointSet;
	m_originalPointSet=NULL;
	if(m_originalColors!=NULL)
		delete[] m_originalColors;
    m_originalColors=NULL;
	if(m_originalNormals!=NULL)
		delete[] m_originalNormals;
	m_originalNormals=NULL; //原始法向
	if(m_filterNormals!=NULL)
		delete[] m_filterNormals;
	m_filterNormals=NULL; //滤波后的法向
	if(m_originalMajorDirection!=NULL)
		delete[] m_originalMajorDirection;
	m_originalMajorDirection=NULL;//原始主方向（对应于最大的特征值）
	if(m_originalMinorDirection!=NULL)
		delete[] m_originalMinorDirection;
	m_originalMinorDirection=NULL;//原始次方向（对应于第二大特征值)
	if(m_gaussCurvature!=NULL)
		delete[] m_gaussCurvature;
	m_gaussCurvature=NULL; //高斯曲率
    if(m_loadVector!=NULL)
		delete[] m_loadVector;
	m_loadVector=NULL;//荷载向量
	if(m_stiffMatrix!=NULL)
		delete[] m_stiffMatrix;
	m_stiffMatrix=NULL;// 刚度矩阵	
	if(m_localOrdinate!=NULL)
		delete[] m_localOrdinate;
	m_localOrdinate=NULL;//             各邻域点的局部坐标
	if(m_localQmatrix!=NULL)
		delete[] m_localQmatrix;
	m_localQmatrix=NULL;//              插值函数矩阵
	if(m_localQTmatrix!=NULL)
		delete[] m_localQTmatrix;
	m_localQTmatrix=NULL;//             插值函数的转置矩阵
	if(m_kernelShapeFunction!=NULL)
		delete[] m_kernelShapeFunction;
	m_kernelShapeFunction=NULL;
	if(m_kernelDerivativeOneShapeFunctionOne!=NULL)
		delete[] m_kernelDerivativeOneShapeFunctionOne;
	m_kernelDerivativeOneShapeFunctionOne=NULL;
	if(m_kernelDerivativeTwoShapeFunctionOne!=NULL)
		delete[] m_kernelDerivativeTwoShapeFunctionOne;
	m_kernelDerivativeTwoShapeFunctionOne=NULL;
	if(m_kernelDerivativeOneShapeFunctionTwo!=NULL)
		delete[] m_kernelDerivativeOneShapeFunctionTwo;
	m_kernelDerivativeOneShapeFunctionTwo=NULL;
	if(m_kernelDerivativeTwoShapeFunctionTwo!=NULL)
	m_kernelDerivativeTwoShapeFunctionTwo=NULL;
	if(m_shapeFunction!=NULL)
		delete[] m_shapeFunction;
	m_shapeFunction=NULL;//                各邻域点的形函数在各采样点的值
	if(m_shapeFunctionDerivativeOne!=NULL)
		delete[] m_shapeFunctionDerivativeOne;
	m_shapeFunctionDerivativeOne=NULL;//   各邻域点的形函数对变量一的倒数在各采样点的值
	if(m_shapeFunctionDerivativeTwo!=NULL)
		delete[] m_shapeFunctionDerivativeTwo;
	m_shapeFunctionDerivativeTwo=NULL;//   各邻域点的形函数对变量二的倒数在各采样点的值
    if(m_determinentMetric!=NULL)
		delete[] m_determinentMetric;
	m_determinentMetric=NULL;//   在采样点处流形标准的代数行列式
	if(m_resultPointSet!=NULL)
		delete[] m_resultPointSet;
	m_resultPointSet=NULL;
	if(m_resultColor!=NULL)
		delete[] m_resultColor;
	m_resultColor=NULL;
}
void MeshlessFilter::ComputeMapKnearest()
{
	float radius=m_radius*m_radius;//存放近邻点到所求点的距离
	std::multimap<float ,int> mapKnearest;//对每一个点建立邻域vector时建立临时map
	int numOfKnearest;//统计一个点的临时邻域大小

	//////////////////////////////////////////////////////////////////////////
	//遍历所有点，并且按照距离由小到大存放存放到map中
	for(int i=0;i<m_numOfPoints;i++){
		numOfKnearest=0;
		for(int j=0;j<m_numOfPoints;j++){
			if(j==i)
				continue;
			float distancePointToPoint;
			distancePointToPoint=(m_originalPointSet[i*3]-m_originalPointSet[j*3])*(m_originalPointSet[i*3]-m_originalPointSet[j*3])
				+(m_originalPointSet[i*3+1]-m_originalPointSet[j*3+1])*(m_originalPointSet[i*3+1]-m_originalPointSet[j*3+1])
				+(m_originalPointSet[i*3+2]-m_originalPointSet[j*3+2])*(m_originalPointSet[i*3+2]-m_originalPointSet[j*3+2]);
			if(distancePointToPoint>radius)
				continue;
			mapKnearest.insert(std::multimap<float,int>::value_type(distancePointToPoint,j));
			numOfKnearest+=1;
			if(numOfKnearest>m_kNearest){
				std::multimap<float,int>::iterator mapIterator = mapKnearest.end();
				mapIterator--;
				float w;
				w=(*mapIterator).first;
				mapKnearest.erase(mapIterator);
				numOfKnearest-=1;
			}
		}
		KnearestField fieldKnearest;//临时结构
		fieldKnearest.m_numOfNearest=numOfKnearest;
		fieldKnearest.m_IfBoundary=FALSE;
		std::multimap<float,int>::iterator mapIterator = mapKnearest.begin();
		for(;mapIterator!=mapKnearest.end();mapIterator++){
			fieldKnearest.m_nearest.push_back((*mapIterator).second);
		}
		m_mapKnearest.insert(std::map<int, KnearestField>::value_type(i,fieldKnearest));

		mapKnearest.clear();
	}	
	//////////////////////////////////////////////////////////////////////////
	//// 测试代码
	CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\knearest.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(i);
			std::vector<int>::iterator vectorIterator=(*mapKnearestIterator).second.m_nearest.begin();
			fprintf(fpout,"%d ",i);
			for(;vectorIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorIterator++){
				fprintf(fpout,"%d ",(*vectorIterator));
			}
			fprintf(fpout,"\n");
		}
		fclose(fpout);
	}
}
void MeshlessFilter::ComputeNormal()
{
	float centroidPosition[3];//重心位置
	float* localVariationVector;//重心到点的位置的向量
	int numOfNearest;

	mat_f8 mat_covaMatrix(3, 3);
	if(m_originalNormals!=NULL)
		delete[] m_originalNormals;
	m_originalNormals=new float[m_numOfPoints*3];
	m_originalMajorDirection=new float[m_numOfPoints*3];
	m_originalMinorDirection=new float[m_numOfPoints*3];
	m_gaussCurvature=new float[m_numOfPoints*2];
	//if(m_meanAreas!=NULL)
	//	delete[] m_meanAreas;
	//m_meanAreas=new float[m_numOfPoints];
	//if(m_meanCurvature!=NULL)
	//	delete[] m_meanCurvature;
	//m_meanCurvature=new float[m_numOfPoints];
	//if(m_resultColor!=NULL)
	//	delete[] m_resultColor;
	//m_resultColor=new float[m_numOfPoints*3];


	for(int i=0;i<m_numOfPoints;i++){
		std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		numOfNearest=(*mapIterator).second.m_numOfNearest;
		for(int j=0;j<3;j++){
			centroidPosition[j]=m_originalPointSet[i*3+j];
		}
		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
			centroidPosition[0]+=m_originalPointSet[(*vectorIterator)*3];
			centroidPosition[1]+=m_originalPointSet[(*vectorIterator)*3+1];
			centroidPosition[2]+=m_originalPointSet[(*vectorIterator)*3+2];
		}
		for(int j=0;j<3;j++){
			centroidPosition[j]/=(numOfNearest+1);
		}
		localVariationVector=new float[(numOfNearest+1)*3];
		localVariationVector[0]=m_originalPointSet[i*3]-centroidPosition[0];
		localVariationVector[1]=m_originalPointSet[i*3+1]-centroidPosition[1];
		localVariationVector[2]=m_originalPointSet[i*3+2]-centroidPosition[2];
		vectorIterator=(*mapIterator).second.m_nearest.begin();
		for(int j=1;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,j++){
			localVariationVector[j*3]=m_originalPointSet[(*vectorIterator)*3]-centroidPosition[0];
			localVariationVector[j*3+1]=m_originalPointSet[(*vectorIterator)*3+1]-centroidPosition[1];
			localVariationVector[j*3+2]=m_originalPointSet[(*vectorIterator)*3+2]-centroidPosition[2];
		}
		//求方差矩阵
		for(int j=0;j<3;j++){
			for(int k=0;k<3;k++){
				mat_covaMatrix(j,k)=0;

				for(int m=0;m<(numOfNearest+1);m++){
					mat_covaMatrix(j,k)=mat_covaMatrix(j,k)+localVariationVector[m*3+j]*localVariationVector[m*3+k];
				}
			}
		}
		vec_f8 eval(3);
		mat_f8 evec(3, 3);
		eigen_symm(eval, evec, mat_covaMatrix);
		m_originalNormals[i*3]=evec(0,2);
		m_originalNormals[i*3+1]=evec(1,2);
		m_originalNormals[i*3+2]=evec(2,2);
		m_originalMajorDirection[i*3]=evec(0,0);
		m_originalMajorDirection[i*3+1]=evec(1,0);
		m_originalMajorDirection[i*3+2]=evec(2,0);
		m_originalMinorDirection[i*3]=evec(0,1);
		m_originalMinorDirection[i*3+1]=evec(1,1);
		m_originalMinorDirection[i*3+2]=evec(2,1);
		m_gaussCurvature[i*2]=eval(0);
		m_gaussCurvature[i*2+1]=eval(1);
		//m_meanAreas[i]=abs(eval(0)*eval(2));
	/*	m_meanCurvature[i]=eval(2)/(eval(0)+eval(1)+eval(2));
		if(m_meanCurvature[i]*4<=0.005){
			m_resultColor[i*3]=0;
			m_resultColor[i*3+1]=0;
			m_resultColor[i*3+2]=0.5+80*m_meanCurvature[i]*4;
		}
		if((m_meanCurvature[i]*4>0.005)&(m_meanCurvature[i]*4<=0.01)){
			m_resultColor[i*3]=0;
			m_resultColor[i*3+1]=0.45+50*(m_meanCurvature[i]*4-0.005);
			m_resultColor[i*3+2]=0;
		}
		if(m_meanCurvature[i]*4>0.01){
			m_resultColor[i*3]=0.15+5*(m_meanCurvature[i]*4-0.01);
			m_resultColor[i*3+1]=0;
			m_resultColor[i*3+2]=0;

		}*/

		delete[] localVariationVector;
	}
	EstimateNormalDirection();

	//////////////////////////////////////////////////////////////////////////
	// 写法向文件
	CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\normal.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f\n",m_originalNormals[i*3],m_originalNormals[i*3+1],m_originalNormals[i*3+2]);

		}
		fclose(fpout);
	}

	//////////////////////////////////////////////////////////////////////////
	// 写主方向文件
	filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\majorDirection.txt";
	
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f\n",m_originalMajorDirection[i*3],m_originalMajorDirection[i*3+1],m_originalMajorDirection[i*3+2]);

		}
		fclose(fpout);
	}

	//////////////////////////////////////////////////////////////////////////
	// 写次方向文件
	filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\minorDirection.txt";

	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f\n",m_originalMinorDirection[i*3],m_originalMinorDirection[i*3+1],m_originalMinorDirection[i*3+2]);

		}
		fclose(fpout);
	}
	

	//filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonShiftmeanCurvature.txt";

	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f\n",m_meanCurvature[i]);
	//	}
	//	fclose(fpout);
	//}
}
void MeshlessFilter::ComputeNearestParameter(int indexPointSet)
{
	//int numOfNearest;
	float localPoint[3], localNormal[3],localMajorDirection[3],localMinorDirection[3];
    std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(indexPointSet);
	std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	m_numOfNeighbors=(*mapIterator).second.m_numOfNearest;
	if(m_localOrdinate!=NULL)
		delete[] m_localOrdinate;
	m_localOrdinate=NULL;
	m_localOrdinate=new float[m_numOfNeighbors*2];
	m_localIntergralRadius=0;
	for(int i=0;i<3;i++){
		localMajorDirection[i]=m_originalMajorDirection[indexPointSet*3+i];
		localMinorDirection[i]=m_originalMinorDirection[indexPointSet*3+i];
	}
	
	for(int i=0;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,i++){
        float vectorNeighborToPoint[3];
		float localRadius;
		for(int j=0;j<3;j++){
			vectorNeighborToPoint[j]=m_resultPointSet[(*vectorIterator)*3+j]-m_resultPointSet[indexPointSet*3+j];
		}
		m_localOrdinate[i*2]=Vector3Vector(vectorNeighborToPoint,localMajorDirection);
		m_localOrdinate[i*2+1]=Vector3Vector(vectorNeighborToPoint,localMinorDirection);    
	
	}
	m_localIntergralRadiusSquare=m_localOrdinate[0]*m_localOrdinate[0]+m_localOrdinate[1]*m_localOrdinate[1];
	m_localIntergralRadius=sqrt(m_localIntergralRadiusSquare);

	m_localWeightRadiusSquare=m_localOrdinate[(m_numOfNeighbors-1)*2]*m_localOrdinate[(m_numOfNeighbors-1)*2]
	         +m_localOrdinate[(m_numOfNeighbors-1)*2+1]*m_localOrdinate[(m_numOfNeighbors-1)*2+1];
    m_localWeightRadius=sqrt(m_localWeightRadiusSquare);
	if(m_localWeightRadius<=m_localIntergralRadius){
		float tempChange;
		tempChange=m_localWeightRadius;
		m_localWeightRadius=m_localIntergralRadius;
		m_localIntergralRadius=m_localWeightRadius;
		tempChange=m_localWeightRadiusSquare;
		m_localWeightRadiusSquare=m_localIntergralRadiusSquare;
		m_localIntergralRadiusSquare=m_localWeightRadiusSquare;
	}


	//////////////////////////////////////////////////////////////////////////
	// 写局部坐标文件
	CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\localOrdinalte.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "a")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		fprintf(fpout,"%d\n",indexPointSet);
		for(int i=0;i<m_numOfNeighbors;i++){
			fprintf(fpout,"%f %f\n",m_localOrdinate[i*2],m_localOrdinate[i*2+1]);

		}
		fclose(fpout);
	}

	
	return;
	

}
void MeshlessFilter::ComputeQmatrix(int indexPointSet)
{
	m_localQmatrix=new float[m_rankOfInterpolate*(m_numOfNeighbors+1)];//按行排列
	m_localQTmatrix=new float[m_rankOfInterpolate*(m_numOfNeighbors+1)];//按列排列
	m_localQmatrix[0]=1;
	for(int i=1;i<m_rankOfInterpolate;i++){
		m_localQmatrix[i]=0;
	}
	for(int i=1;i<m_numOfNeighbors+1;i++){
        m_localQmatrix[i*m_rankOfInterpolate]=1;
		m_localQmatrix[i*m_rankOfInterpolate+1]=m_localOrdinate[(i-1)*2];
		m_localQmatrix[i*m_rankOfInterpolate+2]=m_localOrdinate[(i-1)*2+1];
		m_localQmatrix[i*m_rankOfInterpolate+3]=m_localOrdinate[(i-1)*2]*m_localOrdinate[(i-1)*2];
		m_localQmatrix[i*m_rankOfInterpolate+4]=m_localOrdinate[(i-1)*2]*m_localOrdinate[(i-1)*2+1];
		m_localQmatrix[i*m_rankOfInterpolate+5]=m_localOrdinate[(i-1)*2+1]*m_localOrdinate[(i-1)*2+1];

	}
	for(int i=0;i<m_numOfNeighbors+1;i++){
		for(int j=0;j<m_rankOfInterpolate;j++){
			m_localQTmatrix[j*(m_numOfNeighbors+1)+i]=m_localQmatrix[i*m_rankOfInterpolate+j];
		}
	}
	//////////////////////////////////////////////////////////////////////////
	// 写文件
	CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\localQmatrix.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "a")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		fprintf(fpout,"%d\n",indexPointSet);
		for(int i=0;i<m_numOfNeighbors+1;i++){
			for(int j=0;j<m_rankOfInterpolate;j++){
				fprintf(fpout,"%f ",m_localQmatrix[i*m_rankOfInterpolate+j]);
			}
			fprintf(fpout," \n");            
		}
		fclose(fpout);
	}
	return;    
}
void MeshlessFilter::ComputeMatrixRelatedSampling(int indexPointSet)
{
	/*
	*	函数功能： 计算与采样点有关的矩阵
	*  变量说明： 
	*  int indexPointSet        索引点
	*  float* m_kernelShapeFunction       
	*  float* m_kernelDerivativeOneShapeFunctionOne
	*  float* m_kernelDerivativeTwoShapeFunctionOne
	*  float* m_kernelDerivativeOneShapeFunctionTwo
	*  float* m_kernelDerivativeTwoShapeFunctionTwo
	*/
	//////////////////////////////////////////////////////////////////////////
	// 分配内存空间
	if(m_kernelShapeFunction!=NULL)
		delete[] m_kernelShapeFunction;
	m_kernelShapeFunction=new float[m_rankOfInterpolate*m_numOfInterpolatePoint*(m_numOfNeighbors+1)];
	if(m_kernelDerivativeOneShapeFunctionOne!=NULL)
		delete[] m_kernelDerivativeOneShapeFunctionOne;
	m_kernelDerivativeOneShapeFunctionOne=new float[m_rankOfInterpolate*m_numOfInterpolatePoint*(m_numOfNeighbors+1)];
	if(m_kernelDerivativeTwoShapeFunctionOne!=NULL)
		delete[] m_kernelDerivativeTwoShapeFunctionOne;
	m_kernelDerivativeTwoShapeFunctionOne=new float[m_rankOfInterpolate*m_numOfInterpolatePoint*(m_numOfNeighbors+1)];
	if(m_kernelDerivativeOneShapeFunctionTwo!=NULL)
		delete[] m_kernelDerivativeOneShapeFunctionTwo;
	m_kernelDerivativeOneShapeFunctionTwo=new float[m_rankOfInterpolate*m_numOfInterpolatePoint*(m_numOfNeighbors+1)];
	if(m_kernelDerivativeTwoShapeFunctionTwo!=NULL)
		delete[] m_kernelDerivativeTwoShapeFunctionTwo;
	m_kernelDerivativeTwoShapeFunctionTwo=new float[m_rankOfInterpolate*m_numOfInterpolatePoint*(m_numOfNeighbors+1)];
	//////////////////////////////////////////////////////////////////////////
	// 开始循环计算
	for(int i=0;i<m_numOfInterpolatePoint;i++){
		//////////////////////////////////////////////////////////////////////////
		// 定义临时变量
		float* weightVector=new float[m_numOfNeighbors+1];  // 加权矩阵
		float* weightVectorDerivativeOne=new float[m_numOfNeighbors+1];  //加权矩阵对第一个变量的导数
		float* weightVectorDerivativeTwo=new float[m_numOfNeighbors+1];  //加权矩阵对第二个变量的导数
		float* QTweight=new float[m_rankOfInterpolate*(m_numOfNeighbors+1)];
		float* QTweightDerivativeOne=new float[m_rankOfInterpolate*(m_numOfNeighbors+1)];
		float* QTweightDerivativeTwo=new float[m_rankOfInterpolate*(m_numOfNeighbors+1)];
		float* QTweightQ=new float[m_rankOfInterpolate*m_rankOfInterpolate];
		float* QTweightQDerivativeOne=new float[m_rankOfInterpolate*m_rankOfInterpolate];
		float* QTweightQDerivativeTwo=new float[m_rankOfInterpolate*m_rankOfInterpolate];
		float* QTweightQInver=new float[m_rankOfInterpolate*m_rankOfInterpolate];

		//////////////////////////////////////////////////////////////////////////
		// 计算加权矩阵及其导数阵
		float localInterpolatePoint[2];
		localInterpolatePoint[0]=m_localInterpolatePoint[i][0]*m_localIntergralRadius;
		localInterpolatePoint[1]=m_localInterpolatePoint[i][1]*m_localIntergralRadius;

		weightVector[0]=exp(-0.5*(localInterpolatePoint[0]*localInterpolatePoint[0]
		                        +localInterpolatePoint[1]*localInterpolatePoint[1])/m_localWeightRadiusSquare);
		weightVectorDerivativeOne[0]=weightVector[0]*(-1)*localInterpolatePoint[0]/m_localWeightRadiusSquare;
		weightVectorDerivativeTwo[0]=weightVector[0]*(-1)*localInterpolatePoint[1]/m_localWeightRadiusSquare;

		for(int j=0;j<m_numOfNeighbors;j++){
			float localOrdinateMinus[2];
			localOrdinateMinus[0]=m_localOrdinate[j*2]-localInterpolatePoint[0];
			localOrdinateMinus[1]=m_localOrdinate[j*2+1]-localInterpolatePoint[1];
			weightVector[j+1]=exp(-0.5*(localOrdinateMinus[0]*localOrdinateMinus[0]
							+localOrdinateMinus[1]*localOrdinateMinus[1])/m_localWeightRadiusSquare);
			weightVectorDerivativeOne[j+1]=weightVector[j+1]*localOrdinateMinus[0]/m_localWeightRadiusSquare;
			weightVectorDerivativeTwo[j+1]=weightVector[j+1]*localOrdinateMinus[1]/m_localWeightRadiusSquare;
		}
		//////////////////////////////////////////////////////////////////////////
		//  计算矩阵 QTweight  QTweightDerivativeOne  QTweightDerivativeTwo
		for(int j=0;j<m_rankOfInterpolate;j++){
			for(int k=0;k<m_numOfNeighbors+1;k++){
				QTweight[j*(m_numOfNeighbors+1)+k]=m_localQTmatrix[j*(m_numOfNeighbors+1)+k]*weightVector[k];
				QTweightDerivativeOne[j*(m_numOfNeighbors+1)+k]=m_localQTmatrix[j*(m_numOfNeighbors+1)+k]*weightVectorDerivativeOne[k];
				QTweightDerivativeTwo[j*(m_numOfNeighbors+1)+k]=m_localQTmatrix[j*(m_numOfNeighbors+1)+k]*weightVectorDerivativeTwo[k];

			}
		}
		//////////////////////////////////////////////////////////////////////////
		// 用 mat_f8 定义 矩阵 并赋值
		mat_f8 mat_localQmatrix((m_numOfNeighbors+1),m_rankOfInterpolate);
		mat_f8 mat_localQTweight(m_rankOfInterpolate,(m_numOfNeighbors+1));
		mat_f8 mat_localQTweightDeOne(m_rankOfInterpolate,(m_numOfNeighbors+1));
		mat_f8 mat_localQTweightDeTwo(m_rankOfInterpolate,(m_numOfNeighbors+1));

		for(int j=0;j<m_numOfNeighbors+1;j++){
			for(int k=0;k<m_rankOfInterpolate;k++){
				mat_localQmatrix(j,k)=m_localQmatrix[j*m_rankOfInterpolate+k];
				mat_localQTweight(k,j)=QTweight[k*(m_numOfNeighbors+1)+j];
				mat_localQTweightDeOne(k,j)=QTweightDerivativeOne[k*(m_numOfNeighbors+1)+j];
				mat_localQTweightDeTwo(k,j)=QTweightDerivativeTwo[k*(m_numOfNeighbors+1)+j];
			}
		}

		//////////////////////////////////////////////////////////////////////////
		// 计算矩阵  QTweightQ  QTweightQDerivativeOne  QTweightQDerivativeTwo
		mat_f8 mat_localQTweightQ(m_rankOfInterpolate,m_rankOfInterpolate);
		mat_f8 mat_localQTweightDeOneQ(m_rankOfInterpolate,m_rankOfInterpolate);
		mat_f8 mat_localQTweightDeTwoQ(m_rankOfInterpolate,m_rankOfInterpolate);
		mat_localQTweightQ=mat_localQTweight*mat_localQmatrix;
		mat_localQTweightDeOneQ=mat_localQTweightDeOne*mat_localQmatrix;
		mat_localQTweightDeTwoQ=mat_localQTweightDeTwo*mat_localQmatrix;
		//////////////////////////////////////////////////////////////////////////
		// 计算矩阵  QTweightQDInver
		mat_f8 mat_QTweightQInver(m_rankOfInterpolate,m_rankOfInterpolate);
		mat_QTweightQInver=mat_localQTweightQ;
		if(!matrix_inverse(mat_QTweightQInver)){
			//MessageBox("求逆矩阵出错");            
		}
		//////////////////////////////////////////////////////////////////////////
		// 计算形函数的核矩阵
		mat_f8 mat_shapeMatrixOne(m_rankOfInterpolate,m_numOfNeighbors+1);
		mat_f8 mat_shapeMatrixTwo(m_rankOfInterpolate,m_numOfNeighbors+1);
		mat_f8 mat_shapeMatrixThree(m_rankOfInterpolate,m_numOfNeighbors+1);
		mat_f8 mat_shapeMatrixFour(m_rankOfInterpolate,m_numOfNeighbors+1);
		mat_f8 mat_shapeMatrixFive(m_rankOfInterpolate,m_numOfNeighbors+1);
		mat_shapeMatrixOne=mat_QTweightQInver*mat_localQTweight;
		mat_shapeMatrixTwo=mat_QTweightQInver*mat_localQTweightDeOne;
		mat_shapeMatrixThree=mat_QTweightQInver*mat_localQTweightDeTwo;
		mat_shapeMatrixFour=mat_QTweightQInver*mat_localQTweightDeOneQ;
		mat_shapeMatrixFour=mat_shapeMatrixFour*mat_shapeMatrixOne;
		
		mat_shapeMatrixFive=mat_QTweightQInver*mat_localQTweightDeTwoQ;
		mat_shapeMatrixFive=mat_shapeMatrixFive*mat_shapeMatrixOne;
		//mat_shapeMatrixFive=mat_shapeMatrixFive*mat_localQTweight;
		//////////////////////////////////////////////////////////////////////////
		// 赋值给形函数矩阵
		int sizeNeighbor=m_numOfNeighbors+1;
		int sizeOfMatrix=sizeNeighbor*m_rankOfInterpolate;
		for(int j=0;j<m_rankOfInterpolate;j++){
			for(int k=0;k<sizeNeighbor;k++){
				m_kernelShapeFunction[i*sizeOfMatrix+j*sizeNeighbor+k]=mat_shapeMatrixOne(j,k);
				m_kernelDerivativeOneShapeFunctionOne[i*sizeOfMatrix+j*sizeNeighbor+k]=mat_shapeMatrixTwo(j,k);
				m_kernelDerivativeTwoShapeFunctionOne[i*sizeOfMatrix+j*sizeNeighbor+k]=mat_shapeMatrixThree(j,k);
				m_kernelDerivativeOneShapeFunctionTwo[i*sizeOfMatrix+j*sizeNeighbor+k]=-mat_shapeMatrixFour(j,k);
				m_kernelDerivativeTwoShapeFunctionTwo[i*sizeOfMatrix+j*sizeNeighbor+k]=-mat_shapeMatrixFive(j,k);
			}
		}
		//////////////////////////////////////////////////////////////////////////
		//释放临时变量空间
		delete[] weightVector;  // 加权矩阵
		delete[] weightVectorDerivativeOne;  //加权矩阵对第一个变量的导数
		delete[] weightVectorDerivativeTwo;  //加权矩阵对第二个变量的导数
		delete[] QTweight;
		delete[] QTweightDerivativeOne;
		delete[] QTweightDerivativeTwo;
		delete[] QTweightQ;
		delete[] QTweightQDerivativeOne;
		delete[] QTweightQDerivativeTwo;
		delete[] QTweightQInver;
	}
	//////////////////////////////////////////////////////////////////////////
	// 写文件 kershapeFunciton
	CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\kernelShapeFunction.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "a")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		fprintf(fpout,"%d\n",indexPointSet);
		fprintf(fpout,"/////////////////////////////////");
		for(int i=0;i<m_numOfInterpolatePoint;i++){
			for(int j=0;j<m_rankOfInterpolate;j++){
				for( int k=0;k<m_numOfNeighbors+1;k++){
					fprintf(fpout,"%f ",m_kernelShapeFunction[i*m_rankOfInterpolate*(m_numOfNeighbors+1)+j*(m_numOfNeighbors+1)+k]);
				}
				fprintf(fpout,"\n");  
			}
			fprintf(fpout,"//////////////////////////");
			fprintf(fpout,"\n");            
		}
		fclose(fpout);
	}

	//////////////////////////////////////////////////////////////////////////
	// 写文件 m_kernelDerivativeOneShapeFunctionOne
	filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\kernelDerivativeOneShapeFunctionOne.txt";

	if((fpout = fopen(filename_pw, "a")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		fprintf(fpout,"%d\n",indexPointSet);
		fprintf(fpout,"/////////////////////////////////");
			for(int i=0;i<m_numOfInterpolatePoint;i++){
				for(int j=0;j<m_rankOfInterpolate;j++){
					for( int k=0;k<m_numOfNeighbors+1;k++){
						fprintf(fpout,"%f ",m_kernelDerivativeOneShapeFunctionOne[i*m_rankOfInterpolate*(m_numOfNeighbors+1)+j*(m_numOfNeighbors+1)+k]);
					}
					fprintf(fpout,"\n");  
				}
				fprintf(fpout,"//////////////////////////");
					fprintf(fpout,"\n");            
			}
			fclose(fpout);
	}
	//////////////////////////////////////////////////////////////////////////
	/// 写文件 m_kernelDerivativeTwoShapeFunctionOne
	filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\kernelDerivativeTwoShapeFunctionOne.txt";

	if((fpout = fopen(filename_pw, "a")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		fprintf(fpout,"%d\n",indexPointSet);
		fprintf(fpout,"/////////////////////////////////");
			for(int i=0;i<m_numOfInterpolatePoint;i++){
				for(int j=0;j<m_rankOfInterpolate;j++){
					for( int k=0;k<m_numOfNeighbors+1;k++){
						fprintf(fpout,"%f ", m_kernelDerivativeTwoShapeFunctionOne[i*m_rankOfInterpolate*(m_numOfNeighbors+1)+j*(m_numOfNeighbors+1)+k]);
					}
					fprintf(fpout,"\n");  
				}
				fprintf(fpout,"//////////////////////////");
					fprintf(fpout,"\n");            
			}
			fclose(fpout);
	}

	//////////////////////////////////////////////////////////////////////////
	// 写文件 m_kernelDerivativeOneShapeFunctionTwo
	filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\kernelDerivativeOneShapeFunctionTwo.txt";

	if((fpout = fopen(filename_pw, "a")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		fprintf(fpout,"%d\n",indexPointSet);
		fprintf(fpout,"/////////////////////////////////");
			for(int i=0;i<m_numOfInterpolatePoint;i++){
				for(int j=0;j<m_rankOfInterpolate;j++){
					for( int k=0;k<m_numOfNeighbors+1;k++){
						fprintf(fpout,"%f ", m_kernelDerivativeOneShapeFunctionTwo[i*m_rankOfInterpolate*(m_numOfNeighbors+1)+j*(m_numOfNeighbors+1)+k]);
					}
					fprintf(fpout,"\n");  
				}
				fprintf(fpout,"//////////////////////////");
					fprintf(fpout,"\n");            
			}
			fclose(fpout);
	}

	//////////////////////////////////////////////////////////////////////////
	// 写文件 m_kernelDerivativeTwoShapeFunctionTwo

	filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\kernelDerivativeTwoShapeFunctionTwo.txt";

	if((fpout = fopen(filename_pw, "a")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		fprintf(fpout,"%d\n",indexPointSet);
		fprintf(fpout,"/////////////////////////////////");
			for(int i=0;i<m_numOfInterpolatePoint;i++){
				for(int j=0;j<m_rankOfInterpolate;j++){
					for( int k=0;k<m_numOfNeighbors+1;k++){
						fprintf(fpout,"%f ", m_kernelDerivativeTwoShapeFunctionTwo[i*m_rankOfInterpolate*(m_numOfNeighbors+1)+j*(m_numOfNeighbors+1)+k]);
					}
					fprintf(fpout,"\n");  
				}
				fprintf(fpout,"//////////////////////////");
					fprintf(fpout,"\n");            
			}
			fclose(fpout);
	}
	return;
}
void MeshlessFilter::ComputeShapeFunction(int indexPointSet)
{
	/*
	*	函数功能：计算形函数及其导数
	*  变量说明：
	*  int indexPointSet        索引点
	*  float* m_shapeFunction                各邻域点的形函数在各采样点的值
	*  float* m_shapeFunctionDerivativeOne   各邻域点的形函数对变量一的倒数在各采样点的值
	*  float* m_shapeFunctionDerivativeTwo   各邻域点的形函数对变量二的倒数在各采样点的值
	*/
	int sizeNeighbor=m_numOfNeighbors+1;  //改变邻域大小
	int sizeOfMatrix=sizeNeighbor*m_rankOfInterpolate; //形函数核矩阵的大小
	//////////////////////////////////////////////////////////////////////////
	// 分配内存
	m_shapeFunction=new float[sizeNeighbor*m_numOfInterpolatePoint]; //一行为一个邻域点的形函数在不同的积分点的值
	m_shapeFunctionDerivativeOne=new float[sizeNeighbor*m_numOfInterpolatePoint];
	m_shapeFunctionDerivativeTwo=new float[sizeNeighbor*m_numOfInterpolatePoint];
	//////////////////////////////////////////////////////////////////////////
	// 计算在不同的积分点插值函数及其导数的值
	float* localInterFunction;
	localInterFunction=new float[m_numOfInterpolatePoint*m_rankOfInterpolate];//  一行为每一个积分点处多项式各部分的值
	float* localInterFunctionDeOne;
	localInterFunctionDeOne=new float[m_numOfInterpolatePoint*m_rankOfInterpolate];
	float* localInterFunctionDeTwo;
	localInterFunctionDeTwo=new float[m_numOfInterpolatePoint*m_rankOfInterpolate];
	for(int i=0;i<m_numOfInterpolatePoint;i++){
		localInterFunction[i*m_rankOfInterpolate]=1;
		localInterFunction[i*m_rankOfInterpolate+1]=m_localInterpolatePoint[i][0]*m_localIntergralRadius;
		localInterFunction[i*m_rankOfInterpolate+2]=m_localInterpolatePoint[i][1]*m_localIntergralRadius;
		localInterFunction[i*m_rankOfInterpolate+3]=localInterFunction[i*m_rankOfInterpolate+1]*localInterFunction[i*m_rankOfInterpolate+1];
		localInterFunction[i*m_rankOfInterpolate+4]=localInterFunction[i*m_rankOfInterpolate+1]*localInterFunction[i*m_rankOfInterpolate+2];
		localInterFunction[i*m_rankOfInterpolate+5]=localInterFunction[i*m_rankOfInterpolate+2]*localInterFunction[i*m_rankOfInterpolate+2];
		localInterFunctionDeOne[i*m_rankOfInterpolate]=0;
		localInterFunctionDeOne[i*m_rankOfInterpolate+1]=1;
		localInterFunctionDeOne[i*m_rankOfInterpolate+2]=0;
		localInterFunctionDeOne[i*m_rankOfInterpolate+3]=2*localInterFunction[i*m_rankOfInterpolate+1];
		localInterFunctionDeOne[i*m_rankOfInterpolate+4]=localInterFunction[i*m_rankOfInterpolate+2];
		localInterFunctionDeOne[i*m_rankOfInterpolate+5]=0;
		localInterFunctionDeTwo[i*m_rankOfInterpolate]=0;
		localInterFunctionDeTwo[i*m_rankOfInterpolate+1]=0;
		localInterFunctionDeTwo[i*m_rankOfInterpolate+2]=1;
		localInterFunctionDeTwo[i*m_rankOfInterpolate+3]=0;
		localInterFunctionDeTwo[i*m_rankOfInterpolate+4]=localInterFunction[i*m_rankOfInterpolate+1];
		localInterFunctionDeTwo[i*m_rankOfInterpolate+5]=2*localInterFunction[i*m_rankOfInterpolate+2];
	}
	//////////////////////////////////////////////////////////////////////////
	// 循环计算形函数
	for(int i=0;i<sizeNeighbor;i++){
		for(int j=0;j<m_numOfInterpolatePoint;j++){
            m_shapeFunction[i*m_numOfInterpolatePoint+j]=0;
			m_shapeFunctionDerivativeOne[i*m_numOfInterpolatePoint+j]=0;
			m_shapeFunctionDerivativeTwo[i*m_numOfInterpolatePoint+j]=0;
			for(int k=0;k<m_rankOfInterpolate;k++){
				m_shapeFunction[i*m_numOfInterpolatePoint+j]+=localInterFunction[j*m_rankOfInterpolate+k]*m_kernelShapeFunction[j*sizeOfMatrix+k*sizeNeighbor+i];
				m_shapeFunctionDerivativeOne[i*m_numOfInterpolatePoint+j]+=
					localInterFunctionDeOne[j*m_rankOfInterpolate+k]*m_kernelShapeFunction[j*sizeOfMatrix+k*sizeNeighbor+i]
					+localInterFunction[j*m_rankOfInterpolate+k]*m_kernelDerivativeOneShapeFunctionOne[j*sizeOfMatrix+k*sizeNeighbor+i]
					+localInterFunction[j*m_rankOfInterpolate+k]*m_kernelDerivativeOneShapeFunctionTwo[j*sizeOfMatrix+k*sizeNeighbor+i];
				m_shapeFunctionDerivativeTwo[i*m_numOfInterpolatePoint+j]+=
					localInterFunctionDeTwo[j*m_rankOfInterpolate+k]*m_kernelShapeFunction[j*sizeOfMatrix+k*sizeNeighbor+i]
					+localInterFunction[j*m_rankOfInterpolate+k]*m_kernelDerivativeTwoShapeFunctionOne[j*sizeOfMatrix+k*sizeNeighbor+i]
					+localInterFunction[j*m_rankOfInterpolate+k]*m_kernelDerivativeTwoShapeFunctionTwo[j*sizeOfMatrix+k*sizeNeighbor+i];
            }
		}        
	}
	delete[] localInterFunctionDeOne;
	delete[] localInterFunction;	
	delete[]  localInterFunctionDeTwo;

	//////////////////////////////////////////////////////////////////////////
	/// 写文件 shapefunction
	CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\shapeFunction.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "a")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		fprintf(fpout,"%d\n",indexPointSet);
		fprintf(fpout,"/////////////////////////////////");
		for(int i=0;i<m_numOfNeighbors+1;i++){
			for(int j=0;j<m_numOfInterpolatePoint;j++){				
					fprintf(fpout,"%f ",m_shapeFunction[i*m_numOfInterpolatePoint+j]);				 
			}
			fprintf(fpout,"//////////////////////////");
			fprintf(fpout,"\n");            
		}
		fclose(fpout);
	}

	//////////////////////////////////////////////////////////////////////////
	// 写文件 m_shapeFunctionDerivativeOne
    filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\shapeFunctionDerivativeOne.txt";

	if((fpout = fopen(filename_pw, "a")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		fprintf(fpout,"%d\n",indexPointSet);
		fprintf(fpout,"/////////////////////////////////");
		for(int i=0;i<m_numOfNeighbors+1;i++){
			for(int j=0;j<m_numOfInterpolatePoint;j++){				
				fprintf(fpout,"%f ",m_shapeFunctionDerivativeOne[i*m_numOfInterpolatePoint+j]);				 
			}
			fprintf(fpout,"//////////////////////////");
			fprintf(fpout,"\n");            
		}
		fclose(fpout);
	}

	//////////////////////////////////////////////////////////////////////////
	// 写文件 m_shapeFunctionDerivativeTwo
	filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\shapeFunctionDerivativeTwo.txt";

	if((fpout = fopen(filename_pw, "a")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		fprintf(fpout,"%d\n",indexPointSet);
		fprintf(fpout,"/////////////////////////////////");
		for(int i=0;i<m_numOfNeighbors+1;i++){
			for(int j=0;j<m_numOfInterpolatePoint;j++){				
				fprintf(fpout,"%f ",m_shapeFunctionDerivativeTwo[i*m_numOfInterpolatePoint+j]);				 
			}
			fprintf(fpout,"//////////////////////////");
			fprintf(fpout,"\n");            
		}
		fclose(fpout);
	}

	return;
}
void MeshlessFilter::ComputeManifoldMetric(int indexPointSet)
{
	/*
	*	函数功能： 计算流形标准的代数行列式
	*  变量说明：
	*  int indexPointSet        索引点
	*  float* m_determinentMetric   在采样点处流形标准的代数行列式
	*/
	//////////////////////////////////////////////////////////////////////////
	// 分配内存
	m_determinentMetric=new float[m_numOfInterpolatePoint];
	m_metricMatrix=new float[m_numOfInterpolatePoint*4];
	m_metricMatrixInve=new float[m_numOfInterpolatePoint*4];
	//////////////////////////////////////////////////////////////////////////
	// 循环计算
	for(int i=0;i<m_numOfInterpolatePoint;i++){
		float tempVectorOne[3],tempVectorTwo[3];// 参考3.5
		for(int j=0;j<3;j++){
			tempVectorOne[j]=m_shapeFunctionDerivativeOne[i]*m_resultPointSet[indexPointSet*3+j];
			tempVectorTwo[j]=m_shapeFunctionDerivativeTwo[i]*m_resultPointSet[indexPointSet*3+j];			
		}
		std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(indexPointSet);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		for(int j=1;j<m_numOfNeighbors+1;j++){
			for(int k=0;k<3;k++){
				tempVectorOne[k]+=m_shapeFunctionDerivativeOne[j*m_numOfInterpolatePoint+i]*m_resultPointSet[(*vectorIterator)*3+k];
				tempVectorTwo[k]+=m_shapeFunctionDerivativeTwo[j*m_numOfInterpolatePoint+i]*m_resultPointSet[(*vectorIterator)*3+k];
			}
			vectorIterator++;
		}
		float localMetric[4];
		localMetric[0]=Vector3Vector(tempVectorOne,tempVectorOne);
		localMetric[1]=Vector3Vector(tempVectorOne,tempVectorTwo);
		localMetric[2]=localMetric[1];
		localMetric[3]=Vector3Vector(tempVectorTwo,tempVectorTwo);
		for(int j=0;j<4;j++){
			m_metricMatrix[i*4+j]=localMetric[j];
		}
		m_determinentMetric[i]=localMetric[0]*localMetric[3]-localMetric[1]*localMetric[2];
		m_metricMatrixInve[i*4]=localMetric[3]/m_determinentMetric[i];
		m_metricMatrixInve[i*4+1]=-localMetric[2]/m_determinentMetric[i];
		m_metricMatrixInve[i*4+2]=-localMetric[1]/m_determinentMetric[i];
		m_metricMatrixInve[i*4+3]=localMetric[0]/m_determinentMetric[i];

	}

	//////////////////////////////////////////////////////////////////////////
	/// 写文件 shapefunction
	CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\deteMetric.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "a")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		fprintf(fpout,"%d\n",indexPointSet);
		fprintf(fpout,"/////////////////////////////////");
		for(int i=0;i<m_numOfInterpolatePoint;i++){
					
				fprintf(fpout,"%f ",m_determinentMetric[i]);				 
	
			      
		}
		fprintf(fpout,"//////////////////////////");
		fprintf(fpout,"\n");      
		fclose(fpout);
	}
	return;
}
void MeshlessFilter::ComputeStiffMatrix(int indexPointSet)
{
	//////////////////////////////////////////////////////////////////////////
	//首先计算 对角线上的元素
	float namdaLoadConstant=m_namda*m_loadConstant;
	float elementOne=0;
	float elementTwo=0;
	float elementStiffMatrix;
	//////////////////////////////////////////////////////////////////////////
	// 求积分
	for(int i=0;i<m_numOfInterpolatePoint;i++){
		elementOne+=m_shapeFunction[i]*m_TestFunction[i]*sqrt(m_determinentMetric[i])*m_gaussItegral[i];
		elementTwo+=m_shapeFunctionDerivativeOne[i]*m_TestFunctionDerivativeOne[i]+m_shapeFunctionDerivativeTwo[i]*m_TestFunctionDerivativeTwo[i];
	}
	//////////////////////////////////////////////////////////////////////////
	// 刚度矩阵元素
	elementStiffMatrix=elementOne+m_namda*elementTwo;
	m_stiffMatrix[indexPointSet].m_elementStiffMarix.push_back(elementStiffMatrix);
	for(int i=0;i<3;i++){
		m_loadVector[indexPointSet*3+i]=namdaLoadConstant*elementOne*m_originalPointSet[indexPointSet*3+i]
		        +(1-namdaLoadConstant)*elementOne*m_resultPointSet[indexPointSet*3+i];

	}


	//////////////////////////////////////////////////////////////////////////
	//计算非对角线上的元素	计算荷载向量
	std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(indexPointSet);
	std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	for(int i=1;i<m_numOfNeighbors+1;i++){
		elementOne=0;
		elementTwo=0;
		//////////////////////////////////////////////////////////////////////////
		// 求积分
		for(int j=0;j<m_numOfInterpolatePoint;j++){
            elementOne+=m_shapeFunction[i*m_numOfInterpolatePoint+j]*m_TestFunction[j]*sqrt(m_determinentMetric[j])*m_gaussItegral[j];
			elementTwo+=m_shapeFunctionDerivativeOne[i*m_numOfInterpolatePoint+j]*m_TestFunctionDerivativeOne[j]+m_shapeFunctionDerivativeTwo[i*m_numOfInterpolatePoint+j]*m_TestFunctionDerivativeTwo[j];
		}
		elementStiffMatrix=elementOne+m_namda*elementTwo;
		m_stiffMatrix[indexPointSet].m_elementStiffMarix.push_back(elementStiffMatrix);
		for(int k=0;k<3;k++){
			m_loadVector[indexPointSet*3+k]+=namdaLoadConstant*elementOne*m_originalPointSet[3*(*vectorIterator)+k]
			+(1-namdaLoadConstant)*elementOne*m_resultPointSet[3*(*vectorIterator)+k];
		}
		vectorIterator++;		
		
	}
	//////////////////////////////////////////////////////////////////////////
	/// 写文件 刚度矩阵
	CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\stiffMatrix.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "a")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		fprintf(fpout,"%d\n",indexPointSet);
		fprintf(fpout,"/////////////////////////////////");
		std::vector<float>::iterator iteratorStiffMatrix=m_stiffMatrix[indexPointSet].m_elementStiffMarix.begin();
		for(;iteratorStiffMatrix!=m_stiffMatrix[indexPointSet].m_elementStiffMarix.end();iteratorStiffMatrix++){
            fprintf(fpout,"%f ",(*iteratorStiffMatrix));
		}		
		fprintf(fpout,"//////////////////////////");
		fprintf(fpout,"\n");      
		fclose(fpout);
	}
	return;
}
void MeshlessFilter::EstimateNormalDirection()
{
	//首先确定z值最大的那一点的法向
	int maxZpoint=0;
	for(int i=0;i<m_numOfPoints;i++){
		if(m_originalPointSet[i*3+2]>m_originalPointSet[maxZpoint*3+2])
			maxZpoint=i;

	}
	if(m_originalNormals[maxZpoint*3+2]<0){
		for(int i=0;i<3;i++){
			m_originalNormals[maxZpoint*3+i]=-m_originalNormals[maxZpoint*3+i];
		}
	}

	int* parent=new int[m_numOfPoints];//存放父节点

	std::multimap<float,int> valueKeySeries;//根据点的代价索引的map
	std::map<int,float> pointKeySeries;//根据点索引的map

	//初始化map
	for(int i=0;i<maxZpoint;i++){
		pointKeySeries.insert(std::map<int,float>::value_type(i,i+10));
		valueKeySeries.insert(std::map<float,int>::value_type(i+10,i));
	}
	for(int i=maxZpoint+1;i<m_numOfPoints;i++){
		pointKeySeries.insert(std::map<int,float>::value_type(i,i+10));
		valueKeySeries.insert(std::map<float,int>::value_type(i+10,i));
	}
	pointKeySeries.insert(std::map<int,float>::value_type(maxZpoint,0));
	valueKeySeries.insert(std::map<float,int>::value_type(0,maxZpoint));

	std::map<int,float>::iterator mapIterator = pointKeySeries.begin();
	std::multimap<float,int>::iterator multimapIterator=valueKeySeries.begin();

	while(multimapIterator!=valueKeySeries.end()){
		int pointExtract;
		if((*multimapIterator).first>1)
			break;
		pointExtract=(*multimapIterator).second;
		valueKeySeries.erase(multimapIterator);
		mapIterator=pointKeySeries.find(pointExtract);
		pointKeySeries.erase(mapIterator);
		std::map<int, KnearestField>::iterator mapNearestIterator=m_mapKnearest.find(pointExtract);
		std::vector<int>::iterator vectorIterator=(*mapNearestIterator).second.m_nearest.begin();


		for(;vectorIterator!=(*mapNearestIterator).second.m_nearest.end();vectorIterator++){
			int neighborPoint=(*vectorIterator);
			mapIterator=pointKeySeries.find(neighborPoint);
			if(mapIterator!=pointKeySeries.end()){
				multimapIterator=valueKeySeries.find((*mapIterator).second);
				if(neighborPoint!=(*multimapIterator).second){
					continue;
				}
				float relativeCost,normalDotProduct;
				float pointExtractNormal[3],neighborPointNormal[3];
				for(int j=0;j<3;j++){
					pointExtractNormal[j]=m_originalNormals[pointExtract*3+j];
					neighborPointNormal[j]=m_originalNormals[neighborPoint*3+j];
				}
				normalDotProduct=Vector3Vector(pointExtractNormal,neighborPointNormal);
				if(normalDotProduct<0){
					normalDotProduct=-normalDotProduct;
					for(int j=0;j<3;j++){
						m_originalNormals[neighborPoint*3+j]=-m_originalNormals[neighborPoint*3+j];
					}
				}
				relativeCost=1-normalDotProduct;
				if(relativeCost<(*mapIterator).second){
					parent[neighborPoint]=pointExtract;
					(*mapIterator).second=relativeCost;
					valueKeySeries.erase(multimapIterator);
					valueKeySeries.insert(std::map<float,int>::value_type(relativeCost,neighborPoint));
				}
			}
		}
		mapIterator=pointKeySeries.begin();
		multimapIterator=valueKeySeries.begin();
	}

	/*if(mapIterator!=pointKeySeries.end()){
	int numOfPointsLeave;
	float* pointSetLeave;
	numOfPointsLeave=pointKeySeries.size();
	pointSetLeave=new float[numOfPointsLeave*3];
	int i=0;
	while(mapIterator!=pointKeySeries.end()){
	int tempPointIndex;
	tempPointIndex=(*mapIterator).first;
	for(int j=0;j<3;j++){
	pointSetLeave[i*3+j]=m_originalPointSet[tempPointIndex*3+j];
	}			
	mapIterator++;
	i++;
	}

	pointKeySeries.clear();
	valueKeySeries.clear();
	delete[] parent;
	parent=NULL;

	EstimateNormalDirection(numOfPointsLeave,pointSetLeave,numOfEachPointNeighbor,allPointsNeighbor,allPointsNormal);
	delete[] pointSetLeave;
	}*/
	pointKeySeries.clear();
	valueKeySeries.clear();
	if(parent!=NULL)
		delete[] parent;
	return;
}
void MeshlessFilter::CalculateLinearSystem()
{
	
	//////////////////////////////////////////////////////////////////////////
	// 循环开始
	float threshold=m_ebusainu+1;
	while(threshold>=m_ebusainu){
		threshold=0.0;
		for(int i=0;i<m_numOfPoints;i++){
			std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
			std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
			std::vector<float>::iterator iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
			iteratorStiffMatrix++;
			float tempSum[3];
			for(int j=0;j<3;j++){
				tempSum[j]=0;
			}
			for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,iteratorStiffMatrix++){
				for(int j=0;j<3;j++){
					tempSum[j]+=(*iteratorStiffMatrix)*m_resultPointSet[(*vectorIterator)*3+j];
				}				
			}
			iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
			float localPoint[3];
			for(int j=0;j<3;j++){
				localPoint[j]=(m_loadVector[i*3+j]-tempSum[j])/(*iteratorStiffMatrix);
			}
			for(int j=0;j<3;j++){
				float localThreshold;
				localThreshold=abs(localPoint[j]-m_resultPointSet[i*3+j])/(1.0+abs(localPoint[j]));
				if(localThreshold>threshold)
					threshold=localThreshold;
				m_resultPointSet[i*3+j]=localPoint[j];
			}            
		}
	}
	return;
}