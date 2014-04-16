#include "stdafx.h"
#include "OneFreedomMeshless.h"
#include ".\mathlib\mathlib.h"
#include "Mcube.h"
#include "Matrix.h"
#include <math.h>
using namespace MATHLIB;

const int OneFreedomMeshless::m_numOfIntegralPoint=40;
const int OneFreedomMeshless::m_rankOfInterpolate=3;

const double OneFreedomMeshless::m_ordinateAndWeightCoff[120]={ 
	/*  
	*   第一列为半径为1时候的x坐标
	*   第二列为半径为1时候的y坐标
	*   第三列为两个一维高斯积分的加权值的乘积
	*/
	
    
	   -0.046546,-0.005837,0.023984,
		-0.037660,-0.027969,0.052688,
		-0.003759,-0.046759,0.074326,
		0.039334,-0.025562,0.085930,
		0.039334,0.025562,0.085930,
		-0.003759,0.046759,0.074326,
		-0.037660,0.027969,0.052688,
		-0.046546,0.005837,0.023984,
		-0.228972,-0.028714,0.048451,
		-0.185262,-0.137588,0.106438,
		-0.018490,-0.230023,0.150149,
		0.193496,-0.125745,0.173591,
		0.193496,0.125745,0.173591,
		-0.018490,0.230023,0.150149,
		-0.185262,0.137588,0.106438,
		-0.228972,0.028714,0.048451,
		-0.496114,-0.062215,0.057588,
		-0.401409,-0.298113,0.126510,
		-0.040063,-0.498392,0.178464,
		0.419249,-0.272453,0.206327,
		0.419249,0.272453,0.206327,
		-0.040063,0.498392,0.178464,
		-0.401409,0.298113,0.126510,
		-0.496114,0.062215,0.057588,
		-0.763256,-0.095716,0.048451,
		-0.617555,-0.458637,0.106438,
		-0.061636,-0.766761,0.150149,
		0.645001,-0.419161,0.173591,
		0.645001,0.419161,0.173591,
		-0.061636,0.766761,0.150149,
		-0.617555,0.458637,0.106438,
		-0.763256,0.095716,0.048451,
		-0.046546,-0.005837,0.023984,
		-0.037660,-0.027969,0.052688,
		-0.003759,-0.046759,0.074326,
		0.039334,-0.025562,0.085930,
		0.039334,0.025562,0.085930,
		-0.003759,0.046759,0.074326,
		-0.037660,0.027969,0.052688,
		-0.046546,0.005837,0.023984
};


//const double OneFreedomMeshless::m_localInterpolatePoint[9][2]={
//	{-0.7745966692,-0.7745966692},
//	{0,-0.7745966692},
//	{0.7745966692,-0.7745966692},
//	{-0.7745966692,0},
//	{0,0},
//	{0.7745966692,0},
//	{-0.7745966692,0.7745966692},
//	{0,0.7745966692},
//	{0.7745966692,0.7745966692}
//};
//
//const double OneFreedomMeshless::m_TestFunction[9]={
//	0.00646952576012,
//		0.09626480766995,
//		0.00646952576012,
//		0.09626480766995,
//		1.43239451226107,
//		0.09626480766995,
//		0.00646952576012,
//		0.09626480766995,
//		0.00646952576012
//};
//const double OneFreedomMeshless::m_TestFunctionDerivativeOne[9]={
//	0.04510145794586,
//		0,
//		-0.04510145794586,
//		0.67109759444089,
//		0,
//		-0.67109759444089,
//		0.04510145794586,
//		0,
//		-0.04510145794586
//};
//const double OneFreedomMeshless::m_TestFunctionDerivativeTwo[9]={
//	0.04510145794586, 
//		0.67109759444089,
//		0.04510145794586,
//		0,
//		0,
//		0,
//		-0.04510145794586,
//		-0.67109759444089,
//		-0.04510145794586
//};
//const double OneFreedomMeshless::m_gaussItegral[9]={
//	0.30864197535802469136,
//		0.49382716053950617284,
//		0.30864197535802469136,
//		0.49382716053950617284,
//		0.79012345680987654321,
//		0.49382716053950617284,
//		0.30864197535802469136,
//		0.49382716053950617284,
//		0.30864197535802469136
//};

OneFreedomMeshless::OneFreedomMeshless(void)
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
	m_localIntegralPoints=NULL;
	m_testFunction=NULL;
	m_testFunctionDerivativeOne=NULL;
	m_testFunctionDerivativeTwo=NULL;
}

OneFreedomMeshless::~OneFreedomMeshless(void)
{
	DeleteMeshlessFilter();
}

void  OneFreedomMeshless::GetMeshlessFilter(int numOfPointSet,float* pointSet,int kNearest, double radius, int ifFilterNormal,int interactiveTime,double meanShiftHposition,double meanShiftHnormal,double namda,	double loadConstant,	double ebusainu)
{
	/*
	函数功能：  接口函数，计算点集滤波
	变量说明：
	int numOfPointSet                        点的数量
	double* pointSet							 原始点集
	int kNearest                             k邻域的数量
	double radius                             k邻域搜索的半径
	int ifFilterNormal                       是否对法向进行滤波
	int interactiveTime                      滤波的交互次数
	double meanShiftHposition                法向滤波时的位置方差
	double meanShiftHnormal                  法向滤波时的法向方差
	double namda;                             时间步长
	double loadConstant;                      荷载常数
	double ebusainu;                          循环停止的条件；
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
	m_originalPointSet=new double[m_numOfPoints*3];
	m_resultPointSet=new double[m_numOfPoints*3];
	m_stiffMatrix=new StiffMatrix[m_numOfPoints];
	m_massMatrix=new StiffMatrix[m_numOfPoints];
	m_loadVector=new double[m_numOfPoints*3];
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
			if(j==1){
				int m;
				m=0;
			}
			ComputeNearestParameter(j);
			ComputeTestFunction(j);
			ComputeQmatrix(j);
			ComputeMatrixRelatedSampling(j);
			ComputeShapeFunction(j);
			ComputeManifoldMetric(j);
		//	ComputeGradientShapeFunction();
			
			ComputeStiffMatrix(j);

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
			if(m_metricMatrix!=NULL)
				delete[] m_metricMatrix;
			m_metricMatrix=NULL;
			if(m_metricMatrixInve!=NULL)
				delete[] m_metricMatrixInve;
			m_metricMatrixInve=NULL;
		
			if(m_localOrdinateVectorOne!=NULL)
				delete[] m_localOrdinateVectorOne;
			m_localOrdinateVectorOne=NULL;
			if(m_localOrdinateVectorTwo!=NULL)
				delete[] m_localOrdinateVectorTwo;
			m_localOrdinateVectorTwo=NULL;
			/*if(m_gradientShapeFunction!=NULL)
				delete[] m_gradientShapeFunction;
			m_gradientShapeFunction=NULL;
			if(m_gradientTestFunction!=NULL)
				delete[] m_gradientTestFunction;
			m_gradientTestFunction=NULL;*/
			if(m_testFunction!=NULL)
				delete[] m_testFunction;
			m_testFunction=NULL;
			if(m_testFunctionDerivativeOne!=NULL)
				delete[] m_testFunctionDerivativeOne;
			m_testFunctionDerivativeOne=NULL;
			if(m_testFunctionDerivativeTwo!=NULL)
				delete[] m_testFunctionDerivativeTwo;
			m_testFunctionDerivativeTwo=NULL;
			if(m_localIntegralPoints!=NULL)
				delete[] m_localIntegralPoints;
			m_localIntegralPoints=NULL;	

			//////////////////////////////////////////////////////////////////////////
			/*
			 *	释放内存
			 */
			
			///写文件for matlab
			//CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\matlabMatrix.txt";
			//FILE *fpout;
			//if((fpout = fopen(filename_pw, "a")) == NULL)
			//{
			//	int dkjkd;
			//	//MessageBox("can't open the file!");
			//}
			//else
			//{
			//	/*fprintf(fpout,"%d\n",indexPointSet);
			//	fprintf(fpout,"/////////////////////////////////");*/
			//	std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[j].m_elementStiffMarix.begin();
			//	fprintf(fpout,"%d %d %f\n",j+1,j+1,(*iteratorStiffMatrix));
			//	std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(j);
			//	std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
   //             for(iteratorStiffMatrix++;iteratorStiffMatrix!=m_stiffMatrix[j].m_elementStiffMarix.end();iteratorStiffMatrix++,vectorIterator++){
			//		fprintf(fpout,"%d %d %f\n",j+1,(*vectorIterator)+1,(*iteratorStiffMatrix));
			//	}		
			//	
			//	fclose(fpout);
			//}
		}
		//////////////////////////////////////////////////////////////////////////
		// 解线性系统
		AssembleStiffMatrix();
		CalculateLinearSystem();
		CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\matlabMatrix.txt";
		FILE *fpout;
		if((fpout = fopen(filename_pw, "w")) == NULL)
		{
			int dkjkd;
			//MessageBox("can't open the file!");
		}
		else{
			for(int j=0;j<m_numOfPoints;j++){
				std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[j].m_elementStiffMarix.begin();
				fprintf(fpout,"%d %d %f\n",j+1,j+1,(*iteratorStiffMatrix));
				std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(j);
				std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
				for(iteratorStiffMatrix++;iteratorStiffMatrix!=m_stiffMatrix[j].m_elementStiffMarix.end();iteratorStiffMatrix++,vectorIterator++){
					fprintf(fpout,"%d %d %f\n",j+1,(*vectorIterator)+1,(*iteratorStiffMatrix));
				}					
			}
			fclose(fpout);
		}
		
			
	//	CalculateLinearSystemGauss();
		//////////////////////////////////////////////////////////////////////////
		// 释放内存
		for(int j=0;j<m_numOfPoints;j++){
			m_stiffMatrix[j].m_elementStiffMarix.clear();
			m_massMatrix[j].m_elementStiffMarix.clear();
		}
	
		
	}
    //////////////////////////////////////////////////////////////////////////
    // 写文件 loadVector FOR matlab
	CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\matlabLoadVector.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f\n",m_loadVector[i*3],m_loadVector[i*3+1],m_loadVector[i*3+2]);

		}		
		fclose(fpout);
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
void OneFreedomMeshless:: DeleteMeshlessFilter()
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
void OneFreedomMeshless::ComputeMapKnearest()
{
	double radius=m_radius*m_radius;//存放近邻点到所求点的距离
	std::multimap<double ,int> mapKnearest;//对每一个点建立邻域vector时建立临时map
	int numOfKnearest;//统计一个点的临时邻域大小

	//////////////////////////////////////////////////////////////////////////
	//遍历所有点，并且按照距离由小到大存放存放到map中
	for(int i=0;i<m_numOfPoints;i++){
		numOfKnearest=0;
		for(int j=0;j<m_numOfPoints;j++){
			if(j==i)
				continue;
			double distancePointToPoint;
			distancePointToPoint=(m_originalPointSet[i*3]-m_originalPointSet[j*3])*(m_originalPointSet[i*3]-m_originalPointSet[j*3])
				+(m_originalPointSet[i*3+1]-m_originalPointSet[j*3+1])*(m_originalPointSet[i*3+1]-m_originalPointSet[j*3+1])
				+(m_originalPointSet[i*3+2]-m_originalPointSet[j*3+2])*(m_originalPointSet[i*3+2]-m_originalPointSet[j*3+2]);
			if(distancePointToPoint>radius)
				continue;
			mapKnearest.insert(std::multimap<double,int>::value_type(distancePointToPoint,j));
			numOfKnearest+=1;
			if(numOfKnearest>m_kNearest){
				std::multimap<double,int>::iterator mapIterator = mapKnearest.end();
				mapIterator--;
				double w;
				w=(*mapIterator).first;
				mapKnearest.erase(mapIterator);
				numOfKnearest-=1;
			}
		}
		KnearestField fieldKnearest;//临时结构
		fieldKnearest.m_numOfNearest=numOfKnearest;
		fieldKnearest.m_IfBoundary=FALSE;
		std::multimap<double,int>::iterator mapIterator = mapKnearest.begin();
		for(;mapIterator!=mapKnearest.end();mapIterator++){
			fieldKnearest.m_nearest.push_back((*mapIterator).second);
		}
		m_mapKnearest.insert(std::map<int, KnearestField>::value_type(i,fieldKnearest));

		mapKnearest.clear();
	}	
	//////////////////////////////////////////////////////////////////////////
	//// 测试代码
	//CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\knearest.txt";
	//FILE *fpout;
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(i);
	//		std::vector<int>::iterator vectorIterator=(*mapKnearestIterator).second.m_nearest.begin();
	//		fprintf(fpout,"%d ",i);
	//		for(;vectorIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorIterator++){
	//			fprintf(fpout,"%d ",(*vectorIterator));
	//		}
	//		fprintf(fpout,"\n");
	//	}
	//	fclose(fpout);
	//}
}
void OneFreedomMeshless::ComputeNormal()
{
	double centroidPosition[3];//重心位置
	double* localVariationVector;//重心到点的位置的向量
	int numOfNearest;

	mat_f8 mat_covaMatrix(3, 3);
	if(m_originalNormals!=NULL)
		delete[] m_originalNormals;
	m_originalNormals=new double[m_numOfPoints*3];
	m_originalMajorDirection=new double[m_numOfPoints*3];
	m_originalMinorDirection=new double[m_numOfPoints*3];
	m_gaussCurvature=new double[m_numOfPoints*2];
	//if(m_meanAreas!=NULL)
	//	delete[] m_meanAreas;
	//m_meanAreas=new double[m_numOfPoints];
	//if(m_meanCurvature!=NULL)
	//	delete[] m_meanCurvature;
	//m_meanCurvature=new double[m_numOfPoints];
	//if(m_resultColor!=NULL)
	//	delete[] m_resultColor;
	//m_resultColor=new double[m_numOfPoints*3];


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
		localVariationVector=new double[(numOfNearest+1)*3];
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
	//CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\normal.txt";
	//FILE *fpout;
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f\n",m_originalNormals[i*3],m_originalNormals[i*3+1],m_originalNormals[i*3+2]);

	//	}
	//	fclose(fpout);
	//}

	////////////////////////////////////////////////////////////////////////////
	//// 写主方向文件
	//filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\majorDirection.txt";

	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f\n",m_originalMajorDirection[i*3],m_originalMajorDirection[i*3+1],m_originalMajorDirection[i*3+2]);

	//	}
	//	fclose(fpout);
	//}

	////////////////////////////////////////////////////////////////////////////
	//// 写次方向文件
	//filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\minorDirection.txt";

	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f\n",m_originalMinorDirection[i*3],m_originalMinorDirection[i*3+1],m_originalMinorDirection[i*3+2]);

	//	}
	//	fclose(fpout);
	//}


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
void OneFreedomMeshless::ComputeNearestParameter(int indexPointSet)
{
	//int numOfNearest;
	double localPoint[3], localNormal[3],localMajorDirection[3],localMinorDirection[3];

	
	//if(m_localOrdinate!=NULL)
	//	delete[] m_localOrdinate;
	//m_localOrdinate=NULL;
	//m_localOrdinate=new double[m_numOfNeighbors*3];
	if(m_localIntegralPoints!=NULL)
		delete[] m_localIntegralPoints;
	m_localIntegralPoints=NULL;
	m_localIntegralPoints=new double[m_numOfIntegralPoint*2];
	m_localIntergralRadius=0;
	for(int i=0;i<3;i++){
		localMajorDirection[i]=m_originalMajorDirection[indexPointSet*3+i];
		localMinorDirection[i]=m_originalMinorDirection[indexPointSet*3+i];
		localNormal[i]=m_originalNormals[indexPointSet*3+i];
	}
	std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(indexPointSet);
	std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	m_numOfNeighbors=(*mapIterator).second.m_nearest.size();
	if(m_localOrdinate!=NULL)
		delete[] m_localOrdinate;
	int temp=(*mapIterator).second.m_nearest.size();
	m_localOrdinate=NULL;
	m_localOrdinate=new double[m_numOfNeighbors*3];
	m_localWeightRadiusSquare=0;
	for(int i=0;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,i++){
		double vectorNeighborToPoint[3];
		double localRadius;
		double tempRadiusSquare;
		for(int j=0;j<3;j++){
			vectorNeighborToPoint[j]=m_resultPointSet[(*vectorIterator)*3+j]-m_resultPointSet[indexPointSet*3+j];
		}
		m_localOrdinate[i*3]=Vector3Vector(vectorNeighborToPoint,localMajorDirection);
		m_localOrdinate[i*3+1]=Vector3Vector(vectorNeighborToPoint,localMinorDirection); 
		m_localOrdinate[i*3+2]=Vector3Vector(vectorNeighborToPoint,localNormal);
		tempRadiusSquare=m_localOrdinate[i*3]*m_localOrdinate[i*3]+m_localOrdinate[i*3+1]*m_localOrdinate[i*3+1];
		if(tempRadiusSquare>m_localWeightRadiusSquare)
			m_localWeightRadiusSquare=tempRadiusSquare;

	}
	//m_localIntergralRadiusSquare=m_localOrdinate[0]*m_localOrdinate[0]+m_localOrdinate[1]*m_localOrdinate[1];
	//m_localWeightRadiusSquare=m_localOrdinate[(m_numOfNeighbors-1)*3]*m_localOrdinate[(m_numOfNeighbors-1)*3]
	//	+m_localOrdinate[(m_numOfNeighbors-1)*3+1]*m_localOrdinate[(m_numOfNeighbors-1)*3+1];
	//
	//if(m_localIntergralRadiusSquare>m_localWeightRadiusSquare)
	//	m_localWeightRadiusSquare=m_localIntergralRadiusSquare;
	m_localWeightRadius=sqrt(m_localWeightRadiusSquare);
	

	m_localIntergralRadius=m_localWeightRadius/2;
	m_localWeightRadiusC=m_localWeightRadius/2;
	m_localIntergralRadiusC=m_localIntergralRadius/2;

	m_localIntergralRadiusSquare=m_localIntergralRadius*m_localIntergralRadius;
	m_localWeightRadiusCSquare=m_localWeightRadiusC*m_localWeightRadiusC;
	m_localIntergralRadiusCSquare=m_localIntergralRadiusC*m_localIntergralRadiusC;

	//////////////////////////////////////////////////////////////////////////
	/*
	 *	计算积分点的局部坐标
	 */
	//////////////////////////////////////////////////////////////////////////
	for (int i=0;i<m_numOfIntegralPoint;i++){
		m_localIntegralPoints[i*2]=m_ordinateAndWeightCoff[i*3]*m_localIntergralRadius;
		m_localIntegralPoints[i*2+1]=m_ordinateAndWeightCoff[i*3+1]*m_localIntergralRadius;		
	}

	//////////////////////////////////////////////////////////////////////////
	/*
	 *	写局部坐标文件
	 */
	//CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\localOrdinalte.txt";
	//FILE *fpout;
	//if((fpout = fopen(filename_pw, "a")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	fprintf(fpout,"%d\n",indexPointSet);
	//	for(int i=0;i<m_numOfNeighbors;i++){
	//		fprintf(fpout,"%f %f\n",m_localOrdinate[i*3],m_localOrdinate[i*3+1]);

	//	}
	//	fclose(fpout);
	//}


	return;


}
void OneFreedomMeshless::ComputeQmatrix(int indexPointSet)
{
	m_localQmatrix=new double[m_rankOfInterpolate*(m_numOfNeighbors+1)];//按行排列
	m_localQTmatrix=new double[m_rankOfInterpolate*(m_numOfNeighbors+1)];//按列排列
	m_localQmatrix[0]=1;
	for(int i=1;i<m_rankOfInterpolate;i++){
		m_localQmatrix[i]=0;
	}
	for(int i=1;i<m_numOfNeighbors+1;i++){
		m_localQmatrix[i*m_rankOfInterpolate]=1;
		m_localQmatrix[i*m_rankOfInterpolate+1]=m_localOrdinate[(i-1)*3];
		m_localQmatrix[i*m_rankOfInterpolate+2]=m_localOrdinate[(i-1)*3+1];
		//m_localQmatrix[i*m_rankOfInterpolate+3]=m_localOrdinate[(i-1)*3]*m_localOrdinate[(i-1)*3];
		//m_localQmatrix[i*m_rankOfInterpolate+4]=m_localOrdinate[(i-1)*3]*m_localOrdinate[(i-1)*3+1];
		//m_localQmatrix[i*m_rankOfInterpolate+5]=m_localOrdinate[(i-1)*3+1]*m_localOrdinate[(i-1)*3+1];

	}
	for(int i=0;i<m_numOfNeighbors+1;i++){
		for(int j=0;j<m_rankOfInterpolate;j++){
			m_localQTmatrix[j*(m_numOfNeighbors+1)+i]=m_localQmatrix[i*m_rankOfInterpolate+j];
		}
	}
	//////////////////////////////////////////////////////////////////////////
	// 写文件
	//CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\localQmatrix.txt";
	//FILE *fpout;
	//if((fpout = fopen(filename_pw, "a")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	fprintf(fpout,"%d\n",indexPointSet);
	//	for(int i=0;i<m_numOfNeighbors+1;i++){
	//		for(int j=0;j<m_rankOfInterpolate;j++){
	//			fprintf(fpout,"%f ",m_localQmatrix[i*m_rankOfInterpolate+j]);
	//		}
	//		fprintf(fpout," \n");            
	//	}
	//	fclose(fpout);
	//}
	return;    
}
void OneFreedomMeshless::ComputeMatrixRelatedSampling(int indexPointSet)
{
	//////////////////////////////////////////////////////////////////////////
	/*
	*	函数功能： 计算与采样点有关的矩阵
	*  变量说明： 
	*  int indexPointSet        索引点
	*  double* m_kernelShapeFunction       
	*  double* m_kernelDerivativeOneShapeFunctionOne
	*  double* m_kernelDerivativeTwoShapeFunctionOne
	*  double* m_kernelDerivativeOneShapeFunctionTwo
	*  double* m_kernelDerivativeTwoShapeFunctionTwo
	*/
	//////////////////////////////////////////////////////////////////////////
	// 分配内存空间
	if(m_kernelShapeFunction!=NULL)
		delete[] m_kernelShapeFunction;
	m_kernelShapeFunction=new double[m_rankOfInterpolate*m_numOfIntegralPoint*(m_numOfNeighbors+1)];
	if(m_kernelDerivativeOneShapeFunctionOne!=NULL)
		delete[] m_kernelDerivativeOneShapeFunctionOne;
	m_kernelDerivativeOneShapeFunctionOne=new double[m_rankOfInterpolate*m_numOfIntegralPoint*(m_numOfNeighbors+1)];
	if(m_kernelDerivativeTwoShapeFunctionOne!=NULL)
		delete[] m_kernelDerivativeTwoShapeFunctionOne;
	m_kernelDerivativeTwoShapeFunctionOne=new double[m_rankOfInterpolate*m_numOfIntegralPoint*(m_numOfNeighbors+1)];
	if(m_kernelDerivativeOneShapeFunctionTwo!=NULL)
		delete[] m_kernelDerivativeOneShapeFunctionTwo;
	m_kernelDerivativeOneShapeFunctionTwo=new double[m_rankOfInterpolate*m_numOfIntegralPoint*(m_numOfNeighbors+1)];
	if(m_kernelDerivativeTwoShapeFunctionTwo!=NULL)
		delete[] m_kernelDerivativeTwoShapeFunctionTwo;
	m_kernelDerivativeTwoShapeFunctionTwo=new double[m_rankOfInterpolate*m_numOfIntegralPoint*(m_numOfNeighbors+1)];
	//////////////////////////////////////////////////////////////////////////
	/* 
	 * 计算积分中用到的常量
	 */
	//////////////////////////////////////////////////////////////////////////
	double tempEXP,tempOneMinusEXP,temp;
	temp=-4;
	tempEXP=exp(temp);
	tempOneMinusEXP=1-tempEXP;

	//////////////////////////////////////////////////////////////////////////
	// 开始循环计算
	for(int i=0;i<m_numOfIntegralPoint;i++){
		//////////////////////////////////////////////////////////////////////////
		// 定义临时变量
		double* weightVector=new double[m_numOfNeighbors+1];  // 加权矩阵
		double* weightVectorDerivativeOne=new double[m_numOfNeighbors+1];  //加权矩阵对第一个变量的导数
		double* weightVectorDerivativeTwo=new double[m_numOfNeighbors+1];  //加权矩阵对第二个变量的导数
		double* QTweight=new double[m_rankOfInterpolate*(m_numOfNeighbors+1)];
		double* QTweightDerivativeOne=new double[m_rankOfInterpolate*(m_numOfNeighbors+1)];
		double* QTweightDerivativeTwo=new double[m_rankOfInterpolate*(m_numOfNeighbors+1)];
		double* QTweightQ=new double[m_rankOfInterpolate*m_rankOfInterpolate];
		double* QTweightQDerivativeOne=new double[m_rankOfInterpolate*m_rankOfInterpolate];
		double* QTweightQDerivativeTwo=new double[m_rankOfInterpolate*m_rankOfInterpolate];
		double* QTweightQInver=new double[m_rankOfInterpolate*m_rankOfInterpolate];

		//////////////////////////////////////////////////////////////////////////
		// 计算加权矩阵及其导数阵
		//////////////////////////////////////////////////////////////////////////
	    double tempDistanceSquare;//积分点与邻域点之间的距离的平方
		tempDistanceSquare=m_localIntegralPoints[i*2]*m_localIntegralPoints[i*2]
							+m_localIntegralPoints[i*2+1]*m_localIntegralPoints[i*2+1];
		if(tempDistanceSquare>m_localWeightRadiusSquare){
			weightVector[0]=0;
			weightVectorDerivativeOne[0]=0;
			weightVectorDerivativeTwo[0]=0;
		}
		else{
			double tempEXPDistance;
			tempEXPDistance=exp(-tempDistanceSquare/m_localWeightRadiusCSquare);
			weightVector[0]=(tempEXPDistance-tempEXP)/tempOneMinusEXP;
			weightVectorDerivativeOne[0]=-2*tempEXPDistance*m_localIntegralPoints[i*2]/m_localWeightRadiusCSquare/tempOneMinusEXP;
			weightVectorDerivativeTwo[0]=-2*tempEXPDistance*m_localIntegralPoints[i*2+1]/m_localWeightRadiusCSquare/tempOneMinusEXP;
            
		}

		for(int j=0;j<m_numOfNeighbors;j++){
			double localOrdinateMinus[2];
			localOrdinateMinus[0]=m_localOrdinate[j*3]-m_localIntegralPoints[i*2];
			localOrdinateMinus[1]=m_localOrdinate[j*3+1]-m_localIntegralPoints[i*2+1];
			tempDistanceSquare=localOrdinateMinus[0]*localOrdinateMinus[0]
					+localOrdinateMinus[1]*localOrdinateMinus[1];
			if(tempDistanceSquare>m_localWeightRadiusSquare){
				weightVector[j+1]=0;
				weightVectorDerivativeOne[j+1]=0;
				weightVectorDerivativeTwo[j+1]=0;
			}
			else{
				double tempEXPDistance;
				tempEXPDistance=exp(-tempDistanceSquare/m_localWeightRadiusCSquare);
				weightVector[j+1]=(tempEXPDistance-tempEXP)/tempOneMinusEXP;
				weightVectorDerivativeOne[j+1]=2*tempEXPDistance*localOrdinateMinus[0]/m_localWeightRadiusCSquare/tempOneMinusEXP;
				weightVectorDerivativeTwo[j+1]=2*tempEXPDistance*localOrdinateMinus[1]/m_localWeightRadiusCSquare/tempOneMinusEXP;
			}
		}

		
		//////////////////////////////////////////////////////////////////////////
		/*
		 *	计算矩阵 QTweight  QTweightDerivativeOne  QTweightDerivativeTwo
		 */  
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
	//CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\kernelShapeFunction.txt";
	//FILE *fpout;
	//if((fpout = fopen(filename_pw, "a")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	fprintf(fpout,"%d\n",indexPointSet);
	//	fprintf(fpout,"/////////////////////////////////");
	//	for(int i=0;i<m_numOfIntegralPoint;i++){
	//		for(int j=0;j<m_rankOfInterpolate;j++){
	//			for( int k=0;k<m_numOfNeighbors+1;k++){
	//				fprintf(fpout,"%f ",m_kernelShapeFunction[i*m_rankOfInterpolate*(m_numOfNeighbors+1)+j*(m_numOfNeighbors+1)+k]);
	//			}
	//			fprintf(fpout,"\n");  
	//		}
	//		fprintf(fpout,"//////////////////////////");
	//		fprintf(fpout,"\n");            
	//	}
	//	fclose(fpout);
	//}

	////////////////////////////////////////////////////////////////////////////
	//// 写文件 m_kernelDerivativeOneShapeFunctionOne
	//filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\kernelDerivativeOneShapeFunctionOne.txt";

	//if((fpout = fopen(filename_pw, "a")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	fprintf(fpout,"%d\n",indexPointSet);
	//	fprintf(fpout,"/////////////////////////////////");
	//	for(int i=0;i<m_numOfIntegralPoint;i++){
	//		for(int j=0;j<m_rankOfInterpolate;j++){
	//			for( int k=0;k<m_numOfNeighbors+1;k++){
	//				fprintf(fpout,"%f ",m_kernelDerivativeOneShapeFunctionOne[i*m_rankOfInterpolate*(m_numOfNeighbors+1)+j*(m_numOfNeighbors+1)+k]);
	//			}
	//			fprintf(fpout,"\n");  
	//		}
	//		fprintf(fpout,"//////////////////////////");
	//		fprintf(fpout,"\n");            
	//	}
	//	fclose(fpout);
	//}
	////////////////////////////////////////////////////////////////////////////
	///// 写文件 m_kernelDerivativeTwoShapeFunctionOne
	//filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\kernelDerivativeTwoShapeFunctionOne.txt";

	//if((fpout = fopen(filename_pw, "a")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	fprintf(fpout,"%d\n",indexPointSet);
	//	fprintf(fpout,"/////////////////////////////////");
	//	for(int i=0;i<m_numOfIntegralPoint;i++){
	//		for(int j=0;j<m_rankOfInterpolate;j++){
	//			for( int k=0;k<m_numOfNeighbors+1;k++){
	//				fprintf(fpout,"%f ", m_kernelDerivativeTwoShapeFunctionOne[i*m_rankOfInterpolate*(m_numOfNeighbors+1)+j*(m_numOfNeighbors+1)+k]);
	//			}
	//			fprintf(fpout,"\n");  
	//		}
	//		fprintf(fpout,"//////////////////////////");
	//		fprintf(fpout,"\n");            
	//	}
	//	fclose(fpout);
	//}

	////////////////////////////////////////////////////////////////////////////
	//// 写文件 m_kernelDerivativeOneShapeFunctionTwo
	//filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\kernelDerivativeOneShapeFunctionTwo.txt";

	//if((fpout = fopen(filename_pw, "a")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	fprintf(fpout,"%d\n",indexPointSet);
	//	fprintf(fpout,"/////////////////////////////////");
	//	for(int i=0;i<m_numOfIntegralPoint;i++){
	//		for(int j=0;j<m_rankOfInterpolate;j++){
	//			for( int k=0;k<m_numOfNeighbors+1;k++){
	//				fprintf(fpout,"%f ", m_kernelDerivativeOneShapeFunctionTwo[i*m_rankOfInterpolate*(m_numOfNeighbors+1)+j*(m_numOfNeighbors+1)+k]);
	//			}
	//			fprintf(fpout,"\n");  
	//		}
	//		fprintf(fpout,"//////////////////////////");
	//		fprintf(fpout,"\n");            
	//	}
	//	fclose(fpout);
	//}

	////////////////////////////////////////////////////////////////////////////
	//// 写文件 m_kernelDerivativeTwoShapeFunctionTwo

	//filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\kernelDerivativeTwoShapeFunctionTwo.txt";

	//if((fpout = fopen(filename_pw, "a")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	fprintf(fpout,"%d\n",indexPointSet);
	//	fprintf(fpout,"/////////////////////////////////");
	//	for(int i=0;i<m_numOfIntegralPoint;i++){
	//		for(int j=0;j<m_rankOfInterpolate;j++){
	//			for( int k=0;k<m_numOfNeighbors+1;k++){
	//				fprintf(fpout,"%f ", m_kernelDerivativeTwoShapeFunctionTwo[i*m_rankOfInterpolate*(m_numOfNeighbors+1)+j*(m_numOfNeighbors+1)+k]);
	//			}
	//			fprintf(fpout,"\n");  
	//		}
	//		fprintf(fpout,"//////////////////////////");
	//		fprintf(fpout,"\n");            
	//	}
	//	fclose(fpout);
	//}
	return;
}

void OneFreedomMeshless::ComputeTestFunction(int indexPointSet)
{
	/*
	*	计算测试函数在积分点的函数值及其导数值
	*/
	if(m_testFunction!=NULL)
		delete[] m_testFunction;
	m_testFunction=NULL;
	if(m_testFunctionDerivativeOne!=NULL)
		delete[] m_testFunctionDerivativeOne;
	m_testFunctionDerivativeOne=NULL;
	if(m_testFunctionDerivativeTwo!=NULL)
		delete[] m_testFunctionDerivativeTwo;
	m_testFunctionDerivativeTwo=NULL;
	m_testFunction=new double[m_numOfIntegralPoint];
	m_testFunctionDerivativeOne=new double[m_numOfIntegralPoint];
	m_testFunctionDerivativeTwo=new double[m_numOfIntegralPoint];

	double tempEXP,tempOneMinusEXP;// 计算测试函数时候的一些常数
	double temp=-4;
	tempEXP=exp(temp);
	tempOneMinusEXP=1-tempEXP;

	for(int i=0;i<m_numOfIntegralPoint;i++){
		double tempDistanceSquare;//积分点与原点之间的距离的平方
		tempDistanceSquare=m_localIntegralPoints[i*2]*m_localIntegralPoints[i*2]
				+m_localIntegralPoints[i*2+1]*m_localIntegralPoints[i*2+1];
		if(tempDistanceSquare>m_localIntergralRadiusSquare){ 
			m_testFunction[i]=0;
			m_testFunctionDerivativeOne[i]=0;
			m_testFunctionDerivativeTwo[i]=0;
		
		}
		else{
			double tempEXPDistance;
			tempEXPDistance=exp(-tempDistanceSquare/m_localIntergralRadiusCSquare);
			m_testFunction[i]=(tempEXPDistance-tempEXP)/tempOneMinusEXP;
			m_testFunctionDerivativeOne[i]=-2*tempEXPDistance*m_localIntegralPoints[i*2]/m_localIntergralRadiusCSquare/tempOneMinusEXP;
			m_testFunctionDerivativeTwo[i]=-2*tempEXPDistance*m_localIntegralPoints[i*2+1]/m_localIntergralRadiusCSquare/tempOneMinusEXP;
		}

	}
	return;
}
void OneFreedomMeshless::ComputeShapeFunction(int indexPointSet)
{
	/*
	*	函数功能：计算形函数及其导数
	*  变量说明：
	*  int indexPointSet        索引点
	*  double* m_shapeFunction                各邻域点的形函数在各采样点的值
	*  double* m_shapeFunctionDerivativeOne   各邻域点的形函数对变量一的倒数在各采样点的值
	*  double* m_shapeFunctionDerivativeTwo   各邻域点的形函数对变量二的倒数在各采样点的值
	*/
	int sizeNeighbor=m_numOfNeighbors+1;  //改变邻域大小
	int sizeOfMatrix=sizeNeighbor*m_rankOfInterpolate; //形函数核矩阵的大小
	//////////////////////////////////////////////////////////////////////////
	/*
	 *	分配内存
	 */
	m_shapeFunction=new double[sizeNeighbor*m_numOfIntegralPoint]; //一行为一个邻域点的形函数在不同的积分点的值
	m_shapeFunctionDerivativeOne=new double[sizeNeighbor*m_numOfIntegralPoint];
	m_shapeFunctionDerivativeTwo=new double[sizeNeighbor*m_numOfIntegralPoint];
	//////////////////////////////////////////////////////////////////////////
	// 计算在不同的积分点插值函数及其导数的值
	double* localInterFunction;
	localInterFunction=new double[m_numOfIntegralPoint*m_rankOfInterpolate];//  一行为每一个积分点处多项式各部分的值
	double* localInterFunctionDeOne;
	localInterFunctionDeOne=new double[m_numOfIntegralPoint*m_rankOfInterpolate];
	double* localInterFunctionDeTwo;
	localInterFunctionDeTwo=new double[m_numOfIntegralPoint*m_rankOfInterpolate];
	for(int i=0;i<m_numOfIntegralPoint;i++){
		localInterFunction[i*m_rankOfInterpolate]=1;
		localInterFunction[i*m_rankOfInterpolate+1]=m_localIntegralPoints[i*2];
		localInterFunction[i*m_rankOfInterpolate+2]=m_localIntegralPoints[i*2+1];
		/*localInterFunction[i*m_rankOfInterpolate+3]=localInterFunction[i*m_rankOfInterpolate+1]*localInterFunction[i*m_rankOfInterpolate+1];
		localInterFunction[i*m_rankOfInterpolate+4]=localInterFunction[i*m_rankOfInterpolate+1]*localInterFunction[i*m_rankOfInterpolate+2];
		localInterFunction[i*m_rankOfInterpolate+5]=localInterFunction[i*m_rankOfInterpolate+2]*localInterFunction[i*m_rankOfInterpolate+2];*/
		localInterFunctionDeOne[i*m_rankOfInterpolate]=0;
		localInterFunctionDeOne[i*m_rankOfInterpolate+1]=1;
		localInterFunctionDeOne[i*m_rankOfInterpolate+2]=0;
		/*localInterFunctionDeOne[i*m_rankOfInterpolate+3]=2*localInterFunction[i*m_rankOfInterpolate+1];
		localInterFunctionDeOne[i*m_rankOfInterpolate+4]=localInterFunction[i*m_rankOfInterpolate+2];
		localInterFunctionDeOne[i*m_rankOfInterpolate+5]=0;*/
		localInterFunctionDeTwo[i*m_rankOfInterpolate]=0;
		localInterFunctionDeTwo[i*m_rankOfInterpolate+1]=0;
		localInterFunctionDeTwo[i*m_rankOfInterpolate+2]=1;
		/*localInterFunctionDeTwo[i*m_rankOfInterpolate+3]=0;
		localInterFunctionDeTwo[i*m_rankOfInterpolate+4]=localInterFunction[i*m_rankOfInterpolate+1];
		localInterFunctionDeTwo[i*m_rankOfInterpolate+5]=2*localInterFunction[i*m_rankOfInterpolate+2];*/
	}
	//////////////////////////////////////////////////////////////////////////
	// 循环计算形函数
	for(int i=0;i<sizeNeighbor;i++){
		for(int j=0;j<m_numOfIntegralPoint;j++){
			m_shapeFunction[i*m_numOfIntegralPoint+j]=0;
			m_shapeFunctionDerivativeOne[i*m_numOfIntegralPoint+j]=0;
			m_shapeFunctionDerivativeTwo[i*m_numOfIntegralPoint+j]=0;
			for(int k=0;k<m_rankOfInterpolate;k++){
				m_shapeFunction[i*m_numOfIntegralPoint+j]+=localInterFunction[j*m_rankOfInterpolate+k]*m_kernelShapeFunction[j*sizeOfMatrix+k*sizeNeighbor+i];
				m_shapeFunctionDerivativeOne[i*m_numOfIntegralPoint+j]+=
					localInterFunctionDeOne[j*m_rankOfInterpolate+k]*m_kernelShapeFunction[j*sizeOfMatrix+k*sizeNeighbor+i]
					+localInterFunction[j*m_rankOfInterpolate+k]*m_kernelDerivativeOneShapeFunctionOne[j*sizeOfMatrix+k*sizeNeighbor+i]
					+localInterFunction[j*m_rankOfInterpolate+k]*m_kernelDerivativeOneShapeFunctionTwo[j*sizeOfMatrix+k*sizeNeighbor+i];
					m_shapeFunctionDerivativeTwo[i*m_numOfIntegralPoint+j]+=
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
	//CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\shapeFunction.txt";
	//FILE *fpout;
	//if((fpout = fopen(filename_pw, "a")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	fprintf(fpout,"%d\n",indexPointSet);
	//	fprintf(fpout,"/////////////////////////////////");
	//	for(int i=0;i<m_numOfNeighbors+1;i++){
	//		for(int j=0;j<m_numOfIntegralPoint;j++){				
	//			fprintf(fpout,"%f ",m_shapeFunction[i*m_numOfIntegralPoint+j]);				 
	//		}
	//		fprintf(fpout,"//////////////////////////");
	//		fprintf(fpout,"\n");            
	//	}
	//	fclose(fpout);
	//}

	////////////////////////////////////////////////////////////////////////////
	//// 写文件 m_shapeFunctionDerivativeOne
	//filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\shapeFunctionDerivativeOne.txt";

	//if((fpout = fopen(filename_pw, "a")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	fprintf(fpout,"%d\n",indexPointSet);
	//	fprintf(fpout,"/////////////////////////////////");
	//	for(int i=0;i<m_numOfNeighbors+1;i++){
	//		for(int j=0;j<m_numOfIntegralPoint;j++){				
	//			fprintf(fpout,"%f ",m_shapeFunctionDerivativeOne[i*m_numOfIntegralPoint+j]);				 
	//		}
	//		fprintf(fpout,"//////////////////////////");
	//		fprintf(fpout,"\n");            
	//	}
	//	fclose(fpout);
	//}

	////////////////////////////////////////////////////////////////////////////
	//// 写文件 m_shapeFunctionDerivativeTwo
	//filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\shapeFunctionDerivativeTwo.txt";

	//if((fpout = fopen(filename_pw, "a")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	fprintf(fpout,"%d\n",indexPointSet);
	//	fprintf(fpout,"/////////////////////////////////");
	//	for(int i=0;i<m_numOfNeighbors+1;i++){
	//		for(int j=0;j<m_numOfIntegralPoint;j++){				
	//			fprintf(fpout,"%f ",m_shapeFunctionDerivativeTwo[i*m_numOfIntegralPoint+j]);				 
	//		}
	//		fprintf(fpout,"//////////////////////////");
	//		fprintf(fpout,"\n");            
	//	}
	//	fclose(fpout);
	//}

	return;
}
void OneFreedomMeshless::ComputeManifoldMetric(int indexPointSet)
{
	/*
	*	函数功能： 计算流形标准的代数行列式
	*  变量说明：
	*  int indexPointSet        索引点
	*  double* m_determinentMetric   在采样点处流形标准的代数行列式
	*/
	//////////////////////////////////////////////////////////////////////////
	// 分配内存
	m_determinentMetric=new double[m_numOfIntegralPoint];
	m_metricMatrix=new double[m_numOfIntegralPoint*4];
	m_metricMatrixInve=new double[m_numOfIntegralPoint*4];
	m_localOrdinateVectorOne=new double[m_numOfIntegralPoint*3];
	m_localOrdinateVectorTwo=new double[m_numOfIntegralPoint*3];
	//////////////////////////////////////////////////////////////////////////
	// 循环计算
	for(int i=0;i<m_numOfIntegralPoint;i++){
		double tempVectorOne[3],tempVectorTwo[3];// 参考3.5
		for(int j=0;j<3;j++){
			tempVectorOne[j]=m_shapeFunctionDerivativeOne[i]*m_resultPointSet[indexPointSet*3+j];
			tempVectorTwo[j]=m_shapeFunctionDerivativeTwo[i]*m_resultPointSet[indexPointSet*3+j];			
		}
		std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(indexPointSet);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		for(int j=1;j<m_numOfNeighbors+1;j++){
			for(int k=0;k<3;k++){
				tempVectorOne[k]+=m_shapeFunctionDerivativeOne[j*m_numOfIntegralPoint+i]*m_resultPointSet[(*vectorIterator)*3+k];
				tempVectorTwo[k]+=m_shapeFunctionDerivativeTwo[j*m_numOfIntegralPoint+i]*m_resultPointSet[(*vectorIterator)*3+k];
			}
			vectorIterator++;
		}
		//////////////////////////////////////////////////////////////////////////
		// 计算空间坐标对局部坐标的导数
		for(int j=0;j<3;j++){
			m_localOrdinateVectorOne[i*3+j]=tempVectorOne[j];
			m_localOrdinateVectorTwo[i*3+j]=tempVectorTwo[j];
		}
		double localMetric[4];
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
	//CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\deteMetric.txt";
	//FILE *fpout;
	//if((fpout = fopen(filename_pw, "a")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	fprintf(fpout,"%d\n",indexPointSet);
	//	fprintf(fpout,"/////////////////////////////////");
	//	for(int i=0;i<m_numOfIntegralPoint;i++){

	//		fprintf(fpout,"%f ",m_determinentMetric[i]);				 


	//	}
	//	fprintf(fpout,"//////////////////////////");
	//	fprintf(fpout,"\n");      
	//	fclose(fpout);
	//}
	return;
}
void OneFreedomMeshless::ComputeStiffMatrix(int indexPointSet)
{
	//////////////////////////////////////////////////////////////////////////
	//首先计算 对角线上的元素
	if(indexPointSet==601){
		int k;
		k=0;
	}
	double namdaLoadConstant=m_namda*m_loadConstant;
	double elementOne=0;
	double elementTwo=0;
	double elementStiffMatrix;
	//////////////////////////////////////////////////////////////////////////
	// 求积分
	for(int i=0;i<m_numOfIntegralPoint;i++){
		if(m_determinentMetric[i]<0&abs(m_determinentMetric[i])>0.000001){
			int wwwwwww;
			wwwwwww=0;
		}

		elementOne+=m_shapeFunction[i]*m_testFunction[i]*sqrt(m_determinentMetric[i])*m_ordinateAndWeightCoff[i*3+2];
		//elementTwo+=(m_shapeFunctionDerivativeOne[i]*m_testFunctionDerivativeOne[i]+m_shapeFunctionDerivativeTwo[i]*m_testFunctionDerivativeTwo[i])*m_ordinateAndWeightCoff[i*3+2];
	/*	double tempElementTwo=0;
		for(int j=0;j<3;j++){
			tempElementTwo+=m_gradientShapeFunction[i*3+j]*m_gradientTestFunction[j];
		}*/
		
		double tempElementTwo=0;
		tempElementTwo+=m_shapeFunctionDerivativeOne[i]*m_testFunctionDerivativeOne[i]
				*(m_metricMatrixInve[i*4]*m_metricMatrixInve[i*4]+m_metricMatrixInve[i*4+2]*m_metricMatrixInve[i*4+2]);
		tempElementTwo+=m_shapeFunctionDerivativeOne[i]*m_testFunctionDerivativeTwo[i]
				*(m_metricMatrixInve[i*4]*m_metricMatrixInve[i*4+1]+m_metricMatrixInve[i*4+2]*m_metricMatrixInve[i*4+3]);
		tempElementTwo+=m_shapeFunctionDerivativeTwo[i]*m_testFunctionDerivativeOne[i]
				*(m_metricMatrixInve[i*4+2]*m_metricMatrixInve[i*4]+m_metricMatrixInve[i*4+3]*m_metricMatrixInve[i*4+2]);
		tempElementTwo+=m_shapeFunctionDerivativeTwo[i]*m_testFunctionDerivativeTwo[i]
				*(m_metricMatrixInve[i*4+2]*m_metricMatrixInve[i*4+2]+m_metricMatrixInve[i*4+3]*m_metricMatrixInve[i*4+3]);
		elementTwo+=tempElementTwo*m_ordinateAndWeightCoff[i*3+2];		
	}
	//////////////////////////////////////////////////////////////////////////
	// 刚度矩阵元素
	elementStiffMatrix=elementOne+m_namda*elementTwo;
	m_massMatrix[indexPointSet].m_elementStiffMarix.push_back(elementOne);
	m_stiffMatrix[indexPointSet].m_elementStiffMarix.push_back(m_namda*elementTwo);
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
		for(int j=0;j<m_numOfIntegralPoint;j++){
			elementOne+=m_shapeFunction[i*m_numOfIntegralPoint+j]*m_testFunction[j]*sqrt(m_determinentMetric[j])*m_ordinateAndWeightCoff[j*3+2];
		//	elementTwo+=(m_shapeFunctionDerivativeOne[i*m_numOfIntegralPoint+j]*m_testFunctionDerivativeOne[j]+m_shapeFunctionDerivativeTwo[i*m_numOfIntegralPoint+j]*m_testFunctionDerivativeTwo[j])*m_ordinateAndWeightCoff[j*3+2];
		/*	double tempElementTwo=0;
			for(int k=0;k<3;k++){
				tempElementTwo+=m_gradientShapeFunction[i*m_numOfIntegralPoint*3+j*3+k]*m_gradientTestFunction[j*3+k];
			}*/
		//	elementTwo+=m_shapeFunctionDerivativeOne[i*m_numOfIntegralPoint+j]*m_TestFunctionDerivativeOne[j]+m_shapeFunctionDerivativeTwo[i*m_numOfIntegralPoint+j]*m_TestFunctionDerivativeTwo[j];
			double tempElementTwo=0;
			tempElementTwo+=m_shapeFunctionDerivativeOne[i*m_numOfIntegralPoint+j]*m_testFunctionDerivativeOne[j]
					*(m_metricMatrixInve[j*4]*m_metricMatrixInve[j*4]+m_metricMatrixInve[j*4+2]*m_metricMatrixInve[j*4+2]);
			tempElementTwo+=m_shapeFunctionDerivativeOne[i*m_numOfIntegralPoint+j]*m_testFunctionDerivativeTwo[j]
					*(m_metricMatrixInve[j*4]*m_metricMatrixInve[j*4+1]+m_metricMatrixInve[j*4+2]*m_metricMatrixInve[j*4+3]);
			tempElementTwo+=m_shapeFunctionDerivativeTwo[i*m_numOfIntegralPoint+j]*m_testFunctionDerivativeOne[j]
					*(m_metricMatrixInve[j*4+2]*m_metricMatrixInve[j*4]+m_metricMatrixInve[j*4+3]*m_metricMatrixInve[j*4+2]);
			tempElementTwo+=m_shapeFunctionDerivativeTwo[i*m_numOfIntegralPoint+j]*m_testFunctionDerivativeTwo[j]
					*(m_metricMatrixInve[j*4+2]*m_metricMatrixInve[j*4+2]+m_metricMatrixInve[j*4+3]*m_metricMatrixInve[j*4+3]);
			elementTwo+=tempElementTwo*m_ordinateAndWeightCoff[j*3+2];
		}
		//elementStiffMatrix=elementOne+m_namda*elementTwo;
		elementStiffMatrix=m_namda*elementTwo;
		m_stiffMatrix[indexPointSet].m_elementStiffMarix.push_back(m_namda*elementTwo);
		/*for(int k=0;k<3;k++){
			m_loadVector[indexPointSet*3+k]+=namdaLoadConstant*elementOne*m_originalPointSet[3*(*vectorIterator)+k]
			+(1-namdaLoadConstant)*elementOne*m_resultPointSet[3*(*vectorIterator)+k];
		}*/
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
		std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[indexPointSet].m_elementStiffMarix.begin();
		for(;iteratorStiffMatrix!=m_stiffMatrix[indexPointSet].m_elementStiffMarix.end();iteratorStiffMatrix++){
			fprintf(fpout,"%f ",(*iteratorStiffMatrix));
		}		
		fprintf(fpout,"//////////////////////////");
		fprintf(fpout,"\n");      
		fclose(fpout);
	}
	return;
}
void OneFreedomMeshless::EstimateNormalDirection()
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

	std::multimap<double,int> valueKeySeries;//根据点的代价索引的map
	std::map<int,double> pointKeySeries;//根据点索引的map

	//初始化map
	for(int i=0;i<maxZpoint;i++){
		pointKeySeries.insert(std::map<int,double>::value_type(i,i+10));
		valueKeySeries.insert(std::map<double,int>::value_type(i+10,i));
	}
	for(int i=maxZpoint+1;i<m_numOfPoints;i++){
		pointKeySeries.insert(std::map<int,double>::value_type(i,i+10));
		valueKeySeries.insert(std::map<double,int>::value_type(i+10,i));
	}
	pointKeySeries.insert(std::map<int,double>::value_type(maxZpoint,0));
	valueKeySeries.insert(std::map<double,int>::value_type(0,maxZpoint));

	std::map<int,double>::iterator mapIterator = pointKeySeries.begin();
	std::multimap<double,int>::iterator multimapIterator=valueKeySeries.begin();

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
				double relativeCost,normalDotProduct;
				double pointExtractNormal[3],neighborPointNormal[3];
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
					valueKeySeries.insert(std::map<double,int>::value_type(relativeCost,neighborPoint));
				}
			}
		}
		mapIterator=pointKeySeries.begin();
		multimapIterator=valueKeySeries.begin();
	}

	/*if(mapIterator!=pointKeySeries.end()){
	int numOfPointsLeave;
	double* pointSetLeave;
	numOfPointsLeave=pointKeySeries.size();
	pointSetLeave=new double[numOfPointsLeave*3];
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
void OneFreedomMeshless::CalculateLinearSystem()
{

	//////////////////////////////////////////////////////////////////////////
	// 循环开始
	double threshold=m_ebusainu+1;
	while(threshold>=m_ebusainu){
		threshold=0.0;
		for(int i=0;i<m_numOfPoints;i++){
			std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
			std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
			std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
			iteratorStiffMatrix++;
			double tempSum[3];
			for(int j=0;j<3;j++){
				tempSum[j]=0;
			}
			for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,iteratorStiffMatrix++){
				for(int j=0;j<3;j++){
					tempSum[j]+=(*iteratorStiffMatrix)*m_resultPointSet[(*vectorIterator)*3+j];
				}				
			}
			iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
			double localPoint[3];
			for(int j=0;j<3;j++){
				localPoint[j]=(m_loadVector[i*3+j]-tempSum[j])/(*iteratorStiffMatrix);
				if(localPoint[j]<10000||localPoint[j]>-10000)
					int kkd;
				else{
					int www;
					www=0;
				}
				
			}
			for(int j=0;j<3;j++){
				double localThreshold;
				localThreshold=abs(localPoint[j]-m_resultPointSet[i*3+j]);///(1.0+abs(localPoint[j]));
				if(localThreshold>threshold)
					threshold=localThreshold;
				m_resultPointSet[i*3+j]=localPoint[j];
			}            
		}
		//////////////////////////////////////////////////////////////////////////
		// 调试代码
	threshold=0;
	}
	return;
}
void OneFreedomMeshless::ComputeGradientShapeFunction()
{
	m_gradientShapeFunction=new double[(m_numOfNeighbors+1)*m_numOfIntegralPoint*3];// 存放任一个形函数在所有的积分点的流形梯度
	m_gradientTestFunction=new double[m_numOfIntegralPoint*3];
	for(int i=0;i<m_numOfNeighbors+1;i++){
		for(int j=0;j<m_numOfIntegralPoint;j++){
			 
			double tempSumOne,tempSumTwo;
			tempSumOne=m_metricMatrixInve[j*4]*m_shapeFunctionDerivativeOne[i*m_numOfIntegralPoint+j]
					+m_metricMatrixInve[j*4+1]*m_shapeFunctionDerivativeTwo[i*m_numOfIntegralPoint+j];
			tempSumTwo=m_metricMatrixInve[j*4+2]*m_shapeFunctionDerivativeOne[i*m_numOfIntegralPoint+j]
					+m_metricMatrixInve[j*4+3]*m_shapeFunctionDerivativeTwo[i*m_numOfIntegralPoint+j];
			for(int k=0;k<3;k++){
				m_gradientShapeFunction[i*m_numOfIntegralPoint*3+j*3+k]=tempSumOne*m_localOrdinateVectorOne[j*3+k]
						+tempSumTwo*m_localOrdinateVectorTwo[j*3+k];
			}
			
		}
	}
	for(int i=0;i<m_numOfIntegralPoint;i++){
		double tempSumOne,tempSumTwo;
		tempSumOne=m_metricMatrixInve[i*4]*m_testFunctionDerivativeOne[i]+m_metricMatrixInve[i*4+1]*m_testFunctionDerivativeTwo[i];
		tempSumTwo=m_metricMatrixInve[i*4+2]*m_testFunctionDerivativeOne[i]+m_metricMatrixInve[i*4+3]*m_testFunctionDerivativeTwo[i];
		for(int j=0;j<3;j++){
			m_gradientTestFunction[i*3+j]=tempSumOne*m_localOrdinateVectorOne[i*3+j]
					+tempSumTwo*m_localOrdinateVectorTwo[i*3+j];
		}
	}
	return;
	
}

void OneFreedomMeshless::CalculateLinearSystemGauss()
{
	int* js;
	int i,j,k,is,u,v;
	double d,t;
	double tempMajorElement;
	js=new int[m_numOfPoints];
	struct Knearest{
		std::vector<int> m_nearest;//存放邻域点
	};// 保存点的邻域的结构
	Knearest* tempKnearest;///包括自身点
	tempKnearest=new Knearest[m_numOfPoints];
	//////////////////////////////////////////////////////////////////////////
	// 初始化 邻域结构
	for(i=0;i<m_numOfPoints;i++){
		tempKnearest[i].m_nearest.push_back(i);
		std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
            tempKnearest[i].m_nearest.push_back((*vectorIterator));
		}
	}
	for(k=0;k<m_numOfPoints;k++){// 对角化矩阵
		d=0.0;
		tempMajorElement=0.0;
		for(i=k;i<m_numOfPoints;i++){//选主元
			std::vector<int>::iterator vectorKnearest=tempKnearest[i].m_nearest.begin();
			std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
			for(;vectorKnearest!=tempKnearest[i].m_nearest.end();vectorKnearest++,iteratorStiffMatrix++){
				if((vectorKnearest==tempKnearest[i].m_nearest.end())&(iteratorStiffMatrix!=m_stiffMatrix[i].m_elementStiffMarix.end())){
					int error;
					error=0;
				}
				if((vectorKnearest!=tempKnearest[i].m_nearest.end())&(iteratorStiffMatrix==m_stiffMatrix[i].m_elementStiffMarix.end())){
					int error;
					error=0;
				}
				if((*vectorKnearest)>=k){
					t=abs((*iteratorStiffMatrix));
					if(t>d){
						d=t;
						js[k]=(*vectorKnearest);
						is=i;
						tempMajorElement=(*iteratorStiffMatrix);
					}
				}
			}

		}
		if(is!=k){///交换行
			std::vector<double>::iterator iteratorStiffMatrixOne=m_stiffMatrix[k].m_elementStiffMarix.begin();
			std::vector<double>::iterator iteratorStiffMatrixTwo=m_stiffMatrix[is].m_elementStiffMarix.begin();
			for(;iteratorStiffMatrixTwo!=m_stiffMatrix[is].m_elementStiffMarix.end();iteratorStiffMatrixTwo++,iteratorStiffMatrixOne++){
				//////////////////////////////////////////////////////////////////////////
				// 交换刚度矩阵
				if((iteratorStiffMatrixTwo!=m_stiffMatrix[is].m_elementStiffMarix.end())&(iteratorStiffMatrixOne==m_stiffMatrix[k].m_elementStiffMarix.end())){
					int error;
					error=0;
				}
				if((iteratorStiffMatrixTwo==m_stiffMatrix[is].m_elementStiffMarix.end())&(iteratorStiffMatrixOne!=m_stiffMatrix[k].m_elementStiffMarix.end())){
					int error;
					error=0;
				}
				double tempU,tempV;
				tempU=(*iteratorStiffMatrixOne);
				tempV=(*iteratorStiffMatrixTwo);
				(*iteratorStiffMatrixOne)=tempV;
				(*iteratorStiffMatrixTwo)=tempU;
			}
			std::vector<int>::iterator vectorKnearestOne=tempKnearest[k].m_nearest.begin();
			std::vector<int>::iterator vectorKnearestTwo=tempKnearest[is].m_nearest.begin();
			
			
			for(;vectorKnearestOne!=tempKnearest[k].m_nearest.end();vectorKnearestOne++,vectorKnearestTwo++){
				//////////////////////////////////////////////////////////////////////////
				// 交换 邻域
				if((vectorKnearestOne!=tempKnearest[k].m_nearest.end())&(vectorKnearestTwo==tempKnearest[is].m_nearest.end())){
					int error;
					error=0;
				}
				if((vectorKnearestOne==tempKnearest[k].m_nearest.end())&(vectorKnearestTwo!=tempKnearest[is].m_nearest.end())){
					int error;
					error=0;
				}
				int tempU,tempV;
				tempU=(*vectorKnearestOne);
				tempV=(*vectorKnearestTwo);
				(*vectorKnearestOne)=tempV;
				(*vectorKnearestTwo)=tempU;
			}
			for(int www=0;www<3;www++){
				//////////////////////////////////////////////////////////////////////////
				// 交换常数向量，方程右边的东西
				double tempScalar;
				tempScalar=m_loadVector[is*3+www];
				m_loadVector[is*3+www]=m_loadVector[k*3+www];
				m_loadVector[k+3+www]=tempScalar;
			}
		}

		if(js[k]!=k){//交换列 有错误  在于交换的时候出问题了，冲突了
			for(int www=0;www<m_numOfPoints;www++){
				std::vector<int>::iterator vectorKnearest=tempKnearest[www].m_nearest.begin();
				for(;vectorKnearest!=tempKnearest[www].m_nearest.end();vectorKnearest++){
					if((*vectorKnearest)==js[k]){
						(*vectorKnearest)=k;
						continue;
					}
					if((*vectorKnearest)==k)
						(*vectorKnearest)=js[k];
				}                
			}            
		}
		//////////////////////////////////////////////////////////////////////////
		// 系数矩阵归一化，常数向量归一化
		//std::vector<int>::iterator tempVectorKnearest=tempKnearest[k].m_nearest.begin();
		std::vector<double>::iterator tempiteratorStiffMatrix=m_stiffMatrix[k].m_elementStiffMarix.begin();
	
		for(;tempiteratorStiffMatrix!=m_stiffMatrix[k].m_elementStiffMarix.end();tempiteratorStiffMatrix++){
			if(tempMajorElement==0){
				int error;
				error=0;
			}
			if((*tempiteratorStiffMatrix)!=0.0)
				(*tempiteratorStiffMatrix)/=tempMajorElement;
		}
		for(int www=0;www<3;www++){
			m_loadVector[k*3+www]/=tempMajorElement;
		}
		//////////////////////////////////////////////////////////////////////////
		// 系数矩阵消元
		
	
		for(int www=k+1;www<m_numOfPoints;www++){
			std::vector<int>::iterator iteratorKnearestOne=tempKnearest[www].m_nearest.begin();
			std::vector<double>::iterator iteratorStiffMatrixOne=m_stiffMatrix[www].m_elementStiffMarix.begin();
			double tempScalarUnderMajor;//存放主元正下方的元素
			bool ifIncludeK=false;//判断正下方元素是否为0
			for(;iteratorKnearestOne!=tempKnearest[www].m_nearest.end();iteratorKnearestOne++,iteratorStiffMatrixOne++){
				if(((*iteratorKnearestOne)==k)&((*iteratorStiffMatrixOne)!=0.0)){
					ifIncludeK=true;
					tempScalarUnderMajor=(*iteratorStiffMatrixOne);
					(*iteratorStiffMatrixOne)=0;
					break;
				}
			}
			if(ifIncludeK==true){
				//////////////////////////////////////////////////////////////////////////
				/// 对刚度矩阵消元
				std::vector<double>::iterator tempiteratorStiffMatrixTwo=m_stiffMatrix[k].m_elementStiffMarix.begin();
				std::vector<int>::iterator tempIteratorKnearestTwo=tempKnearest[k].m_nearest.begin();
				for(;tempIteratorKnearestTwo!=tempKnearest[k].m_nearest.end();tempIteratorKnearestTwo++,tempiteratorStiffMatrixTwo++){
					if(((*tempiteratorStiffMatrixTwo)!=0.0)&((*tempIteratorKnearestTwo)>k)){
						//////////////////////////////////////////////////////////////////////////
						/// 对 www 行消元
						iteratorKnearestOne=tempKnearest[www].m_nearest.begin();
						iteratorStiffMatrixOne=m_stiffMatrix[www].m_elementStiffMarix.begin();
						for(;iteratorKnearestOne!=tempKnearest[www].m_nearest.end();iteratorKnearestOne++,iteratorStiffMatrixOne++){
							if((*iteratorKnearestOne)==(*tempIteratorKnearestTwo)){// 查找www行与(*tempIteratorKnearestTwo)同列的元素
								(*iteratorStiffMatrixOne)-=tempScalarUnderMajor*(*tempiteratorStiffMatrixTwo);							
							}
						}

					}
				}
				//////////////////////////////////////////////////////////////////////////
				// 对常数量消元
				for(int wwww=0;wwww<3;wwww++){
					m_loadVector[www*3+wwww]-=tempScalarUnderMajor*m_loadVector[k*3+wwww];
				}

			}
			else
				continue;
		}
	}
	//////////////////////////////////////////////////////////////////////////
	// 解方程
	double* tempResult;
	tempResult=new double[3*m_numOfPoints];
	//////////////////////////////////////////////////////////////////////////
	// 初始化解
	for(i=0;i<m_numOfPoints*3;i++){
		tempResult[i]=0.0;
	}
	//////////////////////////////////////////////////////////////////////////
	// 循环求解
	for(i=(m_numOfPoints-1);i>=0;i--){
		std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
		std::vector<int>::iterator iteratorKnearest=tempKnearest[i].m_nearest.begin();
		for(;iteratorKnearest!=tempKnearest[i].m_nearest.end();iteratorKnearest++,iteratorStiffMatrix++){
			if((*iteratorKnearest)!=i){
				for(j=0;j<3;j++){
					m_loadVector[i*3+j]-=(*iteratorStiffMatrix)*tempResult[i*3+j];
				}
			}
			if((*iteratorKnearest)=i)
				tempMajorElement=(*iteratorStiffMatrix);			
		}
		for(j=0;j<3;j++){
			tempResult[i*3+j]=m_loadVector[i*3+j]/tempMajorElement;
		}
	}

	/// 交换解的次序
	for(i=(m_numOfPoints-1);i>=0;i--){
		if(i!=js[i]){
			for(j=0;j<3;j++){
				t=tempResult[i*3+j];
				tempResult[i*3+j]=tempResult[js[i]*3+j];
				tempResult[js[i]*3+j]=t;
				
			}
		}
	}
	for(i=0;i<m_numOfPoints*3;i++){
		m_resultPointSet[i]=tempResult[i];
	}
	//////////////////////////////////////////////////////////////////////////
	/// 释放内存
	delete[] tempKnearest;
	delete[] js;
	delete[] tempResult;
	return;
}
void OneFreedomMeshless::AssembleStiffMatrix()
{
	//////////////////////////////////////////////////////////////////////////	
	/*
	*	函数功能： 组装刚度矩阵
	*/
	//////////////////////////////////////////////////////////////////////////
	double namdaLoadConstant=m_namda*m_loadConstant;
	for (int i=0;i<m_numOfPoints;i++){
		std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
		std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		double sumOfElement=0;
		//	iteratorMassMatrix++;
		//vectorIterator++;
		for(iteratorStiffMatrix++;iteratorStiffMatrix!=m_stiffMatrix[i].m_elementStiffMarix.end();iteratorStiffMatrix++,vectorIterator++){
			std::map<int, KnearestField>::iterator mapIteratorIner=m_mapKnearest.find((*vectorIterator));
			std::vector<int>::iterator vectorIteratorIner=(*mapIteratorIner).second.m_nearest.begin();
			std::vector<double>::iterator iteratorStiffMatrixIner=m_stiffMatrix[(*vectorIterator)].m_elementStiffMarix.begin();
			bool temp=false;	
		//	vectorIteratorIner++;
			for(iteratorStiffMatrixIner++;iteratorStiffMatrixIner!=m_stiffMatrix[(*vectorIterator)].m_elementStiffMarix.end();iteratorStiffMatrixIner++,vectorIteratorIner++){
				if((*vectorIteratorIner)==i){
					(*iteratorStiffMatrix)=((*iteratorStiffMatrix)+(*iteratorStiffMatrixIner))/2;
					(*iteratorStiffMatrixIner)=(*iteratorStiffMatrix);
					temp=true;
					break;
				}
			}
			if(temp==false){
				(*iteratorStiffMatrix)/=2;
				m_stiffMatrix[(*vectorIterator)].m_elementStiffMarix.push_back((*iteratorStiffMatrix));
				(*mapIteratorIner).second.m_nearest.push_back(i);
			}
			//	sumOfElement+=(*iteratorStiffMatrix);
			//		(*iteratorStiffMatrix)=(*iteratorMassMatrix)+(*iteratorStiffMatrix);
		}
		//iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
		//	iteratorMassMatrix=m_massMatrix[i].m_elementStiffMarix.begin();
		//(*iteratorStiffMatrix)=(*iteratorMassMatrix)-sumOfElement;
	}

	for(int i=0;i<m_numOfPoints;i++){
		std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
		std::vector<double>::iterator iteratorMassMatrix=m_massMatrix[i].m_elementStiffMarix.begin();
		double sumOfElement=0;
		iteratorMassMatrix++;
		for(iteratorStiffMatrix++;iteratorStiffMatrix!=m_stiffMatrix[i].m_elementStiffMarix.end();iteratorStiffMatrix++,iteratorMassMatrix++){
			sumOfElement+=(*iteratorStiffMatrix);
			//(*iteratorStiffMatrix)=(*iteratorMassMatrix)
		}
		iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
		iteratorMassMatrix=m_massMatrix[i].m_elementStiffMarix.begin();
		for(int j=0;j<3;j++){
			m_loadVector[i*3+j]=namdaLoadConstant*(*iteratorMassMatrix)*m_originalPointSet[i*3+j]
			+(1-namdaLoadConstant)*(*iteratorMassMatrix)*m_resultPointSet[i*3+j];

		}
		(*iteratorStiffMatrix)=(*iteratorMassMatrix)-sumOfElement;
	}
	return;
}

