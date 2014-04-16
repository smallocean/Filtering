#include "stdafx.h"
#include ".\rbfmeshless.h"
#include ".\mathlib\mathlib.h"
#include "Mcube.h"
#include "Matrix.h"
#include <math.h>
//#include "engine.h"
using namespace MATHLIB;
const int RBFmeshless::m_numOfIntegralPoint=40;
const double RBFmeshless::m_ordinateAndWeightCoff[120]={ 
	//////////////////////////////////////////////////////////////////////////
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
const double RBFmeshless::m_coffRBFunction[5]={8,40,48,25,5};

RBFmeshless::RBFmeshless(void)
{
	// 初始化所有数据
	m_originalPointSet=NULL;
	m_originalColors=NULL;
	m_numOfPoints=0;
	m_kNearest=0;
	m_originalNormals=NULL;
	m_radius=0;//计算 m_mapKnearest时候的搜索半径
	m_interactiveTime=0; // 滤波的次数
	m_loadVector=NULL;//荷载向量
    m_namda=0;// 时间步长
	//	double m_loadConstant;//荷载常数
    m_ebusainu=0;//循环停止的条件；
	m_numOfNeighbors=0;          //邻域数不包括中心点
    m_localRBFradius=0;//         局部径向基函数半径
	m_localOrdinate=NULL;//             各邻域点的局部坐标
	m_matrixNeighbor=NULL;//            插值矩阵
	m_matrixNeighborInv=NULL;//         插值矩阵的逆阵
	m_localIntegralPoints=NULL;  //积分点的坐标
	m_localIntergralRadius=0;//      局部积分半径
	m_testFunction=NULL;//       测试函数在每一个采样点的值
	m_RBFfunction=NULL;//        每一行为在每一个邻域点的径向基函数在每一个积分点的函数值
    m_shapeFunction=NULL;//                各邻域点的形函数在各采样点的值
	m_determinentMetric=NULL;//   在采样点处流形标准的代数行列式
	m_metricMatrix=NULL;//在采样点处流形标准的矩阵
	m_metricMatrixInve=NULL;////在采样点处流形标准的矩阵的逆
	m_resultPointSet=NULL;
	m_resultColor=NULL;
	nearestTime=0;
	filterTime=0;
	m_RBFfunctionDerOne=NULL;
	m_RBFfunctionDerTwo=NULL;
	m_shapeFunctionDerOne=NULL;
	m_shapeFunctionDerTwo=NULL;
	m_testFunctionDerOne=NULL;
	m_testFunctionDerTwo=NULL;
	m_stiffMatrix=NULL;
	m_massMatrix=NULL;
	m_loadVector=NULL;
	m_originalNormals=NULL;
	m_originalMajorDirection=NULL;
	m_originalMinorDirection=NULL;
	m_meanCurvature=NULL;
	m_resultColor=NULL;
}

RBFmeshless::~RBFmeshless(void)
{
	DeleteMeshlessFilter();
	
}

void RBFmeshless::GetMeshlessFilter(int numOfPointSet,float* pointSet,int kNearest, double radius, int interactiveTime,double namda,	double loadConstant,double ebusainu, int stopN,	int maxIter)
{
	/*
	函数功能：  接口函数，计算点集滤波
	变量说明：
	int numOfPointSet                        点的数量
	double* pointSet							原始点集
	int kNearest                             k邻域的数量
	double radius                             k邻域搜索的半径
	int interactiveTime                      滤波的交互次数
	double namda                             时间步长
	double loadConstant                      为保护体积加的荷载
	double ebusainu                          计算线性方程组时候的停止标准
	*/


	m_numOfPoints=numOfPointSet;
	m_kNearest=kNearest;
	m_radius=radius;
	//m_IfFilterNormal=ifFilterNormal;
	m_interactiveTime=interactiveTime;
	m_namda=namda;
	m_loadConstant=loadConstant;
	m_ebusainu=ebusainu;
	m_itol=stopN;
	m_itmax=maxIter;
	m_tol=ebusainu;
	if(m_originalPointSet!=NULL){
		delete[] m_originalPointSet;
		m_originalPointSet=NULL;
	}
	m_originalPointSet=new double[m_numOfPoints*3];
	if(m_resultPointSet!=NULL){
		delete[] m_resultPointSet;
		m_resultPointSet=NULL;
	}
	m_resultPointSet=new double[m_numOfPoints*3];

	for(int i=0;i<m_numOfPoints*3;i++){
		m_originalPointSet[i]=pointSet[i];
		m_resultPointSet[i]=pointSet[i];
	}

	//////////////////////////////////////////////////////////////////////////
	// 求k 近邻
	//ComputeMapKnearest();
	//////////////////////////////////////////////////////////////////////////
	// 求法向 和两个主方向
	//ComputeNormal();
	//////////////////////////////////////////////////////////////////////////
	//循环求解
	for(int i=0;i<m_interactiveTime;i++){
		//////////////////////////////////////////////////////////////////////////
		// 建立刚度矩阵
		ComputeMapKnearest();
		ComputeNormal();
		if(m_stiffMatrix!=NULL){
			delete[] m_stiffMatrix;
			m_stiffMatrix=NULL;
		}
		m_stiffMatrix=new StiffMatrix[m_numOfPoints];
		if(m_massMatrix!=NULL){
			delete[] m_massMatrix;
			m_massMatrix=NULL;
		}
		m_massMatrix=new StiffMatrix[m_numOfPoints];
		if(m_loadVector!=NULL){
			delete[] m_loadVector;
			m_loadVector=NULL;
		}
		m_loadVector=new double[m_numOfPoints];
		for(int j=0;j<m_numOfPoints;j++){
			
			ComputeNearestParameter(j);
			ComputeTestFunction(j);
			ComputeMatrixNeighbor(j);
			ComputeMatrixRelatedSampling(j);
			ComputeShapeFunction(j);
			ComputeManifoldMetric(j);
			//	ComputeGradientShapeFunction();

			ComputeStiffMatrix(j);

			if(m_localOrdinate!=NULL){
				delete[] m_localOrdinate;
                m_localOrdinate=NULL;
			}
			
			if(m_shapeFunction!=NULL){
				delete[] m_shapeFunction;
                m_shapeFunction=NULL; //一行为一个邻域点的形函数在不同的积分点的值
			}
			if(m_shapeFunctionDerOne!=NULL){
                delete[] m_shapeFunctionDerOne;
			m_shapeFunctionDerOne=NULL;
			if(m_shapeFunctionDerTwo!=NULL)
				delete[] m_shapeFunctionDerTwo;
			m_shapeFunctionDerTwo=NULL;
			if(m_determinentMetric!=NULL)
				delete[] m_determinentMetric;
			m_determinentMetric=NULL;
			if(m_metricMatrix!=NULL)
				delete[] m_metricMatrix;
			m_metricMatrix=NULL;
			if(m_metricMatrixInve!=NULL)
				delete[] m_metricMatrixInve;
			m_metricMatrixInve=NULL;

			//if(m_localOrdinateVectorOne!=NULL)
			//	delete[] m_localOrdinateVectorOne;
			//m_localOrdinateVectorOne=NULL;
			//if(m_localOrdinateVectorTwo!=NULL)
			//	delete[] m_localOrdinateVectorTwo;
			//m_localOrdinateVectorTwo=NULL;
			/*if(m_gradientShapeFunction!=NULL)
			delete[] m_gradientShapeFunction;
			m_gradientShapeFunction=NULL;
			if(m_gradientTestFunction!=NULL)
			delete[] m_gradientTestFunction;
			m_gradientTestFunction=NULL;*/
			if(m_testFunction!=NULL)
				delete[] m_testFunction;
			m_testFunction=NULL;
			if(m_testFunctionDerOne!=NULL)
				delete[] m_testFunctionDerOne;
			m_testFunctionDerOne=NULL;
			if(m_testFunctionDerTwo!=NULL)
				delete[] m_testFunctionDerTwo;
            	m_testFunctionDerTwo=NULL;
			}
			if(m_localIntegralPoints!=NULL){
				delete[] m_localIntegralPoints;
                m_localIntegralPoints=NULL;	
			}
			if(m_matrixNeighbor!=NULL){
				delete[] m_matrixNeighbor;
				m_matrixNeighbor=NULL;
			}
			if(m_matrixNeighborInv!=NULL){
				delete[] m_matrixNeighborInv;
				m_matrixNeighborInv=NULL;
			}
			if(m_RBFfunction!=NULL){
				delete[] m_RBFfunction;
				m_RBFfunction=NULL;
			}
			if(m_RBFfunctionDerOne!=NULL){
				delete[] m_RBFfunctionDerOne;
				m_RBFfunctionDerOne=NULL;
			}
			if(m_RBFfunctionDerTwo!=NULL){
				delete[] m_RBFfunctionDerTwo;
				m_RBFfunctionDerTwo=NULL;
			}
			

			//////////////////////////////////////////////////////////////////////////
			/*
			*	释放内存
			*/

			
		}
		//////////////////////////////////////////////////////////////////////////
		// 解线性系统
		AssembleStiffMatrix();
		AdjustStiffMatrix();
		
			double* b;
			double* x;
			b=new double[m_numOfPoints];
			x=new double[m_numOfPoints];
			for(int k=0;k<m_numOfPoints;k++){
				b[k]=m_loadVector[k];
				x[k]=0;
			}
			int iter;
			double err;
			linbcg(b,x,m_itol,m_tol,m_itmax,iter,  err);
			for(int k=0;k<m_numOfPoints;k++){
				m_resultPointSet[k*3]+=x[k]*m_originalNormals[k*3];
				m_resultPointSet[k*3+1]+=x[k]*m_originalNormals[k*3+1];
				m_resultPointSet[k*3+2]+=x[k]*m_originalNormals[k*3+2];

			}
			delete[] b;
			delete[] x;

	

	//	CalculateLinearSystem();
		//	CalculateLinearSystemGauss();

		//写文件for matlab
		CString filename_pw = "D:\\matlabMatrix.txt";
		FILE *fpout;
		if((fpout = fopen(filename_pw, "w")) == NULL)
		{
			int dkjkd;
			//MessageBox("can't open the file!");
		}
		else
		{
			/*fprintf(fpout,"%d\n",indexPointSet);
			fprintf(fpout,"/////////////////////////////////");*/
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
		filename_pw = "D:\\matlabLoadVector.txt";		
		if((fpout = fopen(filename_pw, "w")) == NULL)
		{
			int dkjkd;
			//MessageBox("can't open the file!");
		}
		else
		{
			for(int i=0;i<m_numOfPoints;i++){
				//for(int j=0;j<3;j++){
					fprintf(fpout,"%f\n",m_loadVector[i]);

				//}	

			}		
			fclose(fpout);
		}

		//////////////////////////////////////////////////////////////////////////
		/*
		 *	matlab 代码
		 */
	 /*   Engine *ep;
		char cmd[256];
		ep=engOpen(NULL);
		engEvalString(ep,"A=zeros(1,10);");
		engEvalString(ep,"load( 'C:\Program Files\MATLAB71\work\matlabMatrix.txt')");
		engEvalString(ep,"A = spconvert(matlabMatrix);");
		engEvalString(ep,"load( 'C:\Program Files\MATLAB71\work\matlabLoadvector.txt')");
		engEvalString(ep,"B=spconvert(matlabLoadvector);");
		engEvalString(ep,"C=A\B;");
		engEvalString(ep,"fid=fopen('C:\Program Files\MATLAB71\work\result.txt','w');");
		for(int k=0;k<m_numOfPoints;k++){
			sprintf(cmd,"fprintf(fid,'%f %f %f\n', C(%d,1),C(%d,2),C(%d,3));",k+1,k+1,k+1);
			engEvalString(ep,cmd);
		}	
		engEvalString(ep,"fclose(fid);");
		engClose(ep);*/
		//////////////////////////////////////////////////////////////////////////
		/*
		 *	读取文件
		 */
		//filename_pw="C:\\Program Files\\MATLAB71\\work\\result.txt";
		//for(int j=0;j<m_numOfPoints*3;j++){
		//	m_resultPointSet[j]=0;
		//}
	
		//if((fpout = fopen(filename_pw, "r")) == NULL)
		//{
		//	int dkjkd;
		//	//MessageBox("can't open the file!");
		//}
		//else
		//{
		//	std::vector<Point> tempPointSets;
		//	while(!feof(fpout)){
		//		Point tempPoint;
		//		float temp1,temp2;
		//		//fscanf(fpout,"%f %f %f %f %f",&tempPoint.vertexPosition[0],&tempPoint.vertexPosition[1],&tempPoint.vertexPosition[2],&temp1,&temp2);
		//		fscanf(fpout,"%f %f %f",&tempPoint.vertexPosition[0],&tempPoint.vertexPosition[1],&tempPoint.vertexPosition[2]);
		//		tempPointSets.push_back(tempPoint);
		//	}
		//	fclose(fpout);
		//	std::vector<Point>::iterator tempPointIterator=tempPointSets.begin();
		//	for(int i=0;tempPointIterator!=tempPointSets.end()&&i<m_numOfPoints;tempPointIterator++,i++){


		//		m_resultPointSet[i*3]=(*tempPointIterator).vertexPosition[0];//*1000;
		//		m_resultPointSet[i*3+1]=(*tempPointIterator).vertexPosition[1];//*1000;
		//		m_resultPointSet[i*3+2]=(*tempPointIterator).vertexPosition[2];//*1000;
		//		/*	m_colors[i*3]=(*tempColorSetsIterator).vertexPosition[0];
		//		m_colors[i*3+1]=(*tempColorSetsIterator).vertexPosition[1];
		//		m_colors[i*3+2]=(*tempColorSetsIterator).vertexPosition[2];*/


		//	}
		//}
		



		//////////////////////////////////////////////////////////////////////////
		// 释放内存
		for(int j=0;j<m_numOfPoints;j++){
			m_stiffMatrix[j].m_elementStiffMarix.clear();
			m_massMatrix[j].m_elementStiffMarix.clear();
		}
		if(m_originalNormals!=NULL){
            delete[] m_originalNormals;
			m_originalNormals=NULL;
		}
		if(m_originalMajorDirection!=NULL){
			delete[] m_originalMajorDirection;
			m_originalMajorDirection=NULL;
		}
		if(m_originalMinorDirection!=NULL){
			delete[] m_originalMinorDirection;
			m_originalMinorDirection=NULL;
		}
		m_mapKnearest.clear();

	}
	//////////////////////////////////////////////////////////////////////////
	// 写文件 loadVector FOR matlab
	
	

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


void RBFmeshless::DeleteMeshlessFilter()
{
	/*
	*	函数功能： 释放变量内存
	*/
	if(m_originalPointSet!=NULL)
		delete[] m_originalPointSet;
	m_originalPointSet=NULL;
	if(m_originalColors!=NULL)
		delete[] m_originalColors;
	m_originalColors=NULL;
	if(m_originalNormals!=NULL)
		delete[] m_originalNormals;
	m_originalNormals=NULL;
	if(m_loadVector!=NULL)
		delete[] m_loadVector;
	m_loadVector=NULL;//荷载向量
	if(m_localOrdinate!=NULL)
		delete[] m_loadVector;
	m_localOrdinate=NULL;//             各邻域点的局部坐标
	if(m_matrixNeighbor!=NULL)
		delete[] m_matrixNeighbor;
	m_matrixNeighbor=NULL;//            插值矩阵
	if(m_matrixNeighbor!=NULL)
		delete[] m_matrixNeighbor;
	m_matrixNeighborInv=NULL;//        插值矩阵的逆阵
	if(m_localIntegralPoints!=NULL)
		delete[] m_localIntegralPoints;
	m_localIntegralPoints=NULL;  //积分点的坐标
	if(m_testFunction!=NULL)
		delete[] m_testFunction;
	m_testFunction=NULL;//       测试函数在每一个采样点的值
	if(m_RBFfunction!=NULL)
		delete[] m_RBFfunction;
	m_RBFfunction=NULL;//       每一行为在每一个邻域点的径向基函数在每一个积分点的函数值
	if(m_shapeFunction!=NULL)
		delete[] m_shapeFunction;
	m_shapeFunction=NULL;//                各邻域点的形函数在各采样点的值
	if(m_determinentMetric!=NULL)
		delete[] m_determinentMetric;
	m_determinentMetric=NULL;//   在采样点处流形标准的代数行列式
	if(m_metricMatrix!=NULL)
		delete[] m_metricMatrix;
	m_metricMatrix=NULL;//在采样点处流形标准的矩阵
	if(m_metricMatrixInve!=NULL)
		delete[] m_metricMatrix;
	m_metricMatrixInve=NULL;////在采样点处流形标准的矩阵的逆
	if(m_resultPointSet!=NULL)
		delete[] m_resultPointSet;
	m_resultPointSet=NULL;
	if(m_resultColor!=NULL)
		delete[] m_resultColor;
	m_resultColor=NULL;
}
void RBFmeshless::ComputeMapKnearest()
{
	/*
	*	函数功能：计算每一个点的 k邻域
	*/
	double radius=m_radius*m_radius;//存放近邻点到所求点的距离
	std::multimap<double ,int> mapKnearest;//对每一个点建立邻域vector时建立临时map
	int numOfKnearest;//统计一个点的临时邻域大小
	if(m_mapKnearest.begin()!=m_mapKnearest.end()){
		m_mapKnearest.clear();
	}

	//////////////////////////////////////////////////////////////////////////
	//遍历所有点，并且按照距离由小到大存放存放到map中
	for(int i=0;i<m_numOfPoints;i++){
		numOfKnearest=0;
		for(int j=0;j<m_numOfPoints;j++){
			if(j==i)
				continue;
			double distancePointToPoint;
			distancePointToPoint=(m_resultPointSet[i*3]-m_resultPointSet[j*3])*(m_resultPointSet[i*3]-m_resultPointSet[j*3])
				+(m_resultPointSet[i*3+1]-m_resultPointSet[j*3+1])*(m_resultPointSet[i*3+1]-m_resultPointSet[j*3+1])
				+(m_resultPointSet[i*3+2]-m_resultPointSet[j*3+2])*(m_resultPointSet[i*3+2]-m_resultPointSet[j*3+2]);
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
	return;
}
void RBFmeshless::ComputeNormal()
{
	/*
	*	函数功能：计算每一个点的法向
	*/
	double centroidPosition[3];//重心位置
	double* localVariationVector;//重心到点的位置的向量
	int numOfNearest;

	mat_f8 mat_covaMatrix(3, 3);
	if(m_originalNormals!=NULL)
		delete[] m_originalNormals;
	if(m_originalNormals!=NULL){
		delete[] m_originalNormals;
		m_originalNormals=NULL;
	}
	if(m_originalMajorDirection!=NULL){
		delete[] m_originalMajorDirection;
		m_originalMajorDirection=NULL;
	}
	if(m_originalMinorDirection!=NULL){
		delete[] m_originalMinorDirection;
		m_originalMinorDirection=NULL;
	}
	m_originalNormals=new double[m_numOfPoints*3];
	m_originalMajorDirection=new double[m_numOfPoints*3];
	m_originalMinorDirection=new double[m_numOfPoints*3];
	//m_gaussCurvature=new double[m_numOfPoints*2];
	//if(m_meanAreas!=NULL)
	//	delete[] m_meanAreas;
	//m_meanAreas=new double[m_numOfPoints];
	if(m_meanCurvature!=NULL)
		delete[] m_meanCurvature;
	m_meanCurvature=new double[m_numOfPoints];
	if(m_resultColor!=NULL)
		delete[] m_resultColor;
	m_resultColor=new double[m_numOfPoints*3];


	for(int i=0;i<m_numOfPoints;i++){
		std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		numOfNearest=(*mapIterator).second.m_numOfNearest;
		for(int j=0;j<3;j++){
			centroidPosition[j]=m_resultPointSet[i*3+j];
		}
		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
			centroidPosition[0]+=m_resultPointSet[(*vectorIterator)*3];
			centroidPosition[1]+=m_resultPointSet[(*vectorIterator)*3+1];
			centroidPosition[2]+=m_resultPointSet[(*vectorIterator)*3+2];
		}
		for(int j=0;j<3;j++){
			centroidPosition[j]/=(numOfNearest+1);
		}
		localVariationVector=new double[(numOfNearest+1)*3];
		localVariationVector[0]=m_resultPointSet[i*3]-centroidPosition[0];
		localVariationVector[1]=m_resultPointSet[i*3+1]-centroidPosition[1];
		localVariationVector[2]=m_resultPointSet[i*3+2]-centroidPosition[2];
		vectorIterator=(*mapIterator).second.m_nearest.begin();
		for(int j=1;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,j++){
			localVariationVector[j*3]=m_resultPointSet[(*vectorIterator)*3]-centroidPosition[0];
			localVariationVector[j*3+1]=m_resultPointSet[(*vectorIterator)*3+1]-centroidPosition[1];
			localVariationVector[j*3+2]=m_resultPointSet[(*vectorIterator)*3+2]-centroidPosition[2];
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
		//m_gaussCurvature[i*2]=eval(0);
		//m_gaussCurvature[i*2+1]=eval(1);
		//m_meanAreas[i]=abs(eval(0)*eval(2));
			m_meanCurvature[i]=eval(2)/(eval(0)+eval(1)+eval(2));
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
		}

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

void RBFmeshless::ComputeNearestParameter(int indexPointSet)
{
	/*
	*	函数功能： 计算任一点处的积分半径，计算其邻域各点的局部坐标，
	*  变量说明：
	*  int indexPointSet        索引点	
	*  double m_localRBFradius;            径向基函数半径
	*  double* m_localOrdinate             各邻域点的局部坐标
	*  double* m_localIntegralPoint        积分点的坐标
	*  double  m_localIntegralRadius       积分域的半径
	*/
	double localPoint[3], localNormal[3],localMajorDirection[3],localMinorDirection[3];

	m_localRBFradius=0;
	for(int i=0;i<3;i++){
		localMajorDirection[i]=m_originalMajorDirection[indexPointSet*3+i];
		localMinorDirection[i]=m_originalMinorDirection[indexPointSet*3+i];
		localNormal[i]=m_originalNormals[indexPointSet*3+i];
	}
	std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(indexPointSet);
	std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	m_numOfNeighbors=(*mapIterator).second.m_nearest.size()+1;
	if(m_localOrdinate!=NULL){
		delete[] m_localOrdinate;
		m_localOrdinate=NULL;
	}
	m_localOrdinate=new double[m_numOfNeighbors*2];
	m_localOrdinate[0]=0;
	m_localOrdinate[1]=0;

	for(int i=1;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,i++){
		double vectorNeighborToPoint[3];
		double localRadius;
		double tempRadiusSquare;
		for(int j=0;j<3;j++){
			vectorNeighborToPoint[j]=m_resultPointSet[(*vectorIterator)*3+j]-m_resultPointSet[indexPointSet*3+j];
		}
		m_localOrdinate[i*2]=Vector3Vector(vectorNeighborToPoint,localMajorDirection);
		m_localOrdinate[i*2+1]=Vector3Vector(vectorNeighborToPoint,localMinorDirection); 
	//	m_localOrdinate[i*3+2]=Vector3Vector(vectorNeighborToPoint,localNormal);
		tempRadiusSquare=m_localOrdinate[i*2]*m_localOrdinate[i*2]+m_localOrdinate[i*2+1]*m_localOrdinate[i*2+1];
		if(tempRadiusSquare>m_localRBFradius)
			m_localRBFradius=tempRadiusSquare;

	}

	m_localRBFradius=sqrt(m_localRBFradius);


	

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

void RBFmeshless::ComputeTestFunction(int indexPointSet)
{
	/*
	*	计算测试函数在积分点的函数值及其导数值
	*/
	if(m_localIntegralPoints!=NULL){
		delete[] m_localIntegralPoints;
		m_localIntegralPoints=NULL;
	}
	m_localIntegralPoints=new double[m_numOfIntegralPoint*2];
	if(m_testFunction!=NULL){
		delete[] m_testFunction;
		m_testFunction=NULL;
	}
	m_testFunction=new double[m_numOfIntegralPoint];
	if(m_testFunctionDerOne!=NULL){
		delete[] m_testFunctionDerOne;
		m_testFunctionDerOne=NULL;
	}
	m_testFunctionDerOne=new double[m_numOfIntegralPoint];
	if(m_testFunctionDerTwo!=NULL){
		delete[] m_testFunctionDerTwo;
		m_testFunctionDerTwo=NULL;
	}
	m_testFunctionDerTwo=new double[m_numOfIntegralPoint];
	m_localIntergralRadius=m_localRBFradius/2;

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
	*	计算形函数及其导数  选择 c和r相等
	*/
	//////////////////////////////////////////////////////////////////////////
	double tempEXP,tempOneMinusEXP;// 计算测试函数时候的一些常数
	double temp=-1;
	tempEXP=exp(temp);
	tempOneMinusEXP=1-tempEXP;
	double integralRaiusSquare, radiusCSquare;

	integralRaiusSquare=m_localIntergralRadius*m_localIntergralRadius;
	radiusCSquare=integralRaiusSquare;

	for(int i=0;i<m_numOfIntegralPoint;i++){
		double tempDistanceSquare;//积分点与原点之间的距离的平方
		tempDistanceSquare=m_localIntegralPoints[i*2]*m_localIntegralPoints[i*2]
				+m_localIntegralPoints[i*2+1]*m_localIntegralPoints[i*2+1];
		if(tempDistanceSquare>integralRaiusSquare){ 
			m_testFunction[i]=0;
			m_testFunctionDerOne[i]=0;
			m_testFunctionDerTwo[i]=0;

		}
		else{
			double tempEXPDistance;
			tempEXPDistance=exp(-tempDistanceSquare/radiusCSquare);
			m_testFunction[i]=(tempEXPDistance-tempEXP)/tempOneMinusEXP;
			m_testFunctionDerOne[i]=-2*tempEXPDistance*m_localIntegralPoints[i*2]/radiusCSquare/tempOneMinusEXP;
			m_testFunctionDerTwo[i]=-2*tempEXPDistance*m_localIntegralPoints[i*2+1]/radiusCSquare/tempOneMinusEXP;
		}

	}
	return;
}

void RBFmeshless::ComputeMatrixNeighbor(int indexPointSet)
{
	/*
	*	函数功能：计算插值矩阵及其逆阵
	*  变量说明：
	*  double* m_matrixNeighbor; 插值矩阵
	*  double* m_matrixNeighborInv; 插值矩阵的逆阵  
	*/
	if(m_matrixNeighbor!=NULL){
		delete[] m_matrixNeighbor;
		m_matrixNeighbor=NULL;
	}
	m_matrixNeighbor=new double[m_numOfNeighbors*m_numOfNeighbors];
	if(m_matrixNeighborInv!=NULL){
		delete[] m_matrixNeighborInv;
		m_matrixNeighborInv=NULL;
	}
	m_matrixNeighborInv=new double[m_numOfNeighbors*m_numOfNeighbors];

	for(int i=0;i<m_numOfNeighbors;i++){
		// 计算每一行  在同一点不同基函数的值
		for(int j=0;j<m_numOfNeighbors;j++){
			double tempDistance;  // d/r
			double tempVector[2];
			for(int k=0;k<2;k++){
				tempVector[k]=m_localOrdinate[i*2+k]-m_localOrdinate[j*2+k];
			}
			tempDistance=tempVector[0]*tempVector[0]+tempVector[1]*tempVector[1];
			tempDistance=sqrt(tempDistance);
			if(tempDistance>m_localRBFradius){
				m_matrixNeighbor[i*m_numOfNeighbors+j]=0;
			}
			else{
				double oneMinusTempDistance;
				double tempOne,tempTwo;
				tempDistance/=m_localRBFradius;				
				oneMinusTempDistance=1-tempDistance;
				tempOne=oneMinusTempDistance*oneMinusTempDistance*oneMinusTempDistance*oneMinusTempDistance*oneMinusTempDistance;
				tempTwo=m_coffRBFunction[0]+m_coffRBFunction[1]*tempDistance
					 +m_coffRBFunction[2]*tempDistance*tempDistance
					 +m_coffRBFunction[3]*tempDistance*tempDistance*tempDistance
					 +m_coffRBFunction[4]*tempDistance*tempDistance*tempDistance*tempDistance;
				m_matrixNeighbor[i*m_numOfNeighbors+j]=tempOne*tempTwo;
			}
		}
	}


	mat_f8 mat_localMatrixNeighbor(m_numOfNeighbors,m_numOfNeighbors);
	mat_f8 mat_localMatrixNeighborInv(m_numOfNeighbors,m_numOfNeighbors);
	
	for(int i=0;i<m_numOfNeighbors;i++){
		for(int j=0;j<m_numOfNeighbors;j++){
			mat_localMatrixNeighbor(i,j)=m_matrixNeighbor[i*m_numOfNeighbors+j];
		}
	}

	if(!matrix_inverse(mat_localMatrixNeighbor)){
		//MessageBox("求逆矩阵出错");            
	}
	
	for(int i=0;i<m_numOfNeighbors;i++){
		for(int j=0;j<m_numOfNeighbors;j++){
			m_matrixNeighborInv[i*m_numOfNeighbors+j]=mat_localMatrixNeighbor(i,j);
		}
	}

	return;
}

void RBFmeshless::ComputeMatrixRelatedSampling(int indexPointSet)
{
	/*
	*	函数功能： 计算与采样点有关的矩阵 
	*  变量说明： 
	*  int indexPointSet        索引点
	*  double* m_RBFfunction;        每一行为在每一个邻域点的径向基函数在每一个积分点的函数值
	*/
	//////////////////////////////////////////////////////////////////////////
	/*
	 *	分配内存空间
	 */
	if(m_RBFfunction!=NULL){
		delete[] m_RBFfunction;
		m_RBFfunction=NULL;
	}
	m_RBFfunction=new double[m_numOfNeighbors*m_numOfIntegralPoint];
	if(m_RBFfunctionDerOne!=NULL){
		delete[] m_RBFfunctionDerOne;
		m_RBFfunctionDerOne;
	}
	m_RBFfunctionDerOne=new double[m_numOfNeighbors*m_numOfIntegralPoint];
	if(m_RBFfunctionDerTwo!=NULL){
		delete[] m_RBFfunctionDerTwo;
		m_RBFfunctionDerTwo=NULL;
	}
	m_RBFfunctionDerTwo=new double[m_numOfNeighbors*m_numOfIntegralPoint];
	
	for(int i=0;i<m_numOfNeighbors;i++){
		for(int j=0;j<m_numOfIntegralPoint;j++){
			double tempDistanceD;  // d/r
			double tempVector[2];
			for(int k=0;k<2;k++){
				tempVector[k]=m_localIntegralPoints[j*2+k]-m_localOrdinate[i*2+k];
			}
			tempDistanceD=tempVector[0]*tempVector[0]+tempVector[1]*tempVector[1];
			tempDistanceD=sqrt(tempDistanceD);
			
			
			if(tempDistanceD>m_localRBFradius){
				m_RBFfunction[i*m_numOfIntegralPoint+j]=0;
				m_RBFfunctionDerTwo[i*m_numOfIntegralPoint+j]=0;
				m_RBFfunctionDerOne[i*m_numOfIntegralPoint+j]=0;
			}
			else{
				double oneMinusTempDistance, tempDistance;
				double tempOne,tempTwo;
				for(int k=0;k<2;k++){
					tempVector[k]/=tempDistanceD;
				}
				tempDistance=tempDistanceD/m_localRBFradius;		
				oneMinusTempDistance=1-tempDistance;
				tempOne=oneMinusTempDistance*oneMinusTempDistance*oneMinusTempDistance*oneMinusTempDistance*oneMinusTempDistance;
				tempTwo=m_coffRBFunction[0]+m_coffRBFunction[1]*tempDistance
					+m_coffRBFunction[2]*tempDistance*tempDistance
					+m_coffRBFunction[3]*tempDistance*tempDistance*tempDistance
					+m_coffRBFunction[4]*tempDistance*tempDistance*tempDistance*tempDistance;
				m_RBFfunction[i*m_numOfIntegralPoint+j]=tempOne*tempTwo;    
				m_RBFfunctionDerOne[i*m_numOfIntegralPoint+j]=-5*(tempOne/oneMinusTempDistance)*tempTwo*tempVector[0]/m_localRBFradius;
				m_RBFfunctionDerTwo[i*m_numOfIntegralPoint+j]=-5*(tempOne/oneMinusTempDistance)*tempTwo*tempVector[1]/m_localRBFradius;
                tempTwo=40/m_localRBFradius
					+96*tempDistanceD/m_localRBFradius/m_localRBFradius
					+75*tempDistanceD*tempDistanceD/m_localRBFradius/m_localRBFradius/m_localRBFradius
					+20*tempDistanceD*tempDistanceD*tempDistanceD/m_localRBFradius/m_localRBFradius/m_localRBFradius/m_localRBFradius;
				m_RBFfunctionDerOne[i*m_numOfIntegralPoint+j]+=tempOne*tempTwo*tempVector[0];
				m_RBFfunctionDerTwo[i*m_numOfIntegralPoint+j]+=tempOne*tempTwo*tempVector[1];
			}
		}
	}
	return;
}

void RBFmeshless::ComputeShapeFunction(int indexPointSet)
{
	/*
	*	函数功能：计算形函数及其导数
	*  变量说明：
	*  int indexPointSet        索引点
	*  double* m_shapeFunction                各邻域点的形函数在各采样点的值
	*/
	//////////////////////////////////////////////////////////////////////////
	/*
	 *	分配内存
	 */
	//////////////////////////////////////////////////////////////////////////
	if(m_shapeFunction!=NULL){
		delete[] m_shapeFunction;
		m_shapeFunction=NULL;
	}
	m_shapeFunction=new double[m_numOfNeighbors*m_numOfIntegralPoint];
	if(m_shapeFunctionDerOne!=NULL){
		delete[] m_shapeFunctionDerOne;
		m_shapeFunctionDerOne;
	}
	m_shapeFunctionDerOne=new double[m_numOfNeighbors*m_numOfIntegralPoint];
	if(m_shapeFunctionDerTwo!=NULL){
		delete[] m_shapeFunctionDerTwo;
		m_shapeFunctionDerTwo=NULL;
	}
	m_shapeFunctionDerTwo=new double[m_numOfNeighbors*m_numOfIntegralPoint];
	//////////////////////////////////////////////////////////////////////////
	/*
	*	计算形函数及其导数
	*/
	//////////////////////////////////////////////////////////////////////////
	for(int i=0;i<m_numOfNeighbors;i++){
		//////////////////////////////////////////////////////////////////////////
		/*
		*	计算节点i的形函数
		*/
		//////////////////////////////////////////////////////////////////////////
		
		for(int j=0;j<m_numOfIntegralPoint;j++){
			m_shapeFunction[i*m_numOfIntegralPoint+j]=0;
			m_shapeFunctionDerOne[i*m_numOfIntegralPoint+j]=0;
			m_shapeFunctionDerTwo[i*m_numOfIntegralPoint+j]=0;
			for(int k=0;k<m_numOfNeighbors;k++){
				m_shapeFunction[i*m_numOfIntegralPoint+j]+=m_RBFfunction[k*m_numOfIntegralPoint+j]
														*m_matrixNeighborInv[k*m_numOfNeighbors+i];
				m_shapeFunctionDerOne[i*m_numOfIntegralPoint+j]+=m_RBFfunctionDerOne[k*m_numOfIntegralPoint+j]
														*m_matrixNeighborInv[k*m_numOfNeighbors+i];
				m_shapeFunctionDerTwo[i*m_numOfIntegralPoint+j]+=m_RBFfunctionDerTwo[k*m_numOfIntegralPoint+j]
														*m_matrixNeighborInv[k*m_numOfNeighbors+i];	

			}
		}
	}
	return;    
}


void RBFmeshless::ComputeManifoldMetric(int indexPointSet)
{
	/*
	*	函数功能： 计算流形标准的代数行列式
	*  变量说明：
	*  int indexPointSet        索引点
	*  double* m_determinentMetric   在采样点处流形标准的代数行列式
	*/
	//////////////////////////////////////////////////////////////////////////
	/*
	*	分配内存
	*/
	//////////////////////////////////////////////////////////////////////////
	if(m_determinentMetric!=NULL){
		delete[] m_determinentMetric;
		m_determinentMetric=NULL;
	}
	m_determinentMetric=new double[m_numOfIntegralPoint];
	if(m_metricMatrix!=NULL){
		delete[] m_metricMatrix;
		m_metricMatrix=NULL;
	}
	m_metricMatrix=new double[m_numOfIntegralPoint*4];
	if(m_metricMatrixInve!=NULL){
		delete[] m_metricMatrixInve;
		m_metricMatrixInve=NULL;
	}
	m_metricMatrixInve=new double[m_numOfIntegralPoint*4];
	//////////////////////////////////////////////////////////////////////////
	/*
	*	计算流形标准矩阵
	*/
	//////////////////////////////////////////////////////////////////////////
	for(int i=0;i<m_numOfIntegralPoint;i++){
		double tempVectorOne[3],tempVectorTwo[3];// 参考3.5
		for(int j=0;j<3;j++){
			tempVectorOne[j]=m_shapeFunctionDerOne[i]*m_resultPointSet[indexPointSet*3+j];
			tempVectorTwo[j]=m_shapeFunctionDerTwo[i]*m_resultPointSet[indexPointSet*3+j];			
		}
		std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(indexPointSet);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		for(int j=1;j<m_numOfNeighbors;j++,vectorIterator++){
			for(int k=0;k<3;k++){
				tempVectorOne[k]+=m_shapeFunctionDerOne[j*m_numOfIntegralPoint+i]*m_resultPointSet[(*vectorIterator)*3+k];
				tempVectorTwo[k]+=m_shapeFunctionDerTwo[j*m_numOfIntegralPoint+i]*m_resultPointSet[(*vectorIterator)*3+k];
			}			
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
	return;


}

 
void RBFmeshless::ComputeStiffMatrix(int indexPointSet)
{
	/*
	*	函数功能： 计算刚度矩阵和荷载向量
	*  变量说明： 
	*  int indexPointSet                    索引点
	*  double* m_loadVector                  荷载向量
	*  StiffMatrix* m_stiffMatrix           刚度矩阵
	*/
	//////////////////////////////////////////////////////////////////////////
	/*
	 *	//首先计算 对角线上的元素  及质量矩阵
	 */
	//////////////////////////////////////////////////////////////////////////
	
	
	double namdaLoadConstant=m_namda*m_loadConstant;
	double elementOne=0;
	double elementTwo=0;
	double elementStiffMatrix;
	//////////////////////////////////////////////////////////////////////////
	// 求积分
	for(int j=0;j<1;j++){
		elementOne=0;
		for(int i=0;i<m_numOfIntegralPoint;i++){
			
			//	elementOne+=m_shapeFunction[i]*m_testFunction[i]*sqrt(m_determinentMetric[i])*m_ordinateAndWeightCoff[i*3+2];
			elementOne+=m_shapeFunction[j*m_numOfNeighbors+i]*m_testFunction[i]*m_ordinateAndWeightCoff[i*3+2];
			//	elementTwo+=(m_shapeFunctionDerOne[i]*m_testFunctionDerOne[i]+m_shapeFunctionDerTwo[i]*m_testFunctionDerTwo[i])*m_ordinateAndWeightCoff[i*3+2];
			/*	double tempElementTwo=0;
			for(int j=0;j<3;j++){
			tempElementTwo+=m_gradientShapeFunction[i*3+j]*m_gradientTestFunction[j];
			}*/

			/*	double tempElementTwo=0;
			tempElementTwo+=m_shapeFunctionDerOne[i]*m_testFunctionDerOne[i]
			*(m_metricMatrixInve[i*4]*m_metricMatrixInve[i*4]+m_metricMatrixInve[i*4+2]*m_metricMatrixInve[i*4+2]);
			tempElementTwo+=m_shapeFunctionDerOne[i]*m_testFunctionDerTwo[i]
			*(m_metricMatrixInve[i*4]*m_metricMatrixInve[i*4+1]+m_metricMatrixInve[i*4+2]*m_metricMatrixInve[i*4+3]);
			tempElementTwo+=m_shapeFunctionDerTwo[i]*m_testFunctionDerOne[i]
			*(m_metricMatrixInve[i*4+2]*m_metricMatrixInve[i*4]+m_metricMatrixInve[i*4+3]*m_metricMatrixInve[i*4+2]);
			tempElementTwo+=m_shapeFunctionDerTwo[i]*m_testFunctionDerTwo[i]
			*(m_metricMatrixInve[i*4+2]*m_metricMatrixInve[i*4+2]+m_metricMatrixInve[i*4+3]*m_metricMatrixInve[i*4+3]);
			elementTwo+=tempElementTwo*m_ordinateAndWeightCoff[i*3+2];		*/
		}
		m_massMatrix[indexPointSet].m_elementStiffMarix.push_back(elementOne);
		//	elementStiffMatrix=elementOne+m_namda*elementTwo;
	//	m_stiffMatrix[indexPointSet].m_elementStiffMarix.push_back(0);
		//m_loadVector[indexPointSet]=(1-namdaLoadConstant)*elementOne;
		/*for(int i=0;i<3;i++){
		m_loadVector[indexPointSet]+=(namdaLoadConstant*elementOne*m_originalPointSet[indexPointSet*3+i]
		+(1-namdaLoadConstant)*elementOne*m_resultPointSet[indexPointSet*3+i])*m_originalNormals[indexPointSet*3+i];

		}*/
	}
	
	//////////////////////////////////////////////////////////////////////////
	/*
	*	//计算刚度矩阵
	*/
	//////////////////////////////////////////////////////////////////////////

	std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(indexPointSet);
	std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	for(int i=0;i<m_numOfNeighbors;i++){
		//elementOne=0;
		elementTwo=0;
		//////////////////////////////////////////////////////////////////////////
		// 求积分
		for(int j=0;j<m_numOfIntegralPoint;j++){
			//elementOne+=m_shapeFunction[i*m_numOfIntegralPoint+j]*m_testFunction[j]*sqrt(m_determinentMetric[j])*m_ordinateAndWeightCoff[j*3+2];
			//elementOne+=m_shapeFunction[i*m_numOfIntegralPoint+j]*m_testFunction[j]*m_ordinateAndWeightCoff[j*3+2];
			elementTwo+=(m_shapeFunctionDerOne[i*m_numOfIntegralPoint+j]*m_testFunctionDerOne[j]+m_shapeFunctionDerTwo[i*m_numOfIntegralPoint+j]*m_testFunctionDerTwo[j])*m_ordinateAndWeightCoff[j*3+2];
			/*	double tempElementTwo=0;
			for(int k=0;k<3;k++){
			tempElementTwo+=m_gradientShapeFunction[i*m_numOfIntegralPoint*3+j*3+k]*m_gradientTestFunction[j*3+k];
			}*/
			//	elementTwo+=m_shapeFunctionDerivativeOne[i*m_numOfIntegralPoint+j]*m_TestFunctionDerivativeOne[j]+m_shapeFunctionDerivativeTwo[i*m_numOfIntegralPoint+j]*m_TestFunctionDerivativeTwo[j];
			/*double tempElementTwo=0;
			tempElementTwo+=m_shapeFunctionDerOne[i*m_numOfIntegralPoint+j]*m_testFunctionDerOne[j]
					*(m_metricMatrixInve[j*4]*m_metricMatrixInve[j*4]+m_metricMatrixInve[j*4+2]*m_metricMatrixInve[j*4+2]);
			tempElementTwo+=m_shapeFunctionDerOne[i*m_numOfIntegralPoint+j]*m_testFunctionDerTwo[j]
					*(m_metricMatrixInve[j*4]*m_metricMatrixInve[j*4+1]+m_metricMatrixInve[j*4+2]*m_metricMatrixInve[j*4+3]);
			tempElementTwo+=m_shapeFunctionDerTwo[i*m_numOfIntegralPoint+j]*m_testFunctionDerOne[j]
					*(m_metricMatrixInve[j*4+2]*m_metricMatrixInve[j*4]+m_metricMatrixInve[j*4+3]*m_metricMatrixInve[j*4+2]);
			tempElementTwo+=m_shapeFunctionDerTwo[i*m_numOfIntegralPoint+j]*m_testFunctionDerTwo[j]
					*(m_metricMatrixInve[j*4+2]*m_metricMatrixInve[j*4+2]+m_metricMatrixInve[j*4+3]*m_metricMatrixInve[j*4+3]);
			elementTwo+=tempElementTwo*m_ordinateAndWeightCoff[j*3+2];*/
		}
		elementStiffMatrix=m_namda*elementTwo;
		m_stiffMatrix[indexPointSet].m_elementStiffMarix.push_back(elementStiffMatrix);
	/*	for(int k=0;k<3;k++){
			m_loadVector[indexPointSet*3+k]+=namdaLoadConstant*elementOne*m_originalPointSet[3*(*vectorIterator)+k]
					+(1-namdaLoadConstant)*elementOne*m_resultPointSet[3*(*vectorIterator)+k];
		}*/
		vectorIterator++;		

	}
	return;
}


void RBFmeshless::EstimateNormalDirection()
{
	/*
	*	函数功能： 估计法向的准确方向
	*/
	//首先确定z值最大的那一点的法向
	int maxZpoint=0;
	for(int i=0;i<m_numOfPoints;i++){
		if(m_resultPointSet[i*3+2]>m_resultPointSet[maxZpoint*3+2])
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

void RBFmeshless::CalculateLinearSystem()
{
	/*
	*	函数功能： 解线性方程
	*/

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
		//threshold=0;
	}
	return;
}


void  RBFmeshless::ComputeGradientShapeFunction()
{
	/*
	*	函数功能： 计算每一个积分点形函数的梯度和测试函数的梯度
	*  变量说明：
	*  double* m_gradientShapeFunction   形函数的梯度
	*  double* m_gradientTestFunction    测试函数的梯度
	*/

}

void RBFmeshless::CalculateLinearSystemGauss()
{
	/*
	*	函数功能： 高斯消元法解方程；
	*/
}

void RBFmeshless:: AssembleStiffMatrix(){
	//////////////////////////////////////////////////////////////////////////	
	/*
	*	函数功能： 组装刚度矩阵
	*/
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	/*
	 *	计算方程组右边和刚度矩阵，先计算右边，再改刚度矩阵 组合质量矩阵和刚度矩阵
	 */
	for(int i=0;i<m_numOfPoints;i++){
		m_loadVector[i]=0;
		std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
		std::vector<double>::iterator iteratorMassMatrix=m_massMatrix[i].m_elementStiffMarix.begin();
		std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		if((*mapIterator).second.m_nearest.size()<m_kNearest){
			(*iteratorMassMatrix)=1;
			for(iteratorStiffMatrix++,iteratorMassMatrix++;iteratorStiffMatrix!=m_stiffMatrix[i].m_elementStiffMarix.end();iteratorStiffMatrix++,vectorIterator++,iteratorMassMatrix++){
				(*iteratorStiffMatrix)=0;
			}
			m_loadVector[i]=0;
		}
		//(*iteratorStiffMatrix)=((*iteratorStiffMatrix)+(*iteratorMassMatrix));	
		else{
			for(iteratorStiffMatrix++,iteratorMassMatrix++;iteratorStiffMatrix!=m_stiffMatrix[i].m_elementStiffMarix.end();iteratorStiffMatrix++,vectorIterator++,iteratorMassMatrix++){
				double distanceNormal; //两个法向的点乘
				double distanceProjection;  //领域点到切平面的距离
				distanceProjection=(m_resultPointSet[(*vectorIterator)*3]-m_resultPointSet[i*3])*m_originalNormals[i*3]
				+(m_resultPointSet[(*vectorIterator)*3+1]-m_resultPointSet[i*3+1])*m_originalNormals[i*3+1]
				+(m_resultPointSet[(*vectorIterator)*3+2]-m_resultPointSet[i*3+2])*m_originalNormals[i*3+2];
				distanceNormal=m_originalNormals[(*vectorIterator)*3]*m_originalNormals[i*3]
				+m_originalNormals[(*vectorIterator)*3+1]*m_originalNormals[i*3+1]
				+m_originalNormals[(*vectorIterator)*3+2]*m_originalNormals[i*3+2];
				if(distanceNormal<0)
					distanceNormal=0;
				m_loadVector[i]-=(*iteratorStiffMatrix)*distanceProjection*distanceNormal;
				(*iteratorStiffMatrix)=((*iteratorStiffMatrix));	
				/*m_loadVector[i]-=(*iteratorStiffMatrix)*distanceProjection;
				(*iteratorStiffMatrix)=((*iteratorStiffMatrix))/distanceNormal;	*/	
				//	(*iteratorStiffMatrix)=((*iteratorStiffMatrix)+(*iteratorMassMatrix))*distanceNormal;		

			}

		}
		
	}




	//////////////////////////////////////////////////////////////////////////
	/*
	 *	对角化刚度矩阵
	 */

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
		/*for(int j=0;j<3;j++){
			m_loadVector[i*3+j]=namdaLoadConstant*(*iteratorMassMatrix)*m_originalPointSet[i*3+j]
			+(1-namdaLoadConstant)*(*iteratorMassMatrix)*m_resultPointSet[i*3+j];

		}*/
		(*iteratorStiffMatrix)=(*iteratorMassMatrix)-sumOfElement;
	}
	return;
}
//////////////////////////////////////////////////////////////////////////
/*
*	调整刚度矩阵
*/
void RBFmeshless::AdjustStiffMatrix()
{
	for (int i=0;i<m_numOfPoints;i++){
		double majorElement,minorElement;
		std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
		majorElement=abs((*iteratorStiffMatrix));
		minorElement=0;
		for(iteratorStiffMatrix++;iteratorStiffMatrix!=m_stiffMatrix[i].m_elementStiffMarix.end();iteratorStiffMatrix++){
			minorElement+=abs((*iteratorStiffMatrix));
		}
		if(majorElement>minorElement){
			continue;
		}
		else{
			iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
			std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
			std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
			for(iteratorStiffMatrix++;iteratorStiffMatrix!=m_stiffMatrix[i].m_elementStiffMarix.end();iteratorStiffMatrix++,vectorIterator++){
				(*iteratorStiffMatrix)=0;
				std::map<int, KnearestField>::iterator mapIteratorIner=m_mapKnearest.find((*vectorIterator));
				std::vector<int>::iterator vectorIteratorIner=(*mapIteratorIner).second.m_nearest.begin();
				std::vector<double>::iterator iteratorStiffMatrixIner=m_stiffMatrix[(*vectorIterator)].m_elementStiffMarix.begin();
				for(iteratorStiffMatrixIner++;iteratorStiffMatrixIner!=m_stiffMatrix[(*vectorIterator)].m_elementStiffMarix.end();iteratorStiffMatrixIner++,vectorIteratorIner++){
					if((*vectorIteratorIner)==i){
                        (*iteratorStiffMatrixIner)=0;
						break;
					}
				}

			}

		}

		
	}

	return;	
}
//////////////////////////////////////////////////////////////////////////
/*
*	pbcg 求稀疏线性方程组
*/
//////////////////////////////////////////////////////////////////////////

void RBFmeshless::linbcg(double* b,double* x,int itol,double tol,int itmax, int & iter, double &err)
{
	double ak,akden,bk,bkden=1.0,bknum,bnrm,dxnrm,xnrm,zminrm,znrm;
	const double EPS=1.0E-14;
	int j;
	int n=m_numOfPoints;
	double* p=new double[n];
	double* pp=new double[n];
	double* r=new double[n];
	double* rr=new double[n];
	double* z=new double[n];
	double* zz=new double[n];
	iter=0;
	atimes(x,r,0);
	for(j=0;j<n;j++){
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	//atimes(r,rr,0);
	if(itol==1){
		bnrm=snrm(b,itol);
		asolve(r,z,0);
	}
	else if (itol==2){
		asolve(b,z,0);
		bnrm=snrm(z,itol);
		asolve(r,z,0);
	}
	else if(itol==3 || itol==4){
		asolve(b,z,0);
		bnrm=snrm(z,itol);
		asolve(r,z,0);
		znrm=snrm(z,itol);
	}

	while(iter <itmax){
		++iter;
		asolve(rr,zz,1);
		for(bknum=0.0,j=0;j<n;j++) bknum+=z[j]*rr[j];
		if(iter==1){
			for(j=0;j<n;j++){
				p[j]=z[j];
				pp[j]=zz[j];
			}
		}
		else{
			bk=bknum/bkden;
			for(j=0;j<n;j++){
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;
		atimes(p,z,0);
		for(akden=0.0,j=0;j<n;j++) akden+=z[j]*pp[j];
		ak=bknum/akden;
		atimes(pp,zz,1);
		for(j=0;j<n;j++){
			x[j]+=ak*p[j];
			r[j]-=ak*z[j];
			rr[j]-=ak*zz[j];
		}
		asolve(r,z,0);
		if(itol==1)
			err=snrm(r,itol)/bnrm;
		else if(itol==2)
			err=snrm(z,itol)/bnrm;
		else if (itol==3||itol||4){
			zminrm=znrm;
			znrm=snrm(z,itol);
			if(fabs(zminrm-znrm)>EPS*znrm){
				dxnrm=fabs(ak)*snrm(p,itol);
				err=znrm/fabs(zminrm-znrm)*dxnrm;
			}else{
				err=znrm/bnrm;
				continue;
			}
			xnrm=snrm(x,itol);
			if(err<=0.5*xnrm) err/=xnrm;
			else{
				err=znrm/bnrm;
				continue;
			}

		}
		if(err<=tol) break;

	}

}
double RBFmeshless::snrm(double* sx,const int itol)
{
	int i,isamax;
	double ans;
	int n=m_numOfPoints;
	if(itol<=3){
		ans=0.0;
		for(i=0;i<n;i++) ans+=sx[i]*sx[i];
		return sqrt(ans);
	}else{
		isamax=0;
		for(i=0;i<n;i++){
			if(fabs(sx[i])>fabs(sx[isamax])) isamax=i;
		}
		return fabs(sx[isamax]);
	}
}
void RBFmeshless::atimes(double* x,double* r,const int itrnsp)
{
	if(itrnsp==0){
		for (int i=0;i<m_numOfPoints;i++){
			std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
			std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
			std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
			r[i]=(*iteratorStiffMatrix)*x[i];

			for(iteratorStiffMatrix++;iteratorStiffMatrix!=m_stiffMatrix[i].m_elementStiffMarix.end();iteratorStiffMatrix++,vectorIterator++){
				r[i]+=(*iteratorStiffMatrix)*x[(*vectorIterator)];				
			}		
		}
	}
	else{
		for(int i=0;i<m_numOfPoints;i++){
			r[i]=0;
		}
		for (int i=0;i<m_numOfPoints;i++){
			std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
			std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
			std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
			r[i]+=(*iteratorStiffMatrix)*x[i];

			for(iteratorStiffMatrix++;iteratorStiffMatrix!=m_stiffMatrix[i].m_elementStiffMarix.end();iteratorStiffMatrix++,vectorIterator++){
				r[(*vectorIterator)]+=(*iteratorStiffMatrix)*x[i];				
			}		
		}

	}
	return;

}
void RBFmeshless::asolve(double* b,double* x,const int itrnsp)
{
	int i;
	int n=m_numOfPoints;
	for(i=0;i<n;i++){
		std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMarix.begin();
		x[i]=((*iteratorStiffMatrix)!=0.0 ? b[i]/(*iteratorStiffMatrix):b[i]);
	}
	return;
}


