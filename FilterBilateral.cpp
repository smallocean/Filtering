// FilterBilateral.cpp: implementation of the FilterBilateral class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Mcube.h"
#include "FilterBilateral.h"
#include ".\mathlib\mathlib.h"
#include "Matrix.h"
#include <time.h>
//#include "GenerateNoise.h"

using namespace MATHLIB;


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FilterBilateral::FilterBilateral()
{
	m_originalPointSet=NULL;
	m_originalColors=NULL;
	m_originalNormals=NULL;
    m_numOfPoints=0;
	m_resultPointSet=NULL;
	m_resultColors=NULL;
	m_resultNormals=NULL;
    m_filterNormals=NULL;
    m_originalColorGradients[0]=NULL;
	m_originalColorGradients[1]=NULL;
	m_originalColorGradients[2]=NULL;
	m_resultColorGradients[0]=NULL;
	m_resultColorGradients[1]=NULL;
	m_resultColorGradients[2]=NULL;
	m_meanShiftHposition=0.0;
	m_meanShiftHnormal=0.0;
	m_meanShiftStopNormal=0.0;
	m_knearest=0;
	m_radius=0;

}

FilterBilateral::~FilterBilateral()
{
	DeleteFilterBilateral();	
}

void FilterBilateral::GetFilterBilateral(int numOfPoints, float* pointSet, float* colors)
{
	long nearestStart,nearestFinish;
	long filterStart,filterFinish;
	m_numOfPoints=numOfPoints;
	m_originalPointSet=new float[m_numOfPoints*3];
	m_originalColors=new float[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints*3;i++){
		m_originalPointSet[i]=pointSet[i];
		m_originalColors[i]=colors[i];
	}
	m_meanShiftHposition=2;
	m_meanShiftHnormal=0.3;
	m_meanShiftStopNormal=0.01;
	m_knearest=20;
	m_radius=4;
	for(int i=0;i<3;i++){
		m_thresholdColor[i]=0.1;
	}
	//AddPositionGaussianNoise();
	//AddColorGaussianNorse();
	nearestStart=clock();
	ComputeMapKnearest();
	nearestFinish=clock();
 	nearestTime=nearestFinish-nearestStart;
	filterStart=clock();
	CalculateOriginalNormals();
	MeanShiftNormals();
	CalculateColorGradient();
	//CalculateBilateralFilterOne();
	//CalculateBilateralFilterTwo();
	CalculateBilateralFilterThree();
	//CalculateBilateralFilterFour();
	//CalculateBilateralFilterFive();
	filterFinish=clock();
	filterTime=filterFinish-filterStart;
	m_mapKnearest.clear();
	for(int i=0;i<3;i++){
		delete[] m_originalColorGradients[i];
		m_originalColorGradients[i]=NULL;
	}

}

void FilterBilateral::DeleteFilterBilateral()
{
	if(m_originalPointSet!=NULL)
		delete[] m_originalPointSet;
	m_originalPointSet=NULL;
	if(m_originalColors!=NULL)
		delete[] m_originalColors;
	m_originalColors=NULL;
	if(m_originalNormals!=NULL)
		delete[] m_originalNormals;
	m_originalNormals=NULL;
    m_numOfPoints=0;
	if(m_resultPointSet!=NULL)
		delete[] m_resultPointSet;
	m_resultPointSet=NULL;
	if(m_resultColors!=NULL)
		delete[] m_resultColors;
	m_resultColors=NULL;
	if(m_resultNormals!=NULL)
		delete[] m_resultNormals;
	m_resultNormals=NULL;
	if(m_filterNormals!=NULL)
		delete[] m_filterNormals;
    m_filterNormals=NULL;
	for(int i=0;i<3;i++){
		if(m_originalColorGradients[i]!=NULL)
			delete[] m_originalColorGradients[i];
        m_originalColorGradients[i]=NULL;
		if(m_resultColorGradients[i]!=NULL)
			delete[] m_resultColorGradients[i];
		m_resultColorGradients[i]=NULL;
	}
}
void FilterBilateral::ComputeMapKnearest()
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
			if(numOfKnearest>m_knearest){
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
    // 测试代码
	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\knearest.txt";
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
	
	
	//////////////////////////////////////////////////////////////////////////
}

void FilterBilateral::CalculateOriginalNormals()
{
	
	float centroidPosition[3];//重心位置
	float* localVariationVector;//重心到点的位置的向量
	int numOfNearest;

	mat_f8 mat_covaMatrix(3, 3);
	m_originalNormals=new float[m_numOfPoints*3];
	
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
		delete[] localVariationVector;
	}

	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\originalNormal.txt";
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

}

void FilterBilateral::MeanShiftNormals()
{
	float distance_position;
	float distance_normal;
	float position_kernel,range_kernel;
	float sum_kernel;
	float length_normal_vector;
	int num_stop=0;//已经停止的数
	bool* flag_stop=new bool[m_numOfPoints];
	for(int i=0;i<m_numOfPoints;i++){
		flag_stop[i]=false;
	}

    m_filterNormals=new float[m_numOfPoints*3];
	while(num_stop<m_numOfPoints){
		for(int i=0;i<m_numOfPoints;i++){
			if(!flag_stop[i]){
				sum_kernel=0;
				for(int j=0;j<3;j++){
                    m_filterNormals[i*3+j]=0;
				}
			
				std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(i);
				//int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
				std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
				for(;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++){
					distance_position=(m_originalPointSet[i*3]-m_originalPointSet[(*vectorNearestIterator)*3])*(m_originalPointSet[i*3]-m_originalPointSet[(*vectorNearestIterator)*3])
						+(m_originalPointSet[i*3+1]-m_originalPointSet[(*vectorNearestIterator)*3+1])*(m_originalPointSet[i*3+1]-m_originalPointSet[(*vectorNearestIterator)*3+1])
						+(m_originalPointSet[i*3+2]-m_originalPointSet[(*vectorNearestIterator)*3+2])*(m_originalPointSet[i*3+2]-m_originalPointSet[(*vectorNearestIterator)*3+2]);
					distance_normal=(m_originalNormals[i*3]-m_originalNormals[(*vectorNearestIterator)*3])*(m_originalNormals[i*3]-m_originalNormals[(*vectorNearestIterator)*3])
						+(m_originalNormals[i*3+1]-m_originalNormals[(*vectorNearestIterator)*3+1])*(m_originalNormals[i*3+1]-m_originalNormals[(*vectorNearestIterator)*3+1])
						+(m_originalNormals[i*3+2]-m_originalNormals[(*vectorNearestIterator)*3+2])*(m_originalNormals[i*3+2]-m_originalNormals[(*vectorNearestIterator)*3+2]);
					position_kernel=exp(-0.5*distance_position/m_meanShiftHposition);
					range_kernel=exp(-0.5*distance_normal/m_meanShiftHnormal);
					sum_kernel+=position_kernel*range_kernel;
					
					for(int j=0;j<3;j++){
						m_filterNormals[i*3+j]+=m_originalNormals[(*vectorNearestIterator)*3+j]*position_kernel*range_kernel;
					}
				}
				sum_kernel+=exp(0.0)*exp(0.0);
				for(int j=0;j<3;j++){
					m_filterNormals[i*3+j]+=m_originalNormals[i*3+j]*exp(0.0)*exp(0.0);
				}
				for(int j=0;j<3;j++){
					m_filterNormals[i*3+j]/=sum_kernel;
				}
				//归一化 m_originalFilterNormals
				float tempNormal[3];
				for(int j=0;j<3;j++){
					tempNormal[j]=m_filterNormals[i*3+j];
				}
				length_normal_vector=sqrt(Vector3Vector(tempNormal,tempNormal));
				for(int j=0;j<3;j++){
					m_filterNormals[i*3+j]/=length_normal_vector;
				}

				length_normal_vector=((m_filterNormals[i*3]-m_originalNormals[i*3])*(m_filterNormals[i*3]-m_originalNormals[i*3])
					+(m_filterNormals[i*3+1]-m_originalNormals[i*3+1])*(m_filterNormals[i*3+1]-m_originalNormals[i*3+1])
					+(m_filterNormals[i*3+2]-m_originalNormals[i*3+2])*(m_filterNormals[i*3+2]-m_originalNormals[i*3+2]));
				
				if(((m_filterNormals[i*3]-m_originalNormals[i*3])*(m_filterNormals[i*3]-m_originalNormals[i*3])
					+(m_filterNormals[i*3+1]-m_originalNormals[i*3+1])*(m_filterNormals[i*3+1]-m_originalNormals[i*3+1])
					+(m_filterNormals[i*3+2]-m_originalNormals[i*3+2])*(m_filterNormals[i*3+2]-m_originalNormals[i*3+2]))<m_meanShiftStopNormal){
						flag_stop[i]=true;
						num_stop+=1;
					}
				else{
					m_originalNormals[i*3]=m_filterNormals[i*3];
					m_originalNormals[i*3+1]=m_filterNormals[i*3+1];
					m_originalNormals[i*3+2]=m_filterNormals[i*3+2];
				}
				
			}
		}

	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;
	//////////////////////////////////////////////////////////////////////////
	// 测试代码
	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\filterNormal.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f\n",m_filterNormals[i*3],m_filterNormals[i*3+1],m_filterNormals[i*3+2]);

		}
		fclose(fpout);
	}
	//////////////////////////////////////////////////////////////////////////
	
}

void FilterBilateral::CalculateColorGradient()
{
	float vectorOne[3], vectorTwo[3], localNormal[3];
	float coffecientOne,coffecientTwo, cofecientD;// 梯度 
	float coffecientMatrix[3][3], coffecientConstantRed[3],coffecientConstantGreen[3],coffecientConstantBlue[3];
	float localPosition[3];
	float localEXOne,localEXTwo;
	m_originalColorGradients[0]=new float[m_numOfPoints*4];
	m_originalColorGradients[1]=new float[m_numOfPoints*4];
	m_originalColorGradients[2]=new float[m_numOfPoints*4];
	for(int i=0;i<m_numOfPoints;i++){
		//首先求与法向量垂直的单位向量
		localNormal[0]=m_filterNormals[i*3];
		localNormal[1]=m_filterNormals[i*3+1];
		localNormal[2]=m_filterNormals[i*3+2];
		Perpendicular(localNormal,vectorOne);
		Cross3(localNormal, vectorOne,vectorTwo);
		// 计算方程组的系数
		for(int j=0;j<3;j++){
			//coffecientConstant[j]=0;
			for(int k=0;k<3;k++){
				coffecientMatrix[j][k]=0;
			}
		}
		std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		localPosition[0]=m_originalPointSet[i*3];
		localPosition[1]=m_originalPointSet[i*3+1];
		localPosition[2]=m_originalPointSet[i*3+2];
		localEXOne=Vector3Vector(localPosition,vectorOne);
		localEXTwo=Vector3Vector(localPosition,vectorTwo);
		coffecientMatrix[0][0]+=localEXOne*localEXOne;
		coffecientMatrix[0][1]+=localEXOne*localEXTwo;
		coffecientMatrix[0][2]+=localEXOne;
		coffecientMatrix[1][1]+=localEXTwo*localEXTwo;
		coffecientMatrix[1][2]+=localEXTwo;
		coffecientConstantRed[0]=m_originalColors[i*3]*localEXOne;
		coffecientConstantRed[1]=m_originalColors[i*3]*localEXTwo;
		coffecientConstantRed[2]=m_originalColors[i*3];
		coffecientConstantGreen[0]=m_originalColors[i*3+1]*localEXOne;
		coffecientConstantGreen[1]=m_originalColors[i*3+1]*localEXTwo;
		coffecientConstantGreen[2]=m_originalColors[i*3+1];
		coffecientConstantBlue[0]=m_originalColors[i*3+2]*localEXOne;
		coffecientConstantBlue[1]=m_originalColors[i*3+2]*localEXTwo;
		coffecientConstantBlue[2]=m_originalColors[i*3+2];
		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
			
			localPosition[0]=m_originalPointSet[(*vectorIterator)*3];
			localPosition[1]=m_originalPointSet[(*vectorIterator)*3+1];
			localPosition[2]=m_originalPointSet[(*vectorIterator)*3+2];
			localEXOne=Vector3Vector(localPosition,vectorOne);
			localEXTwo=Vector3Vector(localPosition,vectorTwo);
			coffecientMatrix[0][0]+=localEXOne*localEXOne;
			coffecientMatrix[0][1]+=localEXOne*localEXTwo;
			coffecientMatrix[0][2]+=localEXOne;
			coffecientMatrix[1][1]+=localEXTwo*localEXTwo;
			coffecientMatrix[1][2]+=localEXTwo;
			coffecientConstantRed[0]+=m_originalColors[(*vectorIterator)*3]*localEXOne;
			coffecientConstantRed[1]+=m_originalColors[(*vectorIterator)*3]*localEXTwo;
			coffecientConstantRed[2]+=m_originalColors[(*vectorIterator)*3];
			coffecientConstantGreen[0]+=m_originalColors[(*vectorIterator)*3+1]*localEXOne;
			coffecientConstantGreen[1]+=m_originalColors[(*vectorIterator)*3+1]*localEXTwo;
			coffecientConstantGreen[2]+=m_originalColors[(*vectorIterator)*3+1];
			coffecientConstantBlue[0]+=m_originalColors[(*vectorIterator)*3+2]*localEXOne;
			coffecientConstantBlue[1]+=m_originalColors[(*vectorIterator)*3+2]*localEXTwo;
			coffecientConstantBlue[2]+=m_originalColors[(*vectorIterator)*3+2];
		}

		//求方程的解 首先求矩阵的拟,然后求解
		mat_f8 m_matrix(3,3);
		vec_f8 m_vectorRed(3),m_vectorGreen(3),m_vectorBlue(3);
		vec_f8 m_gradientRed(3),m_gradientGreen(3),m_gradientBlue(3);
		coffecientMatrix[1][0]=coffecientMatrix[0][1];
		coffecientMatrix[2][0]=coffecientMatrix[0][2];
		coffecientMatrix[2][1]=coffecientMatrix[1][2];
		coffecientMatrix[2][2]=(*mapIterator).second.m_numOfNearest+1;
		for(int j=0;j<3;j++){
			for(int k=0;k<3;k++){
				m_matrix(j,k)=coffecientMatrix[j][k];
			}
			m_vectorRed(j)=coffecientConstantRed[j];
			m_vectorGreen(j)=coffecientConstantGreen[j];
			m_vectorBlue(j)=coffecientConstantBlue[j];
		}
		if(matrix_inverse(m_matrix)){
			m_gradientRed=m_matrix*m_vectorRed;
			m_gradientGreen=m_matrix*m_vectorGreen;
			m_gradientBlue=m_matrix*m_vectorBlue;
		}
		for(int j=0;j<3;j++){
			m_originalColorGradients[0][i*4+j]=m_gradientRed(0)*vectorOne[j]+m_gradientRed(1)*vectorTwo[j];
			m_originalColorGradients[1][i*4+j]=m_gradientGreen(0)*vectorOne[j]+m_gradientGreen(1)*vectorTwo[j];
			m_originalColorGradients[2][i*4+j]=m_gradientBlue(0)*vectorOne[j]+m_gradientBlue(1)*vectorTwo[j];
		}
		m_originalColorGradients[0][i*4+3]=m_gradientRed(2);
		m_originalColorGradients[1][i*4+3]=m_gradientGreen(2);
		m_originalColorGradients[2][i*4+3]=m_gradientBlue(2);
		
	}

	//////////////////////////////////////////////////////////////////////////
	// 测试代码
	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\GradientRed.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f %f\n",m_originalColorGradients[0][i*4],m_originalColorGradients[0][i*4+1],m_originalColorGradients[0][i*4+2],m_originalColorGradients[0][i*4+3]);

		}
		fclose(fpout);
	}

	 filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\GradientBlue.txt";

	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f %f\n",m_originalColorGradients[2][i*4],m_originalColorGradients[2][i*4+1],m_originalColorGradients[2][i*4+2],m_originalColorGradients[2][i*4+3]);

		}
		fclose(fpout);
	}

	 filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\GradientGreen.txt";
	
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f %f\n",m_originalColorGradients[1][i*4],m_originalColorGradients[1][i*4+1],m_originalColorGradients[1][i*4+2],m_originalColorGradients[1][i*4+3]);

		}
		fclose(fpout);
	}


	//////////////////////////////////////////////////////////////////////////
	

	return;
}

void FilterBilateral::CalculateBilateralFilterOne()
{
	//有shift的线性假定
	m_resultPointSet=new float[m_numOfPoints*3];
	m_resultColors=new float[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints;i++){
		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]=0;
			m_resultColors[i*3+j]=0;
		}
	}
	float sumKernel_position;
	float sumKernel_colors[3];
	for(int i=0;i<m_numOfPoints;i++){
		std::map<int ,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		
		sumKernel_position=0;
		sumKernel_colors[0]=0;
		sumKernel_colors[1]=0;
		sumKernel_colors[2]=0;
		ComputeOneOrderPositionEstimation(i);
		CalculateVariation();
		std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();
		std::vector<float>::iterator localProjectionDistanceIterator=m_localProjectionDistance.begin();
		std::vector<POINTVECTOR3D>::iterator localPositionEstimationIterator=m_localPositionEstimation.begin();
		std::vector<Color>::iterator localColorEstimation=m_localColorEstimation.begin();
		std::vector<float>::iterator localColorDistance[3];
		std::vector<float>::iterator localColorProjectionDistance[3];
		for(int j=0;j<3;j++){
			localColorDistance[j]=m_localColorDistance[j].begin();
			localColorProjectionDistance[j]=m_localColorProjectionDistance[j].begin();
		}
		for(;localPointDistanceIterator!=m_localPointDistance.end();localPointDistanceIterator++){
			float kernel_distance;
			float kernel_color[3];
			/*for(int j=0;j<3;j++){
				m_colorVariation[j]=0.1;
				m_colorEstimationVariation[j]=.1;
			}*/
			kernel_distance=exp(-0.5*(*localPointDistanceIterator)/m_pointVariation)*exp(-0.5*(*localProjectionDistanceIterator)/m_estimationVariation);
			if(((*localColorDistance[0])>m_thresholdColor[0])||((*localColorProjectionDistance[0])>m_thresholdColor[0])){
				kernel_color[0]=0;				
			}
			else{
				kernel_color[0]=kernel_distance
					*exp(-0.5*(*localColorDistance[0])/m_colorVariation[0])
					*exp(-0.5*(*localColorProjectionDistance[0])/m_colorEstimationVariation[0]);
			}
			if(((*localColorDistance[1])>m_thresholdColor[1])||((*localColorProjectionDistance[1])>m_thresholdColor[1])){
				kernel_color[1]=0;				
			}
			else{
				kernel_color[1]=kernel_distance
					*exp(-0.5*(*localColorDistance[1])/m_colorVariation[1])
					*exp(-0.5*(*localColorProjectionDistance[1])/m_colorEstimationVariation[1]);
			}
			if(((*localColorDistance[2])>m_thresholdColor[2])||((*localColorProjectionDistance[2])>m_thresholdColor[2])){
				kernel_color[2]=0;				
			}
			else{
				kernel_color[2]=kernel_distance
					*exp(-0.5*(*localColorDistance[2])/m_colorVariation[2])
					*exp(-0.5*(*localColorProjectionDistance[2])/m_colorEstimationVariation[2]);
			}
			

		
			sumKernel_position+=kernel_distance;
			sumKernel_colors[0]+=kernel_color[0];
			sumKernel_colors[1]+=kernel_color[1];
			sumKernel_colors[2]+=kernel_color[2];
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]+=(*localPositionEstimationIterator).pointVector[j]*kernel_distance;
				m_resultColors[i*3+j]+=(*localColorEstimation).pointColor[j]*kernel_color[j];
			}
			localProjectionDistanceIterator++;
		    localPositionEstimationIterator++;
			localColorEstimation++;
			for(int j=0;j<3;j++){
				localColorDistance[j]++;
				localColorProjectionDistance[j]++;			
			}
		}

		sumKernel_position+=1;
		sumKernel_colors[0]+=1;
		sumKernel_colors[1]+=1;
		sumKernel_colors[2]+=1;
		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]+=m_originalPointSet[i*3+j];
			m_resultPointSet[i*3+j]/=sumKernel_position;
			m_resultColors[i*3+j]+=m_originalColors[i*3+j];
			m_resultColors[i*3+j]/=sumKernel_colors[j];
		}
		m_localPointDistance.clear();
        m_localProjectionDistance.clear();
        m_localPositionEstimation.clear();
        m_localColorEstimation.clear();
		for(int j=0;j<3;j++){
			m_localColorDistance[j].clear();
		    m_localColorProjectionDistance[j].clear();	
		}	
	}
	delete[] m_originalPointSet;
	m_originalPointSet=NULL;
	delete[] m_originalColors;
	m_originalColors=NULL;

    //////////////////////////////////////////////////////////////////////////
    // 测试代码
	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\resultColors.txt";
	FILE *fpout;
	float averageColor[3];
	float variationColor[3];
	for(int j=0;j<3;j++){
		averageColor[j]=0;
		variationColor[j]=0;
	}
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f\n",m_resultColors[i*3],m_resultColors[i*3+1],m_resultColors[i*3+2]);
			averageColor[0]+=m_resultColors[i*3];
			averageColor[1]+=m_resultColors[i*3+1];
			averageColor[2]+=m_resultColors[i*3+2];
		}
		for(int i=0;i<3;i++){
			averageColor[i]/=m_numOfPoints;
		}
		fprintf(fpout,"%f %f %f\n",averageColor[0],averageColor[1],averageColor[2]);
		for(int i=0;i<m_numOfPoints;i++){
			variationColor[0]+=(m_resultColors[i*3]-averageColor[0])*(m_resultColors[i*3]-averageColor[0]);
			variationColor[1]+=(m_resultColors[i*3+1]-averageColor[1])*(m_resultColors[i*3+1]-averageColor[1]);
			variationColor[2]+=(m_resultColors[i*3+2]-averageColor[2])*(m_resultColors[i*3+2]-averageColor[2]);
		}
		for(int i=0;i<3;i++){
			variationColor[i]/=(m_numOfPoints-1);
		}
		fprintf(fpout,"%f %f %f\n",variationColor[0],variationColor[1],variationColor[2]);
		fclose(fpout);
	}
	//////////////////////////////////////////////////////////////////////////
	
}

void FilterBilateral::CalculateBilateralFilterTwo()
{
		//常数假定但加入了一个梯度聚类
	m_resultPointSet=new float[m_numOfPoints*3];
	m_resultColors=new float[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints;i++){
		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]=0;
			m_resultColors[i*3+j]=0;
		}
	}
	float sumKernel_position;
	float sumKernel_colors[3];
	for(int i=0;i<m_numOfPoints;i++){
		std::map<int ,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		
		sumKernel_position=0;
		sumKernel_colors[0]=0;
		sumKernel_colors[1]=0;
		sumKernel_colors[2]=0;
		ComputeOneOrderPositionEstimation(i);
		CalculateVariation();
		std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();
		std::vector<float>::iterator localProjectionDistanceIterator=m_localProjectionDistance.begin();
		std::vector<POINTVECTOR3D>::iterator localPositionEstimationIterator=m_localPositionEstimation.begin();
		std::vector<Color>::iterator localColorEstimation=m_localColorEstimation.begin();
		std::vector<float>::iterator localColorDistance[3];
		std::vector<float>::iterator localColorProjectionDistance[3];
		for(int j=0;j<3;j++){
			localColorDistance[j]=m_localColorDistance[j].begin();
			localColorProjectionDistance[j]=m_localColorProjectionDistance[j].begin();
		}
		for(;localPointDistanceIterator!=m_localPointDistance.end();localPointDistanceIterator++){
			float kernel_distance;
			float kernel_color[3];
			kernel_distance=exp(-0.5*(*localPointDistanceIterator)/m_pointVariation)*exp(-0.5*(*localProjectionDistanceIterator)/m_estimationVariation);
			if(((*localColorProjectionDistance[0])>m_thresholdColor[0])||((*localColorProjectionDistance[0])>m_thresholdColor[0])){
				kernel_color[0]=0;
			}
			else{
				kernel_color[0]=kernel_distance
					*exp(-0.5*(*localColorDistance[0])/m_colorVariation[0])
					*exp(-0.5*(*localColorProjectionDistance[0])/m_colorEstimationVariation[0]);
			}

			if(((*localColorProjectionDistance[1])>m_thresholdColor[1])||((*localColorProjectionDistance[1])>m_thresholdColor[1])){
				kernel_color[1]=0;
			}
			else{
				kernel_color[1]=kernel_distance
					*exp(-0.5*(*localColorDistance[1])/m_colorVariation[1])
					*exp(-0.5*(*localColorProjectionDistance[1])/m_colorEstimationVariation[1]);
			}
			if(((*localColorProjectionDistance[2])>m_thresholdColor[2])||((*localColorProjectionDistance[2])>m_thresholdColor[2])){
				kernel_color[2]=0;
			}
			else{
				kernel_color[2]=kernel_distance
					*exp(-0.5*(*localColorDistance[2])/m_colorVariation[2])
					*exp(-0.5*(*localColorProjectionDistance[2])/m_colorEstimationVariation[2]);
			}
		
			sumKernel_position+=kernel_distance;
			sumKernel_colors[0]+=kernel_color[0];
			sumKernel_colors[1]+=kernel_color[1];
			sumKernel_colors[2]+=kernel_color[2];
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]+=(*localPositionEstimationIterator).pointVector[j]*kernel_distance;
				m_resultColors[i*3+j]+=m_originalColors[(*vectorIterator)*3+j]*kernel_color[j];
			}
			localProjectionDistanceIterator++;
		    localPositionEstimationIterator++;
			localColorEstimation++;
			vectorIterator++;
			for(int j=0;j<3;j++){
				localColorDistance[j]++;
				localColorProjectionDistance[j]++;			
			}
		}

		sumKernel_position+=1;
		sumKernel_colors[0]+=1;
		sumKernel_colors[1]+=1;
		sumKernel_colors[2]+=1;
		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]+=m_originalPointSet[i*3+j];
			m_resultPointSet[i*3+j]/=sumKernel_position;
			m_resultColors[i*3+j]+=m_originalColors[i*3+j];
			m_resultColors[i*3+j]/=sumKernel_colors[j];
		}
		m_localPointDistance.clear();
        m_localProjectionDistance.clear();
        m_localPositionEstimation.clear();
        m_localColorEstimation.clear();
		for(int j=0;j<3;j++){
			m_localColorDistance[j].clear();
		    m_localColorProjectionDistance[j].clear();	
		}	
	}
	delete[] m_originalPointSet;
	delete[] m_originalColors;
	m_originalPointSet=NULL;
	m_originalColors=NULL;


	//////////////////////////////////////////////////////////////////////////
	// 测试代码
	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\resultColorsTwo.txt";

	FILE *fpout;
	float averageColor[3];
	float variationColor[3];
	for(int j=0;j<3;j++){
		averageColor[j]=0;
		variationColor[j]=0;
	}
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f\n",m_resultColors[i*3],m_resultColors[i*3+1],m_resultColors[i*3+2]);
			averageColor[0]+=m_resultColors[i*3];
			averageColor[1]+=m_resultColors[i*3+1];
			averageColor[2]+=m_resultColors[i*3+2];
		}
		for(int i=0;i<3;i++){
			averageColor[i]/=m_numOfPoints;
		}
		fprintf(fpout,"%f %f %f\n",averageColor[0],averageColor[1],averageColor[2]);
		for(int i=0;i<m_numOfPoints;i++){
			variationColor[0]+=(m_resultColors[i*3]-averageColor[0])*(m_resultColors[i*3]-averageColor[0]);
			variationColor[1]+=(m_resultColors[i*3+1]-averageColor[1])*(m_resultColors[i*3+1]-averageColor[1]);
			variationColor[2]+=(m_resultColors[i*3+2]-averageColor[2])*(m_resultColors[i*3+2]-averageColor[2]);
		}
		for(int i=0;i<3;i++){
			variationColor[i]/=(m_numOfPoints-1);
		}
		fprintf(fpout,"%f %f %f\n",variationColor[0],variationColor[1],variationColor[2]);
		fclose(fpout);
	}
	//////////////////////////////////////////////////////////////////////////


}
//////////////////////////////////////////////////////////////////////////
// 各项同性滤波 （颜色的各项同性）仅仅考虑 原始点的距离
void FilterBilateral::CalculateBilateralFilterThree()
{
	m_resultPointSet=new float[m_numOfPoints*3];
	m_resultColors=new float[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints;i++){
		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]=0;
			m_resultColors[i*3+j]=0;
		}
	}
	float sumKernel_position;
	float sumKernel_colors[3];
	for(int i=0;i<m_numOfPoints;i++){
		std::map<int ,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();

		sumKernel_position=0;
		sumKernel_colors[0]=0;
		sumKernel_colors[1]=0;
		sumKernel_colors[2]=0;
		ComputeOneOrderPositionEstimation(i);
		CalculateVariation();
		std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();
		std::vector<float>::iterator localProjectionDistanceIterator=m_localProjectionDistance.begin();
		std::vector<POINTVECTOR3D>::iterator localPositionEstimationIterator=m_localPositionEstimation.begin();
	/*	std::vector<Color>::iterator localColorEstimation=m_localColorEstimation.begin();
		std::vector<float>::iterator localColorDistance[3];
		std::vector<float>::iterator localColorProjectionDistance[3];
		for(int j=0;j<3;j++){
			localColorDistance[j]=m_localColorDistance[j].begin();
			localColorProjectionDistance[j]=m_localColorProjectionDistance[j].begin();
		}*/
		for(;localPointDistanceIterator!=m_localPointDistance.end();localPointDistanceIterator++){
			float kernel_distance;
			float kernel_color[3];
			kernel_distance=exp(-0.5*(*localPointDistanceIterator)/m_pointVariation)*exp(-0.5*(*localProjectionDistanceIterator)/m_estimationVariation);
			kernel_color[0]=kernel_color[1]=kernel_color[2]=exp(-0.5*(*localPointDistanceIterator)/m_pointVariation);
			/*kernel_color[0]=kernel_distance
				*exp(-0.5*(*localColorDistance[0])/m_colorVariation[0])
				*exp(-0.5*(*localColorProjectionDistance[0])/m_colorEstimationVariation[0]);
			kernel_color[1]=kernel_distance
				*exp(-0.5*(*localColorDistance[1])/m_colorVariation[1])
				*exp(-0.5*(*localColorProjectionDistance[1])/m_colorEstimationVariation[1]);
			kernel_color[2]=kernel_distance
				*exp(-0.5*(*localColorDistance[2])/m_colorVariation[2])
				*exp(-0.5*(*localColorProjectionDistance[2])/m_colorEstimationVariation[2]);*/
			sumKernel_position+=kernel_distance;
			sumKernel_colors[0]+=kernel_color[0];
			sumKernel_colors[1]+=kernel_color[1];
			sumKernel_colors[2]+=kernel_color[2];
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]+=(*localPositionEstimationIterator).pointVector[j]*kernel_distance;
				m_resultColors[i*3+j]+=m_originalColors[(*vectorIterator)*3+j]*kernel_color[j];
			}
			localProjectionDistanceIterator++;
			localPositionEstimationIterator++;
			//localColorEstimation++;
			vectorIterator++;
			/*for(int j=0;j<3;j++){
				localColorDistance[j]++;
				localColorProjectionDistance[j]++;			
			}*/
		}

		sumKernel_position+=1;
		sumKernel_colors[0]+=1;
		sumKernel_colors[1]+=1;
		sumKernel_colors[2]+=1;
		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]+=m_originalPointSet[i*3+j];
			m_resultPointSet[i*3+j]/=sumKernel_position;
			m_resultColors[i*3+j]+=m_originalColors[i*3+j];
			m_resultColors[i*3+j]/=sumKernel_colors[j];
		}
		m_localPointDistance.clear();
		m_localProjectionDistance.clear();
		m_localPositionEstimation.clear();
		m_localColorEstimation.clear();
		for(int j=0;j<3;j++){
			m_localColorDistance[j].clear();
			m_localColorProjectionDistance[j].clear();	
		}	
	}
	delete[] m_originalPointSet;
	delete[] m_originalColors;
	m_originalPointSet=NULL;
	m_originalColors=NULL;


	//////////////////////////////////////////////////////////////////////////
	// 测试代码
	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\resultColorsThree.txt";

	FILE *fpout;
	float averageColor[3];
	float variationColor[3];
	for(int j=0;j<3;j++){
		averageColor[j]=0;
		variationColor[j]=0;
	}
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f\n",m_resultColors[i*3],m_resultColors[i*3+1],m_resultColors[i*3+2]);
			averageColor[0]+=m_resultColors[i*3];
			averageColor[1]+=m_resultColors[i*3+1];
			averageColor[2]+=m_resultColors[i*3+2];
		}
		for(int i=0;i<3;i++){
			averageColor[i]/=m_numOfPoints;
		}
		fprintf(fpout,"%f %f %f\n",averageColor[0],averageColor[1],averageColor[2]);
		for(int i=0;i<m_numOfPoints;i++){
			variationColor[0]+=(m_resultColors[i*3]-averageColor[0])*(m_resultColors[i*3]-averageColor[0]);
			variationColor[1]+=(m_resultColors[i*3+1]-averageColor[1])*(m_resultColors[i*3+1]-averageColor[1]);
			variationColor[2]+=(m_resultColors[i*3+2]-averageColor[2])*(m_resultColors[i*3+2]-averageColor[2]);
		}
		for(int i=0;i<3;i++){
			variationColor[i]/=(m_numOfPoints-1);
		}
		fprintf(fpout,"%f %f %f\n",variationColor[0],variationColor[1],variationColor[2]);
		fclose(fpout);
	}
	//////////////////////////////////////////////////////////////////////////


}

void FilterBilateral::ComputeOneOrderPositionEstimation(int indexPoint)
{
	if(indexPoint==1062){
		int k;
		k=0;
	}
	
	std::map<int,KnearestField>::iterator mapIterator=m_mapKnearest.find(indexPoint);
	std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	float localNormal[3];//局部法向
	localNormal[0]=m_filterNormals[indexPoint*3];
	localNormal[1]=m_filterNormals[indexPoint*3+1];
	localNormal[2]=m_filterNormals[indexPoint*3+2];
	for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
	   float distance_point_point;//点到点的距离
	   float distance_point_projection;//点到一次估计平面的距离
	   float localVector_point_point[3];//点到点的向量
	   float localColorDistance[3];//颜色距离
	   
	   POINTVECTOR3D localProjection;//投影点的位置

	   localVector_point_point[0]=m_originalPointSet[(*vectorIterator)*3]-m_originalPointSet[indexPoint*3];
	   localVector_point_point[1]=m_originalPointSet[(*vectorIterator)*3+1]-m_originalPointSet[indexPoint*3+1];
	   localVector_point_point[2]=m_originalPointSet[(*vectorIterator)*3+2]-m_originalPointSet[indexPoint*3+2];

	   distance_point_projection=Vector3Vector(localVector_point_point,localNormal);

	   distance_point_point=(Vector3Vector(localVector_point_point,localVector_point_point));

	   localProjection.pointVector[0]=m_originalPointSet[indexPoint*3]+distance_point_projection*localNormal[0];
	   localProjection.pointVector[1]=m_originalPointSet[indexPoint*3+1]+distance_point_projection*localNormal[1];
	   localProjection.pointVector[2]=m_originalPointSet[indexPoint*3+2]+distance_point_projection*localNormal[2];

	   localColorDistance[0]=((m_originalColors[indexPoint*3]-m_originalColors[(*vectorIterator)*3])*(m_originalColors[indexPoint*3]-m_originalColors[(*vectorIterator)*3]));
       localColorDistance[1]=((m_originalColors[indexPoint*3+1]-m_originalColors[(*vectorIterator)*3+1])*(m_originalColors[indexPoint*3+1]-m_originalColors[(*vectorIterator)*3+1]));
       localColorDistance[2]=((m_originalColors[indexPoint*3+2]-m_originalColors[(*vectorIterator)*3+2])*(m_originalColors[indexPoint*3+2]-m_originalColors[(*vectorIterator)*3+2]));
       
	   m_localProjectionDistance.push_back((distance_point_projection)*distance_point_projection);
	   m_localPointDistance.push_back(distance_point_point);
	   m_localPositionEstimation.push_back(localProjection);
	   m_localColorDistance[0].push_back(localColorDistance[0]);
	   m_localColorDistance[1].push_back(localColorDistance[1]);
       m_localColorDistance[2].push_back(localColorDistance[2]);

	}
	//计算一次颜色估计和颜色投影距离
	ComputeOneOrderColorEstimation(indexPoint);
    return;
}

void FilterBilateral::ComputeOneOrderColorEstimation(int indexPoint)
{

	std::map<int,KnearestField>::iterator mapIterator=m_mapKnearest.find(indexPoint);
	std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	std::vector<POINTVECTOR3D>::iterator positionEstimatorIterator=m_localPositionEstimation.begin();
	for(;positionEstimatorIterator!=m_localPositionEstimation.end();positionEstimatorIterator++,vectorIterator++){
		Color localColorEstimation;
		float localColorEstimationDistance[3];

		for(int i=0;i<3;i++){
			localColorEstimation.pointColor[i]=m_originalColorGradients[i][(*vectorIterator)*4]*(*positionEstimatorIterator).pointVector[0]
			+m_originalColorGradients[i][(*vectorIterator)*4+1]*(*positionEstimatorIterator).pointVector[1]
			+m_originalColorGradients[i][(*vectorIterator)*4+2]*(*positionEstimatorIterator).pointVector[2]
			+m_originalColorGradients[i][(*vectorIterator)*4+3];
			localColorEstimationDistance[i]=((localColorEstimation.pointColor[i]-m_originalColors[indexPoint*3+i])*(localColorEstimation.pointColor[i]-m_originalColors[indexPoint*3+i]));
			m_localColorProjectionDistance[i].push_back(localColorEstimationDistance[i]);
		}
		m_localColorEstimation.push_back(localColorEstimation);
	}
	//////////////////////////////////////////////////////////////////////////
	//测试代码
	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\colorEstimation.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		std::vector<Color>::iterator vectorIterator=m_localColorEstimation.begin();
		for(;vectorIterator!=m_localColorEstimation.end();vectorIterator++){
			fprintf(fpout,"%f %f %f\n",(*vectorIterator).pointColor[0],(*vectorIterator).pointColor[1],(*vectorIterator).pointColor[2]);
		}		
		fclose(fpout);
	}
	return;
}

void FilterBilateral::CalculateVariation()
{
	// 
	float numOfColorOff[3],numOfColorEstimationOff[3];//
	for(int i=0;i<3;i++){
		numOfColorOff[i]=0;
		numOfColorEstimationOff[i]=0;

	}
	
	std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();
	std::vector<float>::iterator localProjectionDistanceIterator=m_localProjectionDistance.begin();
	std::vector<float>::iterator localColorDistanceIterator[3];
	std::vector<float>::iterator localColorProjectionDistanceIterator[3];
	for(int i=0;i<3;i++){
		localColorDistanceIterator[i]=m_localColorDistance[i].begin();
		localColorProjectionDistanceIterator[i]=m_localColorProjectionDistance[i].begin();
	}
	m_pointVariation=0;
	m_estimationVariation=0;
	for(int i=0;i<3;i++){
		m_colorVariation[i]=0;
		m_colorEstimationVariation[i]=0;
	}
	for(;localPointDistanceIterator!=m_localPointDistance.end();localPointDistanceIterator++,localProjectionDistanceIterator++){
		m_pointVariation+=(*localPointDistanceIterator);
	    m_estimationVariation+=(*localProjectionDistanceIterator);
		for(int i=0;i<3;i++){
			if((*localColorDistanceIterator[i])<m_thresholdColor[i]){
				m_colorVariation[i]+=(*localColorDistanceIterator[i]);
			}
			else{
				numOfColorOff[i]+=1;
			}
			if((*localColorProjectionDistanceIterator[i])<m_thresholdColor[i]){
                	m_colorEstimationVariation[i]+=(*localColorProjectionDistanceIterator[i]);
				
			}
			else{
				numOfColorEstimationOff[i]+=1;
			}
			/*m_colorEstimationVariation[i]+=(*localColorProjectionDistanceIterator[i]);
		    m_colorVariation[i]+=(*localColorDistanceIterator[i]);*/
			localColorDistanceIterator[i]++;
			localColorProjectionDistanceIterator[i]++;
		}
	}
	m_pointVariation/=m_localPointDistance.size();
	m_estimationVariation/=m_localPointDistance.size();
	for(int i=0;i<3;i++){
		if(m_localPointDistance.size()==numOfColorOff[i]){
			m_colorVariation[i]=0.01;
		}
		else{
			m_colorVariation[i]/=(m_localPointDistance.size()-numOfColorOff[i]);
			if(m_colorVariation[i]==0)
				m_colorVariation[i]=0.01;
		}		
		

		if(m_localPointDistance.size()==numOfColorEstimationOff[i]){
			m_colorEstimationVariation[i]=0.01;
		}
		else{
			m_colorEstimationVariation[i]/=(m_localPointDistance.size()-numOfColorEstimationOff[i]);
			if(m_colorEstimationVariation[i]==0)
				m_colorEstimationVariation[i]=0.01;
		}
	
	}
	return;
}

void FilterBilateral::AddPositionGaussianNoise()
{
	////求邻近点的最小值

	////加入高斯噪声
	////////////////////////////////////////////////////////////////////////////
	//// 首先加入 均值为1的 
	////m_resultPointSet=new float[m_numOfPoints*3];
	//float noise[3];
	//long idum[3];
	//idum[0]=1;
	//idum[1]=1;
	//idum[2]=1;

	//for(int i=0;i<m_numOfPoints;i++){
	//	for(int j=0;j<3;j++){
	//		noise[j]=gasdev(&idum[j]);
	//		//m_resultPointSet[i*3+j]=m_originalPointSet[i*3+j]+noise[j]*0.05;
	//		m_originalPointSet[i*3+j]=m_originalPointSet[i*3+j]+noise[j]*0.05;

	//	//	m_resultPointSet[i*3+j]=m_originalPointSet[i*3+j];//+noise[j]*0.01;
	//	}        
	//}
	//return;
}

void FilterBilateral::AddColorGaussianNorse()
{
	////把颜色转化为0-1之间的浮点数
	///*for(int i=0;i<m_numOfPoints*3;i++){
	//m_originalFloatColor[i]=m_originalColor[i]/255;
	//}
	//delete [] m_originalColor;*/
	////加入高斯白噪声

	////m_resultColors=new float[m_numOfPoints*3];

	//float noise[3];
	//long idum[3];
	//idum[0]=2;
	//idum[1]=2;
	//idum[2]=2;

	//for(int i=0;i<m_numOfPoints;i++){
	//	for(int j=0;j<3;j++){
	//		noise[j]=gasdev(&idum[0]);
	//		m_originalColors[i*3+j]=m_originalColors[i*3+j]+noise[j]*0.05;

	//		/*if(m_originalColor[i*3+j]<0){
	//		m_originalColor[i*3+j]=0;
	//		}
	//		else
	//		if(m_originalColor[i*3+j]>1){
	//		m_originalColor[i*3+j]=1;

	//		}*/
	//	//	m_resultColors[i*3+j]=m_originalColors[i*3+j];
	//	}        
	//}
	////////////////////////////////////////////////////////////////////////////
	//// 测试代码
	//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\originalColors.txt";
	//FILE *fpout;
	//float averageColor[3],variationColor[3];
	//for(int j=0;j<3;j++){
	//	averageColor[j]=0;
	//	variationColor[j]=0;
	//}
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f\n",m_originalColors[i*3],m_originalColors[i*3+1],m_originalColors[i*3+2]);
	//		averageColor[0]+=m_originalColors[i*3];
	//		averageColor[1]+=m_originalColors[i*3+1];
	//		averageColor[2]+=m_originalColors[i*3+2];
	//	}
	//	for(int i=0;i<3;i++){
	//		averageColor[i]/=m_numOfPoints;
	//	}
	//	fprintf(fpout,"%f %f %f\n",averageColor[0],averageColor[1],averageColor[2]);
	//	for(int i=0;i<m_numOfPoints;i++){
	//		variationColor[0]+=(m_originalColors[i*3]-averageColor[0])*(m_originalColors[i*3]-averageColor[0]);
	//		variationColor[1]+=(m_originalColors[i*3+1]-averageColor[1])*(m_originalColors[i*3+1]-averageColor[1]);
	//		variationColor[2]+=(m_originalColors[i*3+2]-averageColor[2])*(m_originalColors[i*3+2]-averageColor[2]);
	//	}
	//	for(int i=0;i<3;i++){
	//		variationColor[i]/=(m_numOfPoints-1);
	//	}
	//	fprintf(fpout,"%f %f %f\n",variationColor[0],variationColor[1],variationColor[2]);
	//	
	//	fclose(fpout);
	//}
	//////////////////////////////////////////////////////////////////////////
}

//////////////////////////////////////////////////////////////////////////
// 各项同性滤波 （颜色的各项同性）考虑 原始点的距离 和投影点的距离
void FilterBilateral::CalculateBilateralFilterFour()
{
	m_resultPointSet=new float[m_numOfPoints*3];
	m_resultColors=new float[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints;i++){
		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]=0;
			m_resultColors[i*3+j]=0;
		}
	}
	float sumKernel_position;
	float sumKernel_colors[3];
	for(int i=0;i<m_numOfPoints;i++){
		std::map<int ,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();

		sumKernel_position=0;
		sumKernel_colors[0]=0;
		sumKernel_colors[1]=0;
		sumKernel_colors[2]=0;
		ComputeOneOrderPositionEstimation(i);
		CalculateVariation();
		std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();
		std::vector<float>::iterator localProjectionDistanceIterator=m_localProjectionDistance.begin();
		std::vector<POINTVECTOR3D>::iterator localPositionEstimationIterator=m_localPositionEstimation.begin();
		/*std::vector<Color>::iterator localColorEstimation=m_localColorEstimation.begin();
		std::vector<float>::iterator localColorDistance[3];
		std::vector<float>::iterator localColorProjectionDistance[3];
		for(int j=0;j<3;j++){
			localColorDistance[j]=m_localColorDistance[j].begin();
			localColorProjectionDistance[j]=m_localColorProjectionDistance[j].begin();
		}*/
		for(;localPointDistanceIterator!=m_localPointDistance.end();localPointDistanceIterator++){
			float kernel_distance;
			float kernel_color[3];
			kernel_distance=exp(-0.5*(*localPointDistanceIterator)/m_pointVariation)*exp(-0.5*(*localProjectionDistanceIterator)/m_estimationVariation);
			kernel_color[0]=kernel_color[1]=kernel_color[2]=kernel_distance;
			/*kernel_color[0]=kernel_distance
				*exp(-0.5*(*localColorDistance[0])/m_colorVariation[0])
				*exp(-0.5*(*localColorProjectionDistance[0])/m_colorEstimationVariation[0]);
			kernel_color[1]=kernel_distance
				*exp(-0.5*(*localColorDistance[1])/m_colorVariation[1])
				*exp(-0.5*(*localColorProjectionDistance[1])/m_colorEstimationVariation[1]);
			kernel_color[2]=kernel_distance
				*exp(-0.5*(*localColorDistance[2])/m_colorVariation[2])
				*exp(-0.5*(*localColorProjectionDistance[2])/m_colorEstimationVariation[2]);*/
			sumKernel_position+=kernel_distance;
			sumKernel_colors[0]+=kernel_color[0];
			sumKernel_colors[1]+=kernel_color[1];
			sumKernel_colors[2]+=kernel_color[2];
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]+=(*localPositionEstimationIterator).pointVector[j]*kernel_distance;
				m_resultColors[i*3+j]+=m_originalColors[(*vectorIterator)*3+j]*kernel_color[j];
			}
			localProjectionDistanceIterator++;
			localPositionEstimationIterator++;
			//localColorEstimation++;
			vectorIterator++;
			/*for(int j=0;j<3;j++){
				localColorDistance[j]++;
				localColorProjectionDistance[j]++;			
			}*/
		}

		sumKernel_position+=1;
		sumKernel_colors[0]+=1;
		sumKernel_colors[1]+=1;
		sumKernel_colors[2]+=1;
		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]+=m_originalPointSet[i*3+j];
			m_resultPointSet[i*3+j]/=sumKernel_position;
			m_resultColors[i*3+j]+=m_originalColors[i*3+j];
			m_resultColors[i*3+j]/=sumKernel_colors[j];
		}
		m_localPointDistance.clear();
		m_localProjectionDistance.clear();
		m_localPositionEstimation.clear();
		m_localColorEstimation.clear();
		for(int j=0;j<3;j++){
			m_localColorDistance[j].clear();
			m_localColorProjectionDistance[j].clear();	
		}	
	}
	delete[] m_originalPointSet;
	delete[] m_originalColors;
	m_originalPointSet=NULL;
	m_originalColors=NULL;


	//////////////////////////////////////////////////////////////////////////
	// 测试代码
	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\resultColorsFour.txt";

	FILE *fpout;
	float averageColor[3];
	float variationColor[3];
	for(int j=0;j<3;j++){
		averageColor[j]=0;
		variationColor[j]=0;
	}
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f\n",m_resultColors[i*3],m_resultColors[i*3+1],m_resultColors[i*3+2]);
			averageColor[0]+=m_resultColors[i*3];
			averageColor[1]+=m_resultColors[i*3+1];
			averageColor[2]+=m_resultColors[i*3+2];
		}
		for(int i=0;i<3;i++){
			averageColor[i]/=m_numOfPoints;
		}
		fprintf(fpout,"%f %f %f\n",averageColor[0],averageColor[1],averageColor[2]);
		for(int i=0;i<m_numOfPoints;i++){
			variationColor[0]+=(m_resultColors[i*3]-averageColor[0])*(m_resultColors[i*3]-averageColor[0]);
			variationColor[1]+=(m_resultColors[i*3+1]-averageColor[1])*(m_resultColors[i*3+1]-averageColor[1]);
			variationColor[2]+=(m_resultColors[i*3+2]-averageColor[2])*(m_resultColors[i*3+2]-averageColor[2]);
		}
		for(int i=0;i<3;i++){
			variationColor[i]/=(m_numOfPoints-1);
		}
		fprintf(fpout,"%f %f %f\n",variationColor[0],variationColor[1],variationColor[2]);
		fclose(fpout);
	}
	//////////////////////////////////////////////////////////////////////////
	
}
//////////////////////////////////////////////////////////////////////////
// 加入颜色信息 常数假定 但没有考虑梯度信息
void FilterBilateral::CalculateBilateralFilterFive()
{
	m_resultPointSet=new float[m_numOfPoints*3];
	m_resultColors=new float[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints;i++){
		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]=0;
			m_resultColors[i*3+j]=0;
		}
	}
	float sumKernel_position;
	float sumKernel_colors[3];
	for(int i=0;i<m_numOfPoints;i++){
		std::map<int ,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();

		sumKernel_position=0;
		sumKernel_colors[0]=0;
		sumKernel_colors[1]=0;
		sumKernel_colors[2]=0;
		
		ComputeOneOrderPositionEstimation(i);
		
		CalculateVariation();
		std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();
		std::vector<float>::iterator localProjectionDistanceIterator=m_localProjectionDistance.begin();
		std::vector<POINTVECTOR3D>::iterator localPositionEstimationIterator=m_localPositionEstimation.begin();
		std::vector<Color>::iterator localColorEstimation=m_localColorEstimation.begin();
		std::vector<float>::iterator localColorDistance[3];
		std::vector<float>::iterator localColorProjectionDistance[3];
		
		for(int j=0;j<3;j++){
			localColorDistance[j]=m_localColorDistance[j].begin();
			localColorProjectionDistance[j]=m_localColorProjectionDistance[j].begin();
		}
		for(;localPointDistanceIterator!=m_localPointDistance.end();localPointDistanceIterator++){
			float kernel_distance;
			float kernel_color[3];
			kernel_distance=exp(-0.5*(*localPointDistanceIterator)/m_pointVariation)*exp(-0.5*(*localProjectionDistanceIterator)/m_estimationVariation);
			if((*localColorDistance[0])>m_thresholdColor[0]){
				kernel_color[0]=0;
			}
			else{
				kernel_color[0]=kernel_distance
					*exp(-0.5*(*localColorDistance[0])/m_colorVariation[0]);
			}
			if((*localColorDistance[1])>m_thresholdColor[1]){
				kernel_color[1]=0;
			}
			else{
				kernel_color[1]=kernel_distance
					*exp(-0.5*(*localColorDistance[1])/m_colorVariation[1]);
			}
			if((*localColorDistance[2])>m_thresholdColor[2]){
				kernel_color[2]=0;
			}
			else{
				kernel_color[2]=kernel_distance
					*exp(-0.5*(*localColorDistance[2])/m_colorVariation[2]);
			}
			
			sumKernel_position+=kernel_distance;
			sumKernel_colors[0]+=kernel_color[0];
			sumKernel_colors[1]+=kernel_color[1];
			sumKernel_colors[2]+=kernel_color[2];
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]+=(*localPositionEstimationIterator).pointVector[j]*kernel_distance;
				m_resultColors[i*3+j]+=m_originalColors[(*vectorIterator)*3+j]*kernel_color[j];
			}
			localProjectionDistanceIterator++;
			localPositionEstimationIterator++;
			localColorEstimation++;
			vectorIterator++;
			for(int j=0;j<3;j++){
				localColorDistance[j]++;
				localColorProjectionDistance[j]++;			
			}
		}

		sumKernel_position+=1;
		sumKernel_colors[0]+=1;
		sumKernel_colors[1]+=1;
		sumKernel_colors[2]+=1;
		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]+=m_originalPointSet[i*3+j];
			m_resultPointSet[i*3+j]/=sumKernel_position;
			m_resultColors[i*3+j]+=m_originalColors[i*3+j];
			m_resultColors[i*3+j]/=sumKernel_colors[j];
			if(m_resultColors[i*3+j]<0||m_resultColors[i*3+j]>1){
				int w;
				w=0;
			}
			else{
				
			}
		}
		m_localPointDistance.clear();
		m_localProjectionDistance.clear();
		m_localPositionEstimation.clear();
		m_localColorEstimation.clear();
		for(int j=0;j<3;j++){
			m_localColorDistance[j].clear();
			m_localColorProjectionDistance[j].clear();	
		}	
	}
	delete[] m_originalPointSet;
	delete[] m_originalColors;
	m_originalPointSet=NULL;
	m_originalColors=NULL;


	//////////////////////////////////////////////////////////////////////////
	// 测试代码
	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\resultColorsFive.txt";

	FILE *fpout;
	float averageColor[3];
	float variationColor[3];
	for(int j=0;j<3;j++){
		averageColor[j]=0;
		variationColor[j]=0;
	}
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f\n",m_resultColors[i*3],m_resultColors[i*3+1],m_resultColors[i*3+2]);
			averageColor[0]+=m_resultColors[i*3];
			averageColor[1]+=m_resultColors[i*3+1];
			averageColor[2]+=m_resultColors[i*3+2];
		}
		for(int i=0;i<3;i++){
			averageColor[i]/=m_numOfPoints;
		}
		fprintf(fpout,"%f %f %f\n",averageColor[0],averageColor[1],averageColor[2]);
		for(int i=0;i<m_numOfPoints;i++){
			variationColor[0]+=(m_resultColors[i*3]-averageColor[0])*(m_resultColors[i*3]-averageColor[0]);
			variationColor[1]+=(m_resultColors[i*3+1]-averageColor[1])*(m_resultColors[i*3+1]-averageColor[1]);
			variationColor[2]+=(m_resultColors[i*3+2]-averageColor[2])*(m_resultColors[i*3+2]-averageColor[2]);
		}
		for(int i=0;i<3;i++){
			variationColor[i]/=(m_numOfPoints-1);
		}
		fprintf(fpout,"%f %f %f\n",variationColor[0],variationColor[1],variationColor[2]);
		fclose(fpout);
	}
	//////////////////////////////////////////////////////////////////////////
}
void FilterBilateral::CalculateBilateralFilterSix()
{
	
	m_pointVariation=0.2;
	m_resultPointSet=new float[m_numOfPoints*3];
	m_resultColors=new float[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints;i++){
		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]=0;
		}
	}
	float sumKernel_position;
	float sumKernel_colors[3];
	for(int i=0;i<m_numOfPoints;i++){
		std::map<int ,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();

		sumKernel_position=0;
		ComputeOneOrderPositionEstimation(i);
		std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();
		std::vector<POINTVECTOR3D>::iterator localPositionEstimationIterator=m_localPositionEstimation.begin();
	
		for(;localPointDistanceIterator!=m_localPointDistance.end();localPointDistanceIterator++){
			float kernel_distance;
			kernel_distance=exp(-0.5*(*localPointDistanceIterator)/m_pointVariation);
		
		
	

			sumKernel_position+=kernel_distance;

			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]+=(*localPositionEstimationIterator).pointVector[j]*kernel_distance;
				
			}
		
			localPositionEstimationIterator++;
			
		
		}

		sumKernel_position+=1;
	
		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]+=m_originalPointSet[i*3+j];
			m_resultPointSet[i*3+j]/=sumKernel_position;
		
		}
		m_localPointDistance.clear();
	
		m_localPositionEstimation.clear();
		m_localColorEstimation.clear();
	
	}
	delete[] m_originalPointSet;
	delete[] m_originalColors;
	m_originalPointSet=NULL;
	m_originalColors=NULL;


	//////////////////////////////////////////////////////////////////////////
	// 测试代码
	//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\resultColorsFive.txt";

	//FILE *fpout;
	//float averageColor[3];
	//float variationColor[3];
	//for(int j=0;j<3;j++){
	//	averageColor[j]=0;
	//	variationColor[j]=0;
	//}
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f\n",m_resultColors[i*3],m_resultColors[i*3+1],m_resultColors[i*3+2]);
	//		averageColor[0]+=m_resultColors[i*3];
	//		averageColor[1]+=m_resultColors[i*3+1];
	//		averageColor[2]+=m_resultColors[i*3+2];
	//	}
	//	for(int i=0;i<3;i++){
	//		averageColor[i]/=m_numOfPoints;
	//	}
	//	fprintf(fpout,"%f %f %f\n",averageColor[0],averageColor[1],averageColor[2]);
	//	for(int i=0;i<m_numOfPoints;i++){
	//		variationColor[0]+=(m_resultColors[i*3]-averageColor[0])*(m_resultColors[i*3]-averageColor[0]);
	//		variationColor[1]+=(m_resultColors[i*3+1]-averageColor[1])*(m_resultColors[i*3+1]-averageColor[1]);
	//		variationColor[2]+=(m_resultColors[i*3+2]-averageColor[2])*(m_resultColors[i*3+2]-averageColor[2]);
	//	}
	//	for(int i=0;i<3;i++){
	//		variationColor[i]/=(m_numOfPoints-1);
	//	}
	//	fprintf(fpout,"%f %f %f\n",variationColor[0],variationColor[1],variationColor[2]);
	//	fclose(fpout);
	//}
	//////////////////////////////////////////////////////////////////////////
}
void FilterBilateral::CalculateBilateralFilterSeven()
{
	m_resultPointSet=new float[m_numOfPoints*3];
	m_resultColors=new float[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints;i++){
		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]=0;
			m_resultColors[i*3+j]=0;
		}
	}
	float sumKernel_position;
	float sumKernel_colors[3];
	for(int i=0;i<m_numOfPoints;i++){
		std::map<int ,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();

		sumKernel_position=0;
		sumKernel_colors[0]=0;
		sumKernel_colors[1]=0;
		sumKernel_colors[2]=0;
		if(i==6513){
			int k;
			k=0;
		}
		ComputeOneOrderPositionEstimation(i);

		CalculateVariation();
		std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();
		std::vector<float>::iterator localProjectionDistanceIterator=m_localProjectionDistance.begin();
		std::vector<POINTVECTOR3D>::iterator localPositionEstimationIterator=m_localPositionEstimation.begin();
		std::vector<Color>::iterator localColorEstimation=m_localColorEstimation.begin();
		std::vector<float>::iterator localColorDistance[3];
		std::vector<float>::iterator localColorProjectionDistance[3];

		for(int j=0;j<3;j++){
			localColorDistance[j]=m_localColorDistance[j].begin();
			localColorProjectionDistance[j]=m_localColorProjectionDistance[j].begin();
		}
		for(;localPointDistanceIterator!=m_localPointDistance.end();localPointDistanceIterator++){
			float kernel_distance;
			float kernel_color[3];
			kernel_distance=exp(-0.5*(*localPointDistanceIterator)/m_pointVariation)*exp(-0.5*(*localProjectionDistanceIterator)/m_estimationVariation);
			if((*localColorDistance[0])>m_thresholdColor[0]){
				kernel_color[0]=0;
			}
			else{
				kernel_color[0]=kernel_distance
					*exp(-0.5*(*localColorDistance[0])/m_colorVariation[0]);
			}
			if((*localColorDistance[1])>m_thresholdColor[1]){
				kernel_color[1]=0;
			}
			else{
				kernel_color[1]=kernel_distance
					*exp(-0.5*(*localColorDistance[1])/m_colorVariation[1]);
			}
			if((*localColorDistance[2])>m_thresholdColor[2]){
				kernel_color[2]=0;
			}
			else{
				kernel_color[2]=kernel_distance
					*exp(-0.5*(*localColorDistance[2])/m_colorVariation[2]);
			}

			sumKernel_position+=kernel_distance;
			sumKernel_colors[0]+=kernel_color[0];
			sumKernel_colors[1]+=kernel_color[1];
			sumKernel_colors[2]+=kernel_color[2];
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]+=(*localPositionEstimationIterator).pointVector[j]*kernel_distance;
				m_resultColors[i*3+j]+=m_originalColors[(*vectorIterator)*3+j]*kernel_color[j];
			}
			localProjectionDistanceIterator++;
			localPositionEstimationIterator++;
			localColorEstimation++;
			vectorIterator++;
			for(int j=0;j<3;j++){
				localColorDistance[j]++;
				localColorProjectionDistance[j]++;			
			}
		}

		sumKernel_position+=1;
		sumKernel_colors[0]+=1;
		sumKernel_colors[1]+=1;
		sumKernel_colors[2]+=1;
		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]+=m_originalPointSet[i*3+j];
			m_resultPointSet[i*3+j]/=sumKernel_position;
			m_resultColors[i*3+j]+=m_originalColors[i*3+j];
			m_resultColors[i*3+j]/=sumKernel_colors[j];
			if(m_resultColors[i*3+j]<0||m_resultColors[i*3+j]>1){
				int w;
				w=0;
			}
			else{

			}
		}
		m_localPointDistance.clear();
		m_localProjectionDistance.clear();
		m_localPositionEstimation.clear();
		m_localColorEstimation.clear();
		for(int j=0;j<3;j++){
			m_localColorDistance[j].clear();
			m_localColorProjectionDistance[j].clear();	
		}	
	}
	delete[] m_originalPointSet;
	delete[] m_originalColors;
	m_originalPointSet=NULL;
	m_originalColors=NULL;


	//////////////////////////////////////////////////////////////////////////
	// 测试代码
	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\resultColorsFive.txt";

	FILE *fpout;
	float averageColor[3];
	float variationColor[3];
	for(int j=0;j<3;j++){
		averageColor[j]=0;
		variationColor[j]=0;
	}
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f\n",m_resultColors[i*3],m_resultColors[i*3+1],m_resultColors[i*3+2]);
			averageColor[0]+=m_resultColors[i*3];
			averageColor[1]+=m_resultColors[i*3+1];
			averageColor[2]+=m_resultColors[i*3+2];
		}
		for(int i=0;i<3;i++){
			averageColor[i]/=m_numOfPoints;
		}
		fprintf(fpout,"%f %f %f\n",averageColor[0],averageColor[1],averageColor[2]);
		for(int i=0;i<m_numOfPoints;i++){
			variationColor[0]+=(m_resultColors[i*3]-averageColor[0])*(m_resultColors[i*3]-averageColor[0]);
			variationColor[1]+=(m_resultColors[i*3+1]-averageColor[1])*(m_resultColors[i*3+1]-averageColor[1]);
			variationColor[2]+=(m_resultColors[i*3+2]-averageColor[2])*(m_resultColors[i*3+2]-averageColor[2]);
		}
		for(int i=0;i<3;i++){
			variationColor[i]/=(m_numOfPoints-1);
		}
		fprintf(fpout,"%f %f %f\n",variationColor[0],variationColor[1],variationColor[2]);
		fclose(fpout);
	}
	//////////////////////////////////////////////////////////////////////////
}
void FilterBilateral::CalculateBilateralFilterEight()
{
    m_pointVariation=0.2;
	m_resultPointSet=new float[m_numOfPoints*3];

	for(int i=0;i<m_numOfPoints;i++){
		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]=0;
	
		}
	}
	float sumKernel_position;

	for(int i=0;i<m_numOfPoints;i++){
		std::map<int ,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		sumKernel_position=0;
		ComputeOneOrderPositionEstimation(i);
		//CalculateVariation();
		std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();
		std::vector<float>::iterator localProjectionDistanceIterator=m_localProjectionDistance.begin();
		std::vector<POINTVECTOR3D>::iterator localPositionEstimationIterator=m_localPositionEstimation.begin();
	
		for(;localPointDistanceIterator!=m_localPointDistance.end();localPointDistanceIterator++){
			float kernel_distance;
		
			kernel_distance=exp(-0.5*(*localPointDistanceIterator)/m_pointVariation)*exp(-0.5*(*localProjectionDistanceIterator)/m_estimationVariation);
		
			

			sumKernel_position+=kernel_distance;
	
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]+=(*localPositionEstimationIterator).pointVector[j]*kernel_distance;
			}
			localProjectionDistanceIterator++;
			localPositionEstimationIterator++;
		
			vectorIterator++;
	
		}

		sumKernel_position+=1;

		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]+=m_originalPointSet[i*3+j];
			m_resultPointSet[i*3+j]/=sumKernel_position;
			
		}
		m_localPointDistance.clear();
		m_localProjectionDistance.clear();
		m_localPositionEstimation.clear();

	
	}
	delete[] m_originalPointSet;
	delete[] m_originalColors;
	m_originalPointSet=NULL;
	m_originalColors=NULL;


	//////////////////////////////////////////////////////////////////////////
	// 测试代码
	//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\resultColorsFive.txt";

	//FILE *fpout;
	//float averageColor[3];
	//float variationColor[3];
	//for(int j=0;j<3;j++){
	//	averageColor[j]=0;
	//	variationColor[j]=0;
	//}
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f\n",m_resultColors[i*3],m_resultColors[i*3+1],m_resultColors[i*3+2]);
	//		averageColor[0]+=m_resultColors[i*3];
	//		averageColor[1]+=m_resultColors[i*3+1];
	//		averageColor[2]+=m_resultColors[i*3+2];
	//	}
	//	for(int i=0;i<3;i++){
	//		averageColor[i]/=m_numOfPoints;
	//	}
	//	fprintf(fpout,"%f %f %f\n",averageColor[0],averageColor[1],averageColor[2]);
	//	for(int i=0;i<m_numOfPoints;i++){
	//		variationColor[0]+=(m_resultColors[i*3]-averageColor[0])*(m_resultColors[i*3]-averageColor[0]);
	//		variationColor[1]+=(m_resultColors[i*3+1]-averageColor[1])*(m_resultColors[i*3+1]-averageColor[1]);
	//		variationColor[2]+=(m_resultColors[i*3+2]-averageColor[2])*(m_resultColors[i*3+2]-averageColor[2]);
	//	}
	//	for(int i=0;i<3;i++){
	//		variationColor[i]/=(m_numOfPoints-1);
	//	}
	//	fprintf(fpout,"%f %f %f\n",variationColor[0],variationColor[1],variationColor[2]);
	//	fclose(fpout);
	//}
	//////////////////////////////////////////////////////////////////////////
}
void FilterBilateral::CalculateColorGradientFair()
{
	float vectorOne[3], vectorTwo[3], localNormal[3];
	float coffecientOne,coffecientTwo, cofecientD;// 梯度 
	float coffecientMatrix[3][3], coffecientConstantRed[3],coffecientConstantGreen[3],coffecientConstantBlue[3];
	float localPosition[3];
	float localEXOne,localEXTwo;
	m_originalColorGradients[0]=new float[m_numOfPoints*4];
	m_originalColorGradients[1]=new float[m_numOfPoints*4];
	m_originalColorGradients[2]=new float[m_numOfPoints*4];
	for(int i=0;i<m_numOfPoints;i++){
		//首先求与法向量垂直的单位向量
		localNormal[0]=m_filterNormals[i*3];
		localNormal[1]=m_filterNormals[i*3+1];
		localNormal[2]=m_filterNormals[i*3+2];
		Perpendicular(localNormal,vectorOne);
		Cross3(localNormal, vectorOne,vectorTwo);
		// 计算方程组的系数
		for(int j=0;j<3;j++){
			//coffecientConstant[j]=0;
			for(int k=0;k<3;k++){
				coffecientMatrix[j][k]=0;
			}
		}
		std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		localPosition[0]=m_originalPointSet[i*3];
		localPosition[1]=m_originalPointSet[i*3+1];
		localPosition[2]=m_originalPointSet[i*3+2];
		localEXOne=Vector3Vector(localPosition,vectorOne);
		localEXTwo=Vector3Vector(localPosition,vectorTwo);
		coffecientMatrix[0][0]+=localEXOne*localEXOne;
		coffecientMatrix[0][1]+=localEXOne*localEXTwo;
		coffecientMatrix[0][2]+=localEXOne;
		coffecientMatrix[1][1]+=localEXTwo*localEXTwo;
		coffecientMatrix[1][2]+=localEXTwo;
		coffecientMatrix[2][2]=1;
		coffecientConstantRed[0]=m_originalColors[i*3]*localEXOne;
		coffecientConstantRed[1]=m_originalColors[i*3]*localEXTwo;
		coffecientConstantRed[2]=m_originalColors[i*3];
		coffecientConstantGreen[0]=m_originalColors[i*3+1]*localEXOne;
		coffecientConstantGreen[1]=m_originalColors[i*3+1]*localEXTwo;
		coffecientConstantGreen[2]=m_originalColors[i*3+1];
		coffecientConstantBlue[0]=m_originalColors[i*3+2]*localEXOne;
		coffecientConstantBlue[1]=m_originalColors[i*3+2]*localEXTwo;
		coffecientConstantBlue[2]=m_originalColors[i*3+2];
		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
            float normal_to_normal_vector[3];
			float variation_normal=0.6;
			float fairing_coffecient;
			localPosition[0]=m_originalPointSet[(*vectorIterator)*3];
			localPosition[1]=m_originalPointSet[(*vectorIterator)*3+1];
			localPosition[2]=m_originalPointSet[(*vectorIterator)*3+2];
			normal_to_normal_vector[0]=m_filterNormals[i*3]-m_filterNormals[(*vectorIterator)*3];
			normal_to_normal_vector[1]=m_filterNormals[i*3+1]-m_filterNormals[(*vectorIterator)*3+1];
			normal_to_normal_vector[2]=m_filterNormals[i*3+2]-m_filterNormals[(*vectorIterator)*3+2];
			fairing_coffecient=exp(-0.5*(normal_to_normal_vector[0]*normal_to_normal_vector[0]
							+normal_to_normal_vector[1]*normal_to_normal_vector[1]
							+normal_to_normal_vector[2]*normal_to_normal_vector[2])/variation_normal);
			localEXOne=Vector3Vector(localPosition,vectorOne);
			localEXTwo=Vector3Vector(localPosition,vectorTwo);
			coffecientMatrix[0][0]+=localEXOne*localEXOne*fairing_coffecient;
			coffecientMatrix[0][1]+=localEXOne*localEXTwo*fairing_coffecient;
			coffecientMatrix[0][2]+=localEXOne*fairing_coffecient;
			coffecientMatrix[1][1]+=localEXTwo*localEXTwo*fairing_coffecient;
			coffecientMatrix[1][2]+=localEXTwo*fairing_coffecient;
			coffecientMatrix[2][2]+=fairing_coffecient;
			coffecientConstantRed[0]+=m_originalColors[(*vectorIterator)*3]*localEXOne*fairing_coffecient;
			coffecientConstantRed[1]+=m_originalColors[(*vectorIterator)*3]*localEXTwo*fairing_coffecient;
			coffecientConstantRed[2]+=m_originalColors[(*vectorIterator)*3]*fairing_coffecient;
			coffecientConstantGreen[0]+=m_originalColors[(*vectorIterator)*3+1]*localEXOne*fairing_coffecient;
			coffecientConstantGreen[1]+=m_originalColors[(*vectorIterator)*3+1]*localEXTwo*fairing_coffecient;
			coffecientConstantGreen[2]+=m_originalColors[(*vectorIterator)*3+1]*fairing_coffecient;
			coffecientConstantBlue[0]+=m_originalColors[(*vectorIterator)*3+2]*localEXOne*fairing_coffecient;
			coffecientConstantBlue[1]+=m_originalColors[(*vectorIterator)*3+2]*localEXTwo*fairing_coffecient;
			coffecientConstantBlue[2]+=m_originalColors[(*vectorIterator)*3+2]*fairing_coffecient;
		}

		//求方程的解 首先求矩阵的拟,然后求解
		mat_f8 m_matrix(3,3);
		vec_f8 m_vectorRed(3),m_vectorGreen(3),m_vectorBlue(3);
		vec_f8 m_gradientRed(3),m_gradientGreen(3),m_gradientBlue(3);
		coffecientMatrix[1][0]=coffecientMatrix[0][1];
		coffecientMatrix[2][0]=coffecientMatrix[0][2];
		coffecientMatrix[2][1]=coffecientMatrix[1][2];
//		coffecientMatrix[2][2]=((*mapIterator).second.m_numOfNearest+1)*fairing_coffecient;
		for(int j=0;j<3;j++){
			for(int k=0;k<3;k++){
				m_matrix(j,k)=coffecientMatrix[j][k];
			}
			m_vectorRed(j)=coffecientConstantRed[j];
			m_vectorGreen(j)=coffecientConstantGreen[j];
			m_vectorBlue(j)=coffecientConstantBlue[j];
		}
		if(matrix_inverse(m_matrix)){
			m_gradientRed=m_matrix*m_vectorRed;
			m_gradientGreen=m_matrix*m_vectorGreen;
			m_gradientBlue=m_matrix*m_vectorBlue;
		}
		for(int j=0;j<3;j++){
			m_originalColorGradients[0][i*4+j]=m_gradientRed(0)*vectorOne[j]+m_gradientRed(1)*vectorTwo[j];
			m_originalColorGradients[1][i*4+j]=m_gradientGreen(0)*vectorOne[j]+m_gradientGreen(1)*vectorTwo[j];
			m_originalColorGradients[2][i*4+j]=m_gradientBlue(0)*vectorOne[j]+m_gradientBlue(1)*vectorTwo[j];
		}
		m_originalColorGradients[0][i*4+3]=m_gradientRed(2);
		m_originalColorGradients[1][i*4+3]=m_gradientGreen(2);
		m_originalColorGradients[2][i*4+3]=m_gradientBlue(2);

	}

	//////////////////////////////////////////////////////////////////////////
	// 测试代码
	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\GradientRed.txt";
	FILE *fpout;
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f %f\n",m_originalColorGradients[0][i*4],m_originalColorGradients[0][i*4+1],m_originalColorGradients[0][i*4+2],m_originalColorGradients[0][i*4+3]);

		}
		fclose(fpout);
	}

	filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\GradientBlue.txt";

	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f %f\n",m_originalColorGradients[2][i*4],m_originalColorGradients[2][i*4+1],m_originalColorGradients[2][i*4+2],m_originalColorGradients[2][i*4+3]);

		}
		fclose(fpout);
	}

	filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\GradientGreen.txt";

	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f %f\n",m_originalColorGradients[1][i*4],m_originalColorGradients[1][i*4+1],m_originalColorGradients[1][i*4+2],m_originalColorGradients[1][i*4+3]);

		}
		fclose(fpout);
	}


	//////////////////////////////////////////////////////////////////////////


	return;
}
