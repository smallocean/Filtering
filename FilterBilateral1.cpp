// FilterBilateral.cpp: implementation of the FilterBilateral class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Mcube.h"
#include "FilterBilateral1.h"
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




FilterBilateral1::FilterBilateral1()

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
	m_gradientKernelWide=0;
	m_thresholdColorGradient=0;

}

FilterBilateral1::~FilterBilateral1()
{
	DeleteFilterBilateral();	
}

void FilterBilateral1::GetFilterBilateral(int numOfPoints, float* pointSet,float* colors,float meanLength,float variationF,float variationG1,float variationG2,float variationH1,float variationH2,float variationH3,int knearest,int meanShift,int interactiveTimes,int indexFunction,float variationNstop,float functionVariation,float thresholdOfColor,float functionGradientWide,float functionGradientThreshold,float thresholdDistanceNormal,float thresholdDistanceGradientFunction)
{
	long nearestStart,nearestFinish;
	long filterStart,filterFinish;
	m_numOfPoints=numOfPoints;
	m_originalPointSet=new float[m_numOfPoints*3];
	m_resultColors=new float[m_numOfPoints*3];
	//m_originalColors=new float[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints*3;i++){
		m_originalPointSet[i]=pointSet[i];
	//	m_originalColors[i]=colors[i];
	}
	//////////////////////////////////////////////////////////////////////////
	// 参数设置	
	m_pointVariation=(variationF*meanLength)*(variationF*meanLength);
	m_estimationVariation=(variationG1*meanLength)*(variationG1*meanLength);
	m_meanShiftHposition=variationH1*meanLength*variationH1*meanLength;
	m_meanShiftHnormal=variationH2;
	m_meanShiftStopNormal=variationNstop;
	m_knearest=knearest;
	m_thresholdDistanceNormal=thresholdDistanceNormal;
	m_thresholdDistanceGradientFunction=thresholdDistanceGradientFunction;
	m_radius=10;
	for (int i=0;i<3;i++){
		m_colorVariation[i]=variationG2;
		m_thresholdColor[i]=thresholdOfColor;
		m_colorEstimationVariation[i]=functionVariation;//
	}
	m_functionWeight=variationH3;
	m_gradientKernelWide=functionGradientWide;
	m_thresholdColorGradient=functionGradientThreshold;
	
	//////////////////////////////////////////////////////////////////////////	
	nearestStart=clock();
	ComputeMapKnearest();
	nearestFinish=clock();
 	nearestTime=nearestFinish-nearestStart;
	filterStart=clock();
	CalculateOriginalNormals();
	if(meanShift==1){
		MeanShiftNormals();
	}
	else{
		m_filterNormals=new float[m_numOfPoints*3];
		for(int i=0;i<m_numOfPoints*3;i++){
			m_filterNormals[i]=m_originalNormals[i];
		}
	}


	//MidianFilter();
	if(indexFunction==1){
	
			CalculateBilateralFilterOne(interactiveTimes);

	
	}
	if(indexFunction==2){
		
			CalculateBilateralFilterTwo(interactiveTimes);
	
		
	}
	if(indexFunction==3){
		CalculateBilateralFilterThree(interactiveTimes);		
	}
	if(indexFunction==4){
		
			CalculateBilateralFilterFour(interactiveTimes);
		
		
	}
	if(indexFunction==5){
		
			CalculateBilateralFilterFive(interactiveTimes);

		
	}
	if(indexFunction==6){
		
			CalculateBilateralFilterSix(interactiveTimes);		
		
	}
	if(indexFunction==7){
	
			CalculateBilateralFilterSeven(interactiveTimes);
		
		
	}
	if(indexFunction==8){	
			CalculateBilateralFilterEight(interactiveTimes);
			
	}
	if(indexFunction==10){	
		CalculateBilateralFilterTen(interactiveTimes);

	}
	if(indexFunction==11){	
		CalculateBilateralFilterEleven(interactiveTimes);

	}
	if(indexFunction==12){	
		CalculateBilateralFilterTwelve(interactiveTimes);

	}
	
	if(indexFunction==0){
		if(m_resultPointSet!=NULL)
			delete[] m_resultPointSet;
        m_resultPointSet=new float[m_numOfPoints*3];
		for(int i=0;i<m_numOfPoints*3;i++){
			m_resultPointSet[i]=m_originalPointSet[i];
		}

	}

	filterFinish=clock();
	filterTime=filterFinish-filterStart;
	if(m_originalPointSet!=NULL)
		delete[] m_originalPointSet;
	m_originalPointSet=new float[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints*3;i++){
		m_originalPointSet[i]=m_resultPointSet[i];
	}
	
	

	CalculateOriginalNormals();
	
	delete[] m_originalPointSet;
	m_originalPointSet=NULL;
//	m_mapKnearest.clear();


}

void FilterBilateral1::DeleteFilterBilateral()
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
void FilterBilateral1::ComputeMapKnearest()
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
	//CString filename_pw = "D:\\qin\\sourcecode\\Shiftknearest.txt";
	//FILE *fpout;
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	m_meanRadius=0;
	//	for(int i=0;i<m_numOfPoints;i++){
	//		std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(i);
	//		std::vector<int>::iterator vectorIterator=(*mapKnearestIterator).second.m_nearest.begin();
	//		fprintf(fpout,"%d ",i);
	//		float tempDistance=0;
	//		for(;vectorIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorIterator++){
	//			fprintf(fpout,"%d ",(*vectorIterator));
	//			tempDistance+=sqrt((m_originalPointSet[i*3]-m_originalPointSet[(*vectorIterator)*3])*(m_originalPointSet[i*3]-m_originalPointSet[(*vectorIterator)*3])
	//				+(m_originalPointSet[i*3+1]-m_originalPointSet[(*vectorIterator)*3+1])*(m_originalPointSet[i*3+1]-m_originalPointSet[(*vectorIterator)*3+1])
	//				+(m_originalPointSet[i*3+2]-m_originalPointSet[(*vectorIterator)*3+2])*(m_originalPointSet[i*3+2]-m_originalPointSet[(*vectorIterator)*3+2]));

	//		}
	//		m_meanRadius+=tempDistance/m_knearest;
	//		fprintf(fpout,"\n");
	//	}
	//	m_meanRadius/=m_numOfPoints;
	//	fclose(fpout);
	//}
	//
	
	//////////////////////////////////////////////////////////////////////////
}

void FilterBilateral1::CalculateOriginalNormals()
{
	
	float centroidPosition[3];//重心位置
	float* localVariationVector;//重心到点的位置的向量
	float* m_meanAreas;
	float* m_meanCurvature;
	m_meanAreas=new float[m_numOfPoints];
	m_meanCurvature=new float[m_numOfPoints];
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
		m_meanAreas[i]=abs(eval(0)*eval(2));
		m_meanCurvature[i]=eval(2)/(eval(0)+eval(1)+eval(2));
		if(m_meanCurvature[i]*4<=0.01){
			m_resultColors[i*3]=0;
			m_resultColors[i*3+1]=0;
			m_resultColors[i*3+2]=0.2+70*m_meanCurvature[i]*4;
		}
		if((m_meanCurvature[i]*4)>0.01&(m_meanCurvature[i]*4)<0.05){
			m_resultColors[i*3]=0.2+20*(m_meanCurvature[i]*4-0.01);
			m_resultColors[i*3+1]=0;
			m_resultColors[i*3+2]=0.2+20*(m_meanCurvature[i]*4-0.01);
		}
		if((m_meanCurvature[i]*4>0.05)&(m_meanCurvature[i]*4<=0.1)){
			m_resultColors[i*3]=0;
			m_resultColors[i*3+1]=0.45+5*(m_meanCurvature[i]*4-0.05);
			m_resultColors[i*3+2]=0;
		}
		if((m_meanCurvature[i]*4>0.1)&(m_meanCurvature[i]*4<=0.2)){
			m_resultColors[i*3]=0.45+5*(m_meanCurvature[i]*4-0.1);
			m_resultColors[i*3+1]=0.45+5*(m_meanCurvature[i]*4-0.1);
			m_resultColors[i*3+2]=0;
		}
		if(m_meanCurvature[i]*4>0.2){
			m_resultColors[i*3]=0.4+5*(m_meanCurvature[i]*4-0.1);
			m_resultColors[i*3+1]=0;
			m_resultColors[i*3+2]=0;
		}

		delete[] localVariationVector;
	}
	delete[] m_meanAreas;
	delete[] m_meanCurvature;
	EstimateNormalDirection();
	

	CString filename_pw = "D:\\qin\\sourcecode\\ShiftoriginalNormal.txt";
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

void FilterBilateral1::MeanShiftNormals()
{
	float distance_position;
	float distance_normal;
	float position_kernel,range_kernel;
	float sum_kernel;
	float length_normal_vector;
	int num_stop=0;//已经停止的数
	float tempNormal[3];
	bool* flag_stop=new bool[m_numOfPoints];
	for(int i=0;i<m_numOfPoints;i++){
		flag_stop[i]=false;
	}

    m_filterNormals=new float[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints*3;i++){
		m_filterNormals[i]=m_originalNormals[i];
	}
	while(num_stop<m_numOfPoints){
		for(int i=0;i<m_numOfPoints;i++){
			if(!flag_stop[i]){
				sum_kernel=0;
				for(int j=0;j<3;j++){
                    tempNormal[j]=0;
				}			
				std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(i);
				//int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
				std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
				for(;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++){
					distance_position=(m_originalPointSet[i*3]-m_originalPointSet[(*vectorNearestIterator)*3])*(m_originalPointSet[i*3]-m_originalPointSet[(*vectorNearestIterator)*3])
						+(m_originalPointSet[i*3+1]-m_originalPointSet[(*vectorNearestIterator)*3+1])*(m_originalPointSet[i*3+1]-m_originalPointSet[(*vectorNearestIterator)*3+1])
						+(m_originalPointSet[i*3+2]-m_originalPointSet[(*vectorNearestIterator)*3+2])*(m_originalPointSet[i*3+2]-m_originalPointSet[(*vectorNearestIterator)*3+2]);
					distance_normal=(m_filterNormals[i*3]-m_originalNormals[(*vectorNearestIterator)*3])*(m_filterNormals[i*3]-m_originalNormals[(*vectorNearestIterator)*3])
						+(m_filterNormals[i*3+1]-m_originalNormals[(*vectorNearestIterator)*3+1])*(m_filterNormals[i*3+1]-m_originalNormals[(*vectorNearestIterator)*3+1])
						+(m_filterNormals[i*3+2]-m_originalNormals[(*vectorNearestIterator)*3+2])*(m_filterNormals[i*3+2]-m_originalNormals[(*vectorNearestIterator)*3+2]);
					position_kernel=exp(-0.5*distance_position/m_meanShiftHposition);
					range_kernel=exp(-0.5*distance_normal/m_meanShiftHnormal);
					sum_kernel+=position_kernel*range_kernel;
					
					for(int j=0;j<3;j++){
						tempNormal[j]+=m_originalNormals[(*vectorNearestIterator)*3+j]*position_kernel*range_kernel;
					}
				}
				distance_normal=(m_filterNormals[i*3]-m_originalNormals[i*3])*(m_filterNormals[i*3]-m_originalNormals[i*3])
					+(m_filterNormals[i*3+1]-m_originalNormals[i*3+1])*(m_filterNormals[i*3+1]-m_originalNormals[i*3+1])
					+(m_filterNormals[i*3+2]-m_originalNormals[i*3+2])*(m_filterNormals[i*3+2]-m_originalNormals[i*3+2]);
				range_kernel=exp(-0.5*distance_normal/m_meanShiftHnormal);
				sum_kernel+=range_kernel;
				
				for(int j=0;j<3;j++){
					tempNormal[j]+=m_originalNormals[i*3+j]*range_kernel;
					tempNormal[j]/=sum_kernel;
				}
				length_normal_vector=sqrt(Vector3Vector(tempNormal,tempNormal));
				for(int j=0;j<3;j++){
					tempNormal[j]/=length_normal_vector;
					//m_filterNormals[i*3+j]=tempNormal[j];
				}								
				if(((m_filterNormals[i*3]-tempNormal[0])*(m_filterNormals[i*3]-tempNormal[0])
					+(m_filterNormals[i*3+1]-tempNormal[1])*(m_filterNormals[i*3+1]-tempNormal[1])
					+(m_filterNormals[i*3+2]-tempNormal[2])*(m_filterNormals[i*3+2]-tempNormal[2]))<m_meanShiftStopNormal){
						flag_stop[i]=true;
						num_stop+=1;
				}
					for(int j=0;j<3;j++){
						m_filterNormals[i*3+j]=tempNormal[j];
					//	m_originalNormals[i*3+j]=tempNormal[j];
					}
			}
		}

	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;
	////////////////////////////////////////////////////////////////////////////
	//// 测试代码
	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\ShiftfilterNormal.txt";
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

void FilterBilateral1::CalculateColorGradient()
{
	float vectorOne[3], vectorTwo[3], localNormal[3];
	float coffecientOne,coffecientTwo, cofecientD;// 梯度 
	float coffecientMatrix[3][3], coffecientConstantRed[3],coffecientConstantGreen[3],coffecientConstantBlue[3];
	float localPosition[3];
	float localWeight;//求颜色梯度时候的加权值
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
		coffecientMatrix[2][2]+=1;
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
			localWeight=exp(sqrt((m_filterNormals[i*3]-m_filterNormals[(*vectorIterator)*3])*(m_filterNormals[i*3]-m_filterNormals[(*vectorIterator)*3])
				+(m_filterNormals[i*3+1]-m_filterNormals[(*vectorIterator)*3+1])*(m_filterNormals[i*3+1]-m_filterNormals[(*vectorIterator)*3+1])
			    +(m_filterNormals[i*3+2]-m_filterNormals[(*vectorIterator)*3+2])*(m_filterNormals[i*3+2]-m_filterNormals[(*vectorIterator)*3+2]))/m_functionWeight);
			coffecientMatrix[0][0]+=localEXOne*localEXOne*localWeight;
			coffecientMatrix[0][1]+=localEXOne*localEXTwo*localWeight;
			coffecientMatrix[0][2]+=localEXOne*localWeight;
			coffecientMatrix[1][1]+=localEXTwo*localEXTwo*localWeight;
			coffecientMatrix[1][2]+=localEXTwo*localWeight;
			coffecientMatrix[2][2]+=localWeight;
			coffecientConstantRed[0]+=m_originalColors[(*vectorIterator)*3]*localEXOne*localWeight;
			coffecientConstantRed[1]+=m_originalColors[(*vectorIterator)*3]*localEXTwo*localWeight;
			coffecientConstantRed[2]+=m_originalColors[(*vectorIterator)*3]*localWeight;
			coffecientConstantGreen[0]+=m_originalColors[(*vectorIterator)*3+1]*localEXOne*localWeight;
			coffecientConstantGreen[1]+=m_originalColors[(*vectorIterator)*3+1]*localEXTwo*localWeight;
			coffecientConstantGreen[2]+=m_originalColors[(*vectorIterator)*3+1]*localWeight;
			coffecientConstantBlue[0]+=m_originalColors[(*vectorIterator)*3+2]*localEXOne*localWeight;
			coffecientConstantBlue[1]+=m_originalColors[(*vectorIterator)*3+2]*localEXTwo*localWeight;
			coffecientConstantBlue[2]+=m_originalColors[(*vectorIterator)*3+2]*localWeight;
		}

		//求方程的解 首先求矩阵的拟,然后求解
		mat_f8 m_matrix(3,3);
		vec_f8 m_vectorRed(3),m_vectorGreen(3),m_vectorBlue(3);
		vec_f8 m_gradientRed(3),m_gradientGreen(3),m_gradientBlue(3);
		coffecientMatrix[1][0]=coffecientMatrix[0][1];
		coffecientMatrix[2][0]=coffecientMatrix[0][2];
		coffecientMatrix[2][1]=coffecientMatrix[1][2];
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
	//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\GradientRed.txt";
	//FILE *fpout;
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f %f\n",m_originalColorGradients[0][i*4],m_originalColorGradients[0][i*4+1],m_originalColorGradients[0][i*4+2],m_originalColorGradients[0][i*4+3]);

	//	}
	//	fclose(fpout);
	//}

	// filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\GradientBlue.txt";

	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f %f\n",m_originalColorGradients[2][i*4],m_originalColorGradients[2][i*4+1],m_originalColorGradients[2][i*4+2],m_originalColorGradients[2][i*4+3]);

	//	}
	//	fclose(fpout);
	//}

	// filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\GradientGreen.txt";
	//
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f %f\n",m_originalColorGradients[1][i*4],m_originalColorGradients[1][i*4+1],m_originalColorGradients[1][i*4+2],m_originalColorGradients[1][i*4+3]);

	//	}
	//	fclose(fpout);
	//}


	//////////////////////////////////////////////////////////////////////////
	

	return;
}

void FilterBilateral1::CalculateBilateralFilterOne(int interactiveTimes)
{
	//有shift的线性假定
	m_resultPointSet=new float[m_numOfPoints*3];
	m_resultColors=new float[m_numOfPoints*3];
	for(int www=0;www<interactiveTimes;www++){
		CalculateColorGradient();
		MeanShiftColorGradient();
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
				
				if((m_filterNormals[i*3]*m_filterNormals[(*vectorIterator)*3]
					+m_filterNormals[i*3+1]*m_filterNormals[(*vectorIterator)*3+1]
					+m_filterNormals[i*3+2]*m_filterNormals[(*vectorIterator)*3+2])<m_thresholdDistanceNormal)
				continue;
				float kernel_distance;
				float kernel_color[3];
				/*for(int j=0;j<3;j++){
				m_colorVariation[j]=0.1;
				m_colorEstimationVariation[j]=.1;
				}*/
				kernel_distance=exp(-0.5*(*localPointDistanceIterator)/m_pointVariation)*exp(-0.5*(*localProjectionDistanceIterator)/m_estimationVariation);
				if((m_originalColorGradients[0][i*4]-m_originalColorGradients[0][(*vectorIterator)*4])*(m_originalColorGradients[0][i*4]-m_originalColorGradients[0][(*vectorIterator)*4])
					+(m_originalColorGradients[0][i*4+1]-m_originalColorGradients[0][(*vectorIterator)*4+1])*(m_originalColorGradients[0][i*4+1]-m_originalColorGradients[0][(*vectorIterator)*4+1])
					+(m_originalColorGradients[0][i*4+2]-m_originalColorGradients[0][(*vectorIterator)*4+2])*(m_originalColorGradients[0][i*4+2]-m_originalColorGradients[0][(*vectorIterator)*4+2])
					+(m_originalColorGradients[0][i*4+2]-m_originalColorGradients[0][(*vectorIterator)*4+2])*(m_originalColorGradients[0][i*4+2]-m_originalColorGradients[0][(*vectorIterator)*4+2])>m_thresholdDistanceGradientFunction){
						kernel_color[0]=0;
					}
				else{
					if(((*localColorDistance[0])>m_thresholdColor[0])||((*localColorProjectionDistance[0])>m_thresholdColor[0])){
						kernel_color[0]=0;				
					}
					else{
						kernel_color[0]=kernel_distance
							//*exp(-0.5*(*localColorDistance[0])/m_colorVariation[0])
							*exp(-0.5*(*localColorProjectionDistance[0])/m_colorEstimationVariation[0]);
					}
				}

				if((m_originalColorGradients[1][i*4]-m_originalColorGradients[1][(*vectorIterator)*4])*(m_originalColorGradients[1][i*4]-m_originalColorGradients[1][(*vectorIterator)*4])
					+(m_originalColorGradients[1][i*4+1]-m_originalColorGradients[1][(*vectorIterator)*4+1])*(m_originalColorGradients[1][i*4+1]-m_originalColorGradients[1][(*vectorIterator)*4+1])
					+(m_originalColorGradients[1][i*4+2]-m_originalColorGradients[1][(*vectorIterator)*4+2])*(m_originalColorGradients[1][i*4+2]-m_originalColorGradients[1][(*vectorIterator)*4+2])
					+(m_originalColorGradients[1][i*4+2]-m_originalColorGradients[1][(*vectorIterator)*4+2])*(m_originalColorGradients[1][i*4+2]-m_originalColorGradients[1][(*vectorIterator)*4+2])>m_thresholdDistanceGradientFunction){
						kernel_color[1]=0;
					}		
				else{
					if(((*localColorDistance[1])>m_thresholdColor[1])||((*localColorProjectionDistance[1])>m_thresholdColor[1])){
						kernel_color[1]=0;				
					}
					else{
						kernel_color[1]=kernel_distance
							//*exp(-0.5*(*localColorDistance[1])/m_colorVariation[1])
							*exp(-0.5*(*localColorProjectionDistance[1])/m_colorEstimationVariation[1]);
					}
				}

				if((m_originalColorGradients[2][i*4]-m_originalColorGradients[2][(*vectorIterator)*4])*(m_originalColorGradients[2][i*4]-m_originalColorGradients[2][(*vectorIterator)*4])
					+(m_originalColorGradients[2][i*4+1]-m_originalColorGradients[2][(*vectorIterator)*4+1])*(m_originalColorGradients[2][i*4+1]-m_originalColorGradients[2][(*vectorIterator)*4+1])
					+(m_originalColorGradients[2][i*4+2]-m_originalColorGradients[2][(*vectorIterator)*4+2])*(m_originalColorGradients[2][i*4+2]-m_originalColorGradients[2][(*vectorIterator)*4+2])
					+(m_originalColorGradients[2][i*4+2]-m_originalColorGradients[2][(*vectorIterator)*4+2])*(m_originalColorGradients[2][i*4+2]-m_originalColorGradients[2][(*vectorIterator)*4+2])>m_thresholdDistanceGradientFunction){
						kernel_color[2]=0;
					}	
				else{
					if(((*localColorDistance[2])>m_thresholdColor[2])||((*localColorProjectionDistance[2])>m_thresholdColor[2])){
						kernel_color[2]=0;				
					}
					else{
						kernel_color[2]=kernel_distance
							//*exp(-0.5*(*localColorDistance[2])/m_colorVariation[2])
							*exp(-0.5*(*localColorProjectionDistance[2])/m_colorEstimationVariation[2]);
					}

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
		for(int k=0;k<m_numOfPoints*3;k++){
			m_originalColors[k]=m_resultColors[k];
		}
	}
	
	delete[] m_originalPointSet;
	m_originalPointSet=NULL;
	delete[] m_originalColors;
	m_originalColors=NULL;

    //////////////////////////////////////////////////////////////////////////
    // 测试代码
	//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\resultColors.txt";
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

void FilterBilateral1::CalculateBilateralFilterTwo(int interactiveTimes)
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
			if((m_filterNormals[i*3]*m_filterNormals[(*vectorIterator)*3]
				+m_filterNormals[i*3+1]*m_filterNormals[(*vectorIterator)*3+1]
				+m_filterNormals[i*3+2]*m_filterNormals[(*vectorIterator)*3+2])<m_thresholdDistanceNormal)
				continue;
			float kernel_distance;
			float kernel_color[3];
			kernel_distance=exp(-0.5*(*localPointDistanceIterator)/m_pointVariation)*exp(-0.5*(*localProjectionDistanceIterator)/m_estimationVariation);
			if(((*localColorProjectionDistance[0])>m_thresholdColor[0])||((*localColorProjectionDistance[0])>m_thresholdColor[0])){
				kernel_color[0]=0;
			}
			else{
				kernel_color[0]=kernel_distance
					*exp(-0.5*(*localColorDistance[0])/m_colorVariation[0]);
				//	*exp(-0.5*(*localColorProjectionDistance[0])/m_colorEstimationVariation[0]);
			}

			if(((*localColorProjectionDistance[1])>m_thresholdColor[1])||((*localColorProjectionDistance[1])>m_thresholdColor[1])){
				kernel_color[1]=0;
			}
			else{
				kernel_color[1]=kernel_distance
					*exp(-0.5*(*localColorDistance[1])/m_colorVariation[1]);
					//*exp(-0.5*(*localColorProjectionDistance[1])/m_colorEstimationVariation[1]);
			}
			if(((*localColorProjectionDistance[2])>m_thresholdColor[2])||((*localColorProjectionDistance[2])>m_thresholdColor[2])){
				kernel_color[2]=0;
			}
			else{
				kernel_color[2]=kernel_distance
					*exp(-0.5*(*localColorDistance[2])/m_colorVariation[2]);
					//*exp(-0.5*(*localColorProjectionDistance[2])/m_colorEstimationVariation[2]);
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
	//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\resultColorsTwo.txt";

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
//////////////////////////////////////////////////////////////////////////
// 各项同性滤波 （颜色的各项同性）仅仅考虑 原始点的距离
void FilterBilateral1::CalculateBilateralFilterThree(int interactiveTimes)
{
	m_resultPointSet=new float[m_numOfPoints*3];
	m_resultColors=new float[m_numOfPoints*3];
	for(int www=0;www<interactiveTimes;www++){
		CalculateColorGradient();
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

	}
	
	delete[] m_originalPointSet;
	delete[] m_originalColors;
	m_originalPointSet=NULL;
	m_originalColors=NULL;


	//////////////////////////////////////////////////////////////////////////
	// 测试代码
	//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\resultColorsThree.txt";

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

void FilterBilateral1::ComputeOneOrderPositionEstimation(int indexPoint)
{

	std::map<int,KnearestField>::iterator mapIterator=m_mapKnearest.find(indexPoint);
	std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	float localNormal[3];//局部法向
	/*localNormal[0]=m_filterNormals[indexPoint*3];
	localNormal[1]=m_filterNormals[indexPoint*3+1];
	localNormal[2]=m_filterNormals[indexPoint*3+2];*/
	for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
	   float distance_point_point;//点到点的距离
	   float distance_point_projection;//点到一次估计平面的距离
	   float localVector_point_point[3];//点到点的向量
	   float localColorDistance[3];//颜色距离
	   POINTVECTOR3D localProjection;//投影点的位置
	   localNormal[0]=m_filterNormals[(*vectorIterator)*3];
	   localNormal[1]=m_filterNormals[(*vectorIterator)*3+1];
	   localNormal[2]=m_filterNormals[(*vectorIterator)*3+2];


	   localVector_point_point[0]=m_originalPointSet[(*vectorIterator)*3]-m_originalPointSet[indexPoint*3];
	   localVector_point_point[1]=m_originalPointSet[(*vectorIterator)*3+1]-m_originalPointSet[indexPoint*3+1];
	   localVector_point_point[2]=m_originalPointSet[(*vectorIterator)*3+2]-m_originalPointSet[indexPoint*3+2];

	   distance_point_projection=Vector3Vector(localVector_point_point,localNormal);

	   distance_point_point=(Vector3Vector(localVector_point_point,localVector_point_point));

	   localProjection.pointVector[0]=m_originalPointSet[indexPoint*3]+distance_point_projection*localNormal[0];
	   localProjection.pointVector[1]=m_originalPointSet[indexPoint*3+1]+distance_point_projection*localNormal[1];
	   localProjection.pointVector[2]=m_originalPointSet[indexPoint*3+2]+distance_point_projection*localNormal[2];
	/*   localColorDistance[0]=((m_originalColors[indexPoint*3]-m_originalColors[(*vectorIterator)*3])*(m_originalColors[indexPoint*3]-m_originalColors[(*vectorIterator)*3]));
	   localColorDistance[1]=((m_originalColors[indexPoint*3+1]-m_originalColors[(*vectorIterator)*3+1])*(m_originalColors[indexPoint*3+1]-m_originalColors[(*vectorIterator)*3+1]));
	   localColorDistance[2]=((m_originalColors[indexPoint*3+2]-m_originalColors[(*vectorIterator)*3+2])*(m_originalColors[indexPoint*3+2]-m_originalColors[(*vectorIterator)*3+2]));*/


	 
	   m_localProjectionDistance.push_back((distance_point_projection)*distance_point_projection);
	   m_localPointDistance.push_back(distance_point_point);
	   m_localPositionEstimation.push_back(localProjection);
	  /* m_localColorDistance[0].push_back(localColorDistance[0]);
	   m_localColorDistance[1].push_back(localColorDistance[1]);
	   m_localColorDistance[2].push_back(localColorDistance[2]);*/


	}
	//计算一次颜色估计和颜色投影距离
	//ComputeOneOrderColorEstimation(indexPoint);

    return;
}

void FilterBilateral1::ComputeOneOrderColorEstimation(int indexPoint)
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

void FilterBilateral1::CalculateVariation()
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
	if(m_pointVariation==0){
		m_pointVariation=1;
	}
	m_estimationVariation/=m_localPointDistance.size();
	if(m_estimationVariation==0){
		m_estimationVariation==1;
	}
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

void FilterBilateral1::AddPositionGaussianNoise()
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

void FilterBilateral1::AddColorGaussianNorse()
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
void FilterBilateral1::CalculateBilateralFilterFour(int interactiveTimes)
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
	//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\resultColorsFour.txt";

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
//////////////////////////////////////////////////////////////////////////
// 加入颜色信息 常数假定 但没有考虑梯度信息
void FilterBilateral1::CalculateBilateralFilterFive(int interactiveTimes)
{
	m_resultPointSet=new float[m_numOfPoints*3];
	m_resultColors=new float[m_numOfPoints*3];
	for(int www=0;www<interactiveTimes;www++){
		CalculateColorGradient();
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
		for(int k=0;k<m_numOfPoints*3;k++){
			m_originalColors[k]=m_resultColors[k];
		}
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
void FilterBilateral1::CalculateBilateralFilterSix(int interactiveTimes)
{
	m_resultPointSet=new float[m_numOfPoints*3];

	for(int www=0;www<interactiveTimes;www++)
	{
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

			std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();
			std::vector<float>::iterator localProjectionDistanceIterator=m_localProjectionDistance.begin();
			std::vector<POINTVECTOR3D>::iterator localPositionEstimationIterator=m_localPositionEstimation.begin();
			for(;localPointDistanceIterator!=m_localPointDistance.end();localPointDistanceIterator++){
				float kernel_distance;
				kernel_distance=exp(-0.5*(*localPointDistanceIterator)/m_pointVariation);
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
		for(int i=0;i<m_numOfPoints*3;i++){
			m_originalPointSet[i]=m_resultPointSet[i];
		}
	}
	delete[] m_originalPointSet;
	m_originalPointSet=NULL;	
}
void FilterBilateral1::CalculateBilateralFilterSeven(int interactiveTimes)
{

	m_resultPointSet=new float[m_numOfPoints*3];
	
	for(int www=0;www<interactiveTimes;www++)
	{
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

			std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();
			std::vector<float>::iterator localProjectionDistanceIterator=m_localProjectionDistance.begin();
			std::vector<POINTVECTOR3D>::iterator localPositionEstimationIterator=m_localPositionEstimation.begin();
			for(;localPointDistanceIterator!=m_localPointDistance.end();localPointDistanceIterator++){
				if((m_filterNormals[i*3]*m_filterNormals[(*vectorIterator)*3]
					+m_filterNormals[i*3+1]*m_filterNormals[(*vectorIterator)*3+1]
					+m_filterNormals[i*3+2]*m_filterNormals[(*vectorIterator)*3+2])<m_thresholdDistanceNormal)
				continue;
				float kernel_distance;
				kernel_distance=exp(-0.5*(*localPointDistanceIterator)/m_pointVariation)*exp(-0.5*(*localProjectionDistanceIterator)/m_estimationVariation);
				sumKernel_position+=kernel_distance;
				for(int j=0;j<3;j++){
					m_resultPointSet[i*3+j]+=(*localPositionEstimationIterator).pointVector[j]*kernel_distance;
				}
				localProjectionDistanceIterator++;
				localPositionEstimationIterator++;
				vectorIterator++;

				// test
				if((m_resultPointSet[i*3]>-9999999)&(m_resultPointSet[i*3]<9999999)){
					int n=0;
					n=2;
				}
				else{
					int n=0;
					n=2;
				}
					


			}
			sumKernel_position+=1;
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]+=m_originalPointSet[i*3+j];
				m_resultPointSet[i*3+j]/=sumKernel_position;
			}

			//////////////////////////////////////////////////////////////////////////
			// 调整滤波后的点在法向位置移动
			float tempdistance=0;
			tempdistance=((m_resultPointSet[i*3]-m_originalPointSet[i*3])*m_filterNormals[i*3]
				+(m_resultPointSet[i*3+1]-m_originalPointSet[i*3+1])*m_filterNormals[i*3+1]
				+(m_resultPointSet[i*3+2]-m_originalPointSet[i*3+2])*m_filterNormals[i*3+2]);
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]=m_originalPointSet[i*3+j]+tempdistance*m_filterNormals[i*3+j];
			}



			m_localPointDistance.clear();
			m_localProjectionDistance.clear();
			m_localPositionEstimation.clear();
		}
		for(int i=0;i<m_numOfPoints*3;i++){
			m_originalPointSet[i]=m_resultPointSet[i];
		}
	}	
		
	delete[] m_originalPointSet;

	m_originalPointSet=NULL;



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
void FilterBilateral1::CalculateBilateralFilterEight(int interactiveTimes)// diffusion
{
  
	m_resultPointSet=new float[m_numOfPoints*3];
	float nanmuda=0.1;
	for(int www=0;www<interactiveTimes;www++){
		for(int i=0;i<m_numOfPoints;i++){
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]=0;
			}
		}	

		for(int i=0;i<m_numOfPoints;i++){
			float weightcoffecient=0;
			std::map<int ,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
			std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
			ComputeOneOrderPositionEstimation(i);
			std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();

			for(;localPointDistanceIterator!=m_localPointDistance.end();localPointDistanceIterator++){

				
				float kernel_distance;
				if((*localPointDistanceIterator)<0.000001){
					kernel_distance=0;
				}
				else{
					kernel_distance=1/sqrt(*localPointDistanceIterator);
				}
				
				weightcoffecient+=kernel_distance;	

				for(int j=0;j<3;j++){
					m_resultPointSet[i*3+j]+=kernel_distance*(m_originalPointSet[(*vectorIterator)*3+j]-m_originalPointSet[i*3+j]);
				}
				vectorIterator++;
			}


			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]=m_originalPointSet[i*3+j]+m_resultPointSet[i*3+j]*nanmuda/weightcoffecient;						
			}
			m_localPointDistance.clear();
			m_localProjectionDistance.clear();
			m_localPositionEstimation.clear();	
		}
		//保护体积
		/*float* tempPointSet;
		tempPointSet=new float[m_numOfPoints*3];
		for(int i=0;i<m_numOfPoints;i++){
		std::map<int ,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		float tempVector[3];
		tempVector[0]=(m_originalPointSet[i*3]-m_resultPointSet[i*3])/(*mapIterator).second.m_nearest.size();
		tempVector[1]=(m_originalPointSet[i*3+1]-m_resultPointSet[i*3+1])/(*mapIterator).second.m_nearest.size();
		tempVector[2]=(m_originalPointSet[i*3+2]-m_resultPointSet[i*3+2])/(*mapIterator).second.m_nearest.size();
		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
		tempPointSet[(*vectorIterator)*3]+=tempVector[0];
		tempPointSet[(*vectorIterator)*3+1]+=tempVector[1];
		tempPointSet[(*vectorIterator)*3+2]+=tempVector[2];
		}			
		}*/
		for(int i=0;i<m_numOfPoints*3;i++){
			m_originalPointSet[i]=m_resultPointSet[i];//+tempPointSet[i];
		}
		//	delete[] tempPointSet;
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
void FilterBilateral1::CalculateColorGradientFair()
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

void FilterBilateral1::MidianFilter()
{ 
	float* tempColors=new float[m_numOfPoints*3];
	for (int i=0;i<m_numOfPoints;i++){
		std::map<int ,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		std::multimap<float,int> tempMidianVector[3];
		for(int k=0;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,k++){
			tempMidianVector[0].insert(std::multimap<float,int>::value_type(m_originalColors[(*vectorIterator)*3],(*vectorIterator)));
			tempMidianVector[1].insert(std::multimap<float,int>::value_type(m_originalColors[(*vectorIterator)*3+1],(*vectorIterator)));
			tempMidianVector[2].insert(std::multimap<float,int>::value_type(m_originalColors[(*vectorIterator)*3+2],(*vectorIterator)));
		}
		tempMidianVector[0].insert(std::multimap<float,int>::value_type(m_originalColors[i*3],i));
		tempMidianVector[1].insert(std::multimap<float,int>::value_type(m_originalColors[i*3+1],i));
		tempMidianVector[2].insert(std::multimap<float,int>::value_type(m_originalColors[i*3+2],i));
		for(int k=0;k<3;k++){
            std::multimap<float,int>::iterator multimapIterator=tempMidianVector[k].begin();
			for(int w=0;w<(int)(*mapIterator).second.m_nearest.size()/2;w++)
				multimapIterator++;
			tempColors[i*3+k]=(*multimapIterator).first;
		}
	}
	for(int i=0;i<m_numOfPoints;i++){
		for(int j=0;j<3;j++){
			m_originalColors[i*3+j]=tempColors[i*3+j];
		}
	}
	delete[] tempColors;
}
void FilterBilateral1::MeanShiftColorGradient()
{
	float distance_position;
	float distance_color_gradient[3];
	float* m_resultGradient[3];
	float distance_normal;
	float position_kernel,range_kernel;
	float sum_kernel[3];
	float tempGradient[3][4];
	float length_normal_vector;
	int num_stop[3];//已经停止的数
	bool* flag_stop[3];
	for(int i=0;i<3;i++){
		num_stop[i]=0;
		flag_stop[i]=new bool[m_numOfPoints];
		m_resultGradient[i]=new float[m_numOfPoints*4];
	}
	for(int i=0;i<m_numOfPoints;i++){
		for(int j=0;j<3;j++){
			flag_stop[j][i]=false;
			for(int k=0;k<4;k++){
				m_resultGradient[j][i*4+k]=m_originalColorGradients[j][i*4+k];
			}
		}	
		
	}

	
	while(num_stop[0]<m_numOfPoints||num_stop[1]<m_numOfPoints||num_stop[2]<m_numOfPoints){
		for(int i=0;i<m_numOfPoints;i++){
			if((!flag_stop[0][i])||(!flag_stop[1][i])||(!flag_stop[2][i])){
				sum_kernel[0]=0;
				sum_kernel[1]=0;
				sum_kernel[2]=0;

				for(int j=0;j<4;j++){				
                        tempGradient[0][j]=0;				
						tempGradient[1][j]=0;					
						tempGradient[2][j]=0;
				}

				std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(i);
				//int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
				std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
				for(;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++){
					distance_position=(m_originalPointSet[i*3]-m_originalPointSet[(*vectorNearestIterator)*3])*(m_originalPointSet[i*3]-m_originalPointSet[(*vectorNearestIterator)*3])
						+(m_originalPointSet[i*3+1]-m_originalPointSet[(*vectorNearestIterator)*3+1])*(m_originalPointSet[i*3+1]-m_originalPointSet[(*vectorNearestIterator)*3+1])
						+(m_originalPointSet[i*3+2]-m_originalPointSet[(*vectorNearestIterator)*3+2])*(m_originalPointSet[i*3+2]-m_originalPointSet[(*vectorNearestIterator)*3+2]);
					distance_normal=(m_filterNormals[i*3]-m_filterNormals[(*vectorNearestIterator)*3])*(m_filterNormals[i*3]-m_filterNormals[(*vectorNearestIterator)*3])
						+(m_filterNormals[i*3+1]-m_filterNormals[(*vectorNearestIterator)*3+1])*(m_filterNormals[i*3+1]-m_filterNormals[(*vectorNearestIterator)*3+1])
						+(m_filterNormals[i*3+2]-m_filterNormals[(*vectorNearestIterator)*3+2])*(m_filterNormals[i*3+2]-m_filterNormals[(*vectorNearestIterator)*3+2]);
					position_kernel=exp(-0.5*distance_position/m_meanShiftHposition);
					range_kernel=exp(-0.5*distance_normal/m_meanShiftHnormal);
					range_kernel*=position_kernel;
					if(!flag_stop[0][i]){
                        float gradient_kernel;
						distance_color_gradient[0]=(m_resultGradient[0][i*4]-m_originalColorGradients[0][(*vectorNearestIterator)*4])*(m_resultGradient[0][i*4]-m_originalColorGradients[0][(*vectorNearestIterator)*4])
							+(m_resultGradient[0][i*4+1]-m_originalColorGradients[0][(*vectorNearestIterator)*4+1])*(m_resultGradient[0][i*4+1]-m_originalColorGradients[0][(*vectorNearestIterator)*4+1])
							+(m_resultGradient[0][i*4+2]-m_originalColorGradients[0][(*vectorNearestIterator)*4+2])*(m_resultGradient[0][i*4+2]-m_originalColorGradients[0][(*vectorNearestIterator)*4+2])
							+(m_resultGradient[0][i*4+3]-m_originalColorGradients[0][(*vectorNearestIterator)*4+3])*(m_resultGradient[0][i*4+3]-m_originalColorGradients[0][(*vectorNearestIterator)*4+3]);
					    gradient_kernel=exp(-0.5*distance_color_gradient[0]/m_gradientKernelWide);
						sum_kernel[0]+=range_kernel*gradient_kernel;
						for(int j=0;j<4;j++){						
							tempGradient[0][j]+=m_originalColorGradients[0][(*vectorNearestIterator)*4+j]*gradient_kernel*range_kernel;
						}

					}
					if(!flag_stop[1][i]){
						float gradient_kernel;
						distance_color_gradient[1]=(m_resultGradient[1][i*4]-m_originalColorGradients[1][(*vectorNearestIterator)*4])*(m_resultGradient[1][i*4]-m_originalColorGradients[1][(*vectorNearestIterator)*4])
							+(m_resultGradient[1][i*4+1]-m_originalColorGradients[1][(*vectorNearestIterator)*4+1])*(m_resultGradient[1][i*4+1]-m_originalColorGradients[1][(*vectorNearestIterator)*4+1])
							+(m_resultGradient[1][i*4+2]-m_originalColorGradients[1][(*vectorNearestIterator)*4+2])*(m_resultGradient[1][i*4+2]-m_originalColorGradients[1][(*vectorNearestIterator)*4+2])
							+(m_resultGradient[1][i*4+3]-m_originalColorGradients[1][(*vectorNearestIterator)*4+3])*(m_resultGradient[1][i*4+3]-m_originalColorGradients[1][(*vectorNearestIterator)*4+3]);
						gradient_kernel=exp(-0.5*distance_color_gradient[1]/m_gradientKernelWide);
						sum_kernel[1]+=range_kernel*gradient_kernel;
						for(int j=0;j<4;j++){
							tempGradient[1][j]+=m_originalColorGradients[1][(*vectorNearestIterator)*4+j]*gradient_kernel*range_kernel;
						}
					}
					if(!flag_stop[2][i]){
						float gradient_kernel;
						distance_color_gradient[2]=(m_resultGradient[2][i*4]-m_originalColorGradients[2][(*vectorNearestIterator)*4])*(m_resultGradient[2][i*4]-m_originalColorGradients[2][(*vectorNearestIterator)*4])
							+(m_resultGradient[2][i*4+1]-m_originalColorGradients[2][(*vectorNearestIterator)*4+1])*(m_resultGradient[2][i*4+1]-m_originalColorGradients[2][(*vectorNearestIterator)*4+1])
							+(m_resultGradient[2][i*4+2]-m_originalColorGradients[2][(*vectorNearestIterator)*4+2])*(m_resultGradient[2][i*4+2]-m_originalColorGradients[2][(*vectorNearestIterator)*4+2])
							+(m_resultGradient[2][i*4+3]-m_originalColorGradients[2][(*vectorNearestIterator)*4+3])*(m_resultGradient[2][i*4+3]-m_originalColorGradients[2][(*vectorNearestIterator)*4+3]);
						gradient_kernel=exp(-0.5*distance_color_gradient[2]/m_gradientKernelWide);
						sum_kernel[2]+=range_kernel*gradient_kernel;
						for(int j=0;j<4;j++){
							tempGradient[2][j]+=m_originalColorGradients[2][(*vectorNearestIterator)*4+j]*gradient_kernel*range_kernel;
						}
					}
					
				}
				if(!flag_stop[0][i]){
					float gradient_kernel;
					distance_color_gradient[0]=(m_resultGradient[0][i*4]-m_originalColorGradients[0][i*4])*(m_resultGradient[0][i*4]-m_originalColorGradients[0][i*4])
						+(m_resultGradient[0][i*4+1]-m_originalColorGradients[0][i*4+1])*(m_resultGradient[0][i*4+1]-m_originalColorGradients[0][i*4+1])
						+(m_resultGradient[0][i*4+2]-m_originalColorGradients[0][i*4+2])*(m_resultGradient[0][i*4+2]-m_originalColorGradients[0][i*4+2])
						+(m_resultGradient[0][i*4+3]-m_originalColorGradients[0][i*4+3])*(m_resultGradient[0][i*4+3]-m_originalColorGradients[0][i*4+3]);
					gradient_kernel=exp(-0.5*distance_color_gradient[0]/m_gradientKernelWide);

					sum_kernel[0]+=gradient_kernel;
					for(int j=0;j<4;j++){
						tempGradient[0][j]+=m_originalColorGradients[0][i*4+j]*gradient_kernel;
						tempGradient[0][j]/=sum_kernel[0];
					}
					distance_color_gradient[0]=(m_resultGradient[0][0]-tempGradient[0][0])*(m_resultGradient[0][i*4]-tempGradient[0][0])
						+(m_resultGradient[0][i*4+1]-tempGradient[0][1])*(m_resultGradient[0][i*4+1]-tempGradient[0][1])
						+(m_resultGradient[0][i*4+2]-tempGradient[0][2])*(m_resultGradient[0][i*4+2]-tempGradient[0][2])
						+(m_resultGradient[0][i*4+3]-tempGradient[0][3])*(m_resultGradient[0][i*4+3]-tempGradient[0][3]);
					m_resultGradient[0][i*4]=tempGradient[0][0];
					m_resultGradient[0][i*4+1]=tempGradient[0][1];
					m_resultGradient[0][i*4+2]=tempGradient[0][2];
					m_resultGradient[0][i*4+3]=tempGradient[0][3];

					if(distance_color_gradient[0]<m_thresholdColorGradient){
						flag_stop[0][i]=true;
						num_stop[0]+=1;
					}
				}
			

				if(!flag_stop[1][i]){

					float gradient_kernel;
					distance_color_gradient[1]=(m_resultGradient[1][i*4]-m_originalColorGradients[1][i*4])*(m_resultGradient[1][i*4]-m_originalColorGradients[1][i*4])
						+(m_resultGradient[1][i*4+1]-m_originalColorGradients[1][i*4+1])*(m_resultGradient[1][i*4+1]-m_originalColorGradients[1][i*4+1])
						+(m_resultGradient[1][i*4+2]-m_originalColorGradients[1][i*4+2])*(m_resultGradient[1][i*4+2]-m_originalColorGradients[1][i*4+2])
						+(m_resultGradient[1][i*4+3]-m_originalColorGradients[1][i*4+3])*(m_resultGradient[1][i*4+3]-m_originalColorGradients[1][i*4+3]);
					gradient_kernel=exp(-0.5*distance_color_gradient[1]/m_gradientKernelWide);

					sum_kernel[1]+=gradient_kernel;
					for(int j=0;j<4;j++){
						tempGradient[1][j]+=m_originalColorGradients[1][i*4+j]*gradient_kernel;
						tempGradient[1][j]/=sum_kernel[1];
					}
					distance_color_gradient[1]=(m_resultGradient[1][i*4]-tempGradient[1][0])*(m_resultGradient[1][i*4]-tempGradient[1][0])
						+(m_resultGradient[1][i*4+1]-tempGradient[1][1])*(m_resultGradient[1][i*4+1]-tempGradient[1][1])
						+(m_resultGradient[1][i*4+2]-tempGradient[1][2])*(m_resultGradient[1][i*4+2]-tempGradient[1][2])
						+(m_resultGradient[1][i*4+3]-tempGradient[1][3])*(m_resultGradient[1][i*4+3]-tempGradient[1][3]);
					m_resultGradient[1][i*4]=tempGradient[1][0];
					m_resultGradient[1][i*4+1]=tempGradient[1][1];
					m_resultGradient[1][i*4+2]=tempGradient[1][2];
					m_resultGradient[1][i*4+3]=tempGradient[1][3];
					if(distance_color_gradient[1]<m_thresholdColorGradient){
						flag_stop[1][i]=true;
						num_stop[1]+=1;
					}				

				}

				if(!flag_stop[2][i]){

					float gradient_kernel;
					distance_color_gradient[2]=(m_resultGradient[2][i*4]-m_originalColorGradients[2][i*4])*(m_resultGradient[2][i*4]-m_originalColorGradients[2][i*4])
						+(m_resultGradient[2][i*4+1]-m_originalColorGradients[2][i*4+1])*(m_resultGradient[2][i*4+1]-m_originalColorGradients[2][i*4+1])
						+(m_resultGradient[2][i*4+2]-m_originalColorGradients[2][i*4+2])*(m_resultGradient[2][i*4+2]-m_originalColorGradients[2][i*4+2])
						+(m_resultGradient[2][i*4+3]-m_originalColorGradients[2][i*4+3])*(m_resultGradient[2][i*4+3]-m_originalColorGradients[2][i*4+3]);
					gradient_kernel=exp(-0.5*distance_color_gradient[2]/m_gradientKernelWide);

					sum_kernel[2]+=gradient_kernel;
					for(int j=0;j<4;j++){
						tempGradient[2][j]+=m_originalColorGradients[2][i*4+j]*gradient_kernel;
						tempGradient[2][j]/=sum_kernel[2];
					}
					distance_color_gradient[2]=(m_resultGradient[2][i*4]-tempGradient[2][0])*(m_resultGradient[2][i*4]-tempGradient[2][0])
						+(m_resultGradient[2][i*4+1]-tempGradient[2][1])*(m_resultGradient[2][i*4+1]-tempGradient[2][1])
						+(m_resultGradient[2][i*4+2]-tempGradient[2][2])*(m_resultGradient[2][i*4+2]-tempGradient[2][2])
						+(m_resultGradient[2][i*4+3]-tempGradient[2][3])*(m_resultGradient[2][i*4+3]-tempGradient[2][3]);
					m_resultGradient[2][i*4]=tempGradient[2][0];
					m_resultGradient[2][i*4+1]=tempGradient[2][1];
					m_resultGradient[2][i*4+2]=tempGradient[2][2];
					m_resultGradient[2][i*4+3]=tempGradient[2][3];

					if(distance_color_gradient[2]<m_thresholdColorGradient){
						flag_stop[2][i]=true;
						num_stop[2]+=1;
					}
				
				}
				
			}
		}

	}

	for(int i=0;i<3;i++){
		for(int j=0;j<m_numOfPoints*4;j++){
			m_originalColorGradients[i][j]=m_resultGradient[i][j];
		}
	}
	for (int i=0;i<3;i++){
		delete[] flag_stop[i];
	
		delete[] m_resultGradient[i];
	}
}
void FilterBilateral1:: EstimateNormalDirection()
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
		valueKeySeries.insert(std::map<float,int,lessFloat>::value_type(i+10,i));
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
		if((*multimapIterator).first>1){
			maxZpoint=(*multimapIterator).second;
			for(;multimapIterator!=valueKeySeries.end();multimapIterator++){
				if(m_originalPointSet[(*multimapIterator).second*3+2]>m_originalPointSet[maxZpoint*3+2])
					maxZpoint=(*multimapIterator).second;				
			}
			if(m_originalNormals[maxZpoint*3+2]<0){
				for(int i=0;i<3;i++){
					m_originalNormals[maxZpoint*3+i]=-m_originalNormals[maxZpoint*3+i];
				}
			}
			mapIterator=pointKeySeries.find(maxZpoint);
			(*mapIterator).second=0;
			multimapIterator=valueKeySeries.begin();
			for(;multimapIterator!=valueKeySeries.end();multimapIterator++){
				if((*multimapIterator).second==maxZpoint){

					valueKeySeries.erase(multimapIterator);
					valueKeySeries.insert(std::map<float,int>::value_type(0,maxZpoint));
					break;
				}
			}			

		}
		else{

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
					bool tempBool=false;
					multimapIterator=valueKeySeries.find((*mapIterator).second);
					if(neighborPoint==(*multimapIterator).second){
						tempBool=true;
					}

					std::multimap<float,int>::iterator  tempMultimapIteratorFront=multimapIterator;
					std::multimap<float,int>::iterator  tempMultimapIteratorBack=multimapIterator;

					while(!tempBool){
						if(tempMultimapIteratorFront!=valueKeySeries.begin()){
							tempMultimapIteratorFront--;
							if(neighborPoint==(*tempMultimapIteratorFront).second){
								tempBool=true;
								multimapIterator=tempMultimapIteratorFront;
							}

						}
						else{
							if(neighborPoint==(*tempMultimapIteratorFront).second){
								tempBool=true;
								multimapIterator=tempMultimapIteratorFront;
							}
						}

						if(tempMultimapIteratorBack!=valueKeySeries.end()){
							tempMultimapIteratorBack++;
							if(neighborPoint==(*tempMultimapIteratorBack).second){
								tempBool=true;
								multimapIterator=tempMultimapIteratorBack;
							}
						}
					}				

					if(tempBool==false){
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


void FilterBilateral1::CalculateBilateralFilterTen(int interactiveTimes)
{

	m_resultPointSet=new float[m_numOfPoints*3];

	/*
	 *	计算重心
	 */
	float* m_centerOfPoints=new float[m_numOfPoints*3];

	for(int i=0;i<m_numOfPoints;i++){
		std::map<int ,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		int numOfNearest=(*mapIterator).second.m_nearest.size();
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		float tempCenter[3];
		for(int j=0;j<3;j++){
			tempCenter[j]=m_originalPointSet[i*3+j];;
		}

		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
			for(int j=0;j<3;j++){
				tempCenter[j]+=m_originalPointSet[(*vectorIterator)*3+j];
			}
		}
		for(int j=0;j<3;j++){
			m_centerOfPoints[i*3+j]=tempCenter[j]/(numOfNearest+1);
		}
	}



	

	for(int www=0;www<interactiveTimes;www++)
	{
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

			float tempVector[3];//要滤波的点到邻域点的重心的向量
			float tempDistancePointToPoint; //要滤波的点到邻域点的重心的距离
			float tempDistancePointToPlane; //要滤波的点到邻域点到重心所在平面的距离
			float tempEstimationPosition[3]; //估计点的位置
			float kernel_distance;
			
			for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
				if((m_filterNormals[i*3]*m_filterNormals[(*vectorIterator)*3]
				+m_filterNormals[i*3+1]*m_filterNormals[(*vectorIterator)*3+1]
				+m_filterNormals[i*3+2]*m_filterNormals[(*vectorIterator)*3+2])<m_thresholdDistanceNormal)
					continue;			

				for(int j=0;j<3;j++){
					tempVector[j]=m_originalPointSet[i*3+j]-m_centerOfPoints[(*vectorIterator)*3+j];
				}

				tempDistancePointToPoint=sqrt(tempVector[0]*tempVector[0]
												+tempVector[1]*tempVector[1]
												+tempVector[2]*tempVector[2]);
				tempDistancePointToPlane=(tempVector[0]*m_filterNormals[(*vectorIterator)*3+0]
											+tempVector[1]*m_filterNormals[(*vectorIterator)*3+1]
											+tempVector[2]*m_filterNormals[(*vectorIterator)*3+2]);
				for(int j=0;j<3;j++){
					tempEstimationPosition[j]=m_originalPointSet[i*3+j]-tempDistancePointToPlane*m_filterNormals[(*vectorIterator)*3+j];
				}
				//tempDistancePointToPlane=abs(tempDistancePointToPlane);


				kernel_distance=exp(-0.5*tempDistancePointToPoint*tempDistancePointToPoint/m_pointVariation)*exp(-0.5*tempDistancePointToPlane*tempDistancePointToPlane/m_estimationVariation);
				sumKernel_position+=kernel_distance;
				for(int j=0;j<3;j++){
					m_resultPointSet[i*3+j]+=tempEstimationPosition[j]*kernel_distance;
				}
		
			//	vectorIterator++;

			}


			
			sumKernel_position+=1;
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]+=m_originalPointSet[i*3+j];
				m_resultPointSet[i*3+j]/=sumKernel_position;
			}

			/*
			 *	计算所滤波点的位置的重心
			 */
		/*	for(int j=0;j<3;j++){
				tempVector[j]=m_originalPointSet[i*3+j]-m_centerOfPoints[i*3+j];
			}

			tempDistancePointToPoint=sqrt(tempVector[0]*tempVector[0]
										+tempVector[1]*tempVector[1]
										+tempVector[2]*tempVector[2]);
			tempDistancePointToPlane=(tempVector[0]*m_filterNormals[i*3+0]
									+tempVector[1]*m_filterNormals[i*3+1]
									+tempVector[2]*m_filterNormals[i*3+2]);
			for(int j=0;j<3;j++){
				tempEstimationPosition[j]=m_originalPointSet[i*3+j]-tempDistancePointToPlane*m_filterNormals[i*3+j];
			}

			kernel_distance=exp(-0.5*tempDistancePointToPoint*tempDistancePointToPoint/m_pointVariation)*exp(-0.5*tempDistancePointToPlane*tempDistancePointToPlane/m_estimationVariation);
			sumKernel_position+=kernel_distance;
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]+=tempEstimationPosition[j]*kernel_distance;
				m_resultPointSet[i*3+j]/=sumKernel_position;
			}*/



			//////////////////////////////////////////////////////////////////////////
			// 调整滤波后的点在法向位置移动
			float tempdistance=0;
			tempdistance=((m_resultPointSet[i*3]-m_originalPointSet[i*3])*m_filterNormals[i*3]
			+(m_resultPointSet[i*3+1]-m_originalPointSet[i*3+1])*m_filterNormals[i*3+1]
			+(m_resultPointSet[i*3+2]-m_originalPointSet[i*3+2])*m_filterNormals[i*3+2]);
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]=m_originalPointSet[i*3+j]+tempdistance*m_filterNormals[i*3+j];
			}



	
		}
		for(int i=0;i<m_numOfPoints*3;i++){
			m_originalPointSet[i]=m_resultPointSet[i];
		}
	}	

	delete[] m_originalPointSet;
	delete[] m_centerOfPoints;

	m_originalPointSet=NULL;



	
}
void FilterBilateral1::CalculateBilateralFilterEleven(int interactiveTimes)
{

	m_resultPointSet=new float[m_numOfPoints*3];

	/*
	*	计算重心
	*/
	float* m_centerOfPoints=new float[m_numOfPoints*3];

	for(int i=0;i<m_numOfPoints;i++){
		std::map<int ,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		int numOfNearest=(*mapIterator).second.m_nearest.size();
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		float tempCenter[3];
		for(int j=0;j<3;j++){
			tempCenter[j]=m_originalPointSet[i*3+j];;
		}

		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
			for(int j=0;j<3;j++){
				tempCenter[j]+=m_originalPointSet[(*vectorIterator)*3+j];
			}
		}
		for(int j=0;j<3;j++){
			m_centerOfPoints[i*3+j]=tempCenter[j]/(numOfNearest+1);
		}
	}





	for(int www=0;www<interactiveTimes;www++)
	{
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

			float tempVector[3];//要滤波的点到邻域点的重心的向量
			float tempDistancePointToPoint; //要滤波的点到邻域点的重心的距离
			float tempDistancePointToPlane; //要滤波的点到邻域点到重心所在平面的距离
			float tempEstimationPosition[3]; //估计点的位置
			float kernel_distance;

			for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
				if((m_filterNormals[i*3]*m_filterNormals[(*vectorIterator)*3]
				+m_filterNormals[i*3+1]*m_filterNormals[(*vectorIterator)*3+1]
				+m_filterNormals[i*3+2]*m_filterNormals[(*vectorIterator)*3+2])<m_thresholdDistanceNormal)
					continue;			

				for(int j=0;j<3;j++){
					tempVector[j]=m_originalPointSet[i*3+j]-m_centerOfPoints[(*vectorIterator)*3+j];
				}

				tempDistancePointToPoint=sqrt(tempVector[0]*tempVector[0]
				+tempVector[1]*tempVector[1]
				+tempVector[2]*tempVector[2]);
				tempDistancePointToPlane=(tempVector[0]*m_filterNormals[(*vectorIterator)*3+0]
				+tempVector[1]*m_filterNormals[(*vectorIterator)*3+1]
				+tempVector[2]*m_filterNormals[(*vectorIterator)*3+2]);
				for(int j=0;j<3;j++){
					tempEstimationPosition[j]=m_originalPointSet[i*3+j]-tempDistancePointToPlane*m_filterNormals[(*vectorIterator)*3+j];
				}
				//tempDistancePointToPlane=abs(tempDistancePointToPlane);


				kernel_distance=exp(-0.5*tempDistancePointToPoint*tempDistancePointToPoint/m_pointVariation)*exp(-0.5*tempDistancePointToPlane*tempDistancePointToPlane/m_estimationVariation);
				sumKernel_position+=kernel_distance;
				for(int j=0;j<3;j++){
					m_resultPointSet[i*3+j]+=tempEstimationPosition[j]*kernel_distance;
				}

			//	vectorIterator++;

			}



		/*	sumKernel_position+=1;
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]+=m_originalPointSet[i*3+j];
				m_resultPointSet[i*3+j]/=sumKernel_position;
			}*/

			/*
			*	计算所滤波点的位置的重心
			*/
			for(int j=0;j<3;j++){
				tempVector[j]=m_originalPointSet[i*3+j]-m_centerOfPoints[i*3+j];
			}

			tempDistancePointToPoint=sqrt(tempVector[0]*tempVector[0]
			+tempVector[1]*tempVector[1]
			+tempVector[2]*tempVector[2]);
			tempDistancePointToPlane=(tempVector[0]*m_filterNormals[i*3+0]
			+tempVector[1]*m_filterNormals[i*3+1]
			+tempVector[2]*m_filterNormals[i*3+2]);
			for(int j=0;j<3;j++){
				tempEstimationPosition[j]=m_originalPointSet[i*3+j]-tempDistancePointToPlane*m_filterNormals[i*3+j];
			}

			kernel_distance=exp(-0.5*tempDistancePointToPoint*tempDistancePointToPoint/m_pointVariation)*exp(-0.5*tempDistancePointToPlane*tempDistancePointToPlane/m_estimationVariation);
			sumKernel_position+=kernel_distance;
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]+=tempEstimationPosition[j]*kernel_distance;
				m_resultPointSet[i*3+j]/=sumKernel_position;
			}



			//////////////////////////////////////////////////////////////////////////
			// 调整滤波后的点在法向位置移动
			float tempdistance=0;
			tempdistance=((m_resultPointSet[i*3]-m_originalPointSet[i*3])*m_filterNormals[i*3]
			+(m_resultPointSet[i*3+1]-m_originalPointSet[i*3+1])*m_filterNormals[i*3+1]
			+(m_resultPointSet[i*3+2]-m_originalPointSet[i*3+2])*m_filterNormals[i*3+2]);
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]=m_originalPointSet[i*3+j]+tempdistance*m_filterNormals[i*3+j];
			}




		}
		for(int i=0;i<m_numOfPoints*3;i++){
			m_originalPointSet[i]=m_resultPointSet[i];
		}
	}	

	delete[] m_originalPointSet;
	delete[] m_centerOfPoints;

	m_originalPointSet=NULL;




}
void FilterBilateral1::CalculateBilateralFilterTwelve(int interactiveTimes)
{

	m_resultPointSet=new float[m_numOfPoints*3];

	/*
	*	计算重心
	*/
	float* m_centerOfPoints=new float[m_numOfPoints*3];

	for(int i=0;i<m_numOfPoints;i++){
		std::map<int ,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		int numOfNearest=(*mapIterator).second.m_nearest.size();
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		float tempCenter[3];
		for(int j=0;j<3;j++){
			tempCenter[j]=m_originalPointSet[i*3+j];;
		}

		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
			for(int j=0;j<3;j++){
				tempCenter[j]+=m_originalPointSet[(*vectorIterator)*3+j];
			}
		}
		for(int j=0;j<3;j++){
			m_centerOfPoints[i*3+j]=tempCenter[j]/(numOfNearest+1);
		}
	}





	for(int www=0;www<interactiveTimes;www++)
	{
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

			float tempVector[3];//要滤波的点到邻域点的重心的向量
			float tempDistancePointToPoint; //要滤波的点到邻域点的重心的距离
			float tempDistancePointToPlane; //要滤波的点到邻域点到重心所在平面的距离
			float tempEstimationPosition[3]; //估计点的位置
			float kernel_distance;

			for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
				if((m_filterNormals[i*3]*m_filterNormals[(*vectorIterator)*3]
				+m_filterNormals[i*3+1]*m_filterNormals[(*vectorIterator)*3+1]
				+m_filterNormals[i*3+2]*m_filterNormals[(*vectorIterator)*3+2])<m_thresholdDistanceNormal)
					continue;			

				for(int j=0;j<3;j++){
					tempVector[j]=m_originalPointSet[i*3+j]-m_centerOfPoints[(*vectorIterator)*3+j];
				}

				tempDistancePointToPoint=sqrt(tempVector[0]*tempVector[0]
				+tempVector[1]*tempVector[1]
				+tempVector[2]*tempVector[2]);
				tempDistancePointToPlane=(tempVector[0]*m_filterNormals[(*vectorIterator)*3+0]
				+tempVector[1]*m_filterNormals[(*vectorIterator)*3+1]
				+tempVector[2]*m_filterNormals[(*vectorIterator)*3+2]);
				for(int j=0;j<3;j++){
					tempEstimationPosition[j]=m_originalPointSet[i*3+j]-tempDistancePointToPlane*m_filterNormals[(*vectorIterator)*3+j];
				}
				//tempDistancePointToPlane=abs(tempDistancePointToPlane);

				for(int j=0;j<3;j++){
					tempVector[j]=m_originalPointSet[i*3+j]-m_originalPointSet[(*vectorIterator)*3+j];

				}
				tempDistancePointToPoint=sqrt(tempVector[0]*tempVector[0]
				+tempVector[1]*tempVector[1]
				+tempVector[2]*tempVector[2]);



				kernel_distance=exp(-0.5*tempDistancePointToPoint*tempDistancePointToPoint/m_pointVariation)*exp(-0.5*tempDistancePointToPlane*tempDistancePointToPlane/m_estimationVariation);
				sumKernel_position+=kernel_distance;
				for(int j=0;j<3;j++){
					m_resultPointSet[i*3+j]+=tempEstimationPosition[j]*kernel_distance;
				}

				//	vectorIterator++;

			}



			/*	sumKernel_position+=1;
			for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]+=m_originalPointSet[i*3+j];
			m_resultPointSet[i*3+j]/=sumKernel_position;
			}*/

			/*
			*	计算所滤波点的位置的重心
			*/
			for(int j=0;j<3;j++){
				tempVector[j]=m_originalPointSet[i*3+j]-m_centerOfPoints[i*3+j];
			}

			tempDistancePointToPoint=sqrt(tempVector[0]*tempVector[0]
			+tempVector[1]*tempVector[1]
			+tempVector[2]*tempVector[2]);
			tempDistancePointToPlane=(tempVector[0]*m_filterNormals[i*3+0]
			+tempVector[1]*m_filterNormals[i*3+1]
			+tempVector[2]*m_filterNormals[i*3+2]);
			for(int j=0;j<3;j++){
				tempEstimationPosition[j]=m_originalPointSet[i*3+j]-tempDistancePointToPlane*m_filterNormals[i*3+j];
			}

			//kernel_distance=exp(-0.5*tempDistancePointToPoint*tempDistancePointToPoint/m_pointVariation)*exp(-0.5*tempDistancePointToPlane*tempDistancePointToPlane/m_estimationVariation);
			kernel_distance=exp(-0.5*tempDistancePointToPlane*tempDistancePointToPlane/m_estimationVariation);
			sumKernel_position+=kernel_distance;
			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]+=tempEstimationPosition[j]*kernel_distance;
				m_resultPointSet[i*3+j]/=sumKernel_position;
			}



			//////////////////////////////////////////////////////////////////////////
			// 调整滤波后的点在法向位置移动
			float tempdistance=0;
			tempdistance=((m_resultPointSet[i*3]-m_originalPointSet[i*3])*m_filterNormals[i*3]
			+(m_resultPointSet[i*3+1]-m_originalPointSet[i*3+1])*m_filterNormals[i*3+1]
			+(m_resultPointSet[i*3+2]-m_originalPointSet[i*3+2])*m_filterNormals[i*3+2]);
			for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]=m_originalPointSet[i*3+j]+tempdistance*m_filterNormals[i*3+j];
			}




		}
		for(int i=0;i<m_numOfPoints*3;i++){
			m_originalPointSet[i]=m_resultPointSet[i];
		}
	}	

	delete[] m_originalPointSet;
	delete[] m_centerOfPoints;

	m_originalPointSet=NULL;




}