#include "stdafx.h"
#include ".\nonshiftbilateralfilter.h"
#include "Mcube.h"
#include ".\mathlib\mathlib.h"
#include "Matrix.h"
#include <time.h>
//#include "GenerateNoise.h"

using namespace MATHLIB;

NonShiftBilateralFilter::NonShiftBilateralFilter(void)
{
	m_originalPointSet=NULL;
	m_originalColors=NULL;
	m_meanAreas=NULL;
	m_minPosition=0;
	m_maxPosition=0;
	m_normOriginalPointSetAndColor=NULL;
	m_numOfPoints=0;
	m_kNearest=0;
	m_IfFilterNormal=0;
	m_ifVolumePreserve=0;
	m_indexFunction=0;
	m_meanShiftStopNormal=0;
	m_originalNormals=NULL;
	m_originalNormalColor=NULL;	
	m_filterNormals=NULL;
	m_resultNormalColor=NULL;	
	m_resultPointSetAndColor=NULL;
	m_volumeChange=NULL;
	m_radius=0;
	m_interactiveTime=0;
	m_meanShiftHnormal=0;
	m_meanShiftHposition=0;
	m_variationTangent=0;
	m_variationNormal=0;	
	m_numDimension=0;
	m_numDimensionColor=0;
	m_resultPointSet=NULL;
	m_resultColor=NULL;
	m_meanRadius=NULL;
	m_meanCurvature=NULL;
	m_distanceMove=NULL;
}

NonShiftBilateralFilter::~NonShiftBilateralFilter(void)
{
	DeleteNonShiftBilateralFilter();

}

void NonShiftBilateralFilter::GetNonShiftBilateralFilter(int numOfPointSet,float* pointSet,int dimensionColor,float* colorSet,int kNearest, int ifFilterNormal, int ifPreserveVolume, int indexFunction,double thresholdNormal,double variationNormal,double variationTangent,int interactiveTime,double meanShiftHposition,double meanShiftHnormal,float variationColor,float thresholdFunction,float thresholdGradient,float thresholdDistanceNormal,float radius,float timeStep,int tangentOrManifold,int ifNormalWeight,int ifAreaWeight,int ifVariationNormal)
{
	//////////////////////////////////////////////////////////////////////////
	// 初始化数据
	long nearestStart,nearestFinish;
	long filterStart,filterFinish;
	m_numOfPoints=numOfPointSet;
	m_kNearest=kNearest;
	m_numDimension=3+dimensionColor;
	m_IfFilterNormal=ifFilterNormal;
	m_originalPointSet=new double[m_numOfPoints*3];
	m_ifVolumePreserve=ifPreserveVolume;
	m_indexFunction=indexFunction;
	m_meanShiftStopNormal=thresholdNormal;
	m_variationNormal=variationNormal*variationNormal;
	m_variationTangent=variationTangent*variationTangent;
	m_meanShiftHposition=meanShiftHposition*meanShiftHposition;
	m_meanShiftHnormal=meanShiftHnormal;
	m_thresholdGradient=thresholdGradient;
	m_radius=radius;
	m_timeStep=timeStep;
	m_ifVariationNormal=ifVariationNormal;

	m_tangentOrManifold=tangentOrManifold;
	m_ifNormalWeigtht=ifNormalWeight;
	m_ifAreaWeight=ifAreaWeight;

	m_interactiveTime=interactiveTime;
	m_variationColor=variationColor;
	m_thresholdFunction=thresholdFunction;
	m_thresholdDistanceNormal=thresholdDistanceNormal;
	for(int i=0;i<m_numOfPoints*3;i++){
		m_originalPointSet[i]=pointSet[i];
	}
	if(colorSet!=NULL){
		m_originalColors=new double[m_numOfPoints*dimensionColor];
		for(int i=0;i<m_numOfPoints*dimensionColor;i++){
			m_originalColors[i]=colorSet[i];
		}
	}
	nearestStart=clock();
	ComputeMapKnearest();
	nearestFinish=clock();
	nearestTime=nearestFinish-nearestStart;
	filterStart=clock();
	
	if(m_indexFunction==1){
		ComputeFilterOne(m_interactiveTime);

	}
	if(m_indexFunction==2){
		ComputeFilterTwo(m_interactiveTime);
	}
	if(m_indexFunction==3){
		ComputeFilterThree(m_interactiveTime);
	}
	if(m_indexFunction==4){
		ComputeFilterFour(m_interactiveTime);
	}
	if(m_indexFunction==0){
		if(m_resultPointSet!=NULL)
			delete[] m_resultPointSet;
		m_resultPointSet=new double[m_numOfPoints*3];
		for(int i=0;i<m_numOfPoints*3;i++){
			m_resultPointSet[i]=m_originalPointSet[i];
		}
	}
	if(m_indexFunction==5){
		ComputeFilterFive(m_interactiveTime);
	}
	if(m_indexFunction==6){
		ComputeFilterSix(m_interactiveTime);
	}
	if(m_indexFunction==7){
		ComputeFilterSeven(m_interactiveTime);
	}

	if(m_indexFunction==10){
		ComputeFilterTen(m_interactiveTime);
	}
	if(m_indexFunction==21){
		ComputeFilterTwentyone(m_interactiveTime);
	}
	if(m_indexFunction==22){
		ComputeFilterTwentytwo(m_interactiveTime);
	}
	if(m_indexFunction==23){
		ComputeFilterTwentyThree(m_interactiveTime);
	}
	if(m_indexFunction==24){
		ComputeFilterTwentyFour(m_interactiveTime);
	}
	if(m_indexFunction==25){
		ComputeFilterTwentyFive(m_interactiveTime);
	}

	if(m_indexFunction==26){
		ComputeFilterTwentySix(m_interactiveTime);
	}

	if(m_indexFunction==30){
		ComputeFilterThirty(m_interactiveTime);
	}

	ComputeNormal();
	FilterNormal();
	for(int i=0;i<m_numOfPoints*3;i++){
		m_originalNormals[i]=m_filterNormals[i];

	}
	filterFinish=clock();
	filterTime=filterFinish-filterStart;
	return;

}
void NonShiftBilateralFilter:: DeleteNonShiftBilateralFilter()
{
	if(m_originalPointSet!=NULL)
		delete[] m_originalPointSet;
	m_originalPointSet=NULL;
	if(m_originalColors!=NULL)
		delete[] m_originalColors;
	m_originalColors=NULL;
	if(m_normOriginalPointSetAndColor!=NULL)
		delete[] m_normOriginalPointSetAndColor;
	m_normOriginalPointSetAndColor=NULL;

	if(m_originalNormals!=NULL)
		delete[] m_originalNormals;
	m_originalNormals=NULL;
	if(m_originalNormalColor!=NULL)
		delete[] m_originalNormalColor;
	m_originalNormalColor=NULL;	
	if(m_filterNormals!=NULL)
		delete[] m_filterNormals;
	m_filterNormals=NULL;
	if(m_resultNormalColor!=NULL)
		delete[] m_resultNormalColor;
	m_resultNormalColor=NULL;	
	if(m_resultPointSetAndColor!=NULL)
		delete[] m_resultPointSetAndColor;
	m_resultPointSetAndColor=NULL;
	if(m_volumeChange!=NULL)
		delete[] m_volumeChange;
	m_volumeChange=NULL;
	if(m_resultPointSet!=NULL)
		delete[] m_resultPointSet;
	m_resultPointSet=NULL;
	if(m_resultColor!=NULL)
		delete[] m_resultColor;
	m_resultColor=NULL;
	if(m_meanRadius!=NULL)
		delete[] m_meanRadius;
	m_meanRadius=NULL;
	if(m_meanAreas!=NULL)
		delete[] m_meanAreas;
	m_meanAreas=NULL;
	if(m_distanceMove!=NULL)
		delete[] m_distanceMove;
	m_distanceMove=NULL;
	if(m_meanCurvature!=NULL)
		delete[] m_meanCurvature;
	m_meanCurvature=NULL;

}
void NonShiftBilateralFilter:: ComputeMapKnearest()
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
	CString filename_pw = "D:\\qin\\sourcecode\\NonShiftknearest.txt";
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
void NonShiftBilateralFilter::ComputeNormal()
{
	double centroidPosition[3];//重心位置
	double* localVariationVector;//重心到点的位置的向量
	int numOfNearest;
	
	mat_f8 mat_covaMatrix(3, 3);
	if(m_originalNormals!=NULL){
		delete[] m_originalNormals;
		m_originalNormals=NULL;
	}	
	m_originalNormals=new double[m_numOfPoints*3];
	if(m_meanAreas!=NULL){
		delete[] m_meanAreas;
		m_meanAreas=NULL;
	}
	m_meanAreas=new double[m_numOfPoints];
	if(m_meanCurvature!=NULL){
		delete[] m_meanCurvature;
		m_meanCurvature=NULL;
	}
	m_meanCurvature=new double[m_numOfPoints];
	if(m_resultColor!=NULL){
		delete[] m_resultColor;
	    m_resultColor=NULL;
	}
	m_resultColor=new double[m_numOfPoints*3];


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
		m_meanAreas[i]=abs(eval(0)*eval(1));
		m_meanCurvature[i]=eval(2)/(eval(0)+eval(1)+eval(2));
	
		if(m_meanCurvature[i]*4<=0.005){
			m_resultColor[i*3]=0;
			m_resultColor[i*3+1]=0;
			m_resultColor[i*3+2]=0.5+80*m_meanCurvature[i]*4;
		}
		if((m_meanCurvature[i]*4>0.005)&(m_meanCurvature[i]*4<=0.1)){
			m_resultColor[i*3]=0;
			m_resultColor[i*3+1]=0.15+6*(m_meanCurvature[i]*4-0.005);
			m_resultColor[i*3+2]=0;
		}
		if(m_meanCurvature[i]*4>0.1){
			m_resultColor[i*3]=0.10+10*(m_meanCurvature[i]*4-0.1);
			m_resultColor[i*3+1]=0;
			m_resultColor[i*3+2]=0;
		}


		delete[] localVariationVector;
	}
	EstimateNormalDirection();

	//CString filename_pw = "D:\\qin\\experiment\\normal.txt";
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

	//filename_pw = "D:\\qin\\experiment\\NonShiftmeanCurvature.txt";

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
void NonShiftBilateralFilter::ComputeNormalColor()
{
	double* centroidPosition;//重心位置
	centroidPosition=new double[m_numDimension];
	double* localVariationVector;//重心到点的位置的向量
	int numOfNearest;

	mat_f8 mat_covaMatrix(m_numDimension, m_numDimension);
	if(m_originalNormalColor!=NULL)
		delete[] m_originalNormalColor;
	m_originalNormalColor=new double[m_numOfPoints*m_numDimension];

	for(int i=0;i<m_numOfPoints;i++){
		std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		numOfNearest=(*mapIterator).second.m_numOfNearest;
		for(int j=0;j<m_numDimension;j++){
			centroidPosition[j]=m_normOriginalPointSetAndColor[i*m_numDimension+j];
		}
		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
			for(int j=0;j<m_numDimension;j++){
				centroidPosition[j]+=m_normOriginalPointSetAndColor[(*vectorIterator)*m_numDimension+j];

			}
		}
		for(int j=0;j<m_numDimension;j++){
			centroidPosition[j]/=(numOfNearest+1);
		}
		localVariationVector=new double[(numOfNearest+1)*m_numDimension];
		for(int j=0;j<m_numDimension;j++){
			localVariationVector[j]=m_normOriginalPointSetAndColor[i*m_numDimension+j]-centroidPosition[j];
		}
		vectorIterator=(*mapIterator).second.m_nearest.begin();
		for(int j=1;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,j++){
			for(int k=0;k<m_numDimension;k++){
				localVariationVector[j*m_numDimension+k]=m_normOriginalPointSetAndColor[(*vectorIterator)*3]-centroidPosition[k];
			}
		}
		//求方差矩阵
		for(int j=0;j<m_numDimension;j++){
			for(int k=0;k<m_numDimension;k++){
				mat_covaMatrix(j,k)=0;
				for(int m=0;m<(numOfNearest+1);m++){
					mat_covaMatrix(j,k)=mat_covaMatrix(j,k)+localVariationVector[m*m_numDimension+j]*localVariationVector[m*m_numDimension+k];
				}
			}
		}
		vec_f8 eval(m_numDimension);
		mat_f8 evec(m_numDimension, m_numDimension);
		eigen_symm(eval, evec, mat_covaMatrix);
		for(int j=0;j<m_numDimension;j++){
			m_originalNormalColor[i*m_numDimension+j]=evec(j,m_numDimension-1);
		}
		
		delete[] localVariationVector;
	}
}
void NonShiftBilateralFilter::FilterNormal()
{
	double distance_position=0;
	double distance_normal=0;
	double position_kernel=0;
	double range_kernel=0;
	double sum_kernel=0;
	double length_normal_vector=0;
	int num_stop=0;//已经停止的数

	bool* flag_stop=new bool[m_numOfPoints];
	for(int i=0;i<m_numOfPoints;i++){
		flag_stop[i]=false;
	}
	if(m_filterNormals!=NULL){
		delete[] m_filterNormals;
		m_filterNormals=NULL;
	}
	m_filterNormals=new double[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints*3;i++){
		m_filterNormals[i]=m_originalNormals[i];
	}
	while(num_stop<m_numOfPoints){
		for(int i=0;i<m_numOfPoints;i++){
			if(!flag_stop[i]){
				sum_kernel=0;
				double tempNormal[3];
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
					//tempNormal[j]/=length_normal_vector;
					m_filterNormals[i*3+j]=tempNormal[j];
				}

			}
		}

	}
	//delete[] m_originalNormals;
	//m_originalNormals=NULL;
	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonShiftfilterNormal.txt";
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
}
void NonShiftBilateralFilter::FilterNormalColor()
{
	double distance_position;
	double distance_normal;
	double position_kernel,range_kernel;
	double sum_kernel;
	double length_normal_vector;
	int num_stop=0;//已经停止的数

	double* tempNormal=new double[m_numDimension];
	
	bool* flag_stop=new bool[m_numOfPoints];
	for(int i=0;i<m_numOfPoints;i++){
		flag_stop[i]=false;
	}
	if(m_resultNormalColor!=NULL)
		delete[] m_resultNormalColor;
	m_resultNormalColor=new double[m_numOfPoints*m_numDimension];

	for(int i=0;i<m_numOfPoints*m_numDimension;i++){
		m_resultNormalColor[i]=m_originalNormalColor[i];
	}
	while(num_stop<m_numOfPoints){
		for(int i=0;i<m_numOfPoints;i++){
			if(!flag_stop[i]){
				sum_kernel=0;
				for(int j=0;j<m_numDimension;j++){
					tempNormal[j]=0;
				}			
				std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(i);
				//int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
				std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
				
				for(;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++){
					distance_position=0;
					distance_normal=0;
			
					for(int j=0;j<m_numDimension;j++){
						distance_position+=(m_normOriginalPointSetAndColor[i*m_numDimension+j]-m_normOriginalPointSetAndColor[(*vectorNearestIterator)*m_numDimension+j])*(m_normOriginalPointSetAndColor[i*m_numDimension+j]-m_normOriginalPointSetAndColor[(*vectorNearestIterator)*m_numDimension+j]);
						distance_normal+=(m_resultNormalColor[i*m_numDimension+j]-m_originalNormalColor[(*vectorNearestIterator)*m_numDimension+j])*(m_resultNormalColor[i*m_numDimension+j]-m_originalNormalColor[(*vectorNearestIterator)*m_numDimension+j]);
					}
				    position_kernel=exp(-0.5*distance_position/m_meanShiftHposition);
					range_kernel=exp(-0.5*distance_normal/m_meanShiftHnormal);
					sum_kernel+=position_kernel*range_kernel;

					for(int j=0;j<m_numDimension;j++){
						tempNormal[j]+=m_originalNormalColor[(*vectorNearestIterator)*m_numDimension+j]*position_kernel*range_kernel;
					}
				}
				distance_normal=0;
				for(int j=0;j<m_numDimension;j++){
					distance_normal+=(m_resultNormalColor[i*m_numDimension+j]-m_originalNormalColor[i*m_numDimension+j])
						*(m_resultNormalColor[i*m_numDimension+j]-m_originalNormalColor[i*m_numDimension+j]);
				}
			
				range_kernel=exp(-0.5*distance_normal/m_meanShiftHnormal);
				sum_kernel+=range_kernel;

				for(int j=0;j<m_numDimension;j++){
					tempNormal[j]+=m_originalNormalColor[i*3+j]*range_kernel;
					tempNormal[j]/=sum_kernel;
				}
				length_normal_vector=0;
				for(int j=0;j<m_numDimension;j++){
					length_normal_vector+=tempNormal[j]*tempNormal[j];
				}
				length_normal_vector=sqrt(length_normal_vector);
				for(int j=0;j<m_numDimension;j++){
					tempNormal[j]/=length_normal_vector;
					//m_resultNormalColor[i*m_numDimension+j]=tempNormal[j];
				}
				double length_normal=0;
				for(int j=0;j<m_numDimension;j++){
					length_normal+=(m_resultNormalColor[i*m_numDimension+j]-tempNormal[j])
						*(m_resultNormalColor[i*m_numDimension+j]-tempNormal[j]);
					m_resultNormalColor[i*m_numDimension+j]=tempNormal[j];
				}


				if(length_normal<m_meanShiftStopNormal){
						flag_stop[i]=true;
						num_stop+=1;
					}				
			}
		}

	}
	delete[] m_originalNormalColor;
	m_originalNormalColor=NULL;
}
void NonShiftBilateralFilter::NormPositionAndColor()
{
	m_minPosition=(double)m_originalPointSet[0];
	m_maxPosition=(double)m_originalPointSet[0];
	for(int i=1;i<m_numOfPoints*3;i++){
		if(m_minPosition>m_originalPointSet[i])
			m_minPosition=(double)m_originalPointSet[i];
		else
			if(m_maxPosition<m_originalPointSet[i])
				m_maxPosition=(double)m_originalPointSet[i];
	}
	m_normOriginalPointSetAndColor=new double[m_numOfPoints*m_numDimension];
	double gradient=double(1/(double)(m_maxPosition-m_minPosition));
	double cutConstant=double(m_minPosition/(double)(m_minPosition-m_maxPosition));
	for(int i=0;i<m_numOfPoints;i++){
		for(int j=0;j<3;j++){
			m_normOriginalPointSetAndColor[i*m_numDimension+j]=gradient*(double)(m_originalPointSet[i*3+j])+m_minPosition/(m_minPosition-m_maxPosition);
		}
		for(int j=3;j<m_numDimension;j++){
			m_normOriginalPointSetAndColor[i*m_numDimension+j]=(double)m_originalColors[i*3+j-3];
		}		
	}
	return;

}

//////////////////////////////////////////////////////////////////////////
// 仅仅计算位置，且没有volume 保护
void NonShiftBilateralFilter::ComputeFilterOne(int interractiveTime)
{

    m_volumeChange=new double[m_numOfPoints];
	m_resultPointSet=new double[m_numOfPoints*3];
	m_distanceMove=new double[m_numOfPoints];
	ComputeNormal();
	/*
	*	对法向进行滤波
	*/
	if(m_IfFilterNormal==1){
		FilterNormal();
	}
	else{
		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
			m_filterNormals[j]=m_originalNormals[j];
		}
	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;
	
	for(int i=0;i<interractiveTime;i++){
	
		double* distanceProjection;
		double* distanceTangent;
		double* weightTangent;
		double* weightNormal;
		double weightSum;
		double radiusMean;
		double distanceMove;

		/*
		 *	对位置进行滤波
		 */
		for(int j=0;j<m_numOfPoints;j++){
		
			double localNormal[3];
			for(int k=0;k<3;k++){
				localNormal[k]=m_filterNormals[j*3+k];
			}
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(j);
			int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
			std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
			distanceProjection=new double[numOfNearest];
			distanceTangent=new double[numOfNearest];
			weightTangent=new double[numOfNearest];
			weightNormal=new double[numOfNearest];
			int* outFlag=new int[numOfNearest];
			for(int k=0;k<numOfNearest;k++){
				distanceProjection[k]=0.0;
				distanceTangent[k]=0.0;
				weightTangent[k]=0.0;
				weightNormal[k]=0.0;
				outFlag[k]=0;
			}
		
			weightSum=0;
			distanceMove=0;
			double sumDistanceProjection=0;
			double sumDistanceTangent=0;
			//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonDistance.txt";
			//FILE *fpout;
			//if((fpout = fopen(filename_pw, "a")) == NULL)
			//{
			//	int dkjkd;
			//	//MessageBox("can't open the file!");
			//}
			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){			
                /*
                 *	法向的threshold处理
                 */
				if((m_filterNormals[i*3]*m_filterNormals[(*vectorNearestIterator)*3]
					+m_filterNormals[i*3+1]*m_filterNormals[(*vectorNearestIterator)*3+1]
					+m_filterNormals[i*3+2]*m_filterNormals[(*vectorNearestIterator)*3+2])<m_thresholdDistanceNormal){
						outFlag[k]=0;						
						numOfNearest-=1;
						continue;
					}
				else{
					outFlag[k]=1;
				}
				double localVectorPosition[3];
				for(int w=0;w<3;w++){
					localVectorPosition[w]=m_originalPointSet[(*vectorNearestIterator)*3+w]-m_originalPointSet[j*3+w];					
				}
				distanceProjection[k]=Vector3Vector(localVectorPosition,localNormal);
				distanceTangent[k]=ComputeNormOfVector3(localVectorPosition);
				distanceTangent[k]=sqrt(distanceTangent[k]*distanceTangent[k]-distanceProjection[k]*distanceProjection[k]);
				sumDistanceTangent+=distanceTangent[k];
				sumDistanceProjection+=abs(distanceProjection[k]);
				//weightNormal[k]=exp(-0.5*distanceProjection*distanceProjection/m_variationNormal);
				//weightTangent=exp(-0.5*distanceTangent*distanceTangent/m_variationTangent);
				//double tempWeight=weightNormal*weightTangent;
				//distanceMove+=distanceProjection*tempWeight;
				//weightSum+=tempWeight;
			
			
			
				//		fprintf(fpout," %f ",distanceProjection[k]);	
		
			
			}
			//fprintf(fpout,"\n");
		
		//	fclose(fpout);
			sumDistanceTangent/=(numOfNearest+1);
			sumDistanceProjection/=(numOfNearest+1);
			double localVariationNormal,localVarationTangent;
			localVariationNormal=m_variationNormal*sumDistanceProjection*sumDistanceProjection;
			if(localVariationNormal<0.000001)
				localVariationNormal=0.000001;
			localVarationTangent=m_variationTangent*sumDistanceTangent*sumDistanceTangent;
			if(localVarationTangent<0.000001)
				localVarationTangent=0.000001;
			vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();

		
			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if(outFlag[k]==0)
					continue;
				weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
				weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);				
				double tempWeight=weightNormal[k]*weightTangent[k];
				distanceMove+=distanceProjection[k]*tempWeight;
				weightSum+=tempWeight;			
			}

			distanceMove=distanceMove/(weightSum+1);
			m_distanceMove[j]=distanceMove;
			m_volumeChange[j]=distanceMove*m_meanAreas[j];
			for(int k=0;k<3;k++){
				m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];
			}
			delete distanceProjection;
			delete distanceTangent;
			delete weightTangent;
			delete weightNormal;
			delete outFlag;
			
		
	    
		}
	//	RestoreVolume();
		for(int j=0;j<m_numOfPoints*3;j++){
			m_originalPointSet[j]=m_resultPointSet[j];
		}
		
	}
	//	m_mapKnearest.clear();

	//delete[] m_originalPointSet;
	//m_originalPointSet=NULL;
	return;
	
}
void NonShiftBilateralFilter::ComputeFilterTwo(int interractiveTime)
{
	//ComputeNormal();
	//if(m_IfFilterNormal==1){
	//	FilterNormal();
	//}
	//else{

	//	m_filterNormals=new double[m_numOfPoints*3];
	//	for(int i=0;i<m_numOfPoints*3;i++){
	//		m_filterNormals[i]=m_originalNormals[i];
	//	}
	//}
	//delete[] m_originalNormals;
	//m_originalNormals=NULL;
	m_volumeChange=new double[m_numOfPoints];
	m_resultPointSet=new double[m_numOfPoints*3];
	m_distanceMove=new double[m_numOfPoints];
	ComputeNormal();
	if(m_IfFilterNormal==1){
		FilterNormal();
	}
	else{

		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
			m_filterNormals[j]=m_originalNormals[j];
		}
	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;
	
	for(int i=0;i<interractiveTime;i++){
		
		/*ComputeNormal();
		if(m_IfFilterNormal==1){
			FilterNormal();
		}
		else{

			m_filterNormals=new double[m_numOfPoints*3];
			for(int j=0;j<m_numOfPoints*3;j++){
				m_filterNormals[j]=m_originalNormals[j];
			}
		}
		delete[] m_originalNormals;
		m_originalNormals=NULL;*/
	
		double* distanceProjection;
		double* distanceTangent;
		double* weightTangent;
		double* weightNormal;
		double weightSum;
		double radiusMean;
		double distanceMove;
		for(int j=0;j<m_numOfPoints;j++){
			double localNormal[3];
			for(int k=0;k<3;k++){
				localNormal[k]=m_filterNormals[j*3+k];
			}
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(j);
			int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
			std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
			distanceProjection=new double[numOfNearest];
			distanceTangent=new double[numOfNearest];
			weightTangent=new double[numOfNearest];
			weightNormal=new double[numOfNearest];
			int* outFlag=new int[numOfNearest];
			for(int k=0;k<numOfNearest;k++){
				distanceProjection[k]=0.0;
				distanceTangent[k]=0.0;
				weightTangent[k]=0.0;
				weightNormal[k]=0.0;
				outFlag[k]=0;
			}
			weightSum=0;
			distanceMove=0;
			double sumDistanceProjection=0;
			double sumDistanceTangent=0;
			//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonDistance.txt";
			//FILE *fpout;
			//if((fpout = fopen(filename_pw, "a")) == NULL)
			//{
			//	int dkjkd;
			//	//MessageBox("can't open the file!");
			//}
			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if((m_filterNormals[i*3]*m_filterNormals[(*vectorNearestIterator)*3]
					+m_filterNormals[i*3+1]*m_filterNormals[(*vectorNearestIterator)*3+1]
					+m_filterNormals[i*3+2]*m_filterNormals[(*vectorNearestIterator)*3+2])<m_thresholdDistanceNormal){
						outFlag[k]=0;
						
						numOfNearest-=1;continue;
					}
				else{
					outFlag[k]=1;
				}
				double localVectorPosition[3];				
				for(int w=0;w<3;w++){
					localVectorPosition[w]=m_originalPointSet[(*vectorNearestIterator)*3+w]-m_originalPointSet[j*3+w];

				}
				distanceProjection[k]=Vector3Vector(localVectorPosition,localNormal);
				distanceTangent[k]=ComputeNormOfVector3(localVectorPosition);
				distanceTangent[k]=sqrt(distanceTangent[k]*distanceTangent[k]-distanceProjection[k]*distanceProjection[k]);
				sumDistanceTangent+=distanceTangent[k];
				sumDistanceProjection+=abs(distanceProjection[k]);
				//weightNormal[k]=exp(-0.5*distanceProjection*distanceProjection/m_variationNormal);
				//weightTangent=exp(-0.5*distanceTangent*distanceTangent/m_variationTangent);
				//double tempWeight=weightNormal*weightTangent;
				//distanceMove+=distanceProjection*tempWeight;
				//weightSum+=tempWeight;



		//		fprintf(fpout," %f ",distanceProjection[k]);




			}
//			fprintf(fpout,"\n");

		//	fclose(fpout);
			sumDistanceTangent/=(numOfNearest+1);
			sumDistanceProjection/=(numOfNearest+1);
			double localVariationNormal,localVarationTangent;
			localVariationNormal=m_variationNormal*sumDistanceProjection*sumDistanceProjection;
			if(localVariationNormal<0.000001)
				localVariationNormal=0.000001;
			localVarationTangent=m_variationTangent*sumDistanceTangent*sumDistanceTangent;
			if(localVarationTangent<0.000001)
				localVarationTangent=0.000001;
			vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();


			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if(outFlag[k]==0)
					continue;
				weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
				weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);
				

				double tempWeight=weightNormal[k]*weightTangent[k]*m_meanAreas[(*vectorNearestIterator)];
				/*tempWeight*=abs(m_filterNormals[(*vectorNearestIterator)*3]*m_filterNormals[j*3]
						+m_filterNormals[(*vectorNearestIterator)*3+1]*m_filterNormals[j*3+1]
						+m_filterNormals[(*vectorNearestIterator)*3+2]*m_filterNormals[j*3+2]);*/
				distanceMove+=distanceProjection[k]*tempWeight;
				weightSum+=tempWeight;
			}

			distanceMove=distanceMove/(weightSum+m_meanAreas[j]);
			m_distanceMove[j]=distanceMove;
			m_volumeChange[j]=distanceMove*m_meanAreas[j];
			for(int k=0;k<3;k++){
				m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];

			}
			delete distanceProjection;
			delete distanceTangent;
			delete weightTangent;
			delete weightNormal;
			delete outFlag;

		}
	//	RestoreVolume();
		for(int j=0;j<m_numOfPoints*3;j++){
			m_originalPointSet[j]=m_resultPointSet[j];
		}
	}
	//delete[] m_originalPointSet;
	//m_originalPointSet=NULL;
	return;
}
void NonShiftBilateralFilter::ComputeFilterThree(int interractiveTime)
{
	m_volumeChange=new double[m_numOfPoints];
	m_resultPointSet=new double[m_numOfPoints*3];
	m_resultColor=new double[m_numOfPoints*3];
	m_distanceMove=new double[m_numOfPoints];
	ComputeNormal();
	if(m_IfFilterNormal==1){
		FilterNormal();
	}
	else{
		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
			m_filterNormals[j]=m_originalNormals[j];
		}
	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;

	for(int i=0;i<interractiveTime;i++){
		//ComputeNormal();
		//if(m_IfFilterNormal==1){
		//	FilterNormal();
		//}
		//else{
		//	m_filterNormals=new double[m_numOfPoints*3];
		//	for(int j=0;j<m_numOfPoints*3;j++){
		//		m_filterNormals[j]=m_originalNormals[j];
		//	}
		//}
		//delete[] m_originalNormals;
		//m_originalNormals=NULL;
		
		double* distanceProjection;
		double* distanceTangent;
		double* weightTangent;
		double* weightNormal;
		double* weightColor[3];
		double* distanceColor[3];
		double weightSum;
		double weightSumColor[3];
		double radiusMean;
		double distanceMove;
		double distanceMoveColor[3];
		for(int j=0;j<m_numOfPoints;j++){
		
			double localNormal[3];
			for(int k=0;k<3;k++){
				localNormal[k]=m_filterNormals[j*3+k];
			}
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(j);
			int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
			std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
			distanceProjection=new double[numOfNearest];
			distanceTangent=new double[numOfNearest];			
			weightTangent=new double[numOfNearest];
			weightNormal=new double[numOfNearest];
			int* outFlag=new int[numOfNearest];
			for(int k=0;k<numOfNearest;k++){
				distanceProjection[k]=0.0;
				distanceTangent[k]=0.0;
				weightTangent[k]=0.0;
				weightNormal[k]=0.0;
				outFlag[k]=0;
			}
			for(int k=0;k<3;k++){
				distanceColor[k]=new double[numOfNearest];
				weightColor[k]=new double[numOfNearest];
				weightSumColor[k]=0;
				distanceMoveColor[k]=0;
			
			}
			weightSum=0;
			distanceMove=0;
			double sumDistanceProjection=0;
			double sumDistanceTangent=0;
			
			//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonDistance.txt";
			//FILE *fpout;
			//if((fpout = fopen(filename_pw, "a")) == NULL)
			//{
			//	int dkjkd;
			//	//MessageBox("can't open the file!");
			//}
			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if((m_filterNormals[(*vectorNearestIterator)*3]*localNormal[0]
				+m_filterNormals[(*vectorNearestIterator)*3+1]*localNormal[1]
				+m_filterNormals[(*vectorNearestIterator)*3+2]*localNormal[2])<m_thresholdDistanceNormal){
					outFlag[k]=0;
					numOfNearest-=1;
					continue;
				}
				double localVectorPosition[3];				
				for(int w=0;w<3;w++){
					localVectorPosition[w]=m_originalPointSet[(*vectorNearestIterator)*3+w]-m_originalPointSet[j*3+w];
					distanceColor[w][k]=m_originalColors[(*vectorNearestIterator)*3+w]-m_originalColors[j*3+w];

				}
				distanceProjection[k]=Vector3Vector(localVectorPosition,localNormal);
				distanceTangent[k]=ComputeNormOfVector3(localVectorPosition);
				distanceTangent[k]=sqrt(distanceTangent[k]*distanceTangent[k]-distanceProjection[k]*distanceProjection[k]);
				sumDistanceTangent+=distanceTangent[k];
				sumDistanceProjection+=abs(distanceProjection[k]);		
				
				//weightNormal[k]=exp(-0.5*distanceProjection*distanceProjection/m_variationNormal);
				//weightTangent=exp(-0.5*distanceTangent*distanceTangent/m_variationTangent);
				//double tempWeight=weightNormal*weightTangent;
				//distanceMove+=distanceProjection*tempWeight;
				//weightSum+=tempWeight;



				//		fprintf(fpout," %f ",distanceProjection[k]);
			}
			//fprintf(fpout,"\n");

			//	fclose(fpout);
			sumDistanceTangent/=(numOfNearest+1);
			sumDistanceProjection/=(numOfNearest+1);
			double localVariationNormal,localVarationTangent;
			localVariationNormal=m_variationNormal*m_variationNormal*sumDistanceProjection*sumDistanceProjection;
			if(localVariationNormal<0.000001)
				localVariationNormal=0.000001;
			localVarationTangent=m_variationTangent*m_variationTangent*sumDistanceTangent*sumDistanceTangent;
			if(localVarationTangent<0.000001)
				localVarationTangent=0.000001;
			vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();


			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if(outFlag[0]==0){
					continue;
				}				
			
				weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
				weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);
				for(int w=0;w<3;w++){
					if(abs(distanceColor[w][k])>m_thresholdFunction)
						weightColor[w][k]=0;
					else
                        weightColor[w][k]=exp(-0.5*distanceColor[w][k]*distanceColor[w][k]/m_variationColor);
					
				}
				double tempWeight=weightNormal[k]*weightTangent[k];
				double tempWeightColor[3];
				tempWeightColor[0]=tempWeight*weightColor[0][k];
				tempWeightColor[1]=tempWeight*weightColor[1][k];
				tempWeightColor[2]=tempWeight*weightColor[2][k];
				distanceMove+=distanceProjection[k]*tempWeight;
				distanceMoveColor[0]+=distanceColor[0][k]*tempWeightColor[0];
				distanceMoveColor[1]+=distanceColor[1][k]*tempWeightColor[1];
				distanceMoveColor[2]+=distanceColor[2][k]*tempWeightColor[2];

				weightSum+=tempWeight;		
				weightSumColor[0]+=tempWeightColor[0];
				weightSumColor[1]+=tempWeightColor[1];
				weightSumColor[2]+=tempWeightColor[2];
			}

			distanceMove=distanceMove/(weightSum+1);
			distanceMoveColor[0]=distanceMoveColor[0]/(weightSumColor[0]+1);
			distanceMoveColor[1]=distanceMoveColor[1]/(weightSumColor[1]+1);
			distanceMoveColor[2]=distanceMoveColor[2]/(weightSumColor[2]+1);
			m_distanceMove[j]=distanceMove;
			m_volumeChange[j]=distanceMove*m_meanAreas[j];
			for(int k=0;k<3;k++){
				m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];
				m_resultColor[j*3+k]=m_originalColors[j*3+k]+distanceMoveColor[k];
			}
			delete distanceProjection;
			delete distanceTangent;
			delete weightTangent;
			delete weightNormal;
			for(int k=0;k<3;k++){
				delete distanceColor[k];
				delete weightColor[k];
			}
			delete[] outFlag;

		}
	   RestoreVolume();
		for(int j=0;j<m_numOfPoints*3;j++){
			m_originalPointSet[j]=m_resultPointSet[j];
			m_originalColors[j]=m_resultColor[j];
		}
	}
	//	m_mapKnearest.clear();

	//delete[] m_originalPointSet;
	//m_originalPointSet=NULL;
	//delete[] m_originalColors;
	//m_originalColors=NULL;
	return;

}
void NonShiftBilateralFilter::ComputeFilterFour(int interractiveTime)
{
	m_volumeChange=new double[m_numOfPoints];
	m_resultPointSet=new double[m_numOfPoints*3];
	m_resultColor=new double[m_numOfPoints*3];
	m_distanceMove=new double[m_numOfPoints];
	ComputeNormal();
	if(m_IfFilterNormal==1){
		FilterNormal();
	}
	else{
		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
			m_filterNormals[j]=m_originalNormals[j];
		}
	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;

	for(int i=0;i<interractiveTime;i++){
	/*	ComputeNormal();
		if(m_IfFilterNormal==1){
			FilterNormal();
		}
		else{
			m_filterNormals=new double[m_numOfPoints*3];
			for(int j=0;j<m_numOfPoints*3;j++){
				m_filterNormals[j]=m_originalNormals[j];
			}
		}
		delete[] m_originalNormals;
		m_originalNormals=NULL;*/
	
		double* distanceProjection;
		double* distanceTangent;
		double* weightTangent;
		double* weightNormal;
		double* weightColor[3];
		double* distanceColor[3];
		double weightSum;
		double weightSumColor[3];
		double radiusMean;
		double distanceMove;
		double distanceMoveColor[3];
		for(int j=0;j<m_numOfPoints;j++){
			double localNormal[3];
			for(int k=0;k<3;k++){
				localNormal[k]=m_filterNormals[j*3+k];
			}
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(j);
			int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
			std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
			distanceProjection=new double[numOfNearest];
			distanceTangent=new double[numOfNearest];			
			weightTangent=new double[numOfNearest];
			weightNormal=new double[numOfNearest];
			int* outFlag=new int[numOfNearest];		
			for(int k=0;k<numOfNearest;k++){
				distanceProjection[k]=0.0;
				distanceTangent[k]=0.0;
				weightTangent[k]=0.0;
				weightNormal[k]=0.0;
				outFlag[k]=0;
			}
			for(int k=0;k<3;k++){
				distanceColor[k]=new double[numOfNearest];
				weightColor[k]=new double[numOfNearest];
				weightSumColor[k]=0;
				distanceMoveColor[k]=0;

			}
			weightSum=0;
			distanceMove=0;
			double sumDistanceProjection=0;
			double sumDistanceTangent=0;
			//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonDistance.txt";
			//FILE *fpout;
			//if((fpout = fopen(filename_pw, "a")) == NULL)
			//{
			//	int dkjkd;
			//	//MessageBox("can't open the file!");
			//}
			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if((m_filterNormals[(*vectorNearestIterator)*3]*localNormal[0]
				+m_filterNormals[(*vectorNearestIterator)*3+1]*localNormal[1]
				+m_filterNormals[(*vectorNearestIterator)*3+2]*localNormal[2])<m_thresholdDistanceNormal){
					outFlag[k]=0;
					numOfNearest-=1;
					continue;
				}
				else{
					outFlag[k]=1;
				}
				double localVectorPosition[3];				
				for(int w=0;w<3;w++){
					localVectorPosition[w]=m_originalPointSet[(*vectorNearestIterator)*3+w]-m_originalPointSet[j*3+w];
					distanceColor[w][k]=m_originalColors[(*vectorNearestIterator)*3+w]-m_originalColors[j*3+w];
				}
				distanceProjection[k]=Vector3Vector(localVectorPosition,localNormal);
				distanceTangent[k]=ComputeNormOfVector3(localVectorPosition);
				distanceTangent[k]=sqrt(distanceTangent[k]*distanceTangent[k]-distanceProjection[k]*distanceProjection[k]);
				sumDistanceTangent+=distanceTangent[k];
				sumDistanceProjection+=abs(distanceProjection[k]);


				//weightNormal[k]=exp(-0.5*distanceProjection*distanceProjection/m_variationNormal);
				//weightTangent=exp(-0.5*distanceTangent*distanceTangent/m_variationTangent);
				//double tempWeight=weightNormal*weightTangent;
				//distanceMove+=distanceProjection*tempWeight;
				//weightSum+=tempWeight;



				//		fprintf(fpout," %f ",distanceProjection[k]);




			}
			//fprintf(fpout,"\n");

			//	fclose(fpout);
			sumDistanceTangent/=(numOfNearest+1);
			sumDistanceProjection/=(numOfNearest+1);
			double localVariationNormal,localVarationTangent;
			localVariationNormal=m_variationNormal*m_variationNormal*sumDistanceProjection*sumDistanceProjection;
			if(localVariationNormal<0.000001)
				localVariationNormal=0.000001;
			localVarationTangent=m_variationTangent*m_variationTangent*sumDistanceTangent*sumDistanceTangent;
			if(localVarationTangent<0.000001)
				localVarationTangent=0.000001;
			vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();


			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if(outFlag[k]==0)
					continue;
				weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
				weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);
				for(int w=0;w<3;w++){
			//		if(abs(distanceColor[w][k])>m_thresholdFunction)
			//			weightColor[w][k]=0;
			//		else
						weightColor[w][k]=exp(-0.5*distanceColor[w][k]*distanceColor[w][k]/m_variationColor);

				}
				double tempWeight=weightNormal[k]*weightTangent[k];
				double tempWeightColor[3];
				tempWeightColor[0]=tempWeight;//*weightColor[0][k];
				tempWeightColor[1]=tempWeight;//*weightColor[1][k];
				tempWeightColor[2]=tempWeight;//*weightColor[2][k];
				distanceMove+=distanceProjection[k]*tempWeight;
				distanceMoveColor[0]+=distanceColor[0][k]*tempWeightColor[0];
				distanceMoveColor[1]+=distanceColor[1][k]*tempWeightColor[1];
				distanceMoveColor[2]+=distanceColor[2][k]*tempWeightColor[2];

				weightSum+=tempWeight;		
				weightSumColor[0]+=tempWeightColor[0];
				weightSumColor[1]+=tempWeightColor[1];
				weightSumColor[2]+=tempWeightColor[2];
			}

			distanceMove=distanceMove/(weightSum+1);
			distanceMoveColor[0]=distanceMoveColor[0]/(weightSumColor[0]+1);
			distanceMoveColor[1]=distanceMoveColor[1]/(weightSumColor[1]+1);
			distanceMoveColor[2]=distanceMoveColor[2]/(weightSumColor[2]+1);
			m_distanceMove[j]=distanceMove;
			m_volumeChange[j]=distanceMove*m_meanAreas[j];
			for(int k=0;k<3;k++){
				m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];
				m_resultColor[j*3+k]=m_originalColors[j*3+k]+distanceMoveColor[k];
			}
			delete distanceProjection;
			delete distanceTangent;
			delete weightTangent;
			delete weightNormal;
			for(int k=0;k<3;k++){
				delete distanceColor[k];
				delete weightColor[k];
			}
            delete[] outFlag;
		}
		RestoreVolume();
		for(int j=0;j<m_numOfPoints*3;j++){
			m_originalPointSet[j]=m_resultPointSet[j];
			m_originalColors[j]=m_resultColor[j];
		}
	}
	//	m_mapKnearest.clear();

	//delete[] m_originalPointSet;
	//m_originalPointSet=NULL;
	//delete[] m_originalColors;
	//m_originalColors=NULL;
	return;
}
void NonShiftBilateralFilter::RestoreVolume()
{
	//for(int i=0;i<m_numOfPoints;i++){
	//	std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
	//	std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	//	float localVolumeChange=0;
	//	for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
	//		localVolumeChange+=m_volumeChange[(*vectorIterator)];			
	//	}
	//	localVolumeChange+=m_volumeChange[i];
	//	localVolumeChange/=((*mapIterator).second.m_nearest.size()+1);
	//	m_resultPointSet[i*3+0]=m_originalPointSet[i*3+0]-((localVolumeChange-m_volumeChange[i])/m_meanAreas[i]-m_distanceMove[i])*m_filterNormals[i*3+0];
	//	m_resultPointSet[i*3+1]=m_originalPointSet[i*3+1]-((localVolumeChange-m_volumeChange[i])/m_meanAreas[i]-m_distanceMove[i])*m_filterNormals[i*3+1];
	//	m_resultPointSet[i*3+2]=m_originalPointSet[i*3+2]-((localVolumeChange-m_volumeChange[i])/m_meanAreas[i]-m_distanceMove[i])*m_filterNormals[i*3+2];
	//}

	//////////////////////////////////////////////////////////////////////////
	// 全局保护volume
	double sumVolumeChange=0;
	double sumArea=0;
	for(int i=0;i<m_numOfPoints;i++){
		sumVolumeChange+=m_volumeChange[i];
		sumArea+=m_meanAreas[i];		
	}
	double meanDistanceMove=sumVolumeChange/sumArea;
	
	for(int i=0;i<m_numOfPoints;i++){
        m_resultPointSet[i*3+0]+=meanDistanceMove*m_filterNormals[i*3+0];
		m_resultPointSet[i*3+1]+=meanDistanceMove*m_filterNormals[i*3+1];
		m_resultPointSet[i*3+2]+=meanDistanceMove*m_filterNormals[i*3+2];
	}
	return;
	

	
}
void NonShiftBilateralFilter::RestorPositinData()
{
	int colorDimension=m_numDimension-3;
	double maxminusmin=m_maxPosition-m_minPosition;
	if(m_resultPointSet!=NULL)
		delete m_resultPointSet;
	m_resultPointSet=new double[m_numOfPoints*3];
	if(m_resultColor!=NULL)
		delete m_resultColor;
	m_resultColor=new double[m_numOfPoints*(m_numDimension-3)];
	for(int i=0;i<m_numOfPoints;i++){
		m_resultPointSet[i*3+0]=m_resultPointSetAndColor[i*m_numDimension+0]*maxminusmin+m_minPosition;
		m_resultPointSet[i*3+1]=m_resultPointSetAndColor[i*m_numDimension+1]*maxminusmin+m_minPosition;
		m_resultPointSet[i*3+2]=m_resultPointSetAndColor[i*m_numDimension+2]*maxminusmin+m_minPosition;
		for(int j=0;j<colorDimension;j++){
			m_resultColor[i*colorDimension]=m_resultPointSetAndColor[i*m_numDimension+3+j];
        }
	}
	delete[] m_resultPointSetAndColor;
	m_resultPointSetAndColor=NULL;
}

void NonShiftBilateralFilter::CalculateMeanRadius()
{
	m_meanRadius=new double[m_numOfPoints];
	for(int i=0;i<m_numOfPoints;i++){
		
	}


}
void NonShiftBilateralFilter::ComputeFilterFive(int interractiveTime)
{

	m_volumeChange=new double[m_numOfPoints];
	m_resultPointSet=new double[m_numOfPoints*3];
	m_distanceMove=new double[m_numOfPoints];
	ComputeNormal();
	if(m_IfFilterNormal==1){
		FilterNormal();
	}
	else{
		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
			m_filterNormals[j]=m_originalNormals[j];
		}
	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;

	for(int i=0;i<interractiveTime;i++){
		/*ComputeNormal();
		if(m_IfFilterNormal==1){
			FilterNormal();
		}
		else{
			m_filterNormals=new double[m_numOfPoints*3];
			for(int j=0;j<m_numOfPoints*3;j++){
				m_filterNormals[j]=m_originalNormals[j];
			}
		}
		delete[] m_originalNormals;
		m_originalNormals=NULL;*/
		double* distanceProjection;
		double* distanceTangent;
		double* weightTangent;
		double* weightNormal;
		double weightSum;
		double radiusMean;
		double distanceMove;
		for(int j=0;j<m_numOfPoints;j++){
		
			double localNormal[3];
			for(int k=0;k<3;k++){
				localNormal[k]=m_filterNormals[j*3+k];
			}
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(j);
			int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
			std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
			distanceProjection=new double[numOfNearest];
			distanceTangent=new double[numOfNearest];
			weightTangent=new double[numOfNearest];
			weightNormal=new double[numOfNearest];
			int* outFlag=new int[numOfNearest];
			for(int k=0;k<numOfNearest;k++){
				distanceProjection[k]=0.0;
				distanceTangent[k]=0.0;
				weightTangent[k]=0.0;
				weightNormal[k]=0.0;
				outFlag[k]=0;
			}
			weightSum=0;
			distanceMove=0;
			double sumDistanceProjection=0;
			double sumDistanceTangent=0;
			//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonDistance.txt";
			//FILE *fpout;
			//if((fpout = fopen(filename_pw, "a")) == NULL)
			//{
			//	int dkjkd;
			//	//MessageBox("can't open the file!");
			//}
			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if((m_filterNormals[(*vectorNearestIterator)*3]*localNormal[0]
				+m_filterNormals[(*vectorNearestIterator)*3+1]*localNormal[1]
				+m_filterNormals[(*vectorNearestIterator)*3+2]*localNormal[2])<m_thresholdDistanceNormal){
					outFlag[k]=0;
					numOfNearest-=1;
					continue;
				}
				else{
					outFlag[k]=1;
				}
				double localVectorPosition[3];				
				for(int w=0;w<3;w++){
					localVectorPosition[w]=m_originalPointSet[(*vectorNearestIterator)*3+w]-m_originalPointSet[j*3+w];

				}
				distanceProjection[k]=Vector3Vector(localVectorPosition,localNormal);
				distanceTangent[k]=ComputeNormOfVector3(localVectorPosition);
				distanceTangent[k]=sqrt(distanceTangent[k]*distanceTangent[k]-distanceProjection[k]*distanceProjection[k]);
				sumDistanceTangent+=distanceTangent[k];
				sumDistanceProjection+=abs(distanceProjection[k]);
			
			}
				sumDistanceTangent/=(numOfNearest+1);
			sumDistanceProjection/=(numOfNearest+1);
			double localVariationNormal,localVarationTangent;
			localVariationNormal=m_variationNormal*m_variationNormal*sumDistanceProjection*sumDistanceProjection;
			if(localVariationNormal<0.000001)
				localVariationNormal=0.000001;
			localVarationTangent=m_variationTangent*m_variationTangent*sumDistanceTangent*sumDistanceTangent;
			if(localVarationTangent<0.000001)
				localVarationTangent=0.000001;
			vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();


			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if(outFlag[k]==0)
					continue;
				weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
				weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);
				double tempWeight=weightNormal[k]*weightTangent[k];
				distanceMove+=distanceProjection[k]*tempWeight;
				weightSum+=tempWeight;			
			}

			distanceMove=distanceMove/(weightSum+1);
			m_distanceMove[j]=distanceMove;
			m_volumeChange[j]=distanceMove*m_meanAreas[j];
			for(int k=0;k<3;k++){
				m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];

			}
			delete distanceProjection;
			delete distanceTangent;
			delete weightTangent;
			delete weightNormal;
			delete outFlag;

		}
		RestoreVolume();
		for(int j=0;j<m_numOfPoints*3;j++){
			m_originalPointSet[j]=m_resultPointSet[j];
		}
	}
	
	return;

}
void NonShiftBilateralFilter::ComputeFilterSix(int interractiveTime)
{
	//ComputeNormal();
	//if(m_IfFilterNormal==1){
	//	FilterNormal();
	//}
	//else{

	//	m_filterNormals=new double[m_numOfPoints*3];
	//	for(int i=0;i<m_numOfPoints*3;i++){
	//		m_filterNormals[i]=m_originalNormals[i];
	//	}
	//}
	//delete[] m_originalNormals;
	//m_originalNormals=NULL;
	m_volumeChange=new double[m_numOfPoints];
	m_resultPointSet=new double[m_numOfPoints*3];
	m_distanceMove=new double[m_numOfPoints];
	ComputeNormal();
	if(m_IfFilterNormal==1){
		FilterNormal();
	}
	else{

		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
			m_filterNormals[j]=m_originalNormals[j];
		}
	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;
	for(int i=0;i<interractiveTime;i++){

		/*ComputeNormal();
		if(m_IfFilterNormal==1){
			FilterNormal();
		}
		else{

			m_filterNormals=new double[m_numOfPoints*3];
			for(int j=0;j<m_numOfPoints*3;j++){
				m_filterNormals[j]=m_originalNormals[j];
			}
		}
		delete[] m_originalNormals;
		m_originalNormals=NULL;*/

		double* distanceProjection;
		double* distanceTangent;
		double* weightTangent;
		double* weightNormal;
		double weightSum;
		double radiusMean;
		double distanceMove;
		for(int j=0;j<m_numOfPoints;j++){
			double localNormal[3];
			for(int k=0;k<3;k++){
				localNormal[k]=m_filterNormals[j*3+k];
			}
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(j);
			int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
			std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
			distanceProjection=new double[numOfNearest];
			distanceTangent=new double[numOfNearest];
			weightTangent=new double[numOfNearest];
			weightNormal=new double[numOfNearest];
			int* outFlag=new int[numOfNearest];
			for(int k=0;k<numOfNearest;k++){
				distanceProjection[k]=0.0;
				distanceTangent[k]=0.0;
				weightTangent[k]=0.0;
				weightNormal[k]=0.0;
				outFlag[k]=0;
			}
			weightSum=0;
			distanceMove=0;
			double sumDistanceProjection=0;
			double sumDistanceTangent=0;
			//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonDistance.txt";
			//FILE *fpout;
			//if((fpout = fopen(filename_pw, "a")) == NULL)
			//{
			//	int dkjkd;
			//	//MessageBox("can't open the file!");
			//}
			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if((m_filterNormals[(*vectorNearestIterator)*3]*localNormal[0]
				+m_filterNormals[(*vectorNearestIterator)*3+1]*localNormal[1]
				+m_filterNormals[(*vectorNearestIterator)*3+2]*localNormal[2])<m_thresholdDistanceNormal){
					outFlag[k]=0;
					numOfNearest-=1;
					continue;
				}
				double localVectorPosition[3];				
				for(int w=0;w<3;w++){
					localVectorPosition[w]=m_originalPointSet[(*vectorNearestIterator)*3+w]-m_originalPointSet[j*3+w];

				}
				distanceProjection[k]=Vector3Vector(localVectorPosition,localNormal);
				distanceTangent[k]=ComputeNormOfVector3(localVectorPosition);
				distanceTangent[k]=sqrt(distanceTangent[k]*distanceTangent[k]-distanceProjection[k]*distanceProjection[k]);
				sumDistanceTangent+=distanceTangent[k];
				sumDistanceProjection+=abs(distanceProjection[k]);
				//weightNormal[k]=exp(-0.5*distanceProjection*distanceProjection/m_variationNormal);
				//weightTangent=exp(-0.5*distanceTangent*distanceTangent/m_variationTangent);
				//double tempWeight=weightNormal*weightTangent;
				//distanceMove+=distanceProjection*tempWeight;
				//weightSum+=tempWeight;



				//		fprintf(fpout," %f ",distanceProjection[k]);




			}
			//			fprintf(fpout,"\n");

			//	fclose(fpout);
			sumDistanceTangent/=(numOfNearest+1);
			sumDistanceProjection/=(numOfNearest+1);
			double localVariationNormal,localVarationTangent;
			localVariationNormal=m_variationNormal*m_variationNormal*sumDistanceProjection*sumDistanceProjection;
			if(localVariationNormal<0.000001)
				localVariationNormal=0.000001;
			localVarationTangent=m_variationTangent*m_variationTangent*sumDistanceTangent*sumDistanceTangent;
			if(localVarationTangent<0.000001)
				localVarationTangent=0.000001;
			vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();


			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if(outFlag[k]==0)
					continue;
				weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
				weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);
				double tempWeight=weightNormal[k]*weightTangent[k]*m_meanAreas[(*vectorNearestIterator)];
				tempWeight*=abs(m_filterNormals[(*vectorNearestIterator)*3]*m_filterNormals[j*3]
				+m_filterNormals[(*vectorNearestIterator)*3+1]*m_filterNormals[j*3+1]
				+m_filterNormals[(*vectorNearestIterator)*3+2]*m_filterNormals[j*3+2]);
				distanceMove+=distanceProjection[k]*tempWeight;
				weightSum+=tempWeight;			
			}

			distanceMove=distanceMove/(weightSum+m_meanAreas[j]);
			m_distanceMove[j]=distanceMove;
			m_volumeChange[j]=distanceMove*m_meanAreas[j];
			for(int k=0;k<3;k++){
				m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];

			}
			delete distanceProjection;
			delete distanceTangent;
			delete weightTangent;
			delete weightNormal;
			delete outFlag;

		}
		RestoreVolume();
		for(int j=0;j<m_numOfPoints*3;j++){
			m_originalPointSet[j]=m_resultPointSet[j];
		}
	}
	/*delete[] m_originalPointSet;
	m_originalPointSet=NULL;*/
	return;
}

void NonShiftBilateralFilter:: EstimateNormalDirection()
{

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

void NonShiftBilateralFilter::ComputeFilterSeven(int interractiveTime)
{
	m_volumeChange=new double[m_numOfPoints];
	m_resultPointSet=new double[m_numOfPoints*3];
	m_distanceMove=new double[m_numOfPoints];
	ComputeNormal();
	if(m_IfFilterNormal==1){
		FilterNormal();
	}
	else{
		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
			m_filterNormals[j]=m_originalNormals[j];
		}
	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;

	for(int i=0;i<interractiveTime;i++){
		/*ComputeNormal();
		if(m_IfFilterNormal==1){
			FilterNormal();
		}
		else{
			if(m_filterNormals!=NULL){
				delete[] m_filterNormals;
				m_filterNormals=NULL;
			}
			m_filterNormals=new double[m_numOfPoints*3];
			for(int j=0;j<m_numOfPoints*3;j++){
				m_filterNormals[j]=m_originalNormals[j];
			}
		}
		delete[] m_originalNormals;
		m_originalNormals=NULL;*/
		double* distanceProjection;
		double* distanceTangent;
		double* distanceTangentNormal;
		double* weightTangent;
		double* weightNormal;
		double* weightTangentNormal;
		double weightSum;
		double radiusMean;
		double distanceMove;
		for(int j=0;j<m_numOfPoints;j++){

			double localNormal[3];
			for(int k=0;k<3;k++){
				localNormal[k]=m_filterNormals[j*3+k];
			}
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(j);
			int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
			std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
			distanceProjection=new double[numOfNearest];
			distanceTangent=new double[numOfNearest];
			distanceTangentNormal=new double[numOfNearest];
			weightTangent=new double[numOfNearest];
			weightNormal=new double[numOfNearest];
			weightTangentNormal=new double[numOfNearest];
			int* outFlag=new int[numOfNearest];
			for(int k=0;k<numOfNearest;k++){
				distanceProjection[k]=0.0;
				distanceTangent[k]=0.0;
				distanceTangentNormal[k]=0;
				weightTangent[k]=0.0;
				weightNormal[k]=0.0;
				weightTangentNormal[k]=0;
				outFlag[k]=0;
			}
			weightSum=0;
			distanceMove=0;
			double sumDistanceProjection=0;
			double sumDistanceTangent=0;
			double sumDistanceTangentNormal=0;
			//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonDistance.txt";
			//FILE *fpout;
			//if((fpout = fopen(filename_pw, "a")) == NULL)
			//{
			//	int dkjkd;
			//	//MessageBox("can't open the file!");
			//}
			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){			
				if((m_filterNormals[i*3]*m_filterNormals[(*vectorNearestIterator)*3]
				+m_filterNormals[i*3+1]*m_filterNormals[(*vectorNearestIterator)*3+1]
				+m_filterNormals[i*3+2]*m_filterNormals[(*vectorNearestIterator)*3+2])<m_thresholdDistanceNormal){
					outFlag[k]=0;

					numOfNearest-=1;
					continue;
				}
				else{
					outFlag[k]=1;
				}
				double localVectorPosition[3];
				for(int w=0;w<3;w++){
					localVectorPosition[w]=m_originalPointSet[(*vectorNearestIterator)*3+w]-m_originalPointSet[j*3+w];					
				}
				double localNearestNormal[3];
				for(int w=0;w<3;w++){
					localNearestNormal[w]=m_filterNormals[(*vectorNearestIterator)*3+w];
				}
				distanceProjection[k]=Vector3Vector(localVectorPosition,localNormal);
				distanceTangent[k]=ComputeNormOfVector3(localVectorPosition);
				distanceTangentNormal[k]=Vector3Vector(localVectorPosition,localNearestNormal);
				//distanceTangent[k]=sqrt(distanceTangent[k]*distanceTangent[k]-distanceProjection[k]*distanceProjection[k]);
				sumDistanceTangent+=distanceTangent[k];
				sumDistanceProjection+=abs(distanceProjection[k]);
				sumDistanceTangentNormal+=abs(distanceTangentNormal[k]);
				//weightNormal[k]=exp(-0.5*distanceProjection*distanceProjection/m_variationNormal);
				//weightTangent=exp(-0.5*distanceTangent*distanceTangent/m_variationTangent);
				//double tempWeight=weightNormal*weightTangent;
				//distanceMove+=distanceProjection*tempWeight;
				//weightSum+=tempWeight;



				//		fprintf(fpout," %f ",distanceProjection[k]);	


			}
			//fprintf(fpout,"\n");

			//	fclose(fpout);
			sumDistanceTangent/=(numOfNearest+1);
			sumDistanceProjection/=(numOfNearest+1);
			sumDistanceTangentNormal/=(numOfNearest+1);
			double localVariationNormal,localVarationTangent;
			localVariationNormal=m_variationNormal*sumDistanceProjection*sumDistanceProjection;
			if(localVariationNormal<0.000001)
				localVariationNormal=0.000001;
			if(m_variationNormal*sumDistanceTangentNormal*sumDistanceTangentNormal<0.000001){
				localVariationNormal=0.000001;
			}
			localVarationTangent=m_variationTangent*sumDistanceTangent*sumDistanceTangent;
			if(localVarationTangent<0.000001)
				localVarationTangent=0.000001;
			vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();


			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if(outFlag[k]==0)
					continue;
				weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
				weightTangentNormal[k]=exp(-0.5*distanceTangentNormal[k]*distanceTangentNormal[k]/localVariationNormal);
				weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);				
				double tempWeight=weightNormal[k]*weightTangent[k]*weightTangentNormal[k];
				distanceMove+=distanceProjection[k]*tempWeight;
				weightSum+=tempWeight;			
			}

			distanceMove=distanceMove/(weightSum+1);
			m_distanceMove[j]=distanceMove;
			m_volumeChange[j]=distanceMove*m_meanAreas[j];
			for(int k=0;k<3;k++){
				m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];
			}
			delete distanceProjection;
			delete distanceTangent;
			delete[] distanceTangentNormal;
			delete[] weightTangentNormal;
			delete weightTangent;
			delete weightNormal;
			delete outFlag;


		}
		//	RestoreVolume();
		for(int j=0;j<m_numOfPoints*3;j++){
			m_originalPointSet[j]=m_resultPointSet[j];
		}
	}
	//	m_mapKnearest.clear();

	//delete[] m_originalPointSet;
	//m_originalPointSet=NULL;
	return;

}

void NonShiftBilateralFilter::ComputeFilterTen(int interractiveTime)//Pauly 方法， isotropic 滤波
{
	
//	m_volumeChange=new double[m_numOfPoints];
	m_resultPointSet=new double[m_numOfPoints*3];
	for (int i=0;i<interractiveTime;i++){
		for(int j=0;j<m_numOfPoints;j++){
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(j);
			int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
			std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
			double* vectorPoints;
			double* localWeight;
			double localLapVector[3];
			double sumOfWeight=0;
			vectorPoints=new double[numOfNearest*3];
			localWeight=new double[numOfNearest];
			for(int num=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,num++){
				double localDistance=0;
				for(int k=0;k<3;k++){
					vectorPoints[num*3+k]=m_originalPointSet[(*vectorNearestIterator)*3+k]-m_originalPointSet[j*3+k];
				}
				for(int k=0;k<3;k++){
					localDistance+=vectorPoints[num*3+k]*vectorPoints[num*3+k];
				}
				localWeight[num]=1/sqrt(localDistance);
				sumOfWeight+=localWeight[num];
			}
			for(int k=0;k<3;k++){
				localLapVector[k]=0;
			}
			for(int k=0;k<numOfNearest;k++){
				localLapVector[0]+=localWeight[k]*vectorPoints[k*3];
				localLapVector[1]+=localWeight[k]*vectorPoints[k*3+1];
				localLapVector[2]+=localWeight[k]*vectorPoints[k*3+2];
			}
			for(int k=0;k<3;k++){
				localLapVector[k]/=sumOfWeight;
				m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+localLapVector[k]*m_timeStep;
			}
			delete[] vectorPoints;
			delete[] localWeight;			
		}
		for(int j=0;j<m_numOfPoints*3;j++){
			m_originalPointSet[j]=m_resultPointSet[j];
		}
	}
}
void NonShiftBilateralFilter::ComputeFilterElevn(int interractiveTime)// Marx 方法， MLS
{

}
void NonShiftBilateralFilter::ComputeFilterTwentyone(int interractiveTime)//定义在流形上的双边滤波
{
	m_volumeChange=new double[m_numOfPoints];
	m_resultPointSet=new double[m_numOfPoints*3];
	m_distanceMove=new double[m_numOfPoints];
	/*
	*	对法向进行滤波
	*/
	ComputeNormal();
	if(m_IfFilterNormal==1){
		FilterNormal();
	}
	else{
		if(m_filterNormals!=NULL){
			delete m_filterNormals;
			m_filterNormals=NULL;
		}
		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
			m_filterNormals[j]=m_originalNormals[j];
		}
	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;
	for(int i=0;i<interractiveTime;i++){
		
		double* distanceProjection;
		double* distanceTangent;
		double* weightTangent;
		double* weightNormal;
		int* outFlag;
		double weightSum;
		double radiusMean;
		double distanceMove;

		/*
		*	对位置进行滤波
		*/
		for(int j=0;j<m_numOfPoints;j++){

			double localNormal[3];
			for(int k=0;k<3;k++){
				localNormal[k]=m_filterNormals[j*3+k];
			}
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(j);
			int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
			std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
			distanceProjection=new double[numOfNearest];
			distanceTangent=new double[numOfNearest];
			weightTangent=new double[numOfNearest];
			weightNormal=new double[numOfNearest];
			outFlag=new int[numOfNearest];
			for(int www=0;www<numOfNearest;www++){
				distanceProjection[www]=0.0;
				distanceTangent[www]=0.0;
				weightTangent[www]=0.0;
				weightNormal[www]=0.0;
				outFlag[www]=0;
			}
			weightSum=0;
			distanceMove=0;
			double sumDistanceProjection=0;
			double sumDistanceTangent=0;
		
			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){			
				/*
				*	法向的threshold处理
				*/
				if((m_filterNormals[i*3]*m_filterNormals[(*vectorNearestIterator)*3]
				+m_filterNormals[i*3+1]*m_filterNormals[(*vectorNearestIterator)*3+1]
				+m_filterNormals[i*3+2]*m_filterNormals[(*vectorNearestIterator)*3+2])<m_thresholdDistanceNormal){
					outFlag[k]=0;

					numOfNearest-=1;
					continue;
				}
				else{
				
					outFlag[k]=1;
				}
				double localVectorPosition[3];
				for(int w=0;w<3;w++){
					localVectorPosition[w]=m_originalPointSet[(*vectorNearestIterator)*3+w]-m_originalPointSet[j*3+w];					
				}
				distanceProjection[k]=Vector3Vector(localVectorPosition,localNormal);
				distanceTangent[k]=ComputeNormOfVector3(localVectorPosition);
				//distanceTangent[k]=sqrt(distanceTangent[k]*distanceTangent[k]-distanceProjection[k]*distanceProjection[k]);
				sumDistanceTangent+=distanceTangent[k];
				sumDistanceProjection+=abs(distanceProjection[k]);
				//weightNormal[k]=exp(-0.5*distanceProjection*distanceProjection/m_variationNormal);
				//weightTangent=exp(-0.5*distanceTangent*distanceTangent/m_variationTangent);
				//double tempWeight=weightNormal*weightTangent;
				//distanceMove+=distanceProjection*tempWeight;
				//weightSum+=tempWeight;



				//		fprintf(fpout," %f ",distanceProjection[k]);	


			}
			//fprintf(fpout,"\n");

			//	fclose(fpout);
			sumDistanceTangent/=(numOfNearest+1);
			sumDistanceProjection/=(numOfNearest+1);
			double localVariationNormal,localVarationTangent;
			localVariationNormal=m_variationNormal*sumDistanceProjection*sumDistanceProjection;
			if(localVariationNormal<0.000001)
				localVariationNormal=0.000001;
			localVarationTangent=m_variationTangent*sumDistanceTangent*sumDistanceTangent;
			if(localVarationTangent<0.000001)
				localVarationTangent=0.000001;
			vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();


			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if(outFlag[k]==0)
					continue;
				weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
				weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);				
				double tempWeight=weightNormal[k]*weightTangent[k];
				distanceMove+=distanceProjection[k]*tempWeight;
				weightSum+=tempWeight;			
			}

			distanceMove=distanceMove/(weightSum+1);
			m_distanceMove[j]=distanceMove;
			m_volumeChange[j]=distanceMove*m_meanAreas[j];
			for(int k=0;k<3;k++){
				m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];
			}
			delete distanceProjection;
			delete distanceTangent;
			delete weightTangent;
			delete weightNormal;
			delete outFlag;



		}
		//	RestoreVolume();
		for(int j=0;j<m_numOfPoints*3;j++){
			m_originalPointSet[j]=m_resultPointSet[j];
		}

	}
	//	m_mapKnearest.clear();

	//delete[] m_originalPointSet;
	//m_originalPointSet=NULL;
	return;

}
void NonShiftBilateralFilter::ComputeFilterTwentytwo(int interractiveTime)//定义在切空间上的双边滤波
{
	m_volumeChange=new double[m_numOfPoints];
	m_resultPointSet=new double[m_numOfPoints*3];
	m_distanceMove=new double[m_numOfPoints];
	ComputeNormal();
	if(m_IfFilterNormal==1){
		FilterNormal();
	}
	else{
		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
			m_filterNormals[j]=m_originalNormals[j];
		}
	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;

	for(int i=0;i<interractiveTime;i++){
		/*
		*	对法向进行滤波
		*/
	/*	ComputeNormal();
		if(m_IfFilterNormal==1){
			FilterNormal();
		}
		else{
			if(m_filterNormals!=NULL){
				delete m_filterNormals;
				m_filterNormals=NULL;
			}
			m_filterNormals=new double[m_numOfPoints*3];
			for(int j=0;j<m_numOfPoints*3;j++){
				m_filterNormals[j]=m_originalNormals[j];
			}
		}
		delete[] m_originalNormals;
		m_originalNormals=NULL;*/
		double* distanceProjection;
		double* distanceTangent;
		double* weightTangent;
		double* weightNormal;
		double weightSum;
		double radiusMean;
		double distanceMove;

		/*
		*	对位置进行滤波
		*/
		for(int j=0;j<m_numOfPoints;j++){

			double localNormal[3];
			for(int k=0;k<3;k++){
				localNormal[k]=m_filterNormals[j*3+k];
			}
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(j);
			int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
			std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
			distanceProjection=new double[numOfNearest];
			distanceTangent=new double[numOfNearest];
			weightTangent=new double[numOfNearest];
			weightNormal=new double[numOfNearest];
			int* outFlag=new int[numOfNearest];
			for(int k=0;k<numOfNearest;k++){
				distanceProjection[k]=0.0;
				distanceTangent[k]=0.0;
				weightTangent[k]=0.0;
				weightNormal[k]=0.0;
				outFlag[k]=0;
			}
			weightSum=0;
			distanceMove=0;
			double sumDistanceProjection=0;
			double sumDistanceTangent=0;
		
			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){			
				/*
				*	法向的threshold处理
				*/
				if((m_filterNormals[j*3]*m_filterNormals[(*vectorNearestIterator)*3]
				+m_filterNormals[j*3+1]*m_filterNormals[(*vectorNearestIterator)*3+1]
				+m_filterNormals[j*3+2]*m_filterNormals[(*vectorNearestIterator)*3+2])<m_thresholdDistanceNormal){
					outFlag[k]=0;

					numOfNearest-=1;
					continue;
				}
				else{
					outFlag[k]=1;
				}
				double localVectorPosition[3];
				for(int w=0;w<3;w++){
					localVectorPosition[w]=m_originalPointSet[(*vectorNearestIterator)*3+w]-m_originalPointSet[j*3+w];					
				}
				distanceProjection[k]=Vector3Vector(localVectorPosition,localNormal);
				distanceTangent[k]=ComputeNormOfVector3(localVectorPosition);
				distanceTangent[k]=sqrt(distanceTangent[k]*distanceTangent[k]-distanceProjection[k]*distanceProjection[k]);
				sumDistanceTangent+=distanceTangent[k];
				sumDistanceProjection+=abs(distanceProjection[k]);
				//weightNormal[k]=exp(-0.5*distanceProjection*distanceProjection/m_variationNormal);
				//weightTangent=exp(-0.5*distanceTangent*distanceTangent/m_variationTangent);
				//double tempWeight=weightNormal*weightTangent;
				//distanceMove+=distanceProjection*tempWeight;
				//weightSum+=tempWeight;



				//		fprintf(fpout," %f ",distanceProjection[k]);	


			}
			//fprintf(fpout,"\n");

			//	fclose(fpout);
			sumDistanceTangent/=(numOfNearest+1);
			sumDistanceProjection/=(numOfNearest+1);
			double localVariationNormal,localVarationTangent;
			localVariationNormal=m_variationNormal*sumDistanceProjection*sumDistanceProjection;
			if(localVariationNormal<0.000001)
				localVariationNormal=0.000001;
			localVarationTangent=m_variationTangent*sumDistanceTangent*sumDistanceTangent;
			if(localVarationTangent<0.000001)
				localVarationTangent=0.000001;
			vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();


			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if(outFlag[k]==0)
					continue;
				weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
				weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);				
				double tempWeight=weightNormal[k]*weightTangent[k];
				distanceMove+=distanceProjection[k]*tempWeight;
				weightSum+=tempWeight;			
			}

			distanceMove=distanceMove/(weightSum+1);
			m_distanceMove[j]=distanceMove;
			m_volumeChange[j]=distanceMove*m_meanAreas[j];
			for(int k=0;k<3;k++){
				m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];
			}
			delete distanceProjection;
			delete distanceTangent;
			delete weightTangent;
			delete weightNormal;
			delete outFlag;



		}
		//	RestoreVolume();
		for(int j=0;j<m_numOfPoints*3;j++){
			m_originalPointSet[j]=m_resultPointSet[j];
		}

	}
	//	m_mapKnearest.clear();

	//delete[] m_originalPointSet;
	//m_originalPointSet=NULL;
	return;

}
void NonShiftBilateralFilter::ComputeFilterTwentyThree(int interractiveTime)//加上了法向距离项
{
	m_volumeChange=new double[m_numOfPoints];
	m_resultPointSet=new double[m_numOfPoints*3];
	m_distanceMove=new double[m_numOfPoints];

	/*
	*	对法向进行滤波
	*/
	ComputeNormal();
	if(m_IfFilterNormal==1){
		FilterNormal();
	}
	else{
		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
			m_filterNormals[j]=m_originalNormals[j];
		}
	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;

	for(int i=0;i<interractiveTime;i++){
		/*
		*	对法向进行滤波
		*/
		/*ComputeNormal();
		if(m_IfFilterNormal==1){
			FilterNormal();
		}
		else{
			if(m_filterNormals!=NULL){
				delete m_filterNormals;
				m_filterNormals=NULL;
			}
			m_filterNormals=new double[m_numOfPoints*3];
			for(int j=0;j<m_numOfPoints*3;j++){
				m_filterNormals[j]=m_originalNormals[j];
			}
		}
		delete[] m_originalNormals;
		m_originalNormals=NULL;*/
		double* distanceProjection;
		double* distanceTangent;
		double* weightTangent;
		double* weightNormal;
		double* normalDot;
		int* outFlag;
		double weightSum;
		double radiusMean;
		double distanceMove;

		/*
		*	对位置进行滤波
		*/
		for(int j=0;j<m_numOfPoints;j++){

			double localNormal[3];
			for(int k=0;k<3;k++){
				localNormal[k]=m_filterNormals[j*3+k];
			}
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(j);
			int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
			std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
			distanceProjection=new double[numOfNearest];
			distanceTangent=new double[numOfNearest];
			weightTangent=new double[numOfNearest];
			weightNormal=new double[numOfNearest];
			normalDot=new double[numOfNearest];
			for(int www=0;www<numOfNearest;www++){
				distanceProjection[www]=0;
				distanceTangent[www]=0;
				weightNormal[www]=0;
				weightTangent[www]=0;
				normalDot[www]=0;
			}
			
	        outFlag=new int[numOfNearest];
			weightSum=0;
			distanceMove=0;
			double sumDistanceProjection=0;
			double sumDistanceTangent=0;
			//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonDistance.txt";
			//FILE *fpout;
			//if((fpout = fopen(filename_pw, "a")) == NULL)
			//{
			//	int dkjkd;
			//	//MessageBox("can't open the file!");
			//}
			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){			
				/*
				*	法向的threshold处理
				*/
				if((m_filterNormals[j*3]*m_filterNormals[(*vectorNearestIterator)*3]
				+m_filterNormals[j*3+1]*m_filterNormals[(*vectorNearestIterator)*3+1]
				+m_filterNormals[j*3+2]*m_filterNormals[(*vectorNearestIterator)*3+2])<m_thresholdDistanceNormal){
					outFlag[k]=0;

					numOfNearest-=1;
					continue;
				}
				else{
					outFlag[k]=1;
				}
				double localVectorPosition[3];
				for(int w=0;w<3;w++){
					localVectorPosition[w]=m_originalPointSet[(*vectorNearestIterator)*3+w]-m_originalPointSet[j*3+w];					
				}

				double localNearestNormal[3];
				for(int w=0;w<3;w++){
					localNearestNormal[w]=m_filterNormals[(*vectorNearestIterator)*3+w];
				}
				distanceProjection[k]=Vector3Vector(localVectorPosition,localNormal);
				normalDot[k]=abs(Vector3Vector(localNormal,localNearestNormal));
				distanceTangent[k]=ComputeNormOfVector3(localVectorPosition);
				distanceTangent[k]=sqrt(distanceTangent[k]*distanceTangent[k]-distanceProjection[k]*distanceProjection[k]);
				sumDistanceTangent+=distanceTangent[k];
				sumDistanceProjection+=abs(distanceProjection[k]);
				//weightNormal[k]=exp(-0.5*distanceProjection*distanceProjection/m_variationNormal);
				//weightTangent=exp(-0.5*distanceTangent*distanceTangent/m_variationTangent);
				//double tempWeight=weightNormal*weightTangent;
				//distanceMove+=distanceProjection*tempWeight;
				//weightSum+=tempWeight;



				//		fprintf(fpout," %f ",distanceProjection[k]);	


			}
			//fprintf(fpout,"\n");

			//	fclose(fpout);
			sumDistanceTangent/=(numOfNearest+1);
			sumDistanceProjection/=(numOfNearest+1);
			double localVariationNormal,localVarationTangent;
			localVariationNormal=m_variationNormal*sumDistanceProjection*sumDistanceProjection;
			if(localVariationNormal<0.000001)
				localVariationNormal=0.000001;
			localVarationTangent=m_variationTangent*sumDistanceTangent*sumDistanceTangent;
			if(localVarationTangent<0.000001)
				localVarationTangent=0.000001;
			vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();


			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if(outFlag[k]==0)
					continue;
				weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
				weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);
	
				double tempWeight=weightNormal[k]*weightTangent[k]*normalDot[k];

				distanceMove+=distanceProjection[k]*tempWeight;
				weightSum+=tempWeight;			
			}

			distanceMove=distanceMove/(weightSum+1);
			m_distanceMove[j]=distanceMove;
			m_volumeChange[j]=distanceMove*m_meanAreas[j];
			for(int k=0;k<3;k++){
				m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];
			}
			delete distanceProjection;
			delete distanceTangent;
			delete weightTangent;
			delete weightNormal;
			delete outFlag;
			delete normalDot;



		}
		//	RestoreVolume();
		for(int j=0;j<m_numOfPoints*3;j++){
			m_originalPointSet[j]=m_resultPointSet[j];
		}

	}
	//	m_mapKnearest.clear();

	//delete[] m_originalPointSet;
	//m_originalPointSet=NULL;
	return;

}
void NonShiftBilateralFilter::ComputeFilterTwentyFour(int interractiveTime)//定义在切空间上的双边滤波
{
	m_volumeChange=new double[m_numOfPoints];
	m_resultPointSet=new double[m_numOfPoints*3];
	m_distanceMove=new double[m_numOfPoints];
	ComputeNormal();
	if(m_IfFilterNormal==1){
		FilterNormal();
	}
	else{
		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
			m_filterNormals[j]=m_originalNormals[j];
		}
	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;

	for(int i=0;i<interractiveTime;i++){
		/*
		*	对法向进行滤波
		*/
		/*	ComputeNormal();
		if(m_IfFilterNormal==1){
		FilterNormal();
		}
		else{
		if(m_filterNormals!=NULL){
		delete m_filterNormals;
		m_filterNormals=NULL;
		}
		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
		m_filterNormals[j]=m_originalNormals[j];
		}
		}
		delete[] m_originalNormals;
		m_originalNormals=NULL;*/
		double* distanceProjection;
		double* distanceTangent;
		double* weightTangent;
		double* weightNormal;
		double weightSum;
		double radiusMean;
		double distanceMove;

		/*
		*	对位置进行滤波
		*/
		for(int j=0;j<m_numOfPoints;j++){

			double localNormal[3];
			for(int k=0;k<3;k++){
				localNormal[k]=m_filterNormals[j*3+k];
			}
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(j);
			int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
			std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
			distanceProjection=new double[numOfNearest];
			distanceTangent=new double[numOfNearest];
			weightTangent=new double[numOfNearest];
			weightNormal=new double[numOfNearest];
			int* outFlag=new int[numOfNearest];
			for(int k=0;k<numOfNearest;k++){
				distanceProjection[k]=0.0;
				distanceTangent[k]=0.0;
				weightTangent[k]=0.0;
				weightNormal[k]=0.0;
				outFlag[k]=0;
			}
			weightSum=0;
			distanceMove=0;
			double sumDistanceProjection=0;
			double sumDistanceTangent=0;
			//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonDistance.txt";
			//FILE *fpout;
			//if((fpout = fopen(filename_pw, "a")) == NULL)
			//{
			//	int dkjkd;
			//	//MessageBox("can't open the file!");
			//}
			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){			
				/*
				*	法向的threshold处理
				*/
				if((m_filterNormals[j*3]*m_filterNormals[(*vectorNearestIterator)*3]
				+m_filterNormals[j*3+1]*m_filterNormals[(*vectorNearestIterator)*3+1]
				+m_filterNormals[j*3+2]*m_filterNormals[(*vectorNearestIterator)*3+2])<m_thresholdDistanceNormal){
					outFlag[k]=0;

					numOfNearest-=1;
					continue;
				}
				else{
					outFlag[k]=1;
				}
				double localVectorPosition[3];
				for(int w=0;w<3;w++){
					localVectorPosition[w]=m_originalPointSet[(*vectorNearestIterator)*3+w]-m_originalPointSet[j*3+w];					
				}
				distanceProjection[k]=Vector3Vector(localVectorPosition,localNormal);
				distanceTangent[k]=ComputeNormOfVector3(localVectorPosition);
				distanceTangent[k]=sqrt(distanceTangent[k]*distanceTangent[k]-distanceProjection[k]*distanceProjection[k]);
				sumDistanceTangent+=distanceTangent[k];
				sumDistanceProjection+=abs(distanceProjection[k]);
				//weightNormal[k]=exp(-0.5*distanceProjection*distanceProjection/m_variationNormal);
				//weightTangent=exp(-0.5*distanceTangent*distanceTangent/m_variationTangent);
				//double tempWeight=weightNormal*weightTangent;
				//distanceMove+=distanceProjection*tempWeight;
				//weightSum+=tempWeight;



				//		fprintf(fpout," %f ",distanceProjection[k]);	


			}
			//fprintf(fpout,"\n");

			//	fclose(fpout);
			sumDistanceTangent/=(numOfNearest+1);
			sumDistanceProjection/=(numOfNearest+1);
			double localVariationNormal,localVarationTangent;
			localVariationNormal=m_variationNormal*sumDistanceProjection*sumDistanceProjection;
			if(localVariationNormal<0.000001)
				localVariationNormal=0.000001;
			localVarationTangent=m_variationTangent*sumDistanceTangent*sumDistanceTangent;
			if(localVarationTangent<0.000001)
				localVarationTangent=0.000001;
			vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();


			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if(outFlag[k]==0)
					continue;
				weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
				weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);				
				double tempWeight=weightNormal[k]*weightTangent[k];
				distanceMove+=distanceProjection[k]*tempWeight;
				weightSum+=tempWeight;			
			}

			distanceMove=distanceMove/(weightSum+1);
			m_distanceMove[j]=distanceMove;
			m_volumeChange[j]=distanceMove*m_meanAreas[j];
			for(int k=0;k<3;k++){
				m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];
			}
			delete distanceProjection;
			delete distanceTangent;
			delete weightTangent;
			delete weightNormal;
			delete outFlag;



		}
		RestoreVolume();
		for(int j=0;j<m_numOfPoints*3;j++){
			m_originalPointSet[j]=m_resultPointSet[j];
		}

	}
	//	m_mapKnearest.clear();

	//delete[] m_originalPointSet;
	//m_originalPointSet=NULL;
	return;

}
void NonShiftBilateralFilter::ComputeFilterTwentyFive(int interractiveTime)//加上了法向投影项的面积不一致性质
{
	m_volumeChange=new double[m_numOfPoints];
	m_resultPointSet=new double[m_numOfPoints*3];
	m_distanceMove=new double[m_numOfPoints];

	/*
	*	对法向进行滤波
	*/
	ComputeNormal();
	if(m_IfFilterNormal==1){
		FilterNormal();
	}
	else{
		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
			m_filterNormals[j]=m_originalNormals[j];
		}
	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;

	for(int i=0;i<interractiveTime;i++){
		/*
		*	对法向进行滤波
		*/
		/*ComputeNormal();
		if(m_IfFilterNormal==1){
		FilterNormal();
		}
		else{
		if(m_filterNormals!=NULL){
		delete m_filterNormals;
		m_filterNormals=NULL;
		}
		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
		m_filterNormals[j]=m_originalNormals[j];
		}
		}
		delete[] m_originalNormals;
		m_originalNormals=NULL;*/
		double* distanceProjection;
		double* distanceTangent;
		double* weightTangent;
		double* weightNormal;
		double* normalDot;
		int* outFlag;
		double weightSum;
		double radiusMean;
		double distanceMove;

		/*
		*	对位置进行滤波
		*/
		for(int j=0;j<m_numOfPoints;j++){

			double localNormal[3];
			for(int k=0;k<3;k++){
				localNormal[k]=m_filterNormals[j*3+k];
			}
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(j);
			int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
			std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
			distanceProjection=new double[numOfNearest];
			distanceTangent=new double[numOfNearest];
			weightTangent=new double[numOfNearest];
			weightNormal=new double[numOfNearest];
			normalDot=new double[numOfNearest];
			outFlag=new int[numOfNearest];
			for(int www=0;www<numOfNearest;www++){
				distanceProjection[www]=0;
				distanceTangent[www]=0;
				weightNormal[www]=0;
				weightTangent[www]=0;
				normalDot[www]=0;
				outFlag[www]=0;
			}

			
			weightSum=0;
			distanceMove=0;
			double sumDistanceProjection=0;
			double sumDistanceTangent=0;
			//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonDistance.txt";
			//FILE *fpout;
			//if((fpout = fopen(filename_pw, "a")) == NULL)
			//{
			//	int dkjkd;
			//	//MessageBox("can't open the file!");
			//}
			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){			
				/*
				*	法向的threshold处理
				*/
				if((m_filterNormals[j*3]*m_filterNormals[(*vectorNearestIterator)*3]
				+m_filterNormals[j*3+1]*m_filterNormals[(*vectorNearestIterator)*3+1]
				+m_filterNormals[j*3+2]*m_filterNormals[(*vectorNearestIterator)*3+2])<m_thresholdDistanceNormal){
					outFlag[k]=0;

					numOfNearest-=1;
					continue;
				}
				else{
					outFlag[k]=1;
				}
				double localVectorPosition[3];
				for(int w=0;w<3;w++){
					localVectorPosition[w]=m_originalPointSet[(*vectorNearestIterator)*3+w]-m_originalPointSet[j*3+w];					
				}

				double localNearestNormal[3];
				for(int w=0;w<3;w++){
					localNearestNormal[w]=m_filterNormals[(*vectorNearestIterator)*3+w];
				}
				distanceProjection[k]=Vector3Vector(localVectorPosition,localNormal);
				normalDot[k]=abs(Vector3Vector(localNormal,localNearestNormal));
				distanceTangent[k]=ComputeNormOfVector3(localVectorPosition);
				distanceTangent[k]=sqrt(distanceTangent[k]*distanceTangent[k]-distanceProjection[k]*distanceProjection[k]);
				sumDistanceTangent+=distanceTangent[k];
				sumDistanceProjection+=abs(distanceProjection[k]);
				//weightNormal[k]=exp(-0.5*distanceProjection*distanceProjection/m_variationNormal);
				//weightTangent=exp(-0.5*distanceTangent*distanceTangent/m_variationTangent);
				//double tempWeight=weightNormal*weightTangent;
				//distanceMove+=distanceProjection*tempWeight;
				//weightSum+=tempWeight;



				//		fprintf(fpout," %f ",distanceProjection[k]);	


			}
			//fprintf(fpout,"\n");

			//	fclose(fpout);
			sumDistanceTangent/=(numOfNearest+1);
			sumDistanceProjection/=(numOfNearest+1);
			double localVariationNormal,localVarationTangent;
			localVariationNormal=m_variationNormal*sumDistanceProjection*sumDistanceProjection;
			if(localVariationNormal<0.000001)
				localVariationNormal=0.000001;
			localVarationTangent=m_variationTangent*sumDistanceTangent*sumDistanceTangent;
			if(localVarationTangent<0.000001)
				localVarationTangent=0.000001;
			vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();


			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if(outFlag[k]==0)
					continue;
				weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
				weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);

				double tempWeight=weightNormal[k]*weightTangent[k]*normalDot[k];

				distanceMove+=distanceProjection[k]*tempWeight;
				weightSum+=tempWeight;			
			}

			distanceMove=distanceMove/(weightSum+1);
			m_distanceMove[j]=distanceMove;
			m_volumeChange[j]=abs(distanceMove*m_meanAreas[j]);
			for(int k=0;k<3;k++){
				m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];
			}
			delete distanceProjection;
			delete distanceTangent;
			delete weightTangent;
			delete weightNormal;
			delete outFlag;
			delete normalDot;



		}
			RestoreVolume();
		for(int j=0;j<m_numOfPoints*3;j++){
			m_originalPointSet[j]=m_resultPointSet[j];
		}

	}
	//	m_mapKnearest.clear();

	//delete[] m_originalPointSet;
	//m_originalPointSet=NULL;
	return;

}
void NonShiftBilateralFilter::ComputeFilterTwentySix(int interractiveTime)//加上了法向投影项的面积不一致性质
{
	m_volumeChange=new double[m_numOfPoints];
	m_resultPointSet=new double[m_numOfPoints*3];
	m_distanceMove=new double[m_numOfPoints];

	/*
	*	对法向进行滤波
	*/
	ComputeNormal();
	if(m_IfFilterNormal==1){
		FilterNormal();
	}
	else{
		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
			m_filterNormals[j]=m_originalNormals[j];
		}
	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;

	for(int i=0;i<interractiveTime;i++){
		
		double* distanceProjection;
		double* distanceTangent;
		double* weightTangent;
		double* weightNormal;
		double* normalDot;
		int* outFlag;
		double weightSum;
		double radiusMean;
		double distanceMove;

		/*
		*	对位置进行滤波
		*/
		for(int j=0;j<m_numOfPoints;j++){

			double localNormal[3];
			for(int k=0;k<3;k++){
				localNormal[k]=m_filterNormals[j*3+k];
			}
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(j);
			int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
			std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
			distanceProjection=new double[numOfNearest];
			distanceTangent=new double[numOfNearest];
			weightTangent=new double[numOfNearest];
			weightNormal=new double[numOfNearest];
			normalDot=new double[numOfNearest];
			outFlag=new int[numOfNearest];
			for(int www=0;www<numOfNearest;www++){
				distanceProjection[www]=0;
				distanceTangent[www]=0;
				weightNormal[www]=0;
				weightTangent[www]=0;
				normalDot[www]=0;
				outFlag[www]=0;
			}

			
			weightSum=0;
			distanceMove=0;
			double sumDistanceProjection=0;
			double sumDistanceTangent=0;
			//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonDistance.txt";
			//FILE *fpout;
			//if((fpout = fopen(filename_pw, "a")) == NULL)
			//{
			//	int dkjkd;
			//	//MessageBox("can't open the file!");
			//}
			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){			
				/*
				*	法向的threshold处理
				*/
				if((m_filterNormals[j*3]*m_filterNormals[(*vectorNearestIterator)*3]
				+m_filterNormals[j*3+1]*m_filterNormals[(*vectorNearestIterator)*3+1]
				+m_filterNormals[j*3+2]*m_filterNormals[(*vectorNearestIterator)*3+2])<m_thresholdDistanceNormal){
					outFlag[k]=0;

					numOfNearest-=1;
					continue;
				}
				else{
					outFlag[k]=1;
				}
				double localVectorPosition[3];
				for(int w=0;w<3;w++){
					localVectorPosition[w]=m_originalPointSet[(*vectorNearestIterator)*3+w]-m_originalPointSet[j*3+w];					
				}

				double localNearestNormal[3];
				for(int w=0;w<3;w++){
					localNearestNormal[w]=m_filterNormals[(*vectorNearestIterator)*3+w];
				}
				distanceProjection[k]=Vector3Vector(localVectorPosition,localNormal);
				normalDot[k]=abs(Vector3Vector(localNormal,localNearestNormal));
				distanceTangent[k]=ComputeNormOfVector3(localVectorPosition);
				distanceTangent[k]=sqrt(distanceTangent[k]*distanceTangent[k]-distanceProjection[k]*distanceProjection[k]);
				sumDistanceTangent+=distanceTangent[k];
				sumDistanceProjection+=abs(distanceProjection[k]);
				//weightNormal[k]=exp(-0.5*distanceProjection*distanceProjection/m_variationNormal);
				//weightTangent=exp(-0.5*distanceTangent*distanceTangent/m_variationTangent);
				//double tempWeight=weightNormal*weightTangent;
				//distanceMove+=distanceProjection*tempWeight;
				//weightSum+=tempWeight;



				//		fprintf(fpout," %f ",distanceProjection[k]);	


			}
			//fprintf(fpout,"\n");

			//	fclose(fpout);
			sumDistanceTangent/=(numOfNearest+1);
			sumDistanceProjection/=(numOfNearest+1);
			double localVariationNormal,localVarationTangent;
			localVariationNormal=m_variationNormal*sumDistanceProjection*sumDistanceProjection;
			if(localVariationNormal<0.000001)
				localVariationNormal=0.000001;
			localVarationTangent=m_variationTangent*sumDistanceTangent*sumDistanceTangent;
			if(localVarationTangent<0.000001)
				localVarationTangent=0.000001;
			vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();


			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
				if(outFlag[k]==0)
					continue;
				weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
				weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);

				double tempWeight=weightNormal[k]*weightTangent[k]*(m_meanAreas[(*vectorNearestIterator)]/m_kNearest)*normalDot[k];

				distanceMove+=distanceProjection[k]*tempWeight;
				weightSum+=tempWeight;			
			}

			distanceMove=distanceMove/(weightSum+m_meanAreas[j]/m_kNearest);
			m_distanceMove[j]=distanceMove;
			m_volumeChange[j]=distanceMove*m_meanAreas[j];
			for(int k=0;k<3;k++){
				m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];
			}
			delete distanceProjection;
			delete distanceTangent;
			delete weightTangent;
			delete weightNormal;
			delete outFlag;
			delete normalDot;



		}
		//	RestoreVolume();
		for(int j=0;j<m_numOfPoints*3;j++){
			m_originalPointSet[j]=m_resultPointSet[j];
		}

	}
	//	m_mapKnearest.clear();

	//delete[] m_originalPointSet;
	//m_originalPointSet=NULL;
	return;

}
void NonShiftBilateralFilter::ComputeFilterThirty(int interractiveTime)
{
	m_volumeChange=new double[m_numOfPoints];
	m_resultPointSet=new double[m_numOfPoints*3];
	m_distanceMove=new double[m_numOfPoints];

	/*
	*	对法向进行滤波
	*/
	ComputeNormal();
	if(m_IfFilterNormal==1){
		FilterNormal();
	}
	else{
		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
			m_filterNormals[j]=m_originalNormals[j];
		}
	}
	delete[] m_originalNormals;
	m_originalNormals=NULL;

	for(int i=0;i<interractiveTime;i++){
		/*
		*	对法向进行滤波
		*/
		/*ComputeNormal();
		if(m_IfFilterNormal==1){
		FilterNormal();
		}
		else{
		if(m_filterNormals!=NULL){
		delete m_filterNormals;
		m_filterNormals=NULL;
		}
		m_filterNormals=new double[m_numOfPoints*3];
		for(int j=0;j<m_numOfPoints*3;j++){
		m_filterNormals[j]=m_originalNormals[j];
		}
		}
		delete[] m_originalNormals;
		m_originalNormals=NULL;*/
		double* distanceProjection;
		double* distanceTangent;
		double* weightTangent;
		double* weightNormal;
		double* normalDot;
		int* outFlag;
		double weightSum;
		double radiusMean;
		double distanceMove;

		/*
		*	对位置进行滤波
		*/
		for(int j=0;j<m_numOfPoints;j++){

			double localNormal[3];
			for(int k=0;k<3;k++){
				localNormal[k]=m_filterNormals[j*3+k];
			}
			std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(j);
			int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
			std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
			distanceProjection=new double[numOfNearest];
			distanceTangent=new double[numOfNearest];
			weightTangent=new double[numOfNearest];
			weightNormal=new double[numOfNearest];
			normalDot=new double[numOfNearest];
			for(int www=0;www<numOfNearest;www++){
				distanceProjection[www]=0;
				distanceTangent[www]=0;
				weightNormal[www]=0;
				weightTangent[www]=0;
				normalDot[www]=0;
			}

			outFlag=new int[numOfNearest];
			weightSum=0;
			distanceMove=0;
			double sumDistanceProjection=0;
			double sumDistanceTangent=0;
			double sumDistanceNormal=0;
			//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonDistance.txt";
			//FILE *fpout;
			//if((fpout = fopen(filename_pw, "a")) == NULL)
			//{
			//	int dkjkd;
			//	//MessageBox("can't open the file!");
			//}
			for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){			
				/*
				*	法向的threshold处理
				*/
				if((m_filterNormals[j*3]*m_filterNormals[(*vectorNearestIterator)*3]
				+m_filterNormals[j*3+1]*m_filterNormals[(*vectorNearestIterator)*3+1]
				+m_filterNormals[j*3+2]*m_filterNormals[(*vectorNearestIterator)*3+2])<m_thresholdDistanceNormal){
					outFlag[k]=0;

					numOfNearest-=1;
					continue;
				}
				else{
					outFlag[k]=1;
				}
				double localVectorPosition[3];
				for(int w=0;w<3;w++){
					localVectorPosition[w]=m_originalPointSet[(*vectorNearestIterator)*3+w]-m_originalPointSet[j*3+w];					
				}

				double localNearestNormal[3];
				for(int w=0;w<3;w++){
					localNearestNormal[w]=m_filterNormals[(*vectorNearestIterator)*3+w];
				}
				distanceProjection[k]=Vector3Vector(localVectorPosition,localNormal);
				normalDot[k]=abs(Vector3Vector(localNormal,localNearestNormal));
				distanceTangent[k]=ComputeNormOfVector3(localVectorPosition);
				if(m_tangentOrManifold==1){
					distanceTangent[k]=sqrt(distanceTangent[k]*distanceTangent[k]-distanceProjection[k]*distanceProjection[k]);
				}
				sumDistanceTangent+=distanceTangent[k]*distanceTangent[k];
				sumDistanceProjection+=distanceProjection[k]*distanceProjection[k];
				sumDistanceNormal+=(1-normalDot[k])*(1-normalDot[k]);

				/*sumDistanceTangent+=distanceTangent[k];
				sumDistanceProjection+=abs(distanceProjection[k]);*/
				//weightNormal[k]=exp(-0.5*distanceProjection*distanceProjection/m_variationNormal);
				//weightTangent=exp(-0.5*distanceTangent*distanceTangent/m_variationTangent);
				//double tempWeight=weightNormal*weightTangent;
				//distanceMove+=distanceProjection*tempWeight;
				//weightSum+=tempWeight;



				//		fprintf(fpout," %f ",distanceProjection[k]);	


			}
			//fprintf(fpout,"\n");

			//	fclose(fpout);
		/*	sumDistanceTangent/=(numOfNearest+1);
			sumDistanceProjection/=(numOfNearest+1);*/
			double localVariationNormal,localVarationTangent,localVariationNormalDistance;

		//	localVariationNormal=m_variationNormal*sumDistanceProjection*sumDistanceProjection;
			localVariationNormal=sumDistanceProjection/numOfNearest;
			if(localVariationNormal<0.000001)
				localVariationNormal=0.000001;
		//	localVarationTangent=m_variationTangent*sumDistanceTangent*sumDistanceTangent;
			localVarationTangent=sumDistanceTangent/numOfNearest;
			if(localVarationTangent<0.000001)
				localVarationTangent=0.000001;

			localVariationNormalDistance=sumDistanceNormal/numOfNearest;
			if(localVariationNormalDistance<0.000001){
				localVariationNormalDistance=0.000001;
			}

			vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();

			if(m_ifNormalWeigtht==0){
				for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
					if(outFlag[k]==0)
						continue;
					weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
					weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);

					double tempWeight=weightNormal[k]*weightTangent[k];

					distanceMove+=distanceProjection[k]*tempWeight;
					weightSum+=tempWeight;			
				}

				distanceMove=distanceMove/(weightSum+1);
				m_distanceMove[j]=distanceMove;
				m_volumeChange[j]=distanceMove*m_meanAreas[j];
				for(int k=0;k<3;k++){
					m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];
				}

			}
			else{
				if(m_ifAreaWeight==0){
					if(m_ifVariationNormal==0){
						for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
							if(outFlag[k]==0)
								continue;
							weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
							weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);

							double tempWeight=weightNormal[k]*weightTangent[k]*normalDot[k];

							distanceMove+=distanceProjection[k]*tempWeight;
							weightSum+=tempWeight;			
						}

					}
					else{
						for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
							if(outFlag[k]==0)
								continue;
							weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
							weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);
							double tempWeightNormalDistance;
							if ((1-normalDot[k])*(1-normalDot[k])>localVariationNormalDistance) {
								tempWeightNormalDistance=0;
							}
							else
								tempWeightNormalDistance=exp(-0.5*(1-normalDot[k])*(1-normalDot[k])/localVariationNormalDistance/4);
							double tempWeight=weightNormal[k]*weightTangent[k]*tempWeightNormalDistance;

							distanceMove+=distanceProjection[k]*tempWeight;
							weightSum+=tempWeight;			
						}

					}

					

					distanceMove=distanceMove/(weightSum+1);
					m_distanceMove[j]=distanceMove;
					m_volumeChange[j]=distanceMove*m_meanAreas[j];
					for(int k=0;k<3;k++){
						m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];
					}

				}
				else{
					if(m_ifVariationNormal==0){
						for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
							if(outFlag[k]==0)
								continue;
							weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
							weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);
							double tempWeightNormalDistance=exp(-0.5*(1-distanceTangent[k])*(1-distanceTangent[k])/localVarationTangent);

							double tempWeight=weightNormal[k]*weightTangent[k]*(m_meanAreas[(*vectorNearestIterator)]/m_kNearest)*normalDot[k];

							distanceMove+=distanceProjection[k]*tempWeight;
							weightSum+=tempWeight;			
						}

					}
					else{
						for(int k=0;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++,k++){
							if(outFlag[k]==0)
								continue;
							weightNormal[k]=exp(-0.5*distanceProjection[k]*distanceProjection[k]/localVariationNormal);
							weightTangent[k]=exp(-0.5*distanceTangent[k]*distanceTangent[k]/localVarationTangent);
							double tempWeightNormalDistance;
							
							if ((1-normalDot[k])*(1-normalDot[k])>localVariationNormalDistance) {
								tempWeightNormalDistance=0;
							}
							else
								tempWeightNormalDistance=exp(-0.5*(1-normalDot[k])*(1-normalDot[k])/localVariationNormalDistance/4);
							

							double tempWeight=weightNormal[k]*weightTangent[k]*(m_meanAreas[(*vectorNearestIterator)]/m_kNearest)*tempWeightNormalDistance;

							distanceMove+=distanceProjection[k]*tempWeight;
							weightSum+=tempWeight;			
						}

					}
					

					distanceMove=distanceMove/(weightSum+m_meanAreas[j]/m_kNearest);
					m_distanceMove[j]=distanceMove;
					m_volumeChange[j]=distanceMove*m_meanAreas[j];
					for(int k=0;k<3;k++){
						m_resultPointSet[j*3+k]=m_originalPointSet[j*3+k]+distanceMove*localNormal[k];
					}

				}
			}
			
			delete distanceProjection;
			delete distanceTangent;
			delete weightTangent;
			delete weightNormal;
			delete outFlag;
			delete normalDot;



		}
		if(m_ifVolumePreserve==1){
			RestoreVolume();
		}
	
		for(int j=0;j<m_numOfPoints*3;j++){
			m_originalPointSet[j]=m_resultPointSet[j];
		}

	}
	//	m_mapKnearest.clear();

	//delete[] m_originalPointSet;
	//m_originalPointSet=NULL;
	return;

}