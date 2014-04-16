// Bilateral_Filter.cpp: implementation of the Bilateral_Filter class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Mcube.h"
#include "Bilateral_Filter.h"
#include ".\mathlib\mathlib.h"
#include "Matrix.h"

using namespace MATHLIB;


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Bilateral_Filter::Bilateral_Filter()
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
}

Bilateral_Filter::~Bilateral_Filter()
{
	DeleteFilterBilateral();	
}

void Bilateral_Filter::GetFilterBilateral(int numOfPoints, float* pointSet, float* colors)
{
	m_numOfPoints=numOfPoints;
	m_originalPointSet=new float[m_numOfPoints*3];
	m_originalColors=new float[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints;i++){
		m_originalPointSet[i]=pointSet[i];
		m_originalColors[i]=colors[i];
	}
	NormalPosition();
	ComputeMapKnearest();
	CalculateOriginalNormals();
	MeanShiftNormals();
	//CalculateColorGradient();
	CalculateBilateralFilterOne();
	//CalculateBilateralFilterTwo();
	m_mapKnearest.clear();
	//for(int i=0;i<3;i++){
	//	delete[] m_originalColorGradients[i];
	//	m_originalColorGradients[i]=NULL;
	//}
	RecoverPosition();
	return;

}

void Bilateral_Filter::DeleteFilterBilateral()
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
void Bilateral_Filter::ComputeMapKnearest()
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
		//////////////////////////////////////////////////////////////////////////
		// 测试代码
        mapIterator=mapKnearest.begin();
		for(;mapIterator!=mapKnearest.end();mapIterator++){
			float k;
			int g;
			k=(*mapIterator).first;
			g=(*mapIterator).second;
			k++;
			g++;
		}
		//////////////////////////////////////////////////////////////////////////
		mapKnearest.clear();
	}	
    //////////////////////////////////////////////////////////////////////////
    // 测试代码
	std::map<int,KnearestField>::iterator mapTempIterator=m_mapKnearest.begin();
	for(;mapTempIterator!=m_mapKnearest.end();mapTempIterator++){
		std::vector<int>::iterator tempTempIterator=(*mapTempIterator).second.m_nearest.begin();
        
		int k=(*mapTempIterator).first;
		int g=(*tempTempIterator);
		g++;
		k++;
	}	
	//////////////////////////////////////////////////////////////////////////
}

void Bilateral_Filter::CalculateOriginalNormals()
{
	
	float centroidPosition[6];//重心位置
	float* localVariationVector;//重心到点的位置的向量
	int numOfNearest;

	mat_f8 mat_covaMatrix(6, 6);
	m_originalNormals=new float[m_numOfPoints*6];
	
	for(int i=0;i<m_numOfPoints;i++){
		std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		numOfNearest=(*mapIterator).second.m_numOfNearest;
		for(int j=0;j<3;j++){
			centroidPosition[j]=m_originalPointSet[i*3+j];
		}
	    centroidPosition[3]=m_originalColors[i*3];
		centroidPosition[4]=m_originalColors[i*3+1];
		centroidPosition[5]=m_originalColors[i*3+2];
		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
			centroidPosition[0]+=m_originalPointSet[(*vectorIterator)*3];
			centroidPosition[1]+=m_originalPointSet[(*vectorIterator)*3+1];
			centroidPosition[2]+=m_originalPointSet[(*vectorIterator)*3+2];
			centroidPosition[3]+=m_originalColors[(*vectorIterator)*3];
			centroidPosition[4]+=m_originalColors[(*vectorIterator)*3+1];
			centroidPosition[5]+=m_originalColors[(*vectorIterator)*3+2];
		}
		for(int j=0;j<6;j++){
			centroidPosition[j]/=(numOfNearest+1);
		}
		localVariationVector=new float[(numOfNearest+1)*6];
		localVariationVector[0]=m_originalPointSet[i*3]-centroidPosition[0];
		localVariationVector[1]=m_originalPointSet[i*3+1]-centroidPosition[1];
		localVariationVector[2]=m_originalPointSet[i*3+2]-centroidPosition[2];
		localVariationVector[3]=m_originalColors[i*3]-centroidPosition[3];
		localVariationVector[4]=m_originalColors[i*3+1]-centroidPosition[4];
		localVariationVector[5]=m_originalColors[i*3+2]-centroidPosition[5];
		vectorIterator=(*mapIterator).second.m_nearest.begin();
		for(int j=1;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,j++){
			localVariationVector[j*3]=m_originalPointSet[(*vectorIterator)*3]-centroidPosition[0];
			localVariationVector[j*3+1]=m_originalPointSet[(*vectorIterator)*3+1]-centroidPosition[1];
			localVariationVector[j*3+2]=m_originalPointSet[(*vectorIterator)*3+2]-centroidPosition[2];
			localVariationVector[j*3+3]=m_originalColors[(*vectorIterator)*3]-centroidPosition[3];
			localVariationVector[j*3+4]=m_originalColors[(*vectorIterator)*3+1]-centroidPosition[4];
			localVariationVector[j*3+5]=m_originalColors[(*vectorIterator)*3+2]-centroidPosition[5];
		}
		//求方差矩阵
		for(int j=0;j<6;j++){
			for(int k=0;k<6;k++){
				mat_covaMatrix(j,k)=0;
			
				for(int m=0;m<(numOfNearest+1);m++){
					mat_covaMatrix(j,k)=mat_covaMatrix(j,k)+localVariationVector[m*3+j]*localVariationVector[m*3+k];
				}
			}
		}
		vec_f8 eval(6);
		mat_f8 evec(6, 6);
		eigen_symm(eval, evec, mat_covaMatrix);
		m_originalNormals[i*6]=evec(0,5);
		m_originalNormals[i*6+1]=evec(1,5);
		m_originalNormals[i*6+2]=evec(2,5);
		m_originalNormals[i*6+3]=evec(3,5);
		m_originalNormals[i*6+4]=evec(4,5);
		m_originalNormals[i*6+5]=evec(5,5);
		delete[] localVariationVector;
	}
}

void Bilateral_Filter::MeanShiftNormals()
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

    m_filterNormals=new float[m_numOfPoints*6];
	while(num_stop<m_numOfPoints){
		for(int i=0;i<m_numOfPoints;i++){
			if(!flag_stop[i]){
				sum_kernel=0;
				for(int j=0;j<6;j++){
                    m_filterNormals[i*3+j]=0;
				}
			
				std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(i);
				//int numOfNearest=(*mapKnearestIterator).second.m_nearest.size();
				std::vector<int>::iterator vectorNearestIterator=(*mapKnearestIterator).second.m_nearest.begin();
				for(;vectorNearestIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorNearestIterator++){
					distance_position=(m_originalPointSet[i*3]-m_originalPointSet[(*vectorNearestIterator)*3])*(m_originalPointSet[i*3]-m_originalPointSet[(*vectorNearestIterator)*3])
						+(m_originalPointSet[i*3+1]-m_originalPointSet[(*vectorNearestIterator)*3+1])*(m_originalPointSet[i*3+1]-m_originalPointSet[(*vectorNearestIterator)*3+1])
						+(m_originalPointSet[i*3+2]-m_originalPointSet[(*vectorNearestIterator)*3+2])*(m_originalPointSet[i*3+2]-m_originalPointSet[(*vectorNearestIterator)*3+2])
						+(m_originalColors[i*3]-m_originalColors[(*vectorNearestIterator)*3])*(m_originalColors[i*3]-m_originalColors[(*vectorNearestIterator)*3])
						+(m_originalColors[i*3+1]-m_originalColors[(*vectorNearestIterator)*3+1])*(m_originalColors[i*3+1]-m_originalColors[(*vectorNearestIterator)*3+1])
						+(m_originalColors[i*3+2]-m_originalColors[(*vectorNearestIterator)*3+2])*(m_originalColors[i*3+2]-m_originalColors[(*vectorNearestIterator)*3+2]);
						
					distance_normal=(m_originalNormals[i*6]-m_originalNormals[(*vectorNearestIterator)*6])*(m_originalNormals[i*6]-m_originalNormals[(*vectorNearestIterator)*6])
						+(m_originalNormals[i*6+1]-m_originalNormals[(*vectorNearestIterator)*6+1])*(m_originalNormals[i*6+1]-m_originalNormals[(*vectorNearestIterator)*6+1])
						+(m_originalNormals[i*6+2]-m_originalNormals[(*vectorNearestIterator)*6+2])*(m_originalNormals[i*6+2]-m_originalNormals[(*vectorNearestIterator)*6+2])
						+(m_originalNormals[i*6+3]-m_originalNormals[(*vectorNearestIterator)*6+3])*(m_originalNormals[i*6+3]-m_originalNormals[(*vectorNearestIterator)*6+3])
						+(m_originalNormals[i*6+4]-m_originalNormals[(*vectorNearestIterator)*6+4])*(m_originalNormals[i*6+4]-m_originalNormals[(*vectorNearestIterator)*6+4])
						+(m_originalNormals[i*6+5]-m_originalNormals[(*vectorNearestIterator)*6+5])*(m_originalNormals[i*6+5]-m_originalNormals[(*vectorNearestIterator)*6+5]);
				
					position_kernel=exp(-0.5*distance_position/m_meanShiftHposition);
					range_kernel=exp(-0.5*distance_normal/m_meanShiftHnormal);
					sum_kernel+=position_kernel*range_kernel;
					
					for(int j=0;j<6;j++){
						m_filterNormals[i*6+j]+=m_originalNormals[(*vectorNearestIterator)*6+j]*position_kernel*range_kernel;
					}
				}
				sum_kernel+=exp(0.0)*exp(0.0);
				for(int j=0;j<6;j++){
					m_filterNormals[i*6+j]+=m_originalNormals[i*6+j]*exp(0.0)*exp(0.0);
				}
				for(int j=0;j<6;j++){
					m_filterNormals[i*6+j]/=sum_kernel;
				}
				//归一化 m_originalFilterNormals
				float tempNormal[6];
				for(int j=0;j<6;j++){
					tempNormal[j]=m_filterNormals[i*6+j];
				}
				length_normal_vector=sqrt(
					tempNormal[0]*tempNormal[0]
					+tempNormal[1]*tempNormal[1]
					+tempNormal[2]*tempNormal[2]
					+tempNormal[3]*tempNormal[3]
					+tempNormal[4]*tempNormal[4]
					+tempNormal[5]*tempNormal[5]);

				for(int j=0;j<6;j++){
					m_filterNormals[i*6+j]/=length_normal_vector;
				}

				length_normal_vector=(m_filterNormals[i*6]-m_originalNormals[i*6])*(m_filterNormals[i*6]-m_originalNormals[i*6])
					+(m_filterNormals[i*6+1]-m_originalNormals[i*6+1])*(m_filterNormals[i*6+1]-m_originalNormals[i*6+1])
					+(m_filterNormals[i*6+2]-m_originalNormals[i*6+2])*(m_filterNormals[i*6+2]-m_originalNormals[i*6+2])
					+(m_filterNormals[i*6+3]-m_originalNormals[i*6+3])*(m_filterNormals[i*6+3]-m_originalNormals[i*6+3])
					+(m_filterNormals[i*6+4]-m_originalNormals[i*6+4])*(m_filterNormals[i*6+4]-m_originalNormals[i*6+4])
					+(m_filterNormals[i*6+5]-m_originalNormals[i*6+5])*(m_filterNormals[i*6+5]-m_originalNormals[i*6+5]);
				
				if(length_normal_vector<m_meanShiftStopNormal){
						flag_stop[i]=true;
						num_stop+=1;
					}
				else{
					m_originalNormals[i*6]=m_filterNormals[i*6];
					m_originalNormals[i*6+1]=m_filterNormals[i*6+1];
					m_originalNormals[i*6+2]=m_filterNormals[i*6+2];
					m_originalNormals[i*6+3]=m_filterNormals[i*6+3];
					m_originalNormals[i*6+4]=m_filterNormals[i*6+4];
					m_originalNormals[i*6+5]=m_filterNormals[i*6+5];					
				}
				
			}
		}

	}
	delete[] m_originalNormals;
}

void Bilateral_Filter::CalculateColorGradient()
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
		//	coffecientConstant[j]=0;
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
		coffecientMatrix[0][3]+=localEXOne;
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
			coffecientMatrix[0][3]+=localEXOne;
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
			m_originalColorGradients[0][i*4+j]=m_gradientRed[0]*vectorOne[j]+m_gradientRed[1]*vectorTwo[j];
			m_originalColorGradients[1][i*4+j]=m_gradientGreen[0]*vectorOne[j]+m_gradientGreen[1]*vectorTwo[j];
			m_originalColorGradients[2][i*4+j]=m_gradientBlue[0]*vectorOne[j]+m_gradientBlue[1]*vectorTwo[j];
		}
		m_originalColorGradients[0][i*4+3]=m_gradientRed[2];
		m_originalColorGradients[1][i*4+3]=m_gradientGreen[2];
		m_originalColorGradients[2][i*4+3]=m_gradientBlue[2];
		
	}
	return;
}

void Bilateral_Filter::CalculateBilateralFilterOne()
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
//	float sumKernel_colors[3];
	for(int i=0;i<m_numOfPoints;i++){
		std::map<int ,KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		
		sumKernel_position=0;

		ComputeOneOrderPositionEstimation(i);
		CalculateVariation();
		std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();
		std::vector<float>::iterator localProjectionDistanceIterator=m_localProjectionDistance.begin();
		std::vector<POINTVECTOR6D>::iterator localPositionEstimationIterator=m_localPositionEstimation.begin();
	
		for(;localPointDistanceIterator!=m_localPointDistance.end();localPointDistanceIterator++){
			float kernel_distance;
			kernel_distance=exp(-0.5*(*localPointDistanceIterator)/m_pointVariation)*exp(-0.5*(*localProjectionDistanceIterator)/m_estimationVariation);
		
			sumKernel_position+=kernel_distance;

			for(int j=0;j<3;j++){
				m_resultPointSet[i*3+j]+=(*localPositionEstimationIterator).pointVector[j]*kernel_distance;
				
				m_resultColors[i*3+j]+=(*localPositionEstimationIterator).pointVector[3+j]*kernel_distance;
			}
			localProjectionDistanceIterator++;
		    localPositionEstimationIterator++;
			
		}

		sumKernel_position+=1;

		for(int j=0;j<3;j++){
			m_resultPointSet[i*3+j]+=m_originalPointSet[i*3+j];
			m_resultPointSet[i*3+j]/=sumKernel_position;
			m_resultColors[i*3+j]+=m_originalColors[i*3+j];
			m_resultColors[i*3+j]/=sumKernel_position;
		}
		m_localPointDistance.clear();
        m_localProjectionDistance.clear();
        m_localPositionEstimation.clear();
      
	}
	delete[] m_originalPointSet;
	delete[] m_originalColors;

}

void Bilateral_Filter::CalculateBilateralFilterTwo()
{
		//常数假定但加入了一个梯度聚类
	/*m_resultPointSet=new float[m_numOfPoints*3];
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
		std::vector<POINTVECTOR6D>::iterator localPositionEstimationIterator=m_localPositionEstimation.begin();
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
			kernel_color[0]=kernel_distance
				*exp(-0.5*(*localColorDistance[0])/m_colorVariation[0])
				*exp(-0.5*(*localColorProjectionDistance[0])/m_colorEstimationVariation[0]);
			kernel_color[1]=kernel_distance
				*exp(-0.5*(*localColorDistance[1])/m_colorVariation[1])
				*exp(-0.5*(*localColorProjectionDistance[1])/m_colorEstimationVariation[1]);
			kernel_color[2]=kernel_distance
				*exp(-0.5*(*localColorDistance[2])/m_colorVariation[2])
				*exp(-0.5*(*localColorProjectionDistance[2])/m_colorEstimationVariation[2]);
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
	delete[] m_originalColors;*/


}

void Bilateral_Filter::CalculateBilateralFilterThree()
{

}

void Bilateral_Filter::ComputeOneOrderPositionEstimation(int indexPoint)
{
	std::map<int,KnearestField>::iterator mapIterator=m_mapKnearest.find(indexPoint);
	std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	float localNormal[6];//局部法向
	localNormal[0]=m_filterNormals[indexPoint*3];
	localNormal[1]=m_filterNormals[indexPoint*3+1];
	localNormal[2]=m_filterNormals[indexPoint*3+2];
	localNormal[3]=m_filterNormals[indexPoint*3+3];
	localNormal[4]=m_filterNormals[indexPoint*3+4];
	localNormal[5]=m_filterNormals[indexPoint*3+5];

	for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
	   float distance_point_point;//点到点的距离
	   float distance_point_projection;//点到一次估计平面的距离
	   float localVector_point_point[6];//点到点的向量
	   float localColorDistance[3];//颜色距离
	   
	   POINTVECTOR6D localProjection;//投影点的位置

	   localVector_point_point[0]=m_originalPointSet[(*vectorIterator)*3]-m_originalPointSet[indexPoint*3];
	   localVector_point_point[1]=m_originalPointSet[(*vectorIterator)*3+1]-m_originalPointSet[indexPoint*3+1];
	   localVector_point_point[2]=m_originalPointSet[(*vectorIterator)*3+2]-m_originalPointSet[indexPoint*3+2];
	   localVector_point_point[3]=m_originalColors[(*vectorIterator)*3]-m_originalColors[indexPoint*3];
	   localVector_point_point[4]=m_originalColors[(*vectorIterator)*3+1]-m_originalColors[indexPoint*3+1];
	   localVector_point_point[5]=m_originalColors[(*vectorIterator)*3+2]-m_originalColors[indexPoint*3+2];

	
	   distance_point_projection=0;
	   distance_point_point=0;
	   for(int i=0;i<6;i++){
		   distance_point_projection+=localVector_point_point[i]*localNormal[i];
		   distance_point_point+=localVector_point_point[i]*localVector_point_point[i];
	   }

	  

	   localProjection.pointVector[0]=m_originalPointSet[indexPoint*3]+distance_point_projection*localNormal[0];
	   localProjection.pointVector[1]=m_originalPointSet[indexPoint*3+1]+distance_point_projection*localNormal[1];
	   localProjection.pointVector[2]=m_originalPointSet[indexPoint*3+2]+distance_point_projection*localNormal[2];
	   localProjection.pointVector[3]=m_originalColors[indexPoint*3]+distance_point_projection*localNormal[3];
	   localProjection.pointVector[4]=m_originalColors[indexPoint*3+1]+distance_point_projection*localNormal[4];
	   localProjection.pointVector[5]=m_originalColors[indexPoint*3+2]+distance_point_projection*localNormal[5];

	   m_localProjectionDistance.push_back(distance_point_projection*distance_point_projection);
	   m_localPointDistance.push_back(distance_point_point);
	   m_localPositionEstimation.push_back(localProjection);
	
	}
	//计算一次颜色估计和颜色投影距离
//	ComputeOneOrderColorEstimation(indexPoint);
    return;
}

void Bilateral_Filter::ComputeOneOrderColorEstimation(int indexPoint)
{
	/*std::map<int,KnearestField>::iterator mapIterator=m_mapKnearest.find(indexPoint);
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
	return;*/
}

void Bilateral_Filter::CalculateVariation()
{
	// 
	std::vector<float>::iterator localPointDistanceIterator=m_localPointDistance.begin();
	std::vector<float>::iterator localProjectionDistanceIterator=m_localProjectionDistance.begin();

	m_pointVariation=0;
	m_estimationVariation=0;

	for(;localPointDistanceIterator!=m_localPointDistance.end();localPointDistanceIterator++,localProjectionDistanceIterator++){
		m_pointVariation+=(*localPointDistanceIterator);
	    m_estimationVariation+=(*localProjectionDistanceIterator);
	
	}
	m_pointVariation/=m_localPointDistance.size();
	m_estimationVariation/=m_localPointDistance.size();

	return;
}

void Bilateral_Filter::NormalPosition()
{
	float max_minus_min;
	max_minus_min=m_maxPosition-m_minPosition;
	for(int i=0;i<m_numOfPoints*3;i++){
		m_originalPointSet[i]=(m_originalPointSet[i]-m_minPosition)/max_minus_min;
	}
	return;
}

void Bilateral_Filter::RecoverPosition()
{
	float max_minus_min=m_maxPosition-m_minPosition;
	for(int i=0;i<m_numOfPoints*3;i++){
		m_resultPointSet[i]=m_resultPointSet[i]*max_minus_min+m_minPosition;
	}
}


