#include "StdAfx.h"
#include ".\geometricestimation.h"
#include <stdlib.h>
#include <math.h>
#include "mathlib\mathlib.h"
#include <map>
#include "Matrix.h"
using namespace MATHLIB;

GeometricEstimation::GeometricEstimation(int numOfPointSet,float* pointSet,uint8 numOfKnearest,int* allPointsNearestNeighbor,float* allPointsNormal)
{
	m_numOfPoint=numOfPointSet;
	m_numOfNeighbor=numOfKnearest;
	m_pointSet=new float[m_numOfPoint*3];
	m_allPointsNormal=new float[m_numOfPoint*3];
	m_allPointsNeighbor=new int[m_numOfNeighbor*m_numOfPoint];
	for(int i=0;i<m_numOfPoint*3;i++){
		m_pointSet[i]=pointSet[i];
		m_allPointsNormal[i]=allPointsNormal[i];
	}
	for(int i=0;i<m_numOfPoint*m_numOfNeighbor;i++){
		m_allPointsNeighbor[i]=allPointsNearestNeighbor[i];
	}
	m_allPointGeometry=NULL;
	m_allPointsAverageCurvature=NULL;
	m_allPointsGaussCurvature=NULL;
	m_allPointsPricipalDirectionOne=NULL;
	m_allPointsPricipalDirectionTwo=NULL;
	m_allPointsPricipalCurvatureOne=NULL;
	m_allPointsPricipalCurvatureTwo=NULL;

}

GeometricEstimation::GeometricEstimation(int numOfPointSet,float* pointSet)
{
	m_numOfPoint=numOfPointSet;	
	
	m_pointSet=new float[m_numOfPoint*3];
	for(int i=0;i<m_numOfPoint;i++){
		for(int j=0;j<3;j++){
			m_pointSet[i*3+j]=pointSet[i*3+j];
		}
	}
	m_allPointsNeighbor=NULL;
	m_allPointsNormal=NULL;
	m_allPointGeometry=NULL;
	m_allPointsAverageCurvature=NULL;
	m_allPointsGaussCurvature=NULL;
	m_allPointsPricipalDirectionOne=NULL;
	m_allPointsPricipalDirectionTwo=NULL;
	m_allPointsPricipalCurvatureOne=NULL;
	m_allPointsPricipalCurvatureTwo=NULL;

	

}
GeometricEstimation::GeometricEstimation(int numOfPointSet,float* pointSet,UINT8 numOfKnearest,int* allPointsNearestNeighbor)
{
	m_numOfPoint=numOfPointSet;
	m_pointSet=new float[m_numOfPoint*3];
	for(int i=0;i<m_numOfPoint;i++){
		for(int j=0;j<3;j++)
			m_pointSet[i*3+j]=pointSet[i*3+j];
	}
	
	m_numOfNeighbor=numOfKnearest;
	m_allPointsNeighbor=new int[m_numOfNeighbor*m_numOfPoint];
	for(int i=0;i<m_numOfNeighbor*m_numOfPoint;i++){
		m_allPointsNeighbor[i]=allPointsNearestNeighbor[i];
	}
	m_allPointGeometry=NULL;
	m_allPointsNormal=NULL;
	m_allPointsAverageCurvature=NULL;
	m_allPointsGaussCurvature=NULL;
	m_allPointsPricipalDirectionOne=NULL;
	m_allPointsPricipalDirectionTwo=NULL;
	m_allPointsPricipalCurvatureOne=NULL;
	m_allPointsPricipalCurvatureTwo=NULL;
}
GeometricEstimation::~GeometricEstimation(void)
{
	if(m_allPointsNeighbor!=NULL)
		delete[] m_allPointsNeighbor;
	m_allPointsNeighbor=NULL;
	if(m_pointSet!=NULL)
		delete[] m_pointSet;
	m_pointSet=NULL;
	if(m_allPointsNormal!=NULL)
		delete[] m_allPointsNormal;
	m_allPointsNormal=NULL;
	if(m_allPointGeometry!=NULL)
		delete[] m_allPointGeometry;
	m_allPointGeometry=NULL;
    if(m_allPointsAverageCurvature!=NULL)
		delete[] m_allPointsAverageCurvature;
	m_allPointsAverageCurvature=NULL;
	if(m_allPointsGaussCurvature!=NULL)
		delete[] m_allPointsGaussCurvature;	
	m_allPointsGaussCurvature=NULL;
	if(m_allPointsPricipalCurvatureOne!=NULL)
		delete[] m_allPointsPricipalCurvatureOne;
	m_allPointsPricipalCurvatureOne=NULL;
	if(m_allPointsPricipalDirectionOne!=NULL)
		delete[] m_allPointsPricipalDirectionOne;
	m_allPointsPricipalDirectionOne=NULL;

	if(m_allPointsPricipalDirectionTwo!=NULL)
		delete[] m_allPointsPricipalDirectionTwo;
	m_allPointsPricipalDirectionTwo=NULL;
    if(m_allPointsPricipalCurvatureTwo!=NULL)
		delete[] m_allPointsPricipalCurvatureTwo;
	m_allPointsPricipalCurvatureTwo=NULL;
}

void GeometricEstimation::Delete()
{
	if(m_allPointsNeighbor!=NULL)
		delete[] m_allPointsNeighbor;
	m_allPointsNeighbor=NULL;
	if(m_pointSet!=NULL)
		delete[] m_pointSet;
	m_pointSet=NULL;
	if(m_allPointsNormal!=NULL)
		delete[] m_allPointsNormal;
	m_allPointsNormal=NULL;
	if(m_allPointGeometry!=NULL)
		delete[] m_allPointGeometry;
	m_allPointGeometry=NULL;

	if(m_allPointsAverageCurvature!=NULL)
		delete[] m_allPointsAverageCurvature;
	m_allPointsAverageCurvature=NULL;
	if(m_allPointsGaussCurvature!=NULL)
		delete[] m_allPointsGaussCurvature;	
	m_allPointsGaussCurvature=NULL;
	if(m_allPointsPricipalCurvatureOne!=NULL)
		delete[] m_allPointsPricipalCurvatureOne;
	m_allPointsPricipalCurvatureOne=NULL;
	if(m_allPointsPricipalDirectionOne!=NULL)
		delete[] m_allPointsPricipalDirectionOne;
	m_allPointsPricipalDirectionOne=NULL;

	if(m_allPointsPricipalDirectionTwo!=NULL)
		delete[] m_allPointsPricipalDirectionTwo;
	m_allPointsPricipalDirectionTwo=NULL;
	if(m_allPointsPricipalCurvatureTwo!=NULL)
		delete[] m_allPointsPricipalCurvatureTwo;
	m_allPointsPricipalCurvatureTwo=NULL;

}
//////////////////////////////////////////////////////////////////////////
//邻域搜索  numofPoint 点的数量，pointset  点的集合， ID 搜索点的index， knearest  邻域数目， radius  邻域球的半径
void GeometricEstimation::SearchNeighbor(UINT8 knearest,float radius,int* pointNeighbor)
{
	
	//int* pointNeighbor;// 存放 近邻点的ID
	float* knearestDistance;// 存放近邻点到所求点的距离
	radius=radius*radius;
	//pointNeighbor=new int[knearest]; 
	knearestDistance=new float[knearest];	
	//////////////////////////////////////////////////////////////////////////
	// 初始化
	for(int i=0;i<knearest;i++){
		pointNeighbor[i]=-1;
	}

	for(int i=0;i<knearest;i++){
		knearestDistance[i]=999999;
	}
	//////////////////////////////////////////////////////////////////////////
	// 遍历所有点，并且按照距离由小到大存放所取到的点
	
	for(int i=0;i<m_numOfPoint;i++){
		if(i==m_ID)
			continue;
		float distancePointToPoint;
		distancePointToPoint=(m_pointSet[m_ID*3]-m_pointSet[i*3])*(m_pointSet[m_ID*3]-m_pointSet[i*3])
			+(m_pointSet[m_ID*3+1]-m_pointSet[i*3+1])*(m_pointSet[m_ID*3+1]-m_pointSet[i*3+1])
			+(m_pointSet[m_ID*3+2]-m_pointSet[i*3+2])*(m_pointSet[m_ID*3+2]-m_pointSet[i*3+2]);
        if(distancePointToPoint>radius)
			continue;
		// 排序和插入新点过程
		for(int j=0;j<knearest;j++){
			if(distancePointToPoint<knearestDistance[j]){
				for(int k=knearest-1;k>j;k--){
					knearestDistance[k]=knearestDistance[k-1];
					pointNeighbor[k]=pointNeighbor[k-1];
				}
				knearestDistance[j]=distancePointToPoint;
				pointNeighbor[j]=i;
				break;
			}
		}               
	}
	delete[] knearestDistance;
	return ;

}
//确定近邻搜索的半径
float GeometricEstimation::EstimateSearchingRadius()
{   
	#define RAND_MAX m_numOfPoint-1;
	float radius;
	int randomNumber;
	int maxIndexNum=15;//估计近邻半径的抽样统计数
	float* nearestDistance;

	float distancePointToPoint;
	nearestDistance=new float[maxIndexNum];
	for(int i=0;i<maxIndexNum;i++){
		nearestDistance[i]=999999;
	}
	for(int i=0;i<maxIndexNum;i++){
		srand(i);// 产生随机数
		randomNumber=rand();// 产生随机数
		// 对每个点求取最小距离点
		for(int j=0;j<randomNumber;j++){
			distancePointToPoint=(m_pointSet[randomNumber*3]-m_pointSet[j*3])*(m_pointSet[randomNumber*3]-m_pointSet[j*3])
				+(m_pointSet[randomNumber*3+1]-m_pointSet[j*3+1])*(m_pointSet[randomNumber*3+1]-m_pointSet[j*3+1])
				+(m_pointSet[randomNumber*3+2]-m_pointSet[j*3+2])*(m_pointSet[randomNumber*3+2]-m_pointSet[j*3+2]);
			if(distancePointToPoint<nearestDistance[i])
				nearestDistance[i]=distancePointToPoint;
		}
			for(int j=randomNumber+1;j<m_numOfPoint;j++){
			distancePointToPoint=(m_pointSet[randomNumber*3]-m_pointSet[j*3])*(m_pointSet[randomNumber*3]-m_pointSet[j*3])
				+(m_pointSet[randomNumber*3+1]-m_pointSet[j*3+1])*(m_pointSet[randomNumber*3+1]-m_pointSet[j*3+1])
				+(m_pointSet[randomNumber*3+2]-m_pointSet[j*3+2])*(m_pointSet[randomNumber*3+2]-m_pointSet[j*3+2]);
			if(distancePointToPoint<nearestDistance[i])
				nearestDistance[i]=distancePointToPoint;
		}
	}
    float maxDistance=nearestDistance[0];
	for(int i=1;i<maxIndexNum;i++){
		if(maxDistance<nearestDistance[i])
			maxDistance=nearestDistance[i];
	}
	radius=sqrt(maxDistance)*3;
	delete[] nearestDistance;
	return radius;
}
//法向估计
void GeometricEstimation::NormalEstimation()
{   
	//////////////////////////////////////////////////////////////////////////
	// 去掉邻域中为－1 的部分
    uint8 numOfExactNeighbor;
	numOfExactNeighbor=m_numOfNeighbor;

	for(int i=0;i<m_numOfNeighbor;i++){
		if(m_pointNeighbor[i]==-1){
			numOfExactNeighbor=(uint8)i;
			break;
		}
	}
	
	
	//////////////////////////////////////////////////////////////////////////
	// 通过求共生矩阵的最小特征值来确定法向
	float covaMatrix[3][3];// 共生矩阵
	float* PointToPointVectorMatrix;// 点到点的向量矩阵
	float eigenValue[3];// 特征值
	float eigenVector[3];// 法向的特征向量
	float averagePoint[3];//平均点位置


	PointToPointVectorMatrix=new float[(numOfExactNeighbor+1)*3];
	//////////////////////////////////////////////////////////////////////////
	// 计算平均点位置
	for(int i=0;i<3;i++){
		averagePoint[i]=m_pointSet[m_ID*3+i];
	}
	for(int i=0;i<numOfExactNeighbor;i++){
		int k=m_pointNeighbor[i];
		for(int j=0;j<3;j++){
            averagePoint[j]+=m_pointSet[k*3+j];
		}
	}

	for(int i=0;i<3;i++){
		averagePoint[i]=averagePoint[i]/(numOfExactNeighbor+1);
	}

	//////////////////////////////////////////////////////////////////////////
	// 计算向量矩阵
	

	for(int i=0;i<numOfExactNeighbor;i++){
		for(int j=0;j<3;j++){
			int k=m_pointNeighbor[i];
			PointToPointVectorMatrix[i*3+j]=m_pointSet[k*3+j]-averagePoint[j];
		}
	}

	for(int i=0;i<3;i++){
		PointToPointVectorMatrix[numOfExactNeighbor*3+i]=m_pointSet[m_ID*3+i]-averagePoint[i];
	}
	

	//////////////////////////////////////////////////////////////////////////
	// 计算共生矩阵
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			covaMatrix[i][j]=0;
			for(int k=0;k<(numOfExactNeighbor+1);k++){
                covaMatrix[i][j]+=PointToPointVectorMatrix[k*3+i]*PointToPointVectorMatrix[k*3+j];
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// 求特征值
	mat_f8 mat_covaMatrix(3, 3);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			mat_covaMatrix(i, j) = covaMatrix[i][j];
		}
	}
	vec_f8 eval(3);
	mat_f8 evec(3, 3);
	eigen_symm(eval, evec, mat_covaMatrix);
    //////////////////////////////////////////////////////////////////////////
    //求特征值特征向量
	for(int i=0;i<3;i++){
		eigenValue[i]=eval(i);
		//for(int j=0;j<3;j++){
		//	eigenVector[i][j]=evec(i,j);
		//}
	}
	for(int i=0;i<3;i++){
		eigenVector[i]=evec(i,2);
	}

	float sum=0;
	for (int i=0;i<3;i++){
		sum+=eigenVector[i]*eigenVector[i];
	}
	sum=sqrt(sum);
	for(int i=0;i<3;i++){
		m_normal[i]=eigenVector[i]/sum;
	}

	delete[] PointToPointVectorMatrix;
	return;    
}
//////////////////////////////////////////////////////////////////////////
//曲率估计  参考文章 《estimating the tensor of curvature of a surface from a polyhedral approximation》
void GeometricEstimation::CurvatureEstimation()
{
	//////////////////////////////////////////////////////////////////////////
	// 去掉邻域中为－1 的部分
	uint8 numOfExactNeighbor;
	numOfExactNeighbor=m_numOfNeighbor;

	for(int i=0;i<m_numOfNeighbor;i++){
		if(m_pointNeighbor[i]==-1){
			numOfExactNeighbor=(uint8)i;
			break;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// 计算 householder transformation matrix
    float EONE[3],EminusNvi[3],EaddNvi[3],Wvi[3],householderMatrix[3][3];
	float identityMatrix[3][3],WviTWvi[3][3];
	float normOfEminusNvi,normOfEaddNvi;
	EONE[0]=1;EONE[1]=0;EONE[2]=0;
	for(int i=0;i<3;i++)
	{
		EminusNvi[i]=EONE[i]-m_normal[i];
		EaddNvi[i]=EONE[i]+m_normal[i];
	}
	normOfEminusNvi=0;
	normOfEaddNvi=0;
	for(int i=0;i<3;i++)
	{
		normOfEminusNvi+=EminusNvi[i]*EminusNvi[i];
		normOfEaddNvi+=EaddNvi[i]*EaddNvi[i];
	}
	if(normOfEaddNvi>normOfEminusNvi)
	{
		normOfEaddNvi=sqrt(normOfEaddNvi);
		for(int i=0;i<3;i++)
		{
			Wvi[i]=EaddNvi[i]/normOfEaddNvi;
		}
	}
	else
	{
		normOfEminusNvi=sqrt(normOfEminusNvi);
		for(int i=0;i<3;i++)
		{
			Wvi[i]=EminusNvi[i]/normOfEminusNvi;
		}
	}
	Vector3TVector(Wvi,Wvi,WviTWvi);
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			identityMatrix[i][j]=0;
	for(int i=0;i<3;i++)
		identityMatrix[i][i]=1;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			householderMatrix[i][j]=identityMatrix[i][j]-2*WviTWvi[i][j];
		}
	}
	//////////////////////////////////////////////////////////////////////////
	// 估计矩阵 Mvi；
	float MviMatrix[3][3];//矩阵mvi
	float* coffOfWeight;//加权系数
	float directionCurvature;//方向曲率
	float tangentDirectionVector[3];//切向方向向量
	VECTOR3D* directionVector;//方向向量
	float* distanceDirectionVector;// 向量长度
	float* backDistanceDirectionVector;//向量长度的倒数
	float matrixNorTNor[3][3];
	float matrixTanTTan[3][3];
	float tempVector[3];
	float tempMatrix[3][3];
	float tempScalar;

	coffOfWeight=new float[numOfExactNeighbor];
	directionVector=new VECTOR3D[numOfExactNeighbor];
	distanceDirectionVector=new float[numOfExactNeighbor];
	backDistanceDirectionVector=new float[numOfExactNeighbor];

	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			MviMatrix[i][j]=0;
		}
	}

	for(int i=0,k=0;i<numOfExactNeighbor;i++,k++){
        int IDofNeighbor;
		IDofNeighbor=m_pointNeighbor[k];
		for(int j=0;j<3;j++){
			directionVector[i][j]=m_pointSet[IDofNeighbor*3+j]-m_pointSet[m_ID*3+j];
		}
		distanceDirectionVector[i]=ComputeNormOfVector3(directionVector[i]);
		if(distanceDirectionVector[i]==0){
			numOfExactNeighbor-=1;
			i=i-1;
			continue;
		}
		backDistanceDirectionVector[i]=1/ComputeNormOfVector3(directionVector[i]);
	
	}
    float distanceSum=0;// 邻域点到所求点的距离和
	for(int i=0;i<numOfExactNeighbor;i++){
		distanceSum+=backDistanceDirectionVector[i];
	}

	for(int i=0;i<numOfExactNeighbor;i++){
		coffOfWeight[i]=backDistanceDirectionVector[i]/distanceSum;
	}
    
	for(int i=0;i<numOfExactNeighbor;i++){	
		directionCurvature=2*Vector3Vector(m_normal,directionVector[i])/distanceDirectionVector[i];
		Vector3TVector(m_normal,m_normal,matrixNorTNor);
		Matrix3MinusMatrix3(identityMatrix,matrixNorTNor,tempMatrix);
		float tempDirectionVector[3];
		for(int j=0;j<3;j++){
			tempDirectionVector[j]=directionVector[i][j];
		}
		Matrix3MultiVector(tempMatrix,tempDirectionVector,tempVector);
		float distanceTempVector;
		distanceTempVector=ComputeNormOfVector3(tempVector);
		for(int j=0;j<3;j++){
			tangentDirectionVector[j]=tempVector[j]/distanceTempVector;
		}
		Vector3TVector(tangentDirectionVector,tangentDirectionVector,matrixTanTTan);
		tempScalar=coffOfWeight[i]*directionCurvature;
		float tempMviMatrix[3][3];
		Matrix3MultiScalar(matrixTanTTan,tempScalar,tempMviMatrix);
		for(int j=0;j<3;j++){
			for(int k=0;k<3;k++){
				MviMatrix[j][k]+=tempMviMatrix[j][k];
			}
		}
   	}


	//tempScalar=0;
	//for(int i=0;i<m_numOfNeighbor;i++){
	//	tempScalar+=coffOfWeight[i];
	//}
 //   for(int i=0;i<3;i++)
	//	for(int j=0;j<3;j++)
	//		MviMatrix[i][j]/=tempScalar;
	//////////////////////////////////////////////////////////////////////////
	// 计算主曲率
	float QMmatrix[3][3],QMQmatrix[3][3],QtransposeMatrix[3][3];
	float principalCurvature[2];
	Matrix3Transpose(householderMatrix,QtransposeMatrix);
	Matrix3Multi(QtransposeMatrix,MviMatrix,QMmatrix);
	Matrix3Multi(QMmatrix,householderMatrix,QMQmatrix);
	principalCurvature[0]=3*QMQmatrix[1][1]-QMQmatrix[2][2];
	principalCurvature[1]=3*QMQmatrix[2][2]-QMQmatrix[1][1];
	//////////////////////////////////////////////////////////////////////////
	// 计算主方向
	float sine,cosine;
	VECTOR3D principalDirection[2];
	if((QMQmatrix[1][1]*QMQmatrix[1][1]+QMQmatrix[1][2]*QMQmatrix[1][2])==0){
		sine=0;
		cosine=1;
	}
	else{
		sine=QMQmatrix[1][2]/sqrt(QMQmatrix[1][1]*QMQmatrix[1][1]+QMQmatrix[1][2]*QMQmatrix[1][2]);
		cosine=-QMQmatrix[1][1]/sqrt(QMQmatrix[1][1]*QMQmatrix[1][1]+QMQmatrix[1][2]*QMQmatrix[1][2]);

	}
	for(int i=0;i<3;i++){
		principalDirection[0][i]=cosine*householderMatrix[1][i]-sine*householderMatrix[2][i];
		principalDirection[1][i]=sine*householderMatrix[1][i]+cosine*householderMatrix[2][i];
    }
	//////////////////////////////////////////////////////////////////////////
	// 计算高斯曲率和平均曲率,主曲率和主方向
	m_gaussCurvature=principalCurvature[0]*principalCurvature[1];
	m_averageCurvature=(principalCurvature[0]+principalCurvature[1])/2;

	if(abs(principalCurvature[0])>abs(principalCurvature[1])){
		m_principalCurvature[0]=principalCurvature[0];
		m_principalCurvature[1]=principalCurvature[1];
		for(int j=0;j<3;j++){
			m_principalDirection[0][j]=principalDirection[0][j];
			m_principalDirection[1][j]=principalDirection[1][j];
		}
	}
	else{
		m_principalCurvature[0]=principalCurvature[1];
		m_principalCurvature[1]=principalCurvature[0];
		for(int j=0;j<3;j++){
			m_principalDirection[0][j]=principalDirection[1][j];
			m_principalDirection[1][j]=principalDirection[0][j];
		}
	}
    delete[]  coffOfWeight;
	delete[]  directionVector;
	delete[] distanceDirectionVector;
	delete[] backDistanceDirectionVector;

	return;
}
//初始化
void GeometricEstimation::SetInput(uint8 ID)
{
}
//返回法向
void GeometricEstimation::GetNormal(VECTOR3D normal)
{
	normal[0]=m_normal[0];
	normal[1]=m_normal[1];
	normal[2]=m_normal[2];
}
//返回平均曲率
float GeometricEstimation::GetAverageCurvature()
{
	return m_averageCurvature;
}
//返回高斯曲率
float GeometricEstimation::GetGaussCurvature()
{
	NormalEstimation();
	CurvatureEstimation();
    

	return m_gaussCurvature;
}
//返回几何特征（包括上面三项）
void GeometricEstimation::GetPointGeometricCharacter(PointGeometry pointGeometry)
{
	NormalEstimation();
	CurvatureEstimation();
	m_pointGeometry.m_gaussCurvature=m_gaussCurvature;
	m_pointGeometry.m_averageCurvature=m_averageCurvature;
	for(int i=0;i<3;i++){
		m_pointGeometry.m_normal[i]=m_normal[i];
	//	m_pointGeometry.m_pointPosition[i]=m_pointSet[m_ID*3+i];
    }
	return;

}
MZSplatEllipse GeometricEstimation::GetSplatEllipse()
{
	MZSplatEllipse m_MZSplatEllipse;
	NormalEstimation();
	CurvatureEstimation();
	for(int i=0;i<3;i++){
		m_MZSplatEllipse.m_position[i]=m_pointSet[m_ID*3+i];
		m_MZSplatEllipse.m_normal[i]=m_normal[i];
		m_MZSplatEllipse.m_Udirection[i]=m_principalDirection[0][i];
		m_MZSplatEllipse.m_Vdirection[i]=m_principalDirection[1][i];

	}
    //求取邻域的最小值
	m_MZSplatEllipse.m_size[1]=0;
	for(int i=0;i<3;i++){
		m_MZSplatEllipse.m_size[1]+=(m_pointSet[m_ID*3+i]-m_pointSet[m_pointNeighbor[0]*3+i])*(m_pointSet[m_ID*3+i]-m_pointSet[m_pointNeighbor[0]*3+i]);
	}
	m_MZSplatEllipse.m_size[1]=1;//sqrt(m_MZSplatEllipse.m_size[1]);
	m_MZSplatEllipse.m_size[0]=1;//m_MZSplatEllipse.m_size[1]*fabs(m_principalCurvature[0]/m_principalCurvature[1]);
	return m_MZSplatEllipse;        
}


void GeometricEstimation:: EstimateNormalDirection(int numOfPoints,float* pointSet,uint8 numOfEachPointNeighbor,int* allPointsNeighbor,float* allPointsNormal)
{
	//首先确定z值最大的那一点的法向
	int maxZpoint=0;
	for(int i=0;i<numOfPoints;i++){
		if(pointSet[i*3+2]>pointSet[maxZpoint*3+2])
			maxZpoint=i;

	}
	if(allPointsNormal[maxZpoint*3+2]<0){
		for(int i=0;i<3;i++){
			allPointsNormal[maxZpoint*3+i]=-allPointsNormal[maxZpoint*3+i];
		}
	}

	int* parent=new int[numOfPoints];//存放父节点

	std::multimap<float,int> valueKeySeries;//根据点的代价索引的map
	std::map<int,float> pointKeySeries;//根据点索引的map

	//初始化map
	for(int i=0;i<maxZpoint;i++){
		pointKeySeries.insert(std::map<int,float>::value_type(i,i+10));
		valueKeySeries.insert(std::map<float,int>::value_type(i+10,i));
	}
	for(int i=maxZpoint+1;i<numOfPoints;i++){
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
		
		for(uint i=0;i<numOfEachPointNeighbor;i++){
			int neighborPoint=allPointsNeighbor[pointExtract*numOfEachPointNeighbor+i];
			mapIterator=pointKeySeries.find(neighborPoint);
			if(mapIterator!=pointKeySeries.end()){
                multimapIterator=valueKeySeries.find((*mapIterator).second);
				if(neighborPoint!=(*multimapIterator).second){
					continue;
				}
				float relativeCost,normalDotProduct;
				float pointExtractNormal[3],neighborPointNormal[3];
				for(int j=0;j<3;j++){
					pointExtractNormal[j]=allPointsNormal[pointExtract*3+j];
					neighborPointNormal[j]=allPointsNormal[neighborPoint*3+j];
				}
				normalDotProduct=Vector3Vector(pointExtractNormal,neighborPointNormal);
				if(normalDotProduct<0){
					normalDotProduct=-normalDotProduct;
					for(int j=0;j<3;j++){
						allPointsNormal[neighborPoint*3+j]=-allPointsNormal[neighborPoint*3+j];
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

	if(mapIterator!=pointKeySeries.end()){
		int numOfPointsLeave;
		float* pointSetLeave;
		numOfPointsLeave=pointKeySeries.size();
		pointSetLeave=new float[numOfPointsLeave*3];
		int i=0;
		while(mapIterator!=pointKeySeries.end()){
			int tempPointIndex;
			tempPointIndex=(*mapIterator).first;
			for(int j=0;j<3;j++){
				pointSetLeave[i*3+j]=m_pointSet[tempPointIndex*3+j];
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
	}
	pointKeySeries.clear();
	valueKeySeries.clear();
	if(parent!=NULL)
        delete[] parent;
	return;
}

//得到所有点的法向
void GeometricEstimation::GetAllPointsNormalFirst()
{
	if(m_allPointsNormal!=NULL)
		delete m_allPointsNormal;
	m_allPointsNormal=new float[m_numOfPoint*3];
	for(int i=0;i<m_numOfPoint;i++){
		m_ID=i;
		for(int j=0;j<3;j++){
			m_position[j]=m_pointSet[i*3+j];
		}
		m_pointNeighbor=new int[m_numOfNeighbor];
		for(int j=0;j<m_numOfNeighbor;j++){
			m_pointNeighbor[j]=m_allPointsNeighbor[m_numOfNeighbor*i+j];
		}

		NormalEstimation();
		delete[] m_pointNeighbor;
		for(int j=0;j<3;j++){
			m_allPointsNormal[i*3+j]=m_normal[j];
		}
	}
	EstimateNormalDirection(m_numOfPoint,m_pointSet,m_numOfNeighbor,m_allPointsNeighbor,m_allPointsNormal);
	return;
}
void GeometricEstimation::GetAllPointsNeighbor(uint8 knearest,float radius)
{
	m_allPointsNeighbor=new int[m_numOfPoint*knearest];
	for(int i=0;i<m_numOfPoint;i++){
		int* neighborOfPoint=new int[knearest];
		m_ID=i;
		SearchNeighbor(knearest,radius,neighborOfPoint);
		for(int j=0;j<knearest;j++){
			m_allPointsNeighbor[i*knearest+j]=neighborOfPoint[j];
		}
		delete[] neighborOfPoint;
	}
	return;
}

void GeometricEstimation::GetAllPointsGeometricChracter()
{
	m_allPointsPricipalCurvatureOne=new float[m_numOfPoint];
	m_allPointsPricipalCurvatureTwo=new float[m_numOfPoint];
	m_allPointsPricipalDirectionTwo=new float[m_numOfPoint*3];
	m_allPointsPricipalDirectionOne=new float[m_numOfPoint*3];
	for(int i=0;i<m_numOfPoint;i++){
		m_ID=i;
		m_pointNeighbor=new int[m_numOfNeighbor];
		for(int j=0;j<m_numOfNeighbor;j++){
			m_pointNeighbor[j]=m_allPointsNeighbor[i*m_numOfNeighbor+j];
		}
		for(j=0;j<3;j++){
			m_normal[j]=m_allPointsNormal[i*3+j];
		}
		CurvatureEstimation();
		m_allPointsPricipalCurvatureOne[i]=m_principalCurvature[0];
		m_allPointsPricipalCurvatureTwo[i]=m_principalCurvature[1];
		for(int j=0;j<3;j++){
			m_allPointsPricipalDirectionOne[i*3+j]=m_principalDirection[0][j];
			m_allPointsPricipalDirectionTwo[i*3+j]=m_principalDirection[1][j];
		}
		delete[] m_pointNeighbor;		
	}
    
}

void GeometricEstimation::GetAllPointsGaussAverageCurvature()
{
	m_allPointsAverageCurvature=new float[m_numOfPoint];
	m_allPointsGaussCurvature=new float[m_numOfPoint];
	for(int i=0;i<m_numOfPoint;i++){
		m_ID=i;
		m_pointNeighbor=new int[m_numOfNeighbor];
		for(int j=0;j<m_numOfNeighbor;j++){
			m_pointNeighbor[j]=m_allPointsNeighbor[i*m_numOfNeighbor+j];
		}
		for(j=0;j<3;j++){
			m_normal[j]=m_allPointsNormal[i*3+j];
		}
		CurvatureEstimation();
	    m_allPointsAverageCurvature[i]=m_averageCurvature;
		m_allPointsGaussCurvature[i]=m_gaussCurvature;		
		delete[] m_pointNeighbor;		
	}

}

void GeometricEstimation::GetAllPointsNormalSecond()
{
	//if(m_allPointsNormal!=NULL)
	//	delete m_allPointsNormal;
	//m_allPointsNormal=new float[m_numOfPoint*3];
	for(int i=0;i<m_numOfPoint;i++){
		m_ID=i;
		for(int j=0;j<3;j++){
			m_position[j]=m_pointSet[i*3+j];
		}
		m_pointNeighbor=new int[m_numOfNeighbor];
		for(int j=0;j<m_numOfNeighbor;j++){
			m_pointNeighbor[j]=m_allPointsNeighbor[m_numOfNeighbor*i+j];
		}

		NormalEstimation();
		delete[] m_pointNeighbor;

		if((m_allPointsNormal[i*3]*m_normal[0]+m_allPointsNormal[i*3+1]*m_normal[1]+m_allPointsNormal[i*3+2]*m_normal[2])>0){
			for(int j=0;j<3;j++){
				m_allPointsNormal[i*3+j]=m_normal[j];
			}
		}
		else{
			for(int j=0;j<3;j++){
				m_allPointsNormal[i*3+j]=-m_normal[j];
			}
		}		
	}	
	return;
}