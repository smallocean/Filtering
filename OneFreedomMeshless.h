#pragma once
#include "NormalStruct.h"

#if !defined (stiff_Matrix)
#define stiff_Matrix
	struct StiffMatrix{
	std::vector<double> m_elementStiffMarix;//存放刚度矩阵的元素
};// 保存点刚度矩阵的结构
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
	double* m_originalNormals; //原始法向
	double* m_filterNormals; //滤波后的法向
	double* m_originalMajorDirection;//原始主方向（对应于最大的特征值）
	double* m_originalMinorDirection;//原始次方向（对应于第二大特征值）
	double* m_meanCurvature;//平均曲率
	double* m_gaussCurvature; //高斯曲率
	double m_radius;//计算 m_mapKnearest时候的搜索半径
	int m_interactiveTime; // 滤波的次数
	double* m_loadVector;//荷载向量
	StiffMatrix* m_stiffMatrix;// 刚度矩阵
	StiffMatrix* m_massMatrix;// 质量矩阵
	double m_namda;// 时间步长
	double m_loadConstant;//荷载常数
	double m_ebusainu;//循环停止的条件；
protected://与插值函数有关的


protected:  //与邻域有关的参数
	int m_numOfNeighbors;          //邻域数不包括中心点
	double m_localIntergralRadius;//      局部积分半径
	double m_localIntergralRadiusSquare;//局部积分半径的平方
	double m_localWeightRadius;//         局部加权半径
	double m_localWeightRadiusC;//局部加权半径的4倍数；
	double m_localWeightRadiusCSquare;//局部加权半径的4倍数的平方
	double m_localIntergralRadiusC;//局部积分半径的4倍数；
	double m_localIntergralRadiusCSquare;//局部积分半径的4倍数的平方
	double m_localWeightRadiusSquare;//   局部加权半径的平方
	double* m_localOrdinate;//             各邻域点的局部坐标
	double* m_localIntegralPoints;  //积分点的坐标
	double* m_localQmatrix;//              插值函数矩阵
	double* m_localQTmatrix;//             插值函数的转置矩阵
protected:  //与采样有关的参数
	double* m_kernelShapeFunction;
	double* m_kernelDerivativeOneShapeFunctionOne;
	double* m_kernelDerivativeTwoShapeFunctionOne;
	double* m_kernelDerivativeOneShapeFunctionTwo;
	double* m_kernelDerivativeTwoShapeFunctionTwo;
	double* m_testFunction;
	double* m_testFunctionDerivativeOne;
	double* m_testFunctionDerivativeTwo;
protected:  //形函数及其导数
	double* m_shapeFunction;//                各邻域点的形函数在各采样点的值
	double* m_shapeFunctionDerivativeOne;//   各邻域点的形函数对变量一的倒数在各采样点的值
	double* m_shapeFunctionDerivativeTwo;//   各邻域点的形函数对变量二的倒数在各采样点的值
	double* m_localOrdinateVectorOne;//空间坐标对局部坐标的导数；
	double* m_gradientShapeFunction;  // 形函数的梯度
	double* m_gradientTestFunction; //    测试函数的梯度
	double* m_localOrdinateVectorTwo;
protected:
	double* m_determinentMetric;//   在采样点处流形标准的代数行列式
	double* m_metricMatrix;//在采样点处流形标准的矩阵
	double* m_metricMatrixInve;////在采样点处流形标准的矩阵的逆

public:
	double* m_resultPointSet;
	double* m_resultColor;
	long nearestTime;
	long filterTime;
	std::map<int, KnearestField> m_mapKnearest;
public:
	void GetMeshlessFilter(int numOfPointSet,float* pointSet,int kNearest, double radius, int ifFilterNormal,int interactiveTime,double meanShiftHposition,double meanShiftHnormal,	double namda,	double loadConstant,	double ebusainu);
	/*
	函数功能：  接口函数，计算点集滤波
	变量说明：
	int numOfPointSet                        点的数量
	double* pointSet							原始点集
	int kNearest                             k邻域的数量
	double radius                             k邻域搜索的半径
	int ifFilterNormal                       是否对法向进行滤波
	int interactiveTime                      滤波的交互次数
	double meanShiftHposition                法向滤波时的位置方差
	double meanShiftHnormal                  法向滤波时的法向方差
	*/

	void DeleteMeshlessFilter();
	/*
	*	函数功能： 释放变量内存
	*/
protected:

	void ComputeMapKnearest();
	/*
	*	函数功能：计算每一个点的 k邻域
	*/

	void ComputeNormal();
	/*
	*	函数功能：计算每一个点的法向
	*/

	void ComputeNearestParameter(int indexPointSet);
	/*
	*	函数功能： 计算任一点处的积分半径，计算其邻域各点的局部坐标， 计算局部Q矩阵及其转置QT
	*  变量说明：
	*  int indexPointSet        索引点
	*  double m_localIntergralRadius       局部积分半径
	*  double* m_localOrdinate             各邻域点的局部坐标
	*  double* m_localQmatrix              插值函数矩阵
	*  double* m_localQTmatrix             插值函数的转置矩阵
	*/
	void ComputeQmatrix(int indexPointSet);
	/*
	*	函数功能： 计算矩阵Q及其转置
	*/

	void ComputeMatrixRelatedSampling(int indexPointSet);
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

	void ComputeShapeFunction(int indexPointSet);
	/*
	*	函数功能：计算形函数及其导数
	*  变量说明：
	*  int indexPointSet        索引点
	*  double* m_shapeFunction                各邻域点的形函数在各采样点的值
	*  double* m_shapeFunctionDerivativeOne   各邻域点的形函数对变量一的倒数在各采样点的值
	*  double* m_shapeFunctionDerivativeTwo   各邻域点的形函数对变量二的倒数在各采样点的值
	*/

	void ComputeManifoldMetric(int indexPointSet);
	/*
	*	函数功能： 计算流形标准的代数行列式
	*  变量说明：
	*  int indexPointSet        索引点
	*  double* m_determinentMetric   在采样点处流形标准的代数行列式
	*/

	void ComputeStiffMatrix(int indexPointSet);
	/*
	*	函数功能： 计算刚度矩阵和荷载向量
	*  变量说明： 
	*  int indexPointSet                    索引点
	*  double* m_loadVector                  荷载向量
	*  StiffMatrix* m_stiffMatrix           刚度矩阵
	*/

	void EstimateNormalDirection();
	/*
	*	函数功能： 估计法向的准确方向
	*/
	void CalculateLinearSystem();
	/*
	*	函数功能： 解线性方程
	*/

	void  ComputeGradientShapeFunction();
	/*
	 *	函数功能： 计算每一个积分点形函数的梯度和测试函数的梯度
	 *  变量说明：
	 *  double* m_gradientShapeFunction   形函数的梯度
	 *  double* m_gradientTestFunction    测试函数的梯度
	 */

	void CalculateLinearSystemGauss();
	/*
	 *	函数功能： 高斯消元法解方程；
	 */
	void ComputeTestFunction(int indexPointSet);
	/*
	 *	计算测试函数在积分点的函数值及其导数值
	 */
	void AssembleStiffMatrix();
	//////////////////////////////////////////////////////////////////////////	
	/*
	*	函数功能： 组装刚度矩阵
	*/
	//////////////////////////////////////////////////////////////////////////

};

