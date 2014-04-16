#pragma once
#include "Math.h"

class NumericalAnalysis
{
public:
	NumericalAnalysis(void);
	~NumericalAnalysis(void);

private:
	static bool IsChasingMatrix(double* a, double* b, double* c,unsigned int n);//////////////////////////////////
	static void ChasingMatrixDec(double* a, double* b, double* c, unsigned int n);
	static void ChasingSolvingL(double* a, double* b, unsigned int n, double* f);//追赶法解矩阵方程
	static void ChasingSolvingU(double* c, unsigned int n, double* f);
public:
	static bool ChasingAlgorithm(double* a, double* b, double* c, unsigned int n, double* f);/////////////////////

public:
	static double GJMatrixInvert(double* a,unsigned int n);//高斯若当消去法求矩阵的逆阵
private:
	static void ExchangeRow(double* a,unsigned int n,unsigned int rowA,unsigned int rowB);//交换矩阵的两行
	static void ExchangeColume(double* a,unsigned int n,unsigned int columeA,unsigned int columeB);//交换矩阵的两列
	static double GJMatrixStep(double* a,unsigned int n,unsigned int step,unsigned int &Ipk);//高斯若当消去法消去第step列，并返回行列式值的一个因子
	static int ChooseMax(double* a,unsigned int n,unsigned int colume);//选取colume列中从colume行开始的绝对值最大主元，返回值为该主元在列中的位置

public:
	static bool GaussMainDec(double* A,unsigned int n,double* b);//高斯选主元的三角分解法
private:
	static bool GaussMainDecStep(double* A,unsigned int n,unsigned int step,unsigned int &Ipr);//高斯消去法第一到第四步（数值分析书P181页）

public:
	static void MatrixAdd(double* a,double* b,unsigned int m,unsigned int n,double* c);//矩阵求和a+b=c，三个矩阵均为m×n矩阵

public:
	static void MatrixFuncSOR(double* A,double* b,double* x,int n,unsigned int &k,double w,double eps);//逐次超松弛迭代法，矩阵A、方程值b，x为方程的解，k为迭代次数，w为松弛因子，eps为精度要求

public:
	static void JacobiIteration(double* A,double* b,double* x,int n,unsigned int &k,double eps);//雅可比迭代法，矩阵A、方程值b，x为方程的解，k为迭代次数，eps为精度要求

public:
	static double Romberg(double(*fun)(double),double a,double b,int &k,double eps);//龙贝格算法，fun为函数指针，a，b为区间，k为二分次数，eps为精度要求
	static double functionA(double x){return 4.0/(pow(x,4)+1);};
};
