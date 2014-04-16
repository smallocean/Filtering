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
	static void ChasingSolvingL(double* a, double* b, unsigned int n, double* f);//׷�Ϸ�����󷽳�
	static void ChasingSolvingU(double* c, unsigned int n, double* f);
public:
	static bool ChasingAlgorithm(double* a, double* b, double* c, unsigned int n, double* f);/////////////////////

public:
	static double GJMatrixInvert(double* a,unsigned int n);//��˹������ȥ������������
private:
	static void ExchangeRow(double* a,unsigned int n,unsigned int rowA,unsigned int rowB);//�������������
	static void ExchangeColume(double* a,unsigned int n,unsigned int columeA,unsigned int columeB);//�������������
	static double GJMatrixStep(double* a,unsigned int n,unsigned int step,unsigned int &Ipk);//��˹������ȥ����ȥ��step�У�����������ʽֵ��һ������
	static int ChooseMax(double* a,unsigned int n,unsigned int colume);//ѡȡcolume���д�colume�п�ʼ�ľ���ֵ�����Ԫ������ֵΪ����Ԫ�����е�λ��

public:
	static bool GaussMainDec(double* A,unsigned int n,double* b);//��˹ѡ��Ԫ�����Ƿֽⷨ
private:
	static bool GaussMainDecStep(double* A,unsigned int n,unsigned int step,unsigned int &Ipr);//��˹��ȥ����һ�����Ĳ�����ֵ������P181ҳ��

public:
	static void MatrixAdd(double* a,double* b,unsigned int m,unsigned int n,double* c);//�������a+b=c�����������Ϊm��n����

public:
	static void MatrixFuncSOR(double* A,double* b,double* x,int n,unsigned int &k,double w,double eps);//��γ��ɳڵ�����������A������ֵb��xΪ���̵Ľ⣬kΪ����������wΪ�ɳ����ӣ�epsΪ����Ҫ��

public:
	static void JacobiIteration(double* A,double* b,double* x,int n,unsigned int &k,double eps);//�ſɱȵ�����������A������ֵb��xΪ���̵Ľ⣬kΪ����������epsΪ����Ҫ��

public:
	static double Romberg(double(*fun)(double),double a,double b,int &k,double eps);//�������㷨��funΪ����ָ�룬a��bΪ���䣬kΪ���ִ�����epsΪ����Ҫ��
	static double functionA(double x){return 4.0/(pow(x,4)+1);};
};
