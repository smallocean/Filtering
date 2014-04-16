#include "StdAfx.h"
#include "numericalanalysis.h"
#include "Math.h"
#include <list>
#include "stdlib.h"
#include ".\numericalanalysis.h"

using namespace std;

NumericalAnalysis::NumericalAnalysis(void)
{
}

NumericalAnalysis::~NumericalAnalysis(void)
{
}

bool NumericalAnalysis::IsChasingMatrix(double* a, double* b, double* c, unsigned int n)//判断矩阵是否满足追赶法要求
																//a为对角线下面元素(n-1)，b为对角线元素(n)，c为对角线上方元素(n-1)
{
	for(unsigned int i=0;i<n;i++)//对角线元素bn不为零
	{
		if(b[i]==0)
			return false;
	}
	for(unsigned int i=0;i<n-1;i++)//对角线上下元素an、cn不为零
	{
		if(a[i]==0||c[i]==0)
			return false;
	}
	if(!(abs(b[0])>abs(c[0])))//|b1|>|c1|>0
		return false;
	for(unsigned int i=0;i<n-2;i++)//|bi|>=|ai|+|ci| i=2,3...n-1
	{
		if(abs(b[i+1])<abs(a[i])+abs(c[i+1]))
			return false;
	}
	if(abs(b[n-1])<=abs(a[n-2]))//|bn|>|an|>0
		return false;
	return true;
}

void NumericalAnalysis::ChasingMatrixDec(double* a, double* b, double* c, unsigned int n)//矩阵的LU分解
										//分解后alpha放入原来的b中，beta放在原来的c中，gama不变即原来的a
										//a为对角线下面元素(n-1)，b为对角线元素(n)，c为对角线上方元素(n-1)
										//经过分解后，用原来的内存空间存放分解后的alpha、beta、gama
{
	c[0]=c[0]/b[0];
	for(unsigned int i=0;i<n-2;i++)
	{
		b[i+1]=b[i+1]-a[i]*c[i];
		c[i+1]=c[i+1]/b[i+1];
	}
	b[n-1]=b[n-1]-a[n-2]*c[n-2];
}

void NumericalAnalysis::ChasingSolvingL(double* a, double* b, unsigned int n, double* f)//解Ly=f方程
										//a为对角线下面元素(n-1)，b为对角线元素(n)，f为方程的值(n)
{
	f[0]=f[0]/b[0];
	for(unsigned int i=0;i<n-1;i++)
	{
		f[i+1]=(f[i+1]-a[i]*f[i])/b[i+1];
	}
}

void NumericalAnalysis::ChasingSolvingU(double* c, unsigned int n, double* f)//解Ux=y方程，c为对角线上方元素(n-1)，f为方程的值(n)
{
	for(int i=n-2;i>=0;i--)
	{
		f[i]=f[i]-c[i]*f[i+1];
	}
}

bool NumericalAnalysis::ChasingAlgorithm(double* a, double* b, double* c, unsigned int n, double* f)//追赶法
										//a为对角线下面元素(n-1)，b为对角线元素(n)，c为对角线上方元素(n-1)，f为方程的值(n)
{
	if(!IsChasingMatrix(a,b,c,n))//判断矩阵是否满足追赶法要求
		return false;

	ChasingMatrixDec(a,b,c,n);//矩阵的LU分解
	ChasingSolvingL(a,b,n,f);//解Ly=f
	ChasingSolvingU(c,n,f);//解Ux=f 最后的结果存在f中

	return true;
}

double NumericalAnalysis::GJMatrixInvert(double* a,unsigned int n)//高斯若当消去法求矩阵的逆阵
{
	double det=1.0;
	unsigned int *Ip=new unsigned int[n];

	for(unsigned int i=0;i<n;i++)
	{
		Ip[i]=0;
		det*=GJMatrixStep(a,n,i,Ip[i]);//重复n次做1——8步
	}

	for(unsigned int i=0;i<n;i++)///////////////////////////////////////////////////
	{
		unsigned int t=Ip[n-i-1];//数值分析书P178页第七步
		if(t==n-i-1)
			continue;
		else
			ExchangeColume(a,n,t,n-i-1);
	}//////////////////////////////////////////////////////////////////////

	return det;
}

void NumericalAnalysis::ExchangeRow(double* a,unsigned int n,unsigned int rowA,unsigned int rowB)//交换两行
{
	double swap=0.0;
	for(unsigned int i=0;i<n;i++)
	{
		swap=a[rowA*n+i];
		a[rowA*n+i]=a[rowB*n+i];
		a[rowB*n+i]=swap;
	}
}

void NumericalAnalysis::ExchangeColume(double* a,unsigned int n,unsigned int columeA,unsigned int columeB)//交换两列
{
	double swap=0.0;
	for(unsigned int i=0;i<n;i++)
	{
		swap=a[i*n+columeA];
		a[i*n+columeA]=a[i*n+columeB];
		a[i*n+columeB]=swap;
	}
}

double NumericalAnalysis::GJMatrixStep(double* a,unsigned int n,unsigned int step,unsigned int &Ipk)//高斯若当消去法消去第step列，并返回行列式值的一个因子
{
	double h=0.0;
	int rowMax=ChooseMax(a,n,step);//获得绝对值最大元素所在的行
	if(rowMax==-1)
		return 0;//如果绝对值最大为0，则行列式为零，不存在逆矩阵
	Ipk=rowMax;
	double Co=a[rowMax*n+step];//获得该列绝对值最大元素
	double det=Co;

	if(rowMax!=step)
	{
		ExchangeRow(a,n,step,rowMax);//如果对角线上元素不是绝对值最大元素，则交换这两行
		det=-det;
	}

	h=1/Co;///////////////////////////////////////////////////////////////
	for(unsigned int i=0;i<n;i++)
	{
		a[i*n+step]=-a[i*n+step]*h;//数值分析书P178页第六步
	}
	a[step*n+step]=1/Co;//////////////////////////////////////////////////

	for(unsigned int i=0;i<n;i++)//////////////////////////////////////////////////
	{
		for(unsigned int j=0;j<n;j++)
		{
			if(i==step||j==step)//数值分析书P178页第七步
				continue;
			else
			{
				a[i*n+j]=a[i*n+j]+a[i*n+step]*a[step*n+j];
			}
		}
	}/////////////////////////////////////////////////////////////////////

	for(unsigned int i=0;i<n;i++)//////////////////////////////////////////////////
	{
		if(i==step)
			continue;
		else
		{
			a[step*n+i]=a[step*n+i]*h;//数值分析书P178页第八步
		}
	}/////////////////////////////////////////////////////////////////////

	return det;
}

int NumericalAnalysis::ChooseMax(double* a,unsigned int n,unsigned int colume)//选取colume列中从colume行开始的绝对值最大主元，返回值为该主元在列中的位置
{
	int rowMax=-1;
	double max=0.0;

	for(unsigned int i=colume;i<n;i++)
	{
		if(max<abs(a[i*n+colume]))
		{
			max=abs(a[i*n+colume]);
			rowMax=i;
		}
	}

	return rowMax;//如果该列所有元素为0，返回-1，表示行列式为0
}

void NumericalAnalysis::MatrixAdd(double* a,double* b,unsigned int m,unsigned int n,double* c)//矩阵求和a+b=c，三个矩阵均为m×n矩阵
{
	for(unsigned int i=0;i<m;i++)
		for(unsigned int j=0;j<n;j++)
		{
			c[i*m+j]=a[i*m+j]+b[i*m+j];
		}
}

void NumericalAnalysis::MatrixFuncSOR(double* A,double* b,double* x,int n,unsigned int &k,double w,double eps)//逐次超松弛迭代法，矩阵A、方程值b，x为方程的解，k为迭代次数，w为松弛因子，eps为精度要求
{
	k=0;
	double Po=0.0;
	for(int i=0;i<n;i++)
		x[i]=0.0;

	while(true)
	{
		k++;
		Po=0.0;
		for(int i=0;i<n;i++)//数值分析P216页第五步
		{
			double detaXi=0.0;
			double sum=0.0;
			for(int j=0;j<i;j++)
			{
				sum+=A[i*n+j]*x[j];
			}
			for(int j=i;j<n;j++)
			{
				sum+=A[i*n+j]*x[j];
			}
			detaXi=w*(b[i]-sum)/A[i*n+i];//数值分析P216页第五步1式
			if(abs(detaXi)>abs(Po))
				Po=detaXi;//2式
			x[i]+=detaXi;//3式
		}
		if(abs(Po)<abs(eps))//若detaX最大值小于精度要求，则迭代终止
			break;
	}
}

bool NumericalAnalysis::GaussMainDec(double* A,unsigned int n,double* b)//高斯选主元的三角分解法
{
	unsigned int *Ip=new unsigned int[n];

	for(unsigned int r=0;r<n;r++)
	{
		bool bIsOdd=GaussMainDecStep(A,n,r,Ip[r]);
		if(bIsOdd==false)
		{
			return false;
		}
	}//完成PA的LU分解，U保存在A的上三角部分，L保存在下三角部分，排列阵P由Ip记录

	for(unsigned int i=0;i<n;i++)
	{
		if(Ip[i]==i)
			continue;
		else
		{
			double swap=b[i];
			b[i]=b[Ip[i]];
			b[Ip[i]]=swap;
		}
	}//完成Pb，数值分析P181页第五步

	for(unsigned int i=1;i<n;i++)
	{
		double sum=0.0;
		for(unsigned int k=0;k<i;k++)
		{
			sum+=A[i*n+k]*b[k];
		}
		b[i]=b[i]-sum;
	}//数值分析P181第六步

	b[n-1]=b[n-1]/A[(n-1)*n+n-1];
	for(int i=n-2;i>=0;i--)
	{
		double sum=0.0;
		for(unsigned int k=i+1;k<n;k++)
		{
			sum+=A[i*n+k]*b[k];
		}
		b[i]=(b[i]-sum)/A[i*n+i];
	}//数值分析P181第七步
	delete Ip;
	return true;
}

bool NumericalAnalysis::GaussMainDecStep(double* A,unsigned int n,unsigned int step,unsigned int &Ipr)//高斯消去法第一到第四步（数值分析书P181页）
{
	for(unsigned int i=step;i<n;i++)/////////////////////////////////////////////////////////
	{
		double sum=0.0;
		for(unsigned int k=0;k<step;k++)
		{
			sum+=A[i*n+k]*A[k*n+step];//按数值分析P181页第一步，计算第step步时，该列各元素值
		}
		A[i*n+step]=A[i*n+step]-sum;
	}////////////////////////////////////////////////////////////////////////////////////////

	unsigned int nMax=ChooseMax(A,n,step);
	if(nMax==-1)
		return false;
	Ipr=unsigned int(nMax);//选取该列中绝对值最大元素，记录其行号

	ExchangeRow(A,n,Ipr,step);//交换Ipr和step对应的两行元素

	for(unsigned int i=step+1;i<n;i++)///////////////////////////////////////////////////////
	{
		A[i*n+step]=A[i*n+step]/A[step*n+step];
		double sum=0.0;
		for(unsigned int k=0;k<step;k++)//计算A的第step行元素和第step列元素
		{
			sum+=A[step*n+k]*A[k*n+i];
		}
		A[step*n+i]=A[step*n+i]-sum;
	}////////////////////////////////////////////////////////////////////////////////////////
	return true;
}

void NumericalAnalysis::JacobiIteration(double* A,double* b,double* x,int n,unsigned int &k,double eps)//雅可比迭代法，矩阵A、方程值b，x为方程的解，k为迭代次数，eps为精度要求
{
	k=0;
	double Po=0.0;
	double *xAssist=new double[n];
	for(int i=0;i<n;i++)
	{
		x[i]=0.0;
		xAssist[i]=0.0;
	}

	while(true)
	{
		k++;
		Po=0.0;
		for(int i=0;i<n;i++)///////////////////////////////////////////////////
		{
			double sum=0.0;
			for(int j=0;j<n;j++)
			{
				if(j==i)
					continue;//按数值分析P205页式2.3进行迭代
				else
				{
					sum+=A[i*n+j]*xAssist[j];
				}
			}
			x[i]=(b[i]-sum)/A[i*n+i];

			if(abs(x[i]-xAssist[i])>abs(Po))//为计算精度，获得最大的detaX5
				Po=x[i]-xAssist[i];
		}//////////////////////////////////////////////////////////////////////

		if(abs(Po)<abs(eps))//若detaX最大值小于精度要求，则迭代终止
			break;
		else memcpy(xAssist,x,n*sizeof(double));//把结果单元中的值赋给辅助单元进行下一次运算
	}
}

double NumericalAnalysis::Romberg(double(*fun)(double),double a,double b,int &k,double eps)//龙贝格算法，fun为函数指针，a，b为区间，k为二分次数，eps为精度要求
{
	list<double> TNumTable;//用于存储T数表的队列
	list<double>::iterator TNumIter;//队列的指针
	double TActive=(b-a)*(fun(a)+fun(b))/2.0;//计算T数表首项
	TNumTable.push_back(TActive);//将T数表首项存入T数表中
	double h=b-a;
	k=1;
	while(true)
	{
		for(int i=0;i<k+1;i++)
		{
			if(i==0)
			{
				TNumIter=TNumTable.end();
				for(int j=0;j<k;j++)//将指针指向T数表上一行的首项
					TNumIter--;

				double sum=0.0;
				for(int j=0;j<pow(2,k-1);j++)
				{
					sum+=fun((a+b)/pow(2.0,k)+(a+b)/pow(2.0,k-1)*double(j));//数值分析P91公式3.1求和项
				}
				TActive=(*TNumIter)/2.0+h/pow(2.0,k)*sum;//数值分析P91公式3.1
				TNumTable.push_back(TActive);
			}
			else
			{
				TNumIter=TNumTable.end();
				for(int j=0;j<k+1;j++)//将指针指向T数表上一行的该项前面一项
					TNumIter--;
				double TUp=*TNumIter;

				TNumIter=TNumTable.end();
				TNumIter--;//将指针指向T数表这一行的该项前面一项
				double TDown=*TNumIter;

				TActive=(pow(4.0,i)/(pow(4.0,i)-1.0)*(TDown))-(1.0/(pow(4.0,i)-1.0)*TUp);//数值分析P94页公式3.11
				TNumTable.push_back(TActive);
			}
		}
		TNumIter=TNumTable.end();
		TNumIter--;
		double TLast=*TNumIter;//该行最后一项
		for(int i=0;i<k+1;i++)
			TNumIter--;//指向上一行最后一项
		if(abs(TLast-(*TNumIter))<eps)//该行最后一项与上一行最后一项之差如果满足精度要求，则返回该行最后一项
			return TLast;
		k++;//否则k+1，继续计算下一行
	}
}
