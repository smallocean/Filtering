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

bool NumericalAnalysis::IsChasingMatrix(double* a, double* b, double* c, unsigned int n)//�жϾ����Ƿ�����׷�Ϸ�Ҫ��
																//aΪ�Խ�������Ԫ��(n-1)��bΪ�Խ���Ԫ��(n)��cΪ�Խ����Ϸ�Ԫ��(n-1)
{
	for(unsigned int i=0;i<n;i++)//�Խ���Ԫ��bn��Ϊ��
	{
		if(b[i]==0)
			return false;
	}
	for(unsigned int i=0;i<n-1;i++)//�Խ�������Ԫ��an��cn��Ϊ��
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

void NumericalAnalysis::ChasingMatrixDec(double* a, double* b, double* c, unsigned int n)//�����LU�ֽ�
										//�ֽ��alpha����ԭ����b�У�beta����ԭ����c�У�gama���伴ԭ����a
										//aΪ�Խ�������Ԫ��(n-1)��bΪ�Խ���Ԫ��(n)��cΪ�Խ����Ϸ�Ԫ��(n-1)
										//�����ֽ����ԭ�����ڴ�ռ��ŷֽ���alpha��beta��gama
{
	c[0]=c[0]/b[0];
	for(unsigned int i=0;i<n-2;i++)
	{
		b[i+1]=b[i+1]-a[i]*c[i];
		c[i+1]=c[i+1]/b[i+1];
	}
	b[n-1]=b[n-1]-a[n-2]*c[n-2];
}

void NumericalAnalysis::ChasingSolvingL(double* a, double* b, unsigned int n, double* f)//��Ly=f����
										//aΪ�Խ�������Ԫ��(n-1)��bΪ�Խ���Ԫ��(n)��fΪ���̵�ֵ(n)
{
	f[0]=f[0]/b[0];
	for(unsigned int i=0;i<n-1;i++)
	{
		f[i+1]=(f[i+1]-a[i]*f[i])/b[i+1];
	}
}

void NumericalAnalysis::ChasingSolvingU(double* c, unsigned int n, double* f)//��Ux=y���̣�cΪ�Խ����Ϸ�Ԫ��(n-1)��fΪ���̵�ֵ(n)
{
	for(int i=n-2;i>=0;i--)
	{
		f[i]=f[i]-c[i]*f[i+1];
	}
}

bool NumericalAnalysis::ChasingAlgorithm(double* a, double* b, double* c, unsigned int n, double* f)//׷�Ϸ�
										//aΪ�Խ�������Ԫ��(n-1)��bΪ�Խ���Ԫ��(n)��cΪ�Խ����Ϸ�Ԫ��(n-1)��fΪ���̵�ֵ(n)
{
	if(!IsChasingMatrix(a,b,c,n))//�жϾ����Ƿ�����׷�Ϸ�Ҫ��
		return false;

	ChasingMatrixDec(a,b,c,n);//�����LU�ֽ�
	ChasingSolvingL(a,b,n,f);//��Ly=f
	ChasingSolvingU(c,n,f);//��Ux=f ���Ľ������f��

	return true;
}

double NumericalAnalysis::GJMatrixInvert(double* a,unsigned int n)//��˹������ȥ������������
{
	double det=1.0;
	unsigned int *Ip=new unsigned int[n];

	for(unsigned int i=0;i<n;i++)
	{
		Ip[i]=0;
		det*=GJMatrixStep(a,n,i,Ip[i]);//�ظ�n����1����8��
	}

	for(unsigned int i=0;i<n;i++)///////////////////////////////////////////////////
	{
		unsigned int t=Ip[n-i-1];//��ֵ������P178ҳ���߲�
		if(t==n-i-1)
			continue;
		else
			ExchangeColume(a,n,t,n-i-1);
	}//////////////////////////////////////////////////////////////////////

	return det;
}

void NumericalAnalysis::ExchangeRow(double* a,unsigned int n,unsigned int rowA,unsigned int rowB)//��������
{
	double swap=0.0;
	for(unsigned int i=0;i<n;i++)
	{
		swap=a[rowA*n+i];
		a[rowA*n+i]=a[rowB*n+i];
		a[rowB*n+i]=swap;
	}
}

void NumericalAnalysis::ExchangeColume(double* a,unsigned int n,unsigned int columeA,unsigned int columeB)//��������
{
	double swap=0.0;
	for(unsigned int i=0;i<n;i++)
	{
		swap=a[i*n+columeA];
		a[i*n+columeA]=a[i*n+columeB];
		a[i*n+columeB]=swap;
	}
}

double NumericalAnalysis::GJMatrixStep(double* a,unsigned int n,unsigned int step,unsigned int &Ipk)//��˹������ȥ����ȥ��step�У�����������ʽֵ��һ������
{
	double h=0.0;
	int rowMax=ChooseMax(a,n,step);//��þ���ֵ���Ԫ�����ڵ���
	if(rowMax==-1)
		return 0;//�������ֵ���Ϊ0��������ʽΪ�㣬�����������
	Ipk=rowMax;
	double Co=a[rowMax*n+step];//��ø��о���ֵ���Ԫ��
	double det=Co;

	if(rowMax!=step)
	{
		ExchangeRow(a,n,step,rowMax);//����Խ�����Ԫ�ز��Ǿ���ֵ���Ԫ�أ��򽻻�������
		det=-det;
	}

	h=1/Co;///////////////////////////////////////////////////////////////
	for(unsigned int i=0;i<n;i++)
	{
		a[i*n+step]=-a[i*n+step]*h;//��ֵ������P178ҳ������
	}
	a[step*n+step]=1/Co;//////////////////////////////////////////////////

	for(unsigned int i=0;i<n;i++)//////////////////////////////////////////////////
	{
		for(unsigned int j=0;j<n;j++)
		{
			if(i==step||j==step)//��ֵ������P178ҳ���߲�
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
			a[step*n+i]=a[step*n+i]*h;//��ֵ������P178ҳ�ڰ˲�
		}
	}/////////////////////////////////////////////////////////////////////

	return det;
}

int NumericalAnalysis::ChooseMax(double* a,unsigned int n,unsigned int colume)//ѡȡcolume���д�colume�п�ʼ�ľ���ֵ�����Ԫ������ֵΪ����Ԫ�����е�λ��
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

	return rowMax;//�����������Ԫ��Ϊ0������-1����ʾ����ʽΪ0
}

void NumericalAnalysis::MatrixAdd(double* a,double* b,unsigned int m,unsigned int n,double* c)//�������a+b=c�����������Ϊm��n����
{
	for(unsigned int i=0;i<m;i++)
		for(unsigned int j=0;j<n;j++)
		{
			c[i*m+j]=a[i*m+j]+b[i*m+j];
		}
}

void NumericalAnalysis::MatrixFuncSOR(double* A,double* b,double* x,int n,unsigned int &k,double w,double eps)//��γ��ɳڵ�����������A������ֵb��xΪ���̵Ľ⣬kΪ����������wΪ�ɳ����ӣ�epsΪ����Ҫ��
{
	k=0;
	double Po=0.0;
	for(int i=0;i<n;i++)
		x[i]=0.0;

	while(true)
	{
		k++;
		Po=0.0;
		for(int i=0;i<n;i++)//��ֵ����P216ҳ���岽
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
			detaXi=w*(b[i]-sum)/A[i*n+i];//��ֵ����P216ҳ���岽1ʽ
			if(abs(detaXi)>abs(Po))
				Po=detaXi;//2ʽ
			x[i]+=detaXi;//3ʽ
		}
		if(abs(Po)<abs(eps))//��detaX���ֵС�ھ���Ҫ���������ֹ
			break;
	}
}

bool NumericalAnalysis::GaussMainDec(double* A,unsigned int n,double* b)//��˹ѡ��Ԫ�����Ƿֽⷨ
{
	unsigned int *Ip=new unsigned int[n];

	for(unsigned int r=0;r<n;r++)
	{
		bool bIsOdd=GaussMainDecStep(A,n,r,Ip[r]);
		if(bIsOdd==false)
		{
			return false;
		}
	}//���PA��LU�ֽ⣬U������A�������ǲ��֣�L�����������ǲ��֣�������P��Ip��¼

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
	}//���Pb����ֵ����P181ҳ���岽

	for(unsigned int i=1;i<n;i++)
	{
		double sum=0.0;
		for(unsigned int k=0;k<i;k++)
		{
			sum+=A[i*n+k]*b[k];
		}
		b[i]=b[i]-sum;
	}//��ֵ����P181������

	b[n-1]=b[n-1]/A[(n-1)*n+n-1];
	for(int i=n-2;i>=0;i--)
	{
		double sum=0.0;
		for(unsigned int k=i+1;k<n;k++)
		{
			sum+=A[i*n+k]*b[k];
		}
		b[i]=(b[i]-sum)/A[i*n+i];
	}//��ֵ����P181���߲�
	delete Ip;
	return true;
}

bool NumericalAnalysis::GaussMainDecStep(double* A,unsigned int n,unsigned int step,unsigned int &Ipr)//��˹��ȥ����һ�����Ĳ�����ֵ������P181ҳ��
{
	for(unsigned int i=step;i<n;i++)/////////////////////////////////////////////////////////
	{
		double sum=0.0;
		for(unsigned int k=0;k<step;k++)
		{
			sum+=A[i*n+k]*A[k*n+step];//����ֵ����P181ҳ��һ���������step��ʱ�����и�Ԫ��ֵ
		}
		A[i*n+step]=A[i*n+step]-sum;
	}////////////////////////////////////////////////////////////////////////////////////////

	unsigned int nMax=ChooseMax(A,n,step);
	if(nMax==-1)
		return false;
	Ipr=unsigned int(nMax);//ѡȡ�����о���ֵ���Ԫ�أ���¼���к�

	ExchangeRow(A,n,Ipr,step);//����Ipr��step��Ӧ������Ԫ��

	for(unsigned int i=step+1;i<n;i++)///////////////////////////////////////////////////////
	{
		A[i*n+step]=A[i*n+step]/A[step*n+step];
		double sum=0.0;
		for(unsigned int k=0;k<step;k++)//����A�ĵ�step��Ԫ�غ͵�step��Ԫ��
		{
			sum+=A[step*n+k]*A[k*n+i];
		}
		A[step*n+i]=A[step*n+i]-sum;
	}////////////////////////////////////////////////////////////////////////////////////////
	return true;
}

void NumericalAnalysis::JacobiIteration(double* A,double* b,double* x,int n,unsigned int &k,double eps)//�ſɱȵ�����������A������ֵb��xΪ���̵Ľ⣬kΪ����������epsΪ����Ҫ��
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
					continue;//����ֵ����P205ҳʽ2.3���е���
				else
				{
					sum+=A[i*n+j]*xAssist[j];
				}
			}
			x[i]=(b[i]-sum)/A[i*n+i];

			if(abs(x[i]-xAssist[i])>abs(Po))//Ϊ���㾫�ȣ��������detaX5
				Po=x[i]-xAssist[i];
		}//////////////////////////////////////////////////////////////////////

		if(abs(Po)<abs(eps))//��detaX���ֵС�ھ���Ҫ���������ֹ
			break;
		else memcpy(xAssist,x,n*sizeof(double));//�ѽ����Ԫ�е�ֵ����������Ԫ������һ������
	}
}

double NumericalAnalysis::Romberg(double(*fun)(double),double a,double b,int &k,double eps)//�������㷨��funΪ����ָ�룬a��bΪ���䣬kΪ���ִ�����epsΪ����Ҫ��
{
	list<double> TNumTable;//���ڴ洢T����Ķ���
	list<double>::iterator TNumIter;//���е�ָ��
	double TActive=(b-a)*(fun(a)+fun(b))/2.0;//����T��������
	TNumTable.push_back(TActive);//��T�����������T������
	double h=b-a;
	k=1;
	while(true)
	{
		for(int i=0;i<k+1;i++)
		{
			if(i==0)
			{
				TNumIter=TNumTable.end();
				for(int j=0;j<k;j++)//��ָ��ָ��T������һ�е�����
					TNumIter--;

				double sum=0.0;
				for(int j=0;j<pow(2,k-1);j++)
				{
					sum+=fun((a+b)/pow(2.0,k)+(a+b)/pow(2.0,k-1)*double(j));//��ֵ����P91��ʽ3.1�����
				}
				TActive=(*TNumIter)/2.0+h/pow(2.0,k)*sum;//��ֵ����P91��ʽ3.1
				TNumTable.push_back(TActive);
			}
			else
			{
				TNumIter=TNumTable.end();
				for(int j=0;j<k+1;j++)//��ָ��ָ��T������һ�еĸ���ǰ��һ��
					TNumIter--;
				double TUp=*TNumIter;

				TNumIter=TNumTable.end();
				TNumIter--;//��ָ��ָ��T������һ�еĸ���ǰ��һ��
				double TDown=*TNumIter;

				TActive=(pow(4.0,i)/(pow(4.0,i)-1.0)*(TDown))-(1.0/(pow(4.0,i)-1.0)*TUp);//��ֵ����P94ҳ��ʽ3.11
				TNumTable.push_back(TActive);
			}
		}
		TNumIter=TNumTable.end();
		TNumIter--;
		double TLast=*TNumIter;//�������һ��
		for(int i=0;i<k+1;i++)
			TNumIter--;//ָ����һ�����һ��
		if(abs(TLast-(*TNumIter))<eps)//�������һ������һ�����һ��֮��������㾫��Ҫ���򷵻ظ������һ��
			return TLast;
		k++;//����k+1������������һ��
	}
}
