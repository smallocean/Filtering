// DialogIntepolate.cpp : implementation file
//

#include "stdafx.h"
#include "Mcube.h"
#include "DialogIntepolate.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDialogIntepolate dialog


//CDialogIntepolate::CDialogIntepolate(CWnd* pParent /*=NULL*/)
//	: CDialog(CDialogIntepolate::IDD, pParent)
//{
//	//{{AFX_DATA_INIT(CDialogIntepolate)
//	m_X = FALSE;
//	m_Y = FALSE;
//	m_Z = FALSE;
//	//}}AFX_DATA_INIT
//}


void CDialogIntepolate::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialogIntepolate)
	DDX_Check(pDX, IDC_CHECK1, m_X);
	DDX_Check(pDX, IDC_CHECK2, m_Y);
	DDX_Check(pDX, IDC_CHECK3, m_Z);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialogIntepolate, CDialog)
	//{{AFX_MSG_MAP(CDialogIntepolate)
	ON_BN_CLICKED(IDC_CHECK1, OnCheck1)
	ON_BN_CLICKED(IDC_CHECK2, OnCheck2)
	ON_BN_CLICKED(IDC_CHECK3, OnCheck3)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialogIntepolate message handlers

void CDialogIntepolate::OnOK() 
{
	// TODO: Add extra validation here
	
	CDialog::OnOK();
}

void CDialogIntepolate::OnCancel() 
{
	// TODO: Add extra cleanup here
	
	CDialog::OnCancel();
}

void CDialogIntepolate::OnCheck1() 
{
	// TODO: Add your control notification handler code here
	m_X=!m_X;
}

void CDialogIntepolate::OnCheck2() 
{
	// TODO: Add your control notification handler code here
	m_Y=!m_Y;
}

void CDialogIntepolate::OnCheck3() 
{
	// TODO: Add your control notification handler code here
	m_Z=!m_Z;
}
