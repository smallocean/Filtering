#if !defined(AFX_DIALOGINTEPOLATE_H__1B140583_E22C_4782_98B4_3CB17AC68AE5__INCLUDED_)
#define AFX_DIALOGINTEPOLATE_H__1B140583_E22C_4782_98B4_3CB17AC68AE5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// DialogIntepolate.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialogIntepolate dialog

class CDialogIntepolate : public CDialog
{
// Construction
public:
	CDialogIntepolate(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CDialogIntepolate)

	BOOL	m_X;
	BOOL	m_Y;
	BOOL	m_Z;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialogIntepolate)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialogIntepolate)
	virtual void OnOK();
	virtual void OnCancel();
	afx_msg void OnCheck1();
	afx_msg void OnCheck2();
	afx_msg void OnCheck3();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DIALOGINTEPOLATE_H__1B140583_E22C_4782_98B4_3CB17AC68AE5__INCLUDED_)
