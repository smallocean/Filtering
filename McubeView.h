// McubeView.h : interface of the CMcubeView class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_MCUBEVIEW_H__7581149F_D5AB_4CFE_8E75_06975F065593__INCLUDED_)
#define AFX_MCUBEVIEW_H__7581149F_D5AB_4CFE_8E75_06975F065593__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "McubeDoc.h"
#include "GL/gl.h"
#include "GL/glu.h"
//#include "GL/glaux.h"


class CMcubeView : public CView
{
protected: // create from serialization only
	CMcubeView();
	DECLARE_DYNCREATE(CMcubeView)

// Attributes
public:
    
	CMcubeDoc* GetDocument();

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CMcubeView)
	public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	protected:
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);
	//}}AFX_VIRTUAL

// Implementation
public:
	GLfloat m_ClearColorBlue;
	GLfloat m_ClearColorRed;
	GLfloat m_ClearColorGreen;
	CPoint m_LeftDownPos;
	GLint m_GLPixelIndex;
	HGLRC m_hGLContext;
	GLbyte m_LeftButtonDown;
	GLbyte m_RightButtonDown;
	GLbyte m_CursorRotation;
	GLfloat m_xRotation;
    GLfloat m_yRotation;
	GLfloat m_xTraslation;
    GLfloat m_yTraslation;
	GLfloat m_zTraslation;
	GLfloat m_xScaling;	
	GLfloat m_yScaling;	
	
	GLfloat m_zScaling;
	virtual ~CMcubeView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	void makeTexture(void);
	void myinit(void);
	BOOL CreateViewGLContext(HDC hDC);
	BOOL SetWindowPixelFormat(HDC hDC);
	void InitGeometry(void);
	//{{AFX_MSG(CMcubeView)
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnPaint();
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnTimer(UINT nIDEvent);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in McubeView.cpp
inline CMcubeDoc* CMcubeView::GetDocument()
   { return (CMcubeDoc*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_MCUBEVIEW_H__7581149F_D5AB_4CFE_8E75_06975F065593__INCLUDED_)
