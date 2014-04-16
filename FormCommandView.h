#if !defined(AFX_FORMCOMMANDVIEW_H__270915B2_62ED_4285_8544_B948F25C0510__INCLUDED_)
#define AFX_FORMCOMMANDVIEW_H__270915B2_62ED_4285_8544_B948F25C0510__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// FormCommandView.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CFormCommandView form view

#ifndef __AFXEXT_H__
#include <afxext.h>
#endif
#include "afxwin.h"

class CFormCommandView : public CFormView
{
protected:
	CFormCommandView();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CFormCommandView)

// Form Data
public:
	//{{AFX_DATA(CFormCommandView)
	enum { IDD = IDD_FORMVIEW };
	CSliderCtrl	m_SliderScaleY;
	CSliderCtrl	m_SliderScaleZ;
	CSliderCtrl	m_SliderScaleX;
	CStatic	m_ControlColorLightSpecular;
	CStatic	m_ControlColorLightDiffuse;
	CStatic	m_ControlColorLightAmbient;
	CStatic	m_ControlBackColor;
	BOOL	m_Antialias;
	BOOL	m_Smooth;
	BOOL	m_Lighting;
	BOOL	m_VRotate;
	BOOL	m_LinkScale;
	int		m_Model;
	//}}AFX_DATA

// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CFormCommandView)
	public:
	virtual void OnInitialUpdate();
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	CView* GetMcubeView();
	BOOL UpdateScale();
	unsigned long BackColor;
	unsigned long AmbientColor;
	unsigned long SpecularColor;
	unsigned long DiffuseColor;
	virtual ~CFormCommandView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
	//{{AFX_MSG(CFormCommandView)
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnCheckAntialias();
	afx_msg void OnCheckSmooth();
	afx_msg void OnCheckLight();
	afx_msg void OnRadioModel10();
	afx_msg void OnRadioModel11();
	afx_msg void OnRadioModel12();
	afx_msg void OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	afx_msg void OnCheckLinkScale();
	afx_msg void OnPaint();
	afx_msg void OnCheckVRotate();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnNMCustomdrawSlider1(NMHDR *pNMHDR, LRESULT *pResult);

	int m_Isovalue;

		 

	float m_variationF;
	float m_variationG1;
	float m_variationG2;
	float m_variationH1;
	float m_variationH2;
	float m_variationH3;
	int m_K;
	float m_variationNstop;
	int m_iterativeTimes;
	int m_indexFunction;
	int m_meanShift;
	afx_msg void OnBnClickedButtonCoffOk();
	float m_noiseCoff;
	float m_funcitonVariation;
	int m_visualization;
	float m_noiseFuncCoff;
	float m_thresholdColor;
	float m_functionGradientWide;
	float m_functionGradientThreshold;
	float m_thresholdDistanceNormal;
	float m_thresholdDistanceGradientFunction;

	float m_timeStep;
	float m_loadConstant;
	float m_positionError;
	float m_radius;
	int m_stopN;
	int m_maxIter;
	float m_curveThreshold;
protected:
	int m_tangentOrManifold;
public:
	int m_ifNormalWeight;
	int m_ifAreaWeight;
protected:
	int m_ifVolumePreserve;
public:
	int m_ifVariationNormal;
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_FORMCOMMANDVIEW_H__270915B2_62ED_4285_8544_B948F25C0510__INCLUDED_)
