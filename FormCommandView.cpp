// FormCommandView.cpp : implementation file
//

#include "stdafx.h"
#include "Mcube.h"
#include "FormCommandView.h"
#include "McubeView.h"
#include "McubeDoc.h"
#include "Mcube.h"
#include "MainFrm.h"
#include ".\formcommandview.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CFormCommandView

IMPLEMENT_DYNCREATE(CFormCommandView, CFormView)

CFormCommandView::CFormCommandView()
	: CFormView(CFormCommandView::IDD)

	, m_variationF(0)
	, m_variationG1(0)
	, m_variationG2(0)
	, m_variationH1(0)
	, m_variationH2(0)
	, m_variationH3(0)
	, m_K(3)
	, m_variationNstop(0)
	, m_iterativeTimes(0)
	, m_indexFunction(0)
	, m_meanShift(0)
	, m_noiseCoff(0)
	, m_funcitonVariation(0)
	, m_visualization(0)
	, m_noiseFuncCoff(0)
	, m_thresholdColor(0)
	, m_functionGradientWide(0)
	, m_functionGradientThreshold(0)
	, m_thresholdDistanceNormal(0)
	, m_thresholdDistanceGradientFunction(0)

	, m_timeStep(0)
	, m_loadConstant(0)
	, m_positionError(0)
	, m_radius(0)
	, m_stopN(0)
	, m_maxIter(0)
	, m_curveThreshold(0)
	, m_tangentOrManifold(0)
	, m_ifNormalWeight(0)
	, m_ifAreaWeight(0)
	, m_ifVolumePreserve(0)
	, m_ifVariationNormal(0)
{
	//{{AFX_DATA_INIT(CFormCommandView)
	m_Antialias = FALSE;
	m_Smooth = FALSE;
	m_Lighting = FALSE;
	m_VRotate = FALSE;
	m_LinkScale = FALSE;
	m_Model = -1;
	m_Isovalue=0;
	//}}AFX_DATA_INIT
}

CFormCommandView::~CFormCommandView()
{
}

void CFormCommandView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CFormCommandView)
	DDX_Control(pDX, IDC_SLIDER2, m_SliderScaleY);
	DDX_Control(pDX, IDC_SLIDER3, m_SliderScaleZ);
	DDX_Control(pDX, IDC_SLIDER1, m_SliderScaleX);
	DDX_Control(pDX, ID_STATIC_SPECULAR, m_ControlColorLightSpecular);
	DDX_Control(pDX, IDC_STATIC_DIFFUSE, m_ControlColorLightDiffuse);
	DDX_Control(pDX, IDC_STATIC_AMBIENT, m_ControlColorLightAmbient);
	DDX_Control(pDX, IDC_STATIC_BACKGROUND, m_ControlBackColor);
	DDX_Check(pDX, IDC_CHECK1, m_Antialias);
	DDX_Check(pDX, IDC_CHECK2, m_Smooth);
	DDX_Check(pDX, IDC_CHECK3, m_Lighting);
	DDX_Check(pDX, IDC_CHECK4, m_VRotate);
	DDX_Check(pDX, IDC_CHECK5, m_LinkScale);
	DDX_Radio(pDX, IDC_RADIO1, m_Model);
	//}}AFX_DATA_MAP
	DDX_Text(pDX, IDC_EDIT_ISOVALUE, m_Isovalue);
	DDX_Text(pDX, IDC_EDIT_VF, m_variationF);
	DDX_Text(pDX, IDC_EDIT_VG1, m_variationG1);
	DDX_Text(pDX, IDC_EDIT_VG2, m_variationG2);
	DDX_Text(pDX, IDC_EDIT_VH1, m_variationH1);
	DDX_Text(pDX, IDC_EDIT_VH2, m_variationH2);
	DDX_Text(pDX, IDC_EDIT_VH3, m_variationH3);
	DDX_Text(pDX, IDC_EDIT_VK, m_K);
	DDX_Text(pDX, IDC_EDIT_VNSTOP, m_variationNstop);
	DDX_Text(pDX, IDC_EDIT_TIMES, m_iterativeTimes);
	DDX_Text(pDX, IDC_EDIT_FUNCTION, m_indexFunction);
	DDX_Text(pDX, IDC_EDIT_MSHIFT, m_meanShift);
	DDX_Text(pDX, IDC_EDIT_NOISECOFF, m_noiseCoff);
	DDX_Text(pDX, IDC_EDIT_FUNCTION_VARIATION, m_funcitonVariation);
	DDX_Text(pDX, IDC_EDIT_VISUALIZATION, m_visualization);
	DDX_Text(pDX, IDC_EDIT_NOISE_FUNCTION, m_noiseFuncCoff);
	DDX_Text(pDX, IDC_EDIT_THEOLD, m_thresholdColor);
	DDX_Text(pDX, IDC_EDIT_GRADIENT_WIDE, m_functionGradientWide);
	DDX_Text(pDX, IDC_EDIT_GRADIENT_THRESHOLD, m_functionGradientThreshold);
	DDX_Text(pDX, IDC_EDIT_Distance_Gradient, m_thresholdDistanceNormal);
	DDX_Text(pDX, IDC_EDIT_Distance_Gradient_Function, m_thresholdDistanceGradientFunction);

	DDX_Text(pDX, IDC_EDIT_Time_Step, m_timeStep);
	DDX_Text(pDX, IDC_EDIT_Load_Constant, m_loadConstant);
	DDX_Text(pDX, IDC_EDIT_Position_Error, m_positionError);
	DDX_Text(pDX, IDC_EDIT_Radius, m_radius);
	DDX_Text(pDX, IDC_EDIT_Radius2, m_stopN);
	DDX_Text(pDX, IDC_EDIT_Radius3, m_maxIter);
	DDX_Text(pDX, IDC_EDIT_Radius4, m_curveThreshold);
	DDX_Text(pDX, IDC_EDIT_Tangent_Manifold, m_tangentOrManifold);
	DDX_Text(pDX, IDC_EDIT_Normal_Weight, m_ifNormalWeight);
	DDX_Text(pDX, IDC_EDIT_Area_Weight, m_ifAreaWeight);
	DDX_Text(pDX, IDC_EDIT_Volume_Preserve, m_ifVolumePreserve);
	DDX_Text(pDX, IDC_EDIT_Variation_Normal, m_ifVariationNormal);
}


BEGIN_MESSAGE_MAP(CFormCommandView, CFormView)
	//{{AFX_MSG_MAP(CFormCommandView)
	ON_WM_LBUTTONUP()
	ON_BN_CLICKED(IDC_CHECK1, OnCheckAntialias)
	ON_BN_CLICKED(IDC_CHECK2, OnCheckSmooth)
	ON_BN_CLICKED(IDC_CHECK3, OnCheckLight)
	ON_BN_CLICKED(IDC_RADIO1, OnRadioModel10)
	ON_BN_CLICKED(IDC_RADIO2, OnRadioModel11)
	ON_BN_CLICKED(IDC_RADIO3, OnRadioModel12)
	ON_WM_HSCROLL()
	ON_BN_CLICKED(IDC_CHECK5, OnCheckLinkScale)
	ON_WM_PAINT()
	ON_BN_CLICKED(IDC_CHECK4, OnCheckVRotate)
	//}}AFX_MSG_MAP
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER1, OnNMCustomdrawSlider1)
//	ON_EN_CHANGE(IDC_EDIT_ISOVALUE, OnEnChangeEditIsovalue)
ON_BN_CLICKED(IDC_BUTTON_COFF_OK, OnBnClickedButtonCoffOk)

END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CFormCommandView diagnostics

#ifdef _DEBUG
void CFormCommandView::AssertValid() const
{
	CFormView::AssertValid();
}

void CFormCommandView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CFormCommandView message handlers

void CFormCommandView::OnLButtonUp(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	CRect rect;
	CMcubeApp* pApp=(CMcubeApp*)AfxGetApp();
	m_ControlBackColor.GetWindowRect(&rect);
	ScreenToClient(&rect);
	if(rect.PtInRect(point))
	{
		CColorDialog dlg(BackColor);
		if(dlg.DoModal()==IDOK)
		{ 
			BackColor=dlg.GetColor();
			CMcubeView *pView=(CMcubeView*)GetMcubeView();
			pView->m_ClearColorRed=(float)GetRValue(BackColor)/255;
			pView->m_ClearColorGreen=(float)GetGValue(BackColor)/255;
			pView->m_ClearColorBlue=(float)GetBValue(BackColor)/255;
			this->InvalidateRect(&rect,FALSE);
			pView->InvalidateRect(NULL,FALSE);
		}
	}
	m_ControlColorLightAmbient.GetWindowRect(&rect);
	ScreenToClient(&rect);
	if(rect.PtInRect(point))
	{
		CColorDialog dlg(AmbientColor);
		if(dlg.DoModal()==IDOK)
		{
			AmbientColor=dlg.GetColor();
			CMcubeView *pView=(CMcubeView *)GetMcubeView();
			float r=(float)GetRValue(AmbientColor)/255.0f;
			float g=(float)GetGValue(AmbientColor)/255.0f;
			float b=(float)GetBValue(AmbientColor)/255.0f;
			float ambientProperties[]={r,g,b,1.0f};
			glLightfv(GL_LIGHT0,GL_AMBIENT,ambientProperties);
			this->InvalidateRect(&rect,FALSE);
			pView->InvalidateRect(NULL,FALSE);			
		}
	}
	
	m_ControlColorLightSpecular.GetWindowRect(&rect);
	ScreenToClient(&rect);
	if(rect.PtInRect(point))
	{
		CColorDialog dlg(SpecularColor);
		if(dlg.DoModal()==IDOK)
		{
			SpecularColor=dlg.GetColor();
			CMcubeView *pView=(CMcubeView *)GetMcubeView();
			float r=(float)GetRValue(SpecularColor)/255.0f;
			float g=(float)GetGValue(SpecularColor)/255.0f;
			float b=(float)GetBValue(SpecularColor)/255.0f;
			float specularProperties[]={r,g,b,1.0f};
			glLightfv(GL_LIGHT0,GL_SPECULAR,specularProperties);
			this->InvalidateRect(&rect,FALSE);
			pView->InvalidateRect(NULL,FALSE);
			
		}
	}
	
	m_ControlColorLightDiffuse.GetWindowRect(&rect);
	ScreenToClient(&rect);
	if(rect.PtInRect(point))
	{
		CColorDialog dlg(DiffuseColor);
		if(dlg.DoModal()==IDOK)
		{
			DiffuseColor=dlg.GetColor();
			CMcubeView *pView=(CMcubeView *)GetMcubeView();
			float r=(float)GetRValue(DiffuseColor)/255.0f;
			float g=(float)GetGValue(DiffuseColor)/255.0f;
			float b=(float)GetBValue(DiffuseColor)/255.0f;
			float diffuseProperties[]={r,g,b,1.0f};
			glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuseProperties);
			this->InvalidateRect(&rect,FALSE);
			pView->InvalidateRect(NULL,FALSE);
			
		}
	}	
	CFormView::OnLButtonUp(nFlags, point);
}

void CFormCommandView::OnCheckAntialias() 
{
	// TODO: Add your control notification handler code here
	m_Antialias=!m_Antialias;
	if(m_Lighting)
	{
		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
		glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
		glLineWidth(1.5f);
	}
	else
	{
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
		glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
		glLineWidth(1.0f);
	}
	this->GetMcubeView()->InvalidateRect(NULL,FALSE);
}

void CFormCommandView::OnCheckSmooth() 
{
	// TODO: Add your control notification handler code here
	m_Smooth=!m_Smooth;
	if(m_Smooth)
		glEnable(GL_SMOOTH);
	else
		glDisable(GL_FLAT);
	this->GetMcubeView()->InvalidateRect(NULL,FALSE);
}

void CFormCommandView::OnCheckLight() 
{
	// TODO: Add your control notification handler code here
	m_Lighting=!m_Lighting;
	if(m_Lighting)
		glEnable(GL_LIGHTING);
	else
		glDisable(GL_LIGHTING);
	this->GetMcubeView()->InvalidateRect(NULL,FALSE);
}

void CFormCommandView::OnRadioModel10() 
{
	// TODO: Add your control notification handler code here
	glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
	this->GetMcubeView()->InvalidateRect(NULL,FALSE);		
}

void CFormCommandView::OnRadioModel11() 
{
	// TODO: Add your control notification handler code here
	glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
	this->GetMcubeView()->InvalidateRect(NULL,FALSE);
}

void CFormCommandView::OnRadioModel12() 
{
	// TODO: Add your control notification handler code here
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	this->GetMcubeView()->InvalidateRect(NULL,FALSE);
	
}

void CFormCommandView::OnInitialUpdate() 
{
	CFormView::OnInitialUpdate();
	
	// TODO: Add your specialized code here and/or call the base class
	TRACE("Sliders:updating...\n");
	m_SliderScaleX.SetRange(1,100,TRUE);
	m_SliderScaleY.SetRange(1,100,TRUE);
    m_SliderScaleZ.SetRange(1,100,TRUE);
	m_SliderScaleX.SetPos(50);
    m_SliderScaleY.SetPos(50);
    m_SliderScaleZ.SetPos(50);
}

void CFormCommandView::OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar) 
{
	// TODO: Add your message handler code here and/or call default
	UpdateScale();
	GetMcubeView()->InvalidateRect(NULL,FALSE);
	CFormView::OnHScroll(nSBCode, nPos, pScrollBar);
}

BOOL CFormCommandView::UpdateScale()
{
    CMcubeView *pView=(CMcubeView *)GetMcubeView();
	pView->m_xScaling=(float)m_SliderScaleX.GetPos()/50.0f;
    pView->m_yScaling=(float)m_SliderScaleY.GetPos()/50.0f;
	pView->m_zScaling=(float)m_SliderScaleZ.GetPos()/50.0f;
	
	if(m_LinkScale)
	{
		m_SliderScaleY.SetPos(m_SliderScaleX.GetPos());
        m_SliderScaleZ.SetPos(m_SliderScaleX.GetPos());
		pView->m_yScaling=pView->m_zScaling=pView->m_xScaling;
	}
    return TRUE;
}

void CFormCommandView::OnCheckLinkScale() 
{
	// TODO: Add your control notification handler code here
	m_LinkScale=!m_LinkScale;
	if(m_LinkScale)
	{
		CMcubeView *pView=(CMcubeView *)GetMcubeView();
		m_SliderScaleY.SetPos(m_SliderScaleX.GetPos());
        m_SliderScaleZ.SetPos(m_SliderScaleX.GetPos());
		pView->m_yScaling=pView->m_zScaling=pView->m_xScaling;
	}
	m_SliderScaleY.EnableWindow(!m_LinkScale);
    m_SliderScaleY.EnableWindow(!m_LinkScale);
	GetMcubeView()->InvalidateRect(NULL,FALSE);
}

void CFormCommandView::OnPaint() 
{
	CPaintDC dc(this); // device context for painting
	
	// TODO: Add your message handler code here
	CRect rect;
	m_ControlBackColor.GetWindowRect(&rect);
	ScreenToClient(&rect);
	CBrush BrushBack(BackColor);
	dc.FillRect(&rect,&BrushBack);
	
	m_ControlColorLightDiffuse.GetWindowRect(&rect);
	ScreenToClient(&rect);
	CBrush BrushLightAmbient(DiffuseColor);
	dc.FillRect(&rect,&BrushLightAmbient);
	
	m_ControlColorLightAmbient.GetWindowRect(&rect);
	ScreenToClient(&rect);
	CBrush BrushLightDiffuse(AmbientColor);
	dc.FillRect(&rect,&BrushLightDiffuse);
	
	m_ControlColorLightSpecular.GetWindowRect(&rect);
	ScreenToClient(&rect);
	CBrush BrushLightSpecular(SpecularColor);
	dc.FillRect(&rect,&BrushLightSpecular);	
	// Do not call CFormView::OnPaint() for painting messages
}

void CFormCommandView::OnCheckVRotate() 
{
	// TODO: Add your control notification handler code here
	m_VRotate=!m_VRotate;
	CMcubeView *pView=(CMcubeView *)GetMcubeView();
	if(m_VRotate)
		pView->SetTimer(1,10,NULL);
	else
		pView->KillTimer(1);
}

CView* CFormCommandView::GetMcubeView()
{
	CMcubeApp *pApp=(CMcubeApp *)AfxGetApp();
	CMainFrame *pFrame=(CMainFrame *)pApp->m_pMainWnd;
	CView * pView=(CView *)pFrame->m_wndSplitter.GetPane(0,1);
	return pView; 
}


void CFormCommandView::OnNMCustomdrawSlider1(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	// TODO: Add your control notification handler code here
	*pResult = 0;
}



//void CFormCommandView::OnEnChangeEditIsovalue()
//{
//	// TODO:  If this is a RICHEDIT control, the control will not
//	// send this notification unless you override the CFormView::OnInitDialog()
//	// function and call CRichEditCtrl().SetEventMask()
//	// with the ENM_CHANGE flag ORed into the mask.
//
//	// TODO:  Add your control notification handler code here
//
//	HWND hWnd=GetSafeHwnd();
//	HDC hDC=::GetDC(hWnd);
//	CMcubeDoc* pMyDoc=GetDocument();
//	pMyDoc->m_Isolevel=m_IsoValue;
//	return;
//
//}

void CFormCommandView::OnBnClickedButtonCoffOk()
{
	// TODO: Add your control notification handler code here
	HWND hWnd=GetSafeHwnd();
	HDC hDC=::GetDC(hWnd);
	CMcubeDoc* pMyDoc=(CMcubeDoc*)GetDocument();
	
	//::GetDlgItem(hWnd,IDC_EDIT_ISOVALUE)->
	this->UpdateData(true);
	
	pMyDoc->m_variationF=m_variationF;
	pMyDoc->m_variationG1=m_variationG1;
	pMyDoc->m_variationG2=m_variationG2;
	pMyDoc->m_variationH1=m_variationH1;
	pMyDoc->m_variationH2=m_variationH2;
	pMyDoc->m_variationH3=m_variationH3;
	pMyDoc->m_K=m_K;
	pMyDoc->m_indexFunction=m_indexFunction;
	pMyDoc->m_meanShift=m_meanShift;
	pMyDoc->m_variationNstop=m_variationNstop;
	pMyDoc->m_iterativeTimes=m_iterativeTimes;
	pMyDoc->m_noiseCoffecient=m_noiseCoff;
	pMyDoc->m_functionVariation=m_funcitonVariation;
	pMyDoc->m_visualization=m_visualization;
	pMyDoc->m_noiseFuncCoff=m_noiseFuncCoff;
	pMyDoc->m_thresholdColor=m_thresholdColor;
	pMyDoc->m_functionGradientWide=m_functionGradientWide;
	pMyDoc->m_functionGradientThreshold=m_functionGradientThreshold;
	pMyDoc->m_thresholdDistanceNormal=m_thresholdDistanceNormal;
	pMyDoc->m_thresholdDistanceGradientFunction=m_thresholdDistanceGradientFunction;
	pMyDoc->m_timeStep=m_timeStep;
	pMyDoc->m_loadConstant=m_loadConstant;
	pMyDoc->m_positionError=m_positionError;
	pMyDoc->m_radius=m_radius;
	pMyDoc->m_stopN=m_stopN;
	pMyDoc->m_maxIter=m_maxIter;
	pMyDoc->m_curveThreshold=m_curveThreshold;
	pMyDoc->m_tangentOrManifold=m_tangentOrManifold;
	pMyDoc->m_ifNormalWeight=m_ifNormalWeight;
	pMyDoc->m_ifAreaWeight=m_ifAreaWeight;
	pMyDoc->m_ifVolumePreserve=m_ifVolumePreserve;
	pMyDoc->m_ifVariationNormal=m_ifVariationNormal;
	return;
}

