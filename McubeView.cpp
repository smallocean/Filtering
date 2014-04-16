// McubeView.cpp : implementation of the CMcubeView class
//

#include "stdafx.h"
#include "Mcube.h"

#include "McubeDoc.h"
#include "McubeView.h"
#include "GL/gl.h"
#include "GL/glu.h"
//#include "GL/glaux.h"
#include "CIsosurface.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#define  TEXTUREWIDTH 64
#define  TEXTUREHEIGHT 64
GLubyte Texture[TEXTUREWIDTH][TEXTUREHEIGHT][3];


/////////////////////////////////////////////////////////////////////////////
// CMcubeView

IMPLEMENT_DYNCREATE(CMcubeView, CView)

BEGIN_MESSAGE_MAP(CMcubeView, CView)
	//{{AFX_MSG_MAP(CMcubeView)
	ON_WM_CREATE()
	ON_WM_PAINT()
	ON_WM_SIZE()
	ON_WM_LBUTTONUP()
	ON_WM_MOUSEMOVE()
	ON_WM_LBUTTONDOWN()
	ON_WM_TIMER()
	//}}AFX_MSG_MAP
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, CView::OnFilePrintPreview)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CMcubeView construction/destruction

CMcubeView::CMcubeView()
{
	// TODO: add construction code here
	m_hGLContext=NULL;
	m_GLPixelIndex=0;
	
	m_LeftButtonDown=FALSE;
	m_RightButtonDown=FALSE;
	//m_CursorRotation=AfxGetApp()->LoadCursor(IDC_CURSOR_ROTATION);
	CMcubeApp *pApp=(CMcubeApp *)AfxGetApp();
	InitGeometry();
}

CMcubeView::~CMcubeView()
{
}

BOOL CMcubeView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

/////////////////////////////////////////////////////////////////////////////
// CMcubeView drawing

void CMcubeView::OnDraw(CDC* pDC)
{
	CMcubeDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	// TODO: add draw code for native data here
}

/////////////////////////////////////////////////////////////////////////////
// CMcubeView printing

BOOL CMcubeView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CMcubeView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add extra initialization before printing
}

void CMcubeView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}

/////////////////////////////////////////////////////////////////////////////
// CMcubeView diagnostics

#ifdef _DEBUG
void CMcubeView::AssertValid() const
{
	CView::AssertValid();
}

void CMcubeView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CMcubeDoc* CMcubeView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CMcubeDoc)));
	return (CMcubeDoc*)m_pDocument;
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CMcubeView message handlers

void CMcubeView::InitGeometry()
{
	m_xRotation=0.0f;
	m_yRotation=0.0f;
	m_xTraslation=0.0f;
	m_yTraslation=0.0f;
	m_zTraslation=-5.0f;
	m_xScaling=1.0f;
	m_yScaling=1.0f;
	m_zScaling=1.0f;
}

int CMcubeView::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	if (CView::OnCreate(lpCreateStruct) == -1)
		return -1;
	
	// TODO: Add your specialized creation code here
	HWND hWnd=GetSafeHwnd();
	HDC hDC=::GetDC(hWnd);
	
	if(SetWindowPixelFormat(hDC)==FALSE)
		return 0;
	if(CreateViewGLContext(hDC)==FALSE)
		return 0;
	
	glPolygonMode(GL_FRONT,GL_LINE);
	glPolygonMode(GL_BACK,GL_LINE);
	
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_NORMALIZE);
	GLfloat mat_ambient[]={0,0,0,0.15};
	GLfloat mat_specular[]={1.0,1.0,1.0,1.0};
	GLfloat mat_shiness[]={15.0};

	glMaterialfv(GL_FRONT,GL_AMBIENT,mat_ambient);
	glMaterialfv(GL_FRONT,GL_SPECULAR,mat_specular);
	glMaterialfv(GL_FRONT,GL_SHININESS,mat_shiness);

	GLfloat ambientProperties[]={0.7f,0.7f,0.7f,1.0f};
	GLfloat diffuseProperties[]={0.8f,0.8f,0.8f,1.0f};
	GLfloat specularProperties[]={0.0f,0.f,0.f,1.0f};

	
	glClearDepth(1.0);

	glLightfv(GL_LIGHT0,GL_AMBIENT,ambientProperties);
	glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuseProperties);
	glLightfv(GL_LIGHT0,GL_SPECULAR,specularProperties);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,5.0f);

	GLfloat position[]={0.0,0.0,1.,1};

	glLightfv(GL_LIGHT0,GL_POSITION,position);
    glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);
	return 0;
}

BOOL CMcubeView::SetWindowPixelFormat(HDC hDC)
{
    PIXELFORMATDESCRIPTOR pixelDesc;
	pixelDesc.nSize=sizeof(PIXELFORMATDESCRIPTOR);
	pixelDesc.nVersion=1;
	pixelDesc.dwFlags=PFD_DRAW_TO_WINDOW|PFD_SUPPORT_OPENGL|PFD_STEREO_DONTCARE|PFD_DOUBLEBUFFER;
	pixelDesc.iPixelType=PFD_TYPE_RGBA;
	pixelDesc.cColorBits=32;
	pixelDesc.cRedShift=0;
	pixelDesc.cGreenBits=0;
	pixelDesc.cGreenShift=0;
	pixelDesc.cBlueBits=0;
	pixelDesc.cBlueShift=0;
	pixelDesc.cAlphaBits=0;
	pixelDesc.cAlphaShift=0;
	pixelDesc.cAccumBits=64;
	pixelDesc.cAccumRedBits=16;
	pixelDesc.cAccumGreenBits=16;
	pixelDesc.cAccumBlueBits=16;
	pixelDesc.cAccumAlphaBits=0;
	pixelDesc.cDepthBits=32;
	pixelDesc.cDepthBits=32;
	pixelDesc.cAuxBuffers=0;
	pixelDesc.iLayerType=PFD_MAIN_PLANE;
    pixelDesc.bReserved=0;
	pixelDesc.dwLayerMask=0;
	pixelDesc.dwVisibleMask=0;
	pixelDesc.dwDamageMask=0;
	
	m_GLPixelIndex=ChoosePixelFormat(hDC,&pixelDesc);
	if(m_GLPixelIndex==0)
	{
		m_GLPixelIndex=1;
		if(DescribePixelFormat(hDC,m_GLPixelIndex,sizeof(PIXELFORMATDESCRIPTOR),&pixelDesc)==0)
			return FALSE;

	}
	if(!SetPixelFormat(hDC,m_GLPixelIndex,&pixelDesc))
		return FALSE;
	return TRUE;
}

BOOL CMcubeView::CreateViewGLContext(HDC hDC)
{
	m_hGLContext=wglCreateContext(hDC);
	if(m_hGLContext==NULL)
		return FALSE;
	
	if(wglMakeCurrent(hDC,m_hGLContext)==FALSE)
		return FALSE;
	return TRUE;

}

void CMcubeView::OnPaint() 
{
	CPaintDC dc(this); // device context for painting
	
	// TODO: Add your message handler code here

	HWND hWnd=GetSafeHwnd();
	HDC hDC=::GetDC(hWnd);
	CMcubeDoc* pMyDoc=GetDocument();

	wglMakeCurrent(hDC,m_hGLContext);
	  
   glClearColor(m_ClearColorRed,m_ClearColorGreen,m_ClearColorBlue,1.0f);
//	glClearColor(1.0,0.0,1.0,1.0f);

    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
	

	glPushMatrix();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
    glTranslated(m_xTraslation,m_yTraslation,-200);
	GLfloat position[]={-0,-0,-20000,1};
	GLfloat position1[]={0,-20000,-20000,1};
	GLfloat position2[]={0,-20000,-20000,1};
	GLfloat position3[]={-20000,0,-20000,1};
	GLfloat position4[]={-20000,0,-20000,1};

 	glLightfv(GL_LIGHT0,GL_POSITION,position);
	glLightfv(GL_LIGHT1,GL_POSITION,position1);
	glLightfv(GL_LIGHT2,GL_POSITION,position2);
	glLightfv(GL_LIGHT3,GL_POSITION,position3);
	glLightfv(GL_LIGHT4,GL_POSITION,position4);
	glRotatef(m_yRotation,1.0f,0.0f,0.0f);
	glRotatef(m_xRotation,0.0f,1.0f,0.0f);
	
	glScalef(m_xScaling,m_yScaling,m_zScaling);
	//glLightfv(GL_LIGHT0,GL_POSITION,position);

	//glTranslated(-256,-256,-256);
	//glLightfv(GL_LIGHT0,GL_POSITION,position);
	//glLightfv(GL_LIGHT1,GL_POSITION,position1);
	//glLightfv(GL_LIGHT2,GL_POSITION,position2);
	//glLightfv(GL_LIGHT3,GL_POSITION,position3);
	//glLightfv(GL_LIGHT4,GL_POSITION,position4);
	glTranslated(-(pMyDoc->x_max+pMyDoc->x_min)/2,-(pMyDoc->y_max+pMyDoc->y_min)/2,-(pMyDoc->z_max+pMyDoc->z_min)/2);
	
	//glTranslated(-pMyDoc->m_fVolLengthX/2,-pMyDoc->m_fVolLengthY/2,-pMyDoc->m_fVolLengthZ/2);


// 	glDrawBuffer(GL_BACK);
//	glDisable(GL_LIGHTING);

	//if(pMyDoc->m_visualization==0){
	//	glEnable(GL_LIGHTING);
	//	for(int i=0;i<pMyDoc->m_numOfPoints;i++){
	//		glBegin(GL_POINTS);
	//		glNormal3f(pMyDoc->m_normals[i*3]
	//		,pMyDoc->m_normals[i*3+1]
	//		,pMyDoc->m_normals[i*3+2]);
	//		//glColor3f(0.5,0.5,0.5);
	//		glVertex3f(pMyDoc->m_pointSets[i*3]
	//		,pMyDoc->m_pointSets[i*3+1]
	//		,pMyDoc->m_pointSets[i*3+2]);

	//		glEnd();

	//	}

	//}
	//else if (pMyDoc->m_visualization==1){
	//	glDisable(GL_LIGHTING);
	//	for(int i=0;i<pMyDoc->m_numOfPoints;i++){
	//		glBegin(GL_POINTS);
	//		/*glNormal3f(pMyDoc->m_normals[i*3]
	//		,pMyDoc->m_normals[i*3+1]
	//		,pMyDoc->m_normals[i*3+2]);*/
	//		glColor3f(pMyDoc->m_colors[i*3],pMyDoc->m_colors[i*3+1],pMyDoc->m_colors[i*3+2]);
	//		glVertex3f(pMyDoc->m_pointSets[i*3]
	//		,pMyDoc->m_pointSets[i*3+1]
	//		,pMyDoc->m_pointSets[i*3+2]);

	//		glEnd();

	//	}

	//}
	if(pMyDoc->m_visualization==0){
		glEnable(GL_LIGHTING);
		
		for (int i=0;i<pMyDoc->m_numOfTriangles;i++)
		{
			int j=pMyDoc->m_triangles[i*3];

			glBegin(GL_TRIANGLES );
			glNormal3f(pMyDoc->m_normals[j*3]
			,pMyDoc->m_normals[j*3+1]
			,pMyDoc->m_normals[j*3+2]);
			glVertex3f(pMyDoc->m_pointSets[j*3]
			,pMyDoc->m_pointSets[j*3+1]
			,pMyDoc->m_pointSets[j*3+2]);
			j=pMyDoc->m_triangles[i*3+1];
			glNormal3f(pMyDoc->m_normals[j*3]
			,pMyDoc->m_normals[j*3+1]
			,pMyDoc->m_normals[j*3+2]);
			glVertex3f(pMyDoc->m_pointSets[j*3]
			,pMyDoc->m_pointSets[j*3+1]
			,pMyDoc->m_pointSets[j*3+2]);

			j=pMyDoc->m_triangles[i*3+2];
			
			glNormal3f(pMyDoc->m_normals[j*3]
			,pMyDoc->m_normals[j*3+1]
			,pMyDoc->m_normals[j*3+2]);
			glVertex3f(pMyDoc->m_pointSets[j*3]
			,pMyDoc->m_pointSets[j*3+1]
			,pMyDoc->m_pointSets[j*3+2]);
			glEnd();
		}

	}
	else{
		glDisable(GL_LIGHTING);
		for (int i=0;i<pMyDoc->m_numOfTriangles;i++)
		{
			int j=pMyDoc->m_triangles[i*3];

			glBegin(GL_TRIANGLES );
			glColor3f(pMyDoc->m_colors[j*3],pMyDoc->m_colors[j*3+1],pMyDoc->m_colors[j*3+2]);
			glVertex3f(pMyDoc->m_pointSets[j*3]
			,pMyDoc->m_pointSets[j*3+1]
			,pMyDoc->m_pointSets[j*3+2]);
			j=pMyDoc->m_triangles[i*3+1];
			glColor3f(pMyDoc->m_colors[j*3],pMyDoc->m_colors[j*3+1],pMyDoc->m_colors[j*3+2]);
			
			glVertex3f(pMyDoc->m_pointSets[j*3]
			,pMyDoc->m_pointSets[j*3+1]
			,pMyDoc->m_pointSets[j*3+2]);
			j=pMyDoc->m_triangles[i*3+2];
			
			glColor3f(pMyDoc->m_colors[j*3],pMyDoc->m_colors[j*3+1],pMyDoc->m_colors[j*3+2]);
			glVertex3f(pMyDoc->m_pointSets[j*3]
			,pMyDoc->m_pointSets[j*3+1]
			,pMyDoc->m_pointSets[j*3+2]);
			glEnd();
		}
	}	
	glPopMatrix();
//  glReadBuffer(GL_BACK);
//	glCopyPixels(0,0,1024,768,GL_COLOR|GL_DEPTH);
    glFlush();
    SwapBuffers(hDC); 


	// glPopMatrix();
	// SwapBuffers(hDC);
	
	// Do not call CView::OnPaint() for painting messages
	dc.SetBkMode(TRANSPARENT);
	dc.SetTextColor(RGB(0,0,0));
	CString error,maxError, meanError;
	error.Format("%f", pMyDoc->m_maxError);
	maxError = "maxError    "+error ;
	dc.TextOut(20,20,maxError);
	error.Format("%f", pMyDoc->m_meanError);
	meanError= "meanError    "+error;
	dc.TextOut(20,50,meanError);

}

void CMcubeView::OnSize(UINT nType, int cx, int cy) 
{
	CView::OnSize(nType, cx, cy);
	
	// TODO: Add your message handler code here
	CSize size(cx,cy);
	double aspect;
	aspect=(cy=0)?(double)size.cx:(double)size.cx/(double)size.cy;
	
	glViewport(0,0,size.cx,size.cy);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45,aspect,1,15000.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glDrawBuffer(GL_BACK);
	glEnable(GL_DEPTH_TEST);
}

void CMcubeView::OnLButtonUp(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	m_LeftButtonDown=FALSE;
	CView::OnLButtonUp(nFlags, point);
	CView::OnLButtonUp(nFlags, point);
}

void CMcubeView::OnMouseMove(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	if(m_LeftButtonDown)
	{
		m_xRotation-=(float)(m_LeftDownPos.x-point.x)/3.0f;
		m_yRotation-=(float)(m_LeftDownPos.y-point.y)/3.0f;
		m_LeftDownPos=point;
		InvalidateRect(NULL,FALSE);
	}
	CView::OnMouseMove(nFlags, point);
}

void CMcubeView::OnLButtonDown(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	m_LeftButtonDown=TRUE;
	m_LeftDownPos=point;
	CView::OnLButtonDown(nFlags, point);
	CView::OnLButtonDown(nFlags, point);
}

void CMcubeView::OnTimer(UINT nIDEvent) 
{
	// TODO: Add your message handler code here and/or call default
	switch(nIDEvent) {
	case 0:
		break;
	case 1:
		m_yRotation+=5.0f;
        InvalidateRect(NULL,FALSE);
		break;
	default:
		{}
	}
	CView::OnTimer(nIDEvent);
}

void CMcubeView::myinit()
{
   // makeTexture();
   //	glPixelStorei(GL_UNPACK_ALIGNMENT,1);
   
	//glTexImage2D(GL_TEXTURE_2D,0,3,TEXTUREWIDTH,TEXTUREHEIGHT,0,GL_RGB,GL_UNSIGNED_BYTE,&Texture[0][0][0]);
	//	glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
	//	glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
	//	glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
	//	glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	//	glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_DECAL);
	//
	//	glEnable(GL_TEXTURE_2D);
	//	glEnable(GL_TEXTURE_2D);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);


}



void CMcubeView::makeTexture()
{
	int i,j,r,g,b;
	for(i=0;i<TEXTUREWIDTH;i++)
	{
		for(j=0;j<TEXTUREHEIGHT;j++)
		{
			r=i*j%255;
			g=(4*i)%255;
			b=(4*j)%255;
			Texture[i][j][0]=(GLubyte)r;
            Texture[i][j][1]=(GLubyte)g;
			Texture[i][j][2]=(GLubyte)b;		
		}
	}

}