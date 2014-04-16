// McubeDoc.cpp : implementation of the CMcubeDoc class
//

#include "stdafx.h"
#include "Mcube.h"
#include <math.h>
#include "McubeDoc.h"
#include "utils.h"
#include "CIsoSurface.h"
#include "McubeView.h"
#include "DialogIntepolate.h"
#include "Vectors.h"
#include ".\mcubedoc.h"
#include "FormCommandView.h"
#include "GenerateNoise.h"
#include "NonShiftBilateralFilter.h"
//#include "BilateralFilter.h"
#include "FilterBilateral1.h"
#include "Delaunay.h"
#include "MeshlessFilter.h"
#include"OneFreedomMeshless.h"
#include "RBFmeshless.h"
#include "AnisoRBFmeshless.h"
#include "MLSfilter.h"
#include "AnisoRBFfilterSurfaceFunction.h"

 
#include "mesh.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CMcubeDoc

IMPLEMENT_DYNCREATE(CMcubeDoc, CDocument)

BEGIN_MESSAGE_MAP(CMcubeDoc, CDocument)
	//{{AFX_MSG_MAP(CMcubeDoc)
	ON_COMMAND(ID_MCUBE_CALCULATER, OnMcubeCalculater)
	ON_COMMAND(ID_VECTOR_INTEPOLATE, OnVectorIntepolate)

//	ON_COMMAND(ID_MCUBE_NODEMOVE, OnMcubeNodemove)
	//}}AFX_MSG_MAP
	ON_EN_CHANGE(IDC_EDIT_ISOVALUE, OnEnChangeEditIsovalue)
	ON_BN_CLICKED(IDC_BUTTON1, OnBnClickedButton1)
	ON_COMMAND(ID_MESH_SKIN, OnMeshSkin)
	ON_COMMAND(ID_MESH_BONE, OnMeshBone)
	ON_COMMAND(ID_Filter_Bilateral, OnFilterBilateral)
	ON_COMMAND(ID_Filter_Addnoise_Geometry, OnFilterAddnoiseGeometry)
	ON_COMMAND(ID_Filter_Addnoise_Function, OnFilterAddnoiseFunction)
	ON_COMMAND(ID_Filter_Non_Shift, OnFilterNonShift)
	ON_COMMAND(ID_Meshless_Filter, OnMeshlessFilter)
	ON_COMMAND(ID_Filter_Aniso_Meshless, OnFilterAnisoMeshless)
	ON_COMMAND(ID_Filter_MLS, OnFilterMls)
	ON_COMMAND(ID_Aniso_Sur_Fun_Filter, OnAnisoSurFunFilter)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CMcubeDoc construction/destruction

CMcubeDoc::CMcubeDoc()
{
	// TODO: add one-time construction code here
    m_nSizeX=0;
	m_nSizeY=0;
	m_nSizeZ=0;
	m_fCellLengthX=1;
	m_fCellLengthY=1;
	m_fCellLengthZ=1;
	m_visualization=0;
	double a=sqrt(-0.0);
	m_pVol=NULL;
   	m_Isolevel=1600;
	m_colors=NULL;
	temp_normals=NULL;
	//m_nMeshTriangleTwo=NULL;
 //   m_nMeshTriangleOne=NULL;
	float  w= (float)(1 << 1);
	//生成cube
	//generateCube();
	//生成球

	//m_nSizeX=m_nSizeY=m_nSizeZ=20;
	//m_pVol=new float[m_nSizeX*m_nSizeY*m_nSizeZ];	
	//for (int i=0;i<m_nSizeX;i++)
	//	for(int j=0;j<m_nSizeY;j++)
	//		for(int k=0;k<m_nSizeZ;k++)
	//		{
	//			m_pVol[i*m_nSizeZ*m_nSizeY+j*m_nSizeZ+k]=(i-(int)(m_nSizeX/2))*(i-(int)(m_nSizeX/2))+(j-(int)(m_nSizeY/2))*(j-(int)(m_nSizeY/2))+(k-(int)(m_nSizeZ/2))*(k-(int)(m_nSizeZ/2));
	//		}

	//		m_fCellLengthX=m_fCellLengthY=m_fCellLengthZ=1;

	//m_pVol=new float[100*100*100];
	//
    // for (int i=0;i<100;i++)
	//	for(int j=0;j<100;j++)
	//		for(int k=0;k<100;k++)
	//		{
	//			m_pVol[i*100*100+j*100+k]=(i-50)*(i-50)+(j-50)*(j-50)+(k-50)*(k-50);
	//		}
	//m_nSizeX=m_nSizeY=m_nSizeZ=100;
	//m_fCellLengthX=m_fCellLengthY=m_fCellLengthZ=1;
	//::GenerateMesh();

	std::vector<Point> tempPointSets;
	std::vector<Triangle> tempTriangleSets;
	std::vector<Point> tempColorSets;

  //CString filename_pw="D:\\qin\\mesh data\\fandiskThree\\point.txt";
  // 	CString filename_pw="D:\\qin\\mesh data\\fandiskpoints.txt";
	/*CString filename_pw = "D:\\qin\\mesh data\\points.txt";
	CString filename_pw = "D:\\qin\\mesh data\\face2 data\\F0002_Vertex.txt";
	CString filename_pw = "D:\\qin\\mesh data\\image data\\point.txt";*/
	//CString filename_pw="D:\\qin\\mesh data\\dragon_recon\\dragon_recon_vrip_res2_point.txt";
	//CString filename_pw="D:\\qin\\mesh data\\dragon_recon\\dragon_recon_vrip_point.txt";
 	CString filename_pw=".\\mesh data\\bunny\\bunny_vertex.txt";
	//  CString  filename_pw="D:\\qin\\mesh data\\face myself\\face myself.txt";
	//CString  filename_pw="D:\\qin\\mesh data\\face\\point.txt";
	
  //CString filename_pw="D:\\qin\\mesh data\\bunny\\bun_zipper_res2_point.txt";
  //  CString filename_pw="C:\\Program Files\\MATLAB71\\work\\result.txt";

 // CString filename_pw="D:\\qin\\mesh data\\face30\\F0030_Vertex.txt";
	//CString filename_pw="D:\\qin\\mesh data\\max-planck\\point.txt";
	//CString filename_pw="D:\\qin\\mesh data\\hand\\point.txt";

	FILE *fpout;

	if((fpout = fopen(filename_pw, "r")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{
		while(!feof(fpout)){
			Point tempPoint;
			float temp1,temp2;
			char tempChar[10];
			//fscanf(fpout,"%f %f %f %f %f",&tempPoint.vertexPosition[0],&tempPoint.vertexPosition[1],&tempPoint.vertexPosition[2],temp1,temp2);
			//fscanf(fpout,"%f %f %f %s",&tempPoint.vertexPosition[0],&tempPoint.vertexPosition[1],&tempPoint.vertexPosition[2],tempChar);
		    fscanf(fpout,"%f %f %f %f %f",&tempPoint.vertexPosition[0],&tempPoint.vertexPosition[1],&tempPoint.vertexPosition[2],&temp1,&temp2);
			tempPointSets.push_back(tempPoint);
		}
		fclose(fpout);
	}


  /* filename_pw = "D:\\qin\\mesh data\\colors.txt";
	 filename_pw = "D:\\qin\\mesh data\\face data\\M0001_Texture.txt";
	filename_pw = "D:\\qin\\mesh data\\image data\\lena_std_color.txt";*/
	//filename_pw="D:\\qin\\mesh data\\face30\\F0030_Texture.txt";
	//if((fpout = fopen(filename_pw, "r")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	while(!feof(fpout)){
	//		Point tempColor;
	//		fscanf(fpout,"%f %f %f",&tempColor.vertexPosition[0],&tempColor.vertexPosition[1],&tempColor.vertexPosition[2]);
	//		tempColorSets.push_back(tempColor);
	//	}
	//	fclose(fpout);
	//}

  //  filename_pw="D:\\qin\\mesh data\\fandiskThree\\triangle.txt";
	//filename_pw="D:\\qin\\mesh data\\fandiskmesh.txt";
	/*filename_pw = "D:\\qin\\mesh data\\triangles.txt";
	filename_pw = "D:\\qin\\mesh data\\face2 data\\F0002_Triangle.txt";
	filename_pw = "D:\\qin\\mesh data\\image data\\triangle.txt";*/
	//filename_pw="D:\\qin\\mesh data\\dragon_recon\\dragon_recon_vrip_res2_triangle.txt";
	//filename_pw="D:\\qin\\mesh data\\dragon_recon\\dragon_recon_vrip_triangle.txt";
	 //filename_pw="D:\\qin\\mesh data\\bunny\\bun_zipper_res2_triangle.txt";
    filename_pw=".\\mesh data\\bunny\\bunny_triangle.txt";
	//filename_pw="D:\\qin\\mesh data\\face myself.txt";
  // filename_pw="D:\\qin\\mesh data\\face30\\F0030_Triangle.txt";
	//filename_pw="D:\\qin\\mesh data\\max-planck\\triangle.txt";
	//   filename_pw="D:\\qin\\mesh data\\hand\\triangle.txt";
	if((fpout = fopen(filename_pw, "r")) == NULL)
	{
		int dkjkd;
	//	MessageBox("can't open the file!");
	}
	else
	{
		while(!feof(fpout)){
			Triangle tempTriangle;
			char temp;
			//fscanf(fpout,"%c %d %d %d",&temp,&tempTriangle.vertexIndex[0],&tempTriangle.vertexIndex[1],&tempTriangle.vertexIndex[2]);
			fscanf(fpout,"%d %d %d",&tempTriangle.vertexIndex[0],&tempTriangle.vertexIndex[1],&tempTriangle.vertexIndex[2]);
			tempTriangleSets.push_back(tempTriangle);
		}
		fclose(fpout);
	}
	m_numOfPoints=tempPointSets.size();
	m_numOfTriangles=tempTriangleSets.size();
	m_pointSets=new float[m_numOfPoints*3];
	m_triangles=new int[m_numOfTriangles*3];
	m_colors=new float[m_numOfPoints*3];
	std::vector<Point>::iterator tempPointIterator=tempPointSets.begin();
	std::vector<Point>::iterator tempColorSetsIterator=tempColorSets.begin();
//	for(int i=0;tempPointIterator!=tempPointSets.end();tempPointIterator++,tempColorSetsIterator++,i++){
	for(int i=0;tempPointIterator!=tempPointSets.end();tempPointIterator++,i++){
		
		
		m_pointSets[i*3]=(*tempPointIterator).vertexPosition[0]*1000;
		m_pointSets[i*3+1]=(*tempPointIterator).vertexPosition[1]*1000;
		m_pointSets[i*3+2]=(*tempPointIterator).vertexPosition[2]*1000;
	/*	m_colors[i*3]=(*tempColorSetsIterator).vertexPosition[0];
		m_colors[i*3+1]=(*tempColorSetsIterator).vertexPosition[1];
		m_colors[i*3+2]=(*tempColorSetsIterator).vertexPosition[2];*/


	}
	std::vector<Triangle>::iterator tempTriangleIterator=tempTriangleSets.begin();
	for(int i=0;tempTriangleIterator!=tempTriangleSets.end();tempTriangleIterator++,i++){
	
		m_triangles[i*3]=(*tempTriangleIterator).vertexIndex[0];
		m_triangles[i*3+1]=(*tempTriangleIterator).vertexIndex[1];
		m_triangles[i*3+2]=(*tempTriangleIterator).vertexIndex[2];

	}
	//filename_pw="D:\\qin\\mesh data\\agnes-hand\\triangle1.txt";;
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfTriangles;i++){
	//		fprintf(fpout,"%d %d %d\n",m_triangles[i*3]-1,m_triangles[i*3+1]-1,m_triangles[i*3+2]-1);

	//	}
	//	fclose(fpout);
	//}

	tempPointSets.clear();
	/*tempColorSets.clear();
	tempTriangleSets.clear();*/

	

	//for(int i=0;i<m_numOfPoints;i++){
	//	for(int j=0;j<3;j++){
	//		m_colors[i*3+j]/=256;
	//	}
	//}

	x_max=x_min=m_pointSets[0];
	y_max=y_min=m_pointSets[1];
	z_max=z_min=m_pointSets[2];
	for(int i=0;i<m_numOfPoints;i++){
		if(m_pointSets[i*3]>x_max)
			x_max=m_pointSets[i*3];
		else{
			if(m_pointSets[i*3]<x_min)
				x_min=m_pointSets[i*3];
		}
		if(m_pointSets[i*3+1]>y_max)
			y_max=m_pointSets[i*3+1];
		else{
			if(m_pointSets[i*3+1]<y_min)
				y_min=m_pointSets[i*3+1];
		}
		if(m_pointSets[i*3+2]>z_max)
			z_max=m_pointSets[i*3+2];
		else{
			if(m_pointSets[i*3+2]<z_min)
				z_min=m_pointSets[i*3+2];
		}
	}
//	AddPositionGaussianNoise();





	//int m_nNormals = m_numOfPoints;
	int m_nNormals = m_numOfTriangles;

	//m_normals=new float[m_nNormals*3];
	m_normals=new float[m_nNormals*3];

	//filename_pw="D:\\qin\\mesh data\\face myself\\ShiftoriginalNormal.txt";
	//filename_pw="D:\\qin\\mesh data\\face\\ShiftoriginalNormal.txt";
	//filename_pw="D:\\qin\\mesh data\\face31\\ShiftoriginalNormal.txt";



	//if((fpout = fopen(filename_pw, "r")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	int i=0;
	//	while(!feof(fpout)){
	//		Point tempPoint;
	//		float temp1,temp2;
	//		char tempChar[10];
	//		//fscanf(fpout,"%f %f %f %f %f",&tempPoint.vertexPosition[0],&tempPoint.vertexPosition[1],&tempPoint.vertexPosition[2],&temp1,&temp2);
	//		fscanf(fpout,"%f %f %f",&tempPoint.vertexPosition[0],&tempPoint.vertexPosition[1],&tempPoint.vertexPosition[2],tempChar);
	//		m_normals[i*3]=tempPoint.vertexPosition[0];
	//		m_normals[i*3+1]=tempPoint.vertexPosition[1];
	//		m_normals[i*3+2]=tempPoint.vertexPosition[2];
	//		i++;			
	//	}
	//	fclose(fpout);
	//}



	//// Set all normals to 0.
	for (int i = 0; i < m_nNormals*3; i++) {
		m_normals[i]=0;
	}

	// Calculate normals.
	for (int i = 0; i < m_numOfTriangles; i++) {
		VECTOR3D vec1, vec2, normal;
		int id0, id1, id2;
		id0 = m_triangles[i*3];
		id1 = m_triangles[i*3+1];
		id2 = m_triangles[i*3+2];
		vec1[0] = m_pointSets[id1*3]- m_pointSets[id0*3];
		vec1[1] = m_pointSets[id1*3+1] - m_pointSets[id0*3+1];
		vec1[2] = m_pointSets[id1*3+2] - m_pointSets[id0*3+2];
		vec2[0] = m_pointSets[id2*3] - m_pointSets[id0*3];
		vec2[1] = m_pointSets[id2*3+1] - m_pointSets[id0*3+1];
		vec2[2] = m_pointSets[id2*3+2] - m_pointSets[id0*3+2];
		normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
		normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
		normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
		
		m_normals[i*3] += normal[0];
				m_normals[i*3+1] += normal[1];
						m_normals[i*3+2] += normal[2];

		//m_normals[id0*3+0] += normal[0];
		//m_normals[id0*3+1] += normal[1];
		//m_normals[id0*3+2] += normal[2];
		//m_normals[id1*3+0] += normal[0];
		//m_normals[id1*3+1] += normal[1];
		//m_normals[id1*3+2] += normal[2];
		//m_normals[id2*3+0] += normal[0];
		//m_normals[id2*3+1] += normal[1];
		//m_normals[id2*3+2] += normal[2];
	}

	// Normalize normals.
	for (int i = 0; i < m_nNormals; i++) {
		
		float length = sqrt(m_normals[i*3+0]*m_normals[i*3+0] + m_normals[i*3+1]*m_normals[i*3+1] + m_normals[i*3+2]*m_normals[i*3+2]);
		if(length<0.000001){
				m_normals[i*3] = 0;		m_normals[i*3+1] = 0;		m_normals[i*3+2] = 1;
		}
		else{
		m_normals[i*3] /= length;
		m_normals[i*3+1] /= length;
		m_normals[i*3+2] /= length;
		}
	}

	m_originalPointSets=new float[m_numOfPoints*3];
	m_originalOriginalPointSets=new float[m_numOfPoints*3];
	m_originalColors=new float[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints*3;i++){
		m_originalPointSets[i]=m_pointSets[i];
		m_originalOriginalPointSets[i]=m_pointSets[i];
		m_originalColors[i]=m_colors[i];
	}
		filename_pw = "D:\\Normals\\Face-Normal.txt";
		//FILE *fpout;
		if((fpout = fopen(filename_pw, "w")) == NULL)
		{
			int dkjkd;
			//MessageBox("can't open the file!");
		}
		else
		{
			//for(int i=0;i<m_numOfPoints;i++){
			for(int i = 0;i<m_nNormals;i++){
				fprintf(fpout,"%f %f %f\n",m_normals[i*3],m_normals[i*3+1],m_normals[i*3+2]);
	
			}
			fclose(fpout);
		}

//	CalculateBoundaryColorOfVertexs(m_numOfPoints,m_pointSets,m_normals);	
  //  GenerateTriangleFromImage();
//  ::GenerateMesh(m_numOfPoints,m_numOfTriangles,m_pointSets,m_triangles);
}

CMcubeDoc::~CMcubeDoc()
{
	delete m_pVol;
}

BOOL CMcubeDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)

}

/////////////////////////////////////////////////////////////////////////////
// CMcubeDoc serialization

void CMcubeDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}

/////////////////////////////////////////////////////////////////////////////
// CMcubeDoc diagnostics

#ifdef _DEBUG
void CMcubeDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CMcubeDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CMcubeDoc commands

BOOL CMcubeDoc::OnOpenDocument(LPCTSTR lpszPathName) 
{
//	if (!CDocument::OnOpenDocument(lpszPathName))
//		return FALSE;
	std::string name = lpszPathName;
//	m_volume.Read(name);

//	xDim dim = m_volume.GetDim();
//	xSpacing spacing = m_volume.GetSpacing();
//	xDataType type = m_volume.GetDataType();

//	m_nSizeX=dim.x();
//	m_nSizeY=dim.y();
//	m_nSizeZ=dim.z();
//	m_fCellLengthX=spacing.x();
///	m_fCellLengthY=spacing.y();
//	m_fCellLengthZ=spacing.z();
    
//	if(type == UnsignedInt8)
//		m_pVol = (uint8*)(m_volume.GetDataPointer());
    return TRUE;
}



void CMcubeDoc::OnMcubeCalculater() 
{
	// TODO: Add your command handler code here

    m_CIsoSurface.GenerateSurface( m_pVol, m_Isolevel,  m_nSizeX-1, m_nSizeY-1,  m_nSizeZ-1,  m_fCellLengthX,  m_fCellLengthY,m_fCellLengthZ);
    m_CIsoSurface.GetVolumeLengths(m_fVolLengthX, m_fVolLengthY, m_fVolLengthZ);
	m_numOfPoints=m_CIsoSurface.m_nVertices;
	m_numOfTriangles=m_CIsoSurface.m_nTriangles;
	m_pointSets=new float[m_numOfPoints*3];
	m_triangles=new int[m_numOfTriangles*3];
	m_normals=new float[m_numOfPoints*3];
	m_colors=new float[m_numOfPoints*3];
	m_originalColors=new float[m_numOfPoints*3];

	for(int i=0;i<m_numOfPoints;i++){
		for(int j=0;j<3;j++){
			m_pointSets[i*3+j]=m_CIsoSurface.m_ppt3dVertices[i][j];
			m_normals[i*3+j]=m_CIsoSurface.m_pvec3dNormals[i][j];

		}
		if(m_pointSets[i*3]>=0&&m_pointSets[i*3]<20){
			m_colors[i*3]=0;
			m_colors[i*3+1]=0;
			m_colors[i*3+2]=0;
		}
		if(m_pointSets[i*3]>=20&&m_pointSets[i*3]<40){
			m_colors[i*3]=0.5;
			m_colors[i*3+1]=0.5;
			m_colors[i*3+2]=0.5;
		}
		if(m_pointSets[i*3]>=40&&m_pointSets[i*3]<60){
			m_colors[i*3]=1;
			m_colors[i*3+1]=1;
			m_colors[i*3+2]=1;
		}
		if(m_pointSets[i*3]>=60&&m_pointSets[i*3]<80){
			m_colors[i*3]=0.5;
			m_colors[i*3+1]=0.5;
			m_colors[i*3+2]=0.5;
		}
		if(m_pointSets[i*3]>=80&&m_pointSets[i*3]<100){
			m_colors[i*3]=0;
			m_colors[i*3+1]=0;
			m_colors[i*3+2]=0;
		}
		
	}
	for(int i=0;i<m_numOfTriangles*3;i++){
		m_triangles[i]=m_CIsoSurface.m_piTriangleIndices[i];
	}
	m_CIsoSurface.DeleteSurface();
		m_originalPointSets=new float[m_numOfPoints*3];
		m_originalOriginalPointSets=new float[m_numOfPoints*3];
		for(int i=0;i<m_numOfPoints*3;i++){
			m_originalPointSets[i]=m_pointSets[i];
			m_originalOriginalPointSets[i]=m_pointSets[i];
			m_originalColors[i]=m_colors[i];
		}


 //   
	this->UpdateAllViews(NULL);

}


void CMcubeDoc::OnVectorIntepolate() 
{
	// TODO: Add your command handler code here
	
	int i,j,k;
	for( i=1;i<m_nSizeZ;i++)
		for(j=0;j<m_nSizeY-1;j++)
			for(k=0;k<m_nSizeX-1;k++)
	{
		m_pVol[2*i*m_nSizeX*m_nSizeY+j*m_nSizeX+k]=m_pVol[2*i*m_nSizeX*m_nSizeY+j*m_nSizeX+k];
		
	}
	for(i=1;i<m_nSizeZ;i++)
		for(j=0;j<m_nSizeY-1;j++)
			for(k=0;k<m_nSizeX;k++)
			{
				m_pVol[(2*i-1)*m_nSizeX*m_nSizeY+j*m_nSizeX+k]=
					(m_pVol[(2*i-2)*m_nSizeX*m_nSizeY+j*m_nSizeX+k]+m_pVol[2*i*m_nSizeX*m_nSizeY+j*m_nSizeX+k])/2;
			}

	m_nSizeZ=2*m_nSizeZ-1;

	m_fCellLengthZ/=2;
	return;
	
}

void CMcubeDoc::Filter()
{
	int i,j,k;
	for(i=1;i<m_nSizeZ-2;i++)
		for(j=1;j<m_nSizeY-2;j++)
			for(k=1;k<m_nSizeX-2;k++)
				m_pVol[i*m_nSizeX*m_nSizeY+j*m_nSizeX+k]=
                 m_pVol[i*m_nSizeX*m_nSizeY+j*m_nSizeX+k]/4+
				 m_pVol[(i-1)*m_nSizeX*m_nSizeY+j*m_nSizeX+k]/8+
				 m_pVol[(i+1)*m_nSizeX*m_nSizeY+j*m_nSizeX+k]/8+
				 m_pVol[i*m_nSizeX*m_nSizeY+(j-1)*m_nSizeX+k]/8+
				 m_pVol[i*m_nSizeX*m_nSizeY+(j+1)*m_nSizeX+k]/8+
				 m_pVol[i*m_nSizeX*m_nSizeY+j*m_nSizeX+k-1]/8+
                 m_pVol[i*m_nSizeX*m_nSizeY+j*m_nSizeX+k+1]/8;
}

/*
void CMcubeDoc::OnMcubeNodemove() 
{
	// TODO: Add your command handler code here
	int Node[20000][3];

	for (int i=0;i<m_nSizeX;i++)
	{
		for(int j=0;j<20000;j++)
			Node[j][0]=Node[j][1]=Node[j][2]=-1;
		float x=i*m_fCellLengthX;
		int number=0;
		for(int k=0;k<m_CIsoSurface.m_nVertices;k++)
		{
			if (x==m_CIsoSurface.m_ppt3dVertices[k][0])
			{
				Node[number][0]=k;
				Node[number][1]=m_CIsoSurface.m_ppt3dVertices[k][1];
                Node[number][2]=m_CIsoSurface.m_ppt3dVertices[k][2];
				number++;
			}
			for(int m=0;m<number-1;m++)
			{
				for(int n=0;n<m_CIsoSurface.m_nTriangles;n++)
				{
                      if (Node[m][k]==m_CIsoSurface.m_piTriangleIndices[n*3])
					  {
						  if(m_CIsoSurface.m_ppt3dVertices
							  [m_CIsoSurface.m_piTriangleIndices[n*3+1]][0]==x)
						  {
							  int a=Node[m+1][0];
						      int b=Node[m+1][1];
							  int c=Node[m+1][2];
							  for(int w=m+1;w<number-1;m++)
								  if (Node[w][0]==(m_CIsoSurface.m_piTriangleIndices[n*3+1]))
									  break;
							  Node[m+1][0]=Node[w][0];
							  Node[m+1][1]=Node[w][1];
							  Node[m+1][2]=Node[w][2];
							  Node[w][0]=a;
							  Node[w][1]=b;
							  Node[w][2]=c;
						  }
						  else if(m_CIsoSurface.m_ppt3dVertices
							  [m_CIsoSurface.m_piTriangleIndices[n*3+2]][0]==x)
						  {
                              int a=Node[m+1][0];
						      int b=Node[m+1][1];
							  int c=Node[m+1][2];
							  for(int w=m+1;w<number-1;m++)
								  if (Node[w][0]==[m_CIsoSurface.m_piTriangleIndices[n*3+1])
									  break;
							  Node[m+1][0]=Node[w][0];
							  Node[m+1][1]=Node[w][1];
							  Node[m+1][2]=Node[w][2];
							  Node[w][0]=a;
							  Node[w][1]=b;
							  Node[w][2]=c;

							  
				}




				
			}
		}


			   


	}

	return;
    
		
}

//指针函数模板类
//template <class T> class List
//{
//	 private:
//		 List<T> *next;
//	 public:
//		 T data;
//		 //构造函数
//		 List(const T& item,List<T>*ptrnext=NULL);
//		 //在当前节点之后插入指针p所指节点
//		 void InsertAfter(List<T> *p);
//		 //删除当前节点的后继节点并返回指向被删除节点的指针
//		 List<T> *DeleteAfter(void);
//		 //返回指向当前节点的后继节点的指针
//		 List<T> *NextNode(void) const;
//};

//找到切片上的等值
struct NodeOnSlice
{
	int id;
	float x;
	float y;
	float z;
};
NodeSlice* SearchNodeOnlice(POINT3D* m_Vertice, float x, float y)
{
	NodeOnSlice Node[3000];
	for(int i=0;i<3000;i++)
	{
		Node[i].id=-1;
		Node[i].x=-1;
	    Node[i].y=-1;
		Node[i].z=-1;
	}
	
	if(y==-1)
	{
		for(int k=0;k<m_CIsoSurface.m_nVertices;k++)
		{
		   	int number=0;
			if (x==m_CIsoSurface.m_ppt3dVertices[k][0])
			{
					Node[number].id=k;
					Node[number].x=m_CIsoSurface.m_ppt3dVertices[k][0];
					Node[number].y=m_CIsoSurface.m_ppt3dVertices[k][1];
					Node[number].z=m_CIsoSurface.m_ppt3dVertices[k][2];
					number++;
			}
		}
	}
	else
	{
        for(int k=0;k<m_CIsoSurface.m_nVertices;k++)
		{
		   	int number=0;
			if (y==m_CIsoSurface.m_ppt3dVertices[k][0])
			{
					Node[number].id=k;
					Node[number].x=m_CIsoSurface.m_ppt3dVertices[k][0];
					Node[number].y=m_CIsoSurface.m_ppt3dVertices[k][1];
					Node[number].z=m_CIsoSurface.m_ppt3dVertices[k][2];
					number++;
			}
		}
	}
	return Node;
}

*/

void CMcubeDoc::OnEnChangeEditIsovalue()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDocument::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
//	m_Isolevel=CFormCommandView.m_IsoValue;
	return;
}

void CMcubeDoc::OnBnClickedButton1()
{
	// TODO: Add your control notification handler code here
	CView *pView;
	POSITION pos=this->GetFirstViewPosition();
	while(pos!=NULL)
	{
		pView=this->GetNextView(pos);
		if(pView->IsKindOf(RUNTIME_CLASS(CFormCommandView)))
			break;
	}
	CFormCommandView * view=(CFormCommandView*) pView;
	view->UpdateData();
    m_Isolevel=view->m_Isovalue;
	return;
}

void CMcubeDoc::OnMeshSkin()
{
	// TODO: Add your command handler code here
	//m_nMeshTriangleOne= new MeshTriangle;
	//m_nMeshTriangleOne->m_nPoint=m_CIsoSurface.m_nVertices;
	//m_nMeshTriangleOne->m_nTriangles=m_CIsoSurface.m_nTriangles;
	//m_nMeshTriangleOne->m_Normals= new float[m_nMeshTriangleOne->m_nPoint*3];
	//m_nMeshTriangleOne->m_Points=new float[m_nMeshTriangleOne->m_nPoint*3];
	//m_nMeshTriangleOne->m_Triangles=new int[m_nMeshTriangleOne->m_nTriangles*3];
	//for (int i=0;i<m_nMeshTriangleOne->m_nPoint;i++){
	//	m_nMeshTriangleOne->m_Points[i*3]=m_CIsoSurface.m_ppt3dVertices[i][0];
	//	m_nMeshTriangleOne->m_Points[i*3+1]=m_CIsoSurface.m_ppt3dVertices[i][1];
	//	m_nMeshTriangleOne->m_Points[i*3+2]=m_CIsoSurface.m_ppt3dVertices[i][2];
	//	m_nMeshTriangleOne->m_Normals[i*3]=m_CIsoSurface.m_pvec3dNormals[i][0];
	//	m_nMeshTriangleOne->m_Normals[i*3+1]=m_CIsoSurface.m_pvec3dNormals[i][1];
	//	m_nMeshTriangleOne->m_Normals[i*3+2]=m_CIsoSurface.m_pvec3dNormals[i][2];
 //  	}

	//for(int i=0;i<m_nMeshTriangleOne->m_nTriangles*3;i++){
	//	m_nMeshTriangleOne->m_Triangles[i]=m_CIsoSurface.m_piTriangleIndices[i];
 //       
	//}

	//return;
	//

}

void CMcubeDoc::OnMeshBone()
{
	// TODO: Add your command handler code here
	//m_nMeshTriangleTwo= new MeshTriangle;
	//m_nMeshTriangleTwo->m_nPoint=m_CIsoSurface.m_nVertices;
	//m_nMeshTriangleTwo->m_nTriangles=m_CIsoSurface.m_nTriangles;
	//m_nMeshTriangleTwo->m_Normals= new float[m_nMeshTriangleTwo->m_nPoint*3];
	//m_nMeshTriangleTwo->m_Points=new float[m_nMeshTriangleTwo->m_nPoint*3];
	//m_nMeshTriangleTwo->m_Triangles=new int[m_nMeshTriangleTwo->m_nTriangles*3];
	//for (int i=0;i<m_nMeshTriangleOne->m_nPoint;i++){
	//	m_nMeshTriangleTwo->m_Points[i*3]=m_CIsoSurface.m_ppt3dVertices[i][0];
	//	m_nMeshTriangleTwo->m_Points[i*3+1]=m_CIsoSurface.m_ppt3dVertices[i][1];
	//	m_nMeshTriangleTwo->m_Points[i*3+2]=m_CIsoSurface.m_ppt3dVertices[i][2];
	//	m_nMeshTriangleTwo->m_Normals[i*3]=m_CIsoSurface.m_pvec3dNormals[i][0];
	//	m_nMeshTriangleTwo->m_Normals[i*3+1]=m_CIsoSurface.m_pvec3dNormals[i][1];
	//	m_nMeshTriangleTwo->m_Normals[i*3+2]=m_CIsoSurface.m_pvec3dNormals[i][2];
	//}

	//for(int i=0;i<m_nMeshTriangleOne->m_nTriangles*3;i++){
	//	m_nMeshTriangleTwo->m_Triangles[i]=m_CIsoSurface.m_piTriangleIndices[i];

	//}

	//return;


}

void CMcubeDoc::OnFilterBilateral()
{
	// TODO: Add your command handler code here
 //   BilateralFilter m_bilateralFilter;
	//m_bilateralFilter.GetBilateralFilter(m_numOfPoints,m_pointSets,m_numOfTriangles,m_triangles,m_normals,m_colors);
	//for(int i=0;i<m_numOfPoints*3;i++){
	//	m_pointSets[i]=m_bilateralFilter.m_pointSets[i];
	//	m_colors[i]=m_bilateralFilter.m_colors[i];
	//}
	//m_bilateralFilter.DeleteBilateralFilter();
	//this->UpdateAllViews(NULL);
	//for(float www=0.20;www<=0.21;www+=0.05){
	//	m_noiseCoffecient=www;
	//AddPositionGaussianNoise();
	float m_meanlength=1.0;	
	//m_meanlength=CalculateMeanLength();

	//////////////////////////////////////////////////////////////////////////
	//非测试代码
	FilterBilateral1 m_bilateralFilter;
	m_bilateralFilter.GetFilterBilateral(m_numOfPoints,m_originalPointSets,m_colors,m_meanlength,m_variationF,m_variationG1,m_variationG2,m_variationH1,m_variationH2,m_variationH3,m_K,m_meanShift,m_iterativeTimes,m_indexFunction,m_variationNstop,m_functionVariation,m_thresholdColor,m_functionGradientWide,m_functionGradientThreshold,m_thresholdDistanceNormal,m_thresholdDistanceGradientFunction);
	//m_bilateralFilter.GetFilterBilateral(m_numOfPoints,m_originalPointSets,m_meanlength,m_variationF,m_variationF,m_variationG2,m_variationF,m_variationH2,m_variationH3,m_K,m_meanShift,m_iterativeTimes,m_indexFunction,m_variationNstop);
	if(m_colors!=NULL){
		delete[] m_colors;		
	}
	m_colors=new float[m_numOfPoints*3];
	
	for(int i=0;i<m_numOfPoints*3;i++){
		m_pointSets[i]=m_bilateralFilter.m_resultPointSet[i];	
		m_colors[i]=m_bilateralFilter.m_resultColors[i];
		m_normals[i]=m_bilateralFilter.m_originalNormals[i];
	}
	//delete[] m_normals;
	//int m_nNormals = m_numOfPoints;
	//m_normals=new float[m_nNormals*3];

	//// Set all normals to 0.
	//for (int i = 0; i < m_nNormals*3; i++) {
	//	m_normals[i]=0;
	//}

	//// Calculate normals.
	//for (int i = 0; i < m_numOfTriangles; i++) {
	//	VECTOR3D vec1, vec2, normal;
	//	int id0, id1, id2;
	//	id0 = m_triangles[i*3];
	//	id1 = m_triangles[i*3+1];
	//	id2 = m_triangles[i*3+2];
	//	vec1[0] = m_pointSets[id1*3]- m_pointSets[id0*3];
	//	vec1[1] = m_pointSets[id1*3+1] - m_pointSets[id0*3+1];
	//	vec1[2] = m_pointSets[id1*3+2] - m_pointSets[id0*3+2];
	//	vec2[0] = m_pointSets[id2*3] - m_pointSets[id0*3];
	//	vec2[1] = m_pointSets[id2*3+1] - m_pointSets[id0*3+1];
	//	vec2[2] = m_pointSets[id2*3+2] - m_pointSets[id0*3+2];
	//	normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
	//	normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
	//	normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
	//	m_normals[id0*3+0] += normal[0];
	//	m_normals[id0*3+1] += normal[1];
	//	m_normals[id0*3+2] += normal[2];
	//	m_normals[id1*3+0] += normal[0];
	//	m_normals[id1*3+1] += normal[1];
	//	m_normals[id1*3+2] += normal[2];
	//	m_normals[id2*3+0] += normal[0];
	//	m_normals[id2*3+1] += normal[1];
	//	m_normals[id2*3+2] += normal[2];
	//}

	//// Normalize normals.
	//for (i = 0; i < m_nNormals; i++) {
	//	float length = sqrt(m_normals[i*3+0]*m_normals[i*3+0] + m_normals[i*3+1]*m_normals[i*3+1] + m_normals[i*3+2]*m_normals[i*3+2]);
	//	m_normals[i*3] /= length;
	//	m_normals[i*3+1] /= length;
	//	m_normals[i*3+2] /= length;
	//}

//	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\time.txt";
//
//	FILE *fpout;
//
//	if((fpout = fopen(filename_pw, "w")) == NULL)
//	{
//		int dkjkd;
//		return;
//		//MessageBox("can't open the file!");
//	}
//	else
//	{
//		CString tempString1,tempString2;
//		tempString1="ifMeanShift           ";
//		tempString2="whichFunction         ";
//		fprintf(fpout,"%s %d %s %d\n",tempString1,m_meanShift,tempString2,m_indexFunction);
//		tempString1="m_numOFnearset        ";
//		tempString2="m_numOFiteract        ";
//		fprintf(fpout,"%s %d %s %d\n",tempString1,m_K,tempString2,m_iterativeTimes);
//
//		tempString1="m_numOfPoints         ";
//		tempString2="m_numOfTriangles      ";
//		fprintf(fpout,"%s %d %s %d\n",tempString1,m_numOfPoints,tempString2,m_numOfTriangles);
//		tempString1="searchKnearestTime    ";
//		tempString2="filterTime            ";
//		fprintf(fpout,"%s %d %s %d\n",tempString1,m_bilateralFilter.nearestTime,tempString2,m_bilateralFilter.filterTime);
//		tempString1="m_variationF          ";
//		tempString2="m_variationG1         ";
//		fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationF,tempString2,m_variationG1);
//		tempString1="m_variationG2         ";
//		tempString2="m_variationH1         ";
//		fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationG2,tempString2,m_variationH1);
//		tempString1="m_variationH2         ";
//		tempString2="m_variationH3         ";
//		fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationH2,tempString2,m_variationH3);
//
//		float tempMSE=0;
//		for(int i=0;i<m_numOfPoints;i++){
//			float* tempDistance;
//			std::map<int ,KnearestField>::iterator mapIterator=m_bilateralFilter.m_mapKnearest.find(i);
//			std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
//			tempDistance=new float[(*mapIterator).second.m_nearest.size()+1];
//			tempDistance[0]=((m_originalOriginalPointSets[i*3]-m_pointSets[i*3])*(m_originalOriginalPointSets[i*3]-m_pointSets[i*3])
//				+(m_originalOriginalPointSets[i*3+1]-m_pointSets[i*3+1])*(m_originalOriginalPointSets[i*3+1]-m_pointSets[i*3+1])
//				+(m_originalOriginalPointSets[i*3+2]-m_pointSets[i*3+2])*(m_originalOriginalPointSets[i*3+2]-m_pointSets[i*3+2])
//				);           
//			int k=1;
//			for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
//				tempDistance[k]=((m_originalOriginalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])*(m_originalOriginalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])
//					+(m_originalOriginalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])*(m_originalOriginalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])
//					+(m_originalOriginalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])*(m_originalOriginalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])
//					);
//				tempDistance[k]=((m_originalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])*(m_originalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])
//				+(m_originalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])*(m_originalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])
//				+(m_originalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])*(m_originalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])
//				);
//				if(tempDistance[k]<tempDistance[0]){
//					tempDistance[0]=tempDistance[k];
//				}
//				k++;
//			}
//			tempMSE+=tempDistance[0];
//			/*tempMSE+=((m_originalPointSets[i*3]-m_pointSets[i*3])*(m_originalPointSets[i*3]-m_pointSets[i*3])
//				+(m_originalPointSets[i*3+1]-m_pointSets[i*3+1])*(m_originalPointSets[i*3+1]-m_pointSets[i*3+1])
//				+(m_originalPointSets[i*3+2]-m_pointSets[i*3+2])*(m_originalPointSets[i*3+2]-m_pointSets[i*3+2])
//				);  */    
//			delete[] tempDistance;
//
//		}
//		tempMSE=sqrt(tempMSE/(m_numOfPoints*m_meanlength*3.14*0.52517*0.52517*m_meanlength));
//		tempString1="tempMSE               ";
//		tempString2="m_normalStop          ";
//		fprintf(fpout,"%s %f %s %f\n",tempString1,tempMSE,tempString2,m_variationNstop);
//		//fprintf(fpout,"%d %f %f %f\n",m_K,m_variationF,m_variationG1,tempMSE);
//		fclose(fpout);
//	}
////////////////////////////////////////////////////////////////////////////
//
//
//
//	//for(int wwww=6;wwww<45;wwww+=1){
//	////	for(float wwwww=0.15;wwwww<3.1;wwwww+=0.15){
//	//	//	for(float wwwwww=0.5;wwwwww<5.1;wwwwww+=0.5 ){
//	//			FilterBilateral1 m_bilateralFilter;
//	//			if(wwww==6)
//	//				m_variationF=1;
//	//			if(wwww>6&wwww<=18){
// //                   m_variationF=1.5+(double)(wwww-6)/24;
//	//			}
//	//			if(wwww>18&wwww<42)
//	//				m_variationF=2.5+(double)(wwww-18)/48;
//	//			if(wwww>=42)
//	//				m_variationF=3;
//	//			m_variationF/=2;
//	//			m_variationG1=m_variationF;
//	//			m_variationH1=m_variationF;
//	//			m_variationH2=0.6;
//	//			m_variationNstop=0.001;
//	//			m_indexFunction=7;
//	//			m_iterativeTimes=1;
//	//			m_meanShift=1;
//	//			//m_bilateralFilter.GetFilterBilateral(m_numOfPoints,m_originalPointSets,m_meanlength,m_variationF,m_variationG1,m_variationG2,m_variationH1,m_variationH2,m_variationH3,m_K,m_meanShift,m_iterativeTimes,m_indexFunction,m_variationNstop);
//	//			m_bilateralFilter.GetFilterBilateral(m_numOfPoints,m_originalPointSets,m_meanlength,m_variationF, m_variationG1,m_variationG2,m_variationH1,m_variationH2,m_variationH3,wwww,m_meanShift,m_iterativeTimes,m_indexFunction,m_variationNstop);
//	//			for(int i=0;i<m_numOfPoints*3;i++){
//	//				m_pointSets[i]=m_bilateralFilter.m_resultPointSet[i];	
//	//			}	
//
//
//
//	//			//////////////////////////////////////////////////////////////////////////
//	//			// 测试代码
//	//			CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\time.txt";
//
//	//			FILE *fpout;
//
//	//			if((fpout = fopen(filename_pw, "a")) == NULL)
//	//			{
//	//				int dkjkd;
//	//				//MessageBox("can't open the file!");
//	//			}
//	//			else
//	//			{
//	//				/*CString tempString1,tempString2;
//	//				tempString1="ifMeanShift           ";
//	//				tempString2="whichFunction         ";
//	//				fprintf(fpout,"%s %d %s %d\n",tempString1,m_meanShift,tempString2,m_indexFunction);
//	//				tempString1="m_numOFnearset        ";
//	//				tempString2="m_numOFiteract        ";
//	//				fprintf(fpout,"%s %d %s %d\n",tempString1,m_K,tempString2,m_iterativeTimes);
//
//	//				tempString1="m_numOfPoints         ";
//	//				tempString2="m_numOfTriangles      ";
//	//				fprintf(fpout,"%s %d %s %d\n",tempString1,m_numOfPoints,tempString2,m_numOfTriangles);
//	//				tempString1="searchKnearestTime    ";
//	//				tempString2="filterTime            ";
//	//				fprintf(fpout,"%s %d %s %d\n",tempString1,m_bilateralFilter.nearestTime,tempString2,m_bilateralFilter.filterTime);
//	//				tempString1="m_variationF          ";
//	//				tempString2="m_variationG1         ";
//	//				fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationF,tempString2,m_variationG1);
//	//				tempString1="m_variationG2         ";
//	//				tempString2="m_variationH1         ";
//	//				fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationG2,tempString2,m_variationH1);
//	//				tempString1="m_variationH2         ";
//	//				tempString2="m_variationH3         ";
//	//				fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationH2,tempString2,m_variationH3);*/
//
//	//				float tempMSE=0;
//	//				for(int i=0;i<m_numOfPoints;i++){
//	//					float* tempDistance;
//	//					std::map<int ,KnearestField>::iterator mapIterator=m_bilateralFilter.m_mapKnearest.find(i);
//	//					std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
//	//					tempDistance=new float[(*mapIterator).second.m_nearest.size()+1];
//	//					tempDistance[0]=((m_originalOriginalPointSets[i*3]-m_pointSets[i*3])*(m_originalOriginalPointSets[i*3]-m_pointSets[i*3])
//	//						+(m_originalOriginalPointSets[i*3+1]-m_pointSets[i*3+1])*(m_originalOriginalPointSets[i*3+1]-m_pointSets[i*3+1])
//	//						+(m_originalOriginalPointSets[i*3+2]-m_pointSets[i*3+2])*(m_originalOriginalPointSets[i*3+2]-m_pointSets[i*3+2])
//	//						);
//	//					int k=1;
//	//					for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
//	//						tempDistance[k]=((m_originalOriginalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])*(m_originalOriginalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])
//	//							+(m_originalOriginalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])*(m_originalOriginalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])
//	//							+(m_originalOriginalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])*(m_originalOriginalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])
//	//							);
//	//						/*tempDistance[k]=((m_originalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])*(m_originalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])
//	//						+(m_originalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])*(m_originalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])
//	//						+(m_originalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])*(m_originalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])
//	//						);*/
//	//						if(tempDistance[k]<tempDistance[0]){
//	//							tempDistance[0]=tempDistance[k];
//	//						}
//	//						k++;
//	//					}
//	//					tempMSE+=tempDistance[0];
//	//					delete[] tempDistance;
//
//	//				}
//	//				tempMSE=sqrt(tempMSE/(m_numOfPoints*m_meanlength*3.14*0.52517*0.52517*m_meanlength));
//	//				/*tempString1="tempMSE               ";
//	//				tempString2="m_normalStop          ";
//	//				fprintf(fpout,"%s %f %s %f\n",tempString1,tempMSE,tempString2,m_variationNstop);*/
//	//				fprintf(fpout,"%d %f\n",wwww,tempMSE);
//	//				fclose(fpout);
//	//			}
//
//	//			
//
//	//			//////////////////////////////////////////////////////////////////////////
//	//			m_bilateralFilter.DeleteFilterBilateral();
//
//	//		
//	//	}
//	//	}
//	//}
//	filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\point.txt";
//	if((fpout = fopen(filename_pw, "w")) == NULL)
//	{
//		int dkjkd;
//		//MessageBox("can't open the file!");
//	}
//	else
//	{
//
//		for(int i=0;i<m_numOfPoints;i++){
//			fprintf(fpout,"%f %f %f\n",m_pointSets[i*3],m_pointSets[i*3+1],m_pointSets[i*3+2]);		
//		}
//		
//		fclose(fpout);
//	}

    m_bilateralFilter.DeleteFilterBilateral();
	this->UpdateAllViews(NULL);

}

void CMcubeDoc::generateCube()
{
	/*m_pVol=new uint8[6859];
	
	for(int i=0;i<19;i++){
		for(int j=0;j<19;j++){
			int a, b;
			a=max(abs(i-9),abs(j-9));
			for(int k=0;k<19;k++){
				
				
				b=max(a,abs(k-9));

				m_pVol[i*19*19+j*19+k]=b;
			}
		}
	}
	return;*/
}

void CMcubeDoc::AddPositionGaussianNoise()
{
	//求邻近点的最小值

	//加入高斯噪声
	//////////////////////////////////////////////////////////////////////////
	// 首先加入 均值为1的 
	//m_resultPointSet=new float[m_numOfPoints*3];
	float m_meanLength=1.0;
	int m_nNormals=m_numOfPoints;
	m_meanLength=CalculateMeanLength();
	//float noise[3];
	//long idum[3];
	//idum[0]=1;
	//idum[1]=1;
	//idum[2]=1;

	//for(int i=0;i<m_numOfPoints;i++){
	//	for(int j=0;j<3;j++){
	//		noise[j]=gasdev(&idum[j]);
	//		//m_resultPointSet[i*3+j]=m_originalPointSet[i*3+j]+noise[j]*0.05;
	//		//m_pointSets[i*3+j]=m_pointSets[i*3+j]+noise[j]*0.15*m_meanLength;
	//		m_originalPointSets[i*3+j]=m_originalOriginalPointSets[i*3+j]+noise[j]*m_noiseCoffecient*m_meanLength;
	//		m_pointSets[i*3+j]=m_originalPointSets[i*3+j];
	//		//m_pointSets[i*3+j]=m_originalPointSets[i*3+j];
	//		//	m_resultPointSet[i*3+j]=m_originalPointSet[i*3+j];//+noise[j]*0.01;
	//	}        
	//}
	//return;

	//////////////////////////////////////////////////////////////////////////
	/*
	 *	加法向噪声
	 */
	//std::vector<Point> tempNormalsSets;
	float noise;
	long idum=1;
	//CString filename_pw="D:\\qin\\experiment\\normal.txt";
	//FILE *fpout;
	//if((fpout = fopen(filename_pw, "r")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	while(!feof(fpout)){
	//		Point tempNormals;
	//		float temp1,temp2;
	//		//fscanf(fpout,"%f %f %f %f %f",&tempPoint.vertexPosition[0],&tempPoint.vertexPosition[1],&tempPoint.vertexPosition[2],&temp1,&temp2);
	//		fscanf(fpout,"%f %f %f",&tempNormals.vertexPosition[0],&tempNormals.vertexPosition[1],&tempNormals.vertexPosition[2]);
	//		tempNormalsSets.push_back(tempNormals);
	//	}
	//	fclose(fpout);
	//}

	//temp_normals=new float[m_numOfPoints*3];
	//std::vector<Point>::iterator tempNormalIterator=tempNormalsSets.begin();
	//for(int i=0;tempNormalIterator!=tempNormalsSets.end();tempNormalIterator++){


	//	temp_normals[i*3]=(*tempNormalIterator).vertexPosition[0];//*1000;
	//	temp_normals[i*3+1]=(*tempNormalIterator).vertexPosition[1];//*1000;
	//	temp_normals[i*3+2]=(*tempNormalIterator).vertexPosition[2];//*1000;
	//}

	for(int i=0;i<m_numOfPoints;i++){
		noise=gasdev(&idum);
		for(int j=0;j<3;j++){			
			//m_resultPointSet[i*3+j]=m_originalPointSet[i*3+j]+noise[j]*0.05;
			//m_pointSets[i*3+j]=m_pointSets[i*3+j]+noise[j]*0.15*m_meanLength;
			m_originalPointSets[i*3+j]=m_originalOriginalPointSets[i*3+j]+noise*m_normals[i*3+j]*m_noiseCoffecient*m_meanLength;
			m_pointSets[i*3+j]=m_originalPointSets[i*3+j];
			//m_pointSets[i*3+j]=m_originalPointSets[i*3+j];
			//	m_resultPointSet[i*3+j]=m_originalPointSet[i*3+j];//+noise[j]*0.01;
		}        
	}

	// Set all normals to 0.
	for (int i = 0; i < m_nNormals*3; i++) {
		m_normals[i]=0;
	}

	// Calculate normals.
	for (int i = 0; i < m_numOfTriangles; i++) {
		VECTOR3D vec1, vec2, normal;
		int id0, id1, id2;
		id0 = m_triangles[i*3];
		id1 = m_triangles[i*3+1];
		id2 = m_triangles[i*3+2];
		vec1[0] = m_pointSets[id1*3]- m_pointSets[id0*3];
		vec1[1] = m_pointSets[id1*3+1] - m_pointSets[id0*3+1];
		vec1[2] = m_pointSets[id1*3+2] - m_pointSets[id0*3+2];
		vec2[0] = m_pointSets[id2*3] - m_pointSets[id0*3];
		vec2[1] = m_pointSets[id2*3+1] - m_pointSets[id0*3+1];
		vec2[2] = m_pointSets[id2*3+2] - m_pointSets[id0*3+2];
		normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
		normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
		normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
		m_normals[id0*3+0] += normal[0];
		m_normals[id0*3+1] += normal[1];
		m_normals[id0*3+2] += normal[2];
		m_normals[id1*3+0] += normal[0];
		m_normals[id1*3+1] += normal[1];
		m_normals[id1*3+2] += normal[2];
		m_normals[id2*3+0] += normal[0];
		m_normals[id2*3+1] += normal[1];
		m_normals[id2*3+2] += normal[2];
	}

	// Normalize normals.
	for (int i = 0; i < m_nNormals; i++) {
		float length = sqrt(m_normals[i*3+0]*m_normals[i*3+0] + m_normals[i*3+1]*m_normals[i*3+1] + m_normals[i*3+2]*m_normals[i*3+2]);
		m_normals[i*3] /= length;
		m_normals[i*3+1] /= length;
		m_normals[i*3+2] /= length;
	}
	return;
}

void CMcubeDoc::AddColorGaussianNorse()
{
	//把颜色转化为0-1之间的浮点数
	/*for(int i=0;i<m_numOfPoints*3;i++){
	m_originalFloatColor[i]=m_originalColor[i]/255;
	}
	delete [] m_originalColor;*/
	//加入高斯白噪声

	//m_resultColors=new float[m_numOfPoints*3];

	float noise[3];
	long idum[3];
	idum[0]=2;
	idum[1]=2;
	idum[2]=2;

	for(int i=0;i<m_numOfPoints;i++){
		for(int j=0;j<3;j++){
			noise[j]=gasdev(&idum[0]);
			m_colors[i*3+j]=m_originalColors[i*3+j]+noise[j]*m_noiseFuncCoff;
			/*m_originalColors[i*3+j]=m_originalColors[i*3+j]+noise[j]*0.05;
			m_colors[i*3+j]=m_originalColors[i*3+j];*/
			/*if(m_originalColor[i*3+j]<0){
			m_originalColor[i*3+j]=0;
			}
			else
			if(m_originalColor[i*3+j]>1){
			m_originalColor[i*3+j]=1;

			}*/
			//	m_resultColors[i*3+j]=m_colors[i*3+j];
		}        
	}
	//////////////////////////////////////////////////////////////////////////
	// 测试代码
	//CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\originalColors.txt";
	//FILE *fpout;
	//float averageColor[3],variationColor[3];
	//for(int j=0;j<3;j++){
	//	averageColor[j]=0;
	//	variationColor[j]=0;
	//}
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f\n",m_colors[i*3],m_colors[i*3+1],m_colors[i*3+2]);
	//		averageColor[0]+=m_colors[i*3];
	//		averageColor[1]+=m_colors[i*3+1];
	//		averageColor[2]+=m_colors[i*3+2];
	//	}
	//	for(int i=0;i<3;i++){
	//		averageColor[i]/=m_numOfPoints;
	//	}
	//	fprintf(fpout,"%f %f %f\n",averageColor[0],averageColor[1],averageColor[2]);
	//	for(int i=0;i<m_numOfPoints;i++){
	//		variationColor[0]+=(m_colors[i*3]-averageColor[0])*(m_colors[i*3]-averageColor[0]);
	//		variationColor[1]+=(m_colors[i*3+1]-averageColor[1])*(m_colors[i*3+1]-averageColor[1]);
	//		variationColor[2]+=(m_colors[i*3+2]-averageColor[2])*(m_colors[i*3+2]-averageColor[2]);
	//	}
	//	for(int i=0;i<3;i++){
	//		variationColor[i]/=(m_numOfPoints-1);
	//	}
	//	fprintf(fpout,"%f %f %f\n",variationColor[0],variationColor[1],variationColor[2]);

	//	fclose(fpout);
	//}
	//////////////////////////////////////////////////////////////////////////
}

void CMcubeDoc::GenerateTriangleFromImage()
{
	float* pointSet;
	pointSet=new float[512*512*2];
	for(int i=0;i<512;i++){
		for(int j=0;j<512;j++){
            pointSet[(i*512+j)*2]=i;
			pointSet[(i*512+j)*2+1]=j;
	
		}
	}
	//Delaunay类中的一些变量
	vertexSet vSetFront;
	Delaunay delaunay;
	triangleSet tSetFront;
	tIterator tIt;
	//std::vector<float[3]> tempTriangle;
	for(int i=0;i<512;i++){
		for(int j=0;j<512;j++){
			vSetFront.insert(vertex(pointSet[(i*512+j)*2],pointSet[(i*512+j)*2+1]));
		}
	}
//	delaunay.Triangulate(vSetFront,tSetFront);
	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\triangle.txt";

	FILE *fpout;
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else{
		for(tIt = tSetFront.begin();tIt!= tSetFront.end(); tIt++){
			int tempTriangle[3];
			for(int j=0;j<3;j++){	
				tempTriangle[j]=int((tIt->m_Vertices[j]->GetX())*512+tIt->m_Vertices[j]->GetY());				
			}
			fprintf(fpout,"%d %d %d\n",tempTriangle[0],tempTriangle[1],tempTriangle[2]);
		}
		fclose(fpout);
	}

	filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\point.txt";
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else{
		for(int i=0;i<512;i++){
			for(int j=0;j<512;j++){
                fprintf(fpout,"%f %f %f\n",(float)i,(float)j,0);
			}
		}
			
		fclose(fpout);
	}	
}


float CMcubeDoc::CalculateMeanLength()
{
	float meanLength=0;
	for(int i=0;i<m_numOfTriangles;i++){
		meanLength+=sqrt((m_originalPointSets[m_triangles[i*3]*3]-m_originalPointSets[m_triangles[i*3+1]*3])*(m_originalPointSets[m_triangles[i*3]*3]-m_originalPointSets[m_triangles[i*3+1]*3])
			+(m_originalPointSets[m_triangles[i*3]*3+1]-m_originalPointSets[m_triangles[i*3+1]*3+1])*(m_originalPointSets[m_triangles[i*3]*3+1]-m_originalPointSets[m_triangles[i*3+1]*3+1])
			+(m_originalPointSets[m_triangles[i*3]*3+2]-m_originalPointSets[m_triangles[i*3+1]*3+2])*(m_originalPointSets[m_triangles[i*3]*3+2]-m_originalPointSets[m_triangles[i*3+1]*3+2]));
		meanLength+=sqrt((m_originalPointSets[m_triangles[i*3]*3]-m_originalPointSets[m_triangles[i*3+2]*3])*(m_originalPointSets[m_triangles[i*3]*3]-m_originalPointSets[m_triangles[i*3+2]*3])
			+(m_originalPointSets[m_triangles[i*3]*3+1]-m_originalPointSets[m_triangles[i*3+2]*3+1])*(m_originalPointSets[m_triangles[i*3]*3+1]-m_originalPointSets[m_triangles[i*3+2]*3+1])
			+(m_originalPointSets[m_triangles[i*3]*3+2]-m_originalPointSets[m_triangles[i*3+2]*3+2])*(m_originalPointSets[m_triangles[i*3]*3+2]-m_originalPointSets[m_triangles[i*3+2]*3+2]));
		meanLength+=sqrt((m_originalPointSets[m_triangles[i*3+2]*3]-m_originalPointSets[m_triangles[i*3+1]*3])*(m_originalPointSets[m_triangles[i*3+2]*3]-m_originalPointSets[m_triangles[i*3+1]*3])
			+(m_originalPointSets[m_triangles[i*3+2]*3+1]-m_originalPointSets[m_triangles[i*3+1]*3+1])*(m_originalPointSets[m_triangles[i*3+2]*3+1]-m_originalPointSets[m_triangles[i*3+1]*3+1])
			+(m_originalPointSets[m_triangles[i*3+2]*3+2]-m_originalPointSets[m_triangles[i*3+1]*3+2])*(m_originalPointSets[m_triangles[i*3+2]*3+2]-m_originalPointSets[m_triangles[i*3+1]*3+2]));
	}

	meanLength/=(3*m_numOfTriangles);
	return meanLength;
}
void CMcubeDoc::OnFilterAddnoiseGeometry()
{
	// TODO: Add your command handler code here
	AddPositionGaussianNoise();
	//delete[] m_normals;

	//int m_nNormals = m_numOfPoints;
	//m_normals=new float[m_nNormals*3];

	//// Set all normals to 0.
	//for (int i = 0; i < m_nNormals*3; i++) {
	//	m_normals[i]=0;
	//}

	//// Calculate normals.
	//for (int i = 0; i < m_numOfTriangles; i++) {
	//	VECTOR3D vec1, vec2, normal;
	//	int id0, id1, id2;
	//	id0 = m_triangles[i*3];
	//	id1 = m_triangles[i*3+1];
	//	id2 = m_triangles[i*3+2];
	//	vec1[0] = m_pointSets[id1*3]- m_pointSets[id0*3];
	//	vec1[1] = m_pointSets[id1*3+1] - m_pointSets[id0*3+1];
	//	vec1[2] = m_pointSets[id1*3+2] - m_pointSets[id0*3+2];
	//	vec2[0] = m_pointSets[id2*3] - m_pointSets[id0*3];
	//	vec2[1] = m_pointSets[id2*3+1] - m_pointSets[id0*3+1];
	//	vec2[2] = m_pointSets[id2*3+2] - m_pointSets[id0*3+2];
	//	normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
	//	normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
	//	normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
	//	m_normals[id0*3+0] += normal[0];
	//	m_normals[id0*3+1] += normal[1];
	//	m_normals[id0*3+2] += normal[2];
	//	m_normals[id1*3+0] += normal[0];
	//	m_normals[id1*3+1] += normal[1];
	//	m_normals[id1*3+2] += normal[2];
	//	m_normals[id2*3+0] += normal[0];
	//	m_normals[id2*3+1] += normal[1];
	//	m_normals[id2*3+2] += normal[2];
	//}

	//// Normalize normals.
	//for (i = 0; i < m_nNormals; i++) {
	//	if(i==131074){
	//		int d;
	//		d=0;
	//	}
	//	float length = sqrt(m_normals[i*3+0]*m_normals[i*3+0] + m_normals[i*3+1]*m_normals[i*3+1] + m_normals[i*3+2]*m_normals[i*3+2]);
	//	m_normals[i*3] /= length;
	//	m_normals[i*3+1] /= length;
	//	m_normals[i*3+2] /= length;
	//}

	this->UpdateAllViews(NULL);

	
}

void CMcubeDoc::OnFilterAddnoiseFunction()
{
	// TODO: Add your command handler code here
	AddColorGaussianNorse();
}

void CMcubeDoc::OnFilterNonShift()
{
	// TODO: Add your command handler code here

	float m_meanlength=1.0;	
	m_meanlength=CalculateMeanLength();
	int m_nNormals=m_numOfPoints;

	//////////////////////////////////////////////////////////////////////////
	//非测试代码
	//m_variationF=m_variationF*m_meanlength;
	//m_variationG1=m_variationG1*m_meanlength;
	//m_variationG2=m_variationG2*m_meanlength;
	//m_variationH1=m_variationH1*m_meanlength;
	m_dimensionColor=3;
	
	
	NonShiftBilateralFilter m_nonShiftBilateralFilter;
	m_nonShiftBilateralFilter.GetNonShiftBilateralFilter(m_numOfPoints,m_originalPointSets,m_dimensionColor,m_colors,m_K,m_meanShift,m_ifVolumePreserve,m_indexFunction,m_variationNstop,m_variationF,m_variationG1,m_iterativeTimes,m_variationG2,m_variationH1,m_variationH2,m_thresholdColor,m_functionGradientThreshold,m_thresholdDistanceNormal,m_radius*m_meanlength,m_timeStep,m_tangentOrManifold,m_ifNormalWeight,m_ifAreaWeight,m_ifVariationNormal);
	
//	m_bilateralFilter.GetFilterBilateral(m_numOfPoints,m_originalPointSets,m_colors,m_meanlength,m_variationF,m_variationG1,m_variationG2,m_variationH1,m_variationH2,m_variationH3,m_K,m_meanShift,m_iterativeTimes,m_indexFunction,m_variationNstop,m_functionVariation,m_thresholdColor,m_functionGradientWide,m_functionGradientThreshold);
	//m_bilateralFilter.GetFilterBilateral(m_numOfPoints,m_originalPointSets,m_meanlength,m_variationF,m_variationF,m_variationG2,m_variationF,m_variationH2,m_variationH3,m_K,m_meanShift,m_iterativeTimes,m_indexFunction,m_variationNstop);
	for(int i=0;i<m_numOfPoints*3;i++){
		m_pointSets[i]=(float)m_nonShiftBilateralFilter.m_resultPointSet[i];	
	    m_colors[i]=m_nonShiftBilateralFilter.m_resultColor[i];
		m_normals[i]=-m_nonShiftBilateralFilter.m_originalNormals[i];
	}


	// Set all normals to 0.
	for (int i = 0; i < m_nNormals*3; i++) {
		m_normals[i]=0;
	}

	// Calculate normals.
	for (int i = 0; i < m_numOfTriangles; i++) {
		VECTOR3D vec1, vec2, normal;
		int id0, id1, id2;
		id0 = m_triangles[i*3];
		id1 = m_triangles[i*3+1];
		id2 = m_triangles[i*3+2];
		vec1[0] = m_pointSets[id1*3]- m_pointSets[id0*3];
		vec1[1] = m_pointSets[id1*3+1] - m_pointSets[id0*3+1];
		vec1[2] = m_pointSets[id1*3+2] - m_pointSets[id0*3+2];
		vec2[0] = m_pointSets[id2*3] - m_pointSets[id0*3];
		vec2[1] = m_pointSets[id2*3+1] - m_pointSets[id0*3+1];
		vec2[2] = m_pointSets[id2*3+2] - m_pointSets[id0*3+2];
		normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
		normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
		normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
		m_normals[id0*3+0] += normal[0];
		m_normals[id0*3+1] += normal[1];
		m_normals[id0*3+2] += normal[2];
		m_normals[id1*3+0] += normal[0];
		m_normals[id1*3+1] += normal[1];
		m_normals[id1*3+2] += normal[2];
		m_normals[id2*3+0] += normal[0];
		m_normals[id2*3+1] += normal[1];
		m_normals[id2*3+2] += normal[2];
	}

	// Normalize normals.
	for (int i = 0; i < m_nNormals; i++) {
		float length = sqrt(m_normals[i*3+0]*m_normals[i*3+0] + m_normals[i*3+1]*m_normals[i*3+1] + m_normals[i*3+2]*m_normals[i*3+2]);
		m_normals[i*3] /= length;
		m_normals[i*3+1] /= length;
		m_normals[i*3+2] /= length;
	}

	CString filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonShifttime.txt";

	FILE *fpout;

	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	return;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	CString tempString1,tempString2;
	//	tempString1="ifMeanShift           ";
	//	tempString2="whichFunction         ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_meanShift,tempString2,m_indexFunction);
	//	tempString1="m_numOFnearset        ";
	//	tempString2="m_numOFiteract        ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_K,tempString2,m_iterativeTimes);

	//	tempString1="m_numOfPoints         ";
	//	tempString2="m_numOfTriangles      ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_numOfPoints,tempString2,m_numOfTriangles);
	//	tempString1="searchKnearestTime    ";
	//	tempString2="filterTime            ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_nonShiftBilateralFilter.nearestTime,tempString2,m_nonShiftBilateralFilter.filterTime);
	//	tempString1="m_variationF          ";
	//	tempString2="m_variationG1         ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationF,tempString2,m_variationG1);
	//	tempString1="m_variationG2         ";
	//	tempString2="m_variationH1         ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationG2,tempString2,m_variationH1);
	//	tempString1="m_variationH2         ";
	//	tempString2="m_variationH3         ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationH2,tempString2,m_variationH3);

	//	float tempMSE=0;
	//	//if(m_colors!=NULL)
	//	//	delete[] m_colors;
	//	//m_colors=new float[m_numOfPoints*3];
	//	for(int i=0;i<m_numOfPoints;i++){
	//		float* tempDistance;
	//		std::map<int ,KnearestField>::iterator mapIterator=m_nonShiftBilateralFilter.m_mapKnearest.find(i);
	//		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	//		tempDistance=new float[(*mapIterator).second.m_nearest.size()+1];
	//		tempDistance[0]=((m_originalOriginalPointSets[i*3]-m_pointSets[i*3])*(m_originalOriginalPointSets[i*3]-m_pointSets[i*3])
	//			+(m_originalOriginalPointSets[i*3+1]-m_pointSets[i*3+1])*(m_originalOriginalPointSets[i*3+1]-m_pointSets[i*3+1])
	//			+(m_originalOriginalPointSets[i*3+2]-m_pointSets[i*3+2])*(m_originalOriginalPointSets[i*3+2]-m_pointSets[i*3+2])
	//			);           
	//		int k=1;
	//		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
	//			tempDistance[k]=((m_originalOriginalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])*(m_originalOriginalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])
	//				+(m_originalOriginalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])*(m_originalOriginalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])
	//				+(m_originalOriginalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])*(m_originalOriginalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])
	//				);
	//			/*tempDistance[k]=((m_originalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])*(m_originalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])
	//				+(m_originalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])*(m_originalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])
	//				+(m_originalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])*(m_originalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])
	//				);*/
	//			if(tempDistance[k]<tempDistance[0]){
	//				tempDistance[0]=tempDistance[k];
	//			}
	//			k++;
	//		}

	///*		if(tempDistance[0]<=0.05){
	//			m_colors[i*3]=0;
	//			m_colors[i*3+1]=0;
	//			m_colors[i*3+2]=0.5+4*tempDistance[0];
	//		}
	//		if(tempDistance[0]>0.05&tempDistance[0]<=0.15){
	//			m_colors[i*3]=0;
	//			m_colors[i*3+1]=0.5+4*(tempDistance[0]-0.1);
	//			m_colors[i*3+2]=0;
	//		}
	//		if(tempDistance[0]>0.15){
	//			m_colors[i*3]=0.3+0.25*(tempDistance[0]-0.2);
	//			m_colors[i*3+1]=0;
	//			m_colors[i*3+2]=0;
	//		}*/
	//		


	//				
	//		tempMSE+=tempDistance[0];
	//		/*tempMSE+=((m_originalPointSets[i*3]-m_pointSets[i*3])*(m_originalPointSets[i*3]-m_pointSets[i*3])
	//		+(m_originalPointSets[i*3+1]-m_pointSets[i*3+1])*(m_originalPointSets[i*3+1]-m_pointSets[i*3+1])
	//		+(m_originalPointSets[i*3+2]-m_pointSets[i*3+2])*(m_originalPointSets[i*3+2]-m_pointSets[i*3+2])
	//		);  */    
	//		delete[] tempDistance;

	//	}
	//
	//	tempMSE=sqrt(tempMSE/(m_numOfPoints*m_meanlength*3.14*0.52517*0.52517*m_meanlength));
	//	tempString1="tempMSE               ";
	//	tempString2="m_normalStop          ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,tempMSE,tempString2,m_variationNstop);
	//	//fprintf(fpout,"%d %f %f %f\n",m_K,m_variationF,m_variationG1,tempMSE);
	//	fclose(fpout);
	//}  
	filename_pw = "D:\\qin\\experiment\\originalPoint.txt";
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{

		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f\n",m_pointSets[i*3],m_pointSets[i*3+1],m_pointSets[i*3+2]);		
		}
		
		fclose(fpout);
	}

	m_nonShiftBilateralFilter.DeleteNonShiftBilateralFilter();
	ComputeError();
	this->UpdateAllViews(NULL);
}

void CMcubeDoc::OnMeshlessFilter()
{
	// TODO: Add your command handler code here
	float m_meanlength;	
	m_meanlength=CalculateMeanLength();

	//////////////////////////////////////////////////////////////////////////
	//非测试代码
	//m_variationF=m_variationF*m_meanlength;
	//m_variationG1=m_variationG1*m_meanlength;
	//m_variationG2=m_variationG2*m_meanlength;
	//m_variationH1=m_variationH1*m_meanlength;
	m_dimensionColor=3;
    //MeshlessFilter m_meshLessFilter;
	
	//OneFreedomMeshless m_meshLessFilter;
	//m_meshLessFilter.GetMeshlessFilter(m_numOfPoints,m_originalPointSets,m_K,m_radius*m_meanlength,m_meanShift,m_iterativeTimes,0,0,m_timeStep,m_loadConstant,m_positionError*m_meanlength);

	RBFmeshless m_meshLessFilter;
	/*for(int i=0;i<m_numOfPoints*3;i++){
		m_originalPointSets[i]=m_originalPointSets[i]/100000;
	}*/
	m_meshLessFilter.GetMeshlessFilter(m_numOfPoints,m_originalPointSets,m_K,m_radius*m_meanlength,m_iterativeTimes,m_timeStep,m_loadConstant,m_positionError*m_meanlength,m_stopN,m_maxIter);
	//for(int i=0;i<m_numOfPoints*3;i++){
	//	m_originalPointSets[i]=m_originalPointSets[i]*100000;
	//}
	for(int i=0;i<m_numOfPoints*3;i++){
		m_pointSets[i]=(float)m_meshLessFilter.m_resultPointSet[i];	
		//m_colors[i]=m_nonShiftBilateralFilter.m_resultColor[i];
	}
	//if(m_visualization==1){
	//	if(m_colors!=NULL)
	//		delete m_colors;
	//	m_colors=new float[m_numOfPoints*3];
	//	for(int i=0;i<m_numOfPoints*3;i++){
	//		m_colors[i]=m_nonShiftBilateralFilter.m_resultColor[i];			
	//	}
	//}

	delete[] m_normals;
	int m_nNormals = m_numOfPoints;
	m_normals=new float[m_nNormals*3];

	// Set all normals to 0.
	for (int i = 0; i < m_nNormals*3; i++) {
		m_normals[i]=0;
	}

	// Calculate normals.
	for (int i = 0; i < m_numOfTriangles; i++) {
		VECTOR3D vec1, vec2, normal;
		int id0, id1, id2;
		id0 = m_triangles[i*3];
		id1 = m_triangles[i*3+1];
		id2 = m_triangles[i*3+2];
		vec1[0] = m_pointSets[id1*3]- m_pointSets[id0*3];
		vec1[1] = m_pointSets[id1*3+1] - m_pointSets[id0*3+1];
		vec1[2] = m_pointSets[id1*3+2] - m_pointSets[id0*3+2];
		vec2[0] = m_pointSets[id2*3] - m_pointSets[id0*3];
		vec2[1] = m_pointSets[id2*3+1] - m_pointSets[id0*3+1];
		vec2[2] = m_pointSets[id2*3+2] - m_pointSets[id0*3+2];
		normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
		normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
		normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
		m_normals[id0*3+0] += normal[0];
		m_normals[id0*3+1] += normal[1];
		m_normals[id0*3+2] += normal[2];
		m_normals[id1*3+0] += normal[0];
		m_normals[id1*3+1] += normal[1];
		m_normals[id1*3+2] += normal[2];
		m_normals[id2*3+0] += normal[0];
		m_normals[id2*3+1] += normal[1];
		m_normals[id2*3+2] += normal[2];
	}

	// Normalize normals.
	for (int i = 0; i < m_nNormals; i++) {
		float length = sqrt(m_normals[i*3+0]*m_normals[i*3+0] + m_normals[i*3+1]*m_normals[i*3+1] + m_normals[i*3+2]*m_normals[i*3+2]);
		m_normals[i*3] /= length;
		m_normals[i*3+1] /= length;
		m_normals[i*3+2] /= length;
	}



	CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\resultPoints.txt";

	FILE *fpout;
	if((fpout=fopen(filename_pw,"w"))==NULL){
		int dkjkd;

		return;
	}
	else{
		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f\n",m_pointSets[i*3],m_pointSets[i*3+1],m_pointSets[i*3+2]);		

		}
		fclose(fpout);
        
	}

	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	return;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	CString tempString1,tempString2;
	//	tempString1="ifMeanShift           ";
	//	tempString2="whichFunction         ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_meanShift,tempString2,m_indexFunction);
	//	tempString1="m_numOFnearset        ";
	//	tempString2="m_numOFiteract        ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_K,tempString2,m_iterativeTimes);

	//	tempString1="m_numOfPoints         ";
	//	tempString2="m_numOfTriangles      ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_numOfPoints,tempString2,m_numOfTriangles);
	//	tempString1="searchKnearestTime    ";
	//	tempString2="filterTime            ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_nonShiftBilateralFilter.nearestTime,tempString2,m_nonShiftBilateralFilter.filterTime);
	//	tempString1="m_variationF          ";
	//	tempString2="m_variationG1         ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationF,tempString2,m_variationG1);
	//	tempString1="m_variationG2         ";
	//	tempString2="m_variationH1         ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationG2,tempString2,m_variationH1);
	//	tempString1="m_variationH2         ";
	//	tempString2="m_variationH3         ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationH2,tempString2,m_variationH3);

	//	float tempMSE=0;
	//	//if(m_colors!=NULL)
	//	//	delete[] m_colors;
	//	//m_colors=new float[m_numOfPoints*3];
	//	for(int i=0;i<m_numOfPoints;i++){
	//		float* tempDistance;
	//		std::map<int ,KnearestField>::iterator mapIterator=m_nonShiftBilateralFilter.m_mapKnearest.find(i);
	//		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	//		tempDistance=new float[(*mapIterator).second.m_nearest.size()+1];
	//		tempDistance[0]=((m_originalOriginalPointSets[i*3]-m_pointSets[i*3])*(m_originalOriginalPointSets[i*3]-m_pointSets[i*3])
	//			+(m_originalOriginalPointSets[i*3+1]-m_pointSets[i*3+1])*(m_originalOriginalPointSets[i*3+1]-m_pointSets[i*3+1])
	//			+(m_originalOriginalPointSets[i*3+2]-m_pointSets[i*3+2])*(m_originalOriginalPointSets[i*3+2]-m_pointSets[i*3+2])
	//			);           
	//		int k=1;
	//		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
	//			tempDistance[k]=((m_originalOriginalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])*(m_originalOriginalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])
	//				+(m_originalOriginalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])*(m_originalOriginalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])
	//				+(m_originalOriginalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])*(m_originalOriginalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])
	//				);
	//			/*tempDistance[k]=((m_originalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])*(m_originalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])
	//			+(m_originalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])*(m_originalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])
	//			+(m_originalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])*(m_originalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])
	//			);*/
	//			if(tempDistance[k]<tempDistance[0]){
	//				tempDistance[0]=tempDistance[k];
	//			}
	//			k++;
	//		}

	//		/*		if(tempDistance[0]<=0.05){
	//		m_colors[i*3]=0;
	//		m_colors[i*3+1]=0;
	//		m_colors[i*3+2]=0.5+4*tempDistance[0];
	//		}
	//		if(tempDistance[0]>0.05&tempDistance[0]<=0.15){
	//		m_colors[i*3]=0;
	//		m_colors[i*3+1]=0.5+4*(tempDistance[0]-0.1);
	//		m_colors[i*3+2]=0;
	//		}
	//		if(tempDistance[0]>0.15){
	//		m_colors[i*3]=0.3+0.25*(tempDistance[0]-0.2);
	//		m_colors[i*3+1]=0;
	//		m_colors[i*3+2]=0;
	//		}*/




	//		tempMSE+=tempDistance[0];
	//		/*tempMSE+=((m_originalPointSets[i*3]-m_pointSets[i*3])*(m_originalPointSets[i*3]-m_pointSets[i*3])
	//		+(m_originalPointSets[i*3+1]-m_pointSets[i*3+1])*(m_originalPointSets[i*3+1]-m_pointSets[i*3+1])
	//		+(m_originalPointSets[i*3+2]-m_pointSets[i*3+2])*(m_originalPointSets[i*3+2]-m_pointSets[i*3+2])
	//		);  */    
	//		delete[] tempDistance;

	//	}

	//	tempMSE=sqrt(tempMSE/(m_numOfPoints*m_meanlength*3.14*0.52517*0.52517*m_meanlength));
	//	tempString1="tempMSE               ";
	//	tempString2="m_normalStop          ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,tempMSE,tempString2,m_variationNstop);
	//	//fprintf(fpout,"%d %f %f %f\n",m_K,m_variationF,m_variationG1,tempMSE);
	//	fclose(fpout);
	//}  
	//filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonShiftPoint.txt";
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{

	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f\n",m_pointSets[i*3],m_pointSets[i*3+1],m_pointSets[i*3+2]);		
	//	}

	//	fclose(fpout);
	//}

	m_meshLessFilter.DeleteMeshlessFilter();
	this->UpdateAllViews(NULL);
}

void CMcubeDoc::OnFilterAnisoMeshless()
{
	// TODO: Add your command handler code here
	// TODO: Add your command handler code here
	float m_meanlength;	
	m_meanlength=CalculateMeanLength();
	m_dimensionColor=3;
	AnisoRBFmeshless m_meshLessFilter;
	//for(int i=0;i<m_numOfPoints*3;i++){
	//	m_originalPointSets[i]=m_originalPointSets[i]/100000;
	//}
   
	m_meshLessFilter.GetMeshlessFilter(m_numOfPoints,m_originalPointSets,m_K,m_radius*m_meanlength,m_iterativeTimes,m_timeStep,m_loadConstant,m_positionError*m_meanlength,m_stopN,m_maxIter,m_curveThreshold);
	//for(int i=0;i<m_numOfPoints*3;i++){
	//	m_originalPointSets[i]=m_originalPointSets[i]*100000;
	//}
	for(int i=0;i<m_numOfPoints*3;i++){
		m_pointSets[i]=(float)m_meshLessFilter.m_resultPointSet[i];	
		m_colors[i]=m_meshLessFilter.m_resultColor[i];
	}
	//if(m_visualization==1){
	//	if(m_colors!=NULL)
	//		delete m_colors;
	//	m_colors=new float[m_numOfPoints*3];
	//	for(int i=0;i<m_numOfPoints*3;i++){
	//		m_colors[i]=m_nonShiftBilateralFilter.m_resultColor[i];			
	//	}
	//}

	delete[] m_normals;
	int m_nNormals = m_numOfPoints;
	m_normals=new float[m_nNormals*3];

	// Set all normals to 0.
	for (int i = 0; i < m_nNormals*3; i++) {
		m_normals[i]=0;
	}

	// Calculate normals.
	for (int i = 0; i < m_numOfTriangles; i++) {
		VECTOR3D vec1, vec2, normal;
		int id0, id1, id2;
		id0 = m_triangles[i*3];
		id1 = m_triangles[i*3+1];
		id2 = m_triangles[i*3+2];
		vec1[0] = m_pointSets[id1*3]- m_pointSets[id0*3];
		vec1[1] = m_pointSets[id1*3+1] - m_pointSets[id0*3+1];
		vec1[2] = m_pointSets[id1*3+2] - m_pointSets[id0*3+2];
		vec2[0] = m_pointSets[id2*3] - m_pointSets[id0*3];
		vec2[1] = m_pointSets[id2*3+1] - m_pointSets[id0*3+1];
		vec2[2] = m_pointSets[id2*3+2] - m_pointSets[id0*3+2];
		normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
		normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
		normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
		m_normals[id0*3+0] += normal[0];
		m_normals[id0*3+1] += normal[1];
		m_normals[id0*3+2] += normal[2];
		m_normals[id1*3+0] += normal[0];
		m_normals[id1*3+1] += normal[1];
		m_normals[id1*3+2] += normal[2];
		m_normals[id2*3+0] += normal[0];
		m_normals[id2*3+1] += normal[1];
		m_normals[id2*3+2] += normal[2];
	}

	// Normalize normals.
	for (int i = 0; i < m_nNormals; i++) {
		float length = sqrt(m_normals[i*3+0]*m_normals[i*3+0] + m_normals[i*3+1]*m_normals[i*3+1] + m_normals[i*3+2]*m_normals[i*3+2]);
		m_normals[i*3] /= length;
		m_normals[i*3+1] /= length;
		m_normals[i*3+2] /= length;
	}



	CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\time.txt";

	FILE *fpout;
	if((fpout=fopen(filename_pw,"w"))==NULL){
		int dkjkd;
		return;
	}
	else{
		
			fprintf(fpout,"%f %f %f\n",(float)(m_meshLessFilter.nearestTime/1000/m_iterativeTimes),(float)(m_meshLessFilter.filterTime/1000/m_iterativeTimes),m_meanlength);		

		
		fclose(fpout);

	}

	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	return;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	CString tempString1,tempString2;
	//	tempString1="ifMeanShift           ";
	//	tempString2="whichFunction         ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_meanShift,tempString2,m_indexFunction);
	//	tempString1="m_numOFnearset        ";
	//	tempString2="m_numOFiteract        ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_K,tempString2,m_iterativeTimes);

	//	tempString1="m_numOfPoints         ";
	//	tempString2="m_numOfTriangles      ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_numOfPoints,tempString2,m_numOfTriangles);
	//	tempString1="searchKnearestTime    ";
	//	tempString2="filterTime            ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_nonShiftBilateralFilter.nearestTime,tempString2,m_nonShiftBilateralFilter.filterTime);
	//	tempString1="m_variationF          ";
	//	tempString2="m_variationG1         ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationF,tempString2,m_variationG1);
	//	tempString1="m_variationG2         ";
	//	tempString2="m_variationH1         ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationG2,tempString2,m_variationH1);
	//	tempString1="m_variationH2         ";
	//	tempString2="m_variationH3         ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationH2,tempString2,m_variationH3);

	//	float tempMSE=0;
	//	//if(m_colors!=NULL)
	//	//	delete[] m_colors;
	//	//m_colors=new float[m_numOfPoints*3];
	//	for(int i=0;i<m_numOfPoints;i++){
	//		float* tempDistance;
	//		std::map<int ,KnearestField>::iterator mapIterator=m_nonShiftBilateralFilter.m_mapKnearest.find(i);
	//		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	//		tempDistance=new float[(*mapIterator).second.m_nearest.size()+1];
	//		tempDistance[0]=((m_originalOriginalPointSets[i*3]-m_pointSets[i*3])*(m_originalOriginalPointSets[i*3]-m_pointSets[i*3])
	//			+(m_originalOriginalPointSets[i*3+1]-m_pointSets[i*3+1])*(m_originalOriginalPointSets[i*3+1]-m_pointSets[i*3+1])
	//			+(m_originalOriginalPointSets[i*3+2]-m_pointSets[i*3+2])*(m_originalOriginalPointSets[i*3+2]-m_pointSets[i*3+2])
	//			);           
	//		int k=1;
	//		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
	//			tempDistance[k]=((m_originalOriginalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])*(m_originalOriginalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])
	//				+(m_originalOriginalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])*(m_originalOriginalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])
	//				+(m_originalOriginalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])*(m_originalOriginalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])
	//				);
	//			/*tempDistance[k]=((m_originalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])*(m_originalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])
	//			+(m_originalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])*(m_originalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])
	//			+(m_originalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])*(m_originalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])
	//			);*/
	//			if(tempDistance[k]<tempDistance[0]){
	//				tempDistance[0]=tempDistance[k];
	//			}
	//			k++;
	//		}

	//		/*		if(tempDistance[0]<=0.05){
	//		m_colors[i*3]=0;
	//		m_colors[i*3+1]=0;
	//		m_colors[i*3+2]=0.5+4*tempDistance[0];
	//		}
	//		if(tempDistance[0]>0.05&tempDistance[0]<=0.15){
	//		m_colors[i*3]=0;
	//		m_colors[i*3+1]=0.5+4*(tempDistance[0]-0.1);
	//		m_colors[i*3+2]=0;
	//		}
	//		if(tempDistance[0]>0.15){
	//		m_colors[i*3]=0.3+0.25*(tempDistance[0]-0.2);
	//		m_colors[i*3+1]=0;
	//		m_colors[i*3+2]=0;
	//		}*/




	//		tempMSE+=tempDistance[0];
	//		/*tempMSE+=((m_originalPointSets[i*3]-m_pointSets[i*3])*(m_originalPointSets[i*3]-m_pointSets[i*3])
	//		+(m_originalPointSets[i*3+1]-m_pointSets[i*3+1])*(m_originalPointSets[i*3+1]-m_pointSets[i*3+1])
	//		+(m_originalPointSets[i*3+2]-m_pointSets[i*3+2])*(m_originalPointSets[i*3+2]-m_pointSets[i*3+2])
	//		);  */    
	//		delete[] tempDistance;

	//	}

	//	tempMSE=sqrt(tempMSE/(m_numOfPoints*m_meanlength*3.14*0.52517*0.52517*m_meanlength));
	//	tempString1="tempMSE               ";
	//	tempString2="m_normalStop          ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,tempMSE,tempString2,m_variationNstop);
	//	//fprintf(fpout,"%d %f %f %f\n",m_K,m_variationF,m_variationG1,tempMSE);
	//	fclose(fpout);
	//}  
	filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\resultPoint.txt";
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{

		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f\n",m_pointSets[i*3],m_pointSets[i*3+1],m_pointSets[i*3+2]);		
		}

		fclose(fpout);
	}

	m_meshLessFilter.DeleteMeshlessFilter();
	this->UpdateAllViews(NULL);
}

void CMcubeDoc::OnFilterMls()
{
	// TODO: Add your command handler code here
	float m_meanlength=1;	
	//m_meanlength=CalculateMeanLength();
	m_dimensionColor=3;
	MLSfilter m_MLSfilter;
	//for(int i=0;i<m_numOfPoints*3;i++){
	//	m_originalPointSets[i]=m_originalPointSets[i]/100000;
	//}
    m_MLSfilter.GetMLSfilter(m_numOfPoints,m_originalPointSets,m_K,m_radius*m_meanlength,m_iterativeTimes);
	//m_meshLessFilter.GetMeshlessFilter(m_numOfPoints,m_originalPointSets,m_K,m_radius*m_meanlength,m_iterativeTimes,m_timeStep,m_loadConstant,m_positionError*m_meanlength,m_stopN,m_maxIter,m_curveThreshold);
	//for(int i=0;i<m_numOfPoints*3;i++){
	//	m_originalPointSets[i]=m_originalPointSets[i]*100000;
	//}
	if(m_colors!=NULL){
		delete[] m_colors;		
	}
	m_colors=new float[m_numOfPoints*3];
	for(int i=0;i<m_numOfPoints*3;i++){
		m_pointSets[i]=(float)m_MLSfilter.m_resultPointSet[i];	
		m_colors[i]=(float)m_MLSfilter.m_resultColor[i];
		m_normals[i]=(float)m_MLSfilter.m_originalNormals[i];
	}
	//if(m_visualization==1){
	//	if(m_colors!=NULL)
	//		delete m_colors;
	//	m_colors=new float[m_numOfPoints*3];
	//	for(int i=0;i<m_numOfPoints*3;i++){
	//		m_colors[i]=m_nonShiftBilateralFilter.m_resultColor[i];			
	//	}
	//}

	//delete[] m_normals;
	//int m_nNormals = m_numOfPoints;
	//m_normals=new float[m_nNormals*3];

	//// Set all normals to 0.
	//for (int i = 0; i < m_nNormals*3; i++) {
	//	m_normals[i]=0;
	//}

	//// Calculate normals.
	//for (int i = 0; i < m_numOfTriangles; i++) {
	//	VECTOR3D vec1, vec2, normal;
	//	int id0, id1, id2;
	//	id0 = m_triangles[i*3];
	//	id1 = m_triangles[i*3+1];
	//	id2 = m_triangles[i*3+2];
	//	vec1[0] = m_pointSets[id1*3]- m_pointSets[id0*3];
	//	vec1[1] = m_pointSets[id1*3+1] - m_pointSets[id0*3+1];
	//	vec1[2] = m_pointSets[id1*3+2] - m_pointSets[id0*3+2];
	//	vec2[0] = m_pointSets[id2*3] - m_pointSets[id0*3];
	//	vec2[1] = m_pointSets[id2*3+1] - m_pointSets[id0*3+1];
	//	vec2[2] = m_pointSets[id2*3+2] - m_pointSets[id0*3+2];
	//	normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
	//	normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
	//	normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
	//	m_normals[id0*3+0] += normal[0];
	//	m_normals[id0*3+1] += normal[1];
	//	m_normals[id0*3+2] += normal[2];
	//	m_normals[id1*3+0] += normal[0];
	//	m_normals[id1*3+1] += normal[1];
	//	m_normals[id1*3+2] += normal[2];
	//	m_normals[id2*3+0] += normal[0];
	//	m_normals[id2*3+1] += normal[1];
	//	m_normals[id2*3+2] += normal[2];
	//}

	//// Normalize normals.
	//for (i = 0; i < m_nNormals; i++) {
	//	float length = sqrt(m_normals[i*3+0]*m_normals[i*3+0] + m_normals[i*3+1]*m_normals[i*3+1] + m_normals[i*3+2]*m_normals[i*3+2]);
	//	m_normals[i*3] /= length;
	//	m_normals[i*3+1] /= length;
	//	m_normals[i*3+2] /= length;
	//}



	/*CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\time.txt";

	FILE *fpout;
	if((fpout=fopen(filename_pw,"w"))==NULL){
		int dkjkd;
		return;
	}
	else{

		fprintf(fpout,"%f %f\n",m_meshLessFilter.nearestTime,m_meshLessFilter.filterTime);		


		fclose(fpout);

	}*/

	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	return;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	CString tempString1,tempString2;
	//	tempString1="ifMeanShift           ";
	//	tempString2="whichFunction         ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_meanShift,tempString2,m_indexFunction);
	//	tempString1="m_numOFnearset        ";
	//	tempString2="m_numOFiteract        ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_K,tempString2,m_iterativeTimes);

	//	tempString1="m_numOfPoints         ";
	//	tempString2="m_numOfTriangles      ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_numOfPoints,tempString2,m_numOfTriangles);
	//	tempString1="searchKnearestTime    ";
	//	tempString2="filterTime            ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_nonShiftBilateralFilter.nearestTime,tempString2,m_nonShiftBilateralFilter.filterTime);
	//	tempString1="m_variationF          ";
	//	tempString2="m_variationG1         ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationF,tempString2,m_variationG1);
	//	tempString1="m_variationG2         ";
	//	tempString2="m_variationH1         ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationG2,tempString2,m_variationH1);
	//	tempString1="m_variationH2         ";
	//	tempString2="m_variationH3         ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationH2,tempString2,m_variationH3);

	//	float tempMSE=0;
	//	//if(m_colors!=NULL)
	//	//	delete[] m_colors;
	//	//m_colors=new float[m_numOfPoints*3];
	//	for(int i=0;i<m_numOfPoints;i++){
	//		float* tempDistance;
	//		std::map<int ,KnearestField>::iterator mapIterator=m_nonShiftBilateralFilter.m_mapKnearest.find(i);
	//		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	//		tempDistance=new float[(*mapIterator).second.m_nearest.size()+1];
	//		tempDistance[0]=((m_originalOriginalPointSets[i*3]-m_pointSets[i*3])*(m_originalOriginalPointSets[i*3]-m_pointSets[i*3])
	//			+(m_originalOriginalPointSets[i*3+1]-m_pointSets[i*3+1])*(m_originalOriginalPointSets[i*3+1]-m_pointSets[i*3+1])
	//			+(m_originalOriginalPointSets[i*3+2]-m_pointSets[i*3+2])*(m_originalOriginalPointSets[i*3+2]-m_pointSets[i*3+2])
	//			);           
	//		int k=1;
	//		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
	//			tempDistance[k]=((m_originalOriginalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])*(m_originalOriginalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])
	//				+(m_originalOriginalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])*(m_originalOriginalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])
	//				+(m_originalOriginalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])*(m_originalOriginalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])
	//				);
	//			/*tempDistance[k]=((m_originalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])*(m_originalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])
	//			+(m_originalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])*(m_originalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])
	//			+(m_originalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])*(m_originalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])
	//			);*/
	//			if(tempDistance[k]<tempDistance[0]){
	//				tempDistance[0]=tempDistance[k];
	//			}
	//			k++;
	//		}

	//		/*		if(tempDistance[0]<=0.05){
	//		m_colors[i*3]=0;
	//		m_colors[i*3+1]=0;
	//		m_colors[i*3+2]=0.5+4*tempDistance[0];
	//		}
	//		if(tempDistance[0]>0.05&tempDistance[0]<=0.15){
	//		m_colors[i*3]=0;
	//		m_colors[i*3+1]=0.5+4*(tempDistance[0]-0.1);
	//		m_colors[i*3+2]=0;
	//		}
	//		if(tempDistance[0]>0.15){
	//		m_colors[i*3]=0.3+0.25*(tempDistance[0]-0.2);
	//		m_colors[i*3+1]=0;
	//		m_colors[i*3+2]=0;
	//		}*/




	//		tempMSE+=tempDistance[0];
	//		/*tempMSE+=((m_originalPointSets[i*3]-m_pointSets[i*3])*(m_originalPointSets[i*3]-m_pointSets[i*3])
	//		+(m_originalPointSets[i*3+1]-m_pointSets[i*3+1])*(m_originalPointSets[i*3+1]-m_pointSets[i*3+1])
	//		+(m_originalPointSets[i*3+2]-m_pointSets[i*3+2])*(m_originalPointSets[i*3+2]-m_pointSets[i*3+2])
	//		);  */    
	//		delete[] tempDistance;

	//	}

	//	tempMSE=sqrt(tempMSE/(m_numOfPoints*m_meanlength*3.14*0.52517*0.52517*m_meanlength));
	//	tempString1="tempMSE               ";
	//	tempString2="m_normalStop          ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,tempMSE,tempString2,m_variationNstop);
	//	//fprintf(fpout,"%d %f %f %f\n",m_K,m_variationF,m_variationG1,tempMSE);
	//	fclose(fpout);
	//}  
	//filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonShiftPoint.txt";
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{

	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f\n",m_pointSets[i*3],m_pointSets[i*3+1],m_pointSets[i*3+2]);		
	//	}

	//	fclose(fpout);
	//}

	m_MLSfilter.DeleteMLSfilter();
	this->UpdateAllViews(NULL);

}

void CMcubeDoc::OnAnisoSurFunFilter()
{
	// TODO: Add your command handler code here
	/*
	 *	同时滤波表面和函数
	 */
	float m_meanlength;	
	m_meanlength=CalculateMeanLength();
	m_dimensionColor=3;
	AnisoRBFfilterSurfaceFunction m_meshLessFilter;
	m_meshLessFilter.GetMeshlessFilter(m_numOfPoints,m_originalPointSets,m_colors,m_K,m_radius*m_meanlength,m_iterativeTimes,m_timeStep,m_loadConstant,m_positionError*m_meanlength,m_stopN,m_maxIter,m_curveThreshold,m_thresholdColor);
	//for(int i=0;i<m_numOfPoints*3;i++){
	//	m_originalPointSets[i]=m_originalPointSets[i]*100000;
	//}
	for(int i=0;i<m_numOfPoints*3;i++){
		m_pointSets[i]=(float)m_meshLessFilter.m_resultPointSet[i];	
		m_colors[i]=m_meshLessFilter.m_resultColor[i];
	}
	//if(m_visualization==1){
	//	if(m_colors!=NULL)
	//		delete m_colors;
	//	m_colors=new float[m_numOfPoints*3];
	//	for(int i=0;i<m_numOfPoints*3;i++){
	//		m_colors[i]=m_nonShiftBilateralFilter.m_resultColor[i];			
	//	}
	//}

	delete[] m_normals;
	int m_nNormals = m_numOfPoints;
	m_normals=new float[m_nNormals*3];

	// Set all normals to 0.
	for (int i = 0; i < m_nNormals*3; i++) {
		m_normals[i]=0;
	}

	// Calculate normals.
	for (int i = 0; i < m_numOfTriangles; i++) {
		VECTOR3D vec1, vec2, normal;
		int id0, id1, id2;
		id0 = m_triangles[i*3];
		id1 = m_triangles[i*3+1];
		id2 = m_triangles[i*3+2];
		vec1[0] = m_pointSets[id1*3]- m_pointSets[id0*3];
		vec1[1] = m_pointSets[id1*3+1] - m_pointSets[id0*3+1];
		vec1[2] = m_pointSets[id1*3+2] - m_pointSets[id0*3+2];
		vec2[0] = m_pointSets[id2*3] - m_pointSets[id0*3];
		vec2[1] = m_pointSets[id2*3+1] - m_pointSets[id0*3+1];
		vec2[2] = m_pointSets[id2*3+2] - m_pointSets[id0*3+2];
		normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
		normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
		normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
		m_normals[id0*3+0] += normal[0];
		m_normals[id0*3+1] += normal[1];
		m_normals[id0*3+2] += normal[2];
		m_normals[id1*3+0] += normal[0];
		m_normals[id1*3+1] += normal[1];
		m_normals[id1*3+2] += normal[2];
		m_normals[id2*3+0] += normal[0];
		m_normals[id2*3+1] += normal[1];
		m_normals[id2*3+2] += normal[2];
	}

	// Normalize normals.
	for (int i = 0; i < m_nNormals; i++) {
		float length = sqrt(m_normals[i*3+0]*m_normals[i*3+0] + m_normals[i*3+1]*m_normals[i*3+1] + m_normals[i*3+2]*m_normals[i*3+2]);
		m_normals[i*3] /= length;
		m_normals[i*3+1] /= length;
		m_normals[i*3+2] /= length;
	}



	CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\time.txt";

	FILE *fpout;
	if((fpout=fopen(filename_pw,"w"))==NULL){
		int dkjkd;
		return;
	}
	else{
		fprintf(fpout,"%f %f %f\n",(float)(m_meshLessFilter.nearestTime/1000/m_iterativeTimes),(float)(m_meshLessFilter.filterTime/1000/m_iterativeTimes),m_meanlength);		
		fclose(fpout);
	}
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	return;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	CString tempString1,tempString2;
	//	tempString1="ifMeanShift           ";
	//	tempString2="whichFunction         ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_meanShift,tempString2,m_indexFunction);
	//	tempString1="m_numOFnearset        ";
	//	tempString2="m_numOFiteract        ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_K,tempString2,m_iterativeTimes);

	//	tempString1="m_numOfPoints         ";
	//	tempString2="m_numOfTriangles      ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_numOfPoints,tempString2,m_numOfTriangles);
	//	tempString1="searchKnearestTime    ";
	//	tempString2="filterTime            ";
	//	fprintf(fpout,"%s %d %s %d\n",tempString1,m_nonShiftBilateralFilter.nearestTime,tempString2,m_nonShiftBilateralFilter.filterTime);
	//	tempString1="m_variationF          ";
	//	tempString2="m_variationG1         ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationF,tempString2,m_variationG1);
	//	tempString1="m_variationG2         ";
	//	tempString2="m_variationH1         ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationG2,tempString2,m_variationH1);
	//	tempString1="m_variationH2         ";
	//	tempString2="m_variationH3         ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,m_variationH2,tempString2,m_variationH3);

	//	float tempMSE=0;
	//	//if(m_colors!=NULL)
	//	//	delete[] m_colors;
	//	//m_colors=new float[m_numOfPoints*3];
	//	for(int i=0;i<m_numOfPoints;i++){
	//		float* tempDistance;
	//		std::map<int ,KnearestField>::iterator mapIterator=m_nonShiftBilateralFilter.m_mapKnearest.find(i);
	//		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	//		tempDistance=new float[(*mapIterator).second.m_nearest.size()+1];
	//		tempDistance[0]=((m_originalOriginalPointSets[i*3]-m_pointSets[i*3])*(m_originalOriginalPointSets[i*3]-m_pointSets[i*3])
	//			+(m_originalOriginalPointSets[i*3+1]-m_pointSets[i*3+1])*(m_originalOriginalPointSets[i*3+1]-m_pointSets[i*3+1])
	//			+(m_originalOriginalPointSets[i*3+2]-m_pointSets[i*3+2])*(m_originalOriginalPointSets[i*3+2]-m_pointSets[i*3+2])
	//			);           
	//		int k=1;
	//		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
	//			tempDistance[k]=((m_originalOriginalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])*(m_originalOriginalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])
	//				+(m_originalOriginalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])*(m_originalOriginalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])
	//				+(m_originalOriginalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])*(m_originalOriginalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])
	//				);
	//			/*tempDistance[k]=((m_originalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])*(m_originalPointSets[i*3]-m_pointSets[(*vectorIterator)*3])
	//			+(m_originalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])*(m_originalPointSets[i*3+1]-m_pointSets[(*vectorIterator)*3+1])
	//			+(m_originalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])*(m_originalPointSets[i*3+2]-m_pointSets[(*vectorIterator)*3+2])
	//			);*/
	//			if(tempDistance[k]<tempDistance[0]){
	//				tempDistance[0]=tempDistance[k];
	//			}
	//			k++;
	//		}

	//		/*		if(tempDistance[0]<=0.05){
	//		m_colors[i*3]=0;
	//		m_colors[i*3+1]=0;
	//		m_colors[i*3+2]=0.5+4*tempDistance[0];
	//		}
	//		if(tempDistance[0]>0.05&tempDistance[0]<=0.15){
	//		m_colors[i*3]=0;
	//		m_colors[i*3+1]=0.5+4*(tempDistance[0]-0.1);
	//		m_colors[i*3+2]=0;
	//		}
	//		if(tempDistance[0]>0.15){
	//		m_colors[i*3]=0.3+0.25*(tempDistance[0]-0.2);
	//		m_colors[i*3+1]=0;
	//		m_colors[i*3+2]=0;
	//		}*/




	//		tempMSE+=tempDistance[0];
	//		/*tempMSE+=((m_originalPointSets[i*3]-m_pointSets[i*3])*(m_originalPointSets[i*3]-m_pointSets[i*3])
	//		+(m_originalPointSets[i*3+1]-m_pointSets[i*3+1])*(m_originalPointSets[i*3+1]-m_pointSets[i*3+1])
	//		+(m_originalPointSets[i*3+2]-m_pointSets[i*3+2])*(m_originalPointSets[i*3+2]-m_pointSets[i*3+2])
	//		);  */    
	//		delete[] tempDistance;

	//	}

	//	tempMSE=sqrt(tempMSE/(m_numOfPoints*m_meanlength*3.14*0.52517*0.52517*m_meanlength));
	//	tempString1="tempMSE               ";
	//	tempString2="m_normalStop          ";
	//	fprintf(fpout,"%s %f %s %f\n",tempString1,tempMSE,tempString2,m_variationNstop);
	//	//fprintf(fpout,"%d %f %f %f\n",m_K,m_variationF,m_variationG1,tempMSE);
	//	fclose(fpout);
	//}  
	filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\resultPoint.txt";
	if((fpout = fopen(filename_pw, "w")) == NULL)
	{
		int dkjkd;
		//MessageBox("can't open the file!");
	}
	else
	{

		for(int i=0;i<m_numOfPoints;i++){
			fprintf(fpout,"%f %f %f\n",m_pointSets[i*3],m_pointSets[i*3+1],m_pointSets[i*3+2]);		
		}

		fclose(fpout);
	}

	m_meshLessFilter.DeleteMeshlessFilter();
	this->UpdateAllViews(NULL);
}


void CMcubeDoc::ComputeError()
{
    m_maxError=0;
	m_meanError=0;

	for(int i=0;i<m_numOfPoints;i++){
		float error=0;
		for(int j=0;j<3;j++){
			error+=(m_pointSets[i*3+j]-m_originalOriginalPointSets[i*3+j])*m_normals[i*3+j];
		}
		error=abs(error);	
		if(error>m_maxError)
			m_maxError=error;
		m_meanError+=error;
	}
	m_meanError/=m_numOfPoints;
	return;
}