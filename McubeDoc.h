// McubeDoc.h : interface of the CMcubeDoc class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_MCUBEDOC_H__62FE69EC_BC70_4E05_8D54_D407A8AACF6F__INCLUDED_)
#define AFX_MCUBEDOC_H__62FE69EC_BC70_4E05_8D54_D407A8AACF6F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "stdafx.h"
#include "CIsoSurface.h"
//#include "xVolume.h"
#include "ElementBuild.h"




class CMcubeDoc : public CDocument
{
protected: // create from serialization only
	CMcubeDoc();
	DECLARE_DYNCREATE(CMcubeDoc)

// Attributes
public:
  // CIsoSurface<uint8> m_CIsoSurface;
	CIsoSurface<float> m_CIsoSurface;
	//xVolume *GetVolume(){return &m_volume;};
	//xVolume m_volume;
   float m_fVolLengthX;
   float m_fVolLengthY;
   float m_fVolLengthZ;
   int m_Isolevel;
   //三角片信息
   float* m_pointSets;
   float* m_originalPointSets;
   float* m_originalOriginalPointSets;
   float* m_originalColors;
   float* m_colors;
   int* m_triangles;
   float* temp_normals; //加噪声和计算误差时候的法向
   float* m_normals;
   int m_numOfPoints;
   int m_numOfTriangles;
   

   float x_max,x_min,y_max,y_min,z_max,z_min;
// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CMcubeDoc)
	public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
	virtual BOOL OnOpenDocument(LPCTSTR lpszPathName);
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CMcubeDoc();

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

    float* m_pVol;
    //uint8 * m_pVol;
	int m_nSizeZ;
	int m_nSizeY;
	int m_nSizeX;	

	
	

	//MeshTriangle* m_nMeshTriangleOne;
	//MeshTriangle* m_nMeshTriangleTwo;
	
// Generated message map functions
public:
	//////////////////////////////////////////////////////////////////////////
	// 滤波参数
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
	float m_noiseCoffecient;
	float m_functionVariation;
	int m_visualization;
	float m_noiseFuncCoff;
	float m_thresholdColor;
	float m_functionGradientWide;
	float m_functionGradientThreshold;
	int m_dimensionColor;
	float m_thresholdDistanceNormal;
	float m_thresholdDistanceGradientFunction;
	float m_timeStep;
	float m_loadConstant;
	float m_positionError;
	float m_radius;
    int m_stopN;
	int m_maxIter;
	double m_curveThreshold;

	int m_tangentOrManifold;
	int m_ifNormalWeight;
	int m_ifAreaWeight;
	int m_ifVolumePreserve;
	int m_ifVariationNormal;


	float m_maxError; //错误定义为与原始点的距离的法向投影的大小
	float m_meanError;
	//////////////////////////////////////////////////////////////////////////
	
protected:
	void Filter();
	float m_fCellLengthZ;
	float m_fCellLengthY;
	float m_fCellLengthX;
	//{{AFX_MSG(CMcubeDoc)
	afx_msg void OnMcubeCalculater();
	afx_msg void OnVectorIntepolate();
	afx_msg void OnMcubeNodemove();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
protected:
	void generateCube();
	void AddPositionGaussianNoise();
	void AddColorGaussianNorse();
	void GenerateTriangleFromImage();
	float CalculateMeanLength();
	void ComputeError();

public:
	afx_msg void OnEnChangeEditIsovalue();
	afx_msg void OnBnClickedButton1();
	afx_msg void OnMeshSkin();
	afx_msg void OnMeshBone();
	afx_msg void OnFilterBilateral();
	afx_msg void OnFilterAddnoiseGeometry();
	afx_msg void OnFilterAddnoiseFunction();
	afx_msg void OnFilterNonShift();
	afx_msg void OnMeshlessFilter();
	afx_msg void OnFilterAnisoMeshless();
	afx_msg void OnFilterMls();
	afx_msg void OnAnisoSurFunFilter();
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_MCUBEDOC_H__62FE69EC_BC70_4E05_8D54_D407A8AACF6F__INCLUDED_)
