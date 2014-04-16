; CLW file contains information for the MFC ClassWizard

[General Info]
Version=1
LastClass=CMcubeDoc
LastTemplate=CDialog
NewFileInclude1=#include "stdafx.h"
NewFileInclude2=#include "Mcube.h"
LastPage=0

ClassCount=7
Class1=CMcubeApp
Class2=CMcubeDoc
Class3=CMcubeView
Class4=CMainFrame

ResourceCount=3
Resource1=IDD_ABOUTBOX
Resource2=IDR_MAINFRAME
Class5=CAboutDlg
Class6=CFormCommandView
Class7=CDialogIntepolate
Resource3=IDD_FORMVIEW (English (U.S.))

[CLS:CMcubeApp]
Type=0
HeaderFile=Mcube.h
ImplementationFile=Mcube.cpp
Filter=N

[CLS:CMcubeDoc]
Type=0
HeaderFile=McubeDoc.h
ImplementationFile=McubeDoc.cpp
Filter=N
BaseClass=CDocument
VirtualFilter=DC
LastObject=CMcubeDoc

[CLS:CMcubeView]
Type=0
HeaderFile=McubeView.h
ImplementationFile=McubeView.cpp
Filter=C
BaseClass=CView
VirtualFilter=VWC


[CLS:CMainFrame]
Type=0
HeaderFile=MainFrm.h
ImplementationFile=MainFrm.cpp
Filter=T
BaseClass=CFrameWnd
VirtualFilter=fWC
LastObject=CMainFrame




[CLS:CAboutDlg]
Type=0
HeaderFile=Mcube.cpp
ImplementationFile=Mcube.cpp
Filter=D
LastObject=CAboutDlg

[DLG:IDD_ABOUTBOX]
Type=1
Class=CAboutDlg
ControlCount=4
Control1=IDC_STATIC,static,1342177283
Control2=IDC_STATIC,static,1342308480
Control3=IDC_STATIC,static,1342308352
Control4=IDOK,button,1342373889

[MNU:IDR_MAINFRAME]
Type=1
Class=CMainFrame
Command1=ID_FILE_NEW
Command2=ID_FILE_OPEN
Command3=ID_FILE_SAVE
Command4=ID_FILE_SAVE_AS
Command5=ID_FILE_PRINT
Command6=ID_FILE_PRINT_PREVIEW
Command7=ID_FILE_PRINT_SETUP
Command8=ID_FILE_MRU_FILE1
Command9=ID_APP_EXIT
Command10=ID_EDIT_UNDO
Command11=ID_EDIT_CUT
Command12=ID_EDIT_COPY
Command13=ID_EDIT_PASTE
Command14=ID_VIEW_TOOLBAR
Command15=ID_VIEW_STATUS_BAR
Command16=ID_MCUBE_CALCULATER
Command17=ID_VECTOR_INTEPOLATE
Command18=ID_MCUBE_NODEMOVE
Command19=ID_HELP_FINDER
Command20=ID_APP_ABOUT
CommandCount=20

[ACL:IDR_MAINFRAME]
Type=1
Class=CMainFrame
Command1=ID_FILE_NEW
Command2=ID_FILE_OPEN
Command3=ID_FILE_SAVE
Command4=ID_FILE_PRINT
Command5=ID_EDIT_UNDO
Command6=ID_EDIT_CUT
Command7=ID_EDIT_COPY
Command8=ID_EDIT_PASTE
Command9=ID_EDIT_UNDO
Command10=ID_EDIT_CUT
Command11=ID_EDIT_COPY
Command12=ID_EDIT_PASTE
Command13=ID_NEXT_PANE
Command14=ID_PREV_PANE
Command15=ID_CONTEXT_HELP
Command16=ID_HELP
CommandCount=16

[TB:IDR_MAINFRAME]
Type=1
Class=?
Command1=ID_FILE_NEW
Command2=ID_FILE_OPEN
Command3=ID_FILE_SAVE
Command4=ID_EDIT_CUT
Command5=ID_EDIT_COPY
Command6=ID_EDIT_PASTE
Command7=ID_FILE_PRINT
Command8=ID_APP_ABOUT
Command9=ID_CONTEXT_HELP
CommandCount=9

[DLG:IDD_FORMVIEW (English (U.S.))]
Type=1
Class=CFormCommandView
ControlCount=27
Control1=IDC_STATIC,button,1342177287
Control2=IDC_STATIC,static,1342308352
Control3=IDC_STATIC,static,1342308352
Control4=IDC_STATIC_BACKGROUND,static,1342177287
Control5=ID_STATIC_SPECULAR,static,1342177287
Control6=IDC_STATIC_DIFFUSE,static,1342177287
Control7=IDC_STATIC_AMBIENT,static,1342177287
Control8=IDC_STATIC,button,1342177287
Control9=IDC_CHECK1,button,1342242819
Control10=IDC_CHECK2,button,1342242819
Control11=IDC_CHECK3,button,1342242819
Control12=IDC_RADIO1,button,1342308361
Control13=IDC_RADIO2,button,1342177289
Control14=IDC_RADIO3,button,1342177289
Control15=IDC_STATIC,button,1342177287
Control16=IDC_CHECK4,button,1342242819
Control17=IDC_STATIC,button,1342177287
Control18=IDC_CHECK5,button,1342242819
Control19=IDC_SLIDER1,msctls_trackbar32,1342242840
Control20=IDC_STATIC,static,1342308352
Control21=IDC_STATIC,static,1342308352
Control22=IDC_SLIDER2,msctls_trackbar32,1342242840
Control23=IDC_STATIC,static,1342308352
Control24=IDC_SLIDER3,msctls_trackbar32,1342242840
Control25=IDC_STATIC,static,1342308352
Control26=IDC_STATIC,static,1342308352
Control27=IDC_STATIC,static,1342308352

[CLS:CFormCommandView]
Type=0
HeaderFile=FormCommandView.h
ImplementationFile=FormCommandView.cpp
BaseClass=CFormView
Filter=D
LastObject=IDC_SLIDER3
VirtualFilter=VWC

[CLS:CDialogIntepolate]
Type=0
HeaderFile=DialogIntepolate.h
ImplementationFile=DialogIntepolate.cpp
BaseClass=CDialog
Filter=D
VirtualFilter=dWC
LastObject=IDC_CHECK3

