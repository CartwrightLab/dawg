// snout.cpp : Defines the class behaviors for the application.
//

#include "stdafx.h"
#include "snout.h"
#include "snoutDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CSnoutApp

BEGIN_MESSAGE_MAP(CSnoutApp, CWinApp)
	ON_COMMAND(ID_HELP, CWinApp::OnHelp)
END_MESSAGE_MAP()


// CSnoutApp construction

CSnoutApp::CSnoutApp()
{
	// TODO: add construction code here,
	// Place all significant initialization in InitInstance
}


// The one and only CSnoutApp object

CSnoutApp theApp;


// CSnoutApp initialization

BOOL CSnoutApp::InitInstance()
{
	// InitCommonControls() is required on Windows XP if an application
	// manifest specifies use of ComCtl32.dll version 6 or later to enable
	// visual styles.  Otherwise, any window creation will fail.
	InitCommonControls();

	CWinApp::InitInstance();

	AfxEnableControlContainer();


	CSnoutDlg dlg;
	m_pMainWnd = &dlg;
	dlg.DoModal();
	dlg.m_dawgout.Close();
	return FALSE;
}
