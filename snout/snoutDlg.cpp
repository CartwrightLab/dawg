// snoutDlg.cpp : implementation file
//

#include "stdafx.h"
#include "snout.h"
#include "snoutDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

void CDawgOutput::WriteStdOut(LPCSTR pszOutput)
{
	if (m_pEdit != NULL)
	{
		int nSize = m_pEdit->GetWindowTextLength();
		m_pEdit->SetSel(nSize, nSize);
		m_pEdit->ReplaceSel(pszOutput);		// add the message to the end of Edit control
	}
}

void CDawgOutput::WriteStdError(LPCSTR pszError)
{
	WriteStdOut(pszError);
}

// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()


// CSnoutDlg dialog
CSnoutDlg::CSnoutDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CSnoutDlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CSnoutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_FileEditCtrl(pDX, IDC_FILENAME, m_ecFilename, FEC_FILEOPEN);
	DDX_Control(pDX, IDC_OUTPUT, m_ecOutput);
}

BEGIN_MESSAGE_MAP(CSnoutDlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
	ON_BN_CLICKED(IDC_RUN, OnBnClickedRun)
	ON_BN_CLICKED(IDC_SAVE, OnBnClickedSave)
	ON_BN_CLICKED(IDC_EDIT, OnBnClickedEdit)
	ON_STN_DBLCLK(IDC_ABOUT, OnStnDblclickAbout)
END_MESSAGE_MAP()

// CSnoutDlg message handlers

BOOL CSnoutDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	HANDLE handle = ::LoadImage(AfxGetInstanceHandle(),MAKEINTRESOURCE(IDI_FOLDER),
                                IMAGE_ICON,0,0,LR_DEFAULTCOLOR);	
	m_ecFilename.SetButtonImage(handle, PJAI_ICON|PJAI_AUTODELETE);
	m_ecFilename.SetButtonWidth(21);
	m_ecFilename.GetOpenFileName()->lpstrFilter = _T("Fud Files (*.fud)\0*.fud\0All Files (*.*)\0*.*\0");

	m_dawgout.m_pEdit = &m_ecOutput;
	m_font.CreatePointFont(100, "Courier New");
	m_ecOutput.SetFont(&m_font);
	
	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CSnoutDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

void CSnoutDlg::OnPaint() 
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

HCURSOR CSnoutDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

void CSnoutDlg::OnBnClickedRun()
{
	//Retreive filename
	UpdateData();
	TCHAR tcFileName[MAX_PATH];
	m_ecFilename.GetWindowText(tcFileName, MAX_PATH);
	CString csCmd = _T("dawg.exe \"");
	csCmd.Append(tcFileName);
	csCmd.Append(_T("\""));
	m_ecOutput.SetWindowText("");
	m_dawgout.Close();
	m_dawgout.Open(csCmd);
}

void CSnoutDlg::OnBnClickedSave()
{
	CFileDialog dlg(false, _T("poo"), NULL, OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
		_T("Poo Files (*.poo)|*.poo|All Files (*.*)|*.*"), this);

	if(dlg.DoModal() == IDOK)
	{
		CString csOut;
		CStdioFile fout;
		CFileException exc;
		if(!fout.Open(dlg.GetFileName(), CFile::modeCreate|CFile::modeWrite|CFile::typeText, &exc))
		{
			exc.ReportError();
			return;
		}
		m_ecOutput.GetWindowText(csOut);
		fout.WriteString(csOut);
	}
}

void CSnoutDlg::OnBnClickedEdit()
{
	CString csOut, csFile,csTemp;
	CStdioFile fin;
	m_ecFilename.GetWindowText(csFile);
	//ShellExecute(m_hWnd, _T("Edit"), csFile, NULL, NULL, SW_SHOWNORMAL);
	ShellExecute(m_hWnd, _T("Open"), _T("notepad.exe"), csFile, NULL, SW_SHOWNORMAL);
}

void CSnoutDlg::OnStnDblclickAbout()
{
	CAboutDlg dlg;
	dlg.DoModal();	
}