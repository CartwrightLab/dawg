// snoutDlg.h : header file
//
#include "afxwin.h"
#include "fileeditctrl.h"
#include "redir.h"

class CDawgOutput : public CRedirector
{
public:
	CDawgOutput() : m_pEdit(NULL) { }
	virtual ~CDawgOutput() { }
	
	CEdit* m_pEdit;

protected:
	// overrides:
	virtual void WriteStdOut(LPCSTR pszOutput);
	virtual void WriteStdError(LPCSTR pszError);
};

// CSnoutDlg dialog
class CSnoutDlg : public CDialog
{
// Construction
public:
	CSnoutDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	enum { IDD = IDD_SNOUT_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support

// Implementation
public:
	CFileEditCtrl m_ecFilename;
	CEdit m_ecOutput;
	CDawgOutput m_dawgout;
	CFont m_font;

protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	afx_msg void OnBnClickedSave();
	afx_msg void OnBnClickedEdit();
	afx_msg void OnStnDblclickAbout();
	afx_msg void OnBnClickedRun();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	DECLARE_MESSAGE_MAP()
};
