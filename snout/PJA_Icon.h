#if !defined(AFX_PJA_ICON_H__984F0FD2_9725_11D5_B625_E51433507822__INCLUDED_)
#define AFX_PJA_ICON_H__984F0FD2_9725_11D5_B625_E51433507822__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// PJA_Icon.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CPJA_Icon window

class CPJA_Icon : public CStatic
{
// Construction
public:
    CPJA_Icon();

// Attributes
public:

// Operations
public:

// Overrides
    // ClassWizard generated virtual function overrides
    //{{AFX_VIRTUAL(CPJA_Icon)
    public:
    virtual BOOL PreTranslateMessage(MSG* pMsg);
    protected:
    virtual void PreSubclassWindow();
    //}}AFX_VIRTUAL

// Implementation
public:
    void SetLinkCursor(HCURSOR Cursor);
    void SetLink (LPCTSTR theLink);
    void Draw(CDC *pDC, int x = 0, int y = 0);
    HICON GetIconHandle();

    // Generated message map functions
protected:
    CString Link;
    CToolTipCtrl ToolTip;
    HCURSOR hCursor;
    HICON hIcon;
    void SetDefaultCursor();
    //{{AFX_MSG(CPJA_Icon)
    afx_msg void OnPaint();
    afx_msg void OnClicked();
    afx_msg BOOL OnSetCursor(CWnd* pWnd, UINT nHitTest, UINT message);
    //}}AFX_MSG

    DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_PJA_ICON_H__984F0FD2_9725_11D5_B625_E51433507822__INCLUDED_)
