// PJA_Icon.cpp : implementation file
//
// written by PJ Arends
// pja@telus.net
//
// some of this code is 'borrowed' from Chris Maunder's CHyperLink class.
// (can be found at http://www.codeproject.com/miscctrl/hyperlink.asp)


#include "stdafx.h"
#include "PJA_Icon.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/*///////////////////////////////////////////////////////////////////////////

                   How I got the data to make the icon
                          (or .ico file format)

The first step was to create the icon. I used VS resource editor to create a
custom 24 pixel by 19 pixel by 2 colour icon which I saved to a .ico file. I
then opened the file in binary mode, which showed me the bytes in hexadecimal
format.

The first WORD in the file will always be zero, and the second WORD contains
the image type, it is one for icons, and two for cursors. The third WORD
contains the number of icon images found in the file (in my case it was one).
Following these three WORDs there is an array of sixteen byte icon headers,
one header for each icon image in the file. The header format is :

        Width    (BYTE)  - The width of the icon in pixels
        Height   (BYTE)  - The height of the icon in pixels
        Colours  (BYTE)  - The number of colours in the icon
        Reserved (BYTE)  - Alignment byte
        Planes   (WORD)  - The number of planes
        BitCount (WORD)  - The number of bits per pixel
        Size     (DWORD) - The number of data bytes in the icon image
        Offset   (DWORD) - The offset (in bytes) from the start of the .ico file
                           to the start of the data bytes for this icon image

So my .ico file in hexadecimal looked like this :

  00 00|01 00|01 00|18|13|02|00|01 00|01 00|C8 00 00 00|16 00 00 00|28 00 00 00 18...
       |     |     |  |  |  |  |     |     |           |           |
  zero  icon    #  |W  H  C  R  Plane  BitC    Size        Offset  |  Data Bytes...
             images|                                               |
                   |<---------------- Icon Header ---------------->|

The information I needed from the header was the width, height, and size of the icon.
And then, starting 'Offset' bytes into the file, I copied 'Size' bytes into the data buffer.

*////////////////////////////////////////////////////////////////////////////

#define LINK_TEXT   _T("http://www3.telus.net/pja/cppcode.htm")
#define ICON_WIDTH  0x18
#define ICON_HEIGHT 0x13
#define ICON_SIZE   0xC8

/////////////////////////////////////////////////////////////////////////////
// CPJA_Icon

CPJA_Icon::CPJA_Icon()
{
    BYTE data[ICON_SIZE] = {0x28, 0x00, 0x00, 0x00, 0x18, 0x00, 0x00, 0x00, 0x26, 0x00,
                            0x00, 0x00, 0x01, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00,
                            0x4C, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                            0x00, 0x00, 0x00, 0x00, 0xFF, 0xFF, 0xFF, 0x00, 0x00, 0x00,

                            0x00, 0x81, 0x00, 0x00, 0x00, 0x81, 0x00, 0x00, 0x00, 0x4E,
                            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                            0x00, 0x03, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00,
                            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00,
                            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,

                            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
                            0x00, 0x03, 0x00, 0x00, 0x00, 0x04, 0x00, 0x00, 0x00, 0x81,
                            0x00, 0x00, 0x00, 0x00, 0xF8, 0x00, 0x1F, 0x00, 0xE7, 0xFF,
                            0xE7, 0x00, 0xDF, 0xFF, 0xFB, 0x00, 0xB7, 0xEF, 0xFD, 0x00,
                            0xB7, 0xF7, 0xFD, 0x00, 0x77, 0xF7, 0xFE, 0x00, 0x74, 0x77,

                            0x0A, 0x00, 0x73, 0xB6, 0xF6, 0x00, 0x77, 0xB6, 0xF6, 0x00,
                            0x77, 0xB6, 0xF6, 0x00, 0x77, 0xB7, 0x06, 0x00, 0x73, 0xB6,
                            0xF6, 0x00, 0x74, 0x77, 0x0E, 0x00, 0x7F, 0xFF, 0xFE, 0x00,
                            0xBF, 0xFF, 0xFD, 0x00, 0xBF, 0xF7, 0xFD, 0x00, 0xDF, 0xFF,
                            0xFB, 0x00, 0xE7, 0xFF, 0xE7, 0x00, 0xF8, 0x00, 0x1F, 0x00};

    hIcon = ::CreateIconFromResourceEx(data,
                                       ICON_SIZE,
                                       TRUE,
                                       0x00030000,
                                       ICON_WIDTH,
                                       ICON_HEIGHT,
                                       LR_DEFAULTCOLOR);
    hCursor = NULL;
    Link = LINK_TEXT;
}

BEGIN_MESSAGE_MAP(CPJA_Icon, CStatic)
    //{{AFX_MSG_MAP(CPJA_Icon)
    ON_CONTROL_REFLECT(STN_CLICKED, OnClicked)
    ON_WM_PAINT()
    ON_WM_SETCURSOR()
    //}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CPJA_Icon message handlers

void CPJA_Icon::Draw(CDC *pDC, int x /* = 0 */, int y /* = 0 */)
{
    ::DrawIconEx(pDC->m_hDC,
                 x,
                 y,
                 hIcon,
                 0,
                 0,
                 0,
                 NULL,
                 DI_NORMAL);
}

HICON CPJA_Icon::GetIconHandle()
{
    return hIcon;
}

void CPJA_Icon::OnClicked() 
{
    ShellExecute(NULL, _T("open"), Link, NULL, NULL, SW_SHOWNORMAL);
}

void CPJA_Icon::OnPaint() 
{
    CPaintDC dc(this); // device context for painting
    Draw(&dc);
}

BOOL CPJA_Icon::OnSetCursor(CWnd*, UINT, UINT) 
{
    if (hCursor)
    {
        ::SetCursor(hCursor);
        return TRUE;
    }
    return FALSE;
}

void CPJA_Icon::PreSubclassWindow() 
{   // Set SS_NOTIFY style so we get mouse messages
    ModifyStyle (0, SS_NOTIFY);
    SetDefaultCursor();

    // size the window so the icon fills the client area, without
    // moving the bottom right corner of the window
    CRect rc;
    GetClientRect(&rc);
    int width = rc.Width();
    GetWindowRect(&rc);
    GetParent()->ScreenToClient(&rc);
    int border = (rc.Width() - width);
    SetWindowPos(NULL,
                 rc.right - ICON_WIDTH - border,
                 rc.bottom - ICON_HEIGHT - border,
                 ICON_WIDTH + border,
                 ICON_HEIGHT + border,
                 SWP_NOZORDER);

    // create and setup the tooltip
    GetClientRect(&rc);
    ToolTip.Create(this);
    ToolTip.AddTool(this, Link, rc, 1);
    ToolTip.SetFont(GetParent()->GetFont());

    CStatic::PreSubclassWindow();
}

BOOL CPJA_Icon::PreTranslateMessage(MSG* pMsg) 
{   // send messages to the tooltip
    ToolTip.RelayEvent(pMsg);
    return CStatic::PreTranslateMessage(pMsg);
}

void CPJA_Icon::SetDefaultCursor()
{
// The following appeared in Paul DiLascia's Jan 1998 MSJ articles,
// It loads a "hand" cursor from the winhlp32.exe module
    if (hCursor == NULL)                // No cursor handle - load our own
    {
        // Get the windows directory
        CString strWndDir;
        GetWindowsDirectory(strWndDir.GetBuffer(MAX_PATH), MAX_PATH);
        strWndDir.ReleaseBuffer();

        strWndDir += _T("\\winhlp32.exe");
        // This retrieves cursor #106 from winhlp32.exe, which is a hand pointer
        HMODULE hModule = LoadLibrary(strWndDir);
        if (hModule) {
            HCURSOR hHandCursor = ::LoadCursor(hModule, MAKEINTRESOURCE(106));
            if (hHandCursor)
                hCursor = CopyCursor(hHandCursor);
        }
        FreeLibrary(hModule);
    }

}

void CPJA_Icon::SetLink(LPCTSTR theLink)
{
    ASSERT(theLink);
    Link = theLink;
    if (IsWindow(m_hWnd))
        ToolTip.UpdateTipText(Link, this, 1);
}

void CPJA_Icon::SetLinkCursor(HCURSOR Cursor)
{
    hCursor = Cursor;
    if (hCursor == NULL)
        SetDefaultCursor();
}
