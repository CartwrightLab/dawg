/****************************************************************************
PJAImage.cpp : implementation file for the PJAImage class
written by PJ Arends
pja@telus.net

For updates check http://www3.telus.net/pja/PJAImage.htm

-----------------------------------------------------------------------------
This code is provided as is, with no warranty as to it's suitability or usefulness
in any application in which it may be used. This code has not been tested for
UNICODE builds, nor has it been tested on a network ( with UNC paths ).

This code may be used in any way you desire. This file may be redistributed by any
means as long as it is not sold for profit, and providing that this notice and the
authors name are included.

If any bugs are found and fixed, a note to the author explaining the problem and
fix would be nice.
-----------------------------------------------------------------------------
****************************************************************************/

/////////////////////////////////////////////////////////////////////////////
// CPJAImage class

#include "stdafx.h"
#include "pjaimage.h"

#define BRUSHWIDTH  8
#define BRUSHHEIGHT 8

/////////////////////////////////////////////////////////////////////////////
//
//  CPJAImage constructor  (public member function)
//    Initializes member variables
//
//  Parameters :
//    None
//
//  Returns :
//    Nothing
//
/////////////////////////////////////////////////////////////////////////////

CPJAImage::CPJAImage()
{
    m_ImageFlags = 0;
    m_DrawFlags = 0;
    m_hImage = NULL;
    m_size.cx = 0;
    m_size.cy = 0;
    m_TransparentColour = CLR_DEFAULT;
}

/////////////////////////////////////////////////////////////////////////////
//
//  CPJAImage destructor  (public member function)
//    frees the memory held by the image handle
//
//  Parameters :
//    None
//
//  Returns :
//    Nothing
//
/////////////////////////////////////////////////////////////////////////////

CPJAImage::~CPJAImage()
{
    if (m_ImageFlags & PJAI_AUTODELETE)
    {
        if (m_ImageFlags & PJAI_ICON)
            DestroyIcon((HICON)m_hImage);
        else
            DeleteObject((HGDIOBJ)m_hImage);
    }
}

/////////////////////////////////////////////////////////////////////////////
//
//  CPJAImage::DitherBlt  (protected member function)
//    Draws the image on the FromDC as a disabled (grayed) image onto the pToDC
//
//  Parameters :
//    pToDC  [in] - pointer to the DC to draw the bitmap onto
//    x      [in] - the left side of the image on the destination DC
//    y      [in] - the top  of the image on the destination DC
//    w      [in] - the width of the image on the destination DC
//    h      [in] - the height of the image on the destination DC
//    FromDC [in] - The DC containing the bitmap to be grayed
//
//  Returns :
//    Nothing
//
//  Note : modified from code found at http://www.codeguru.com/bitmap/dither_blt.shtml
//         original author Jean-Edouard Lachand-Robert (iamwired@geocities.com)
//
/////////////////////////////////////////////////////////////////////////////

void CPJAImage::DitherBlt(CDC *pToDC, int x, int y, int w, int h, CDC *pFromDC)
{
    CDC MonoDC;
    if (MonoDC.CreateCompatibleDC(pToDC))
    {
        struct {
            BITMAPINFOHEADER bmiHeader; 
            RGBQUAD          bmiColors[2]; 
        } RGBBWBITMAPINFO = { { sizeof(BITMAPINFOHEADER),
            w,
            h,
            1,
            1,
            BI_RGB,
            0,
            0,
            0,
            0,
            0},
        { { 0x00, 0x00, 0x00, 0x00 },
        { 0xFF, 0xFF, 0xFF, 0x00 }}};
        VOID *pbitsBW;
        HBITMAP hbmBW = CreateDIBSection(MonoDC.m_hDC,
            (LPBITMAPINFO)&RGBBWBITMAPINFO,
            DIB_RGB_COLORS,
            &pbitsBW,
            NULL,
            0);
        ASSERT(hbmBW);

        if (hbmBW)
        {
            int SavedMonoDC = MonoDC.SaveDC();
            int SavedpToDC = pToDC->SaveDC();

            // Attach the monochrome DIB section to the MonoDC
            MonoDC.SelectObject(hbmBW);
            
            // BitBlt the bitmap into the monochrome DIB section
            MonoDC.BitBlt(0, 0, w, h, pFromDC, 0, 0, SRCCOPY);
            
            // BitBlt the black bits in the monochrome bitmap into COLOR_3DHILIGHT bits in the destination DC
            // The magic ROP comes from the Charles Petzold's book
            HBRUSH hb = CreateSolidBrush(GetSysColor(COLOR_3DHILIGHT));
            pToDC->SelectObject(hb);
            pToDC->BitBlt(x + 1, y + 1, w, h, &MonoDC, 0, 0, 0x00B8074A);
            
            // BitBlt the black bits in the monochrome bitmap into COLOR_3DSHADOW bits in the destination DC
            hb = CreateSolidBrush(GetSysColor(COLOR_3DSHADOW));
            DeleteObject(pToDC->SelectObject(hb));
            pToDC->BitBlt(x, y, w, h, &MonoDC, 0, 0, 0x00B8074A);
            
            pToDC->RestoreDC(SavedpToDC);
            MonoDC.RestoreDC(SavedMonoDC);
            DeleteObject(hb);
        }
        DeleteObject(hbmBW);
        MonoDC.DeleteDC();
    }
}

/////////////////////////////////////////////////////////////////////////////
//
//  CPJAImage::DrawImage  (public member function)
//    Draws the image (set with the SetImage() function) on the given device 
//    context
//
//  Parameters :
//    pDC       [in] - a pointer to the the device context to draw the image on
//    x         [in] - the left side of the image on the destination DC
//    y         [in] - the top  of the image on the destination DC
//    w         [in] - the width of the image on the destination DC
//    h         [in] - the height of the image on the destination DC
//    DrawFlags [in] - How to draw the image
//
//  Returns :
//    Nothing
//
//  Note :
//    See the PJAImage.h file for a description of the flags used
//
//    If the image is an icon or a transparent bitmap, the brush that is selected
//    into the pDC will be used as the background. The brush can be a pattern
//    or bitmap brush. This code assumes the brush is an 8x8 pixel brush
//    (Set with the BUTTONWIDTH an BUTTONHEIGHT macros).
//
//    The image will be drawn entirely within the rectangle specified by the
//    x, y, w, and h parameters.
//
/////////////////////////////////////////////////////////////////////////////

void CPJAImage::DrawImage(CDC *pDC, int x, int y, int w, int h, DWORD DrawFlags /* = 0 */)
{   // sanity check
    if (!m_hImage)
        return;

    // verify flags
#ifdef _DEBUG
    if (DrawFlags & PJAI_DISABLED)
        ASSERT (!(DrawFlags & PJAI_GRAYSCALE));
    if (DrawFlags & PJAI_CENTERED)
        ASSERT (!(DrawFlags & PJAI_STRETCHED));
#endif

    m_DrawFlags = DrawFlags;

    // Get the background brush for transparent images
    CBrush *BackGround = pDC->GetCurrentBrush();
    BackGround->UnrealizeObject();

    // handle for the grayscale bitmap
    HBITMAP GrayBmp = NULL;

    // set the clip region to the specified rectangle
    CRgn ClipRgn;
    ClipRgn.CreateRectRgn(x, y, x + w, y + h);
    pDC->SelectClipRgn(&ClipRgn);
    ClipRgn.DeleteObject();

    // create memory DC
    CDC memDC;
    memDC.CreateCompatibleDC(pDC);
    int savedmemDC = memDC.SaveDC();
    CBitmap memDCBmp;

    CDC* pOutputDC = &memDC;

    int left = x;
    int top = y;
    int width = m_size.cx;
    int height = m_size.cy;

    if (m_DrawFlags & PJAI_CENTERED)
    {   // center the image on the output rectangle
        left = x + (w / 2) - (m_size.cx / 2);
        top = y + (h / 2) - (m_size.cy / 2);
    }

    // Create a DC and bitmap for the stretched image
    CDC StretchDC;
    int savedStretchDC = 0;
    CBitmap stretchbmp;
    if (m_DrawFlags & PJAI_STRETCHED)
    {   // stretch image to fit output rectangle
        width = w;
        height = h;
    }

    // get the brush origins in case we are using a bitmap or pattern brush
    CPoint BrushOrg;
    CPoint Origin = pDC->GetBrushOrg();
    BrushOrg.x = (BRUSHWIDTH - (left - Origin.x) % BRUSHWIDTH);
    BrushOrg.y = (BRUSHHEIGHT - (top - Origin.y) % BRUSHHEIGHT);

    // Create a DC and bitmap for the transparent image
    CDC TransparentDC;
    int savedTransparentDC = 0;
    CBitmap Transparentbmp;

    if (m_ImageFlags & PJAI_ICON)
    {   // draw the icon onto the memory DC
        HICON TheIcon = (HICON)m_hImage;
        if (m_DrawFlags & PJAI_GRAYSCALE)
        {   // convert the colour icon to grayscale
            ICONINFO iconinfo;
            GetIconInfo(TheIcon, &iconinfo);
            if (iconinfo.hbmColor)
            {
                HBITMAP grayscale = GrayScale(pDC, iconinfo.hbmColor);

                ::DeleteObject(iconinfo.hbmColor);
                iconinfo.hbmColor = grayscale;
                TheIcon = ::CreateIconIndirect(&iconinfo);
                ::DeleteObject(iconinfo.hbmColor);
                ::DeleteObject(iconinfo.hbmMask);
            }
        }

        memDCBmp.CreateCompatibleBitmap(pDC, width, height);
        memDC.SelectObject(&memDCBmp);
        memDC.SetBrushOrg(BrushOrg);
        memDC.SelectObject(BackGround);
        memDC.FillRect(CRect(0, 0, width + 1, height + 1), BackGround);

        ::DrawIconEx(memDC.m_hDC, 0, 0, TheIcon, width, height, 0, NULL, DI_NORMAL);
        
        if (TheIcon != m_hImage)
            ::DestroyIcon(TheIcon);
    }
    else if (m_ImageFlags & PJAI_BITMAP)
    {   // place bitmap image into the memory DC
        memDC.SelectObject((HBITMAP)m_hImage);

        if (m_TransparentColour == CLR_DEFAULT)
            m_TransparentColour = memDC.GetPixel(0, 0);

        if (m_DrawFlags & PJAI_STRETCHED)
        {   // stretch the image
            StretchDC.CreateCompatibleDC(pDC);
            savedStretchDC = StretchDC.SaveDC();
            stretchbmp.CreateCompatibleBitmap(pDC, w, h);
            StretchDC.SelectObject(stretchbmp);
            StretchDC.SetStretchBltMode(COLORONCOLOR);
            StretchDC.StretchBlt(0, 0, width, height, &memDC, 0, 0, m_size.cx, m_size.cy, SRCCOPY);
            pOutputDC = &StretchDC;
        }

        if (m_DrawFlags & PJAI_TRANSPARENT)
        {   // draw the image transparently
            TransparentDC.CreateCompatibleDC(pDC);
            savedTransparentDC = TransparentDC.SaveDC();
            Transparentbmp.CreateCompatibleBitmap(pDC, width, height);
            TransparentDC.SelectObject(&Transparentbmp);
            TransparentDC.SetBrushOrg(BrushOrg);
            TransparentDC.SelectObject(BackGround);
            TransparentDC.FillRect(CRect(0, 0, width + 1, height + 1), BackGround);
            DrawTransparent(&TransparentDC, width, height, pOutputDC);
            pOutputDC = &TransparentDC;
        }
        else if (m_DrawFlags & PJAI_GRAYSCALE)
        {   // convert the image to grayscale
            GrayBmp = GrayScale(pDC, *pOutputDC->GetCurrentBitmap());
            pOutputDC->SelectObject(GrayBmp);
        }
    }
    else
    {
        ASSERT (FALSE);  // m_Flags improperly set (should never get here)
    }

    if (m_DrawFlags & PJAI_DISABLED)  // draw the image disabled
        DitherBlt(pDC, left, top, width, height, pOutputDC);
    else   // draw the image
        pDC->BitBlt(left, top, width, height, pOutputDC, 0, 0, SRCCOPY);

    // clean up after ourselves
    if (savedTransparentDC)
    {
        TransparentDC.RestoreDC(savedTransparentDC);
        TransparentDC.DeleteDC();
    }

    if (savedStretchDC)
    {
        StretchDC.RestoreDC(savedStretchDC);
        StretchDC.DeleteDC();
    }

    memDC.RestoreDC(savedmemDC);
    memDC.DeleteDC();

    if (GrayBmp)
        ::DeleteObject(GrayBmp);
}

/////////////////////////////////////////////////////////////////////////////
//
//  CPJAImage::DrawTransparent  (protected member function)
//    transparently draws the image in the source device context onto the 
//    destination device context 
//
//  Parameters :
//    pToDC         [in] - pointer to the destination device context
//    w             [in] - the width of the image
//    h             [in] - the height of the image
//    pFromDC       [in] - pointer to the source DC containing the bitmap to be drawn
//
//  Returns :
//    Nothing
//
//  Note :
//    Uses the 'True Mask' method
//    Modified from code found at http://www.codeguru.com/bitmap/CISBitmap.shtml
//    original author Paul Reynolds (Paul.Reynolds@cmgroup.co.uk)
// 
/////////////////////////////////////////////////////////////////////////////

void CPJAImage::DrawTransparent(CDC *pToDC, int w, int h, CDC *pFromDC)
{
    // handle for the grayscale image
    HBITMAP Gray = NULL;

    CDC MonoDC;
    MonoDC.CreateCompatibleDC(pToDC);

    int savedToDC = pToDC->SaveDC();
    int savedFromDC = pFromDC->SaveDC();
    int savedMonoDC = MonoDC.SaveDC();

    pToDC->SetBkColor(RGB(255, 255, 255));
    pToDC->SetTextColor(RGB(0, 0, 0));
    pFromDC->SetBkColor(m_TransparentColour);

    // Create the mask
    CBitmap MonoDCbmp;
    MonoDCbmp.CreateBitmap(w, h, 1, 1, NULL);
    MonoDC.SelectObject(&MonoDCbmp);
    MonoDC.BitBlt(0, 0, w, h, pFromDC, 0, 0, SRCCOPY);

    if (m_DrawFlags & PJAI_GRAYSCALE)
    {   // convert the image to grayscale. We do this here just
        // in case the bitmap has a grayscale colour as the 
        // transparent colour
        Gray = GrayScale(pFromDC, *pFromDC->GetCurrentBitmap());
        pFromDC->SelectObject(Gray);
    }

    // draw the transparent bitmap
    pToDC->BitBlt(0, 0, w, h, pFromDC, 0, 0, SRCINVERT);
    pToDC->BitBlt(0, 0, w, h, &MonoDC, 0, 0, SRCAND);
    pToDC->BitBlt(0, 0, w, h, pFromDC, 0, 0, SRCINVERT);

    MonoDC.RestoreDC(savedMonoDC);
    pFromDC->RestoreDC(savedFromDC);
    pToDC->RestoreDC(savedToDC);
    MonoDC.DeleteDC();

    if (Gray)
        ::DeleteObject(Gray);
}

/////////////////////////////////////////////////////////////////////////////
//
//  CPJAImage::GetSize  (public member function)
//    Gets the size of the image in pixels
//
//  Parameters :
//    None
//
//  Returns :
//    a CSize containing the size of the image
//
/////////////////////////////////////////////////////////////////////////////

CSize CPJAImage::GetSize()
{
    return m_size;
}

/////////////////////////////////////////////////////////////////////////////
//
//  CPJAImage::Gray   (protected member function)
//    Gets the grayscale value of the given colour
//
//  Parameters :
//    r [in] - the red colour
//    g [in] - the green colour
//    b [in] - the blue colour
//
//  Returns :
//    the grayscale value
//
/////////////////////////////////////////////////////////////////////////////

int CPJAImage::Gray(int r, int g, int b)
{
    return (b * 11 + g * 59 + r * 30) / 100;
}

/////////////////////////////////////////////////////////////////////////////
//
//  CPJAImage::GrayScale   (protected member function)
//    Creates a DIBSection that is the grayscale equivalant of the given bitmap
//
//  Parameters :
//    pDC     [in] - a pointer to a DC used to build the colour table
//    hBitmap [in] - the bitmap to turn grayscale
//
//  Returns :
//    the HBITMAP of the new grayscale DIBSection
//
//  Note :
//    when you are finished with the HBITMAP returned by this function
//    you have to delete it with a call to ::DeleteObject()
//
//    This function does not change the transparent colour to grayscale.
//
/////////////////////////////////////////////////////////////////////////////

HBITMAP CPJAImage::GrayScale(CDC *pDC, HBITMAP hBitmap)
{
    // get the bitmap's size and colour information
    BITMAP bm;
    ::GetObject(hBitmap, sizeof(BITMAP), &bm);

    // create a DIBSection copy of the original bitmap
    HBITMAP hDib = (HBITMAP)::CopyImage(hBitmap, IMAGE_BITMAP, 0, 0, LR_COPYRETURNORG | LR_CREATEDIBSECTION);
    
    if (bm.bmBitsPixel < 16)
    {   // bitmap has a colour table, so we modify the colour table
        CDC memDC;
        memDC.CreateCompatibleDC(pDC);
        int SavedMemDC = memDC.SaveDC();
        memDC.SelectObject(hDib);
        int nColours = 1 << bm.bmBitsPixel;

        RGBQUAD pal[256];

        // Get the colour table
        ::GetDIBColorTable(memDC.m_hDC, 0, nColours, pal);

        // modify the colour table
        for (int x = 0; x < nColours; x++)
        {
            if (!IsTransparent(pal[x].rgbRed, pal[x].rgbGreen, pal[x].rgbBlue))
            {
                BYTE nGray = (BYTE)Gray(pal[x].rgbRed, pal[x].rgbGreen, pal[x].rgbBlue);
                pal[x].rgbRed = nGray;
                pal[x].rgbGreen = nGray;
                pal[x].rgbBlue = nGray;
            }
        }

        // set the modified colour tab to the DIBSection bitmap
        ::SetDIBColorTable(memDC.m_hDC, 0, nColours, pal);
        
        memDC.RestoreDC(SavedMemDC);
        memDC.DeleteDC();
        return hDib;
    }

    else
    {   // the bitmap does not have a colour table, so we modify the bitmap bits directly
        int Size = bm.bmHeight * bm.bmWidth;
        
        BITMAPINFO bmi;
        bmi.bmiHeader.biSize          = sizeof(BITMAPINFOHEADER);
        bmi.bmiHeader.biHeight        = bm.bmHeight;
        bmi.bmiHeader.biWidth         = bm.bmWidth;
        bmi.bmiHeader.biPlanes        = 1;
        bmi.bmiHeader.biBitCount      = bm.bmBitsPixel;
        bmi.bmiHeader.biCompression   = BI_RGB;
        bmi.bmiHeader.biSizeImage     = ((bm.bmWidth * bm.bmBitsPixel + 31) & (~31)) / 8 * bm.bmHeight;
        bmi.bmiHeader.biXPelsPerMeter = 0;
        bmi.bmiHeader.biYPelsPerMeter = 0;
        bmi.bmiHeader.biClrUsed       = 0;
        bmi.bmiHeader.biClrImportant  = 0;
        
        // Get the bitmaps data bits
        BYTE *pBits = new BYTE[bmi.bmiHeader.biSizeImage];
        VERIFY (::GetDIBits(pDC->m_hDC, hDib, 0, bm.bmHeight, pBits, &bmi, DIB_RGB_COLORS));

        if (bm.bmBitsPixel == 32)
        {
            DWORD *dst=(DWORD *)pBits;
            
            while (Size--)
            {
                if (!IsTransparent(GetBValue(*dst), GetGValue(*dst), GetRValue(*dst)))
                {
                    int nGray = Gray(GetBValue(*dst), GetGValue(*dst), GetRValue(*dst));
                    *dst = (DWORD)RGB(nGray, nGray, nGray);
                }
                dst++;
           }
        }
        
        else if (bm.bmBitsPixel == 24)
        {
            BYTE *dst=(BYTE*)pBits;

            for (int dh = 0; dh < bm.bmHeight; dh++)
            {
                for (int dw = 0; dw < bm.bmWidth; dw++)
                {
                    if (!IsTransparent(dst[2], dst[1], dst[0]))
                    {
                        int nGray = Gray(dst[2], dst[1], dst[0]);
                    
                        dst[0]=(BYTE)nGray;
                        dst[1]=(BYTE)nGray;
                        dst[2]=(BYTE)nGray;
                    }
                    dst += 3;
                }

                // each row is DWORD aligned, so when we reach the end of a row, we
                // have to realign the pointer to point to the start of the next row
                int pos = (int)dst - (int)pBits;
                int rem = pos % 4;
                if (rem)
                    dst += 4 - rem;
            }
        }

        else if (bm.bmBitsPixel == 16)
        {
            WORD *dst=(WORD*)pBits;
            
            while (Size--)
            {
                BYTE b = (BYTE)((*dst)&(0x1F));
                BYTE g = (BYTE)(((*dst)>>5)&(0x1F));
                BYTE r = (BYTE)(((*dst)>>10)&(0x1F));
                
                if (!IsTransparent(r, g, b))
                {
                    int nGray = Gray(r, g, b);
                    *dst = ((WORD)(((BYTE)(nGray)|((WORD)((BYTE)(nGray))<<5))|(((DWORD)(BYTE)(nGray))<<10)));
                }
                dst++;
            }
        }

        // set the modified bitmap data bits to the DIBSection
        ::SetDIBits(pDC->m_hDC, hDib, 0, bm.bmHeight, pBits, &bmi, DIB_RGB_COLORS);
        delete[] pBits;
        return hDib;
    }
}

/////////////////////////////////////////////////////////////////////////////
//
//  CPJAImage::IsTransparent  (protected member function)
//    checks if the given colour values make up the transparent colour for the bitmap
//
//  Parameters :
//    r [in] - the red colour
//    g [in] - the green colour
//    b [in] - the blue colour
//
//  Returns :
//    TRUE if the colour is the transparent colour
//    FALSE if not
//
/////////////////////////////////////////////////////////////////////////////


BOOL CPJAImage::IsTransparent(int r, int g, int b)
{
    if (!(m_DrawFlags & PJAI_TRANSPARENT))
        return FALSE;
    return  (RGB(r, g, b) == m_TransparentColour);
}

/////////////////////////////////////////////////////////////////////////////
//
//  CPJAImage::SetImage  (public member function)
//    Sets the image and image flags
//
//  Parameters :
//    Image [in] - a HANDLE of the image to set (either a HBITMAP or a HICON)
//    Flags [in] - the flags that specify the type of image and how it is drawn
//
//  Returns :
//    TRUE on success
//    FALSE on failure
//
//  Note :
//    See the PJAImage.h file for a description of the flags used
//
/////////////////////////////////////////////////////////////////////////////

BOOL CPJAImage::SetImage(HANDLE Image, DWORD Flags)
{
    if (Image)
    {   // verify flags
        if (!((Flags & PJAI_BITMAP ? 1 : 0) ^ (Flags & PJAI_ICON ? 1 : 0)))
        {
            TRACE (_T("PJAI_Image::SetImage() : Must specify either PJAI_BITMAP or PJAI_ICON"));
            ASSERT (FALSE);
            return FALSE;
        }
    }

    if (m_hImage && m_hImage != Image)
    {
        if (m_ImageFlags & PJAI_AUTODELETE)
        {   // remove the old image
            if (m_ImageFlags & PJAI_ICON)
                DestroyIcon((HICON)m_hImage);
            else
                DeleteObject((HGDIOBJ)m_hImage);
        }
        m_hImage = NULL;
        m_size.cx = 0;
        m_size.cy = 0;
        m_TransparentColour = CLR_DEFAULT;
    }

    if (Image)
    {   // get the image dimensions
        if (Flags & PJAI_BITMAP)
        {
            BITMAP bmp;
            if (GetObject((HBITMAP)Image, sizeof(BITMAP), &bmp))
            {
                m_size.cx = bmp.bmWidth;
                m_size.cy = bmp.bmHeight;
            }
        }
        else if (Flags & PJAI_ICON)
        {
            ICONINFO iconinfo;
            GetIconInfo((HICON)Image, &iconinfo);
            BITMAP bmp;
            if (GetObject(iconinfo.hbmMask, sizeof(BITMAP), &bmp))
            {
                m_size.cx = bmp.bmWidth;
                m_size.cy = iconinfo.hbmColor ? bmp.bmHeight : bmp.bmHeight / 2;
            }
            // prevent a resource leak
            DeleteObject(iconinfo.hbmColor);
            DeleteObject(iconinfo.hbmMask);
        }
    }

    m_hImage = Image;
    m_ImageFlags = Flags;
    return TRUE;
}

/////////////////////////////////////////////////////////////////////////////
//
//  CPJAImage::SetTransparentColour  (public member function)
//    Set the colour to be used as the transparent colour 
//
//  Parameters :
//    clr [in] - the colour to be used as the transparent colour
//
//  Returns :
//    the old transparent colour
//
//  Note :
//    If the colour is CLR_DEFAULT (the default), the colour of the
//    top left pixel (0,0) is used as the transparent colour.
//
/////////////////////////////////////////////////////////////////////////////

COLORREF CPJAImage::SetTransparentColour(COLORREF clr /*= CLR_DEFAULT*/)
{
    COLORREF oldclr = m_TransparentColour;
    m_TransparentColour = clr;
    return oldclr;
}

