/****************************************************************************
PJAImage.h : header file for the PJAImage class
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

#if !defined(AFX_PJAIMAGE_H__F15965B0_B05A_11D4_B625_A1459D96AB20__INCLUDED_)
#define AFX_PJAIMAGE_H__F15965B0_B05A_11D4_B625_A1459D96AB20__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

/////////////////////////////////////////////////////////////////////////////
// CPJAImage class

// this class is used to draw an image (HBITMAP or HICON)

// Flags used by the SetImage() function
//
// Usage :
// PJAI_ICON       - the image is an icon, cannot be used with FEC_BITMAP
// PJAI_BITMAP     - the image is a bitmap, cannot be used with FEC_ICON
// PJAI_AUTODELETE - the image handle is deleted and the memory freed when it is no longer needed
#define PJAI_ICON        0x00000001
#define PJAI_BITMAP      0x00000002
#define PJAI_AUTODELETE  0x00000004

// Flags used by the DrawImage() function
//
// Usage :
// PJAI_TRANSPARENT - used with FEC_BITMAP, the bitmap will be drawn transparently
// PJAI_STRETCHED   - Can not be used with PJAI_CENTERED . The image will be stretched to fit the rectangle in DrawImage()
// PJAI_CENTERED    - Can not be used with PJAI_STRETCHED. The image will be centered on the rectangle
// PJAI_DISABLED    - Can not be used with PJAI_GRAYSCALE. Draws the image as disabled (3D Monochrome)
// PJAI_GRAYSCALE   - Can not be used with PJAI_DISABLED. Draws the image in grayscale
#define PJAI_TRANSPARENT 0x00000010
#define PJAI_STRETCHED   0x00000020
#define PJAI_CENTERED    0x00000040
#define PJAI_DISABLED    0x00000080
#define PJAI_GRAYSCALE   0x00000100

class CPJAImage  
{
public:
    CPJAImage();                        // constructor
    virtual ~CPJAImage();               // destructor

    void DrawImage(CDC *pDC, int x, int y, int w, int h, DWORD DrawFlags = 0); // draws the image
    CSize GetSize();                          // retreives the dimensions of the image
    BOOL SetImage(HANDLE Image, DWORD Flags); // setup the image
    COLORREF SetTransparentColour(COLORREF clr = CLR_DEFAULT); // set the colour to used as the transparent colour

protected:
    void DitherBlt(CDC *pToDC, int x, int y, int w, int h, CDC *pFromDC);  // draw the bitmap grayed
    void DrawTransparent(CDC *pToDC, int w, int h, CDC *pFromDC); // draw the bitmap transparently
    int Gray(int r, int g, int b);
	HBITMAP GrayScale (CDC *pDC, HBITMAP hBitmap);
	BOOL IsTransparent(int r, int g, int b);

private:
    DWORD m_ImageFlags;                 // the control flags
    DWORD m_DrawFlags;                  // the drawing flags
    HANDLE m_hImage;                    // the image
    CSize m_size;                       // the image size
    COLORREF m_TransparentColour;       // the transparent colour
};

#endif