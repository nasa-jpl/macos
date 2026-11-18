/* giza - a scientific plotting library built on cairo
 *
 * Copyright (c) 2010      James Wetter and Daniel Price
 * Copyright (c) 2010-2012 Daniel Price
 *
 * This library is free software; and you are welcome to redistribute
 * it under the terms of the GNU General Public License
 * (GPL, see LICENSE file for details) and the provision that
 * this notice remains intact. If you modify this file, please
 * note section 2a) of the GPLv2 states that:
 *
 *  a) You must cause the modified files to carry prominent notices
 *     stating that you changed the files and the date of any change.
 *
 * This software is distributed "AS IS", with ABSOLUTELY NO WARRANTY.
 * See the GPL for specific language governing rights and limitations.
 *
 * The Original code is the giza plotting library.
 *
 * Contributor(s):
 *      Added Win32 windowed driver - April 2026
 *
 * Description:
 *   Native Win32 windowed display driver for giza.
 *   Uses a cairo image surface (ARGB32) as the backing store.
 *   The window runs a message loop in a separate thread so MACOS
 *   does not block while the window is open.
 *   Cursor/keyboard interaction is handled via shared state protected
 *   by a CRITICAL_SECTION.
 */

#ifdef _WIN32

#include "giza-private.h"
#include "giza-drivers-private.h"
#include "giza-driver-win-private.h"
#include "giza-io-private.h"
#include "giza-transforms-private.h"
#include <giza.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

/* ===========================================================================
 * Constants
 * =========================================================================*/
#define GIZA_WIN_CLASS        "GizaWindow"
#define GIZA_WIN_TITLE        "Giza"
#define GIZA_DEFAULT_WIDTH    800
#define GIZA_DEFAULT_HEIGHT   600
#define GIZA_DEVICE_UNITS_PER_MM    3.7054
#define GIZA_DEVICE_UNITS_PER_PIXEL 1.0
#define GIZA_DEVICE_INTERACTIVE     1
#define GIZA_WIN_MARGIN             20

/* Custom window messages */
#define GIZA_WM_PAINT_REQUEST  (WM_USER + 1)
#define GIZA_WM_CLOSE_REQUEST  (WM_USER + 2)

/* ===========================================================================
 * Per-device state
 * =========================================================================*/
struct GIZA_WinDevice
{
    HWND            hwnd;           /* Window handle                        */
    HANDLE          thread;         /* Message loop thread                  */
    DWORD           threadId;       /* Thread ID                            */
    CRITICAL_SECTION cs;            /* Protects shared state below          */

    /* Keyboard/mouse input — written by WndProc, read by get_key_press    */
    volatile int    gotInput;       /* Non-zero when input is ready         */
    volatile char   inputChar;      /* Character pressed                    */
    volatile double inputX;         /* Cursor X in device coords            */
    volatile double inputY;         /* Cursor Y in device coords            */

    int             width;          /* Window client area width  (pixels)   */
    int             height;         /* Window client area height (pixels)   */
    int             in_use;
};

static struct GIZA_WinDevice WD[GIZA_MAX_DEVICES];
static int giza_win_class_registered = 0;

/* ===========================================================================
 * Forward declarations
 * =========================================================================*/
static LRESULT CALLBACK _giza_wndproc (HWND hwnd, UINT msg,
                                        WPARAM wp, LPARAM lp);
static DWORD WINAPI     _giza_win_thread (LPVOID param);
static void             _giza_win_blit   (int devid);

/* ===========================================================================
 * Helper: blit the cairo image surface to the window via GDI
 * =========================================================================*/
static void
_giza_win_blit (int devid)
{
    if (!WD[devid].hwnd || !Dev[devid].surface)
        return;

    int w = Dev[devid].width;
    int h = Dev[devid].height;

    /* Get raw pixel data from cairo surface */
    cairo_surface_flush (Dev[devid].surface);
    unsigned char *data = cairo_image_surface_get_data (Dev[devid].surface);
    int stride          = cairo_image_surface_get_stride (Dev[devid].surface);

    /* Build a BITMAPINFO for a 32-bit top-down DIB */
    BITMAPINFO bmi;
    ZeroMemory (&bmi, sizeof (bmi));
    bmi.bmiHeader.biSize        = sizeof (BITMAPINFOHEADER);
    bmi.bmiHeader.biWidth       = w;
    bmi.bmiHeader.biHeight      = -h;   /* negative = top-down             */
    bmi.bmiHeader.biPlanes      = 1;
    bmi.bmiHeader.biBitCount    = 32;
    bmi.bmiHeader.biCompression = BI_RGB;

    HDC hdc = GetDC (WD[devid].hwnd);
    if (hdc)
    {
        /* Centre the drawing area inside the margin */
        SetDIBitsToDevice (hdc,
                           GIZA_WIN_MARGIN, GIZA_WIN_MARGIN, /* dest x,y  */
                           w, h,                              /* dest w,h  */
                           0, 0,                              /* src x,y   */
                           0, h,                              /* scan line */
                           data, &bmi, DIB_RGB_COLORS);
        ReleaseDC (WD[devid].hwnd, hdc);
    }
}

/* ===========================================================================
 * Window procedure — runs on the message-loop thread
 * =========================================================================*/
static LRESULT CALLBACK
_giza_wndproc (HWND hwnd, UINT msg, WPARAM wp, LPARAM lp)
{
    /* Find which device owns this window */
    int devid = -1;
    for (int i = 0; i < GIZA_MAX_DEVICES; i++)
        if (WD[i].hwnd == hwnd) { devid = i; break; }

    switch (msg)
    {
    case WM_PAINT:
    {
        PAINTSTRUCT ps;
        BeginPaint (hwnd, &ps);
        if (devid >= 0)
            _giza_win_blit (devid);
        EndPaint (hwnd, &ps);
        return 0;
    }

    case GIZA_WM_PAINT_REQUEST:
        if (devid >= 0)
            _giza_win_blit (devid);
        return 0;

    case WM_KEYDOWN:
    case WM_CHAR:
    {
        if (devid < 0) break;
        /* Convert virtual key / char to a single ASCII character */
        char ch = 0;
        if (msg == WM_CHAR)
            ch = (char) wp;
        else
        {
            /* Handle special keys as single chars similar to X11 */
            switch (wp)
            {
            case VK_RETURN: ch = '\r'; break;
            case VK_ESCAPE: ch = '\x1b'; break;
            case VK_LEFT:   ch = 'D'; break;   /* PGPLOT cursor keys */
            case VK_RIGHT:  ch = 'C'; break;
            case VK_UP:     ch = 'A'; break;
            case VK_DOWN:   ch = 'B'; break;
            default:        ch = 0;   break;
            }
        }
        if (ch != 0)
        {
            /* Record cursor position in device coords */
            POINT pt;
            GetCursorPos (&pt);
            ScreenToClient (hwnd, &pt);

            EnterCriticalSection (&WD[devid].cs);
            WD[devid].inputChar = ch;
            WD[devid].inputX    = (double)(pt.x - GIZA_WIN_MARGIN);
            WD[devid].inputY    = (double)(WD[devid].height
                                           - (pt.y - GIZA_WIN_MARGIN));
            WD[devid].gotInput  = 1;
            LeaveCriticalSection (&WD[devid].cs);
        }
        return 0;
    }

    case WM_LBUTTONDOWN:
    {
        if (devid < 0) break;
        int mx = LOWORD (lp) - GIZA_WIN_MARGIN;
        int my = WD[devid].height - (HIWORD (lp) - GIZA_WIN_MARGIN);
        EnterCriticalSection (&WD[devid].cs);
        WD[devid].inputChar = 'A';   /* left click = 'A' like X11 */
        WD[devid].inputX    = (double) mx;
        WD[devid].inputY    = (double) my;
        WD[devid].gotInput  = 1;
        LeaveCriticalSection (&WD[devid].cs);
        return 0;
    }

    case WM_RBUTTONDOWN:
    {
        if (devid < 0) break;
        int mx = LOWORD (lp) - GIZA_WIN_MARGIN;
        int my = WD[devid].height - (HIWORD (lp) - GIZA_WIN_MARGIN);
        EnterCriticalSection (&WD[devid].cs);
        WD[devid].inputChar = 'X';   /* right click = 'X' */
        WD[devid].inputX    = (double) mx;
        WD[devid].inputY    = (double) my;
        WD[devid].gotInput  = 1;
        LeaveCriticalSection (&WD[devid].cs);
        return 0;
    }

    case GIZA_WM_CLOSE_REQUEST:
        DestroyWindow (hwnd);
        return 0;

    case WM_DESTROY:
        if (devid >= 0)
            WD[devid].hwnd = NULL;
        PostQuitMessage (0);
        return 0;
    }
    return DefWindowProc (hwnd, msg, wp, lp);
}

/* ===========================================================================
 * Message loop thread
 * =========================================================================*/
static DWORD WINAPI
_giza_win_thread (LPVOID param)
{
    int devid = (int)(intptr_t) param;

    /* Register window class once */
    if (!giza_win_class_registered)
    {
        WNDCLASSEX wc;
        ZeroMemory (&wc, sizeof (wc));
        wc.cbSize        = sizeof (wc);
        wc.style         = CS_HREDRAW | CS_VREDRAW;
        wc.lpfnWndProc   = _giza_wndproc;
        wc.hInstance     = GetModuleHandle (NULL);
        wc.hCursor       = LoadCursor (NULL, IDC_ARROW);
        wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
        wc.lpszClassName = GIZA_WIN_CLASS;
        RegisterClassEx (&wc);
        giza_win_class_registered = 1;
    }

    int winW = WD[devid].width  + 2 * GIZA_WIN_MARGIN;
    int winH = WD[devid].height + 2 * GIZA_WIN_MARGIN;

    /* Account for title bar / borders */
    RECT rc = { 0, 0, winW, winH };
    AdjustWindowRect (&rc, WS_OVERLAPPEDWINDOW, FALSE);

    WD[devid].hwnd = CreateWindowEx (
        0,
        GIZA_WIN_CLASS,
        GIZA_WIN_TITLE,
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT,
        rc.right - rc.left,
        rc.bottom - rc.top,
        NULL, NULL,
        GetModuleHandle (NULL),
        NULL);

    if (!WD[devid].hwnd)
    {
        _giza_error ("_giza_win_thread", "Could not create window");
        return 1;
    }

    ShowWindow  (WD[devid].hwnd, SW_SHOW);
    UpdateWindow (WD[devid].hwnd);

    /* Standard Win32 message loop */
    MSG msg;
    while (GetMessage (&msg, NULL, 0, 0))
    {
        TranslateMessage (&msg);
        DispatchMessage  (&msg);
    }
    return 0;
}

/* ===========================================================================
 * Driver functions
 * =========================================================================*/

int
_giza_open_device_win (double width, double height, int units)
{
    /* Initialise device array on first use */
    static int didInit = 0;
    if (!didInit)
    {
        memset (&WD[0], 0, GIZA_MAX_DEVICES * sizeof (struct GIZA_WinDevice));
        didInit = 1;
    }

    if (WD[id].in_use)
    {
        _giza_error ("_giza_open_device_win",
                     "Internal error: WD[%d] still in use", id);
        return -1;
    }

    memset (&WD[id], 0, sizeof (struct GIZA_WinDevice));
    WD[id].in_use = 1;
    InitializeCriticalSection (&WD[id].cs);

    Dev[id].deviceUnitsPermm    = GIZA_DEVICE_UNITS_PER_MM;
    Dev[id].deviceUnitsPerPixel = GIZA_DEVICE_UNITS_PER_PIXEL;
    Dev[id].isInteractive       = GIZA_DEVICE_INTERACTIVE;

    if (width > 0. && height > 0. && units > 0)
        _giza_get_specified_size (width, height, units,
                                  &Dev[id].width, &Dev[id].height);
    else
    {
        Dev[id].width  = GIZA_DEFAULT_WIDTH;
        Dev[id].height = GIZA_DEFAULT_HEIGHT;
    }

    WD[id].width  = Dev[id].width;
    WD[id].height = Dev[id].height;

    /* Create cairo image surface as backing store */
    Dev[id].surface = cairo_image_surface_create (
        CAIRO_FORMAT_ARGB32, Dev[id].width, Dev[id].height);
    if (!Dev[id].surface)
    {
        _giza_error ("_giza_open_device_win",
                     "Could not create cairo image surface");
        return -1;
    }

    /* Launch window thread */
    WD[id].thread = CreateThread (NULL, 0, _giza_win_thread,
                                   (LPVOID)(intptr_t) id, 0,
                                   &WD[id].threadId);
    if (!WD[id].thread)
    {
        _giza_error ("_giza_open_device_win", "Could not create window thread");
        cairo_surface_destroy (Dev[id].surface);
        Dev[id].surface = NULL;
        return -1;
    }

    /* Wait for window to be created (up to 5 seconds) */
    for (int i = 0; i < 500 && !WD[id].hwnd; i++)
        Sleep (10);

    if (!WD[id].hwnd)
    {
        _giza_error ("_giza_open_device_win", "Window did not appear");
        return -1;
    }

    return 0;
}

/* ---------------------------------------------------------------------------
 * Flush: blit the backing surface to the window
 * -------------------------------------------------------------------------*/
void
_giza_flush_device_win (void)
{
    if (WD[id].hwnd)
        PostMessage (WD[id].hwnd, GIZA_WM_PAINT_REQUEST, 0, 0);
}

/* ---------------------------------------------------------------------------
 * Change page: clear the surface and repaint
 * -------------------------------------------------------------------------*/
void
_giza_change_page_win (void)
{
    if (!Dev[id].surface) return;

    cairo_destroy (Dev[id].context);

    cairo_surface_destroy (Dev[id].surface);
    Dev[id].surface = cairo_image_surface_create (
        CAIRO_FORMAT_ARGB32, Dev[id].width, Dev[id].height);
    if (!Dev[id].surface)
    {
        _giza_error ("_giza_change_page_win",
                     "Could not recreate cairo surface");
        return;
    }
    Dev[id].context = cairo_create (Dev[id].surface);
    if (!Dev[id].context)
    {
        _giza_error ("_giza_change_page_win",
                     "Could not recreate cairo context");
        return;
    }

    /* Clear to background colour */
    cairo_set_source_rgb (Dev[id].context, 1.0, 1.0, 1.0);
    cairo_paint (Dev[id].context);

    _giza_flush_device_win ();
}

/* ---------------------------------------------------------------------------
 * Close: destroy window and release resources
 * -------------------------------------------------------------------------*/
void
_giza_close_device_win (void)
{
    if (WD[id].hwnd)
    {
        /* Ask the window to close itself from its own thread */
        PostMessage (WD[id].hwnd, GIZA_WM_CLOSE_REQUEST, 0, 0);

        /* Wait for thread to finish (up to 3 seconds) */
        WaitForSingleObject (WD[id].thread, 3000);
        CloseHandle (WD[id].thread);
        WD[id].thread = NULL;
    }

    if (Dev[id].surface)
    {
        cairo_surface_destroy (Dev[id].surface);
        Dev[id].surface = NULL;
    }

    DeleteCriticalSection (&WD[id].cs);
    WD[id].in_use = 0;
}

/* ---------------------------------------------------------------------------
 * Expand clipping: match XW behaviour (reset to full device)
 * -------------------------------------------------------------------------*/
void
_giza_expand_clipping_win (void)
{
    _giza_set_trans (GIZA_TRANS_IDEN);
    cairo_reset_clip (Dev[id].context);
    cairo_rectangle (Dev[id].context,
                     0., 0.,
                     (double) Dev[id].width,
                     (double) Dev[id].height);
    cairo_clip (Dev[id].context);
}

/* ---------------------------------------------------------------------------
 * Get key press: block until user presses a key or clicks the mouse.
 * x, y returned in world coordinates via the normal giza transform.
 * -------------------------------------------------------------------------*/
void
_giza_get_key_press_win (int mode, int moveCurs, int nanc,
                          const double *xanc, const double *yanc,
                          double *x, double *y, char *ch)
{
    if (!WD[id].hwnd)
    {
        *ch = 'Q';
        return;
    }

    /* Flush any pending drawing first */
    _giza_flush_device_win ();

    /* Bring window to front and give it focus */
    SetForegroundWindow (WD[id].hwnd);
    SetFocus (WD[id].hwnd);

    /* Clear any previous input */
    EnterCriticalSection (&WD[id].cs);
    WD[id].gotInput = 0;
    LeaveCriticalSection (&WD[id].cs);

    /* Spin-wait for input — yields the CPU via Sleep(1) */
    int got = 0;
    while (!got)
    {
        Sleep (10);
        EnterCriticalSection (&WD[id].cs);
        got = WD[id].gotInput;
        if (got)
        {
            *ch = WD[id].inputChar;
            /* Convert device coords back to world coords */
            double dx = WD[id].inputX;
            double dy = WD[id].inputY;
            /* Apply inverse of current normalisation transform */
            _giza_set_trans (GIZA_TRANS_NORM);
            cairo_device_to_user (Dev[id].context, &dx, &dy);
            *x = dx;
            *y = dy;
            WD[id].gotInput = 0;
        }
        LeaveCriticalSection (&WD[id].cs);

        /* If window was closed, return 'Q' */
        if (!WD[id].hwnd)
        {
            *ch = 'Q';
            return;
        }
    }
}

/* ---------------------------------------------------------------------------
 * Init band: rubber-band drawing not yet implemented — return 0 (success)
 * -------------------------------------------------------------------------*/
int
_giza_init_band_win (void)
{
    return 0;
}

/* ---------------------------------------------------------------------------
 * Select: make the given device current (no extra work needed)
 * -------------------------------------------------------------------------*/
int
_giza_select_win (int devid)
{
    return 0;
}

#endif /* _WIN32 */
