/*
 *        Polygon Reduction Demo by Stan Melax (c) 1998
 *  Permission to use any of this code wherever you want is granted..
 *  Although, please do acknowledge authorship if appropriate.
 *
 *  This module contains the window setup code, mouse input, timing
 *  routines, and that sort of stuff.  The interesting modules
 *  to see are bunnygut.cpp and progmesh.cpp.
 *
 *  The windows 95 specific code for this application was taken from
 *  an example of processing mouse events in an OpenGL program using
 *  the Win32 API from the www.opengl.org web site.
 *
 *  Under Project->Settings, Link Options, General Category
 *  Add:
 *        Opengl32.lib glu32.lib winmm.lib
 *  to the Object/Library Modules
 *
 *  You will need have OpenGL libs and include files to compile this
 *  Go to the www.opengl.org web site if you need help with this.
 */


#include <windows.h>    /* must include this before GL/gl.h */
#include <GL/gl.h>              /* OpenGL header file */
#include <GL/glu.h>             /* OpenGL utilities header file */
#include <stdio.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include <time.h>

#include "vector.h"
#include "font.h"

// Functions and Variables from bunny module
extern void       InitModel();
extern void       RenderModel();
extern Vector     model_position;      // position of bunny
extern Quaternion model_orientation;   // orientation of bunny

// Global Variables
float   DeltaT = 0.1f;
float   FPS;
int     Width  = 512;
int     Height = 512;
int     MouseX = 0;
int     MouseY = 0;
Vector  MouseVector;      // 3D direction mouse points
Vector  OldMouseVector;
int     MouseState=0;     // true iff left button down
float   ViewAngle=45.0f;

HDC hDC;                /* device context */
HPALETTE hPalette = 0;  /* custom palette (if needed) */


void CalcFPSDeltaT(){
        static int timeinit=0;
        static int start,start2,current,last;
        static int frame=0, frame2=0;
        if(!timeinit){
                frame=0;
                start=timeGetTime();
                timeinit=1;
        }
        frame++;
        frame2++;
        current=timeGetTime(); // found in winmm.lib
        double dif=(double)(current-start)/CLOCKS_PER_SEC;
        double rv = (dif)? (double)frame/(double)dif:-1.0;
        if(dif>2.0 && frame >10) {
                start  = start2;
                frame  = frame2;
                start2 = timeGetTime();
                frame2 = 0;
        }
        DeltaT = (float)(current-last)/CLOCKS_PER_SEC;
        if(current==last) {
                DeltaT = 0.1f/CLOCKS_PER_SEC;  // it just can't be 0
        }
        // if(DeltaT>1.0) DeltaT=1.0;
        FPS = (float)rv;
        last = current;
}


void ComputeMouseVector(){
        OldMouseVector=MouseVector;
        float spread = (float)tan(ViewAngle/2*3.14/180);
        float y = spread * ((Height-MouseY)-Height/2.0f) /(Height/2.0f);
    float x = spread * (MouseX-Width/2.0f)  /(Height/2.0f);
    Vector v(x ,y,-1);
    // v=UserOrientation *v;
    v=normalise(v);
        MouseVector = v;
}

Quaternion VirtualTrackBall(Vector cop,Vector cor,Vector dir1,Vector dir2) {
        // Implement track ball functionality to spin stuf on the screen
        //  cop   center of projection
        //  cor   center of rotation
        //  dir1  old mouse direction
        //  dir2  new mouse direction
        // pretend there is a sphere around cor.  Then find the points
        // where dir1 and dir2 intersect that sphere.  Find the
        // rotation that takes the first point to the second.
        float m;
        // compute plane
        Vector nrml = cor - cop;
         // since trackball proportional to distance from cop
        float fudgefactor = 1.0f/(magnitude(nrml) * 0.25f);
        nrml = normalise(nrml);
        float dist = -(nrml^cor);
        Vector u= planelineintersection(nrml,dist,cop,cop+dir1);
        u=u-cor;
        u=u*fudgefactor;
        m= magnitude(u);
        if(m>1) {u=u*1.0f/m;}
        else {
                u=u - (nrml * (float)sqrt(1-m*m));
        }
        Vector v= planelineintersection(nrml,dist,cop,cop+dir2);
        v=v-cor;
        v=v*fudgefactor;
        m= magnitude(v);
        if(m>1) {v=v*1.0f/m;}
        else {
                v=v - (nrml * (float)sqrt(1-m*m));
        }
        Vector axis = u*v;
        float angle;
        m=magnitude(axis);
        if(m>1)m=1; // avoid potential floating point error
        Quaternion q(Vector(1.0f,0.0f,0.0f),0.0f);
        if(m>0 && (angle=(float)asin(m))>3.14/180) {
                        axis = normalise(axis);
                        q=Quaternion(axis,angle);
        }
        return q;
}

void SpinIt(){
        // Change the orientation of the bunny according to mouse drag
        Quaternion q=VirtualTrackBall(Vector(0,0,0),model_position,
                                      OldMouseVector,MouseVector);
        model_orientation=q*model_orientation;
}

void Reshape(int width, int height){
        // called initially and when the window changes size
        Width=width;
        Height=height;
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(ViewAngle, (float)width/height, 0.1, 50.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void PrintStats(){
        char buf[1024];buf[0]='\0';
        sprintf(buf,"FPS: %5.2f   ",FPS);
        PostString(buf,0,-1,0);
}

void Display(){
        // main drawing routine - called every frame
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushMatrix();
        glLoadIdentity();
        // camera at default (zero) position and orientation
        RenderModel();
        PrintStats();
        glLoadIdentity();
        RenderStrings();
    glPopMatrix();
    glFlush();
    SwapBuffers(hDC);                   /* nop if singlebuffered */
}


LONG WINAPI WindowProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    static PAINTSTRUCT ps;
    static GLboolean left  = GL_FALSE;  /* left button currently down? */
    static GLboolean right = GL_FALSE;  /* right button currently down? */
    static int omx, omy, mx, my;

    switch(uMsg) {
    case WM_PAINT:
                BeginPaint(hWnd, &ps);
                EndPaint(hWnd, &ps);
                return 0;
    case WM_SIZE:
                Reshape(LOWORD(lParam), HIWORD(lParam));
                PostMessage(hWnd, WM_PAINT, 0, 0);
                return 0;
    case WM_CHAR:
                switch (wParam) {
                        case 27: /* ESC key */
                            PostQuitMessage(0);
                            break;
                }
                return 0;

        case WM_LBUTTONDOWN:
            /* if we don't set the capture we won't get mouse move
           messages when the mouse moves outside the window. */
                SetCapture(hWnd);
                MouseX = LOWORD(lParam);
                MouseY = HIWORD(lParam);
                ComputeMouseVector();
                MouseState = 1;
                return 0;

    case WM_LBUTTONUP:
                MouseX = LOWORD(lParam);
                MouseY = HIWORD(lParam);
            if(MouseX & 1 << 15) MouseX -= (1 << 16);
            if(MouseY & 1 << 15) MouseY -= (1 << 16);
                ComputeMouseVector();
                if(MouseState) SpinIt();
                MouseState=0;
                /* remember to release the capture when we are finished. */
                ReleaseCapture();
                return 0;

    case WM_MOUSEMOVE:
                MouseX = LOWORD(lParam);
                MouseY = HIWORD(lParam);
            /* Win32 is pretty braindead about the x, y position that
               it returns when the mouse is off the left or top edge
               of the window (due to them being unsigned). therefore,
               roll the Win32's 0..2^16 pointer co-ord range to the
               more amenable (and useful) 0..+/-2^15. */
            if(MouseX & 1 << 15) MouseX -= (1 << 16);
            if(MouseY & 1 << 15) MouseY -= (1 << 16);
                ComputeMouseVector();
                if(MouseState) SpinIt();
                return 0;

    case WM_PALETTECHANGED:
                if (hWnd == (HWND)wParam) break;
                /* fall through to WM_QUERYNEWPALETTE */
    case WM_QUERYNEWPALETTE:
                if (hPalette) {
                    UnrealizeObject(hPalette);
                    SelectPalette(hDC, hPalette, FALSE);
                    RealizePalette(hDC);
                    return TRUE;
                }
                return FALSE;

    case WM_CLOSE:
                PostQuitMessage(0);
                return 0;
    }
    return DefWindowProc(hWnd, uMsg, wParam, lParam);
}

HWND CreateOpenGLWindow(char* title)
{
        // make a double-buffered, rgba, opengl window
    int         n, pf;
    HWND        hWnd;
    WNDCLASS    wc;
    LOGPALETTE* lpPal;
    PIXELFORMATDESCRIPTOR pfd;
    static HINSTANCE hInstance = 0;

    /* only register the window class once - use hInstance as a flag. */
    if (!hInstance) {
                hInstance = GetModuleHandle(nullptr);
                wc.style         = CS_OWNDC;
                wc.lpfnWndProc   = (WNDPROC)WindowProc;
                wc.cbClsExtra    = 0;
                wc.cbWndExtra    = 0;
                wc.hInstance     = hInstance;
                wc.hIcon         = LoadIcon(nullptr, IDI_WINLOGO);
                wc.hCursor       = LoadCursor(nullptr, IDC_ARROW);
                wc.hbrBackground = nullptr;
                wc.lpszMenuName  = nullptr;
                wc.lpszClassName = "OpenGL";

                if (!RegisterClass(&wc)) {
                        MessageBox(nullptr, "RegisterClass() failed:  "
                                   "Cannot register window class.",
                                   "Error", MB_OK);
                        return nullptr;
                }
    }

    hWnd = CreateWindow("OpenGL", title, WS_OVERLAPPEDWINDOW |
                        WS_CLIPSIBLINGS | WS_CLIPCHILDREN,
                        0,0,Width,Height, nullptr, nullptr, hInstance, nullptr);

    if (hWnd == nullptr) {
                MessageBox(nullptr,
                           "CreateWindow() failed:  Cannot create a window.",
                           "Error", MB_OK);
                return nullptr;
    }

    hDC = GetDC(hWnd);

    /* there is no guarantee that the contents of the stack that become
       the pfd are zeroed, therefore _make sure_ to clear these bits. */
    memset(&pfd, 0, sizeof(pfd));
    pfd.nSize        = sizeof(pfd);
    pfd.nVersion     = 1;
    pfd.dwFlags      = PFD_DRAW_TO_WINDOW
                     | PFD_SUPPORT_OPENGL
                     | PFD_DOUBLEBUFFER;
    pfd.iPixelType   = PFD_TYPE_RGBA;
    pfd.cDepthBits   = 32;
    pfd.cColorBits   = 32;

    pf = ChoosePixelFormat(hDC, &pfd);
    if (pf == 0) {
                MessageBox(nullptr, "ChoosePixelFormat() failed:  "
                           "Cannot find a suitable pixel format.",
                           "Error", MB_OK);
                return 0;
    }

    if (SetPixelFormat(hDC, pf, &pfd) == FALSE) {
                MessageBox(nullptr, "SetPixelFormat() failed:  "
                           "Cannot set format specified.", "Error", MB_OK);
                return 0;
    }

    DescribePixelFormat(hDC, pf, sizeof(PIXELFORMATDESCRIPTOR), &pfd);

    if (pfd.dwFlags & PFD_NEED_PALETTE ||
                pfd.iPixelType == PFD_TYPE_COLORINDEX) {

                n = 1 << pfd.cColorBits;
                if (n > 256) n = 256;

                lpPal = (LOGPALETTE*)malloc(sizeof(LOGPALETTE) +
                                                sizeof(PALETTEENTRY) * n);
                memset(lpPal, 0, sizeof(LOGPALETTE) + sizeof(PALETTEENTRY) * n);
                lpPal->palVersion = 0x300;
                lpPal->palNumEntries = n;

                GetSystemPaletteEntries(hDC, 0, n, &lpPal->palPalEntry[0]);

                /* if the pixel type is RGBA, then we want to make an RGB ramp,
                   otherwise (color index) set individual colors. */
                if (pfd.iPixelType == PFD_TYPE_RGBA) {
                        int redMask = (1 << pfd.cRedBits) - 1;
                        int greenMask = (1 << pfd.cGreenBits) - 1;
                        int blueMask = (1 << pfd.cBlueBits) - 1;
                        int i;

                        /* fill in the entries with an RGB color ramp. */
                        for (i = 0; i < n; ++i) {
                        lpPal->palPalEntry[i].peRed =
                                (((i >> pfd.cRedShift)   & redMask)   * 255)
                               /redMask;
                        lpPal->palPalEntry[i].peGreen =
                                (((i >> pfd.cGreenShift) & greenMask) * 255)
                               /greenMask;
                        lpPal->palPalEntry[i].peBlue =
                                (((i >> pfd.cBlueShift)  & blueMask)  * 255)
                               /blueMask;
                        lpPal->palPalEntry[i].peFlags = 0;
                        }
                } else {
                        lpPal->palPalEntry[0].peRed = 0;
                        lpPal->palPalEntry[0].peGreen = 0;
                        lpPal->palPalEntry[0].peBlue = 0;
                        lpPal->palPalEntry[0].peFlags = PC_NOCOLLAPSE;
                        lpPal->palPalEntry[1].peRed = 255;
                        lpPal->palPalEntry[1].peGreen = 0;
                        lpPal->palPalEntry[1].peBlue = 0;
                        lpPal->palPalEntry[1].peFlags = PC_NOCOLLAPSE;
                        lpPal->palPalEntry[2].peRed = 0;
                        lpPal->palPalEntry[2].peGreen = 255;
                        lpPal->palPalEntry[2].peBlue = 0;
                        lpPal->palPalEntry[2].peFlags = PC_NOCOLLAPSE;
                        lpPal->palPalEntry[3].peRed = 0;
                        lpPal->palPalEntry[3].peGreen = 0;
                        lpPal->palPalEntry[3].peBlue = 255;
                        lpPal->palPalEntry[3].peFlags = PC_NOCOLLAPSE;
                }

                hPalette = CreatePalette(lpPal);
                if (hPalette) {
                        SelectPalette(hDC, hPalette, FALSE);
                        RealizePalette(hDC);
                }

                free(lpPal);
    }

    ReleaseDC(hDC, hWnd);
    return hWnd;
}

int APIENTRY WinMain(HINSTANCE hCurrentInst, HINSTANCE hPreviousInst,
        LPSTR lpszCmdLine, int nCmdShow)
{
    HGLRC hRC;                          /* opengl context */
    HWND  hWnd;                         /* window */
    MSG   msg;                          /* message */

        // InitModel() initialises some data structures and
        // does the progressive mesh polygon reduction algorithm
        // on the model.
        CalcFPSDeltaT(); // to time the algorithm
        InitModel();
        CalcFPSDeltaT();

        hWnd = CreateOpenGLWindow("bunnylod by Stan Melax");
    if (hWnd == nullptr) exit(1);

    hDC = GetDC(hWnd);
    hRC = wglCreateContext(hDC);
    wglMakeCurrent(hDC, hRC);
    ShowWindow(hWnd, nCmdShow);
        glEnable(GL_DEPTH_TEST);

        PostString("Demo by Stan Melax (c)1998",5,-5,20);
        PostString("Model by Viewpoint Datalabs (c)1996",5,-4,20);
        char buf[128];
        PostString("Mesh Reduction Algorithm (non-optimised)",1,0,5);
        sprintf(buf,"was executed in %5.3f seconds",DeltaT);
        PostString(buf,2,1,6);

    while (1) {
                while(PeekMessage(&msg, hWnd, 0, 0, PM_NOREMOVE)) {
                    if(GetMessage(&msg, hWnd, 0, 0)) {
                                TranslateMessage(&msg);
                                DispatchMessage(&msg);
                    } else {
                                // This 'goto' was in the sample code
                                goto quit;
                    }
                }
                CalcFPSDeltaT();
                Display();
    }

  quit:
    wglMakeCurrent(nullptr, nullptr);
    ReleaseDC(hDC, hWnd);
    wglDeleteContext(hRC);
    DestroyWindow(hWnd);
    if (hPalette) DeleteObject(hPalette);
    return msg.wParam;
}
