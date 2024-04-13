#include <windows.h>
#include <math.h>

char SINTAB[256];
HDC pDC;
HBITMAP old;
HBITMAP ourbitmap;
int *framebuf = (int*)malloc(320*200*4);
void render_effect(int tick, int *framebuf){
    //memset(framebuf, RGB(0, 0, 255), 320*200*4);
    //int k = 0;
    //for(int i = 0; i < 200; ++i)
    //    for(int j = 0; j < 320; ++j, ++k)
    //        framebuf[k] = RGB(0, 0, 255);
    for(int i = 0; i < 200*320; ++i) framebuf[i] = RGB(0, 0, 255);
/*
    int i, j, k;
    tick >>= 2;
    for(k = 0, i = 0; i < 200; i++)
        for(j = 0; j < 320; j++, k++)
            *(framebuf+k)=RGB(SINTAB[(i + tick) & 0xff],
                SINTAB[(j - tick) & 0xff],
                SINTAB[(SINTAB[tick & 0xff] + (k >> 6)) & 0xff]);
*/
}
LRESULT CALLBACK winproc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam){
    HDC hDC;
    PAINTSTRUCT PtStr;
    switch (uMsg) {
    case WM_DESTROY:
        SelectObject(pDC,old);
        DeleteDC(pDC);
        DeleteObject(ourbitmap);
        PostQuitMessage(0);
        KillTimer (hWnd, 1);
        break;
    case WM_CREATE:{
        SetTimer (hWnd, 1, 1, NULL);
        HDC hDC;
        BITMAPINFO bitmapinfo{};
        hDC=CreateCompatibleDC(0);
        bitmapinfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
        bitmapinfo.bmiHeader.biWidth = 320;
        bitmapinfo.bmiHeader.biHeight = -200; /* top-down */
        bitmapinfo.bmiHeader.biPlanes = 1;
        bitmapinfo.bmiHeader.biBitCount = 32;
        bitmapinfo.bmiHeader.biCompression = BI_RGB;
        bitmapinfo.bmiHeader.biClrUsed = 256;
        bitmapinfo.bmiHeader.biClrImportant = 256;
        ourbitmap = CreateDIBSection(hDC, &bitmapinfo, DIB_RGB_COLORS, (void**)&framebuf, 0, 0);
        pDC=CreateCompatibleDC(NULL);
        old=(HBITMAP)SelectObject(pDC, ourbitmap);
        DeleteDC(hDC);
        break;
    }
    case WM_TIMER:
        InvalidateRgn(hWnd,0,0);
        UpdateWindow (hWnd);
        break;
    case WM_PAINT:
        hDC=BeginPaint(hWnd,&PtStr);
        render_effect(GetTickCount(),framebuf);
        BitBlt(hDC, 0, 0, 320, 200, pDC, 0, 0, SRCCOPY);
        EndPaint(hWnd,&PtStr);
        break;
    default:
        return DefWindowProc (hWnd, uMsg, wParam, lParam);
        break;
    }
    return 0;
}
int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR lpszCmdLine, int nCmdShow){
    char name[] = "hehehehaw";
    WNDCLASSEX winclass{};
    HWND hWnd;
    MSG msg;
    int i;
    for (i=0;i<256;i++)
        SINTAB[i]=sin(((i+1)*3.14159265359)/128)*127+128;

    winclass.cbSize = sizeof(WNDCLASSEX);
    winclass.style = CS_DBLCLKS;
    winclass.lpfnWndProc = &winproc;
    winclass.hInstance=hInst;
    winclass.hIcon = LoadIcon(NULL,IDI_WINLOGO);
    winclass.hCursor = LoadCursor(NULL,IDC_NO);
    winclass.lpszClassName = name;
    winclass.hIconSm=NULL;

    if (!RegisterClassEx(&winclass))
        return 0;
    hWnd = CreateWindow(name, name, WS_SYSMENU | WS_CAPTION | WS_BORDER | WS_OVERLAPPED | WS_VISIBLE | WS_MINIMIZEBOX, CW_USEDEFAULT, 0, 320+2, 200+16+2, 0, 0, hInst, 0);
    ShowWindow(hWnd,nCmdShow);
    UpdateWindow(hWnd);
    while(GetMessage(&msg, 0, 0, 0)){
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return (msg.wParam);
}