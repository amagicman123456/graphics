#include <windows.h>
HDC pdc;
HBITMAP old;
HBITMAP ourbitmap;
int width, height, *framebuf;
RECT dimensions;
LRESULT CALLBACK WindowProcessMessages(HWND hwnd, UINT msg, WPARAM w, LPARAM l){
    switch (msg){
        case WM_CREATE:{
            //SetTimer(hwnd, 1, 1, 0);
            HDC hdc;
            BITMAPINFO bitmapinfo{};
            hdc = CreateCompatibleDC(0);
            bitmapinfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
            bitmapinfo.bmiHeader.biWidth = width;
            bitmapinfo.bmiHeader.biHeight = -height; // top down is negative
            bitmapinfo.bmiHeader.biPlanes = 1;
            bitmapinfo.bmiHeader.biBitCount = 32;
            bitmapinfo.bmiHeader.biCompression = BI_RGB;
            bitmapinfo.bmiHeader.biClrUsed = 256;
            bitmapinfo.bmiHeader.biClrImportant = 256;
            ourbitmap = CreateDIBSection(hdc, &bitmapinfo, DIB_RGB_COLORS, (void**)&framebuf, 0, 0);
            pdc = CreateCompatibleDC(0);
            old = (HBITMAP)SelectObject(pdc, ourbitmap);
            DeleteDC(hdc);
            break;
        }
        case WM_SIZE:
            GetWindowRect(GetDesktopWindow(), &dimensions);
            width = dimensions.right;
            height = dimensions.bottom;
            break;
        //case WM_TIMER:
        //    InvalidateRgn(hwnd, 0, 0);
        //    UpdateWindow(hwnd);
        //    break;
        case WM_PAINT:{
            // GetTickCount() returns the tick number, which can be useful for animations
            PAINTSTRUCT ps;
            HDC h = BeginPaint(hwnd, &ps);
            for(int i = 0; i < width*height; ++i) framebuf[i] = RGB(255, 0, 0);
            BitBlt(h, 0, 0, width, height, pdc, 0, 0, SRCCOPY);
            EndPaint(hwnd, &ps);
            break;
        }
        case WM_DESTROY:
            SelectObject(pdc, old);
            DeleteDC(pdc);
            DeleteObject(ourbitmap);
            KillTimer(hwnd, 1);
            exit(0);
        default:
            return DefWindowProc(hwnd, msg, w, l);
    }
    return 0;
}
int WINAPI WinMain(HINSTANCE currentInstance, HINSTANCE previousInstance, PSTR cmdLine, int cmdCount){
    GetWindowRect(GetDesktopWindow(), &dimensions);
    width = dimensions.right;
    height = dimensions.bottom;
    // edit the buffer in WM_PAINT
    framebuf = new int[width * height];
    char name[] = "hehehehaw";

    WNDCLASS wc{};
    wc.hInstance = currentInstance;
    wc.lpszClassName = name;
    wc.hCursor = LoadCursor(0, IDC_ARROW);
    wc.lpfnWndProc = WindowProcessMessages;
    RegisterClass(&wc);
    CreateWindow(name, name, WS_OVERLAPPEDWINDOW | WS_VISIBLE | WS_MAXIMIZE, CW_USEDEFAULT, CW_USEDEFAULT, width+2, height+16+2, 0, 0, 0, 0);
    MSG msg{};

    while(GetMessage(&msg, 0, 0, 0)){
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return 0;
}