#include <windows.h>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
HDC pdc;
HBITMAP old;
HBITMAP ourbitmap;
int width, height, *framebuf;
//RECT dimensions;
struct point3d{
    double x, y, z;
    point3d(double a, double b, double c) : x(a), y(b), z(c){}
    double x2d(double fov){
        return x / (z * tan(fov * 3.1415926592653589793 / 360)) * width + width / 2;
    }
    double y2d(double fov){
        return height - (y / (z * tan(fov * 3.1415926592653589793 / 360)) * height + height / 2);
    }
    void e(double fov){
        cout << this->x2d(fov) << ' ' << this->y2d(fov) << '\n';
    }
};
void draw_line(int x1, int y1, int x2, int y2) {
    int dx = abs(x2 - x1),
        dy = abs(y2 - y1),
        sx = (x1 < x2) - (x1 >= x2),// ? 1 : -1,
        sy = (y1 < y2) - (y1 >= y2),// ? 1 : -1,
        err = dx - dy, e2;
    while(x1 != x2 || y1 != y2){
        if(x1 < width && x1 > -1 && y1 < height && y1 > -1) framebuf[x1 + y1 * width] = 0;
        e2 = err << 1;
        if(e2 > -dy) err -= dy, x1 += sx;
        if(e2 < dx) err += dx, y1 += sy;
    }
}
void draw_line(point3d a, point3d b, double fov){
    draw_line(a.x2d(fov), a.y2d(fov), b.x2d(fov), b.y2d(fov));
}
LRESULT CALLBACK WindowProcessMessages(HWND hwnd, UINT msg, WPARAM w, LPARAM l){
    switch(msg){
        case WM_CREATE:{
            SetTimer(hwnd, 1, 1, 0);
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
        case WM_SIZE:{
            width = LOWORD(l);
            height = HIWORD(l);
            SelectObject(pdc, old);
            DeleteDC(pdc);
            DeleteObject(ourbitmap);
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
        case WM_TIMER:
            InvalidateRgn(hwnd, 0, 0);
            UpdateWindow(hwnd);
            break;
        case WM_PAINT:{
            // GetTickCount() returns the tick number, which can be useful for animations
            PAINTSTRUCT ps;
            HDC h = BeginPaint(hwnd, &ps);
            for(int i = 0; i < width*height; ++i)
                framebuf[i] = RGB(255, 0, 0);
            //draw_line(0, 0, width, height);
            //constexpr double fov = 90, a = tan(fov * 3.1415926592653589793 / 360);
            //draw_line(
            //    -5 / (100 * a) * width + width / 2,
            //    height - (-10 / (100 * a) * height + height / 2),
            //    5 / (50 * a) * width + width / 2,
            //    height - (10 / (50 * a) * height + height / 2)
            //);
            {
                point3d a = {50, 50, 1000.5},
                        b = {50, -50, 900.5},
                        c = {50, 50, 900.5},
                        d = {50, -50, 1000.5},
                        e = {-50, 50, 1000.5},
                        f = {-50, -50, 900.5},
                        g = {-50, 50, 900.5},
                        h = {-50, -50, 1000.5};
                double fov = 45;
                draw_line(a, c, fov);
                draw_line(a, d, fov);
                draw_line(a, e, fov);
                draw_line(f, g, fov);
                draw_line(f, h, fov);
                draw_line(f, b, fov);
                draw_line(b, c, fov);
                draw_line(d, h, fov);
                draw_line(d, b, fov);
                draw_line(h, e, fov);
                draw_line(e, g, fov);
                draw_line(g, c, fov);
            }
            /*
            for(point3d (&c)[4] : a){
                //for(point3d& i : c) i.z += 10;
                for(int i = 0; i < 3; ++i){
                    cout << c[i].x2d(90) << ' ' << c[i].y2d(90) << ' ' << c[i + 1].x2d(90) << ' ' << c[i + 1].y2d(90) << '\n';
                    draw_line(
                        c[i].x2d(90),
                        c[i].y2d(90),
                        c[i + 1].x2d(90),
                        c[i + 1].y2d(90)
                    );
                }
            }
            */
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
    //GetWindowRect(GetDesktopWindow(), &dimensions);
    //width = dimensions.right;
    //height = dimensions.bottom;
    width = 600;
    height = 600;

    // edit the buffer in WM_PAINT
    // framebuf = new int[width * height];
    char name[] = "hehehehaw";

    WNDCLASS wc{};
    wc.hInstance = currentInstance;
    wc.lpszClassName = name;
    wc.hCursor = LoadCursor(0, IDC_ARROW);
    wc.lpfnWndProc = WindowProcessMessages;
    RegisterClass(&wc);
    CreateWindow(name, name, WS_OVERLAPPEDWINDOW | WS_VISIBLE/* | WS_MAXIMIZE*/, CW_USEDEFAULT, CW_USEDEFAULT, width+/*2*/16, height+/*16+2*/39, 0, 0, 0, 0);
    MSG msg{};

    while(GetMessage(&msg, 0, 0, 0)){
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return 0;
}