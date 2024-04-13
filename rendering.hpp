#ifndef rendering_hpp
#define rendering_hpp
#ifndef _WINDOWS_H
#include <windows.h>
#endif
#include <functional>
#include <cstdlib>
std::function<void()> _render, _resize = [](){};
int _width, _height, _window_width, _window_height, _offset, _tick_count, *_framebuf;
uint32_t rgb(uint8_t red, uint8_t green, uint8_t blue){
    uint8_t a[4] = {blue, green, red, 0};
    uint32_t b;
    memcpy(&b, a, 4);
    return b;
    //return *(uint32_t*)a;
}
LRESULT CALLBACK WindowProcessMessages(HWND hwnd, UINT msg, WPARAM w, LPARAM l){
    static HDC pdc;
    static HBITMAP old;
    static HBITMAP ourbitmap;
    switch(msg){
        case WM_GETMINMAXINFO:{
            LPMINMAXINFO e = (LPMINMAXINFO)l;
            e->ptMinTrackSize.x = 6 + _offset + 1;
            e->ptMinTrackSize.y = 29 + _offset + 1;
            break;
        }
        case WM_CREATE:{
            SetTimer(hwnd, 1, 1, 0);
            HDC hdc;
            BITMAPINFO bitmapinfo{};
            hdc = CreateCompatibleDC(0);
            bitmapinfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
            bitmapinfo.bmiHeader.biWidth = _width;
            bitmapinfo.bmiHeader.biHeight = -_height; // top down is negative
            bitmapinfo.bmiHeader.biPlanes = 1;
            bitmapinfo.bmiHeader.biBitCount = 32;
            bitmapinfo.bmiHeader.biCompression = BI_RGB;
            bitmapinfo.bmiHeader.biClrUsed = 256;
            bitmapinfo.bmiHeader.biClrImportant = 256;
            ourbitmap = CreateDIBSection(hdc, &bitmapinfo, DIB_RGB_COLORS, (void**)&_framebuf, 0, 0);
            pdc = CreateCompatibleDC(0);
            old = (HBITMAP)SelectObject(pdc, ourbitmap);
            DeleteDC(hdc);
            break;
        }
        case WM_SIZE:{
            _width = LOWORD(l);
            _height = HIWORD(l);
            _resize();
            SelectObject(pdc, old);
            DeleteDC(pdc);
            DeleteObject(ourbitmap);
            HDC hdc;
            BITMAPINFO bitmapinfo{};
            hdc = CreateCompatibleDC(0);
            bitmapinfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
            bitmapinfo.bmiHeader.biWidth = _width;
            bitmapinfo.bmiHeader.biHeight = -_height; // top down is negative
            bitmapinfo.bmiHeader.biPlanes = 1;
            bitmapinfo.bmiHeader.biBitCount = 32;
            bitmapinfo.bmiHeader.biCompression = BI_RGB;
            bitmapinfo.bmiHeader.biClrUsed = 256;
            bitmapinfo.bmiHeader.biClrImportant = 256;
            ourbitmap = CreateDIBSection(hdc, &bitmapinfo, DIB_RGB_COLORS, (void**)&_framebuf, 0, 0);
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
            _tick_count = GetTickCount();
            PAINTSTRUCT ps;
            HDC h = BeginPaint(hwnd, &ps);
            _render();
            BitBlt(h, 0, 0, _width, _height, pdc, 0, 0, SRCCOPY);
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
class screen{
public:
    RECT dimensions, work_area;
    screen(){
        GetWindowRect(GetDesktopWindow(), &dimensions);
        SystemParametersInfo(SPI_GETWORKAREA, 0, &work_area, 0);
    }
    int width() const{
        return dimensions.right;
    }
    int width_area() const{
        return work_area.right;
    }
    int height() const{
        return dimensions.bottom;
    }
    int height_area() const{
        return work_area.bottom;
    }
};
class window{
public:
    int &width = _width, &height = _height, &window_width = _window_width, &window_height = _window_height, &offset = _offset, *&pixels = _framebuf;
    enum{
        invisible = 1,
        maximize,
        no_resize = 4
    };
    window(const char* name, int w, int h, int x = CW_USEDEFAULT, int y = CW_USEDEFAULT, int flags = 0){
        _width = w, _height = h;
        WNDCLASS wc{};
        wc.lpszClassName = name;
        wc.hCursor = LoadCursor(0, IDC_ARROW);
        wc.lpfnWndProc = WindowProcessMessages;
        RegisterClass(&wc);
        bool e = !(flags & no_resize);
        offset = e * 10;
        int resize_flags = WS_THICKFRAME | WS_MINIMIZEBOX | WS_MAXIMIZEBOX;
        //wm_size assigns _width and _height to maximized width and height, but safety first
        if(flags & maximize){
            screen s;
            window_width = s.width_area(), window_height = s.height_area();
            width = window_width - 6 - offset, height = window_height - 29 - offset;
        }else window_width = width + 6 + offset, window_height = height + 29 + offset;
        //window_width = width + 6 + offset, window_height = height + 29 + offset;
        CreateWindow(name, name, WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | (e * resize_flags) | (!(flags & invisible) * WS_VISIBLE) | (!!(flags & maximize) * WS_MAXIMIZE), x - (((uint32_t)x != CW_USEDEFAULT) << 3), y, window_width, window_height, 0, 0, 0, 0);
    }
    /*
    void draw_line(int x1, int y1, int x2, int y2, uint32_t color = 0) {
        int dx = abs(x2 - x1),
            dy = abs(y2 - y1),
            sx = (x1 < x2) - (x1 >= x2),
            sy = (y1 < y2) - (y1 >= y2),
            err = dx - dy, e2;
        while(x1 != x2 || y1 != y2){
            if(x1 < width && x1 > -1 && y1 < height && y1 > -1) pixels[x1 + y1 * width] = color;
            e2 = err << 1;
            if(e2 > -dy) err -= dy, x1 += sx;
            if(e2 < dx) err += dx, y1 += sy;
        }
    }
    */
    int tick_count(){
        return _tick_count;
    }
    void render(std::function<void()> func){
        _render = func;
    }
    void on_resize(std::function<void()> func){
        _resize = func;
    }
    static inline void loop(){
        MSG msg{};
        while(GetMessage(&msg, 0, 0, 0)){
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }
};
#endif