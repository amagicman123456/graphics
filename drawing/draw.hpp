#ifndef draw_hpp
    #define draw_hpp
    #ifndef _WINDOWS_H
        #include <windows.h>
    #endif
    class window{
    private:
        HWND hwindow;
    public:
        int height, width;
        window(){}
        window(const char* name, int w, int h, HINSTANCE currentInstance, WNDPROC WindowProcessMessages, bool resizable) : width(w), height(h){
            WNDCLASS wc{};
            wc.hInstance = currentInstance;
            wc.lpszClassName = name;
            wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
            wc.hbrBackground = CreateSolidBrush(RGB(0, 0, 0));
            wc.lpfnWndProc = WindowProcessMessages;
            RegisterClass(&wc);
            hwindow = CreateWindow(name, name, resizable ? WS_OVERLAPPEDWINDOW | WS_VISIBLE : WS_OVERLAPPEDWINDOW & ~WS_MAXIMIZEBOX & ~WS_THICKFRAME | WS_VISIBLE, CW_USEDEFAULT, CW_USEDEFAULT, w, h, 0, 0, 0, 0);
        }
        operator HWND(){
            return hwindow;
        }
    };
    class scoped_hdc{
    private:
        HWND hwnd;
        bool mem;
        HDC hdc;
    public:
        scoped_hdc(HDC h) : hdc(h), mem(1){}
        scoped_hdc(HWND h) : hwnd(h), hdc(GetDC(h)), mem(0){}
        scoped_hdc(HWND hd, HDC hc) : hwnd(hd), hdc(hc), mem(0){}
        ~scoped_hdc(){
            mem ? DeleteDC(hdc) : ReleaseDC(hwnd, hdc);
        }
        operator HDC(){
            return hdc;
        }
    };
#endif