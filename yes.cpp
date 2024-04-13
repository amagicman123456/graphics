#include <windows.h>
#include <functional>
#include <cmath>
#include <iostream>
#include <process.h>
#include <vector>
#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
int f = 0;
void fps(void*){
    while(1){
        Sleep(1000);
        std::cout << f << '\n';
        f = 0;
    }
}
int width = 1200, height = 600, tick_count, *framebuf;
double horizontal_fov = 120, vertical_fov = atan(tan(height / 2.0) * height / width);
double square(double x){
    return x * x;
}
double greater(double a, double b){
    return a > b ? a : b;
}
double abs_val(double a){
    return a < 0 ? -a : a;
}
double dot_product(double ax, double ay, double az, double bx, double by, double bz){
    return ax * bx + ay * by + az * bz;
}
struct point{
    double x, y, z;
    point(double e, double f, double g) : x(e), y(f), z(g){}
    point operator-(point f){
        return point(x - f.x, y - f.y, z - f.z);
    }
    point& operator-=(point g){
        x -= g.x;
        y -= g.y;
        z -= g.z;
        return *this;
    }
};
typedef point vector;
struct sphere{
    double radius;
    point center;
    sphere(double r, point c) : radius(r), center(c){}
    bool hit(vector u, double& distance){
        double magnitude = sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
        ux /= magnitude, uy /= magnitude, uz /= magnitude;
        double dot = dot_product(u.x, u.y, u.z, center.x, center.y, center.z);
        double determinant = square(dot) - center.x * center.x - center.y * center.y - center.z * center.z + radius * radius;
        distance = -dot - sqrt(determinant);
        return determinant >= 0;
    }
};
struct polygon{
    std::vector<point> points{};
    polygon(auto... l) try{
        if(sizeof...(l) < 3) throw;
        points = {l..};
        a = points[1] - points[0], b = points[2] - points[0], c = cross_product(a, b);
        k = c.x * points[0].x + c.y * points[1].y + c.z * points[2].z;
    }catch(...){std::cout << "error: nunber of points to polygon's constructor must be greater than two\n";}
    bool hit(vector u, double& distance){
        double e = u.x * c.x + u.y * c.y + u.z * c.z;
        double t = k / e;
        if(!e || (uz * t < 0)) return false;
        point intersection(u.x * t, u.y * t, u.z * t);
        distance = sqrt(intersection.x * intersection.x + intersection.y * intersection.y + intersection.z * intersection.z);
        double greatest = greater(greater(c1.z, c2.z), c3.z);
        auto stretch = [&greatest](point& p){
            double factor = greatest / p.z;
            p.x *= factor;
            p.y *= factor;
            p.z = greatest;
        };
        for(point& i : points) stretch(i);
        stretch(intersection);

        #define tri_specialization
        #define first_algorithm

        #ifdef tri_specialization
        if(n == 3){
            #ifdef first_algorithm
                if(points[2].y == points[0].y){
                    point temp = points[2];
                    points[2] = points[1];
                    points[1] = temp;
                }
                double s1 = points[2].y - points[0].y,
                       s2 = points[2].x - points[0].x,
                       s3 = points[1].y - points[0].y,
                       s4 = p.y - points[0].y;
                double w1 = (points[0].x * s1 + s4 * s2 - p.x * s1) / (s3 * s2 - (polygon[1].x - polygon[0].x) * s1),
                       w2 = (s4 - w1 * s3) / s1;
                return w1 >= 0 && w2 >= 0 && (w1 + w2) <= 1;
            #else
                auto sign = [&](point p1, point p2, point p3) -> double{
                    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
                };
                double d1 = sign(p, points[0], points[1]),
                       d2 = sign(p, points[1], points[2]),
                       d3 = sign(p, points[2], points[0]);
                bool neg = d1 < 0 || d2 < 0 || d3 < 0,
                     pos = d1 > 0 || d2 > 0 || d3 > 0;
                return !(neg && pos);
            #endif
        }
        #endif
        bool inside = false;
        for(int i = 0, j = n - 1; i < n; j = i++)
            if(((points[i].y > p.y) != (points[j].y > p.y)) &&
                (p.x < (points[j].x - points[i].x) * (p.y - points[i].y) / (points[j].y - points[i].y) + points[i].x))
                    inside = !inside;
                    return inside;
    }
private:
    vector a, b, c;
    double k;
};
vector cross_product(vector a, vector b){
    return vector(a.y * b.z - b.y * a.z, a.z * b.x - b.z * a.x, a.x * b.y - b.x * a.y);
}
std::function<void()> render =
[&](){
    ++f;
    double z = width / (2 * tan(horizontal_fov / 2)); //z_sq = z * z;
    for(int j = height - 1; j >= 0; --j){
        for(int i = 0; i < width; ++i){
            //double ux = -width / 2.0 + i, uy = height / 2.0 - j, uz = z;
            vector v(-width / 2.0 + i, height / 2.0 - j, z);
            /*
            double magnitude = sqrt(ux * ux + uy * uy + z_sq);
            ux /= magnitude, uy /= magnitude, uz /= magnitude;
            double circle_x = 100, circle_y = 41.5, circle_z = 2000, radius = sqrt(100000);
            double dot = dot_product(ux, uy, uz, circle_x, circle_y, circle_z);
            double determinant = square(dot) - circle_x * circle_x - circle_y * circle_y - circle_z * circle_z + radius * radius;
            bool hit_sphere = determinant >= 0;
            double distance_to_sphere = -dot - sqrt(determinant);
            */
            sphere s(sqrt(100000), point(100, 41.5, 2000));
            double distance_to_sphere = 0;
            bool hit_sphere = s.hit(v, std::ref(distance_to_sphere));
            
            #define TRI
            /*
            //point c1(-250, 100, 2000), c2(250, 100, 2000), c3(250, -100, 2000);
            //point c1(-100, 150, 2000), c2(250, 100, 2000), c3(250, -150, 2000);
            point c1(-250, -300, 2400), c2(250, 200, 1500), c3(350, -100, 1000);
            #ifdef QUAD
                point c4(-250, -100, 2000);
            #endif
            vector a = c2 - c1, b = c3 - c1, c = cross_product(a, b);
            //ax + by + cz = k
            double k = c.x * c1.x + c.y * c2.y + c.z * c3.z;
            double e = ux * c.x + uy * c.y + uz * c.z;
            double t = k / e;
            if(!e || (uz * t < 0)) framebuf[j * width + i] = RGB(255, 0, 0);
            else{
                point intersection(ux * t, uy * t, uz * t);
                double distance_to_plane = sqrt(intersection.x * intersection.x + intersection.y * intersection.y + intersection.z * intersection.z);
                double greatest = greater(greater(c1.z, c2.z), c3.z);
                auto stretch = [&greatest](point& p){
                    double factor = greatest / p.z;
                    p.x *= factor;
                    p.y *= factor;
                    p.z = greatest;
                };
                stretch(c1);
                stretch(c2);
                stretch(c3);
                #ifdef QUAD
                    stretch(c4);
                #endif
                stretch(intersection);
                auto is_inside = [](point p, auto ...points) -> bool{
                    constexpr int n = sizeof...(points);
                    point polygon[n] = {points...};

                    #define tri_specialization
                    #define first_algorithm

                    #ifdef tri_specialization
                        if(n == 3){
                            #ifdef first_algorithm
                                if(polygon[2].y == polygon[0].y){
                                    point temp = polygon[2];
                                    polygon[2] = polygon[1];
                                    polygon[1] = temp;
                                }
                                double s1 = polygon[2].y - polygon[0].y,
                                       s2 = polygon[2].x - polygon[0].x,
                                       s3 = polygon[1].y - polygon[0].y,
                                       s4 = p.y - polygon[0].y;
                                double w1 = (polygon[0].x * s1 + s4 * s2 - p.x * s1) / (s3 * s2 - (polygon[1].x - polygon[0].x) * s1),
                                       w2 = (s4 - w1 * s3) / s1;
                                return w1 >= 0 && w2 >= 0 && (w1 + w2) <= 1;
                            #else
                                auto sign = [&](point p1, point p2, point p3) -> double{
                                    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
                                };
                                double d1 = sign(p, polygon[0], polygon[1]),
                                       d2 = sign(p, polygon[1], polygon[2]),
                                       d3 = sign(p, polygon[2], polygon[0]);
                                bool neg = d1 < 0 || d2 < 0 || d3 < 0,
                                     pos = d1 > 0 || d2 > 0 || d3 > 0;
                                return !(neg && pos);
                            #endif
                        }
                    #endif

                    bool inside = false;
                    for(int i = 0, j = n - 1; i < n; j = i++)
                        if(((polygon[i].y > p.y) != (polygon[j].y > p.y)) &&
                            (p.x < (polygon[j].x - polygon[i].x) * (p.y - polygon[i].y) / (polygon[j].y - polygon[i].y) + polygon[i].x))
                            inside = !inside;
                    return inside;
                };
                /*
                #ifdef TRI
                    bool hit_triangle = is_inside(intersection, c1, c2, c3);
                    if(hit_sphere || hit_triangle) framebuf[j * width + i] = RGB(0, 0, 255);
                    else framebuf[j * width + i] = RGB(255, 0, 0);
                #endif
                #ifdef QUAD
                    bool hit_quadrilateral = is_inside(intersection, c1, c2, c3, c4);
                    if(hit_sphere || hit_quadrilateral) framebuf[j * width + i] = RGB(0, 0, 255);
                    else framebuf[j * width + i] = RGB(255, 0, 0);
                #endif
                *//*
                bool hit_polygon = is_inside(intersection, c1, c2, c3
                #ifdef QUAD
                    ,c4
                #endif
                );
                */
                polygon p(point(-250, -300, 2400), point(250, 200, 1500), point(350, -100, 1000)
                #ifdef QUAD
                    , point(-250, -100, 2000)
                );
                double distance_to_polygon;
                bool hit_polygon = p.hit(v, std::ref(distance_to_polygon));
                // for multiple objects create a list of objects the ray hit
                // with either inheritance or std::any or smth
                // and find which of the distances is least
                if(hit_sphere || hit_polygon) framebuf[j * width + i] = RGB(0, 0, 255);
                else framebuf[j * width + i] = RGB(255, 0, 0);
            }
        }
    }
}
, resize =
[](){

};
LRESULT CALLBACK WindowProcessMessages(HWND hwnd, UINT msg, WPARAM w, LPARAM l){
    static HDC pdc;
    static HBITMAP old;
    static HBITMAP ourbitmap;
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
            resize();
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
            tick_count = GetTickCount();
            PAINTSTRUCT ps;
            HDC h = BeginPaint(hwnd, &ps);
            render();
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
int main(){
    _beginthread(fps, 0, 0);
    char name[] = "omg";
    WNDCLASS wc{};
    wc.lpszClassName = name;
    wc.hCursor = LoadCursor(0, IDC_ARROW);
    wc.lpfnWndProc = WindowProcessMessages;
    RegisterClass(&wc);
    int a = width + 16, b = height + 39;
    CreateWindow(name, name, WS_OVERLAPPEDWINDOW | WS_CAPTION | WS_SYSMENU | WS_VISIBLE, CW_USEDEFAULT, CW_USEDEFAULT, a, b, 0, 0, 0, 0);
    MSG msg{};
    while(GetMessage(&msg, 0, 0, 0)){
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
}