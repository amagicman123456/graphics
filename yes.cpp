#include <windows.h>
#include <process.h>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
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
    point() : x(0), y(0), z(0){}
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
typedef uint32_t color;
struct object{
    color clr;
    virtual std::pair<bool, double> hit() = 0;
};
typedef point vector;
struct sphere : object{
    double radius;
    point center;
    color clr;
    sphere(color cl, double r, point c) : radius(r), center(c), clr(cl){}
    std::pair<bool, double> hit(vector u) override{
        double magnitude = sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
        u.x /= magnitude, u.y /= magnitude, u.z /= magnitude;
        double dot = dot_product(u.x, u.y, u.z, center.x, center.y, center.z);
        double determinant = square(dot) - center.x * center.x - center.y * center.y - center.z * center.z + radius * radius;
        return std::pair<bool, double>(determinant >= 0, -dot - sqrt(determinant));
    }
};
vector cross_product(vector a, vector b){
    return vector(a.y * b.z - b.y * a.z, a.z * b.x - b.z * a.x, a.x * b.y - b.x * a.y);
}
struct polygon : object{
    std::vector<point> points{};
    color clr;
    polygon(color cl, auto... l) try : clr(cl){
        if(sizeof...(l) < 3) throw;
        points = {(point(l))...};
        a = points[1] - points[0], b = points[2] - points[0], c = cross_product(a, b);
        k = c.x * points[0].x + c.y * points[1].y + c.z * points[2].z;
    }catch(...){std::cout << "error: number of points to polygon's constructor must be greater than two\n";}
    void stretch(double& greatest, point& p){
        double factor = greatest / p.z;
        p.x *= factor;
        p.y *= factor;
        p.z = greatest;
    };
    std::pair<bool, double> hit(vector u) override {
        double e = u.x * c.x + u.y * c.y + u.z * c.z;
        double t = k / e;
        if(!e || (u.z * t < 0)) return std::pair<bool, double>(false, 0);
        point intersection(u.x * t, u.y * t, u.z * t);
        double distance = sqrt(intersection.x * intersection.x + intersection.y * intersection.y + intersection.z * intersection.z);
        double greatest = greater(greater(points[0].z, points[1].z), points[2].z);

        for(point& i : points) stretch(greatest, i);
        stretch(greatest, intersection);

        #define tri_specialization
        #define first_algorithm

        #ifdef tri_specialization
        if(points.size() == 3){
            #ifdef first_algorithm
                if(points[2].y == points[0].y){
                    point temp = points[2];
                    points[2] = points[1];
                    points[1] = temp;
                }
                double s1 = points[2].y - points[0].y,
                       s2 = points[2].x - points[0].x,
                       s3 = points[1].y - points[0].y,
                       s4 = intersection.y - points[0].y;
                double w1 = (points[0].x * s1 + s4 * s2 - intersection.x * s1) / (s3 * s2 - (points[1].x - points[0].x) * s1),
                       w2 = (s4 - w1 * s3) / s1;
                return std::pair<bool, double>(w1 >= 0 && w2 >= 0 && (w1 + w2) <= 1, distance);
            #else
                auto sign = [&](point p1, point p2, point p3) -> double{
                    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
                };
                double d1 = sign(p, points[0], points[1]),
                       d2 = sign(p, points[1], points[2]),
                       d3 = sign(p, points[2], points[0]);
                bool neg = d1 < 0 || d2 < 0 || d3 < 0,
                     pos = d1 > 0 || d2 > 0 || d3 > 0;
                return std::pair<bool, double>(!(neg && pos), distance);
            #endif
        }
        #endif
        bool inside = false;
        for(int i = 0, j = points.size() - 1; i < points.size(); j = i++)
            if(((points[i].y > intersection.y) != (points[j].y > intersection.y)) &&
                (intersection.x < (points[j].x - points[i].x) * (intersection.y - points[i].y) / (points[j].y - points[i].y) + points[i].x))
                    inside = !inside;
        return std::pair<bool, double>(inside, distance);
    }
private:
    vector a, b, c;
    double k;
};
sphere s(RGB(0, 0, 255), sqrt(100000), point(100, 41.5, 2000));
point a(-250, -300, 2400), b(250, 200, 1500), c(350, -100, 1000)
#ifdef QUAD
    , d(-250, -100, 2000)
#endif
;
polygon p(RGB(0, 255, 0), a, b, c
#ifdef QUAD
    , d
#endif
);
std::vector<object*> world = {&s, &p};
std::function<void()> render =
[&](){
    ++f;
    double z = width / (2 * tan(horizontal_fov / 2));
    for(int j = height - 1; j >= 0; --j){
        for(int i = 0; i < width; ++i){
            vector v(-width / 2.0 + i, height / 2.0 - j, z);

            //std::pair<bool, double> hit_s = s.hit(v);
            //bool hit_sphere = hit_s.first;
            //double distance_to_sphere = hit_s.second;

            //std::pair<bool, double> hit_p = p.hit(v);
            //bool hit_polygon = hit_p.first;
            //double distance_to_polygon = hit_p.second;

            // for multiple objects create a list of objects the ray hit
            // with either inheritance or std::any or smth
            // and find which of the distances is least
            bool hit_nothing = true;
            object* o = *std::min_element(world.begin(), world.end(), [&hit_nothing](object* a, object *b){
                // bagel supremacy
                std::pair<bool, double> hit_a = a.hit(v), hit_b = b.hit(v);
                /*
                if(unlikely(hit_a.first && hit_b.first)){
                    hit_nothing = false;
                    return hit_a.second < hit_b.second;
                }
                if(unlikely(hit_a.first)){
                    hit_nothing = false;
                    return true;
                }
                //if(unlikely(hit_b.first)) hit_nothing = false;
                hit_nothing = !hit_b.first;
                return false;
                */
                if(unlikely(hit_a.first)){
                    hit_nothing = false;
                    if(unlikely(hit_b.first)) return hit_a.second < hit_b.second;
                    return true;
                }
                //if(unlikely(hit_b.first)) hit_nothing = false;
                hit_nothing = !hit_b.first;
                return false;
            });
            if(likely(hit_nothing)) framebuf[j * width + i] = RGB(255, 0, 0);
            else framebuf[j * width + i] = o->color;

            //if(unlikely(hit_sphere || hit_polygon)) framebuf[j * width + i] = RGB(0, 0, 255);
            //else framebuf[j * width + i] = RGB(255, 0, 0);
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
