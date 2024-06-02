#include <windows.h>
#include <process.h>
#include <functional>
#include <algorithm>
#include <iostream>
#include <atomic>
//#include <memory>
#include <vector>
#include <cmath>
#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
std::atomic<int> f(0);
std::atomic<double> yaw_angle_radians(0);
std::atomic<double> pitch_angle_radians(0);
std::atomic<double> roll_angle_radians(0);
void fps(void*){
    while(1){
        Sleep(1000);
        std::cout << f << '\n';
        f = 0;

        //yaw_angle_radians = yaw_angle_radians + 0.175;
        //pitch_angle_radians = (!pitch_angle_radians) * 0.175;
        pitch_angle_radians = pitch_angle_radians + 0.175;
        //roll_angle_radians = (!roll_angle_radians) * 0.35;
    }
}
int width = 1200, height = 600, tick_count, *framebuf;
double horizontal_fov = 120, vertical_fov = atan(tan(height / 2.0) * height / width);
inline double square(double x){
    return x * x;
}
inline double greater(double a, double b){
    return a > b ? a : b;
}
inline double abs_val(double a){
    return a < 0 ? -a : a;
}
inline double dot_product(double ax, double ay, double az, double bx, double by, double bz){
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
    point operator*(double h){
        return point(x * h, y * h, z * h);
    }
};
typedef uint32_t color;
typedef point vector;
vector cross_product(vector a, vector b){
    return vector(a.y * b.z - b.y * a.z, a.z * b.x - b.z * a.x, a.x * b.y - b.x * a.y);
}
inline double dot_product(vector a, vector b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
struct object{
    virtual std::pair<bool, double> hit(vector u) = 0;
    virtual bool is_in_frustum(std::vector<std::vector<point>> plane) = 0;
    virtual void set_color(color c){clr = c;}
    color clr;
};
struct sphere : public object{
    double radius;
    point center;
    sphere(color cl, double r, point c) : radius(r), center(c){set_color(cl);}
    virtual void set_color(color c){clr = c;}
    std::pair<bool, double> hit(vector u) override{
        //this will be in render instead
        //todo: do the same thing for u.x and u.y
        //todo: also fast reject if sphere is not in field of view
        //here for funny
        if((u.z > 0 && center.z < radius) || (u.z < 0 && center.z > radius)) return std::pair<bool, double>(false, 0);
        //magnitude(cross(A - B, C - B)) 
        //divided by magnitude(C - B)
        double magnitude_squared = u.x * u.x + u.y * u.y + u.z * u.z;
        double dot = dot_product(u.x, u.y, u.z, center.x, center.y, center.z);
        double determinant = dot * dot - magnitude_squared * (center.x * center.x + center.y * center.y + center.z * center.z - radius * radius);
        // if determinant < 0 you can't sqrt it
        if(determinant < 0) return std::pair<bool, double>(false, 0);
        return std::pair<bool, double>(true, (-dot - sqrt(determinant)) / magnitude_squared);
    }
    bool is_in_frustum(std::vector<std::vector<point>> plane) override{
        //FOR A LINE
        //magnitude(cross(A - B, C - B)) 
        //divided by magnitude(C - B)
        //so let j = C - B, sqrt(cross(A - B, j)^2 / (j.x * j.x + j.y * j.y + j.z * j.z))
        
        //if( < radius * radius)

        for(int c = 0; c < (int)plane.size(); ++c){
            /*
            for(point& i : plane[c]) i -= center;

            vector a = plane[c][1] - plane[c][0], b = plane[c][2] - plane[c][0];
            // literally just the cross product btw
            double i = a.y * b.z - a.z * b.y,
                   j = a.x * b.z - a.z * b.x,
                   k = a.x * b.y - a.y * b.x;
            double d = i * -plane[c][0].x + j * -plane[c][0].y + k * -plane[c][0].z;

            //todo: its always printing the same result
            std::cout << d / sqrt(i * i + j * j + k * k) << '\n';
            if(d / sqrt(i * i + j * j + k * k) < radius) return true;
            */

            vector a = plane[c][1] - plane[c][0], b = plane[c][2] - plane[c][0], cross = cross_product(a, b);

            double cross_mag = sqrt(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z);
            cross.x /= cross_mag, cross.y /= cross_mag, cross.z /= cross_mag;

            //double distance = dot_product(cross, center);
            //std::cout << c << ": " << distance << '\n';

            //if its positive its facing the same direction as the cross product
            //if its negative its facing the opposite
            //right hand rule for cross product btw

            //erm probably remove this
			//double t = -cross.x * center.x - cross.y * center.y - cross.z * center.z;
			//if(t * cross.x - center.x < radius){std::cout << "tru moo flat\n"; return true;}

			double distance = dot_product(cross, center);
			//if(distance < 0 || distance < radius){std::cout << "tru moo\n"; return true;}
            if(abs(distance) < radius){std::cout << "tru moo: " << distance << '\n'; return true;}
        }
        return false;

        //bool in = d / (i * i + j * j + k * k) < radius;
        //return d / (i * i + j * j + k * k) < radius; //for now
        //return in;
        //return true;
    }
};
struct polygon : public object{
    std::vector<point> points{};
    polygon(color cl, auto... l) try /*: clr(cl)*/{
        set_color(cl);
        if(sizeof...(l) < 3) throw;
        points = {(point(l))...};
        a = points[1] - points[0], b = points[2] - points[0], c = cross_product(a, b);
        k = c.x * points[0].x + c.y * points[1].y + c.z * points[2].z;
    }catch(...){std::cout << "error: number of points to polygon's constructor must be greater than two\n";}
    virtual void set_color(color c){clr = c;}
    void stretch(double& greatest, point& p){
        double factor = greatest / p.z;
        p.x *= factor;
        p.y *= factor;
        p.z = greatest;
    };
    std::pair<bool, double> hit(vector u) override{
        //this will be in render instead
        //todo: do the same thing for u.x and u.y
        //todo: also fast reject if polygon is not in field of view
        //here for funny
        for(point& i : points)
            if((u.z < 0 && i.z < 0) || (u.z > 0 && i.z > 0)) goto polygon_start;
        return std::pair<bool, double>(false, 0);
        polygon_start:

        double e = u.x * c.x + u.y * c.y + u.z * c.z;
        double t = k / e;
        if(likely(!e || (u.z * t < 0))) return std::pair<bool, double>(false, 0);
        point intersection(u.x * t, u.y * t, u.z * t);

        double distance = sqrt(intersection.x * intersection.x + intersection.y * intersection.y + intersection.z * intersection.z);
        double greatest = greater(greater(points[0].z, points[1].z), points[2].z);

        for(point& i : points) stretch(greatest, i);
        stretch(greatest, intersection);

        //#define tri_specialization
        //#define first_algorithm

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
                     pos = d1 > 0 || d2 > 0 || d3 > 0,
                     ret = !(neg && pos);
                return std::pair<bool, double>(ret, distance);
            #endif
        }
        #endif
        bool inside = false;
        for(int i = 0, j = points.size() - 1; i < (int)points.size(); j = i++)
            if(((points[i].y > intersection.y) != (points[j].y > intersection.y)) &&
                (intersection.x < (points[j].x - points[i].x) * (intersection.y - points[i].y) / (points[j].y - points[i].y) + points[i].x))
                    inside = !inside;
        return std::pair<bool, double>(inside, distance);
    }
    bool is_in_frustum(std::vector<std::vector<point>> plane) override{
        return true; // for now
    }
private:
    vector a, b, c;
    double k;
};

struct doughnut : public object{
    point center;
    double minor_radius, major_radius;
    double pitch_rad, roll_rad;
    doughnut(color cl, double minor_r, double major_r, point c, double pitch = 0, double roll = 0) :
        center(c), minor_radius(minor_r), major_radius(major_r), pitch_rad(pitch), roll_rad(roll), epsilon(major_r * major_r - minor_r * minor_r){set_color(cl);}
    virtual void set_color(color c){clr = c;}
    std::pair<bool, double> hit(vector u) override{
        u -= center;
        //todo: rotate u by the matrix which is the product of the negative pitch_rad matrix and the negative roll_rad matrix

        //this will be in render
        //todo: then instantly return false if donut is not in field of view
        double magnitude = sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
        double big = 10;
        u.x /= (magnitude * big), u.y /= (magnitude * big), u.z /= (magnitude * big);
        long double a, b, c, d, e;
        {
            long double _a = -center.x, _b = u.x, _c = -center.y, _d = u.y, _e = -center.z, _f = u.z;

            #define SIMPLIFY_DOUGHNUT

            #ifdef SIMPLIFY_DOUGHNUT
                a = square(_b * _b + _d * _d + _f * _f);
                b = 4*(_a*_b*_b*_b + _a*_b*_d*_d + _b*_b*_c*_d + _c*_d*_d*_d + _e*_f*_f*_f + _a*_b*_f*_f + _e*_b*_b*_f + _c*_d*_f*_f + _e*_d*_d*_f);
                c = 2*(3*_a*_a*_b*_b + _a*_a*_d*_d + 4*_a*_b*_c*_d + _b*_b*_c*_c + 3*_c*_c*_d*_d * 3*_e*_e*_f*_f + _a*_a*_f*_f + 4*_e*_a*_b*_f + _e*_e*_b*_b + _c*_c*_f*_f + 4*_e*_c*_d*_f + _e*_e*_d*_d + epsilon*_b*_b + epsilon*_d*_d + epsilon*_f*_f - 2*major_radius*major_radius*_b*_b - 2*major_radius*major_radius*_d*_d);
                d = 4*(_a*_a*_a*_b + _a*_a*_c*_d + _a*_b*_c*_c + _c*_c*_c*_d + _e*_e*_e*_f + _e*_e*_a*_b + _e*_c*_c*_f + _e*_e*_c*_d + epsilon*_a*_b + epsilon*_c*_d + _e*epsilon*_f - 2*major_radius*major_radius*_a*_b - 2*major_radius*major_radius*_c*_d);
                e = _a*_a*_a*_a + 2*_a*_a*_c*_c + _c*_c*_c*_c + _e*_e*_e*_e + epsilon*epsilon + 2*_e*_e*_a*_a + 2*_e*_e*_c*_c + 2*epsilon*_a*_a + 2*epsilon*_c*_c + 2*_e*_e*epsilon - 4*major_radius*major_radius*_a*_a - 4*major_radius*major_radius*_c*_c;
            #else
                a = _b*_b*_b*_b + 2*_b*_b*_d*_d + _d*_d*_d*_d + _f*_f*_f*_f + 2*_b*_b*_f*_f + 2*_d*_d*_f*_f;
                b = 4*_a*_b*_b*_b + 4*_a*_b*_d*_d + 4*_b*_b*_c*_d + 4*_c*_d*_d*_d + 4*_e*_f*_f*_f + 4*_a*_b*_f*_f + 4*_e*_b*_b*_f + 4*_c*_d*_f*_f + 4*_e*_d*_d*_f;
                c = 6*_a*_a*_b*_b + 2*_a*_a*_d*_d + 8*_a*_b*_c*_d + 2*_b*_b*_c*_c + 6*_c*_c*_d*_d * 6*_e*_e*_f*_f + 2*_a*_a*_f*_f + 8*_e*_a*_b*_f + 2*_e*_e*_b*_b + 2*_c*_c*_f*_f + 8*_e*_c*_d*_f + 2*_e*_e*_d*_d + 2*epsilon*_b*_b + 2*epsilon*_d*_d + 2*epsilon*_f*_f - 4*major_radius*major_radius*_b*_b - 4*major_radius*major_radius*_d*_d;
                d = 4*_a*_a*_a*_b + 4*_a*_a*_c*_d + 4*_a*_b*_c*_c + 4*_c*_c*_c*_d + 4*_e*_e*_e*_f + 4*_e*_e*_a*_b + 4*_e*_c*_c*_f + 4*_e*_e*_c*_d + 4*epsilon*_a*_b + 4*epsilon*_c*_d + 4*_e*epsilon*_f - 8*major_radius*major_radius*_a*_b - 8*major_radius*major_radius*_c*_d;
                e = _a*_a*_a*_a + 2*_a*_a*_c*_c + _c*_c*_c*_c + _e*_e*_e*_e + epsilon*epsilon + 2*_e*_e*_a*_a + 2*_e*_e*_c*_c + 2*epsilon*_a*_a + 2*epsilon*_c*_c + 2*_e*_e*epsilon - 4*major_radius*major_radius*_a*_a - 4*major_radius*major_radius*_c*_c;
            #endif
        }
        b /= a;
        c /= a;
        d /= a;
        e /= a;
        //std::cout << a << ' ' << b << ' ' << c << ' ' << d << ' ' << e << ' ' << '\n';
        long double discriminant = 256 * e * e * e - 192 * b * d * e * e - 128 * c * c * e * e + 144 * c * d * d * e - 27 * d * d * d *d + 144 * b * b * c * e * e - 6 * b * b * d * d * e - 80 * b * c * c * d * e + 18 * b * c * d * d * d + 16 * c * c * c * c * e - 4 * c * c * c * d * d - 27 * b * b * b * b * e * e + 18 * b * b * b * c * d * e - 4 * b * b * b * d * d * d - 4 * b * b * c * c * c * e + b * b * c * c * d * d;
        long double p = 8 * c - 3 * b * b;
        long double g = 64 * e - 16 * c * c + 16 * b * b * c - 16 * b * d - 3 * b * b * b * b;
        //todo: find the distance instead of just putting 100
        if(discriminant < 0) return std::pair<bool, double>(true, 100);
        if(!discriminant) return std::pair<bool, double>(!(!g && p > 0 && b * b * b + 8 * d - 4 * b * c), 100);
        //if(discriminant > 0)
        //return std::pair<bool, double>(!(p > 0 || g > 0), 100);
        return std::pair<bool, double>(p < 0 && g < 0, 100);
    }
    bool is_in_frustum(std::vector<std::vector<point>> plane) override{
        return true; //for now
    }
private:
    double epsilon;
};

sphere s(RGB(0, 0, 255), sqrt(100000), point(100, 41.5, 2000));
//doughnut d(RGB(0, 255, 0), 200, 300, point(100, 100, 5000));
doughnut d(RGB(255, 0, 255), 60, 90, point(0, 0, 0));
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
std::vector<object*> world = {&s, &p, &d};
double z = width / (2 * tan(horizontal_fov / 2));
std::function<void()> render =
[](){
    //todo: use the gpu for calculations
    //magnitude(cross(A - B, C - B)) 
    //divided by magnitude(C - B)
	std::vector<object*> can_hit = world; // copy
    point bounding[4]{
        point(-width / 2.0, height / 2.0, z), //upper left
        point(width / 2.0, height / 2.0, z), //upper right
        point(-width / 2.0, -height / 2.0, z), //bottom left
        point(width / 2.0, -height / 2.0, z) //bottom right
    };
    //todo: rotate the bounding points by yaw pitch and roll, during initialization or after

    for(point& i : bounding){
        double &vx = i.x, &vy = i.y, &vz = i.z;
        if(yaw_angle_radians){
            double cos_yar = cos(yaw_angle_radians), sin_yar = sin(yaw_angle_radians), sin_yar_vx = sin_yar * vx;
            vx = cos_yar * vx + sin_yar * vz;
            vz = -sin_yar_vx + cos_yar * vz;
        }
        if(pitch_angle_radians){
            //double temp = vx;
            //vx = cos(pitch_angle_radians) * vx + sin(pitch_angle_radians) * vz;
            //vz = -sin(pitch_angle_radians) * temp + cos(pitch_angle_radians) * vz;
            //vy += vz * sin(pitch_angle_radians);

            double cos_par = cos(pitch_angle_radians), sin_par = sin(pitch_angle_radians), sin_par_vy = sin_par * vy;
            vy = cos_par * vy - sin_par * vz;
            vz = sin_par_vy + cos_par * vz;

            //todo: update pitch to something instead of this
            //vy += vz * sin(pitch_angle_radians);
        }
        if(roll_angle_radians){
            //double temp = vy;
            //vy = cos(roll_angle_radians) * vy - sin(roll_angle_radians) * vz;
            //vz = sin(roll_angle_radians) * temp + cos(roll_angle_radians) * vz;
            double cos_rar = cos(roll_angle_radians), sin_rar = sin(roll_angle_radians), sin_rar_vx = sin_rar * vx;
            vx = cos_rar * vx - sin_rar * vy;
            vy = sin_rar_vx + cos_rar * vy;
        }
    }
    constexpr double scaler = 0.5; // any number that isn't 1, maybe just change it to 0 idk
    std::vector<std::vector<point>> plane{
        {bounding[0], bounding[2], bounding[0] * scaler}, //left
        {bounding[1], bounding[3], bounding[1] * scaler}, //right
        {bounding[0], bounding[1], bounding[0] * scaler}, //top
        {bounding[2], bounding[3], bounding[2] * scaler}, //bottom
    };
    //todo: in is_in_frustum() for every object, find the distance from either the center or all points to each plane, and check if any part of the shape is inside
    #if __cplusplus >= 202002L
        std::erase_if(can_hit, [&plane](object* i){return !i->is_in_frustum(plane);});
    #else
        can_hit.erase(std::remove_if(can_hit.begin(), can_hit.end(), [&plane](object* i){return !i->is_in_frustum(plane);}), can_hit.end());
    #endif

    ++f;
    //todo: start from upper left instead
    //todo: increment vx, vy, and vz by the appropriate amount starting from the top left
    for(int j = height - 1; j >= 0; --j){
        for(int i = 0; i < width; ++i){
            double vx = -width / 2.0 + i;
            double vy = height / 2.0 - j;
            double vz = z;
            //todo: its kinda slow ngl
            //todo: remove yaw pitch and roll here
            if(yaw_angle_radians){
                double cos_yar = cos(yaw_angle_radians), sin_yar = sin(yaw_angle_radians), sin_yar_vx = sin_yar * vx;
                vx = cos_yar * vx + sin_yar * vz;
                vz = -sin_yar_vx + cos_yar * vz;
            }
            if(pitch_angle_radians){
                //double temp = vx;
                //vx = cos(pitch_angle_radians) * vx + sin(pitch_angle_radians) * vz;
                //vz = -sin(pitch_angle_radians) * temp + cos(pitch_angle_radians) * vz;
                //vy += vz * sin(pitch_angle_radians);

                double cos_par = cos(pitch_angle_radians), sin_par = sin(pitch_angle_radians), sin_par_vy = sin_par * vy;
                vy = cos_par * vy - sin_par * vz;
                vz = sin_par_vy + cos_par * vz;

                //todo: update pitch to something instead of this
                //vy += vz * sin(pitch_angle_radians);
            }
            if(roll_angle_radians){
                //double temp = vy;
                //vy = cos(roll_angle_radians) * vy - sin(roll_angle_radians) * vz;
                //vz = sin(roll_angle_radians) * temp + cos(roll_angle_radians) * vz;
                double cos_rar = cos(roll_angle_radians), sin_rar = sin(roll_angle_radians), sin_rar_vx = sin_rar * vx;
                vx = cos_rar * vx - sin_rar * vy;
                vy = sin_rar_vx + cos_rar * vy;
            }
            vector v(vx, vy, vz);

            bool hit_nothing = true, small_change = false;
            object* smallest = can_hit[0];
            std::pair<bool, double> hit_b = smallest->hit(v);
            //todo: likely and unlikely might be unnecessary
            for(int w = 1; w < (int)can_hit.size(); ++w){
                //if(comp(hit_nothing, v, world[i], smallest)) smallest = world[i];
                std::pair<bool, double> hit_a = can_hit[w]->hit(v);
                if(small_change) hit_b = smallest->hit(v), small_change = false;
                //std::pair<bool, double> hit_a(0, 0), hit_b(0, 0);

                if(unlikely(hit_a.first)){
                    hit_nothing = false;
                    if(unlikely(hit_b.first)){
                        if(hit_a.second < hit_b.second)
                            smallest = can_hit[w], small_change = true;
                    }
                    else smallest = can_hit[w], small_change = true;
                }
                else hit_nothing = !hit_b.first;
            }
            if(likely(hit_nothing)) framebuf[j * width + i] = RGB(255, 0, 0);
            else framebuf[j * width + i] = smallest->clr;
        }
    }
}
, resize =
[](){
    z = width / (2 * tan(horizontal_fov / 2));
    vertical_fov = atan(tan(height / 2.0) * height / width);
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
    //std::ios_base::sync_with_stdio(false);
    //std::cout.tie(0);
    //std::cin.tie(0);
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
