#ifndef ray_tracing_hpp
#define ray_tracing_hpp
#include <memory>
#include <vector>
#include <limits>
#include <cmath>
constexpr double infinity = std::numeric_limits<double>::infinity();
class vec3{
public:
    double x, y, z;
    vec3() : x(0), y(0), z(0){}
    vec3(double a, double b, double c) : x(a), y(b), z(c){};
    vec3 operator-() const{
        return vec3{-x, -y, -z};
    }
    vec3& operator+=(const vec3& v){
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    vec3& operator*=(double t){
        x *= t;
        y *= t;
        z *= t;
        return *this;
    }
    vec3& operator/=(double t){
        return *this *= 1 / t;
    }
    double length() const{
        return sqrt(length_squared());
    }
    double length_squared() const{
        return x * x + y * y + z * z;
    }
};
using point3 = vec3;
inline vec3 operator+(const vec3& a, const vec3& b){
    return vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}
inline vec3 operator+(const vec3& a, double t){
    return vec3{a.x + t, a.y + t, a.z + t};
}
inline vec3 operator+(double t, const vec3& a){
    return vec3{a.x + t, a.y + t, a.z + t};
}
inline vec3 operator-(const vec3& a, const vec3& b){
    return vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}
inline vec3 operator*(const vec3& a, const vec3& b){
    return vec3{a.x * b.x, a.y * b.y, a.z * b.z};
}
inline vec3 operator*(double t, const vec3& v){
    return vec3{t * v.x, t * v.y, t * v.z};
}
inline vec3 operator*(const vec3& v, double t){
    return vec3{t * v.x, t * v.y, t * v.z};
}
inline vec3 operator/(const vec3& v, double t){
    return (1 / t) * v;
}
inline double dot(const vec3& a, const vec3& b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline vec3 cross(const vec3& a, const vec3& b){
    return vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline vec3 unit_vector(const vec3& v){
    return v / v.length();
}
class ray{
public:
    point3 origin;
    vec3 direction;
    ray(){}
    ray(const point3& o, const vec3& d) : origin(o), direction(d){}
    point3 at(double t) const{
        return origin + t * direction;
    }
};
class hit_record{
public:
    uint32_t color;
    point3 p;
    vec3 normal;
    double t;
    bool front;
    void set_face_normal(const ray& r, const vec3& out){
        front = dot(r.direction, out) < 0;
        normal = (front - !front) * out;
        //normal = front ? out : -out;
    }
};
class object{
public:
    uint32_t color;

    virtual bool hit(const ray& r, double ray_tmin, double ray_tmax, hit_record& rec) const = 0;
    virtual ~object() = default;
};
using object_ptr = std::shared_ptr<object>;
class sphere : public object{
public:
    uint32_t color;

    point3 center;
    double radius;
    sphere(){}
    sphere(point3 c, double r, uint32_t clr) : center(c), radius(r), color(clr){}
    static std::shared_ptr<sphere> create(point3 c, double r, uint32_t clr){
        return std::make_shared<sphere>(c, r, clr);
    }
    bool hit(const ray& r, double ray_tmin, double ray_tmax, hit_record& rec) const override{
        vec3 o = r.origin - center;
        double a = r.direction.length_squared(),
               half_b = dot(o, r.direction),
               c = o.length_squared() - radius * radius,
               discriminant = half_b * half_b - a * c;
        if(discriminant < 0) return 0;
        double sqrt_d = sqrt(discriminant),
               root = (-half_b - sqrt_d) / a;
        if(root <= ray_tmin || ray_tmax <= root){
            root = (-half_b + sqrt_d) / a;
            if(root <= ray_tmin || ray_tmax <= root) return 0;
        }

        rec.t = root;
        rec.p = r.at(rec.t);
        rec.set_face_normal(r, vec3{(rec.p - center) / radius});
        rec.color = color;

        return 1;
    }
};
using sphere_ptr = std::shared_ptr<sphere>;
class plane : public object{
public:
    uint32_t color;

    bool hit(const ray& r, double ray_tmin, double ray_tmax, hit_record& rec) const override{

    }
};
class object_list : public object{
public:
    std::vector<object_ptr> objects;
    object_list(){}
    object_list(object_ptr a){
        objects.emplace_back(a);
    }
    void clear(){
        objects.clear();
    }
    void add(object_ptr a){
        objects.emplace_back(a);
    }
    bool hit(const ray& r, double ray_tmin, double ray_tmax, hit_record& rec) const override{
        hit_record h;
        bool has_hit = 0;
        double closest = ray_tmax;
        for(const object_ptr& o : objects){
            if(o->hit(r, ray_tmin, closest, h)){
                has_hit = 1;
                closest = h.t;
                rec = h;
            }
        }
        return has_hit;
    }
    object_ptr& operator[](int i){
        return objects[i];
    }
    object_ptr operator[](int i) const{
        return objects[i];
    }
};
class camera{
public:
    double viewport_height, viewport_width, offset_x, offset_y;
    window w;
    point3 c;
    camera(window e) : w(e), c(0, 0, 0), viewport_height(2){
        viewport_width = ((w.width) << 1) / (double)w.height;
        offset_x = viewport_width / w.width;
        offset_y = viewport_height / w.height;
        w.on_resize([&](){
            viewport_height = 2, viewport_width = ((w.width) << 1) / (double)w.height,
            offset_x = viewport_width / w.width, offset_y = viewport_height / w.height;
        });
    }
    void render(const object_list& world){
        static vec3 current{-viewport_width / 2, 1, -1};
        for(int i = 0, e = w.height; --e; current.y -= offset_y){
            for(int f = 0; f < w.width; ++f, ++i, current.x += offset_x){
                ray r{c, current};
                hit_record h;
                w.pixels[i] = world.hit(r, 0, infinity, h) ? h.color : rgb(255, 255, 255);
            }
            current.x = -viewport_width / 2;
        }
        current.y = 1;
    }
};
#endif