#include "rendering.hpp"
#include "ray_tracing.hpp"
#include <process.h>
#include <iostream>
int f = 0;
void fps(void*){
    while(1){
        Sleep(1000);
        std::cout << f << '\n';
        f = 0;
    }
}
int main(){
    screen s;
    int width = 800, height = 600;
    window w{"hehehehaw", width, height, (s.width() >> 1) - width / 2, 50};
    //double viewport_height = 2, viewport_width = ((w.width) << 1) / (double)w.height,
    //       offset_x = viewport_width / w.width, offset_y = viewport_height / w.height;
    //point3 camera{0, 0, 0};
    //vec3 current;
    //current.z = -1;
    //w.on_resize([&](){
    //    viewport_height = 2, viewport_width = ((w.width) << 1) / (double)w.height,
    //    offset_x = viewport_width / w.width, offset_y = viewport_height / w.height;
    //});
    object_list world;
    world.add(sphere::create(point3{0, 0, -1}, 0.5, rgb(0, 0, 255)));
    world.add(sphere::create(point3{0, 0, -1}, 0.5, rgb(255, 0, 0)));
    uint8_t counter = 0;
    camera c(w);
    _beginthread(fps, 0, 0);
    w.render([&](){
        ++counter, ++f;
        //current.y = 1;
        /*
        auto hit_sphere = [&](const point3& center, double radius, const ray& r) -> double{
            vec3 o = r.origin - center;
            double a = r.direction.length_squared(),
                   half_b = dot(o, r.direction),
                   c = o.length_squared() - radius * radius,
                   discriminant = half_b * half_b - a * c;
            if(discriminant < 0) return -1;
            return (-half_b - sqrt(discriminant)) / a;
        };
        double g = 0.75 - (w.tick_count() & 1);
        for(int i = 0, e = w.height; --e; current.y -= offset_y){
            current.x = -viewport_width / 2;
            for(int f = 0; f < w.width; ++f, ++i, current.x += offset_x){
                ray r{camera, current};
                w.pixels[i] = hit_sphere(point3{-g, 0, -1}, 0.5, r) > 0.0
                           || hit_sphere(point3{g, 0, -1}, 0.5, r) > 0.0 ? rgb(255, 0, 0) : rgb(255, 255, 255);
            }
        }
        */
        //sphere center_sphere(point3{0, 0, -1}, 0.5);
        double g = (counter & 1) * 0.75;
        std::dynamic_pointer_cast<sphere>(world[0])->center.x = g;
        std::dynamic_pointer_cast<sphere>(world[1])->center.x = -g;
        //for(int i = 0, e = w.height; --e; current.y -= offset_y){
        //    current.x = -viewport_width / 2;
        //    for(int f = 0; f < w.width; ++f, ++i, current.x += offset_x){
        //        ray r{camera, current};
        //        hit_record h;
        //        w.pixels[i] = world.hit(r, 0, infinity, h) ? rgb(255, 0, 0) : rgb(255, 255, 255);
        //    }
        //}
        c.render(world);
    });
    window::loop();
    return 0;
}