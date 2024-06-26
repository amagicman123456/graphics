#include <iostream>
#include <cmath>
#include <vector>
struct point{
    double x, y, z;
    point() : x(0), y(0), z(0){}
    point(double e, double f, double g) : x(e), y(f), z(g){}
	inline point operator-(){
		return point(-x, -y, -z);
	}
    inline point operator-(point f){
        return point(x - f.x, y - f.y, z - f.z);
    }
    inline point& operator-=(point g){
        x -= g.x;
        y -= g.y;
        z -= g.z;
        return *this;
    }
    inline point operator*(double h){
        return point(x * h, y * h, z * h);
    }
	inline point operator/(double i){
		return point(x / i, y / i, z / i);
	}
};
typedef std::vector<point> plane;
//typedef std::vector<plane> plane_array;
typedef point vector;
inline vector cross_product(vector a, vector b){
    return vector(a.y * b.z - b.y * a.z, a.z * b.x - b.z * a.x, a.x * b.y - b.x * a.y);
}
inline double dot_product(vector a, vector b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
double plane_distance(plane bound, point p){ 
	vector a = bound[1] - bound[0], b = bound[2] - bound[0], cross = cross_product(a, b);
	double cross_mag = sqrt(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z);
	cross.x /= cross_mag, cross.y /= cross_mag, cross.z /= cross_mag;
	return dot_product(cross, p);
}
int main(){
	//plane p{point(0, 0, 0), point(0, 0, -50), point(0, -50, 0)};
	constexpr int width_px = 1200, height_px = 600;
	constexpr double ratio = height / (double)width;
	constexpr double width = 1, height = ratio;
	constexpr double z = 0.497903;
    plane p{
        point(-width / 2, -height / 2, z),
        point(width / 2, -height / 2, z),
        point(0, 0, 0)
    };
    //std::cout << plane_distance(p, point(50, 0, 0)) << '\n';
    std::cout << plane_distance(p, point(-240, -320, 2400)) << '\n';
    return 0;
}
