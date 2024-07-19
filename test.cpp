#include <iostream>
#include <cmath>
#include <vector>
/*
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
*/
struct point{
    double x, y, z;
	bool is_normalized = 0;
    point() : x(0), y(0), z(0){}
    point(double e, double f, double g) : x(e), y(f), z(g){}
	double magnitude_squared() const{
		return x * x + y * y + z * z;
	}
	double magnitude() const{
		return sqrt(this->magnitude_squared());
	}
	void set_normalized(const bool a){is_normalized = a;}
	void normalize(){
		set_normalized(true);
		const double magnitude = this->magnitude();
		this->x /= magnitude, this->y /= magnitude, this->z /= magnitude;
	}
	inline point operator-() const{
		return point(-x, -y, -z);
	}
    inline point operator-(const point& f) const{
        return point(x - f.x, y - f.y, z - f.z);
    }
    inline point& operator-=(const point& g){
        x -= g.x;
        y -= g.y;
        z -= g.z;
        return *this;
    }
	inline point& operator+=(const point& o){
		x += o.x;
		y += o.y;
		z += o.z;
		return *this;
	}
    inline point operator*(const double h) const{
        return point(x * h, y * h, z * h);
    }
	inline point operator/(const double i) const{
		return point(x / i, y / i, z / i);
	}
	friend std::ostream& operator<<(std::ostream& o, const point& b){
		o << '{' << b.x << ", " << b.y << ", " << b.z << '}';
		return o;
	}
	void rotate(const double y, /*const*/ double p, const double r){
		//test
		p = -p;
		//todo: cache previous values of sin(angle) and cos(angle), maybe save the number rounded to the thousandths place, and round the next one to the same place before looking it up
		//double sin_yar = sin(y), cos_yar = cos(y),
		//	   sin_par = sin(p), cos_par = cos(p),
		//	   sin_rar = sin(r), cos_rar = cos(r);
		double sin_yar, cos_yar, sin_par, cos_par, sin_rar, cos_rar;
		if(y) sin_yar = sin(y), cos_yar = cos(y);
		else sin_yar = 0, cos_yar = 1;
		if(p) sin_par = sin(p), cos_par = cos(p);
		else sin_par = 0, cos_par = 1;
		if(r) sin_rar = sin(r), cos_rar = cos(r);
		else sin_rar = 0, cos_rar = 1;
		const double x_copy = this->x, y_copy = this->y;
		this->x = (cos_yar * cos_rar + sin_yar * sin_par * sin_rar) * this->x
				+ (sin_yar * sin_par * cos_rar - cos_yar * sin_rar) * this->y
				+ (sin_yar * cos_par) * this->z;
		this->y = (cos_par * sin_rar) * x_copy
				+ (cos_par * cos_rar) * this->y
				- sin_par * this->z;
		this->z = (cos_yar * sin_par * sin_rar - sin_yar * cos_rar) * x_copy
				+ (sin_yar * sin_rar + cos_yar * sin_par * cos_rar) * y_copy
				+ (cos_yar * cos_par) * this->z;
	}
	void rotate_yaw(const double y){
		//todo: cache previous values of sin(angle) and cos(angle)
		if(y){
			double cos_yar = cos(y), sin_yar = sin(y), sin_yar_x = sin_yar * this->x;
			this->x = cos_yar * this->x + sin_yar * this->z;
			this->z = -sin_yar_x + cos_yar * this->z;
		}
	}
	void rotate_pitch(/*const*/ double p){
		//todo: cache previous values of sin(angle) and cos(angle)
		if(p){
			//test
			p = -p;
			//todo: check if pitch is actually correct (make it so positive always migrates to top instead of going in circles)
			double cos_par = cos(p), sin_par = sin(p), sin_par_y = sin_par * this->y;
			std::cout << "point " << *this << '\n';
			std::cout << "cos_par: " << cos_par << "; sin_par: " << sin_par << '\n';
			this->y = cos_par * this->y - sin_par * this->z;
			std::cout << "this->y: " << this->y << '\n';
			this->z = sin_par_y + cos_par * this->z;
		}
	}
	void rotate_roll(const double r){
		//todo: cache previous values of sin(angle) and cos(angle)
		if(r){
			//todo: check if roll is actually correct
			double cos_rar = cos(r), sin_rar = sin(r), sin_rar_x = sin_rar * this->x;
			this->x = cos_rar * this->x - sin_rar * this->y;
			this->y = sin_rar_x + cos_rar * this->y;
		}
	}
	/* if needed
	#if __cplusplus >= 202002L
	auto operator<=>(const point& other) const{
		//this is not correct
		if(x == other.x){
			if(y == other.y){
				if(z == other.z)
					return 0;
				if(z < other.z) return -1;
				return 1;
			}
			if(y < other.y) return -1;
			return 1;
		}
		if(x < other.x) return -1;
		return 1;
	}
	#else
	bool operator==(const point& other) const{
		return x == other.x && y == other.y && z == other.z;
	}
	#if __cplusplus < 202002L
	bool operator!=(const point& other) const{
		return !this->operator==(other);
	}
	#endif
	//implement >, <, >=, <= if needed
	#endif
	*/
};
/*
std::ostream& operator<<(std::ostream& o, const point& b){
	o << '{' << b.x << ", " << b.y << ", " << b.z << '}';
	return o;
}
*/
typedef std::vector<point> plane;
//typedef std::vector<plane> plane_array;
typedef point vector;
inline vector cross_product(const vector& a, const vector& b){
    return vector(a.y * b.z - b.y * a.z, a.z * b.x - b.z * a.x, a.x * b.y - b.x * a.y);
}
inline double dot_product(const vector& a, const vector& b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
bool on_plane_normal_side(const plane& bound, const point& p, int plane_number = 0){
	vector a = bound[1] - bound[0], b = bound[2] - bound[0], cross = cross_product(a, b);
	double A = -cross.x * bound[0].x, B = -cross.y * bound[0].y, C = -cross.z * bound[0].z;
	//todo: idk if this is correct
	return A * p.x + B * p.y + C * p.z <= 0;
}
double plane_distance(const plane& bound, const point& p, int plane_number = 0){ 
	vector a = bound[1] - bound[0], b = bound[2] - bound[0], cross = cross_product(a, b);
	std::cout << plane_number << " cross: " << cross << '\n';
	double cross_mag = cross.magnitude();//sqrt(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z);
	cross.x /= cross_mag, cross.y /= cross_mag, cross.z /= cross_mag;
	std::cout << "normalized cross: " << cross << '\n';
	std::cout << "bound[0]: " << bound[2] << '\n';
	return dot_product(cross, p);
	//vector v = bound[0] - p;
	//return dot_product(v, cross);
	//return (dot_product(cross, p) - dot_product(cross, bound[0])) / cross_mag;
}
int main(){
	//plane p{point(0, 0, 0), point(0, 0, -50), point(0, -50, 0)};
	constexpr int width_px = 1200, height_px = 600;
	constexpr double ratio = height_px / (double)width_px;
	constexpr double width = 1, height = ratio;
	//constexpr double z = 0.497903;
	constexpr double convert_horizontal_fov_to_radians = 0.0175;
	constexpr double horizontal_fov_degrees = 90;
    constexpr double z = width / (2 * tan(horizontal_fov_degrees * convert_horizontal_fov_to_radians / 2));
	/*
	plane p{
        point(-width / 2, -height / 2, z),
        point(width / 2, -height / 2, z),
        point(0, 0, 0)
    };
	*/
	constexpr double yaw_angle_radians = 0, pitch_angle_radians = -0.6, roll_angle_radians = 0;
	point bounding[4]{
		point(-width / 2.0, height / 2.0, z), //upper left
		point(width / 2.0, height / 2.0, z), //upper right
		point(-width / 2.0, -height / 2.0, z), //bottom left
		point(width / 2.0, -height / 2.0, z) //bottom right
	};
	std::cout << "bounding:\n";
	for(auto& i : bounding) std::cout << i << '\n';
	for(point& i : bounding)
		//ROTATE BY NEGATIVE pitch_angle_radians
		//i.rotate(yaw_angle_radians, -pitch_angle_radians, roll_angle_radians);
		i.rotate_pitch(-pitch_angle_radians);
	std::cout << "after rotation:\n";
	for(auto& i : bounding) std::cout << i << '\n';
	/*
	plane p{
		bounding[3], 
		bounding[2], 
		bounding[3] * 2
	};
	*/
	std::cout << "bottom y: " << bounding[2].y << '\n';
	std::cout << "top y: " << bounding[0].y << '\n';
	std::vector<plane> planes = {
		{bounding[0], bounding[2], point(0, 0, 0)}, //left
		{bounding[3], bounding[1], point(0, 0, 0)}, //right
		{bounding[1], bounding[0], point(0, 0, 0)}, //top
		{bounding[3], bounding[2], bounding[3] * 2}, //todo: fix bottom
	};
	//for(int p = 0; auto& i : planes[2])
	//	i.y = -i.y, std::cout << "plane[" << p << "]: " << i << '\n', ++p;
	//constexpr double _x = 350, _y = -533.579, _z = 2809.41;
	constexpr double _x = 0, _y = -1400, _z = 2400;
	point pt = point(_x, _y, _z);
	point base = point(0, 0, 0);
	std::cout << "pt: " << pt << '\n';
	for(int i = 0; const auto& p : planes) std::cout << i << ":\n" << std::boolalpha << /*plane_distance(p, point(_x, _y, _z), i)*/ on_plane_normal_side(p, pt, i) << '\n', ++i;
	std::cout << "base: " << base << '\n';
	for(int i = 0; const auto& p : planes) std::cout << i << ":\n" << std::boolalpha << on_plane_normal_side(p, base, i) << '\n', ++i;
	//std::cout << plane_distance(p, point(350, -974.321, 2506.49)) << '\n';
	//std::cout << plane_distance(p, point(350, -533.579, 2809.41)) << '\n';
	//std::cout << plane_distance(p, point(50, 0, 0)) << '\n';
    //std::cout << plane_distance(p, point(-240, -320, 2400)) << '\n';
    return 0;
}
