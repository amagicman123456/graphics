#include <windowsx.h>
#include <windows.h>
#include <process.h>
#include <functional>
#include <algorithm>
#include <iostream>
#include <atomic>
#include <vector>
#include <cmath>
#if __cplusplus < 201703L
	#define M_PI 3.14159265358979323846
#endif
#ifdef SOUND
	#include <memory>
	#if __cplusplus >= 201703L
		#include <filesystem>
		namespace fs = std::filesystem;
	#elif __cplusplus >= 201402L
		#include <experimental/filesystem>
		namespace fs = std::experimental::filesystem;
	#endif
#endif
#ifdef IMAGE
#include <fstream>
struct image{
	image() : width(0), row_bytes(0), height(0), data(nullptr){}
	image(image&& other){
		width = other.width;
		row_bytes = other.row_bytes;
		height = other.height;
		data = std::move(other.data);
	}
	image(uint32_t a, int32_t b, std::unique_ptr<uint8_t[]>&& c) : width(a), row_bytes((3 * a + 3) & -4) /* round up to multiple of 4 */, height(b), data(std::move(c)){}
	uint32_t width;
	uint32_t row_bytes;
	int32_t height;
	std::unique_ptr<uint8_t[]> data;
	uint32_t color_at(uint32_t x, uint32_t y) const{
		uint32_t index = /*(height - y - 1)*/y * row_bytes + 3 * x;
		//return (uint32_t)data[index] << 24 |
		//	   (uint32_t)data[index + 1] << 16 |
		//	   (uint32_t)data[index + 2] << 8;
		return (uint32_t)data[index] << 16 |
			   (uint32_t)data[index + 1] << 8 |
			   (uint32_t)data[index + 2];
	}
	image& operator=(image&& other){
		width = other.width;
		row_bytes = other.row_bytes;
		height = other.height;
		data = std::move(other.data);
		return *this;
	}
	//make iterators if wanted
};
/*
#if __cplusplus >= 201703L
#define read_rgb_image_from_buffer(buf) []() constexpr -> image{\
	uint32_t width = *reinterpret_cast<uint32_t*>(buf + 18)\
	int32_t height = *reinterpret_cast<int32_t*>(buf + 22);\
	bool sign = height > 0;\
	uint32_t size = 3 * width * (sign ? height : -height);\
	std::unique_ptr<uint8_t[]> data = std::make_unique<uint8_t[]>(size);\
	bitmap.seekg(*reinterpret_cast<uint32_t*>(buf + 10), std::ios_base::beg);\
	bitmap.read(reinterpret_cast<char*>(data.get()), size);\
	if(sign)\
		for(uint32_t i = 0; i < size; i += 3){\
			uint8_t temp = data[i];\
			data[i] = data[i + 2];\
			data[i + 2] = temp;\
		}\
	return image(width, height, std::move(data));\
}()
#else
#define read_rgb_image_from_buffer(buf) []() -> image{\
	uint32_t width = *reinterpret_cast<uint32_t*>(buf + 18)\
	int32_t height = *reinterpret_cast<int32_t*>(buf + 22);\
	bool sign = height > 0;\
	uint32_t size = 3 * width * (sign ? height : -height);\
	std::unique_ptr<uint8_t[]> data = std::make_unique<uint8_t[]>(size);\
	bitmap.seekg(*reinterpret_cast<uint32_t*>(buf + 10), std::ios_base::beg);\
	bitmap.read(reinterpret_cast<char*>(data.get()), size);\
	if(sign)\
		for(uint32_t i = 0; i < size; i += 3){\
			uint8_t temp = data[i];\
			data[i] = data[i + 2];\
			data[i + 2] = temp;\
		}\
	return image(width, height, std::move(data));\
}()
*/
image read_rgb_image(const char* path){
	std::ifstream bitmap{path, std::ios::in | std::ios::binary};
	if(!bitmap) return image(0, 0, std::unique_ptr<uint8_t[]>(nullptr));
	char buf[26];
	bitmap.read(buf, 26);
	uint32_t width = *reinterpret_cast<uint32_t*>(buf + 18);
	int32_t height = *reinterpret_cast<int32_t*>(buf + 22);
	bool sign = height > 0;
	uint32_t size = 3 * width * (sign ? height : -height);
	std::unique_ptr<uint8_t[]> data = std::make_unique<uint8_t[]>(size);
	bitmap.seekg(*reinterpret_cast<uint32_t*>(buf + 10), std::ios_base::beg);
	bitmap.read(reinterpret_cast<char*>(data.get()), size);
	if(sign)
		for(uint32_t i = 0; i < size; i += 3){
			uint8_t temp = data[i];
			data[i] = data[i + 2];
			data[i + 2] = temp;
		}
	return image(width, height, std::move(data));
}
#endif
#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
std::atomic<int> framerate(0);
std::atomic<double> yaw_angle_radians(0);
std::atomic<double> pitch_angle_radians(0);
std::atomic<double> roll_angle_radians(0);
void fps(void*){
    while(1){
        Sleep(1000);
        std::cout << framerate << '\n';
        framerate = 0;
    }
}
int width_px = 1200, height_px = 600;//, tick_count;
uint32_t *framebuf;
double width = 1, height = height_px / (double)width_px,
       pixel_inc = width / width_px;
#ifndef horizontal_fov
	#define horizontal_fov 90
#endif
double /*horizontal_fov = 90,*/ vertical_fov = atan(tan(height_px / 2.0) * height_px / (double)width_px);
inline double square(double x){
    return x * x;
}
inline double greater(double a, double b){
    return a > b ? a : b;
}
inline double abs_val(double a){
    return a < 0 ? -a : a;
}
struct point{
    double x, y, z;
    point() : x(0), y(0), z(0){}
    point(double e, double f, double g) : x(e), y(f), z(g){}
	inline point operator-() const{
		return point(-x, -y, -z);
	}
    inline point operator-(point f) const{
        return point(x - f.x, y - f.y, z - f.z);
    }
    inline point& operator-=(point g){
        x -= g.x;
        y -= g.y;
        z -= g.z;
        return *this;
    }
	inline point& operator+=(point o){
		x += o.x;
		y += o.y;
		z += o.z;
		return *this;
	}
    inline point operator*(double h) const{
        return point(x * h, y * h, z * h);
    }
	inline point operator/(double i) const{
		return point(x / i, y / i, z / i);
	}
	void rotate_yaw(double y){
		if(y){
			double cos_yar = cos(y), sin_yar = sin(y), sin_yar_x = sin_yar * this->x;
			this->x = cos_yar * this->x + sin_yar * this->z;
			this->z = -sin_yar_x + cos_yar * this->z;
		}
	}
	void rotate_pitch(double p){
		if(p){
			//todo: check if pitch is actually correct (make it so positive always migrates to top instead of going in circles)
			double cos_par = cos(p), sin_par = sin(p), sin_par_y = sin_par * this->y;
			this->y = cos_par * this->y - sin_par * this->z;
			this->z = sin_par_y + cos_par * this->z;
		}
	}
	void rotate_roll(double r){
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
typedef uint32_t color;
typedef point vector;
inline vector cross_product(vector a, vector b){
    return vector(a.y * b.z - b.y * a.z, a.z * b.x - b.z * a.x, a.x * b.y - b.x * a.y);
}
inline double dot_product(vector a, vector b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
typedef std::vector<point> plane;
typedef std::vector<plane> plane_array;
#define no_color 4278190080
struct hit_val{
	bool first;
	double second;
	color third;
	hit_val(bool b, double d, color c = no_color/*use unused byte*/) : first(b), second(d), third(c){}
};
struct object{
    virtual /*std::pair<bool, double>*/ hit_val hit(vector u) const = 0;
    virtual bool is_in_frustum(const plane_array& plane) const = 0;
    virtual void set_color(const color c){clr = c;}
    color clr;
	#ifdef IMAGE
	virtual void set_image(image&& i){img = std::move(i);}
	image img;
	#endif
	virtual void set_name(const char* str){name = str;}
	const char *class_name, *name;
	double yaw_rad, pitch_rad, roll_rad;
	std::function<void(const vector&)> pressed = [](const vector&){};
	virtual void on_hit(const std::function<void(const vector&)>& input_function){pressed = input_function;}
	virtual void set_rotation(const double y, const double p, const double r){yaw_rad = y, pitch_rad = p, roll_rad = r;}
	virtual void set_yaw(const double y){yaw_rad = y;}
	virtual void set_pitch(const double p){pitch_rad = p;}
	virtual void set_roll(const double r){roll_rad = r;}
};
//typedef object* object_pointer; //maybe change it to unique_ptr or shared_ptr idk
typedef object* object_pointer;
double plane_distance(const plane& bound, const point& p){ 
	vector a = bound[1] - bound[0], b = bound[2] - bound[0], cross = cross_product(a, b);
	double cross_mag = sqrt(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z);
	cross.x /= cross_mag, cross.y /= cross_mag, cross.z /= cross_mag;
	return dot_product(cross, p);
}
struct sphere : public object{
    double radius;
    point center;
	//#ifdef IMAGE
	//double image_yaw_rad, image_pitch_rad, image_roll_rad;
	//#endif
    sphere(color cl, double r, point c, const char* n = "the default sphere"
	#ifdef IMAGE
	, image&& i = read_rgb_image(""), double yr = 0, double pr  = 0, double rr = 0
	#endif 
	) : radius(r), center(c)
	//#ifdef IMAGE
	//, image_yaw_rad(yaw_rad), image_pitch_rad(pitch_rad), image_roll_rad(roll_rad)
	//#endif
	//, yaw_rad(y), pitch_rad(p), roll_rad(r)
	{set_color(cl), set_class_name(), set_name(n)//, set_rotation(yr, pr, rr)
	#ifdef IMAGE
	, set_image(std::move(i)), set_rotation(yr, pr, rr)
	#endif
	;
	}
    virtual void set_color(const color c){clr = c;}
	virtual void set_name(const char* str){name = str;}
	#ifdef IMAGE
	virtual void set_image(image&& i){img = std::move(i);}
	#endif
	virtual void set_class_name(){class_name = "sphere";}
    virtual void on_hit(const std::function<void(const vector&)>& input_function){pressed = input_function;}
	virtual void set_rotation(const double y, const double p, const double r){yaw_rad = y, pitch_rad = p, roll_rad = r;}
	virtual void set_yaw(const double y){yaw_rad = y;}
	virtual void set_pitch(const double p){pitch_rad = p;}
	virtual void set_roll(const double r){roll_rad = r;}
	/*std::pair<bool, double>*/ hit_val hit(vector u) const override{
        //todo: do the same thing for u.x and u.y
		if((u.z > 0 && center.z < radius) || (u.z < 0 && center.z > radius)) return /*std::pair<bool, double>*/ hit_val(false, 0);
		/*	
        double magnitude_squared = u.x * u.x + u.y * u.y + u.z * u.z;
        //double dot = dot_product(u.x, u.y, u.z, center.x, center.y, center.z);
		double dot = dot_product(u, -center);
		double determinant = dot * dot - magnitude_squared * (center.x * center.x + center.y * center.y + center.z * center.z - radius * radius);
        // if determinant < 0 you can't sqrt it
        if(determinant < 0) return std::pair<bool, double>(false, 0);
		//std::cout << (-dot + sqrt(determinant)) / magnitude_squared << '\n';
		//todo: provide correct distance
		//std::cout << "-dot: " << -dot << "\n-sqrt determinant: " << -sqrt(determinant) << "\nmagnitude squared: " << magnitude_squared << '\n';
        return std::pair<bool, double>(true, (-dot - sqrt(determinant)) / magnitude_squared);
		*/
		//todo: its the correct difference but normalizing requires a sqrt
		#if 1
		double magnitude = sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
		u.x /= magnitude, u.y /= magnitude, u.z /= magnitude;
		double dot = dot_product(u, -center);
		double determinant = dot * dot - (center.x * center.x + center.y * center.y + center.z * center.z - radius * radius);
		#else
		//todo: doesnt work for non-unit vectors ig
		double a = dot_product(u, u), b = 2 * dot_product(u, -center), c = dot_product(center, center) - radius * radius;
		double determinant = b * b - 4 * a * c;
		#endif
		if(determinant < 0) return /*std::pair<bool, double>*/ hit_val(false, 0);
		double distance = -dot - sqrt(determinant);
		//double distance = (-b - sqrt(determinant)) / (2 * a);
		//std::cout << distance << '\n';
		pressed(u);
		#ifndef IMAGE
		return /*std::pair<bool, double>*/ hit_val(true, distance);
		#else
		if(!img.width || !img.height) return hit_val(true, distance);
		vector hit_spot = (u * distance) - center;
		hit_spot.rotate_yaw(yaw_rad);
		hit_spot.rotate_pitch(pitch_rad);
		hit_spot.rotate_roll(roll_rad);
		#if 0
		if(/*image_*/yaw_rad){
			double cos_yar = cos(/*image_*/yaw_rad), sin_yar = sin(/*image_*/yaw_rad), sin_yar_x = sin_yar * hit_spot.x;
			hit_spot.x = cos_yar * hit_spot.x + sin_yar * hit_spot.z;
			hit_spot.z = -sin_yar_x + cos_yar * hit_spot.z;
		}
		if(/*image_*/pitch_rad){
			//todo: check if pitch is actually correct (make it so positive always migrates to top instead of going in circles)
			double cos_par = cos(/*image_*/pitch_rad), sin_par = sin(/*image_*/pitch_rad), sin_par_y = sin_par * hit_spot.y;
			hit_spot.y = cos_par * hit_spot.y - sin_par * hit_spot.z;
			hit_spot.z = sin_par_y + cos_par * hit_spot.z;
		}
		if(/*image_*/roll_rad){
			//todo: check if roll is actually correct
			double cos_rar = cos(/*image_*/roll_rad), sin_rar = sin(/*image_*/roll_rad), sin_rar_x = sin_rar * hit_spot.x;
			hit_spot.x = cos_rar * hit_spot.x - sin_rar * hit_spot.y;
			hit_spot.y = sin_rar_x + cos_rar * hit_spot.y;
		}
		#endif
		double hit_spot_magnitude = sqrt(hit_spot.x * hit_spot.x + hit_spot.y * hit_spot.y + hit_spot.z * hit_spot.z);
		hit_spot.x /= hit_spot_magnitude, hit_spot.y /= hit_spot_magnitude, hit_spot.z /= hit_spot_magnitude;
		{
			double u = 0.5 + atan2(hit_spot.z, hit_spot.x) / (2 * M_PI);
			double v = 0.5 + asin(hit_spot.y) / M_PI;
			long x = lrint(u * img.width), y = lrint(v * img.height);
			color hit_color = img.color_at(x, y);
			return hit_val(true, distance, hit_color);
		}
		#endif
	}
    bool is_in_frustum(const plane_array& plane) const override{
		#if 0
		//todo: i guess change vector(0, -1, 0) to the correct vector for all pitch_angle_radians too (up and down)
		// fabricate a back side
		vector left_vector = plane[2][0] - plane[2][1], // vector pointing left from top right to top left
			   projection = left_vector * (dot_product(plane[2][1], left_vector) / dot_product(left_vector, left_vector)), // point on back side
			   plane_normal = cross_product(left_vector /* should it be projection? */, vector(0, -1, 0)); // calculate normal facing back
		//double plane_normal_mag = sqrt(plane_normal.x * plane_normal.x + plane_normal.y * plane_normal.y + plane_normal.z * plane_normal.z);
        //plane_normal.x /= plane_normal_mag, plane_normal.y /= plane_normal_mag, plane_normal.z /= plane_normal_mag;
		double distance = dot_product(plane_normal, center);
		//todo: might not be exactly when its on the other side idk
		if(distance > 0){std::cout << "nah fam\n"; return false;}
		#endif
		#if 0
		for(/*int c = 0; c < (int)plane.size(); ++c*/ std::vector<point>& bound : plane){
            vector a = bound/*plane[c]*/[1] - bound/*plane[c]*/[0], b = bound/*plane[c]*/[2] - bound/*plane[c]*/[0], cross = cross_product(a, b);

            double cross_mag = sqrt(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z);
            cross.x /= cross_mag, cross.y /= cross_mag, cross.z /= cross_mag;

            //if its positive its facing the same direction as the cross product
            //if its negative its facing the opposite

			double distance = dot_product(cross, center);
            if(abs_val(distance) < radius){/*std::cout << "tru moo: " << distance << '\n';*/ return true;}
        }
		
        return false;
		#else
    	//todo: fix returning true when sphere isnt in frustum 
        //corner case: sphere with distance less than radius / sqrt(2) or radius * 0.707106871 to (top or bottom) and (left or right) and vise versa ex. near bottom and left is bottom-left
        //middle case: sphere near a side and is bounded by its two adjacent sides ex. sphere near bottom and bounded by left and right
		//assuming plane[0] is left, plane[1] is right, plane[2] is top, plane[3] is bottom
		double left_distance = plane_distance(plane[0], center), right_distance = plane_distance(plane[1], center),
			   top_distance = plane_distance(plane[2], center), bottom_distance = plane_distance(plane[3], center);
		double corner_distance = radius * 0.707;
		//assuming > 0 is pointing in for some reason idk
		#if 1
		return(
			(left_distance > 0 && right_distance > 0 && top_distance > 0 && bottom_distance > 0)
			||
			//todo: combine left and right
			#if 1
			(left_distance < 0 && left_distance > -radius && ((/*middle*/top_distance > 0 && bottom_distance > 0) || (/*corner*/top_distance > -corner_distance || bottom_distance > -corner_distance)))
			||
			(right_distance < 0 && right_distance > -radius && ((/*middle*/top_distance > 0 && bottom_distance > 0) || (/*corner*/top_distance > -corner_distance || bottom_distance > -corner_distance)))
			#else
			(((left_distance < 0 && left_distance > -radius) || (right_distance < 0 && right_distance > -radius)) && ((/*middle*/top_distance > 0 && bottom_distance > 0) || (/*corner*/top_distance > -corner_distance || bottom_distance > -corner_distance)))
			#endif
			||
			//todo: combine top and bottom
			#if 1
			(top_distance < 0 && top_distance > -radius && ((/*middle*/left_distance > 0 && right_distance > 0) || (/*corner*/left_distance > -corner_distance || right_distance > -corner_distance)))
			||
			(bottom_distance < 0 && bottom_distance > -radius && ((/*middle*/left_distance > 0 && right_distance > 0) || (/*corner*/left_distance > -corner_distance || right_distance > -corner_distance)))
			#else
			(((top_distance < 0 && top_distance > -radius) || (bottom_distance < 0 && bottom_distance > -radius)) && ((/*middle*/left_distance > 0 && right_distance > 0) || (/*corner*/left_distance > -corner_distance || right_distance > -corner_distance)))
			#endif
		);
		#else
			bool a = (
				#if 1
				(left_distance > 0 && right_distance > 0 && top_distance > 0 && bottom_distance > 0)
				||
				(right_distance < 0 && right_distance > -radius && ((/*middle*/top_distance > 0 && bottom_distance > 0) || (/*corner*/top_distance > -corner_distance || bottom_distance > -corner_distance)))
				#else // doesnt work here for some reason
				(((left_distance < 0 && left_distance > -radius) || (right_distance < 0 && right_distance > -radius)) && ((/*middle*/top_distance > 0 && bottom_distance > 0) || (/*corner*/top_distance > -corner_distance || bottom_distance > -corner_distance)))
				#endif
				||
				//todo: combine top and bottom
				#if 1
				(top_distance < 0 && top_distance > -radius && ((/*middle*/left_distance > 0 && right_distance > 0) || (/*corner*/left_distance > -corner_distance || right_distance > -corner_distance)))
				||
				(bottom_distance < 0 && bottom_distance > -radius && ((/*middle*/left_distance > 0 && right_distance > 0) || (/*corner*/left_distance > -corner_distance || right_distance > -corner_distance)))
				#else // doesnt work here for some reason
				(((top_distance < 0 && top_distance > -radius) || (bottom_distance < 0 && bottom_distance > -radius)) && ((/*middle*/left_distance > 0 && right_distance > 0) || (/*corner*/left_distance > -corner_distance || right_distance > -corner_distance)))
				#endif
			);
			std::cout << std::boolalpha << a << '\n';
			return a;
		#endif
		#endif
	}
};
struct polygon : public object{
    std::vector<point> points/*{}*/;
	#ifdef IMAGE
    polygon(color cl, const char* n, image&& i, auto... l) try : points{(point(l))...}{
		set_color(cl);
		set_class_name();
		set_name(n);
		set_image(std::move(i));
		if(sizeof...(l) < 3) throw;
	}catch(...){std::cout << "error: number of points to polygon's constructor must be greater than two\n";}
	#endif
	polygon(color cl, const char* n, auto... l) try /*: clr(cl)*/ : points{(point(l))...}{
        set_color(cl);
		set_class_name();
		set_name(n);
		#ifdef IMAGE
		set_image(std::move(read_rgb_image("")));
		#endif
        if(sizeof...(l) < 3) throw;
        //points = {(point(l))...};
        //a = points[1] - points[0], b = points[2] - points[0], c = cross_product(a, b);
        //k = c.x * points[0].x + c.y * points[1].y + c.z * points[2].z;
    }catch(...){std::cout << "error: number of points to polygon's constructor must be greater than two\n";}
	polygon(color cl, auto... l) try /*: clr(cl)*/ : points{(point(l))...}{
        set_color(cl);
		set_class_name();
		set_name("the default polygon");
        #ifdef IMAGE
		set_image(std::move(read_rgb_image("")));
		#endif
		if(sizeof...(l) < 3) throw;
        //points = {(point(l))...};
        //a = points[1] - points[0], b = points[2] - points[0], c = cross_product(a, b);
        //k = c.x * points[0].x + c.y * points[1].y + c.z * points[2].z;
    }catch(...){std::cout << "error: number of points to polygon's constructor must be greater than two\n";}
	virtual void set_color(const color c){clr = c;}
    virtual void set_name(const char* str){name = str;}
	#ifdef IMAGE
	virtual void set_image(image&& i){img = std::move(i);}
	#endif
	virtual void set_class_name(){class_name = "polygon";}
	virtual void on_hit(const std::function<void(const vector&)>& input_function){pressed = input_function;}
	point calculate_centroid() const{
		static point centroid;
		if(points.size() == 3){
			centroid.x = (points[0].x + points[1].x + points[2].x) / 3, 
			centroid.y = (points[0].y + points[1].y + points[2].y) / 3,
			centroid.z = (points[0].z + points[1].z + points[2].z) / 3;
		}else{/*todo: calculate*/}
		return centroid;
	}
	virtual void set_rotation(const double y, const double p, const double r){
		point centroid = calculate_centroid();
		#if 0
		#if __cplusplus >= 201703L
		[[maybe_unused]]
		#endif
		static bool once = [&](){
			if(is_triangle) 
				centroid.x = (points[0].x + points[1].x + points[2].x) / 3, 
				centroid.y = (points[0].y + points[1].y + points[2].y) / 3,
				centroid.z = greatest;
			else{/*todo: calculate*/}
		};
		#if __cplusplus < 201703L
		(void)once; //for no warning
		#endif
		#endif
		yaw_rad = y, pitch_rad = p, roll_rad = r;
		for(point& i : points){
			i -= centroid;
			i.rotate_yaw(y), i.rotate_pitch(p), i.rotate_roll(r);
			i += centroid;
		}
	}
	virtual void set_yaw(const double y){
		point centroid = calculate_centroid();
		yaw_rad = y;
		for(point& i : points){
			i -= centroid;
			i.rotate_yaw(y);
			i += centroid;
		}
	}
	virtual void set_pitch(const double p){
		point centroid = calculate_centroid();
		pitch_rad = p;
		for(point& i : points){
			i -= centroid;
			i.rotate_pitch(p);
			i += centroid;
		}
	}
	virtual void set_roll(const double r){
		point centroid = calculate_centroid();
		roll_rad = r;
		for(point& i : points){
			i -= centroid;
			i.rotate_roll(r);
			i += centroid;
		}
	}
	void stretch(const double& greatest, point& p) const{
		//fix to accomodate for yar par and rar
        double factor = greatest / p.z;
        p.x *= factor;
        p.y *= factor;
        p.z = greatest;
    };
    /*std::pair<bool, double>*/ hit_val hit(vector u) const override{
		
		static double yar_snap = yaw_angle_radians, par_snap = pitch_angle_radians, rar_snap = roll_angle_radians;
		static double y_snap = yaw_rad, p_snap = pitch_rad, r_snap = roll_rad;
		static std::vector<point> vertices;
		static bool is_triangle;
		static int call_num = 0, calls_per_frame = width_px * height_px;
		static vector a, b, c;
		static double k;
		#ifdef IMAGE
		static point centroid;
		#endif
		#if __cplusplus >= 201703L
		[[maybe_unused]]
		#endif
		static bool once = [&](){
            vertices = points;
			is_triangle = (points.size() == 3);
			a = points[1] - points[0], b = points[2] - points[0], c = cross_product(a, b);
        	//k = c.x * points[0].x + c.y * points[1].y + c.z * points[2].z;
			k = c.x * points[0].x + c.y * points[0].y + c.z * points[0].z;
			double greatest = greater(greater(points[0].z, points[1].z), points[2].z);
			for(point& i : vertices) stretch(greatest, i);
			#ifdef IMAGE
			if(is_triangle) 
				centroid.x = (vertices[0].x + vertices[1].x + vertices[2].x) / 3, 
				centroid.y = (vertices[0].y + vertices[1].y + vertices[2].y) / 3,
				centroid.z = greatest;
			else{/*todo: calculate*/}
			#endif
			return true;
		}();
		#if __cplusplus < 201703L
		(void)once; //for no warning
		#endif

		//todo: do the same thing for u.x and u.y
		for(const point& i : points) //maybe change to vertices
            if((u.z < 0 && i.z < 0) || (u.z > 0 && i.z > 0)) goto polygon_start;
       	return /*std::pair<bool, double>*/ hit_val(false, 0);
        polygon_start:
		#if 0
		static double yar_snap = yaw_angle_radians, par_snap = pitch_angle_radians, rar_snap = roll_angle_radians;
       	static double y_snap = yaw_rad, p_snap = pitch_rad, r_snap = roll_rad;
 		static std::vector<point> vertices;
		static bool is_triangle;
		static int call_num = 0, calls_per_frame = width_px * height_px;
		static vector a, b, c;
		static double k;
		#ifdef IMAGE
		static point centroid;
		#endif
		#if __cplusplus >= 201703L
		[[maybe_unused]]
		#endif
		static bool once = [&](){
            vertices = points;
			is_triangle = (points.size() == 3);
			a = points[1] - points[0], b = points[2] - points[0], c = cross_product(a, b);
        	//k = c.x * points[0].x + c.y * points[1].y + c.z * points[2].z;
			k = c.x * points[0].x + c.y * points[0].y + c.z * points[0].z;
			double greatest = greater(greater(points[0].z, points[1].z), points[2].z);
			for(point& i : vertices) stretch(greatest, i);
			#ifdef IMAGE
			if(is_triangle) 
				centroid.x = (vertices[0].x + vertices[1].x + vertices[2].x) / 3, 
				centroid.y = (vertices[0].y + vertices[1].y + vertices[2].y) / 3,
				centroid.z = greatest;
			else{/*calculate*/}
			#endif
			return true;
		}();
		#if __cplusplus < 201703L
		(void)once; //for no warning
		#endif
		#endif
        double e = u.x * c.x + u.y * c.y + u.z * c.z;
        double t = k / e;
        if(likely(!e || (u.z * t < 0))) return /*std::pair<bool, double>*/ hit_val(false, 0);
        point intersection(u.x * t, u.y * t, u.z * t);
		//std::cout << "t: " << t << '\n';
		//std::cout << "intersection: " << intersection.x << ' ' << intersection.y << ' ' << intersection.z << '\n';
        double distance = sqrt(intersection.x * intersection.x + intersection.y * intersection.y + intersection.z * intersection.z);
		double greatest = greater(greater(points[0].z, points[1].z), points[2].z);

		//return [&greatest, &distance, &intersection, this](std::vector<point>& points){
		/*
		static double yar_snap = yaw_angle_radians, par_snap = pitch_angle_radians, rar_snap = roll_angle_radians;
        static std::vector<point> vertices;
		static int call_num = 0, calls_per_frame = width_px * height_px;;
		#if __cplusplus >= 201703L
		[[maybe_unused]]
		#endif
		static bool once = [&](){
            vertices = points;
			a = points[1] - points[0], b = points[2] - points[0], c = cross_product(a, b);
        	k = c.x * points[0].x + c.y * points[1].y + c.z * points[2].z;
			double greatest = greater(greater(points[0].z, points[1].z), points[2].z);
			for(point& i : vertices) stretch(greatest, i);
			return true;
		}();
		#if __cplusplus < 201703L
		(void)once; //for no warning
		#endif
        */
		if(likely(call_num < calls_per_frame)){
        	stretch(greatest, intersection);
			++call_num;
		}else{
            call_num = 0;
			if(vertices.size() != points.size()) goto update;
			for(int i = 0; i < (int)vertices.size(); ++i)
				if(vertices[i].x != points[i].x || vertices[i].y != points[i].y) goto update;
			if(calls_per_frame != width_px * height_px)
                calls_per_frame = width_px * height_px;
			if(y_snap != yaw_rad || p_snap != pitch_rad || r_snap != roll_rad){
				y_snap = yaw_rad, p_snap = pitch_rad, r_snap = roll_rad;
				goto update;
			}
			if(yar_snap != yaw_angle_radians){
                yar_snap = yaw_angle_radians; 
                goto update;
            }
            if(par_snap != pitch_angle_radians){
                par_snap = pitch_angle_radians;
                goto update;
            }
            if(rar_snap != roll_angle_radians){
                rar_snap = roll_angle_radians;
                update:
                vertices = points;
				is_triangle = (points.size() == 3);
				a = points[1] - points[0], b = points[2] - points[0], c = cross_product(a, b);
        		//k = c.x * points[0].x + c.y * points[1].y + c.z * points[2].z;
                k = c.x * points[0].x + c.y * points[0].y + c.z * points[0].z;
				double greatest = greater(greater(points[0].z, points[1].z), points[2].z);
    			for(point& i : vertices) stretch(greatest, i);
            	#ifdef IMAGE
				if(is_triangle) 
					centroid.x = (vertices[0].x + vertices[1].x + vertices[2].x) / 3, 
					centroid.y = (vertices[0].y + vertices[1].y + vertices[2].y) / 3,
					centroid.z = greatest;
				else{/*calculate*/}
				#endif
				stretch(greatest, intersection);
            }
			/*
			if(y_snap != yaw_rad || p_snap != pitch_rad || r_snap != roll_rad){
				y_snap = yaw_rad, p_snap = pitch_rad, r_snap = roll_rad;
				for(point& i : points)
					i.rotate_yaw(yaw_rad), i.rotate_pitch(pitch_rad), i.rotate_roll(roll_rad);
			}
			*/
        }
        //#define tri_specialization
        //#define first_algorithm

        #ifdef tri_specialization
        if(/*vertices.size() == 3*/ is_triangle){
            #ifdef first_algorithm
				int one = 1, two = 2;
                if(vertices[2].y == vertices[0].y){
                    //point temp = vertices[2];
                    //vertices[2] = vertices[1];
                    //vertices[1] = temp;
					one = 2;
					two = 1;
                }
                double s1 = vertices[two].y - vertices[0].y,
                       s2 = vertices[two].x - vertices[0].x,
                       s3 = vertices[one].y - vertices[0].y,
                       s4 = intersection.y - vertices[0].y;
                double w1 = (vertices[0].x * s1 + s4 * s2 - intersection.x * s1) / (s3 * s2 - (vertices[one].x - vertices[0].x) * s1),
                       w2 = (s4 - w1 * s3) / s1;
				#ifndef IMAGE
				bool hit_polygon = w1 >= 0 && w2 >= 0 && (w1 + w2) <= 1;
				if(hit_polygon) pressed(u);
				return hit_val(hit_polygon, distance);
                //return /*std::pair<bool, double>*/ hit_val(w1 >= 0 && w2 >= 0 && (w1 + w2) <= 1, distance);
            	#else
				//do
				#endif
			#else
                auto sign = [](point p1, point p2, point p3) -> double{
                    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
                };
                double d1 = sign(intersection, vertices[0], vertices[1]),
                       d2 = sign(intersection, vertices[1], vertices[2]),
                       d3 = sign(intersection, vertices[2], vertices[0]);
                bool neg = d1 < 0 || d2 < 0 || d3 < 0,
                     pos = d1 > 0 || d2 > 0 || d3 > 0,
                     ret = !(neg && pos);
				if(ret) pressed(u);
				#ifndef IMAGE
				return /*std::pair<bool, double>*/ hit_val(ret, distance);
            	#else
				//do
				#endif
			#endif
        }
        #endif
        bool inside = false;
        for(int i = 0, j = vertices.size() - 1; i < (int)vertices.size(); j = i++)
            if(((vertices[i].y > intersection.y) != (vertices[j].y > intersection.y)) &&
                (intersection.x < (vertices[j].x - vertices[i].x) * (intersection.y - vertices[i].y) / (vertices[j].y - vertices[i].y) + vertices[i].x))
                    inside = !inside;
		if(inside) pressed(u);
		#ifndef IMAGE
		return /*std::pair<bool, double>*/ hit_val(inside, distance);
		#else
		if(!is_triangle) return hit_val(inside, distance);
		if(!img.width || !img.height) return hit_val(inside, distance);
		//int x = img.width / 2 + intersection.x - centroid.x, y = img.height / 2 + intersection.y - centroid.y;
	 	double triangle_width = greater(greater(abs_val(vertices[1].x - vertices[0].x), abs_val(vertices[2].x - vertices[0].x)), abs_val(vertices[2].x - vertices[1].x)), triangle_height = greater(greater(abs_val(vertices[1].y - vertices[0].y), abs_val(vertices[2].y - vertices[0].y)), abs_val(vertices[2].y - vertices[1].y)); 
		double x = 0.5 + (intersection.x - centroid.x) / /*700 triangle width*/triangle_width, y = 0.5 + (intersection.y - centroid.y) / /*triangle height*/triangle_height;
		//color hit_color = img.color_at(x, y);
		if(x < 0 || x > 1 || y < 0 || y > 1) return hit_val(inside, distance, RGB(255, 255, 255));
		color hit_color = img.color_at(x * img.width, y * img.height);
		return hit_val(inside, distance, hit_color);
		#endif
		//}(const_cast<std::vector<point>&>(points));
	}
    bool is_in_frustum(const plane_array& plane) const override{
		for(const point& i : points){
			/*
			vector left_vector = plane[2][0] - plane[2][1], // vector pointing left from top right to top left
				   projection = left_vector * (dot_product(plane[2][1], left_vector) / dot_product(left_vector, left_vector)), // point on back side
				   plane_normal = cross_product(left_vector, vector(0, -1, 0)); // calculate normal facing back
			double distance = dot_product(plane_normal, i);
			if(distance > 0){this_point = false; break;}
			*/
			/*
			for(std::vector<point>& bound : plane){
				vector a = bound[1] - bound[0], b = bound[2] - bound[0], cross = cross_product(a, b);
				//double cross_mag = sqrt(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z);
            	//cross.x /= cross_mag, cross.y /= cross_mag, cross.z /= cross_mag;
				double distance = dot_product(cross, i);
				//assuming < 0 is pointing out for some reason
				if(distance < 0){
					this_point = false; 
					break;
				}
			}
			if(this_point){std::cout << "polygon in frustum real\n"; return true;}
			*/
			//assuming again
			double left_distance = plane_distance(plane[0], i), right_distance = plane_distance(plane[1], i),
				   top_distance = plane_distance(plane[2], i), bottom_distance = plane_distance(plane[3], i);
        	//std::cout << "point " << i.x << ' ' << i.y << ' ' << i.z << 
			//			 " distances: " << left_distance << ' ' << right_distance << ' ' << top_distance << ' ' << bottom_distance << '\n';
			if(left_distance > 0 && right_distance > 0 && top_distance > 0 && bottom_distance > 0) return true;
		}
        //return false;
		//todo: moving left twice and down once makes the polygon disappear :(
		//todo: when going down the distance to bottom becomes negative for some reason
        //return this->hit(plane[0][0]).first || this->hit(plane[0][1]).first || this->hit(plane[1][0]).first || this->hit(plane[1][1]).first; 
    	bool corner_in_polygon = this->hit(plane[0][0]).first || this->hit(plane[0][1]).first || this->hit(plane[1][0]).first || this->hit(plane[1][1]).first; 
		//if(!corner_in_polygon) std::cout << "nooooooooooooo\n";
		return corner_in_polygon;
	}
//private:
    //vector a, b, c;
    //double k;
};

struct doughnut : public object{
    point center;
    double minor_radius, major_radius;
    //double yaw_rad, pitch_rad;
    doughnut(color cl, double minor_r, double major_r, point c, double yaw = 0, double pitch = 0) :
        center(c), minor_radius(minor_r), major_radius(major_r), /*yaw_rad(yaw), pitch_rad(pitch),*/ epsilon(major_r * major_r - minor_r * minor_r){set_color(cl), set_class_name(), set_name("the default doughnut"), set_yaw(yaw), set_pitch(pitch);}
    doughnut(color cl, double minor_r, double major_r, point c, const char* n, double yaw = 0, double pitch = 0) :
        center(c), minor_radius(minor_r), major_radius(major_r), /*yaw_rad(yaw), pitch_rad(pitch),*/ epsilon(major_r * major_r - minor_r * minor_r){set_color(cl), set_class_name(), set_name(n); set_yaw(yaw), set_pitch(pitch);}
	virtual void set_color(const color c){clr = c;}
    virtual void set_name(const char* str){name = str;}
	virtual void set_class_name(){class_name = "doughnut";}
	virtual void on_hit(const std::function<void(const vector&)>& input_function){pressed = input_function;}
	virtual void set_rotation(const double y, const double p, const double r){yaw_rad = y, pitch_rad = p, roll_rad = r;}
	virtual void set_yaw(const double y){yaw_rad = y;}
	virtual void set_pitch(const double p){pitch_rad = p;}
	virtual void set_roll(const double r){roll_rad = r;}
	/*std::pair<bool, double>*/ hit_val hit(vector u) const override{
        //u -= center; // todo: with a center other than the origin the donut does weird things
		//already accounted for in formula but it still does weird things
		#if 0
		if(yaw_rad){
			double cos_yar = cos(-yaw_rad), sin_yar = sin(-yaw_rad), sin_yar_x = sin_yar * u.x;
			u.x = cos_yar * u.x + sin_yar * u.z;
			u.z = -sin_yar_x + cos_yar * u.z;
		}
		if(pitch_rad){
			//todo: check if pitch is actually correct (make it so positive always migrates to top instead of going in circles)
			double cos_par = cos(-pitch_rad), sin_par = sin(-pitch_rad), sin_par_y = sin_par * u.y;
			u.y = cos_par * u.y - sin_par * u.z;
			u.z = sin_par_y + cos_par * u.z;
		}
		#endif
		u.rotate_yaw(yaw_rad);
		u.rotate_pitch(pitch_rad);
        double magnitude = sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
        constexpr double scaler = 1;
        u.x /= (magnitude * scaler), u.y /= (magnitude * scaler), u.z /= (magnitude * scaler);
        long double a, b, c, d, e;
        {
            long double _a = -center.x, _b = u.x, _c = -center.y, _d = u.y, _e = -center.z, _f = u.z;

            #define simplify_donut

            #ifdef simplify_donut
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
        long double discriminant = 256 * e * e * e - 192 * b * d * e * e - 128 * c * c * e * e + 144 * c * d * d * e - 27 * d * d * d *d + 144 * b * b * c * e * e - 6 * b * b * d * d * e - 80 * b * c * c * d * e + 18 * b * c * d * d * d + 16 * c * c * c * c * e - 4 * c * c * c * d * d - 27 * b * b * b * b * e * e + 18 * b * b * b * c * d * e - 4 * b * b * b * d * d * d - 4 * b * b * c * c * c * e + b * b * c * c * d * d;
        long double p = 8 * c - 3 * b * b;
        long double g = 64 * e - 16 * c * c + 16 * b * b * c - 16 * b * d - 3 * b * b * b * b;
        //todo: find the distance instead of just putting 100
        if(discriminant < 0){pressed(u); return /*std::pair<bool, double>*/ hit_val(true, 100);}
        if(!discriminant){
			bool hit_doughnut = !(!g && p > 0 && b * b * b + 8 * d - 4 * b * c);
			if(hit_doughnut) pressed(u);
			return hit_val(hit_doughnut, 100);
			//return /*std::pair<bool, double>*/ hit_val(!(!g && p > 0 && b * b * b + 8 * d - 4 * b * c), 100);
		}
		//if(discriminant > 0)
        //return std::pair<bool, double>(!(p > 0 || g > 0), 100);
		bool hit_doughnut = p < 0 && g < 0;
       	if(hit_doughnut) pressed(u);
		return hit_val(hit_doughnut, 100);
		//return /*std::pair<bool, double>*/ hit_val(p < 0 && g < 0, 100);
    }
    bool is_in_frustum(const plane_array& plane) const override{
		(void)plane; //for no warning
        return true; //for now
    }
private:
    double epsilon;
};
#if 0
#ifndef IMAGE
sphere s(RGB(0, 0, 255), sqrt(100000), point(100, 41.5, 2000), "the big red sphere");
#else
sphere s(no_color, sqrt(100000), point(0, 50, 2000), "the big blue earth", read_rgb_image("images/earth.bmp"));
#endif
//doughnut d(RGB(0, 255, 0), 200, 300, point(100, 100, 5000));
doughnut d(RGB(255, 0, 255), 90, 1500, point(0, 0, 0), "the laggy doughnut", 1.57, 1);
#ifndef IMAGE
point p1(-250, -300, 2400), p2(250, 200, 1500), p3(350, -100, 1000)
#ifdef QUAD
    , p4(-250, -100, 2000)
#endif
;
#else
point p1(-350, -165, 2400), p2(350, -165, 2400), p3(0, /*-1000*/-1400, 2400)
#ifdef QUAD
	, p4(0, 1000, 24000)
#endif
;
#endif
polygon p(
#ifndef IMAGE
RGB(0, 255, 0)
#else
no_color
#endif
,
#ifdef QUAD
"the funny quadrilateral"
#else
"the funny triangle"
#endif
#ifdef IMAGE
, read_rgb_image("images/cone.bmp")
#endif
, p1, p2, p3
#ifdef QUAD
    , p4
#endif
);
std::vector<object_pointer> world = {&s, &p
#ifdef NO_DOUGHNUT
	};
#else
, &d};
#endif
#endif
std::vector<std::unique_ptr<object>> world;
constexpr double convert_horizontal_fov_to_radians = //maybe just make it a bool and ternary the 0.0175 in calculating for z
#ifdef in_radians
//1 / 0.0175 //to cancel out the 0.0175
1
#else
//1
0.0175
#endif
;
double z = width / (2 * tan(horizontal_fov * convert_horizontal_fov_to_radians /* * 0.0175*/ / 2));
typedef std::function<void(std::vector<object_pointer>&, plane_array&, bool&)> plane_setup_type;
template<typename T> struct lambda_destructor{
public:
    T &lambda, after;
    lambda_destructor(T& l, T a) : lambda(l), after(a){}
    ~lambda_destructor(){lambda = after;}
};
plane_setup_type plane_setup;
std::function<void()> render =
[](){
	++framerate;
    //todo: use the gpu for calculations
	std::vector<object_pointer> can_hit;
	
	static bool yaw_rotated = 0;
	static double yaw_total = 0;
	constexpr int object_num_lower_bound = 0, object_num_upper_bound = 1; //inclusive not exclusive
	if(!yaw_rotated){
		yaw_total += 0.0523599;
		for(int object_num = object_num_lower_bound; object_num <= object_num_upper_bound; ++object_num)
			world[object_num]->set_yaw(yaw_total);
		yaw_rotated = abs(fmod(yaw_total, 2 * M_PI)) <= 0.01;
	}else{
		static bool pitch_rotated = 0;
		static double pitch_total = 0;
		if(!pitch_rotated){
			pitch_total += 0.0523599;
			for(int object_num = object_num_lower_bound; object_num <= object_num_upper_bound; ++object_num)
				world[object_num]->set_pitch(pitch_total);
			pitch_rotated = abs(fmod(pitch_total, 2 * M_PI)) <= 0.01;
		}
		else yaw_rotated = 0, yaw_total = 0, pitch_rotated = 0, pitch_total = 0;
	}
	
	plane_array plane;
    bool ypr;
    //for clicks and stuff
    plane_setup_type create_plane = [](std::vector<object_pointer>& can_hit, plane_array& plane, bool& ypr){
        //can_hit = world; // copy
		//this is where segfault
		if(can_hit.size() != world.size()) can_hit.resize(world.size());
		for(int i = 0; i < (int)world.size(); ++i)
			can_hit[i] = &*world[i];
        point bounding[4]{
            point(-width / 2.0, height / 2.0, z), //upper left
            point(width / 2.0, height / 2.0, z), //upper right
            point(-width / 2.0, -height / 2.0, z), //bottom left
            point(width / 2.0, -height / 2.0, z) //bottom right
        };
        ypr = /*yaw_angle_radians*/ abs_val(fmod(yaw_angle_radians, 2 * M_PI)) > 0.01 ||
    			   /*pitch_angle_radians*/ abs_val(fmod(pitch_angle_radians, 2 * M_PI)) > 0.01 ||
    			   /*roll_angle_radians*/ abs_val(fmod(roll_angle_radians, 2 * M_PI)) > 0.01;
    	if(ypr){
    		for(point& i : bounding){
				/*
    			if(yaw_angle_radians){
    				double cos_yar = cos(yaw_angle_radians), sin_yar = sin(yaw_angle_radians), sin_yar_x = sin_yar * i.x;
    				i.x = cos_yar * i.x + sin_yar * i.z;
    				i.z = -sin_yar_x + cos_yar * i.z;
    			}
    			if(pitch_angle_radians){
    				//todo: check if pitch is actually correct (make it so positive always migrates to top instead of going in circles)
    				double cos_par = cos(pitch_angle_radians), sin_par = sin(pitch_angle_radians), sin_par_y = sin_par * i.y;
    				i.y = cos_par * i.y - sin_par * i.z;
    				i.z = sin_par_y + cos_par * i.z;
    			}
    			if(roll_angle_radians){
    				//todo: check if roll is actually correct
    				double cos_rar = cos(roll_angle_radians), sin_rar = sin(roll_angle_radians), sin_rar_x = sin_rar * i.x;
    				i.x = cos_rar * i.x - sin_rar * i.y;
    				i.y = sin_rar_x + cos_rar * i.y;
    			}
				*/
				i.rotate_yaw(yaw_angle_radians);
				i.rotate_pitch(pitch_angle_radians);
				i.rotate_roll(roll_angle_radians);
    		}
    	}
    	//todo: i guess change vector(0, -1, 0) to the correct vector for all pitch_angle_radians too (up and down)
    	// fabricate a back side
		//todo: make the backside the computer screen
    	vector left_vector = bounding[0] - bounding[1], // vector pointing left from top right to top left
    		   projection = left_vector * (dot_product(bounding[0], left_vector) / dot_product(left_vector, left_vector)); // point on back side
    
        plane = {
            {bounding[0], bounding[2], point(0, 0, 0)}, //left
    		{bounding[3], bounding[1], point(0, 0, 0)}, //right
    		{bounding[1], bounding[0], point(0, 0, 0)}, //top
            {bounding[2], bounding[3], point(0, 0, 0)}, //bottom
    		{point(0, 0, 0), projection, point(0, -1, 0) /*fix*/}, //back test	
    		//todo: maybe add front plane
        };
    	#if __cplusplus >= 202002L
            std::erase_if(can_hit, [&plane](object_pointer i){return !i->is_in_frustum(plane);});
        #else
            can_hit.erase(std::remove_if(can_hit.begin(), can_hit.end(), [&plane](object_pointer i){return !i->is_in_frustum(plane);}), can_hit.end());
        #endif
	};
    create_plane(can_hit, plane, ypr);
    plane_setup = [&](std::vector<object_pointer>& c, plane_array& p, bool& yp){c = can_hit, p = plane, yp = ypr;};
    lambda_destructor<plane_setup_type>(plane_setup, create_plane);
	if(!can_hit.size()){
		for(int i = 0; i < width_px * height_px; ++i) framebuf[i] = RGB(255, 255, 255);
		return;
	}
	if(ypr){
		vector start = plane[0][0]; //start at top left
		vector row_inc = (plane[2][0] - plane[2][1]) / width_px, column_dec = (plane[0][0] - plane[0][1]) / height_px;
		for(int jpx = height_px - 1; jpx >= 0; start -= column_dec, --jpx){
			double row_vx = start.x, row_vy = start.y, row_vz = start.z;
			int prev = jpx * width_px;
			for(int ipx = 0; ipx < width_px; row_vx += row_inc.x, row_vy += row_inc.y, row_vz += row_inc.z, ++ipx){
				int index = prev + ipx;
				vector v(row_vx, -row_vy, row_vz);
				object* smallest = can_hit[0];
				/*std::pair<bool, double>*/ hit_val hit_b = smallest->hit(v);
				bool hit_nothing = !hit_b.first;//, small_change = false;
				//todo: likely and unlikely might be unnecessary
				for(int w = 1; w < (int)can_hit.size(); ++w){
					//if(comp(hit_nothing, v, world[i], smallest)) smallest = world[i];
					/*std::pair<bool, double>*/ hit_val hit_a = can_hit[w]->hit(v);				
					//todo: why small change?
					//if(small_change) hit_b = smallest->hit(v), small_change = false;

					if(unlikely(hit_a.first)){
						hit_nothing = false;
						if(unlikely(hit_b.first)){
							if(hit_a.second < hit_b.second)
								smallest = can_hit[w], hit_b = smallest->hit(v);//small_change = true;
						}
						else smallest = can_hit[w], hit_b = smallest->hit(v);//small_change = true;
					}
					else hit_nothing = !hit_b.first;
				}
				if(likely(hit_nothing)) framebuf[index] = RGB(255, 255, 255);//RGB(255, 0, 0);
				else 
				#ifndef IMAGE
					framebuf[index] = smallest->clr;
				#else
					framebuf[index] = (hit_b.third != no_color ? hit_b.third : smallest->clr);
				#endif
			}
		}
	}else{
		int ipx, jpx = height_px - 1;
		for(double j = height/*- 1*/; j >= 0 && jpx >= 0; /*--j*/j -= pixel_inc, --jpx){
			ipx = 0;
			for(double i = 0; i < width /* 1 */ && ipx < width_px; /*++i*/i += pixel_inc, ++ipx){
				int index = jpx * width_px + ipx;
				vector v(-width / 2 + i, height / 2 - j, z);

				object* smallest = can_hit[0];
				/*std::pair<bool, double>*/ hit_val hit_b = smallest->hit(v);
				bool hit_nothing = !hit_b.first;//, small_change = false;
				//todo: likely and unlikely might be unnecessary
				for(int w = 1; w < (int)can_hit.size(); ++w){
					//if(comp(hit_nothing, v, world[i], smallest)) smallest = world[i];
					/*std::pair<bool, double>*/ hit_val hit_a = can_hit[w]->hit(v);
					//if(small_change) hit_b = smallest->hit(v), small_change = false;

					if(unlikely(hit_a.first)){
						hit_nothing = false;
						if(unlikely(hit_b.first)){
							if(hit_a.second < hit_b.second)
								smallest = can_hit[w], hit_b = smallest->hit(v);//small_change = true;
						}
						else smallest = can_hit[w], hit_b = smallest->hit(v);//small_change = true;
					}
					else hit_nothing = !hit_b.first;
				}
				if(likely(hit_nothing)) framebuf[index] = RGB(255, 255, 255);//RGB(255, 0, 0);
				else 
				#ifndef IMAGE
					framebuf[index] = smallest->clr;
				#else
					framebuf[index] = (hit_b.third != no_color ? hit_b.third : smallest->clr);
				#endif

			}
		}
	}
}
, resize =
[](){
    z = width / (2 * tan(horizontal_fov * convert_horizontal_fov_to_radians /* * 0.0175*/ / 2));
    vertical_fov = atan(tan(height_px / 2.0) * height_px / (double)width_px);
    height = height_px / (double)width_px;
    pixel_inc = width / width_px;
};
#ifdef SOUND
std::vector<std::unique_ptr<char[]>> sounds;
void play_sounds(void*){
	for(int i = 0; i < (int)sounds.size(); ++i)
		PlaySound(TEXT(sounds[i].get()), NULL, SND_FILENAME | SND_SYNC);
	sounds.clear();
	return;
}
#endif
LRESULT CALLBACK WindowProcessMessages(HWND hwnd, UINT msg, WPARAM w, LPARAM l){
    static HDC pdc;
    static HBITMAP old;
    static HBITMAP bitmap;
    switch(msg){
        case WM_CREATE:{
            SetTimer(hwnd, 1, 1, 0);
            HDC hdc;
            BITMAPINFO bitmapinfo{};
            hdc = CreateCompatibleDC(0);
            bitmapinfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
            bitmapinfo.bmiHeader.biWidth = width_px;
            bitmapinfo.bmiHeader.biHeight = -height_px; // top down is negative
            bitmapinfo.bmiHeader.biPlanes = 1;
            bitmapinfo.bmiHeader.biBitCount = 32;
            bitmapinfo.bmiHeader.biCompression = BI_RGB;
            bitmapinfo.bmiHeader.biClrUsed = 256;
            bitmapinfo.bmiHeader.biClrImportant = 256;
            bitmap = CreateDIBSection(hdc, &bitmapinfo, DIB_RGB_COLORS, (void**)&framebuf, 0, 0);
            pdc = CreateCompatibleDC(0);
            old = (HBITMAP)SelectObject(pdc, bitmap);
            DeleteDC(hdc);
            break;
        }
        case WM_SIZE:{
            width_px = LOWORD(l);
            height_px = HIWORD(l);
            resize();
            SelectObject(pdc, old);
            DeleteDC(pdc);
            DeleteObject(bitmap);
            HDC hdc;
            BITMAPINFO bitmapinfo{};
            hdc = CreateCompatibleDC(0);
            bitmapinfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
            bitmapinfo.bmiHeader.biWidth = width_px;
            bitmapinfo.bmiHeader.biHeight = -height_px; // top down is negative
            bitmapinfo.bmiHeader.biPlanes = 1;
            bitmapinfo.bmiHeader.biBitCount = 32;
            bitmapinfo.bmiHeader.biCompression = BI_RGB;
            bitmapinfo.bmiHeader.biClrUsed = 256;
            bitmapinfo.bmiHeader.biClrImportant = 256;
            bitmap = CreateDIBSection(hdc, &bitmapinfo, DIB_RGB_COLORS, (void**)&framebuf, 0, 0);
            pdc = CreateCompatibleDC(0);
            old = (HBITMAP)SelectObject(pdc, bitmap);
            DeleteDC(hdc);
            break;
        }
        case WM_KEYDOWN:
			switch(w){
				case 'W': case VK_UP:
					pitch_angle_radians = pitch_angle_radians + 0.15;
					break;
				case 'S': case VK_DOWN:
					pitch_angle_radians = pitch_angle_radians - 0.15;
					break;
				case 'A': case VK_LEFT:
					yaw_angle_radians = yaw_angle_radians - 0.15;
					break;
				case 'D': case VK_RIGHT:
					yaw_angle_radians = yaw_angle_radians + 0.15;
					break;
			}
        	break;
		case WM_LBUTTONDOWN:{
			int x = GET_X_LPARAM(l), y = GET_Y_LPARAM(l);
            std::vector<object_pointer> can_hit;
            plane_array plane;
            bool ypr;
            plane_setup(can_hit, plane, ypr);
			vector v;
			if(ypr){
				vector start = plane[0][0]; //start at top left
				vector row_inc = (plane[2][0] - plane[2][1]) / width_px, column_dec = (plane[0][0] - plane[0][1]) / height_px;
				start -= column_dec * (height_px - y);
				start.x += row_inc.x * x;
				start.y += row_inc.y * x;
				start.z += row_inc.z * x;
				v = vector(start.x, -start.y, start.z);
			}
			else
				v = vector(-width / 2 + x * pixel_inc, height / 2 - y * pixel_inc, z);
			object* smallest = can_hit[0];
			/*std::pair<bool, double>*/ hit_val hit_b = smallest->hit(v);
			bool hit_nothing = !hit_b.first;//, small_change = false;
			//todo: likely and unlikely might be unnecessary
			for(int w = 1; w < (int)can_hit.size(); ++w){
				//if(comp(hit_nothing, v, world[i], smallest)) smallest = world[i];
				/*std::pair<bool, double>*/ hit_val hit_a = can_hit[w]->hit(v);				
				//if(small_change) hit_b = smallest->hit(v), small_change = false;
				
				if(unlikely(hit_a.first)){
					hit_nothing = false;
					if(unlikely(hit_b.first)){
						if(hit_a.second < hit_b.second)
							smallest = can_hit[w], hit_b = smallest->hit(v);//small_change = true;
					}
					else smallest = can_hit[w], hit_b = smallest->hit(v);//small_change = true;
				}
				else hit_nothing = !hit_b.first;
			}
			if(likely(hit_nothing)) std::cout << "you clicked the vast emptiness of space, devoid of any shred of liveliness and hope...\n";
			else{
				#ifdef SOUND
				//static char sound_name[256];
				//strncpy(sound_name, (std::string("sound/") + std::string(smallest->class_name) + std::string(".wav")).c_str(), 256);
				bool playing = sounds.size();
				std::string path = std::string("sound/") + std::string(smallest->class_name);
				for(const auto& file : fs::directory_iterator(path)){
					//strncpy(sound_name, file.path().string().c_str(), 255);
					//sounds.emplace_back(std::unique_ptr<char[]>(new char[256]));
					sounds.emplace_back(std::make_unique<char[]>(256));
					strncpy(sounds.back().get(), file.path().string().c_str(), 255);
					//PlaySound(TEXT(sound_name), NULL, SND_FILENAME | SND_ASYNC);
				}
				if(!playing) _beginthread(play_sounds, 0, 0);
				#endif
				std::cout << "you clicked on a " << smallest->class_name << " and its name is \'" << smallest->name << "\'!\n";
			}
			break;
		}
		case WM_TIMER:
            InvalidateRgn(hwnd, 0, 0);
            UpdateWindow(hwnd);
            break;
        case WM_PAINT:{
            //tick_count = GetTickCount();
            PAINTSTRUCT ps;
            HDC h = BeginPaint(hwnd, &ps);
            render();
            BitBlt(h, 0, 0, width_px, height_px, pdc, 0, 0, SRCCOPY);
            EndPaint(hwnd, &ps);
            break;
        }
        case WM_DESTROY:
            SelectObject(pdc, old);
            DeleteDC(pdc);
            DeleteObject(bitmap);
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
	#ifdef IMAGE
	std::unique_ptr<sphere> s = std::make_unique<sphere>(no_color, sqrt(100000), point(0, 50, 2000), "the big blue earth", read_rgb_image("images/earth.bmp"));
	#else
	std::unique_ptr<sphere> s = std::make_unique<sphere>(RGB(0, 0, 255), sqrt(100000), point(0, 50, 2000), "the red sphere");
	#endif
	//s->on_hit([](const vector& v){
	//	std::cout << "hi\n";
	//});
	#ifndef NO_DOUGHNUT
	std::unique_ptr<doughnut> d = std::make_unique<doughnut>(RGB(255, 0, 255), 90, 1500, point(0, 0, 0), "the laggy doughnut", 1.57, 1);
	#endif
	//point p1(-250, -300, 2400), p2(250, 200, 1500), p3(350, -100, 1000);
	point p1(-350, -165, 2400), p2(350, -165, 2400), p3(0, /*-1000*/-1400, 2400);
	#ifdef IMAGE
	std::unique_ptr<polygon> p = std::make_unique<polygon>(no_color, "the funny triangle", read_rgb_image("images/cone.bmp"), p1, p2, p3);
	#else
	std::unique_ptr<polygon> p = std::make_unique<polygon>(RGB(0, 255, 0), "the funny triangle", p1, p2, p3);
	#endif
	world.emplace_back(std::move(s));
	#ifndef NO_DOUGHNUT
	world.emplace_back(std::move(d));
	#endif
	world.emplace_back(std::move(p));
	_beginthread(fps, 0, 0);
	char name[] = "omg";
	WNDCLASS wc{};
	wc.lpszClassName = name;
	wc.hCursor = LoadCursor(0, IDC_ARROW);
	wc.lpfnWndProc = WindowProcessMessages;
	RegisterClass(&wc);
	int a = width_px + 16, b = height_px + 39;
	CreateWindow(name, name, WS_OVERLAPPEDWINDOW | WS_CAPTION | WS_SYSMENU | WS_VISIBLE, CW_USEDEFAULT, CW_USEDEFAULT, a, b, 0, 0, 0, 0);
	MSG msg{};
	while(GetMessage(&msg, 0, 0, 0)){
	    TranslateMessage(&msg);
	    DispatchMessage(&msg);
	}
}
