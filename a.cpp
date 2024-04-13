#include <iostream>
#include <cmath>
using namespace std;
int main(){
    constexpr double fov = 90;
    constexpr int x = -5, y = -10, z = 50, width = 1366, height = 768;
    constexpr double a = (x / (z * tan(fov * 3.14159265 / 360) * 2) + 0.5) * width,
                     b = (y / (z * tan(fov * 3.14159265 / 360) * 2) + 0.5) * height;
    cout << -5 / (50 * tan(90 * 3.14159265 / 360)) * height + height / 2 << '\n';
    cout << a << ' ' << b << '\n';
}