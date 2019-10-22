//
// Created by a on 2019/10/20.
//

#include <iostream>

using namespace std;

//三维doublex向量
struct Vec3d {
    double x, y, z;

    Vec3d() {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }

    Vec3d(double dx, double dy, double dz) {
        x = dx;
        y = dy;
        z = dz;
    }

    void Set(double dx, double dy, double dz) {
        x = dx;
        y = dy;
        z = dz;
    }
};

//计算三点成面的法向量
bool Sign_of_Point(const Vec3d &p1, const Vec3d &p2, const Vec3d &p3, const Vec3d &point) {
    //p1(n1,n2,n3)
    //平面方程 na * (x – n1) + nb * (y – n2) + nc * (z – n3) = 0 ;
    double a = (p2.y - p1.y) * (p3.z - p1.z) - (p2.z - p1.z) * (p3.y - p1.y);
    double b = (p2.z - p1.z) * (p3.x - p1.x) - (p2.x - p1.x) * (p3.z - p1.z);
    double c = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
    double val_point = a*(point.x-p1.x)+b*(point.y-p1.y)+c*(point.z-p1.z);
    return val_point > 0;
}
int main(){
    Vec3d v1(1.0,0.0,0.0);
    Vec3d v2(0.0,1.0,0.0);
    Vec3d v3(1.0,1.0,1.0);
    Vec3d v4(1.0,1.0,10.0);
    Vec3d v5(0.0,0.0,0.0);
    if(Sign_of_Point(v1,v2,v3,v4)==Sign_of_Point(v1,v2,v3,v5)){
        cout<<"on the same side"<<endl;
    } else{
        cout<<"not ont the same side"<<endl;
    }
    //cout<<Sign_of_Point(v1,v2,v3,v4)<<endl;
    return 0;
}