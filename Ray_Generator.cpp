//
// Created by a on 2019/10/27.
//
#include <iostream>
#include <eigen3/Eigen/Geometry>
#include "Ray_File_Func.h"

using namespace std;

inline Eigen::Matrix3d rotation_Matrix(const Eigen::Ref<Eigen::Vector4d> &quater_ori);

inline double
Distance_Point_Face(const Eigen::Ref<Eigen::Vector4d> &face_para, const Eigen::Ref<Eigen::Vector3d> &point);

int main() {
    Ray_File_Func test("n100-id1.tess", "n100-id1.stpoly", "n100-id1.ori");
    cout << rotation_Matrix(test.poly_ori()[1]) << endl;
    return 0;
}

inline Eigen::Matrix3d rotation_Matrix(const Eigen::Ref<Eigen::Vector4d> &quater_ori) {
    return Eigen::Matrix3d(Eigen::Quaterniond(quater_ori(0, 0),
                                              quater_ori(0, 0),
                                              quater_ori(0, 0),
                                              quater_ori(0, 0)));
}

inline double
Distance_Point_Face(const Eigen::Ref<Eigen::Vector4d> &face_para, const Eigen::Ref<Eigen::Vector3d> &point) {
    return (face_para(0, 0) * point(0, 0) + face_para(0, 1) * point(0, 1) + face_para(0, 2) * point(0, 2) +
            face_para(0, 3)) / sqrt(face_para(0, 0) * face_para(0, 0) + face_para(0, 1) * face_para(0, 1) +
                                    face_para(0, 2) * face_para(0, 2));
}