//
// Created by a on 2019/10/25.
//

#include "generator_algo.h"

bool Check_In_Plane(const Ref<MatrixXd> &para, const Ref<MatrixXd> &point) {
    double result = para(0, 1) * point(0, 0) +
                    para(0, 2) * point(0, 1) +
                    para(0, 3) * point(0, 2) -
                    para(0, 0);
    return abs(result) <= EPSILON;
}


vector<string> split(const string &str, const string &pattern) {
    vector<string> result;
    int position_start = str.find_first_not_of(pattern);
    string new_str = str.substr(position_start, str.size());
    split(result, new_str, is_any_of(pattern), token_compress_on);
    return result;
}

MatrixXd
Plane_Para(const Ref<const MatrixXd> &point_1, const Ref<const MatrixXd> &point_2, const Ref<const MatrixXd> &point_3) {

    double a =
            (point_2(0, 1) - point_1(0, 1)) * (point_3(0, 2) - point_1(0, 2)) -
            (point_3(0, 1) - point_1(0, 1)) * (point_2(0, 2) - point_1(0, 2));
    double b =
            (point_2(0, 2) - point_1(0, 2)) * (point_3(0, 0) - point_1(0, 0)) -
            (point_3(0, 2) - point_1(0, 2)) * (point_2(0, 0) - point_1(0, 0));
    double c =
            (point_2(0, 0) - point_1(0, 0)) * (point_3(0, 1) - point_1(0, 1)) -
            (point_2(0, 1) - point_1(0, 1)) * (point_3(0, 0) - point_1(0, 0));
    double length = sqrt(a * a + b * b + c * c);

    double d = a / length * point_1(0, 0) +
               b / length * point_1(0, 1) +
               c / length * point_1(0, 2);
    MatrixXd para(1, 4);
    para << d, a / length, b / length, c / length;
    return para;
}

int Is_Point_In_Poly(const Ref<const MatrixXd> &para, const Ref<const MatrixXd> &point) {
    double distance = para(0, 1) * point(0, 0) +
                      para(0, 2) * point(0, 1) +
                      para(0, 3) * point(0, 2) -
                      para(0, 0);
    if (distance >= 0)
        return 1;
    else
        return -1;
}


