//
// Created by a on 2019/10/25.
//

#ifndef POLYCRY_GENE_GENERATOR_ALGO_H
#define POLYCRY_GENE_GENERATOR_ALGO_H


#include <string>
#include <cmath>
#include <vector>
#include <eigen3/Eigen/SVD>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost;
using namespace Eigen;

const double EPSILON = 0.0000001;

bool Check_In_Plane(const Ref<MatrixXd> &para, const Ref<MatrixXd> &point);

vector<string> split(const string &str, const string &pattern);

MatrixXd
Plane_Para(const Ref<const MatrixXd> &point_1, const Ref<const MatrixXd> &point_2, const Ref<const MatrixXd> &point_3);

int Is_Point_In_Poly(const Ref<const MatrixXd> &para, const Ref<const MatrixXd> &point);


#endif //POLYCRY_GENE_GENERATOR_ALGO_H
