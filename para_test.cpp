//
// Created by a on 2019/10/24.
//

#include <iostream>
#include <eigen3/Eigen/SVD>
#include <string>
#include <fstream>
#include <vector>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost;
using namespace Eigen;

vector<string> split(const string &str, const string &pattern);

MatrixXd
Plane_Para(const Ref<const MatrixXd> &point_1, const Ref<const MatrixXd> &point_2, const Ref<const MatrixXd> &point_3);

int main() {
    int face_number = 594;
    int vertex_number = 448;
    string skip_line_text;
    string line_text;

    ifstream read_TESS;
    string vertex_array;
    MatrixXd vertex_cartesian(vertex_number, 3);
    read_TESS.open("n100-id1.tess", ios::in);
    while (!read_TESS.eof()) {
        getline(read_TESS, line_text);
        if (line_text == " **vertex") {
            getline(read_TESS, skip_line_text);
            for (int number_vertex = 0; number_vertex < vertex_number; number_vertex++) {
                getline(read_TESS, vertex_array);
                vector<string> vertex_array_split = split(vertex_array, " ");
                vertex_cartesian(number_vertex, 0) = stod(vertex_array_split[1]);
                vertex_cartesian(number_vertex, 1) = stod(vertex_array_split[2]);
                vertex_cartesian(number_vertex, 2) = stod(vertex_array_split[3]);
            }
        }
    }
    read_TESS.close();

    string face_vertex_text;
    MatrixXi face_vertex_array(face_number, 3);
    read_TESS.open("n100-id1.tess", ios::in);
    while (!read_TESS.eof()) {
        getline(read_TESS, line_text);
        if (line_text == " **face") {
            getline(read_TESS, skip_line_text);
            for (int number_face = 0; number_face < face_number; ++number_face) {
                getline(read_TESS, face_vertex_text);
                getline(read_TESS, skip_line_text);
                getline(read_TESS, skip_line_text);
                getline(read_TESS, skip_line_text);
                vector<string> face_vertex_text_split = split(face_vertex_text, " ");
                face_vertex_array(number_face, 0) = stoi(face_vertex_text_split[2]) - 1;
                face_vertex_array(number_face, 1) = stoi(face_vertex_text_split[3]) - 1;
                face_vertex_array(number_face, 2) = stoi(face_vertex_text_split[4]) - 1;
            }
        }
    }
    read_TESS.close();
    vector<Matrix<double, 1, 4>> face_equation_parameters;
    for (int i_fep = 0; i_fep < face_number; ++i_fep) {
        face_equation_parameters.emplace_back(
                Plane_Para(vertex_cartesian.row(face_vertex_array(i_fep, 0)),
                           vertex_cartesian.row(face_vertex_array(i_fep, 1)),
                           vertex_cartesian.row(face_vertex_array(i_fep, 2))));
        //cout << i_fep << " ";
        //cout << face_equation_parameters(i_fep) << endl;
    }
    cout << endl;
    for (int j = 0; j < 594; ++j) {
        cout << j << " " << face_equation_parameters[j] << endl;
    }

    return 0;
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
    MatrixXd para(1, 4);
    para(0, 1) =
            (point_2(0, 1) - point_1(0, 1)) * (point_3(0, 2) - point_1(0, 2)) -
            (point_3(0, 1) - point_1(0, 1)) * (point_2(0, 2) - point_1(0, 2));
    para(0, 2) =
            (point_2(0, 2) - point_1(0, 2)) * (point_3(0, 0) - point_1(0, 0)) -
            (point_3(0, 2) - point_1(0, 2)) * (point_2(0, 0) - point_1(0, 0));
    para(0, 3) =
            (point_2(0, 0) - point_1(0, 0)) * (point_3(0, 1) - point_1(0, 1)) -
            (point_2(0, 1) - point_1(0, 1)) * (point_3(0, 0) - point_1(0, 0));
    para(0, 0) =
            para(0, 0) * point_1(0, 0) +
            para(0, 1) * point_1(0, 1) +
            para(0, 2) * point_1(0, 2);
    return para;
}