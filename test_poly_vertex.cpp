//
// Created by a on 2019/10/20.
//
//polyhedron_faces无法使用eigen矩阵，因为长度相同。其他数据存储方式皆为eigen矩阵。
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Geometry>
#include "generator_algo.h"


int seed_number;
int vertex_number;
int edge_number;
int face_number;
int polyhedron_number;


int main() {
    //读取.tess文件
    ifstream read_TESS;
    string line_text;
    string cell_text;
    string vertex_text;
    string edge_text;
    string face_text;
    string polyhedron_text;
    string skip_line_text;
    read_TESS.open("n100-id1.tess", ios::in);
    while (getline(read_TESS, line_text)) {
        if (line_text == " **cell") {
            getline(read_TESS, cell_text);
            seed_number = stoi(cell_text, nullptr, 10);
        }
        if (line_text == " **vertex") {
            getline(read_TESS, vertex_text);
            vertex_number = stoi(vertex_text, nullptr, 10);
        }
        if (line_text == " **edge") {
            getline(read_TESS, edge_text);
            edge_number = stoi(edge_text, nullptr, 10);
        }
        if (line_text == " **face") {
            getline(read_TESS, face_text);
            face_number = stoi(face_text, nullptr, 10);
        }
        if (line_text == " **polyhedron") {
            getline(read_TESS, polyhedron_text);
            polyhedron_number = stoi(polyhedron_text, nullptr, 10);
        }
    }
    read_TESS.close();

    string cell_array;
    MatrixXd cell_centroid(seed_number, 3);
    read_TESS.open("n100-id1.tess", ios::in);
    while (!read_TESS.eof()) {
        getline(read_TESS, line_text);
        if (line_text == "  *seed") {
            for (int number_seed = 0; number_seed < seed_number; number_seed++) {
                getline(read_TESS, cell_array);
                vector<string> cell_array_split = split(cell_array, " ");
                cell_centroid(number_seed, 0) = stod(cell_array_split[1]);
                cell_centroid(number_seed, 1) = stod(cell_array_split[2]);
                cell_centroid(number_seed, 2) = stod(cell_array_split[3]);
            }
        }
    }

    read_TESS.close();
    string vertex_array;
    vector<Vector3d> vertex_points;
    read_TESS.open("n100-id1.tess", ios::in);
    while (!read_TESS.eof()) {
        getline(read_TESS, line_text);
        if (line_text == " **vertex") {
            getline(read_TESS, skip_line_text);
            for (int number_vertex = 0; number_vertex < vertex_number; number_vertex++) {
                getline(read_TESS, vertex_array);
                vector<string> vertex_array_split = split(vertex_array, " ");
                Vector3d vertex_back;
                vertex_back << stod(vertex_array_split[1]),stod(vertex_array_split[2]),stod(vertex_array_split[3]);
                vertex_points.emplace_back(vertex_back);

            }
        }
    }
    read_TESS.close();
    for (int j = 0; j < 100; ++j) {
        cout<<vertex_points[j].transpose()<<" "<<endl;
    }
}
