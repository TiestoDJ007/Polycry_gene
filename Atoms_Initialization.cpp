//
// Created by a on 2019/10/16.
//
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;
//初始化晶体常数矩阵
//初始化中心点
//镶嵌最大距离
const double max_radius = 12.8;

//两点距离原型
double Points_Distance(vector<double> &p1, vector<double> &p2);

double Inner_Cross(vector<double> &p1, vector<int> &p2);

//template<class T>
//int get_array_len(T &array) { return sizeof(array) / sizeof(array[0]); }

int main() {
    //初始化晶体常数矩阵
    vector<vector<double >> lattice_constant_array = {{0.0, 0.0, 0.0},
                                                      {0.5, 0.5, 0.0},
                                                      {0.5, 0.0, 0.5},
                                                      {0.0, 0.5, 0.5}};
    //初始化中心点
    vector<double> centroid_cartesian = {1.0,
                                         2.0,
                                         3.0};
    int lattice_array_size = lattice_constant_array.size();
    vector<vector<double >> cell_array = {{1, 0, 0},
                                          {0, 1, 0},
                                          {0, 0, 1}};
    //计算能包含住球体的正方体一个顶点
    double cubic_vertex_x = centroid_cartesian[0] - max_radius;
    double cubic_vertex_y = centroid_cartesian[1] - max_radius;
    double cubic_vertex_z = centroid_cartesian[2] - max_radius;
    vector<double> cubic_vertex = {cubic_vertex_x, cubic_vertex_y, cubic_vertex_z};
    //计算正方体中的原子排列,再由正方体生成球形晶粒
    int cubic_size_element = (int) max_radius * 2 / 1;
    vector<int> cubic_size = {cubic_size_element, cubic_size_element, cubic_size_element};
    vector<vector<double>> cubic_atoms;
    for (int number_x = 0; number_x < cubic_size[0]; ++number_x) {
        for (int number_y = 0; number_y < cubic_size[1]; ++number_y) {
            for (int number_z = 0; number_z < cubic_size[2]; ++number_z) {
                vector<int> base_position = {number_x, number_y, number_z};
                vector<double> cart_position = {Inner_Cross(cell_array[0], base_position),
                                                Inner_Cross(cell_array[1], base_position),
                                                Inner_Cross(cell_array[2], base_position)};
                for (int number_lattice_array = 0; number_lattice_array < lattice_array_size; ++number_lattice_array) {
                    vector<double> atom_position = {
                            cart_position[0] + lattice_constant_array[number_lattice_array][0] + cubic_vertex_x,
                            cart_position[1] + lattice_constant_array[number_lattice_array][1] + cubic_vertex_y,
                            cart_position[2] + lattice_constant_array[number_lattice_array][2] + cubic_vertex_z};
                    double compare_val = Points_Distance(atom_position, centroid_cartesian) - max_radius;
                    if (compare_val > 0) {
                        cubic_atoms.push_back(atom_position);
                    }
                }
            }
        }
    }
    ofstream outfile;
    outfile.open("structure.dat");
    outfile<<"Crystalline Cu atoms"<<'\n'<<'\n';
    outfile<<(int)cubic_atoms.size()<<" atoms"<<'\n';
    outfile<<"1"<<" atom types"<<'\n';
    outfile.close();
    //
    return 0;
}

//两点之间距离函数
double Points_Distance(vector<double> &p1, vector<double> &p2) {
    double Distance = sqrt(pow(p2[0] - p1[0], 2) + pow(p2[1] - p1[1], 2) + pow(p2[2] - p1[2], 2));
    return Distance;
}

double Inner_Cross(vector<double> &p1, vector<int> &p2) {
    return (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]);
}

