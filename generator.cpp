//
// Created by a on 2019/10/20.
//
//polyhedron_faces无法使用eigen矩阵，因为长度相同。其他数据存储方式皆为eigen矩阵。
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <cmath>
//#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>

using namespace std;
using namespace boost;
using namespace Eigen;

int seed_number;
int vertex_number;
int edge_number;
int face_number;
int polyhedron_number;


vector<string> split(const string &str, const string &pattern);


int Is_Point_In_Poly(const Ref<const MatrixXd> para, const Ref<const MatrixXd> point);


int main() {
    //读取.tess文件
    ifstream read_TESS;
    read_TESS.open("n100-id1.tess", ios::in);
    string line_text;
    string cell_text;
    string vertex_text;
    string edge_text;
    string face_text;
    string polyhedron_text;
    string skip_line_text;
    int test = 0;
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
        test = test + 1;
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
                //cell_centroid(number_seed,0) = centroid_x;
                //{centroid_x, centroid_y, centroid_z};
            }
        }
    }
    //cout << cell_centroid << endl;
    read_TESS.close();
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
    //cout << vertex_cartesian[1][1] << endl;


    string edge_array_text;
    MatrixXi edge_array(edge_number, 2);
    read_TESS.open("n100-id1.tess", ios::in);
    while (!read_TESS.eof()) {
        getline(read_TESS, line_text);
        if (line_text == " **edge") {
            getline(read_TESS, skip_line_text);
            for (int number_edge = 0; number_edge < edge_number; number_edge++) {
                getline(read_TESS, edge_array_text);
                vector<string> edge_text_split = split(edge_array_text, " ");
                edge_array(number_edge, 0) = stoi(edge_text_split[1]);
                edge_array(number_edge, 1) = stoi(edge_text_split[2]);
            }
        }
    }
    read_TESS.close();
    //cout << edge_array[1][1] << endl;


    string face_vertex_text;
    string face_equation_parameter_text;
    MatrixXi face_vertex_array(face_number, 3);
    MatrixXd face_equation_parameters(face_number, 4);
    read_TESS.open("n100-id1.tess", ios::in);
    while (!read_TESS.eof()) {
        getline(read_TESS, line_text);
        if (line_text == " **face") {
            getline(read_TESS, skip_line_text);
            for (int number_face = 0; number_face < face_number; ++number_face) {
                getline(read_TESS, face_vertex_text);
                getline(read_TESS, skip_line_text);
                getline(read_TESS, face_equation_parameter_text);
                getline(read_TESS, skip_line_text);
                vector<string> face_vertex_text_split = split(face_vertex_text, " ");
                face_vertex_array(number_face, 0) = stoi(face_vertex_text_split[1]);
                face_vertex_array(number_face, 1) = stoi(face_vertex_text_split[2]);
                face_vertex_array(number_face, 2) = stoi(face_vertex_text_split[3]);
                vector<string> face_eq_para_split = split(face_equation_parameter_text, " ");
                face_equation_parameters(number_face, 0) = stod(face_eq_para_split[0]);
                face_equation_parameters(number_face, 1) = stod(face_eq_para_split[1]);
                face_equation_parameters(number_face, 2) = stod(face_eq_para_split[2]);
                face_equation_parameters(number_face, 3) = stod(face_eq_para_split[3]);
            }
        }
    }
    read_TESS.close();
    string polyhedron_face_text;
    vector<vector<int >> polyhedron_faces;
    read_TESS.open("n100-id1.tess", ios::in);
    while (!read_TESS.eof()) {
        getline(read_TESS, line_text);
        if (line_text == " **polyhedron") {
            getline(read_TESS, skip_line_text);
            for (int number_polyhedron = 0; number_polyhedron < polyhedron_number; ++number_polyhedron) {
                getline(read_TESS, polyhedron_face_text);
                vector<string> polyhedron_text_split = split(polyhedron_face_text, " ");
                vector<int> face_list;
                int polyhedron_faces_number = polyhedron_text_split.size();
                for (int number_polyhedron_face = 1;
                     number_polyhedron_face < polyhedron_faces_number; ++number_polyhedron_face) {
                    int face_list_element = abs(stoi(polyhedron_text_split[number_polyhedron_face]));
                    face_list.push_back(face_list_element);
                }
                polyhedron_faces.push_back(face_list);
            }
        }
    }
    read_TESS.close();
    //cout << polyhedron_faces[1][1] << endl;
    //读取最大半径
    string radius_text;
    ifstream read_STPOLY;
    read_STPOLY.open("n100-id1.stpoly", ios::in);
    MatrixXd radius_MAX_array(1, seed_number);
    for (int number_radius = 0; number_radius < seed_number; ++number_radius) {
        getline(read_STPOLY, radius_text);
        radius_MAX_array(0, number_radius) = stod(radius_text);
    }
    cout.precision(15);
    //cout<<fixed<<radius_MAX_array<<endl;
    //读取文件完毕


    //cell_centroid包含各个镶嵌的中心点
    //vertex_cartesian包含各个镶嵌的顶点
    //edge_array包含每条边所对应的两点
    //face_vertex_array包含每个面所对应的顶点
    //face_equation_parameters包含每个面所对应的方程参数
    //polyhedron_faces包含镶嵌所对应的面
    //##################################################
    //所有数组皆为c++排列方式，即从0-N，而不是从1-N
    //需要在重新组合时，将对应数值加1
    //##################################################

    //fcc单位晶格
    MatrixXd lattice_constant_array(4, 3);
    lattice_constant_array << 0.0, 0.0, 0.0,
            0.5, 0.5, 0.0,
            0.5, 0.0, 0.5,
            0.0, 0.5, 0.5;
    //晶格取向
    MatrixXd cell_martix(3, 3);
    cell_martix << 1, 0, 0,
            0, 1, 0,
            0, 0, 1;


    //对区域 “0” 进行实验
    //选取镶嵌0的区域
    int choosing_polyhedron = 50;
    double lattice_parameter = 3.6149;
    vector<int> poly_chosen = polyhedron_faces[choosing_polyhedron];
    double radius = radius_MAX_array(0,choosing_polyhedron);

    MatrixXd centroid_points(1, 3);
    centroid_points << cell_centroid(choosing_polyhedron, 0),
            cell_centroid(choosing_polyhedron, 1),
            cell_centroid(choosing_polyhedron, 2);
    //cout << lattice_constant_array.row(1) << endl;
    MatrixXd cubic_vertex_array(1, 3);
    cubic_vertex_array << centroid_points(0) - radius, centroid_points(1) - radius, centroid_points(2) - radius;
    //MatrixXd cubic_size(1,3);
    //cubic_size<<(int)radius*2/1,(int)radius*2/1,(int)radius*2/1;
    int cubic_size = (int) radius * 2 / lattice_parameter;
    int cubic_atoms_number = 0;
    MatrixXd cubic_atoms(cubic_size * cubic_size * cubic_size * 4, 3);
    for (int number_x = 0; number_x < cubic_size; ++number_x) {
        for (int number_y = 0; number_y < cubic_size; ++number_y) {
            for (int number_z = 0; number_z < cubic_size; ++number_z) {
                MatrixXd base_position(1, 3);
                base_position << number_x, number_y, number_z;
                Vector3d cart_position = cell_martix * base_position.transpose();
                //cout<<cart_position.transpose()<<endl;
                for (int number_lattice_array = 0; number_lattice_array < 4; ++number_lattice_array) {
                    MatrixXd atom_position(1, 3);
                    atom_position = cart_position.transpose() + lattice_constant_array.row(number_lattice_array) +
                                    cubic_vertex_array;
                    cubic_atoms(cubic_atoms_number, 0) = atom_position(0, 0)*lattice_parameter;
                    cubic_atoms(cubic_atoms_number, 1) = atom_position(0, 1)*lattice_parameter;
                    cubic_atoms(cubic_atoms_number, 2) = atom_position(0, 2)*lattice_parameter;
                    cubic_atoms_number = cubic_atoms_number + 1;
                    //cout << atom_position << endl;
                }
            }
        }
    }
    //构建镶嵌区域内原子

    //double a;
    //a = Is_Point_In_Poly(face_equation_parameters.row(0), cubic_atoms.row(0));
    //cout << a << endl;
    //对区域一进行实验
    int face_tot_number = polyhedron_faces[choosing_polyhedron].size();
    cout<<face_tot_number<<endl;
    MatrixXi centroid_signal(1, face_tot_number);
    for (int number_face = 0; number_face < face_tot_number; ++number_face) {
        centroid_signal(0, number_face) = Is_Point_In_Poly(
                face_equation_parameters.row(polyhedron_faces[choosing_polyhedron][number_face] - 1),
                cell_centroid.row(choosing_polyhedron));
    }
    cout << centroid_signal << endl;
    //使用判断原子和中心点在所有面的同一方向
    //MatrixXd print_atom;
    //Matrix<float, Dynamic, 3> xxx;
    //print_atom.row(0) <<0,0,0;
    int pre_atoms_tot_number = cubic_atoms.rows();
    int atoms_number = 0;
    vector<int> atoms_number_array;
    for (int number_pre_atom = 0; number_pre_atom < pre_atoms_tot_number; ++number_pre_atom) {
        MatrixXi atom_signal(1, face_tot_number);
        for (int number_face = 0; number_face < face_tot_number; ++number_face) {
            atom_signal(0, number_face) = Is_Point_In_Poly(
                    face_equation_parameters.row(polyhedron_faces[choosing_polyhedron][number_face] - 1),
                    cubic_atoms.row(number_pre_atom));
        }
        if (atom_signal == centroid_signal) {
            atoms_number_array.push_back(number_pre_atom);
            atoms_number++;
        }
    }
    cout << pre_atoms_tot_number << endl;
    cout << atoms_number << endl;
    cout<<cubic_vertex_array(0,0)<<endl;
    ofstream outdata;
    outdata.precision(6);
    outdata.open("test.dat", ios::out);
    outdata << "Crystalline Cu atoms\n\n";
    outdata << atoms_number << " atoms\n";
    outdata << "1 atom types\n";
    outdata << fixed<<cubic_vertex_array(0, 0) << " " << cubic_vertex_array(0, 0) + radius * 2 << " xlo xhi\n";
    outdata << cubic_vertex_array(0, 1) << " " << cubic_vertex_array(0, 1) + radius * 2 << " ylo yhi\n";
    outdata << cubic_vertex_array(0, 2) << " " << cubic_vertex_array(0, 2) + radius * 2 << " zlo zhi\n";
    outdata << "\n";
    outdata << "Atoms\n\n";
    for (int number_atom = 0; number_atom < atoms_number; ++number_atom) {
        outdata << number_atom + 1 << " 1 " << cubic_atoms.row(number_atom)<<"\n";
    }
    outdata.close();
    return 0;


}

vector<string> split(const string &str, const string &pattern) {
    vector<string> result;
    int position_start = str.find_first_not_of(pattern);
    string new_str = str.substr(position_start, str.size());
    split(result, new_str, is_any_of(pattern), token_compress_on);
    return result;
}

int Is_Point_In_Poly(const Ref<const MatrixXd> para, const Ref<const MatrixXd> point) {
    double distance = para(0, 0) * point(0, 0) +
                      para(0, 1) * point(0, 1) +
                      para(0, 2) * point(0, 2) -
                      para(0, 3);
    if (distance > 0)
        return 1;
    else
        return -1;
}
