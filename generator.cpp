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
                vertex_back << stod(vertex_array_split[1]), stod(vertex_array_split[2]), stod(vertex_array_split[3]);
                vertex_points.emplace_back(vertex_back);

            }
        }
    }
    read_TESS.close();

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

    string face_vertex_text;
    string face_equation_parameter_text;
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
                Plane_Para(vertex_points[face_vertex_array(i_fep, 0)].transpose(),
                           vertex_points[face_vertex_array(i_fep, 1)].transpose(),
                           vertex_points[face_vertex_array(i_fep, 2)].transpose()));
    }

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
                for (int number_polyhedron_face = 2;
                     number_polyhedron_face < polyhedron_faces_number; ++number_polyhedron_face) {
                    int face_list_element = abs(stoi(polyhedron_text_split[number_polyhedron_face]));
                    face_list.push_back(face_list_element);
                }
                polyhedron_faces.push_back(face_list);
            }
        }
    }
    read_TESS.close();
    //读取最大半径并且读取多面体顶点组合
    string stpoly_file_string;
    ifstream read_STPOLY;
    read_STPOLY.open("n100-id1.stpoly", ios::in);
    MatrixXd radius_MAX_array(1, seed_number);
    vector<vector<int >> polyhedron_vertex;
    for (int number_line = 0; number_line < seed_number; ++number_line) {
        getline(read_STPOLY, stpoly_file_string);
        vector<string> stpoly_file_split = split(stpoly_file_string, " ");
        radius_MAX_array(0, number_line) = stod(stpoly_file_split[0]);
        int vertex_num = (int) stpoly_file_split.size();
        vector<int> vertex_push_back;
        for (int i_vertex = 1; i_vertex < vertex_num; ++i_vertex) {
            vertex_push_back.emplace_back(stoi(stpoly_file_split[i_vertex]) - 1); //减1是为了与读取的顶点容器编号相匹配
        }
        polyhedron_vertex.emplace_back(vertex_push_back);
    }
    read_STPOLY.close();
    //######读取四元数##############
    string quaternion_text;
    ifstream read_ORI;
    read_ORI.open("n100-id1.ori", ios::in);
    MatrixXd quaternion_array(seed_number, 4);
    for (int number_ori = 0; number_ori < seed_number; ++number_ori) {
        getline(read_ORI, quaternion_text);
        vector<string> quaternion_test_split = split(quaternion_text, " ");
        quaternion_array.row(number_ori) << stod(quaternion_test_split[0]), stod(quaternion_test_split[1]), stod(
                quaternion_test_split[2]), stod(quaternion_test_split[3]);
    }
    read_ORI.close();
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
    //Matrix<double, 1, 3> cell_cen_test;
    //cell_cen_test << 57, 56, 460;
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
    //#########设置晶格常数
    double lattice_parameter = 0.8;
    //#############初始化实际原子位置容器
    vector<vector<double >> atoms_position;
    //开始添加原子
    int cout_number = 1;
    for (int number_poly = 18; number_poly < 19; ++number_poly) {

        //#########初始化此镶嵌所对应的面
        vector<int> poly_chosen = polyhedron_faces[number_poly];
        //########选取此晶粒的等效半径
        double radius = radius_MAX_array(0, number_poly);
        //#########计算所需立方体晶粒的分割次数
        int cubic_size = 2 * radius / lattice_parameter;
        //##########计算旋转矩阵
        Matrix3d rotation_Matrix;
        Quaterniond quaternion(quaternion_array(number_poly, 0),
                               quaternion_array(number_poly, 1),
                               quaternion_array(number_poly, 2),
                               quaternion_array(number_poly, 3));
        rotation_Matrix = Matrix3d(quaternion);
        vector<MatrixXd> cubic_atoms;
        for (int number_x = -cubic_size; number_x < cubic_size + 1; ++number_x) {
            for (int number_y = -cubic_size; number_y < cubic_size + 1; ++number_y) {
                for (int number_z = -cubic_size; number_z < cubic_size + 1; ++number_z) {
                    //#######循环到所在的分格点
                    MatrixXd base_position(1, 3);
                    base_position << number_x, number_y, number_z;
                    //#######换算成预生成晶格所在坐标
                    MatrixXd cart_position(1, 3);
                    cart_position = (cell_martix * base_position.transpose()).transpose();
                    //##########生成实际预生成原子位置
                    for (int number_lattice_array = 0; number_lattice_array < 4; ++number_lattice_array) {
                        MatrixXd atom_position(1, 3);
                        atom_position =
                                (cart_position + lattice_constant_array.row(number_lattice_array)) * lattice_parameter;
                        //################真实原子坐标
                        cubic_atoms.emplace_back((rotation_Matrix * atom_position.transpose()).transpose() +
                                                 cell_centroid.row(number_poly));
                    }
                }
            }
        }
        //预生成原子总个数
        int pre_atoms_tot = cubic_atoms.size();
        //cout << pre_atoms_tot << endl;
        //所选多面体面总数
        int poly_face_tot = polyhedron_faces[number_poly].size();
        //所选镶嵌顶点数
        int poly_vertex_tot = polyhedron_vertex[number_poly].size();
        //构建镶嵌顶点坐标组
        MatrixXd poly_chosen_vertex(polyhedron_vertex[number_poly].size(), 3);
        for (int number_poly_vertex = 0; number_poly_vertex < poly_vertex_tot; ++number_poly_vertex) {
            poly_chosen_vertex.row(number_poly_vertex)
                    << vertex_points[polyhedron_vertex[number_poly][number_poly_vertex]].transpose();
        }
        //选取多面体极值
        double min_x_poly = poly_chosen_vertex.col(0).minCoeff();
        double max_x_poly = poly_chosen_vertex.col(0).maxCoeff();
        double min_y_poly = poly_chosen_vertex.col(1).minCoeff();
        double max_y_poly = poly_chosen_vertex.col(1).maxCoeff();
        double min_z_poly = poly_chosen_vertex.col(2).minCoeff();
        double max_z_poly = poly_chosen_vertex.col(2).maxCoeff();
        //cout << "min_x_poly = " << min_x_poly << "\n";
        //cout << "max_x_poly = " << max_x_poly << "\n";
        //cout << "min_y_poly = " << min_y_poly << "\n";
        //cout << "max_y_poly = " << max_y_poly << "\n";
        //cout << "min_z_poly = " << min_z_poly << "\n";
        //cout << "max_z_poly = " << max_z_poly << endl;
        //构建镶嵌区域内原子
        for (int number_atom = 0; number_atom < pre_atoms_tot; ++number_atom) {
            MatrixXd atom;
            atom = cubic_atoms[number_atom];
            //cout<<atom<<endl;
            if (min_x_poly < atom(0, 0) && atom(0, 0) < max_x_poly &&
                min_y_poly < atom(0, 1) && atom(0, 1) < max_y_poly &&
                min_z_poly < atom(0, 2) && atom(0, 2) < max_z_poly) {
                int cross_number = 0;
                for (int number_face = 0; number_face < poly_face_tot; ++number_face) {
                    if(Check_In_Plane(face_equation_parameters[poly_chosen[number_face]-1],atom)){
                        vector<double> back ={(double)number_poly,atom(0,0),atom(0,1),atom(0,2)};
                        atoms_position.emplace_back(back);
                    }
                }
            }
        }

    }

    ofstream outdata;
    outdata.precision(6);
    outdata.open("test.dat", ios::out);
    outdata << "Crystalline Cu atoms\n\n";
    outdata << atoms_position.size() << " atoms\n";
    outdata << seed_number << " atom types\n";
    outdata << fixed;
    outdata << "0" << " " << "100" << " xlo xhi\n";
    outdata << "0" << " " << "100" << " ylo yhi\n";
    outdata << "0" << " " << "100" << " zlo zhi\n";
    outdata << "\n";
    outdata << "Atoms\n\n";
    for (int number_atom = 0; number_atom < atoms_position.size(); ++number_atom) {
        outdata << number_atom + 1 << " " << (int) atoms_position[number_atom][0];
        for (int number_out = 1; number_out < 4; ++number_out) {
            outdata << " " << atoms_position[number_atom][number_out];
        }
        outdata << endl;
    }
    outdata.close();

    return 0;
}

