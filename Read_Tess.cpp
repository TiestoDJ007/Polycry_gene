//
// Created by a on 2019/10/14.
//
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <cmath>

using namespace std;
using namespace boost;

int seed_number;
int vertex_number;
int edge_number;
int face_number;
int polyhedron_number;

vector<string> split(const string &str, const string &pattern);

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
    cout << "test is " << test << '\n'
         << "seed_number is " << seed_number << '\n'
         << "vertex_number is " << vertex_number << '\n'
         << "edge_number is " << edge_number << '\n'
         << "face_number is " << face_number << '\n'
         << "polyhedron_number is " << polyhedron_number << endl;
    read_TESS.close();

    string cell_array;
    vector<vector<double >> cell_centroid(seed_number, vector<double>(3));
    read_TESS.open("n100-id1.tess", ios::in);
    while (!read_TESS.eof()) {
        getline(read_TESS, line_text);
        if (line_text == "  *seed") {
            for (int number_seed = 0; number_seed < seed_number; number_seed++) {
                getline(read_TESS, cell_array);
                vector<string> cell_array_split = split(cell_array, " ");
                double centroid_x = stod(cell_array_split[1]);
                double centroid_y = stod(cell_array_split[2]);
                double centroid_z = stod(cell_array_split[3]);
                cell_centroid[number_seed] = {centroid_x, centroid_y, centroid_z};
            }
        }
    }
    cout << cell_centroid[1][2] << endl;
    read_TESS.close();

    string vertex_array;
    vector<vector<double >> vertex_cartesian(vertex_number, vector<double>(3));
    read_TESS.open("n100-id1.tess", ios::in);
    while (!read_TESS.eof()) {
        getline(read_TESS, line_text);
        if (line_text == " **vertex") {
            getline(read_TESS, skip_line_text);
            for (int number_vertex = 0; number_vertex < vertex_number; number_vertex++) {
                getline(read_TESS, vertex_array);
                vector<string> vertex_array_split = split(vertex_array, " ");
                double vertex_x = stod(vertex_array_split[1]);
                double vertex_y = stod(vertex_array_split[2]);
                double vertex_z = stod(vertex_array_split[3]);
                vertex_cartesian[number_vertex] = {vertex_x, vertex_y, vertex_z};
            }
            cout << "1" << endl;
        }
    }
    read_TESS.close();
    cout << vertex_cartesian[1][1] << endl;


    string edge_array_text;
    vector<vector<int >> edge_array(edge_number, vector<int>(2));
    read_TESS.open("n100-id1.tess", ios::in);
    while (!read_TESS.eof()) {
        getline(read_TESS, line_text);
        if (line_text == " **edge") {
            getline(read_TESS, skip_line_text);
            for (int number_edge = 0; number_edge < edge_number; number_edge++) {
                getline(read_TESS, edge_array_text);
                vector<string> edge_text_split = split(edge_array_text, " ");
                int edge_coordination_point_1 = stoi(edge_text_split[1]);
                int edge_coordination_point_2 = stoi(edge_text_split[2]);
                edge_array[number_edge] = {edge_coordination_point_1, edge_coordination_point_2};
            }
        }
    }
    read_TESS.close();
    cout << edge_array[1][1] << endl;


    string face_vertex_text;
    string face_equation_parameter_text;
    vector<vector<int >> face_vertex_array(face_number, vector<int>(3));
    vector<vector<double>> face_equation_parameters(face_number, vector<double>(4));
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
                int face_vertex_point_1 = stoi(face_vertex_text_split[1]);
                int face_vertex_point_2 = stoi(face_vertex_text_split[2]);
                int face_vertex_point_3 = stoi(face_vertex_text_split[3]);
                face_vertex_array[number_face] = {face_vertex_point_1, face_vertex_point_2, face_vertex_point_3};
                vector<string> face_eq_para_split = split(face_equation_parameter_text, " ");
                double face_eq_para_a = stod(face_eq_para_split[0]);
                double face_eq_para_b = stod(face_eq_para_split[1]);
                double face_eq_para_c = stod(face_eq_para_split[2]);
                double face_eq_para_d = stod(face_eq_para_split[3]);
                face_equation_parameters[number_face] = {face_eq_para_a, face_eq_para_b, face_eq_para_c,
                                                         face_eq_para_d};
            }
        }
    }
    read_TESS.close();
    cout << face_vertex_array[1][1] << " " << face_equation_parameters[1][1] << endl;

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
    cout << polyhedron_faces[1][1] << endl;
    return 0;

}

vector<string> split(const string &str, const string &pattern) {
    vector<string> result;
    int position_start = str.find_first_not_of(pattern);
    string new_str = str.substr(position_start, str.size());
    split(result, new_str, is_any_of(pattern), token_compress_on);
    return result;
}
