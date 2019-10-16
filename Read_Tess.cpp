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
//int polyhedron_faces_number;

//vector<vector<double >>cell_centroid(cell_number, vector<double>(3));
vector<string> split(const string &str, const string &pattern);

int main() {
    //读取.tess文件
    ifstream read_TESS;
    read_TESS.open("n100-id1.tess", ios::in);
    string line_text;
    //string line_text_1;
    string cell_text;
    string vertex_text;
    string edge_text;
    string face_text;
    string face_equation_parameter_text;
    string face_vertex_text;
    string polyhedron_text;
    string polyhedron_face_text;
    string skip_line_text;
    string cell_array;
    string vertex_array;
    int test = 0;
    while (getline(read_TESS, line_text)) {
        //getline(read_TESS, line_text);
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


    vector<vector<double >> cell_centroid(seed_number, vector<double>(3));
    read_TESS.open("n100-id1.tess", ios::in);
    while (!read_TESS.eof()) {
        getline(read_TESS, line_text);
        if (line_text == "  *seed") {
            for (int number_seed = 0; number_seed < seed_number; number_seed++) {
                getline(read_TESS, cell_array);
                //cout << cell_array << endl;
                vector<string> cell_array_split = split(cell_array, " ");
                double centroid_x = stod(cell_array_split[1]);
                double centroid_y = stod(cell_array_split[2]);
                double centroid_z = stod(cell_array_split[3]);
                cell_centroid[number_seed] = {centroid_x, centroid_y, centroid_z};
                //cout << " " << centroid_x << " " << centroid_y << " " << centroid_z << endl;
            }
        }
    }
    cout << cell_centroid[1][2] << endl;
    read_TESS.close();

    vector<vector<double >> vertex_cartesian(vertex_number, vector<double>(3));
    read_TESS.open("n100-id1.tess", ios::in);
    while (!read_TESS.eof()) {
        if (line_text.compare(" **vertex") == 0) {
            getline(read_TESS, skip_line_text);
            for (int number_vertex = 0; number_vertex < vertex_number; number_vertex++) {
                getline(read_TESS, vertex_array);
                cout << vertex_array << endl;
                vector<string> vertex_array_split = split(vertex_array, " ");
                double vertex_x = stod(vertex_array_split[1]);
                double vertex_y = stod(vertex_array_split[2]);
                double vertex_z = stod(vertex_array_split[3]);
                vertex_cartesian[number_vertex] = {vertex_x, vertex_y, vertex_z};
                cout << vertex_array_split[0] << " " << vertex_array_split[1] << " " << vertex_array_split[4] << endl;
            }
        }
    }
    read_TESS.close();


    return 0;

}

vector<string> split(const string &str, const string &pattern) {
    vector<string> result;
    int position_start = str.find_first_not_of(pattern);
    string new_str = str.substr(position_start, str.size());
    split(result, new_str, is_any_of(pattern), token_compress_on);
    return result;
}
