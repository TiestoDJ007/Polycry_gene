//
// Created by a on 2019/10/27.
//
#include <cmath>
#include "Ray_File_Func.h"
#include <boost/algorithm/string.hpp>

std::ifstream read_file;

Ray_File_Func::Ray_File_Func() {
    cell_tot_ = 0;
    vertex_tot_ = 0;
    face_tot_ = 0;
    polyhedron_tot_ = 0;
    file_name_tess = "none_file_name";
    file_name_stopoly = "none_file_name";
    file_name_ori = "none_file_name";
}

Ray_File_Func::Ray_File_Func(const std::string &file_name_tess_, const std::string &file_name_stopoly_,
                             const std::string &file_name_ori_) {
    file_name_tess = file_name_tess_;
    file_name_stopoly = file_name_stopoly_;
    file_name_ori = file_name_ori_;
    read_file.open(file_name_tess, std::ios::in);
    while (getline(read_file, line_text)) {
        if (line_text == " **cell") {
            getline(read_file, line_text);
            cell_tot_ = stoi(line_text, nullptr, 10);
        }
        if (line_text == " **vertex") {
            getline(read_file, line_text);
            vertex_tot_ = stoi(line_text, nullptr, 10);
        }
        if (line_text == " **edge") {
            getline(read_file, line_text);
            face_tot_ = stoi(line_text, nullptr, 10);
        }
        if (line_text == " **face") {
            getline(read_file, line_text);
            face_tot_ = stoi(line_text, nullptr, 10);
        }
        if (line_text == " **polyhedron") {
            getline(read_file, line_text);
            polyhedron_tot_ = stoi(line_text, nullptr, 10);
        }
    }
    read_file.close();
}

Ray_File_Func::~Ray_File_Func() {

}

int Ray_File_Func::cell_tot() {
    return cell_tot_;
}

int Ray_File_Func::vertex_tot() {
    return vertex_tot_;
}

int Ray_File_Func::face_tot() {
    return face_tot_;
}

int Ray_File_Func::polyhedron_tot() {
    return polyhedron_tot_;
}

std::vector<Eigen::Vector3d> Ray_File_Func::cell_position() {
    std::vector<Eigen::Vector3d> return_data;
    read_file.open(file_name_tess, std::ios::in);
    while (!read_file.eof()) {
        getline(read_file, line_text);
        if (line_text == "  *seed") {
            for (int number_seed = 0; number_seed < cell_tot_; number_seed++) {
                getline(read_file, line_text);
                std::vector<std::string> line_text_split = split_string(line_text, " ");
                Eigen::Vector3d back_vector;
                back_vector << stod(line_text_split[1]),
                        stod(line_text_split[2]),
                        stod(line_text_split[3]);
                return_data.emplace_back(back_vector);
            }
        }
    }
    read_file.close();
    return return_data;
}

std::vector<Eigen::Vector3d> Ray_File_Func::vertex_position() {
    std::vector<Eigen::Vector3d> return_data;
    read_file.open(file_name_tess, std::ios::in);
    while (!read_file.eof()) {
        getline(read_file, line_text);
        if (line_text == " **vertex") {
            getline(read_file, line_text);
            for (int i_vertex = 0; i_vertex < vertex_tot_; ++i_vertex) {
                getline(read_file, line_text);
                std::vector<std::string> line_text_split = split_string(line_text, " ");
                Eigen::Vector3d back_vector;
                back_vector << stod(line_text_split[1]),
                        stod(line_text_split[2]),
                        stod(line_text_split[3]);
                return_data.emplace_back(back_vector);
            }
        }
    }
    read_file.close();
    return return_data;
}

std::vector<std::vector<int >> Ray_File_Func::face_vertices() {
    std::vector<std::vector<int >> return_data;
    read_file.open(file_name_tess, std::ios::in);
    while (!read_file.eof()) {
        getline(read_file, line_text);
        if (line_text == " **face") {
            getline(read_file, line_text);
            for (int i_face = 0; i_face < face_tot_; ++i_face) {
                getline(read_file, line_text);
                std::vector<std::string> line_text_split = split_string(line_text, " ");
                std::vector<int> line_back;
                for (int i_vertex = 2; i_vertex < line_text_split.size(); ++i_vertex) {
                    line_back.emplace_back(std::stoi(line_text_split[i_vertex]) - 1);
                }
                return_data.emplace_back(line_back);
                getline(read_file, line_text);
                getline(read_file, line_text);
                getline(read_file, line_text);
            }
        }
    }
    read_file.close();
    return return_data;
}

std::vector<std::vector<int >> Ray_File_Func::poly_faces() {
    std::vector<std::vector<int >> return_data;
    read_file.open(file_name_tess, std::ios::in);
    while (!read_file.eof()) {
        getline(read_file, line_text);
        if (line_text == " **polyhedron") {
            getline(read_file, line_text);
            for (int i_poly = 0; i_poly < polyhedron_tot_; ++i_poly) {
                getline(read_file, line_text);
                std::vector<std::string> line_text_split = split_string(line_text, " ");
                std::vector<int> line_back;
                for (int i_face = 2; i_face < line_text_split.size(); ++i_face) {
                    line_back.emplace_back(std::abs(std::stoi(line_text_split[i_face]) - 1));
                }
                return_data.emplace_back(line_back);
            }
        }
    }
    read_file.close();
    return return_data;
}

std::vector<double> Ray_File_Func::poly_eqradius() {
    std::vector<double> return_data;
    read_file.open(file_name_stopoly, std::ios::in);
    for (int i_radius = 0; i_radius < cell_tot_; ++i_radius) {
        getline(read_file, line_text);
        return_data.emplace_back(std::stod(split_string(line_text, " ")[0]));
    }
    read_file.close();
    return return_data;
}

std::vector<std::vector<int >> Ray_File_Func::poly_vertices() {
    std::vector<std::vector<int >> return_data;
    read_file.open(file_name_stopoly, std::ios::in);
    for (int i_poly = 0; i_poly < cell_tot_; ++i_poly) {
        getline(read_file, line_text);
        std::vector<std::string> line_text_split = split_string(line_text, " ");
        std::vector<int> line_back;
        for (int i_ver = 1; i_ver < line_text_split.size(); ++i_ver) {
            line_back.emplace_back(std::stoi(line_text_split[i_ver]) - 1);
        }
        return_data.emplace_back(line_back);
    }
    read_file.close();
    return return_data;
}

std::vector<Eigen::Vector4d> Ray_File_Func::poly_ori() {
    std::vector<Eigen::Vector4d> return_data;
    read_file.open(file_name_ori, std::ios::in);
    for (int i_poly = 0; i_poly < cell_tot_; ++i_poly) {
        getline(read_file, line_text);
        std::vector<std::string> line_text_split = split_string(line_text, " ");
        Eigen::Vector4d back_vector;
        back_vector << std::stod(line_text_split[0]),
                std::stod(line_text_split[1]),
                std::stod(line_text_split[2]),
                std::stod(line_text_split[3]);
        return_data.emplace_back(back_vector);
    }
    read_file.close();
    return return_data;
}

inline std::vector<std::string> split_string(const std::string &str, const std::string &pattern) {
    std::vector<std::string> result_str;
    int position_start = str.find_first_not_of(pattern);
    std::string new_str = str.substr(position_start, str.size());
    boost::split(result_str, new_str, boost::is_any_of(pattern), boost::token_compress_on);
    return result_str;
}

