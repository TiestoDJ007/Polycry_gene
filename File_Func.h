//
// Created by a on 2019/10/27.
//

#ifndef POLYCRY_GENE_FILE_FUNC_H
#define POLYCRY_GENE_FILE_FUNC_H

#include <vector>
#include <string>
#include <fstream>
#include <eigen3/Eigen/SVD>


class File_Func {
private:
    std::ifstream read_file;
    int cell_tot_;
    int vertex_tot_;
    int face_tot_;
    int polyhedron_tot_;
    std::string line_text;
    std::string file_name_tess;
    std::string file_name_stopoly;
    std::string file_name_ori;
public:
    File_Func();

    File_Func(const std::string &file_name_tess_, const std::string &file_name_stopoly_,
              const std::string &file_name_ori_);

    ~File_Func();

    int cell_tot();

    std::vector<Eigen::Vector3d> cell_position();

    std::vector<Eigen::Vector3d> vertex_position();

    std::vector<std::vector<int >> face_vertices();

    std::vector<std::vector<int >> poly_faces();

    std::vector<double> poly_eqradius();

    std::vector<Eigen::Vector4d> poly_ori();

    std::vector<std::vector<int >> poly_vertices();

};

inline std::vector<std::string> split_string(const std::string &str, const std::string &pattern);

#endif //POLYCRY_GENE_FILE_FUNC_H
