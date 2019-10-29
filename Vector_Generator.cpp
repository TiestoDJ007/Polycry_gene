//neper -T -n 100 /
//      -id 1     /
//      -reg 1    /
//      -morpho "diameq:dirac(1),sphericity:lognormal(0.145,0.03,1-x)"  /
//      -domain "cube(200,200,200)" /
//      -statpoly "radeq,verlist"   /
//      -oridescriptor q
#include <iostream>
#include <eigen3/Eigen/Geometry>
#include "File_Func.h"

using namespace std;
using namespace Eigen;

inline double
Distance_Point_Face(const Eigen::Ref<Eigen::Vector4d> &face_para, const Eigen::Ref<Eigen::Vector3d> &point);

inline Vector4d Face_Para(const Ref<Vector3d> &point_1, const Ref<Vector3d> &point_2, const Ref<Vector3d> &point_3);

inline bool Point_Signal(const Ref<Vector4d> &face_para, const Ref<Vector3d> &point);

int main() {

    clock_t begin = clock();

    double lattice_parameter = 3.6149;

    MatrixXd lattice_constant_array(4, 3);
    lattice_constant_array << 0.0, 0.0, 0.0,
            0.5, 0.5, 0.0,
            0.5, 0.0, 0.5,
            0.0, 0.5, 0.5;

    Matrix3d cell_martix;
    cell_martix << 1, 0, 0,
            0, 1, 0,
            0, 0, 1;

    File_Func voro_info("n100-id1.tess", "n100-id1.stpoly", "n100-id1.ori");

    int cell_tot = voro_info.cell_tot();
    vector<Vector3d> cell_position = voro_info.cell_position();
    vector<Vector3d> vertex_position = voro_info.vertex_position();
    vector<vector<int >> face_vertices = voro_info.face_vertices();
    cout<<face_vertices[1][0]<<endl;
    vector<vector<int >> poly_faces = voro_info.poly_faces();
    cout<<poly_faces[1][0]<<endl;
    vector<double> poly_eqradius = voro_info.poly_eqradius();
    vector<vector<int >> poly_vertices = voro_info.poly_vertices();
    cout<<poly_vertices[1][0]<<endl;
    vector<Vector4d> poly_ori = voro_info.poly_ori();

    vector<MatrixXd> face_vertex;
    for (int i_face = 0; i_face < face_vertices.size(); ++i_face) {
        MatrixXd back_Matrix(face_vertices[i_face].size(), 3);
        for (int i_vex = 0; i_vex < face_vertices[i_face].size(); ++i_vex) {
            back_Matrix.row(i_vex) = vertex_position[face_vertices[i_face][i_vex]].transpose();
        }
        face_vertex.emplace_back(back_Matrix);
    }

    vector<Vector4d> face_para;
    for (int i_face = 0; i_face < face_vertices.size(); ++i_face) {
        face_para.emplace_back(Face_Para(vertex_position[face_vertices[i_face][0]],
                                         vertex_position[face_vertices[i_face][1]],
                                         vertex_position[face_vertices[i_face][2]]));
    }
    cout<<face_para[5].transpose()<<endl;

    vector<vector<double >> atoms_outfile;

    for (int i_poly = 0; i_poly < 0; ++i_poly) {

        int cubic_size = 2 * poly_eqradius[i_poly] / lattice_parameter;

        MatrixXd poly_vertex(poly_vertices[i_poly].size(), 3);
        for (int i_vertex = 0; i_vertex < poly_vertices[i_poly].size(); ++i_vertex) {
            poly_vertex(i_vertex, 0) = vertex_position[poly_vertices[i_poly][i_vertex]](0, 0);
            poly_vertex(i_vertex, 1) = vertex_position[poly_vertices[i_poly][i_vertex]](1, 0);
            poly_vertex(i_vertex, 2) = vertex_position[poly_vertices[i_poly][i_vertex]](2, 0);
        }

        double min_x = poly_vertex.col(0).minCoeff();
        double min_y = poly_vertex.col(1).minCoeff();
        double min_z = poly_vertex.col(2).minCoeff();
        double max_x = poly_vertex.col(0).maxCoeff();
        double max_y = poly_vertex.col(1).maxCoeff();
        double max_z = poly_vertex.col(2).maxCoeff();

        Matrix3d rotation_Matrix;
        rotation_Matrix = Matrix3d(Quaterniond(poly_ori[i_poly]));

        vector<Vector3d> pre_cell_atoms;
        for (int i_x = -cubic_size; i_x < cubic_size + 1; ++i_x) {
            for (int i_y = -cubic_size; i_y < cubic_size + 1; ++i_y) {
                for (int i_z = -cubic_size; i_z < cubic_size + 1; ++i_z) {
                    Vector3d base_point;
                    base_point << i_x, i_y, i_z;
                    Vector3d cart_point;
                    cart_point = cell_martix * base_point;
                    for (int i_lattice = 0; i_lattice < 4; ++i_lattice) {
                        Vector3d origin_atom;
                        origin_atom =
                                rotation_Matrix * (cart_point + lattice_constant_array.row(i_lattice).transpose()) *
                                lattice_parameter;
                        Vector3d return_atom = cell_position[i_poly] + origin_atom;
                        if (return_atom(0, 0) > min_x && return_atom(0, 0) < max_x &&
                            return_atom(1, 0) > min_y && return_atom(1, 0) < max_y &&
                            return_atom(2, 0) > min_z && return_atom(2, 0) < max_z) {
                            pre_cell_atoms.emplace_back(return_atom);
                        }
                    }
                }
            }
        }

        VectorXi Signal_point(poly_faces[i_poly].size());
        for (int i_face = 0; i_face < poly_faces[i_poly].size(); ++i_face) {
            Signal_point.row(i_face) << Point_Signal(face_para[poly_faces[i_poly][i_face]],
                                                     cell_position[i_poly]);
        }
        cout << pre_cell_atoms.size() << endl;
        cout << Signal_point.transpose() << endl;

        for (int i_atom = 0; i_atom < pre_cell_atoms.size(); ++i_atom) {
            VectorXi Atom_Signal(poly_faces[i_poly].size());
            for (int i_face = 0; i_face < poly_faces[i_poly].size(); ++i_face) {
                Atom_Signal.row(i_face) << Point_Signal(face_para[poly_faces[i_poly][i_face]],
                                                        pre_cell_atoms[i_atom]);
            }
            if (Atom_Signal == Signal_point)
                atoms_outfile.emplace_back(vector<double>{(double) i_poly + 1,
                                                          pre_cell_atoms[i_atom](0, 0),
                                                          pre_cell_atoms[i_atom](1, 0),
                                                          pre_cell_atoms[i_atom](2, 0)});
        }
    }
    cout << atoms_outfile.size() << endl;
    cout << (clock() - begin) * 1.0 / CLOCKS_PER_SEC * 1000 << " ms" << endl;

    ofstream outdata;
    outdata.precision(6);
    outdata.open("test_1.dat", ios::out);
    outdata << "Crystalline Cu atoms\n\n";
    outdata << atoms_outfile.size() << " atoms\n";
    outdata << cell_tot << " atom types\n";
    outdata << fixed;
    outdata << "0" << " " << "200" << " xlo xhi\n";
    outdata << "0" << " " << "200" << " ylo yhi\n";
    outdata << "0" << " " << "200" << " zlo zhi\n";
    outdata << "\n";
    outdata << "Atoms\n\n";
    for (int number_atom = 0; number_atom < atoms_outfile.size(); ++number_atom) {
        outdata << number_atom + 1 << " " << (int) atoms_outfile[number_atom][0];
        for (int number_out = 1; number_out < 4; ++number_out) {
            outdata << " " << atoms_outfile[number_atom][number_out];
        }
        outdata << endl;
    }
    outdata.close();
    return 0;
}


inline double
Distance_Point_Face(const Eigen::Ref<Eigen::Vector4d> &face_para, const Eigen::Ref<Eigen::Vector3d> &point) {
    return (face_para(0, 0) * point(0, 0) +
            face_para(1, 0) * point(1, 0) +
            face_para(2, 0) * point(2, 0) +
            face_para(3, 0)) /
           sqrt(face_para(0, 0) * face_para(0, 0) +
                face_para(1, 0) * face_para(1, 0) +
                face_para(2, 0) * face_para(2, 0));
}

inline Vector4d Face_Para(const Ref<Vector3d> &p_1, const Ref<Vector3d> &p_2, const Ref<Vector3d> &p_3) {
    Vector4d data_return;
    double a = (p_2(1, 0) - p_1(1, 0)) * (p_3(2, 0) - p_1(2, 0)) -
               (p_3(1, 0) - p_1(1, 0)) * (p_2(2, 0) - p_1(2, 0));
    double b = (p_2(2, 0) - p_1(2, 0)) * (p_3(0, 0) - p_1(0, 0)) -
               (p_3(2, 0) - p_1(2, 0)) * (p_2(0, 0) - p_1(0, 0));
    double c = (p_2(0, 0) - p_1(0, 0)) * (p_3(1, 0) - p_1(1, 0)) -
               (p_3(0, 0) - p_1(0, 0)) * (p_2(1, 0) - p_1(1, 0));
    double norm = sqrt(a * a + b * b + c * c);
    double d = -a * p_1(0, 0) - b * p_1(1, 0) - c * p_1(2, 0);
    data_return << a / norm, b / norm, c / norm, d / norm;
    return data_return;
}

inline bool Point_Signal(const Ref<Vector4d> &face_para, const Ref<Vector3d> &point) {
    return (face_para(0, 0) * point(0, 0) +
            face_para(1, 0) * point(1, 0) +
            face_para(2, 0) * point(2, 0) +
            face_para(3, 0)) > 0;
}