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

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

using namespace std;
using namespace Eigen;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 3, bg::cs::cartesian> Point;
typedef std::pair<Point, int> Point_Index;
typedef bg::model::box<Point> Box;

inline Vector4d Face_Para(const Ref<Vector3d> &point_1, const Ref<Vector3d> &point_2, const Ref<Vector3d> &point_3);

inline bool Point_Signal(const Ref<Vector4d> &face_para, const Ref<Vector3d> &point);

const double lattice_parameter = 3.6149;

const double min_distance = 0.425 * lattice_parameter;

int main() {

    clock_t begin = clock();


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
    vector<vector<int >> poly_faces = voro_info.poly_faces();
    vector<double> poly_eqradius = voro_info.poly_eqradius();
    vector<vector<int >> poly_vertices = voro_info.poly_vertices();
    vector<Vector4d> poly_ori = voro_info.poly_ori();

    vector<Vector4d> face_para;
    for (int i_face = 0; i_face < face_vertices.size(); ++i_face) {
        face_para.emplace_back(Face_Para(vertex_position[face_vertices[i_face][0]],
                                         vertex_position[face_vertices[i_face][1]],
                                         vertex_position[face_vertices[i_face][2]]));
    }

    vector<Point_Index> atoms_outfile;

    for (int i_poly = 0; i_poly < cell_tot; ++i_poly) {

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

        for (int i_atom = 0; i_atom < pre_cell_atoms.size(); ++i_atom) {
            VectorXi Atom_Signal(poly_faces[i_poly].size());
            for (int i_face = 0; i_face < poly_faces[i_poly].size(); ++i_face) {
                Atom_Signal.row(i_face) << Point_Signal(face_para[poly_faces[i_poly][i_face]],
                                                        pre_cell_atoms[i_atom]);
            }
            if (Atom_Signal == Signal_point)
                atoms_outfile.emplace_back(make_pair(Point(
                        pre_cell_atoms[i_atom](0, 0),
                        pre_cell_atoms[i_atom](1, 0),
                        pre_cell_atoms[i_atom](2, 0)), i_poly + 1));
        }
    }
    clock_t Pre_Gener = clock();

    cout << "Pre Generator Done" << '\n';
    cout << "Pre Generator Time: " << (Pre_Gener - begin) * 1.0 / CLOCKS_PER_SEC * 1000 << " ms" << '\n';
    cout << "Number of Pre Generator Atoms: " << atoms_outfile.size() << '\n';
    bgi::rtree<Point_Index, bgi::quadratic<16>> Rtree;
    for (int i_atom = 0; i_atom < atoms_outfile.size(); ++i_atom) {
        Rtree.insert(make_pair(atoms_outfile[i_atom].first, i_atom));
    }

    vector<vector<int>> atom_index;
    for (int i_atom = 0; i_atom < atoms_outfile.size(); ++i_atom) {
        //不能用box方法，肯定会发生误判，box中两点最长为sqrt(3)*radius,肯定超过判断条件。
        //换做KNN筛选方法。
        vector<Point_Index> nearest_atom;
        Rtree.query(bgi::nearest(atoms_outfile[i_atom].first, 5), std::back_inserter(nearest_atom));
        for (int i_near = 0; i_near < nearest_atom.size(); ++i_near) {
            if (bg::distance(nearest_atom[i_near].first, atoms_outfile[i_atom].first) < min_distance &&
                bg::distance(nearest_atom[i_near].first, atoms_outfile[i_atom].first) > 0.0) {
                vector<int> back_index{i_atom, nearest_atom[i_near].second};
                sort(back_index.begin(), back_index.end());
                atom_index.emplace_back(back_index);
            }
        }
    }

    clock_t check_pair = clock();
    cout << "Un-Check Pair Number: " << atom_index.size() << '\n';
    cout << "Select Atoms Time: " << (check_pair - Pre_Gener) * 1.0 / CLOCKS_PER_SEC * 1000 << " ms" << '\n';

    for (int i = 0; i < atom_index.size(); ++i) {
        vector<int> initial_val = atom_index[i];
        for (int j = 0; j < atom_index.size(); ++j) {
            if (initial_val == atom_index[j] && i != j)
                atom_index.erase(atom_index.begin() + j);
        }
    }

    vector<int> Delete_Index;
    int back_index = -1;
    for (int i_index = 0; i_index < atom_index.size(); ++i_index) {
        if (back_index != atom_index[i_index][0])
            Delete_Index.emplace_back(atom_index[i_index][0]);
        back_index = atom_index[i_index][0];
    }
    sort(Delete_Index.begin(),Delete_Index.end());
    cout << "Number of Deleted Atoms: " << Delete_Index.size() << '\n';
    clock_t check_out = clock();

    cout << "Check Atom-Couple Time: " << (check_out - check_pair) * 1.0 / CLOCKS_PER_SEC * 1000 << " ms" << '\n';

    for (int i_delete = Delete_Index.size() - 1; i_delete >= 0; i_delete = i_delete - 1) {
        atoms_outfile.erase(atoms_outfile.begin() + Delete_Index[i_delete]);
    }


    clock_t starting_write = clock();
    ofstream outdata;
    outdata.precision(6);
    outdata.open("Selected_Atoms.dat", ios::out);
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
        outdata << number_atom + 1 << " " << atoms_outfile[number_atom].second;
        outdata << " " << atoms_outfile[number_atom].first.get<0>() << " " << atoms_outfile[number_atom].first.get<1>()
                << " " << atoms_outfile[number_atom].first.get<2>() << "\n";
    }
    outdata << endl;
    outdata.close();

    clock_t writing_data = clock();
    cout << "Writing File Time: " << (writing_data - starting_write) * 1.0 / CLOCKS_PER_SEC * 1000 << " ms" << '\n';
    cout << "Total Time: " << (clock() - begin) * 1.0 / CLOCKS_PER_SEC * 1000 << " ms" << endl;

    return 0;
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
    return ((double) face_para(0, 0) * point(0, 0) +
            face_para(1, 0) * point(1, 0) +
            face_para(2, 0) * point(2, 0) +
            face_para(3, 0)) > 0.0;
}