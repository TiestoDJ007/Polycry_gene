//
// Created by a on 2019/10/20.
//
#include <iostream>
#include <list>
#include <vector>
#include <eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
using namespace std;
using namespace Eigen;
int main(){
    MatrixXi a=MatrixXi::Ones(5,5);
    MatrixXi b=MatrixXi::Zero(5,5);
    MatrixXi c=MatrixXi::Zero(5,5);
    vector<MatrixXi> xxx;
    xxx.emplace_back(a);
    xxx.emplace_back(b);
    cout<<xxx[1]<<endl;
    //cout<<a<<endl;
}
