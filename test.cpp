//
// Created by a on 2019/10/20.
//
#include <iostream>
#include<eigen3/Eigen/Dense>
using namespace std;
using namespace Eigen;
int main(){
    MatrixXi a=MatrixXi::Ones(5,5);
    MatrixXi b=MatrixXi::Zero(5,5);
    MatrixXi c=MatrixXi::Zero(5,5);
    if(a==b)
        cout<<"yes"<<endl;
    else
        cout<<"no"<<endl;
    //cout<<a<<endl;
}
