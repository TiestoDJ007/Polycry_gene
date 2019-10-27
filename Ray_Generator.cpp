//
// Created by a on 2019/10/27.
//
#include <iostream>
#include "Ray_File_Func.h"
using namespace std;
int main(){
    Ray_File_Func test("n100-id1.tess", "n100-id1.stpoly", "n100-id1.ori");
    cout<<test.poly_ori()[0]<<endl;
    return 0;
}
