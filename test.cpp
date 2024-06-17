#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include "face.h"


int main(int argc, char const *argv[])
{
    ostringstream node_out;
    ostringstream triangle_out;
    int num_nodes = 0;
    int num_triangles = 0;


    NonuniformFace face_test("x",0,7,7,2,2,2,1);
    face_test.genNodes();
    face_test.shiftNodes(-7,-7,-7);
    face_test.genTriangles();
    face_test.genConnectorsBottom(7);
    face_test.genConnectorsTop(7);
    face_test.genOutput(num_nodes,num_triangles,node_out,triangle_out,true);

    cout << num_nodes << '\t' << 3 << '\t' << 0 << '\t' << 0 << endl;
    
    cout << node_out.str();

    cout << triangle_out.str();

    return 0;
}
