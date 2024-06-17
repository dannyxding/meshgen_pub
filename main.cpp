#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <limits>
#include <cmath>
#include <unordered_map>
#include <queue>
#include "face.h"
using namespace std;

// Static Config (DIR)
static std::filesystem::path OFF_DIR = "./off";
static std::filesystem::path POLY_DIR = "./poly";
static std::filesystem::path MESH_DIR = "./mesh";
static std::filesystem::path TMSMESH_DIR = "./tool_surface/tmsmesh";
static std::filesystem::path TETGEN_DIR = "./tool_mesh/tetgen";

// DEBUG
static int PROTEIN_LABEL = 2;
static int MEMBRANE_LABEL = 3;
static int SOLVENT_LABEL = 1;


struct MeshNode {
    double x, y, z;
    MeshNode(double x_in, double y_in , double z_in) :  x(x_in), y(y_in), z(z_in) {}
    MeshNode() :  x(0.0), y(0.0), z(0.0) {}
};

// For SeararchAlg
struct MeshTetragonSearch {
    // Definition of label here: -1 - Undefined, 0/1 - might be membrane or solvent region yet to be determined, 2 - Protein region
    int node1, node2, node3, node4, neigh1, neigh2, neigh3, neigh4, label;
    bool visited;

    MeshTetragonSearch(int node1_in, int node2_in, int node3_in, int node4_in, int neigh1_in = -1, int neigh2_in = -1, int neigh3_in = -1, int neigh4_in = -1, int label_in = -1) 
    : node1(node1_in), node2(node2_in), node3(node3_in), node4(node4_in), neigh1(neigh1_in), neigh2(neigh2_in), neigh3(neigh3_in), neigh4(neigh4_in), visited(false) {
        if(label_in == PROTEIN_LABEL){
            label = 2;
        }else{
            label = -1;
        }
    }
    MeshTetragonSearch() : node1(0), node2(0), node3(0), node4(0), neigh1(-1), neigh2(-1), neigh3(-1), neigh4(-1), label(-1), visited(false) {}
};


struct MeshTetragon{
    int node1, node2, node3, node4, label;
    MeshTetragon(int node1_in, int node2_in, int node3_in, int node4_in, int label_in) : node1(node1_in), node2(node2_in), node3(node3_in), node4(node4_in), label(label_in) {}
    MeshTetragon() : node1(0), node2(0), node3(0), node4(0), label(-1) {}

    void print_tet(ostringstream& out_s, int n){
        out_s << n << '\t' << node1 << '\t' << node2 <<  '\t' << node3 << '\t' << node4 << '\t' << label << endl;
    }
};


// Main Data (cetnral region used for search)
// central_ndoes refer to central area in XYZ three direction that we are gonna use for 
std::unordered_map<int, MeshNode> central_nodes;
std::unordered_map<int, MeshNode> middle_nodes;
std::unordered_map<int, MeshNode> memsurface_nodes;
std::unordered_map<int, MeshTetragonSearch> central_tetragon;
// Keep tracking if the label 0 group touched the boundary (belong to the Membrane region..)

bool node_is_central(int n_node){
    return central_nodes.find(n_node) != central_nodes.end();
}

bool node_is_middle(int n_node){
     return middle_nodes.find(n_node) != middle_nodes.end();
}

bool node_is_memsurface(int n_node){
    return memsurface_nodes.find(n_node) != memsurface_nodes.end();
}

bool tet_is_middle(MeshTetragon tet){
    return node_is_middle(tet.node1) && node_is_middle(tet.node2) && node_is_middle(tet.node3) && node_is_middle(tet.node4);
}

bool tet_is_central(MeshTetragon tet){
    return tet_is_middle(tet) &&  (node_is_central(tet.node1) || node_is_central(tet.node2) || node_is_central(tet.node3) || node_is_central(tet.node4));
}

bool tet_is_middle(MeshTetragonSearch tet){
    return node_is_middle(tet.node1) && node_is_middle(tet.node2) && node_is_middle(tet.node3) && node_is_middle(tet.node4);
}

bool tet_is_central(MeshTetragonSearch tet){
    return tet_is_middle(tet) &&  (node_is_central(tet.node1) || node_is_central(tet.node2) || node_is_central(tet.node3) || node_is_central(tet.node4));
}

/* Helper Function for BST
 * Given index of a neighbour triangle, retrurn ture if can keep searching...
 * n_curr_tri must be an index of valid VISITED CENTRAL triangle
 * Makesure used only for LABEL:0   update the touch_boundary variable
*/ 
bool real_neighbor(int n_neighbor_tri, bool& touch_boundary){
    MeshTetragonSearch neighbor = central_tetragon[n_neighbor_tri];
    if(!tet_is_middle(neighbor)){
        return false;
    }
    if(n_neighbor_tri == -1 || !tet_is_central(neighbor)){
        // This Neighbor not exist / isn't central...
        touch_boundary = true;
        return false;
    }
    return (!neighbor.visited);
}



/**
 * Helper Functin for BST
 * n_tretragon must be an index of valid (UNVISITED) CENTRAL triangle
*/
void visit_triangle(int n_tretragon, int label){
    central_tetragon[n_tretragon].visited = true;
    central_tetragon[n_tretragon].label  = label;

    //cout << n_tretragon << endl;
    //cout << n_tretragon << '\t' << central_tetragon[n_tretragon].node1 << '\t' << central_tetragon[n_tretragon].node2 << '\t' << central_tetragon[n_tretragon].node3 << '\t' << central_tetragon[n_tretragon].node4 << endl;
}


string double_to_string(double d){
    ostringstream out;
    out << d;
    std::string str = out.str();

    str.erase(str.find_last_not_of('0') + 1, std::string::npos);

    if (!str.empty() && str.back() == '.') {
        str.pop_back();
    }

    return str;
}

bool is_in_range(double num, double lower, double upper) {
    return num >= lower && num <= upper;
}


// Write a nodes to the output stream (in .poly format)
void add_node(ostringstream& nodes_out, int& num_nodes, double x, double y, double z){
    nodes_out <<  (num_nodes + 1) << '\t' << x << '\t' << y << '\t' << z << '\n';
    num_nodes++;
}

// Add a facet to the .poly file with a single triangle..
void add_facet(ostringstream& triangle_out, int n_a, int n_b, int n_c){
    triangle_out << "1\n";
    triangle_out << "3\t" << n_a << '\t' << n_b << '\t' << n_c << '\n';
}

/**
 * Generate the PROTEIN surface..
 * Return the ifstream of the .off file from the pqr file indicated, use TMSmesh to generate that .off file if such .off file doesn't exists yet.
 * pqr_file_name is the file name and we assume such file is placed under the root directory
 * 
*/
ifstream get_off_file(string pqr_file_name, string str_tms_d, string str_tms_e){
     //Check if pqr file exists
    ifstream pqr_file(pqr_file_name);
    if(pqr_file.fail()){
        cout << "Could not find the pqr file: " << pqr_file_name << endl;
        exit(1);
    }
    pqr_file.close();

    // Check if corresponding off file exists, if it do, we can skip the TMSmesh part
    string off_file_name = pqr_file_name + '-' + str_tms_d + '_' + str_tms_e + ".off_modify.off";
    std::filesystem::path off_file_path = OFF_DIR / off_file_name;
    ifstream off_file(off_file_path);
    if(off_file.fail()){
        // Use TMSmesh (TMSMesh2.1) to generate the protein surface mesh (.off)
        cout << "Running TMSmesh...\n";
        string t_cmd = TMSMESH_DIR.string() +  "/TMSmesh2.1 "  + pqr_file_name + ' ' + str_tms_d + ' ' + str_tms_e;
        // cout << t_cmd << endl;
        int result = std::system(t_cmd.c_str());
        if(result!=0){
            cout << "Failed to run TMSMesh2.1, check configuration" << endl;
            exit(1);
        }
        // Move .off file to the OFF_DIR
        std::filesystem::rename(off_file_name, off_file_path);
        off_file.open(off_file_path); // This time it should succeed.
        if(off_file.fail()){
            cout << "TMSMesh run successfully but failed to oepn the file: " << off_file_path << endl;
            exit(1); 
        }
    }

    cout << ".off file opened successfully" << endl;
    return off_file;

}


/**
 * Determine the boundary of the box domain through reading the protein .off file, figure out the boundary, add ETA and neccessary padding..
*/
void determine_boundary(ifstream& protein_off_file, int ETA_X, int ETA_Y, int ETA_Z, int H_M, int H_S, int Z_1, int Z_2, int Connecting_len,
 int& L_x1, int& L_x2, int& L_y1, int& L_y2, int& L_z1, int& L_z2){
        // Read the off_file once to find the boundry of nodes' coord
    string line;
    getline(protein_off_file, line); // reads "OFF"
    getline(protein_off_file, line); // read first line..
    istringstream s_line(line);
    int protein_num_nodes, protein_num_triangles;
    s_line >> protein_num_nodes >> protein_num_triangles;
    // Iterate through the nodes' list
    float x_max = std::numeric_limits<float>::lowest();
    float y_max = std::numeric_limits<float>::lowest();
    float z_max = std::numeric_limits<float>::lowest();
    float x_min = std::numeric_limits<float>::max();
    float y_min = std::numeric_limits<float>::max();
    float z_min = std::numeric_limits<float>::max();
    for(size_t i = 0; i < protein_num_nodes; i++){
        getline(protein_off_file, line);
        istringstream stm_line(line);
        double x, y ,z;
        stm_line >> x >> y >> z;    
        x_max = (x > x_max)? x : x_max;
        y_max = (y > y_max)? y : y_max;
        z_max = (z > z_max)? z : z_max;
        x_min = (x < x_min)? x : x_min;
        y_min = (y < y_min)? y : y_min;
        z_min = (z < z_min)? z : z_min;
    }

    // cout << "Finds the boundary of protein meshes: " << x_max << " " << x_min << " " << y_max << " " << y_min << " " << z_max << " " << z_min << endl; 
    // Set the boundary of the box domain... Apply some padding, want to make sure the surfact box's boundary is some interger and that it's length are multiple of theh H_S(large mesh length size..)
    L_x1 = std::floor(x_min) - ETA_X;
    L_y1 = std::floor(y_min) - ETA_Y;
 
    L_x2 = L_x1 + std::ceil((x_max + ETA_X - L_x1)/H_S) * H_S;
    L_y2 = L_y1 + std::ceil((y_max + ETA_Y - L_y1)/H_S) * H_S;



    // Padding for Z direction...
    if((Z_2 - Z_1) % H_M != 0){
        cout << "Improper value for Membrane_boundary_Z_2 and  Membrane_boundary_Z_1, the width of membrane part should be multiple of H_M" << endl;
        exit(1);
    }
    L_z1 = std::floor(z_min) - ETA_Z;
    L_z2 = std::ceil(z_max) + ETA_Z;
    L_z1 = Z_1 - Connecting_len - ((Z_1 - Connecting_len - L_z1)/H_S)*H_S;
    L_z2 = Z_2 + Connecting_len + ((L_z2 - Z_2 - Connecting_len)/H_S)*H_S;

    cout << "Determined the boundary of box domain: " << L_x1 << " " << L_x2 << " " << L_y1 << " " << L_y2 << " " << L_z1 << " " << L_z2 << endl; 

    protein_off_file.clear();
    protein_off_file.seekg(0, std::ios::beg);
    // Reset file stream to the begeining

    return;
}


void generate_surface_mesh(ostringstream& node_out, ostringstream& triangle_out, int& num_nodes, int& num_triangles, int Membrane_boundary_Z1, int Membrane_boundary_Z2, int Connecting_len, int H_M, int H_S,
                           int& L_x1, int& L_x2, int& L_y1, int& L_y2, int& L_z1, int& L_z2){
    int x_len = L_x2 - L_x1;
    int y_len = L_y2 - L_y1;
    int z_len = L_z2 - L_z1;
    int x_shift = L_x1;
    int y_shift = L_y1;
    int z_shift = L_z1;


    if(Connecting_len % H_M != 0){
        cout << "Improper value for Connecting_len, it should be multiple of H_M" << endl;
        exit(1);
    }
    //Margin_len = (z_len / 2) - Connecting_len - Membrane_boundary_Z;
    
    int zDimNum_T = (L_z2 - Connecting_len - Membrane_boundary_Z2) / H_S; 
    int zDimNum_B = (Membrane_boundary_Z1 - Connecting_len - L_z1) / H_S;

    cout << "Divided the height (" << z_len << ") into several regions;  Membrane: " << -Membrane_boundary_Z1  << '+' << Membrane_boundary_Z2  << " Connecting(each): " << Connecting_len << "  Margin: " << zDimNum_T * H_S << '+' << zDimNum_B*H_S << endl;


    Face face_top(0, "z", z_len, x_len/H_S, y_len/H_S, 0, H_S);
    face_top.genNodes();
    face_top.shiftNodes(x_shift, y_shift, z_shift);
    face_top.genTriangles();
    face_top.genOutput(num_nodes, num_triangles, node_out, triangle_out, true, true);

    Face face_bottom(0, "z", 0, x_len/H_S, y_len/H_S, 0, H_S);
    face_bottom.genNodes();
    face_bottom.shiftNodes(x_shift, y_shift, z_shift);
    face_bottom.genTriangles();
    face_bottom.genOutput(num_nodes, num_triangles, node_out, triangle_out, true, true);

    NonuniformFace face_front("x",0,x_len/H_S, y_len/H_S, zDimNum_T, zDimNum_B, (Membrane_boundary_Z2-Membrane_boundary_Z1)/H_M, Connecting_len, H_M);
    face_front.genNodes();
    face_front.shiftNodes(x_shift, y_shift, z_shift);
    face_front.genTrianglesAndConnectors(y_len/H_S);
    face_front.genOutput(num_nodes, num_triangles, node_out, triangle_out, true);

    NonuniformFace face_back("x",x_len,x_len/H_S, y_len/H_S, zDimNum_T, zDimNum_B, (Membrane_boundary_Z2-Membrane_boundary_Z1)/H_M, Connecting_len, H_M);
    face_back.genNodes();
    face_back.shiftNodes(x_shift, y_shift, z_shift);
    face_back.genTrianglesAndConnectors(y_len/H_S);
    face_back.genOutput(num_nodes, num_triangles, node_out, triangle_out, true);

    NonuniformFace face_left("y",0,x_len/H_S, y_len/H_S,zDimNum_T, zDimNum_B, (Membrane_boundary_Z2-Membrane_boundary_Z1)/H_M, Connecting_len, H_M);
    face_left.genNodes();
    face_left.shiftNodes(x_shift, y_shift, z_shift);
    face_left.genTrianglesAndConnectors(x_len/H_S);
    face_left.genOutput(num_nodes, num_triangles, node_out, triangle_out, true);

    NonuniformFace face_right("y",y_len,x_len/H_S, y_len/H_S, zDimNum_T, zDimNum_B, (Membrane_boundary_Z2-Membrane_boundary_Z1)/H_M, Connecting_len, H_M);
    face_right.genNodes();
    face_right.shiftNodes(x_shift, y_shift, z_shift);
    face_right.genTrianglesAndConnectors(x_len/H_S);
    face_right.genOutput(num_nodes, num_triangles, node_out, triangle_out, true);
}



void add_grids_nodes(ifstream& protein_off_file, ostringstream& node_out, int& num_nodes,
                    int L_x1, int L_x2, int L_y1, int L_y2, double H_G, int Membrane_boundary_Z1, int Membrane_boundary_Z2){
    // Add crossing nodes..
    int x_len = L_x2 - L_x1;
    int y_len = L_y2 - L_y1;
    //int H_G = H_M; //The length of grids of crossing nodes. Here assume to be the same as H_M
    int num_X = x_len / H_G;
    int num_Y = y_len / H_G;

    double vertical_min_max_Z1[num_X+1][2];
    double horizontal_min_max_Z1[num_Y+1][2];
    double vertical_min_max_Z2[num_X+1][2];
    double horizontal_min_max_Z2[num_Y+1][2];
    for (size_t i = 0; i < num_X+1 ; i++)
    {
        vertical_min_max_Z1[i][0] = std::numeric_limits<float>::max();
        vertical_min_max_Z1[i][1] = std::numeric_limits<float>::lowest();  
        vertical_min_max_Z2[i][0] = std::numeric_limits<float>::max();
        vertical_min_max_Z2[i][1] = std::numeric_limits<float>::lowest();
    }
    for (size_t i = 0; i < num_Y+1 ; i++)
    {
        horizontal_min_max_Z1[i][0] = std::numeric_limits<float>::max();
        horizontal_min_max_Z1[i][1] = std::numeric_limits<float>::lowest();
        horizontal_min_max_Z2[i][0] = std::numeric_limits<float>::max();
        horizontal_min_max_Z2[i][1] = std::numeric_limits<float>::lowest();
    }
    
    string line;
    getline(protein_off_file, line); // Reads "OFF"
    getline(protein_off_file, line); // read first line..
    istringstream s_line(line);
    int protein_num_nodes, protein_num_triangles;
    s_line >> protein_num_nodes >> protein_num_triangles;

    int Z_1l = Membrane_boundary_Z1 - H_G;
    int Z_1u = Membrane_boundary_Z1 + H_G;
    int Z_2l = Membrane_boundary_Z2 - H_G;
    int Z_2u = Membrane_boundary_Z2 + H_G;

    // Iterate through the nodes' list
    for(size_t i = 0; i < protein_num_nodes; i++){
        getline(protein_off_file, line);
        istringstream stm_line(line);
        double x, y ,z;
        stm_line >> x >> y >> z; 

        if(z >= Z_1l && z <= Z_1u){
            int x_i = (x + H_G/2 - L_x1) / H_G;
            int y_i = (y + H_G/2 - L_y1) / H_G;

            horizontal_min_max_Z1[x_i][0] = (y < horizontal_min_max_Z1[x_i][0])? y : horizontal_min_max_Z1[x_i][0];
            horizontal_min_max_Z1[x_i][1] = (y > horizontal_min_max_Z1[x_i][1])? y : horizontal_min_max_Z1[x_i][1];

            vertical_min_max_Z1[y_i][0] = (x < vertical_min_max_Z1[y_i][0])? x : vertical_min_max_Z1[y_i][0];
            vertical_min_max_Z1[y_i][1] = (x > vertical_min_max_Z1[y_i][1])? x : vertical_min_max_Z1[y_i][1];
        }

        if(z >= Z_2l && z <= Z_2u){
            int x_i = (x + H_G/2 - L_x1) / H_G;
            int y_i = (y + H_G/2 - L_y1) / H_G;


            horizontal_min_max_Z2[x_i][0] = (y < horizontal_min_max_Z2[x_i][0])? y : horizontal_min_max_Z2[x_i][0];
            horizontal_min_max_Z2[x_i][1] = (y > horizontal_min_max_Z2[x_i][1])? y : horizontal_min_max_Z2[x_i][1];

            vertical_min_max_Z2[y_i][0] = (x < vertical_min_max_Z2[y_i][0])? x : vertical_min_max_Z2[y_i][0];
            vertical_min_max_Z2[y_i][1] = (x > vertical_min_max_Z2[y_i][1])? x : vertical_min_max_Z2[y_i][1];
        }
    }


    for(int i = 0; i <= num_X; i++){
        for(int j = 0; j <= num_Y; j++){
            
            int x = L_x1 + H_G * i;
            int y = L_y1 + H_G * j;

            if(y < horizontal_min_max_Z1[i][0] || y > horizontal_min_max_Z1[i][1] || x < vertical_min_max_Z1[j][0] || x > vertical_min_max_Z1[j][1]){
                add_node(node_out, num_nodes, x, y, Membrane_boundary_Z1);
            }

            if(y < horizontal_min_max_Z2[i][0] || y > horizontal_min_max_Z2[i][1] || x < vertical_min_max_Z2[j][0] || x > vertical_min_max_Z2[j][1]){
                add_node(node_out, num_nodes, x, y, Membrane_boundary_Z2);
            }
        }
    }

    protein_off_file.clear();
    protein_off_file.seekg(0, std::ios::beg);
    // Reset file stream to the begeining


}


// MSMS NanoShaper


/**
 * @return number of triangles in the protein, this number should be added to the final output as number of facets.
*/
int add_protein_surface(ifstream& protein_off_file, ostringstream& node_out, ostringstream& triangle_out, int& num_nodes){
    // Add protein surfaces as facets, reads from the protein_off_file
    ostringstream test;
    // First, read nodes and add to the poly file..
    int num_nodes_offset = num_nodes + 1; //No. in .poly file (p) = No. in .off file (o) + num_nodes_offset
    
    string line;
    getline(protein_off_file, line); // Reads "OFF"
    getline(protein_off_file, line); // read first line..
    istringstream s_line(line);
    int protein_num_nodes, protein_num_triangles;
    s_line >> protein_num_nodes >> protein_num_triangles;

    for(size_t i = 0; i < protein_num_nodes; i++){
        getline(protein_off_file, line);
        istringstream stm_line(line);
        double x, y ,z;
        stm_line >> x >> y >> z; 
        add_node(node_out,num_nodes,x,y,z);
        //add_node(test,num_nodes,x,y,z);
    }

    for(size_t i = 0; i < protein_num_triangles; i++){
        getline(protein_off_file, line);
        istringstream stm_line(line);
        int x, n_a, n_b, n_c;
        stm_line >> x >> n_a >> n_b >> n_c;
        add_facet(triangle_out, n_a + num_nodes_offset, n_b + num_nodes_offset, n_c + num_nodes_offset); 
        //add_facet(test, n_a + num_nodes_offset, n_b + num_nodes_offset, n_c + num_nodes_offset); 
    }
    return protein_num_triangles;
}


/*
    Run tmsmesh from the .poly file at POLY_DIR and move the result into MESH_DIR
*/
void run_mesh_tmsmesh(string Tetgen_flag, string surface_poly_file_name, string surface_poly_file_name_stem){

    cout << "Running TMSMesh: " << endl;
    
    string t_cmd = (TETGEN_DIR/"tetgen").string() + ' ' + Tetgen_flag + ' ' + (POLY_DIR/surface_poly_file_name).string();
    cout << t_cmd << endl;
    int result = std::system(t_cmd.c_str());
    if(result!=0){
        cout << "Failed to run Tetgen, check configuration" << endl;
        exit(1);
    }
    // Move all files into MESH_DIR
    string surface_poly_file_name_node = surface_poly_file_name_stem + ".1.node";
    string surface_poly_file_name_ele = surface_poly_file_name_stem + ".1.ele";
    string surface_poly_file_name_face = surface_poly_file_name_stem + ".1.face";
    string surface_poly_file_name_edge = surface_poly_file_name_stem + ".1.edge";
    string surface_poly_file_name_neigh = surface_poly_file_name_stem + ".1.neigh";
    if (filesystem::exists( MESH_DIR / surface_poly_file_name_node)) filesystem::remove( MESH_DIR / surface_poly_file_name_node);
    if (filesystem::exists( MESH_DIR / surface_poly_file_name_ele)) filesystem::remove( MESH_DIR / surface_poly_file_name_ele);
    if (filesystem::exists( MESH_DIR / surface_poly_file_name_face)) filesystem::remove( MESH_DIR / surface_poly_file_name_face);
    if (filesystem::exists( MESH_DIR / surface_poly_file_name_edge)) filesystem::remove( MESH_DIR / surface_poly_file_name_edge);
    if (filesystem::exists( MESH_DIR / surface_poly_file_name_neigh)) filesystem::remove( MESH_DIR / surface_poly_file_name_neigh);
    std::filesystem::rename(POLY_DIR / surface_poly_file_name_node, MESH_DIR / surface_poly_file_name_node);
    std::filesystem::rename(POLY_DIR / surface_poly_file_name_ele, MESH_DIR / surface_poly_file_name_ele);
    std::filesystem::rename(POLY_DIR / surface_poly_file_name_face, MESH_DIR / surface_poly_file_name_face);
    std::filesystem::rename(POLY_DIR / surface_poly_file_name_edge, MESH_DIR / surface_poly_file_name_edge);
    std::filesystem::rename(POLY_DIR / surface_poly_file_name_neigh, MESH_DIR / surface_poly_file_name_neigh);

    cout << "Mesh Generated by TMSMesh and moved to: " << MESH_DIR.string() << endl;

}



void write_poly_file(std::ofstream& surface_poly_file, ostringstream& node_out, ostringstream& triangle_out ,int num_nodes, int protein_num_triangles){
    // Write .poly file
    surface_poly_file << "# .POLY File..." << endl;
    surface_poly_file << num_nodes << '\t' << 3 << '\t' << 0 << '\t' << 0 << endl;
    surface_poly_file << node_out.str();
    surface_poly_file << 6 + protein_num_triangles << '\t' << 0 << endl;
    surface_poly_file << triangle_out.str();
    surface_poly_file << 0 << endl;
    surface_poly_file << 0 << endl;
}



void print_run_time(std::clock_t  start, string task){
    double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    cout << "(" << task << ") finished in run time: " << duration << endl;
}

/**
 * Args:
 *  1 - protein's file (.)
 * 
*/
int main(int argc, char const *argv[])
{   

    // Config..
    string Tetgen_flag = "-pqAnQa1";

    // Argument.. e.x. ....
    int ETA_X = 20;
    int ETA_Y = 20;
    int ETA_Z = 20;
    int H_M = 1; // small mesh size (middle..) 
    int H_S = 2 * H_M; // large mesh size (should be some multiple of H_M)

    int Membrane_boundary_Z = 7; // This value (Mb) s.t. (Mb|H_M)  The Z coord of membrane should be (-Mb,Mb)
    // Membrane_boundary Z1 < 0, Membrane Boundary Z2 > 0
    int Membrane_boundary_Z1 = -11;
    int Membrane_boundary_Z2 = 6;
    int Connecting_len = 2; // (C) s.t. (C|H_M) The length of a connecting region
    //int Margin_len; // (M) s.t. (M|H_S) Will calculate later. This is measured bt H_S as unit
    // s.t. (z_len = 2 * (M+C+Mb))

    double tms_d = 0.4;
    //double tms_c = 0.4;
    double tms_e = 0.9;
    // 0.4 0.9 0.9

    string str_tms_d = double_to_string(tms_d);
    string str_tms_e = double_to_string(tms_e);


    // Parse Command line;
    if(argc != 2){
        cout << "Invalid input\n";
        exit(1);
    }
    string protein_pqr_file_name = argv[1];
    std::size_t last_dot_pos = protein_pqr_file_name.rfind('.');
    string file_stem = protein_pqr_file_name.substr(0, last_dot_pos);
    string protein_off_file_name =  protein_pqr_file_name + '-' + str_tms_d + '_' + str_tms_e + ".off_modify.off";

   
    // Generate or Retrieve the .off file of the protein
    ifstream protein_off_file = get_off_file(protein_pqr_file_name, str_tms_d, str_tms_e);

    // Determine the boundary of the box domain through reading the protein file
    int L_x1, L_x2, L_y1, L_y2, L_z1, L_z2;
    determine_boundary(protein_off_file, ETA_X, ETA_Y, ETA_Z,H_M, H_S,Membrane_boundary_Z1 , Membrane_boundary_Z2, Connecting_len, L_x1, L_x2, L_y1, L_y2, L_z1, L_z2);
    //protein_off_file.close();
 
    ostringstream node_out;
    ostringstream triangle_out;
    int num_nodes = 0;
    int num_triangles = 0;

    int x_len = L_x2 - L_x1;
    int y_len = L_y2 - L_y1;
    int z_len = L_z2 - L_z1;
    int x_shift = L_x1;
    int y_shift = L_y1;
    int z_shift = L_z1;    

    string surface_poly_file_name = "surface-" + file_stem + ".poly";
    string surface_poly_file_name_stem = "surface-" + file_stem;

    
    //SKIP FOR DEBUG
    
    // Genereate surface mesh of the box domain in the .POLY form..
    generate_surface_mesh(node_out, triangle_out, num_nodes, num_triangles, Membrane_boundary_Z1, Membrane_boundary_Z2, Connecting_len, H_M, H_S, L_x1, L_x2, L_y1, L_y2, L_z1, L_z2);
    // Add grids nodes to the Membrane surface...
    // H_Test = 0.6
    double H_Test = 0.6;
    add_grids_nodes(protein_off_file,node_out,num_nodes,L_x1, L_x2, L_y1, L_y2, H_M, Membrane_boundary_Z1, Membrane_boundary_Z2);
    // Add protein surfaces as facets.
    int protein_num_triangles = add_protein_surface(protein_off_file,node_out,triangle_out, num_nodes);
    // Write the .poly file


    std::ofstream surface_poly_file(POLY_DIR / surface_poly_file_name);
    if (!surface_poly_file) {
        std::cerr << "Unable to write to the .poly file";
        exit(1);
    }
    write_poly_file(surface_poly_file, node_out, triangle_out, num_nodes, protein_num_triangles);
    cout << ".poly file of the box surface saved to " << (POLY_DIR / surface_poly_file_name) << endl;

    // NEXT - GENERATE Mesh.
    run_mesh_tmsmesh(Tetgen_flag, surface_poly_file_name, surface_poly_file_name_stem);
    
    

    

    // Open 3 mesh files...
    ifstream mesh_nodes_file(MESH_DIR / (surface_poly_file_name_stem + ".1.node"));
    if(mesh_nodes_file.fail()){
        cout << "Failed to open the file: " << MESH_DIR / (surface_poly_file_name_stem + ".1.node") << endl;
        exit(1);
    }
    ifstream mesh_ele_file(MESH_DIR / (surface_poly_file_name_stem + ".1.ele"));
    if(mesh_ele_file.fail()){
        cout << "Failed to open the file: " << MESH_DIR / (surface_poly_file_name_stem + ".1.ele") << endl;
        exit(1);
    }
       ifstream mesh_neigh_file(MESH_DIR / (surface_poly_file_name_stem + ".1.neigh"));
    if(mesh_neigh_file.fail()){
        cout << "Failed to open the file: " << MESH_DIR / (surface_poly_file_name_stem + ".1.neigh") << endl;
        exit(1);
    }

    cout << "DEBUG - A" << endl;

    // NEXT - Locate Central Region
    int central_X_min = L_x1 + ETA_X - 2*H_S;
    int central_X_max = L_x2 - ETA_X + 2*H_S;
    int central_Y_min = L_y1 + ETA_Y - 2*H_S;
    int central_Y_max = L_y2 - ETA_Y + 2*H_S;

    
    std::clock_t bfs_start_time = std::clock();

     cout << "DEBUG - B" << endl;

    // Read Nodes..
    string line;
    getline(mesh_nodes_file, line); // Reads First line (n 3 0 0)
    istringstream s_line(line);
    int num_mesh_nodes;
    s_line >> num_mesh_nodes;
    double x,y,z;
    int n;
    for (size_t i = 0; i < num_mesh_nodes; i++)
    {
        getline(mesh_nodes_file, line);
        if (line[0] == '#'){
            //Skip comments lines.
            i--;
            continue;
        }
        istringstream s_line(line);
        s_line >> n >> x >> y >> z;

        if(is_in_range(z,Membrane_boundary_Z1, Membrane_boundary_Z2)){
            middle_nodes[n] = MeshNode(x,y,z);
        }
        if(is_in_range(x,central_X_min,central_X_max) && is_in_range(y,central_Y_min,central_Y_max) && is_in_range(z,Membrane_boundary_Z1,Membrane_boundary_Z2)){
            central_nodes[n] = MeshNode(x,y,z);
        }
        if(z == Membrane_boundary_Z1 || z == Membrane_boundary_Z2){
            memsurface_nodes[n] = MeshNode(x,y,z);
        }
    }
    // Read Tetragons (.ele and .neigh at the same time..)
    getline(mesh_ele_file,line);
    getline(mesh_neigh_file,line);
    s_line.str(line);
    int num_mesh_triangles;
    s_line >> num_mesh_triangles;
    // Note: all_tetragon's index also start from 1, leave the first all_tetragon[0] as blank
    MeshTetragon all_tetragon[num_mesh_triangles+1];


     cout << "DEBUG - C" << endl;


    int n_node_a, n_node_b, n_node_c, n_node_d, label;
    int n_neigh_a, n_neigh_b, n_neigh_c, n_neigh_d;
    for (size_t i = 0; i < num_mesh_triangles; i++)
    {
        // n = i+1
        getline(mesh_ele_file, line);
        istringstream s_line(line);
        s_line >> n >> n_node_a >> n_node_b >> n_node_c >> n_node_d >> label;
        MeshTetragon new_tet(n_node_a, n_node_b, n_node_c, n_node_d, label);
        all_tetragon[n] = new_tet;
        getline(mesh_neigh_file, line);
        // Update Membrane.. (Label 3)
        if(tet_is_middle(new_tet) && new_tet.label != PROTEIN_LABEL){
            all_tetragon[n].label = MEMBRANE_LABEL;
        }        
        if(tet_is_central(new_tet) && new_tet.label != PROTEIN_LABEL){
            istringstream s_line(line);
            s_line >> n >> n_neigh_a >> n_neigh_b >> n_neigh_c >> n_neigh_d;
            central_tetragon[n] = MeshTetragonSearch(n_node_a,n_node_b,n_node_c,n_node_d,n_neigh_a,n_neigh_b,n_neigh_c,n_neigh_d,label);
        }
    }

    cout << "(BFS Search) A Total of " << central_tetragon.size() << " Tetragons in serach area" << endl;



    // NEXT - Start Search
    // TRIAL1: BFS
    // For testing
    std::vector<MeshTetragonSearch> triangles_0;
    double x_min_0, x_max_0, y_min_0, y_max_0;

    std::queue<int> queue;
    // Find the first start triangle (not protein) and let it be labled 0 and start first serach...
    for (auto it = central_tetragon.begin(); it != central_tetragon.end(); ++it) {
        if(it->second.label == -1){
            it->second.label = 0;
            it->second.visited = true;
            queue.push(it->first);
            break;
        }   
    }
    bool label_0_touch_boundary = false;


    // Start BFS
    while(!queue.empty()){
        int n_tretragon = queue.front();
        queue.pop();
        MeshTetragonSearch curr_triangle = central_tetragon[n_tretragon];
        if(real_neighbor(curr_triangle.neigh1, label_0_touch_boundary)){
            visit_triangle(curr_triangle.neigh1,curr_triangle.label);
            queue.push(curr_triangle.neigh1);
        }
        if(real_neighbor(curr_triangle.neigh2, label_0_touch_boundary)){
            visit_triangle(curr_triangle.neigh2,curr_triangle.label);
            queue.push(curr_triangle.neigh2);
        }
        if(real_neighbor(curr_triangle.neigh3, label_0_touch_boundary)){
            visit_triangle(curr_triangle.neigh3,curr_triangle.label);
            queue.push(curr_triangle.neigh3);
        }
        if(real_neighbor(curr_triangle.neigh4, label_0_touch_boundary)){
            visit_triangle(curr_triangle.neigh4,curr_triangle.label);
            queue.push(curr_triangle.neigh4);
        }
    }
    int count = 0;
    if(label_0_touch_boundary){
        // Label 0 is membrane, Label -1 is solvent, change all tet with label -1 into solvent.
       for (std::unordered_map<int, MeshTetragonSearch>::iterator it = central_tetragon.begin(); it != central_tetragon.end(); it++) {
            if(it->second.label != 0){
                all_tetragon[it->first].label = SOLVENT_LABEL;
                count++;
            }
        }
    }else{
        // Label 0 is solvent, Label -1 is membrane, change all tet with label 0 into solvent
        for (std::unordered_map<int, MeshTetragonSearch>::iterator it = central_tetragon.begin(); it != central_tetragon.end(); it++) {
            if(it->second.label == 0){
                all_tetragon[it->first].label = SOLVENT_LABEL;
                count++;
            }
        }  
    }
    print_run_time(bfs_start_time, "BFS Search");
    cout << "(BFS Search) Found " << count << " Tetragon as solvent..." << endl;


    // (REWRITE) the .ele file based on the updated label...
    ostringstream new_ele_file_s;
    new_ele_file_s << num_mesh_triangles << '\t' << 4 << '\t' << 1 << endl;
    int all_tetragon_size = sizeof(all_tetragon) / sizeof( all_tetragon[0]);
    for(int i = 1; i < all_tetragon_size; i++){
        all_tetragon[i].print_tet(new_ele_file_s,i);
    }

    ofstream mesh_ele_file_rewrite(MESH_DIR / (surface_poly_file_name_stem + ".1.ele"));
    mesh_ele_file_rewrite << new_ele_file_s.str();
    

    /*
    // DEBUG A - {PRINT CENTRAL NODES}
    for(auto it = central_nodes.begin(); it != central_nodes.end(); it++){
        //cout << it->first << '\t' << it->second.x <<  '\t' << it->second.y <<  '\t' << it->second.z <<  endl;
    }
    // DEBUG B - {PRINT CEMTRAL TET}
    for(auto it = central_tetragon.begin(); it != central_tetragon.end(); it++){
        //cout << it->first << '\t' << it->second.node1 <<  '\t' << it->second.node2 <<  '\t' << it->second.node3 <<   '\t' << it->second.node4 << endl;
    }
    */

    // surface_poly_file.close();
    protein_off_file.close();

    return 0;
}



