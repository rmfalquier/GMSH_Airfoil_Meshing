// TODO: Better comments
// TODO: Make more input based
// TODO: Functions
// TODO: Make more general

// GMSH Needs to export the 3D mesh in Version2 ASCII Format

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

int main(){
    // AIRFOILE FILE NAME
    std::string in_airfoil_file_name{"./.AIRFOILS_IN/NACA2412.txt"};

    // STORAGE CONTAINERS
    std::string line{};
    std::vector<std::vector<double>> airfoil_pts;

    // OPEN AND READ AIRFOIL FILE TO SAVE DATA IN THE AIRFOIL POINT VECTOR
    std::ifstream in_file;
    in_file.open(in_airfoil_file_name);
    if (!in_file) {
        std::cerr << "Problem opening airfoil file" << std::endl;
        return 1;
    }else{
        while (std::getline(in_file, line)) {
            std::istringstream test_stream {line};
            double test_dbl;
            if(test_stream >> test_dbl){
                // CREATE TEMPORARY STORAGE FOR COORDINATES
                double temp_x_coord {};
                double temp_y_coord {};
                std::istringstream coord_stream {line};
                std::vector<double> pt_xy_coords;

                // READ COORDINATES FROM STREAM, PUSH TO VECTOR IN DESIRED ORDER, ADD TO AIRFOIL VECTOR
                coord_stream >> temp_x_coord >> temp_y_coord;
                pt_xy_coords.push_back(temp_x_coord);
                pt_xy_coords.push_back(temp_y_coord);
                airfoil_pts.push_back(pt_xy_coords);
            }
        }
    }

    // CREATE AND WRITE OUTPUT FILE
    std::string out_airfoil_file_name {"./.AIRFOILS_OUT/NACA2412.geo"};
    std::ofstream out_file{out_airfoil_file_name};
    if (!out_file) {
        std::cerr << "Error opening output file" << std::endl;
        return 1;
    }
    else{
        // WRITE AIRFOIL POINTS
        int point_number {0}; //GMSH Formatting Starts at 1
        for(auto xy_pt:airfoil_pts){
            ++point_number;
            out_file << "Point(" << point_number << ") = {" <<
                xy_pt.front() << ", " << xy_pt.back() << 
                ", 0.0, 1.0};" << std::endl; 
        }

        // WRITE BOX POINTS - CLOCKWISE
        out_file << "Point(" << point_number+1 << ") = {-13.5, -14.0, 0.0, 1.0};" << std::endl;
        out_file << "Point(" << point_number+2 << ") = {14.5, -14.0, 0.0, 1.0};" << std::endl; 
        out_file << "Point(" << point_number+3 << ") = {14.5, 14.0, 0.0, 1.0};" << std::endl; 
        out_file << "Point(" << point_number+4 << ") = {-13.5, 14.0, 0.0, 1.0};" << std::endl;  
        point_number += 4;

        // WRITE AIRFOIL LINES
        for (int i = 1; i <= (point_number-4); i++){
            if(i==(point_number-4)){
                // out_file << "Line(" << i << ") = {" << i << ", " << 1 << "};" << std::endl; 
                out_file << "Line(" << (point_number+i) << ") = {" << i << ", " << 1 << "};" << std::endl; 
            }
            else{
                // out_file << "Line(" << i << ") = {" << i << ", " << i+1 << "};" << std::endl; 
                out_file << "Line(" << (point_number+i) << ") = {" << i << ", " << i+1 << "};" << std::endl; 
            }
        }

        // WRITE BOX LINES
        for (int i = (point_number-3); i <= point_number; i++){
            if(i==point_number){
                // out_file << "Line(" << i << ") = {" << i << ", " << (point_number-3) << "};" << std::endl; 
                out_file << "Line(" << (point_number+i) << ") = {" << i << ", " << (point_number-3) << "};" << std::endl; 
            }
            else{
                // out_file << "Line(" << i << ") = {" << i << ", " << i+1 << "};" << std::endl; 
                out_file << "Line(" << (point_number+i) << ") = {" << i << ", " << i+1 << "};" << std::endl; 
            }
        }
        
        // UPDATE POINT NUMBER
        point_number *= 2;

        // WRITE CURVE LOOPS
        out_file << "Curve Loop(" << (point_number+1) << ") = {";
        for (int i = (point_number-3); i <= point_number; i++){
            if(i==point_number){
                out_file << i << "};" << std::endl;
            }else{
                out_file << i << ", ";
            }
        }

        out_file << "Curve Loop(" << (point_number+2) << ") = {";
        for (int i = ((point_number/2)+1); i <= (point_number-4); i++){
            if(i==(point_number-4)){
                out_file << i << "};" << std::endl;
            }else{
                out_file << i << ", ";
            }
        }

        // WRITE PLANE SURFACE
        out_file << "Plane Surface(" << (point_number+3) << ") = {" 
            << (point_number+1) << ", " << (point_number+2) << "};" << std::endl; 

        // WRITE EXTRUDE
        out_file 
         << "Extrude {0.0, 0.0, 2.0} {" << std::endl 
         << "   Surface{" << (point_number+3) << "};" << std::endl 
         << "   Layers{1};" << std::endl 
         << "   Recombine;" << std::endl //Recombine is for quadrilateral versus triangular mesh
         << "}" << std::endl; 

        // PHYSICAL GROUPS 
        int foil_pts {(point_number-8)/2};
        // Per above definition
        int right {point_number+3};

        // First Surface after last line is defined is the bottom of the extrude:
        // (Jump from last line + remaining extruded face lines + number of airfoil lines + 
        // number of bottom conecting vertices + 2 for indexing)
        int bottom {point_number + 5 + 3 + foil_pts + 2 + 2}; 

        // Next Surfaces clockwise around the extrude add 4 to index
        int outlet {bottom + 4};
        int top {outlet + 4};
        int inlet {top + 4};

        // Right: 
        // One face per airfoil point where each face adds 4 to index + 1 for last index
        int left {inlet + (foil_pts*4) + 1};

        // Writing:
        out_file << "Physical Surface(\"Sides\") = {" << right << ", " << left << "};" << std::endl;
        out_file << "Physical Surface(\"Bottom\") = {" << bottom << "};" << std:: endl;
        out_file << "Physical Surface(\"Outlet\") = {" << outlet << "};" << std::endl;
        out_file << "Physical Surface(\"Top\") = {" << top << "};" << std::endl;
        out_file << "Physical Surface(\"Inlet\") = {" << inlet << "};" << std::endl;
        out_file << "Physical Surface(\"Foil\") = {";
        for (int i = 1; i<=foil_pts; i++){
            if (i < foil_pts){
                out_file << (inlet + (i*4)) << ", ";
            }else{
                out_file << (inlet + (i*4)) << "};" << std::endl;
            }
        }

        // VOLUME
        out_file << "Physical Volume(\"Vol\") = {1};" << std::endl;  

        /*
        // TRANSFINITE SURFACE: RIGHT
        out_file << "Transfinite Surface {" << right << "} = {";
        for (int i = 1; i <= (foil_pts+4); i++){
            if (i<(foil_pts+4)){
                out_file << i << ", ";
            }else {
                out_file << i << "};" << std::endl;
            }
        }

        // TRANSFINITE SURFACE: LEFT
        out_file << "Transfinite Surface {" << left << "} = {"
            << ((point_number/2)+1) << ", " << ((point_number/2)+2) << ", " 
            << ((point_number/2)+6) << ", " << ((point_number/2)+10) << ", ";
        int left_clockwise_start {((point_number/2)+18)+((foil_pts-2)*4)};
        for (int i = left_clockwise_start; i > ((point_number/2)+18); i-=4){
                out_file << i << ", ";
        }
        out_file << ((point_number/2)+18) << ", " << ((point_number/2)+17) << "};" << std::endl;
        */

        // RECOMBINE SURFACES
        out_file << "Recombine Surface {" << left << "};" << std::endl;
        out_file << "Recombine Surface {" << right << "};" << std::endl;
    }

    // CLOSE FILES
    in_file.close();
    out_file.close();

    return 0;
}