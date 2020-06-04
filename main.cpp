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
        // WRITE POINTS
        int point_number {0}; //GMSH Formatting Starts at 1
        for(auto xy_pt:airfoil_pts){
            ++point_number;
            out_file << "Point(" << point_number << ") = {" <<
                xy_pt.front() << ", " << xy_pt.back() << 
                ", 0.0, 1.0};" << std::endl; 
        }

        // WRITE LINES
        for (int i = 1; i <= point_number; i++){
            if(i==point_number){
                out_file << "Line(" << i << ") = {" << i << ", " << 1 << "};" << std::endl; 
            }
            else{
                out_file << "Line(" << i << ") = {" << i << ", " << i+1 << "};" << std::endl; 
            }
        }
    }

    // ! LEFT OFF AT GMSH PART4

    // CLOSE FILES
    in_file.close();
    out_file.close();

    return 0;
}