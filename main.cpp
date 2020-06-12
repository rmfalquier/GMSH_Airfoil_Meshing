// TODO: Better comments
// TODO: Make more input based
// TODO: Make more function based
// TODO: Experiment with splines instead of lines: BSpline ( expression ) = { expression-list }; and Spline ( expression ) = { expression-list };
// TODO: Need to better understand transfinite definitions

//! Left off at the definition of transverse curve loops for surface definitions, remember to use RH rule pointing away from surfaces and into holes

//! GMSH Needs to export the 3D mesh in Version2 ASCII Format

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

int main(){
    // AIRFOILE FILE NAME
    // std::string in_airfoil_file_name{"./.AIRFOILS_IN/NACA2412.txt"};
    // std::string out_airfoil_file_name {"./.AIRFOILS_OUT/NACA2412.geo"};
    std::string in_airfoil_file_name{"./.AIRFOILS_IN/riot_ctr_prof.dat"};
    std::string out_airfoil_file_name {"./.AIRFOILS_OUT/RIOT_CTR_PROF.geo"};

    // MESH PARAMETERS
    int extrude_dist {2};

    // STORAGE CONTAINERS
    std::string line{};
    std::vector<std::vector<double>> airfoil_pts;

    // FILE PARAMETERS
    std::ifstream in_file;
    in_file.open(in_airfoil_file_name);
    std::ofstream out_file{out_airfoil_file_name};

    // OPEN AND READ AIRFOIL FILE TO SAVE DATA IN THE AIRFOIL POINT VECTOR
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
    if (!out_file) {
        std::cerr << "Error opening output file" << std::endl;
        return 1;
    }
    else{
        // DEFINE POINTS
        int point_ct {1}; // GMSH Formatting Starts at 1

        std::vector<int> af_side_one_points{};
        std::vector<int> af_side_two_points{};
        std::vector<int> box_side_one_points{};
        std::vector<int> box_side_two_points{};

        // Write Airfoil Points - Clockwise
        for(auto xy_pt:airfoil_pts){
            out_file << "Point(" << point_ct << ") = {" <<
                xy_pt.front() << ", " << xy_pt.back() << 
                ", 0.0, 1.0};" << std::endl; 

            out_file << "Point(" << (point_ct + airfoil_pts.size()) << ") = {" 
                << xy_pt.front() << ", " << xy_pt.back() << ", " 
                << extrude_dist << ", 1.0};" << std::endl; 

            af_side_one_points.push_back(point_ct);    
            af_side_two_points.push_back(point_ct + airfoil_pts.size());
            point_ct++;
        }
        point_ct += airfoil_pts.size();

        // Write Side One Box Points - Clockwise
        out_file << "Point(" << point_ct << ") = {-13.5, -14.0, 0.0, 1.0};" << std::endl;
        box_side_one_points.push_back(point_ct);
        point_ct++;
        out_file << "Point(" << point_ct << ") = {14.5, -14.0, 0.0, 1.0};" << std::endl; 
        box_side_one_points.push_back(point_ct);
        point_ct++;
        out_file << "Point(" << point_ct << ") = {14.5, 14.0, 0.0, 1.0};" << std::endl; 
        box_side_one_points.push_back(point_ct);
        point_ct++;
        out_file << "Point(" << point_ct << ") = {-13.5, 14.0, 0.0, 1.0};" << std::endl;  
        box_side_one_points.push_back(point_ct);
        point_ct++;

        // Write Side Two Box Points - Clockwise
        out_file << "Point(" << point_ct << ") = {-13.5, -14.0, "<< extrude_dist << ", 1.0};" << std::endl;
        box_side_two_points.push_back(point_ct);
        point_ct++;
        out_file << "Point(" << point_ct << ") = {14.5, -14.0, "<< extrude_dist << ", 1.0};" << std::endl;
        box_side_two_points.push_back(point_ct);
        point_ct++;
        out_file << "Point(" << point_ct << ") = {14.5, 14.0, "<< extrude_dist << ", 1.0};" << std::endl;
        box_side_two_points.push_back(point_ct);
        point_ct++;
        out_file << "Point(" << point_ct << ") = {-13.5, 14.0, "<< extrude_dist << ", 1.0};" << std::endl; 
        box_side_two_points.push_back(point_ct);
        point_ct++;

        // DEFINE LINES
        int line_ct {1}; // GMSH Formatting Starts at 1
        
        std::vector<int> af_side_one_lines{};
        std::vector<int> af_side_two_lines{};
        std::vector<int> af_tranv_lines{};

        std::vector<int> box_side_one_lines{};
        std::vector<int> box_side_two_lines{};
        std::vector<int> box_tranv_lines{};

        // Write Side One Airfoil Lines - Clockwise
        for(auto pt_num:af_side_one_points){
            if (pt_num==af_side_one_points.back()){
                out_file << "Line(" << (line_ct) << ") = {" << pt_num << ", " << af_side_one_points.front() << "};" << std::endl;
            }else {
                out_file << "Line(" << (line_ct) << ") = {" << pt_num << ", " << pt_num+1 << "};" << std::endl;
            }
            af_side_one_lines.push_back(line_ct);
            line_ct++; 
        }

        // Write Side Two Airfoil Lines - Clockwise
        for(auto pt_num:af_side_two_points){
            if (pt_num==af_side_two_points.back()){
                out_file << "Line(" << (line_ct) << ") = {" << pt_num << ", " << af_side_two_points.front() << "};" << std::endl;
            }else {
                out_file << "Line(" << (line_ct) << ") = {" << pt_num << ", " << pt_num+1 << "};" << std::endl;
            }
            af_side_two_lines.push_back(line_ct);
            line_ct++; 
        }

        // Write Transverse Airfoil Lines - Clockwise
        for (long long unsigned int i = 0; i<af_side_one_points.size(); i++){
            out_file << "Line(" << (line_ct) << ") = {" << af_side_one_points.at(i) << ", " << af_side_two_points.at(i) << "};" << std::endl;
            af_tranv_lines.push_back(line_ct);
            line_ct++; 
        }

        // Write Side One Box Lines - Clockwise
        for(auto pt_num:box_side_one_points){
            if (pt_num==box_side_one_points.back()){
                out_file << "Line(" << (line_ct) << ") = {" << pt_num << ", " << box_side_one_points.front() << "};" << std::endl;
            }else {
                out_file << "Line(" << (line_ct) << ") = {" << pt_num << ", " << pt_num+1 << "};" << std::endl;
            }
            box_side_one_lines.push_back(line_ct);
            line_ct++; 
        }

        // Write Side Two Airfoil Lines - Clockwise
        for(auto pt_num:box_side_two_points){
            if (pt_num==box_side_two_points.back()){
                out_file << "Line(" << (line_ct) << ") = {" << pt_num << ", " << box_side_two_points.front() << "};" << std::endl;
            }else {
                out_file << "Line(" << (line_ct) << ") = {" << pt_num << ", " << pt_num+1 << "};" << std::endl;
            }
            box_side_two_lines.push_back(line_ct);
            line_ct++; 
        }

        // Write Transverse Airfoil Lines - Clockwise
        for (long long unsigned int i = 0; i<box_side_one_points.size(); i++){
            out_file << "Line(" << (line_ct) << ") = {" << box_side_one_points.at(i) << ", " << box_side_two_points.at(i) << "};" << std::endl;
            box_tranv_lines.push_back(line_ct);
            line_ct++; 
        }

        // DEFINE CURVE LOOPS
        int curve_ct {1};

        int af_side_one_curve{};
        int af_side_two_curve{};
        // std::vector<int> af_tranv_curves{};

        int box_side_one_curve{};
        int box_side_two_curve{};
        // std::vector<int> box_tranv_curves{};

        // Airfoil Side One Curve Loop - Hole
        out_file << "Curve Loop(" << (curve_ct) << ") = {";
        for (auto ln_num:af_side_one_lines){
            if(ln_num==af_side_one_lines.back()){
                out_file << ln_num << "};" << std::endl;
            }else{
                out_file << ln_num << ", ";
            }
        }
        af_side_one_curve = curve_ct;
        curve_ct++;

        // Airfoil Side Two Curve Loop - Hole
        out_file << "Curve Loop(" << (curve_ct) << ") = {";
        for (int ln_num = af_side_two_lines.back(); ln_num>=af_side_two_lines.front(); ln_num--){
            if(ln_num==af_side_two_lines.front()){
                out_file << ln_num << "};" << std::endl;
            }else{
                out_file << ln_num << ", ";
            }
        }
        af_side_two_curve = curve_ct;
        curve_ct++;

        // Box Side One Curve Loop - Face
        out_file << "Curve Loop(" << (curve_ct) << ") = {";
        for (int ln_num = box_side_one_lines.back(); ln_num>=box_side_one_lines.front(); ln_num--){
            if(ln_num==box_side_one_lines.front()){
                out_file << ln_num << "};" << std::endl;
            }else{
                out_file << ln_num << ", ";
            }
        }
        box_side_one_curve = curve_ct;
        curve_ct++;

        // Box Side Two Curve Loop - Face
        out_file << "Curve Loop(" << (curve_ct) << ") = {";
        for (auto ln_num:box_side_two_lines){
            if(ln_num==box_side_two_lines.back()){
                out_file << ln_num << "};" << std::endl;
            }else{
                out_file << ln_num << ", ";
            }
        }
        box_side_two_curve = curve_ct;
        curve_ct++;

        // DEFINE PLANE SURFACES
        int surface_ct {1};

        out_file << "Plane Surface(" << (surface_ct) << ") = {" 
            << (box_side_one_curve) << ", " << (af_side_one_curve) << "};" << std::endl; 
        surface_ct++;

        out_file << "Plane Surface(" << (surface_ct) << ") = {" 
            << (box_side_two_curve) << ", " << (af_side_two_curve) << "};" << std::endl; 
        surface_ct++;    

        // Define Surface Loops 
        // Define Volumes
        // Define Physical Surfaces
        // Define Physical Volumes
        // Recombine External Surfaces 
        //? Boundary Layer 

        /*

        // PHYSICAL GROUPS 
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

        // RECOMBINE SURFACES
        out_file << "Recombine Surface {" << left << "};" << std::endl;
        out_file << "Recombine Surface {" << right << "};" << std::endl;
        */
    }

    // CLOSE FILES
    in_file.close();
    out_file.close();

    return 0;
}