// TODO: Better comments
// TODO: Make more input based
// TODO: Make more function based
// TODO: Experiment with splines instead of lines: BSpline ( expression ) = { expression-list }; and Spline ( expression ) = { expression-list };
// TODO: Need to better understand transfinite definitions

// ! THIS VERSION DOES NOT PRODUCE WHAT WE NEED. IT IS ELEGANT IN THAT IT ATTEMPTS TO DEFINE EVERYTHING FROM SCRATCH BUT DOES NOT PRODUCE WHAT WE WANT

// ! GMSH Needs to export the 3D mesh in Version2 ASCII Format

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
    int extrude_dist {1};

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

        // Check for identical first and last coordinates, modify accordingly
        if(airfoil_pts.front()==airfoil_pts.back()){
            airfoil_pts.pop_back();
        }

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

        // Write Side Two Box Lines - Clockwise
        for(auto pt_num:box_side_two_points){
            if (pt_num==box_side_two_points.back()){
                out_file << "Line(" << (line_ct) << ") = {" << pt_num << ", " << box_side_two_points.front() << "};" << std::endl;
            }else {
                out_file << "Line(" << (line_ct) << ") = {" << pt_num << ", " << pt_num+1 << "};" << std::endl;
            }
            box_side_two_lines.push_back(line_ct);
            line_ct++; 
        }

        // Write Transverse Box Lines - Clockwise
        for (long long unsigned int i = 0; i<box_side_one_points.size(); i++){
            out_file << "Line(" << (line_ct) << ") = {" << box_side_one_points.at(i) << ", " << box_side_two_points.at(i) << "};" << std::endl;
            box_tranv_lines.push_back(line_ct);
            line_ct++; 
        }

        // DEFINE CURVE LOOPS
        int curve_ct {1};

        int af_side_one_curve{};
        int af_side_two_curve{};
        std::vector<int> af_tranv_curves{};

        int box_side_one_curve{};
        int box_side_two_curve{};
        std::vector<int> box_tranv_curves{};

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

        // Airfoil Transverse Curve Loops - RH Rule pointing away from faces
        for(long long unsigned int i = 0; i<af_tranv_lines.size(); i++){
            out_file << "Curve Loop(" << (curve_ct) << ") = {";
            if (af_tranv_lines.at(i)==af_tranv_lines.back()){
                out_file << -1*af_side_one_lines.at(i) << ", "
                         << af_side_two_lines.at(i) << ", "
                         << -1*af_tranv_lines.front() << ", "
                         << af_tranv_lines.at(i) << "};" << std::endl;
            }
            else{
                out_file << -1*af_side_one_lines.at(i) << ", "
                         << af_side_two_lines.at(i) << ", "
                         << af_tranv_lines.at(i) << ", "
                         << -1*af_tranv_lines.at(i+1) << "};" << std::endl;
            }
            af_tranv_curves.push_back(curve_ct);
            curve_ct++;
        }

        // Box Side One Curve Loop - Face
        out_file << "Curve Loop(" << (curve_ct) << ") = {";
        for (int ln_num = box_side_one_lines.back(); ln_num>=box_side_one_lines.front(); ln_num--){
            if(ln_num==box_side_one_lines.front()){
                out_file << -ln_num << "};" << std::endl;
            }else{
                out_file << -ln_num << ", ";
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

        // Box Transverse Curve Loops
        for (long long unsigned int i = 0; i<box_tranv_lines.size(); i++){
            out_file << "Curve Loop(" << (curve_ct) << ") = {";
            if(box_tranv_lines.at(i)==box_tranv_lines.back()){
                out_file 
                         << box_side_one_lines.at(i) << ", "
                         << -1*box_side_two_lines.at(i) << ", "
                         << box_tranv_lines.front() << ", "
                         << -1*box_tranv_lines.at(i) << "};" << std::endl;
            }else{
                out_file 
                         << box_side_one_lines.at(i) << ", "
                         << -1*box_side_two_lines.at(i) << ", "
                         << -1*box_tranv_lines.at(i) << ", "
                         << box_tranv_lines.at(i+1) << "};" << std::endl;
            }
            box_tranv_curves.push_back(curve_ct);
            curve_ct++;
        }

        // DEFINE PLANE SURFACES
        int surface_ct {1};

        std::vector<int> af_tranv_psurfaces{};

        int side_one_psurface{};
        int side_two_psurface{};
        int bottom_psurface{};
        int outlet_psurface{};
        int top_psurface{};
        int inlet_psurface{};

        // Side One Plane Surface
        out_file << "Plane Surface(" << (surface_ct) << ") = {" 
                 << (box_side_one_curve) << ", " << (af_side_one_curve) << "};" << std::endl; 
        side_one_psurface = surface_ct;
        surface_ct++;

        // Side Two Plane Surface
        out_file << "Plane Surface(" << (surface_ct) << ") = {" 
                 << (box_side_two_curve) << ", " << (af_side_two_curve) << "};" << std::endl; 
        side_two_psurface = surface_ct;
        surface_ct++;    

        // Airfoil Plane Surfaces
        for (auto tranv_loop:af_tranv_curves){
            out_file << "Plane Surface(" << (surface_ct) << ") = {" 
                     << (tranv_loop) << "};" << std::endl; 
            af_tranv_psurfaces.push_back(surface_ct);
            surface_ct++;  
        }

        // Bottom Plane Surface
        out_file << "Plane Surface(" << (surface_ct) << ") = {" 
                 << (box_tranv_curves.at(0)) << "};" << std::endl; 
        bottom_psurface = surface_ct;
        surface_ct++;

        // Outlet Plane Surface
        out_file << "Plane Surface(" << (surface_ct) << ") = {" 
                 << (box_tranv_curves.at(1)) << "};" << std::endl; 
        outlet_psurface = surface_ct;
        surface_ct++;

        // Top Plane Surface
        out_file << "Plane Surface(" << (surface_ct) << ") = {" 
                 << (box_tranv_curves.at(2)) << "};" << std::endl; 
        top_psurface = surface_ct;
        surface_ct++;

        // Inlet Plane Surface
        out_file << "Plane Surface(" << (surface_ct) << ") = {" 
                 << (box_tranv_curves.at(3)) << "};" << std::endl; 
        inlet_psurface = surface_ct;
        surface_ct++;

        // DEFINE SURFACE LOOPS - Surface Loop ( expression ) = { expression-list } < Using Sewing >;
        int surface_loop_ct {1};

        int box_surface_loop{};
        int af_surface_loop{};

        //Box Surface Loop
        out_file << "Surface Loop (" << surface_loop_ct << ") = {"
                 << inlet_psurface << ", "
                 << side_one_psurface << ", "
                 << outlet_psurface << ", "
                 << side_two_psurface << ", "
                 << bottom_psurface << ", "
                 << top_psurface << "};" << std::endl;
        box_surface_loop = surface_loop_ct;
        surface_loop_ct++;

        // Airfoil Surface Loop
        out_file << "Surface Loop (" << surface_loop_ct << ") = {";
        for (auto af_surf:af_tranv_psurfaces){
            if(af_surf == af_tranv_psurfaces.back()){
                out_file << af_surf << "};" << std::endl;
            }else{
                out_file << af_surf << ", ";
            }
        }
        af_surface_loop = surface_loop_ct;
        surface_loop_ct++;

        // DEFINE VOLUMES - Volume ( expression ) = { expression-list };
        int volume_ct {1};
        
        int mesh_vol {};

        out_file << "Volume (" << volume_ct << ") = {"
                 << box_surface_loop << ", "
                 << af_surface_loop << "};" << std::endl;

        mesh_vol = volume_ct;
        volume_ct++; 

        // DEFINE PHYSICAL SURFACES
        out_file << "Physical Surface(\"Inlet\") = {" << inlet_psurface << "};" << std::endl;
        out_file << "Physical Surface(\"Bottom\") = {" << bottom_psurface << "};" << std:: endl;
        out_file << "Physical Surface(\"Outlet\") = {" << outlet_psurface << "};" << std::endl;
        out_file << "Physical Surface(\"Top\") = {" << top_psurface << "};" << std::endl;
        out_file << "Physical Surface(\"Sides\") = {" << side_one_psurface << ", " << side_two_psurface << "};" << std::endl;
        out_file << "Physical Surface(\"Foil\") = {";
        for (auto af_surf:af_tranv_psurfaces){
            if (af_surf == af_tranv_psurfaces.back()){
                out_file << af_surf << "};" << std::endl;
            }else{
                out_file << af_surf << ", ";
            }
        }

        // DEFINE PHYSICAL VOLUMES
        out_file << "Physical Volume(\"Vol\") = {" << mesh_vol << "};" << std::endl;  

        // DEFINE TRANSFINITE LINES
        out_file << "Transfinite Curve {"; 
                for (auto tranv_line:af_tranv_lines){
                    out_file << tranv_line << ", ";
                }
                for (auto tranv_line:box_tranv_lines){
                    if(tranv_line == box_tranv_lines.back()){
                        out_file << tranv_line << "} = 1 Using Progression 1;" << std::endl;
                    }else{
                        out_file << tranv_line << ", ";
                    }
                }

        // RECOMBINE EXTERNAL SURFACES 
        out_file << "Recombine Surface {" << side_one_psurface << "};" << std::endl;
        out_file << "Recombine Surface {" << side_two_psurface << "};" << std::endl;

        //? DEFINE BOUNDARY LAYER


    }

    // CLOSE FILES
    in_file.close();
    out_file.close();

    return 0;
}