// TODO: Make more input based
// TODO: Make more function based
// TODO: Experiment with splines instead of lines: BSpline ( expression ) = { expression-list }; and Spline ( expression ) = { expression-list };

// ! GMSH Needs to export the 3D mesh in Version2 ASCII Format

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

int main(){
    // AIRFOILE FILE NAME
    std::string in_airfoil_file_name{"./.AIRFOILS_IN/NACA2412.txt"};
    std::string out_airfoil_file_name {"./.AIRFOILS_OUT/NACA2412.geo"};
    // std::string in_airfoil_file_name{"./.AIRFOILS_IN/riot_ctr_prof.dat"};
    // std::string out_airfoil_file_name {"./.AIRFOILS_OUT/RIOT_CTR_PROF.geo"};
    // std::string in_airfoil_file_name{"./.AIRFOILS_IN/riot_ctr_braked_100.dat"};
    // std::string out_airfoil_file_name {"./.AIRFOILS_OUT/RIOT_CTR_BRAKED_PROF.geo"};

    // MESH PARAMETERS
    double square_side {30.0};
    double grid_size {square_side/20};
    double mesh_thickness {square_side/10};

    // BOUNDARY LAYER PARAMETERS
    float anisomax {270.0}; // Threshold angle for creating a mesh fan in the boundary layer
    int quads {1}; // Generate recombined elements in the boundary layer
    float hfar {0.03}; // Element size far from the wall
    float hwall_n {0.001}; // Mesh Size Normal to the The Wall
    float ratio {1.1}; // Size Ratio Between Two Successive Layers
    float thickness {0.05}; // must be bigger than hwall_n

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
        // MAKE ANY NECCESSARY MODS TO AIRFOIL POINTS AND FLAG THE RIGHT CONDITIONS
        // Check for identical first and last airfoil coordinates, modify accordingly
        if(airfoil_pts.front()==airfoil_pts.back()){
            airfoil_pts.pop_back();
        }

        // Check for odd number of lines
        std::vector<double> xy_pt_last {airfoil_pts.back()};
        std::vector<double> xy_pt_first {airfoil_pts.front()};
        if((airfoil_pts.size()%2)!=0){
            std::vector<double> new_pt{};

            double new_x {xy_pt_last.front() + ((xy_pt_first.front() - xy_pt_last.front())/2)};
            double new_y {xy_pt_last.back() + ((xy_pt_first.back() - xy_pt_last.back())/2)};
            new_pt.push_back(new_x);
            new_pt.push_back(new_y);
            airfoil_pts.push_back(new_pt);
        }

        // Check for fan conditions
        bool double_fan {false};
        if(xy_pt_first.front() == xy_pt_last.front()){
            double_fan =true;
        }     
        
        // DEFINE POINTS
        int point_ct {1}; // GMSH Formatting Starts at 1

        std::vector<int> box_side_points{};
        std::vector<int> af_side_points{};

        // Write Box Points - Counterclockwise from bottom left corner
        out_file << "Point(" << point_ct << ") = {"
                 << -((square_side/2)-0.5) << ", " << -(square_side/2) << ", 0.0, " 
                 << grid_size << "};" << std::endl;
        box_side_points.push_back(point_ct);
        point_ct++;

        out_file << "Point(" << point_ct << ") = {"
                 << -((square_side/2)-0.5) << ", " << (square_side/2) << ", 0.0, " 
                 << grid_size << "};" << std::endl;
        box_side_points.push_back(point_ct);
        point_ct++;

        out_file << "Point(" << point_ct << ") = {"
                 << ((square_side/2)+0.5) << ", " << (square_side/2) << ", 0.0, " 
                 << grid_size << "};" << std::endl;
        box_side_points.push_back(point_ct);
        point_ct++;

        out_file << "Point(" << point_ct << ") = {"
                 << ((square_side/2)+0.5) << ", " << -(square_side/2) << ", 0.0, " 
                 << grid_size << "};" << std::endl;
        box_side_points.push_back(point_ct);
        point_ct++;

        // Write Airfoil Points - Clockwise
        for(auto xy_pt:airfoil_pts){
            out_file << "Point(" << point_ct << ") = {" 
                     << xy_pt.front() << ", " << xy_pt.back() << ", 0.0, " 
                     << grid_size << "};" << std::endl; 

            af_side_points.push_back(point_ct);    
            point_ct++;
        }

        // DEFINE LINES
        int line_ct {1}; // GMSH Formatting Starts at 1
        
        std::vector<int> af_side_lines{};
        std::vector<int> box_side_line{};


        // Write Side Box Lines - Counterclockwise
        for(auto pt_num:box_side_points){
            if (pt_num == box_side_points.back()){
                out_file << "Line(" << (line_ct) << ") = {" << pt_num << ", " << box_side_points.front() << "};" << std::endl;
            }else {
                out_file << "Line(" << (line_ct) << ") = {" << pt_num << ", " << pt_num+1 << "};" << std::endl;
            }
            box_side_line.push_back(line_ct);
            line_ct++; 
        }

        // Write Side Airfoil Lines - Clockwise
        for(auto pt_num:af_side_points){
            if (pt_num == af_side_points.back()){
                out_file << "Line(" << (line_ct) << ") = {" << pt_num << ", " << af_side_points.front() << "};" << std::endl;
            }else {
                out_file << "Line(" << (line_ct) << ") = {" << pt_num << ", " << pt_num+1 << "};" << std::endl;
            }
            af_side_lines.push_back(line_ct);
            line_ct++; 
        }

        // DEFINE CURVE LOOPS
        int curve_ct {line_ct}; // GMSH numbers curves as a new line

        int af_side_curve{};
        int box_side_curve{};

        // Box Side Curve Loop - Face
        out_file << "Curve Loop(" << (curve_ct) << ") = {";
        for (auto ln_num:box_side_line){
            if(ln_num == box_side_line.back()){
                out_file << ln_num << "};" << std::endl;
            }else{
                out_file << ln_num << ", ";
            }
        }
        box_side_curve = curve_ct;
        curve_ct++;

        // Airfoil Side Curve Loop - Hole
        out_file << "Curve Loop(" << (curve_ct) << ") = {";
        for (auto ln_num:af_side_lines){
            if(ln_num == af_side_lines.back()){
                out_file << ln_num << "};" << std::endl;
            }else{
                out_file << ln_num << ", ";
            }
        }
        af_side_curve = curve_ct;
        curve_ct++;

        // DEFINE PLANE SURFACE
        int surface_ct {curve_ct};

        int side_psurface{};

        out_file << "Plane Surface(" << (surface_ct) << ") = {" 
                 << (box_side_curve) << ", " << (af_side_curve) << "};" << std::endl; 
        side_psurface = surface_ct;
        surface_ct++;

        // DEFINE BOUNDARY LAYER
        // Field[1] = BoundaryLayer; // hwall * ratio^(dist/hwall)
        // Field[1].AnisoMax = 1.0; // Threshold angle for creating a mesh fan in the boundary layer
        // Field[1].EdgesList = {};  // EdgesList for 2D: Tags of curves in the geometric model for which a boundary layer is needed  
        // Field[1].FacesList = {};  // FacesList for 3D: Tags of faces in the geometric model for which a boundary layer is needed  
        // Field[1].ExcludedFaceList = {}; // Tags of surfaces in the geometric model where the boundary layer should not be applied
        // Field[1].FanNodesList = {}; // Tags of points in the geometric model for which a fan is created
        // Field[1].IntersectMetrics = 0; // Intersect metrics of all faces
        // Field[1].NodesList = {}; //Tags of points in the geometric model for which a boundary layer ends
        // Field[1].Quads = 1; // Generate recombined elements in the boundary layer
        // Field[1].hfar = 0.03; // Element size far from the wall
        // Field[1].hwall_n = 0.001; // Mesh Size Normal to the The Wall
        // Field[1].hwall_n_nodes = {}; // Mesh Size Normal to the The Wall at nodes (overwrite hwall_n when defined)
        // Field[1].hwall_t = 0.001; // ???
        // Field[1].ratio = 1.1; // Size Ratio Between Two Successive Layers
        // Field[1].thickness = 0.05;  // Maximal thickness of the boundary layer, must be bigger than hwall_n
        // BoundaryLayer Field = 1; // Maximal thickness of the boundary layer
        
        out_file << "Field[1] = BoundaryLayer;" << std::endl;
        out_file << "Field[1].AnisoMax = " << anisomax << ";" << std::endl;
        out_file << "Field[1].EdgesList = {";
        for(auto af_line:af_side_lines){
            if(af_line == af_side_lines.back()){
                out_file << af_line << "};" << std::endl;
            }else{
                out_file << af_line << ", ";
            }
        }
        if(double_fan){
            out_file << "Field[1].FanNodesList = {" 
                 << af_side_points.front() << ", " << af_side_points.back() << "};" << std::endl;
        }
        else{
            out_file << "Field[1].FanNodesList = {" 
                 << af_side_points.front() << "};" << std::endl;
        }
        out_file << "Field[1].Quads = " << quads << ";" << std::endl;
        out_file << "Field[1].hfar = " << hfar << ";" << std::endl;
        out_file << "Field[1].hwall_n = " << hwall_n << ";" << std::endl; 
        out_file << "Field[1].ratio = " << ratio << ";" << std::endl;
        out_file << "Field[1].thickness = " << thickness << ";" << std::endl;
        out_file << "BoundaryLayer Field = 1;" << std::endl;

        // DEFINE EXTRUSION
        out_file << "surfaceVector[] = Extrude {0, 0, " << mesh_thickness << "} {" << std::endl;
        out_file << "   Surface{" << side_psurface << "};" << std::endl;
        out_file << "   Layers{1};" << std::endl;
        out_file << "   Recombine;" << std::endl;
        out_file << "};" << std::endl;

        // DEFINE PHYSICAL SURFACES AND ASSIGN FROM SURFACE VECTOR
        
        // surfaceVector contains in the following order:
        // [0] - front surface (opposed to source surface)
        // [1] - extruded volume
        // [2] - inlet surface (belonging to 1st arg in curve loop defining the box)
        // [3] - top surface (belonging to 2nd arg in curve loop defining the box)
        // [4] - outlet surface (belonging to 3rd arg in curve loop defining the box)
        // [5] - bottom surface (belonging to 4th arg in curve loop defining the box)
        // [6] - first airfoil transverse surface (belonging to 1st arg in curve loop defining airfoil)
        // [...]
        // [6+n] - where n is the number of airfoil lines, last airfoil transverse surface (belonging to n+6th arg in curve loop defining airfoil)

        out_file << "Physical Surface(\"Inlet\") = {surfaceVector[2]};" << std::endl;
        out_file << "Physical Surface(\"Top\") = {surfaceVector[3]};" << std::endl;
        out_file << "Physical Surface(\"Outlet\") = {surfaceVector[4]};" << std::endl;
        out_file << "Physical Surface(\"Bottom\") = {surfaceVector[5]};" << std:: endl;
        out_file << "Physical Surface(\"Sides\") = {surfaceVector[0], " << side_psurface << "};" << std::endl;
        out_file << "Physical Surface(\"Foil\") = {";
        for (long long unsigned int i = 6; i < (af_side_lines.size()+6); i++){
            if (i == ((af_side_lines.size()+6)-1)){
                out_file << "surfaceVector[" << i << "]};" << std::endl;
            }else{
                out_file << "surfaceVector[" << i << "], ";
            }
        }

        // DEFINE PHYSICAL VOLUMES
        out_file << "Physical Volume(\"Vol\") = surfaceVector[1];" << std::endl;  

        //! RECOMBINE EXTERNAL SURFACES - BOUNDARY LAYER DOES NOT LIKE THIS
        // out_file << "Recombine Surface {" << side_psurface << "};" << std::endl;
        // out_file << "Recombine Surface {surfaceVector[0]};" << std::endl;

    }

    // CLOSE FILES
    in_file.close();
    out_file.close();

    return 0;
}