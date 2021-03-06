// TODO: Make more input based
// TODO: Make more function based
// TODO: Experiment with splines instead of lines: BSpline ( expression ) = { expression-list }; and Spline ( expression ) = { expression-list };

// ! GMSH Needs to export the 3D mesh in Version2 ASCII Format

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>

int main(){
    // AIRFOILE FILE NAMES
    std::string in_airfoil_file_name{"./.AIRFOILS_IN/NACA2412.txt"};
    std::string out_airfoil_file_name {"./.AIRFOILS_OUT/NACA2412.geo"};

    // MESH PARAMETERS
    double square_side {30.0};
    double grid_size {square_side/50};
    double mesh_thickness {square_side/10};

    // Y+ ESTIMATION
    double freestream_vel {30};                                                         //[m/s]
    double char_length {1};                                                             //[m]
    double density {1.225};                                                             //[kg/m^3]
    double dyn_visc {0.00001789};                                                       //kg/m.s]
    double kin_visc {dyn_visc/density};                                                 //[m^2/s]
    double Re_num {(density*freestream_vel*char_length)/dyn_visc};                      //[]
    double skin_fric_coeff_over_two {0.0359*std::pow(Re_num,-0.2)};                     //[]
    double fric_vel {std::abs(freestream_vel)*std::sqrt(skin_fric_coeff_over_two)};     //[]
    double y_plus {1};                                                                //[]
    double cell_ctr_height {y_plus*kin_visc/fric_vel};                                  //[m]    

    std::cout << "Kinematic Viscosity: " << kin_visc << std::endl;
    std::cout << "Skin Friction Coeff Cf/2: " << skin_fric_coeff_over_two << std::endl;
    std::cout << "Friction Velocity " << fric_vel << std::endl;
    std::cout << "Reynolds Number: " << Re_num << std::endl;   
    std::cout << "Near wall cell centre height: " << cell_ctr_height << std::endl;   

    // BOUNDARY LAYER PARAMETERS
    double anisomax {90}; // Threshold angle for creating a mesh fan in the boundary layer
    int quads {1}; // Generate recombined elements in the boundary layer
    double hfar {0.1}; // Element size far from the wall
    double hwall_n {cell_ctr_height*2}; // Mesh Size Normal to the The Wall
    double ratio {1.05}; // Size Ratio Between Two Successive Layers
    double thickness {0.95*(square_side/2)}; // must be bigger than hwall_n

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
        // ! NEEDS SOME SERIOUS CLEANUP
        bool double_fan {false};
        if(xy_pt_first.front() == xy_pt_last.front()){
            double_fan =true;

            // Calculate Top Edge Vector
            std::vector<double> temp_A_pt {airfoil_pts.at(0)};
            std::vector<double> temp_B_pt {airfoil_pts.at(1)};
            std::vector<double> top_edge_vec {};

            top_edge_vec.push_back(temp_B_pt.front() - temp_A_pt.front());
            top_edge_vec.push_back(temp_B_pt.back() - temp_A_pt.back());

            // Calculate Mid Point Length           
            double mid_pt_length {(xy_pt_first.back() - xy_pt_last.back())/2};

            // Calculate Perpendicular Vector and Radius from Mid Point
            std::vector<double> temp_C_vec {};
            double temp_y {mid_pt_length+xy_pt_last.back()};
            double temp_x {((top_edge_vec.back()*(xy_pt_first.back()-temp_y))/top_edge_vec.front())+xy_pt_first.front()};

            temp_C_vec.push_back(xy_pt_first.front() - temp_x);
            temp_C_vec.push_back(xy_pt_first.back() - temp_y);

            double radius {std::sqrt((temp_C_vec.front()*temp_C_vec.front()) + (temp_C_vec.back()*temp_C_vec.back()))};

            // Calculate Arc Center
            double arc_center_x {xy_pt_first.front() - radius};
            double arc_center_y {mid_pt_length+xy_pt_last.back()};

            // Make Points
            std::vector <double> temp_new_pt {};

            airfoil_pts.erase(airfoil_pts.begin());
            airfoil_pts.pop_back();

            temp_new_pt.push_back(arc_center_x + temp_C_vec.front());
            temp_new_pt.push_back(arc_center_y + temp_C_vec.back());
            airfoil_pts.insert(airfoil_pts.begin(),temp_new_pt);
            temp_new_pt.clear();

            // ! THE INCREMENT AND FOR LOOP CONDITIONS SHOULD BE LINKED 
            double y_increment {mid_pt_length/40};
            double temp_new_y_pt {-1*temp_C_vec.back()};
            for (size_t i = 0; i<80; i++){
                double temp_new_x_pt {std::sqrt((radius*radius)-(temp_new_y_pt*temp_new_y_pt))};

                temp_new_y_pt = temp_new_y_pt + arc_center_y;
                temp_new_x_pt = temp_new_x_pt + arc_center_x;

                temp_new_pt.push_back(temp_new_x_pt);
                temp_new_pt.push_back(temp_new_y_pt);
                airfoil_pts.push_back(temp_new_pt);
                temp_new_pt.clear();

                temp_new_y_pt = temp_new_y_pt + y_increment;
            }
        }     
        
        // ! UNDER CONSTRUCTION
        // REFINE ASPECT RATIO
        // Variables
        std::vector<std::vector<double>> airfoil_pts_dummy {airfoil_pts};
        airfoil_pts.clear();
        double x1 {};
        double y1 {};
        double x2 {};
        double y2 {};
        double line_x {};
        double line_y {};
        double line_length{};
        double line_AR{};

        // Refinement Loop
        // ! Need to add stuff for last point
        for(auto current_pt:airfoil_pts_dummy){
            if(current_pt == airfoil_pts_dummy.front()){
                x1 = current_pt.front();
                y1 = current_pt.back();
                airfoil_pts.push_back(current_pt);
            }else{
                std::vector<double> temp_pt {};
                
                x2 = current_pt.front();
                y2 = current_pt.back();
                line_x = x2-x1;
                line_y = y2-y1;
                line_length = std::sqrt((line_x*line_x)+(line_y*line_y));
                line_AR = line_length/hwall_n;

                if(line_AR>5){
                    double new_line_length {line_length};
                    double new_line_AR {line_AR};
                    size_t divisor {1};
                    
                    while(new_line_AR>100){
                        divisor = divisor*2;
                        new_line_length = line_length/divisor;
                        new_line_AR = new_line_length/hwall_n;
                    }

                    for(size_t i {1}; i<divisor; i++){
                        temp_pt.push_back(x1 + (i*(line_x/divisor)));
                        temp_pt.push_back(y1 + (i*(line_y/divisor)));
                        airfoil_pts.push_back(temp_pt);
                        temp_pt.clear();
                    }
                    airfoil_pts.push_back(current_pt);
                }else{
                    airfoil_pts.push_back(current_pt);
                }

                // Update Points
                x1 = x2;
                y1 = y2;
            } 
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
            // out_file << "Field[1].FanNodesList = {" 
            //      << af_side_points.front() << ", " << af_side_points.back() << "};" << std::endl;
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