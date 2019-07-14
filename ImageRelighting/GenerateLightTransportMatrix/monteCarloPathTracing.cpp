#include "monteCarloPathTracing.h"


double random0to1()
{
    static double invRAND_MAX = 1.0/RAND_MAX;
    return rand()*invRAND_MAX;
}


//double cube_bbmin[3] = { -285, -285, -285 }, cube_bbmax[3] = { 797, 797, 797 };
double cube_bbmin[3] = { -500,-500,-500 }, cube_bbmax[3] = { 1012,1012,1012 };
//double cube_bbmin[3] = { -717, -717, -717 }, cube_bbmax[3] = { 1229, 1229, 1229 };



double cylinder_center[11][3] = { {170, 0, 120}, {340, 0, 120}, {85, 0, 200}, {255, 0, 200}, 
                                 {425, 0, 200}, {140, 0, 350}, {370, 0, 350}, {180, 0, 440}, {330, 0, 440}, {80, 0, 450}, 
                                 {430, 0, 450} };
double cylinder_radius[11] = { 40, 40, 20, 20, 20, 15, 15, 20, 20, 30, 30,  };
double cylinder_height[11] = { 350, 350, 300, 230, 300, 230, 230, 100, 100, 150, 150 };


double sphere_center[3][3] = { {256, 20, 350}, };
double sphere_radius[3]    = {50};


int FindNearestIntersection(double *org, double *dir, int current_object_index, double *hit_point, int depth)
{

    int    nearest_object_index = -1;
    double nearest_t = 1.0e6+1;
   
    double inv_dir[3] = { 1.0/dir[0], 1.0/dir[1], 1.0/dir[2], };

   
    // intersection test between ray and 6 walls constituting cube maps
    for(int i = 0; i < 6; i++){

        if(i == current_object_index) continue;
        if(depth == 0 && i == 4)      continue; // ray from camera can go through back wall. (otherwise, nothing is visible)

        double t;

        switch(i){
        case 0: t = (cube_bbmax[1]-org[1])*inv_dir[1];  break;
        case 1: t = (cube_bbmin[1]-org[1])*inv_dir[1];  break;
        case 2: t = (cube_bbmin[2]-org[2])*inv_dir[2];  break;
        case 3: t = (cube_bbmax[0]-org[0])*inv_dir[0];  break;
        case 4: t = (cube_bbmax[2]-org[2])*inv_dir[2];  break;
        case 5: t = (cube_bbmin[0]-org[0])*inv_dir[0];  break;
        }


        if(t > 0.0 && t < nearest_t){
            double intersect_1, intersect_2;

            switch(i){
            case 0:
            case 1:
                intersect_1 = org[0]+dir[0]*t;
                intersect_2 = org[2]+dir[2]*t;
                if(cube_bbmin[0] <= intersect_1 && intersect_1 <= cube_bbmax[0] && 
                   cube_bbmin[2] <= intersect_2 && intersect_2 <= cube_bbmax[2]){
                    nearest_t = t;
                    nearest_object_index = i;
                }
                break;
            case 2:
            case 4:
                intersect_1 = org[0]+dir[0]*t;
                intersect_2 = org[1]+dir[1]*t;
                if(cube_bbmin[0] <= intersect_1 && intersect_1 <= cube_bbmax[0] && 
                   cube_bbmin[1] <= intersect_2 && intersect_2 <= cube_bbmax[1]){
                    nearest_t = t;
                    nearest_object_index = i;
                }
                break;
            case 3:
            case 5:
                intersect_1 = org[1]+dir[1]*t;
                intersect_2 = org[2]+dir[2]*t;
                if(cube_bbmin[1] <= intersect_1 && intersect_1 <= cube_bbmax[1] && 
                   cube_bbmin[2] <= intersect_2 && intersect_2 <= cube_bbmax[2]){
                    nearest_t = t;
                    nearest_object_index = i;
                }
                break;
            }

        } // if(t > nearest_t){

    } // for(int i = 1; i <= 6; i++){


    // intersection test between ray and table
    if(current_object_index != 6){
        double t = -org[1]*inv_dir[1];

        if(t > 0.0 && t < nearest_t){
            double intersect_1 = org[0]+dir[0]*t;
            double intersect_2 = org[2]+dir[2]*t;

            if(0.0 <= intersect_1 && intersect_1 <= 512.0 && 0.0 <= intersect_2 && intersect_2 <= 512.0){
                nearest_t = t;
                nearest_object_index = 6;
            }
        }
    } // if(current_object_index != 6){


    // intersection test between ray and cylinders standing along y-axis

    double a = dir[0]*dir[0] + dir[2]*dir[2];
    double inv_a = 1.0/a;

    for(int i = 0; i < 11; i++){

        if(i+7 == current_object_index) continue;
        
        // perform intersecton between ray and circle that is the slice of the cylinder
        double b = 2.0 * (dir[0]*(org[0]-cylinder_center[i][0]) + dir[2]*(org[2]-cylinder_center[i][2]));
        double c = (org[0]-cylinder_center[i][0])*(org[0]-cylinder_center[i][0]) + (org[2]-cylinder_center[i][2])*(org[2]-cylinder_center[i][2]) - cylinder_radius[i]*cylinder_radius[i];

        double discriminant = b*b-4*a*c;

        if(discriminant >= 0.0){
            double t1 = (-b-sqrt(discriminant))*(0.5*inv_a);
            double t2 = (-b+sqrt(discriminant))*(0.5*inv_a);

            if(t1 > 0.0 && t1 < nearest_t){
                double hit_point_y = org[1] + dir[1]*t1;

                if(hit_point_y > cylinder_center[i][1] && hit_point_y < cylinder_center[i][1]+cylinder_height[i]){
                    nearest_t = t1;
                    nearest_object_index = i+7;
                }
            }
            if(t2 > 0.0 && t2 < nearest_t){
                double hit_point_y = org[1] + dir[1]*t2;

                if(hit_point_y > cylinder_center[i][1] && hit_point_y < cylinder_center[i][1]+cylinder_height[i]){
                    nearest_t = t2;
                    nearest_object_index = i+7;
                }
            }

        } // if(discriminant >= 0.0){

        // intersection between ray and cap of cylinder
        double t = (cylinder_height[i] - org[1]) * inv_dir[1]; 

        if(t > 0.0 && t < nearest_t){
            double hit_point_x = org[0] + dir[0]*t;
            double hit_point_z = org[2] + dir[2]*t;

            double diff_to_center_x = hit_point_x-cylinder_center[i][0];
            double diff_to_center_z = hit_point_z-cylinder_center[i][2];

            if(diff_to_center_x*diff_to_center_x + diff_to_center_z*diff_to_center_z < cylinder_radius[i]*cylinder_radius[i]){
                nearest_t = t;
                nearest_object_index = i+7;
            }
        }

    }


    // intersection test between ray and 2 spheres

    a     = DotProduct(dir, dir);
    inv_a = 1.0/a;
    double dot_org_dir = DotProduct(org, dir);
    double dot_org_org = DotProduct(org, org);

    for(int i = 0; i < 1; i++){
        if(i+18 == current_object_index) continue;

        double b = 2.0*(dot_org_dir-DotProduct(dir, sphere_center[i]));
        double c = dot_org_org - 2.0*DotProduct(org, sphere_center[i]) + DotProduct(sphere_center[i], sphere_center[i]) - sphere_radius[i]*sphere_radius[i];

        double discriminant = b*b-4*a*c;

        if(discriminant >= 0.0){
            double t = (-b-sqrt(discriminant))*(0.5*inv_a);

            if(i+18 == current_object_index && fabs(t) < 1.0e-6){
                // ray passes through inside of the sphere
                return -1;
            }else if(i+18 != current_object_index){

                if(t > 0.0 && t < nearest_t){
                    nearest_t = t;
                    nearest_object_index = i+18;
                }
            }
        }

    }


    for(int i = 0; i < 3; i++) hit_point[i] = org[i] + dir[i]*nearest_t;

   
    return nearest_object_index;
}



void GetMaterialProperty(double *hit_point, int object_index, bool isLight, 
                         double *object_normal, double *object_diffuse, double *object_specular, double &object_shininess, 
                         double &probability_for_diffuse_out, double &inv_probability_for_diffuse_out, double &inv_probability_for_specular_out)
{

    if(object_index == 0){
        object_diffuse[0] = 0.0;
        object_diffuse[1] = 0.0; 
        object_diffuse[2] = 0.0;
        object_specular[0] = object_specular[1] = object_specular[2] = 0.0;
        object_shininess = 1.0;
    }else if(object_index == 1){
        object_diffuse[0] = object_diffuse[1] = object_diffuse[2] = 1.0;
        object_specular[0] = object_specular[1] = object_specular[2] = 0.0;
        object_shininess = 1.0;
    }else if(object_index >= 2 && object_index <= 5){

        static double inv_cube_size_x = 1.0 / (cube_bbmax[0]-cube_bbmin[0]);
        static double inv_cube_size_y = 1.0 / (cube_bbmax[1]-cube_bbmin[1]);
        static double inv_cube_size_z = 1.0 / (cube_bbmax[2]-cube_bbmin[2]);

        int pixel_coord_x = 64*(object_index-2);

        switch(object_index){
        case 2: pixel_coord_x += (int)((hit_point[0]-cube_bbmin[0]) * inv_cube_size_x * 64);  break;
        case 3: pixel_coord_x += (int)((hit_point[2]-cube_bbmin[2]) * inv_cube_size_z * 64);  break;
        case 4: pixel_coord_x += (int)((cube_bbmax[0]-hit_point[0]) * inv_cube_size_x * 64);  break;
        case 5: pixel_coord_x += (int)((cube_bbmax[2]-hit_point[2]) * inv_cube_size_z * 64);  break;
        }

        int pixel_coord_y = (int)((hit_point[1]-cube_bbmin[1]) * inv_cube_size_y * 64);

        int pixel_index = (pixel_coord_y*png_width + pixel_coord_x) * 3;

        object_diffuse[0] = png_data[pixel_index+0];
        object_diffuse[1] = png_data[pixel_index+1];
        object_diffuse[2] = png_data[pixel_index+2];

        object_specular[0] = object_specular[1] = object_specular[2] = 0.0;
        object_shininess = 1.0;
    }


    if(object_index == 6){

        object_normal[0]  = 0.0;  object_normal[1] = 1.0;  object_normal[2]  = 0.0; 

        object_diffuse[0]  = object_diffuse[1]  = object_diffuse[2] = 0.9;
        object_specular[0] = object_specular[1] = object_specular[2] = 0.1;
        object_shininess = 10.0;

    }else if(object_index >= 7 && object_index <= 17){

        int id = object_index-7;

        double diff_to_center_x = hit_point[0] - cylinder_center[id][0];
        double diff_to_center_z = hit_point[2] - cylinder_center[id][2];

        if(diff_to_center_x*diff_to_center_x + diff_to_center_z*diff_to_center_z < cylinder_radius[id]*cylinder_radius[id]){
            // ray hits with the cap of cylinder
            object_normal[0]  = 0.0;  object_normal[1]  =  1.0;  object_normal[2]  = 0.0; 
        }else{
            object_normal[0] = hit_point[0]-cylinder_center[id][0];
            object_normal[1] = 0.0;
            object_normal[2] = hit_point[2]-cylinder_center[id][2];

            Normalize(object_normal);
        }

        if(object_index == 12 || object_index == 13){
            object_diffuse[0]  = 0.5;  object_diffuse[1] = 0.3;  object_diffuse[2] = 0.1;
            object_specular[0] = 0.5; object_specular[1] = 0.1; object_specular[2] = 0.1;
            object_shininess = 2.0;
        }else if(object_index == 10){
            object_diffuse[0]  = 0.0;  object_diffuse[1] = 0.3;  object_diffuse[2] = 0.5;
            object_specular[0] = 0.1; object_specular[1] = 0.1; object_specular[2] = 0.5;
            object_shininess = 2.0;
        }else if(object_index == 14 || object_index == 15){
            object_diffuse[0]  = 0.8;  object_diffuse[1] = 0.5;  object_diffuse[2] = 0.0;
            object_specular[0] = 0.2; object_specular[1] = 0.3; object_specular[2] = 0.0;
            object_shininess = 2.0;
        }else{
            object_diffuse[0]  = 0.2;  object_diffuse[1] = 0.2;  object_diffuse[2] = 0.2;
            object_specular[0] = 0.8; object_specular[1] = 0.8; object_specular[2] = 0.8;
            object_shininess = 5000.0;
        }
    }else{
        for(int i = 0; i < 3; i++) object_normal[i] = hit_point[i]-sphere_center[0][i];
        Normalize(object_normal); 

        object_diffuse[0]  = 0.2;  object_diffuse[1] = 0.2;  object_diffuse[2] = 0.2;
        object_specular[0] = 0.8; object_specular[1] = 0.8; object_specular[2] = 0.8;
        object_shininess = 5000.0;
    }

    static bool is_probability_diffuse_computed[19];
    static bool init = true;

    if(init){
        for(int i = 0; i < 19; i++) is_probability_diffuse_computed[i] = false;
        init = false;
    }

    // precompute probability to shoot diffuse ray and specular ray for each material 
    static double probability_for_diffuse[19];
    static double inv_probability_for_diffuse[19], inv_probability_for_specular[19]; 

    if(is_probability_diffuse_computed[object_index] == false){
        double sum_diffuse  =  object_diffuse[0] +  object_diffuse[1] +  object_diffuse[2];
        double sum_specular = object_specular[0] + object_specular[1] + object_specular[2];

        probability_for_diffuse[object_index] = sum_diffuse / (sum_diffuse+sum_specular);

        inv_probability_for_diffuse[object_index]  = (sum_diffuse+sum_specular) / sum_diffuse;
        if(sum_specular > 1.0e-6) inv_probability_for_specular[object_index] = (sum_diffuse+sum_specular) / sum_specular;
        else                      inv_probability_for_specular[object_index] = 0.0;

        is_probability_diffuse_computed[object_index] = true;
    }

    probability_for_diffuse_out      = probability_for_diffuse[object_index];
    inv_probability_for_diffuse_out  = inv_probability_for_diffuse[object_index];
    inv_probability_for_specular_out = inv_probability_for_specular[object_index];
}


void GetPixelColor(int cubemap_id, int *cubemap_pixel_coord, unsigned int *pixel_color, double view_theta, double view_phi)
{
    static double cubemap_theta[6][CUBEMAP_RESOLUTION][CUBEMAP_RESOLUTION];
    static double   cubemap_phi[6][CUBEMAP_RESOLUTION][CUBEMAP_RESOLUTION];

    static bool init = true;

    if(init){
        double inv_cube_size_x = 1.0 / (cube_bbmax[0]-cube_bbmin[0]);
        double inv_cube_size_y = 1.0 / (cube_bbmax[1]-cube_bbmin[1]);
        double inv_cube_size_z = 1.0 / (cube_bbmax[2]-cube_bbmin[2]);

        double inv_cubemap_resolution = 1.0 / (double)CUBEMAP_RESOLUTION;

        for(int i = 0; i < 6; i++){

            for(int j = 0; j < CUBEMAP_RESOLUTION; j++){
                for(int k = 0; k < CUBEMAP_RESOLUTION; k++){

                    double retrieve_vector[3];

                    if(i == 0 || i == 1){
                        retrieve_vector[0] = (j+0.5)*inv_cubemap_resolution*2.0 - 1.0;
                        retrieve_vector[2] = (k+0.5)*inv_cubemap_resolution*2.0 - 1.0;
                        if(i == 0) retrieve_vector[1] =  1.0;
                        else       retrieve_vector[1] = -1.0; 
                    }else if(i == 2 || i == 4){
                        retrieve_vector[0] = (j+0.5)*inv_cubemap_resolution*2.0 - 1.0;
                        retrieve_vector[1] = (k+0.5)*inv_cubemap_resolution*2.0 - 1.0;
                        if(i == 2) retrieve_vector[2] = -1.0;
                        else       retrieve_vector[2] =  1.0;
                    }else{ // i == 3 || i == 5
                        retrieve_vector[2] = (j+0.5)*inv_cubemap_resolution*2.0 - 1.0;
                        retrieve_vector[1] = (k+0.5)*inv_cubemap_resolution*2.0 - 1.0;
                        if(i == 3) retrieve_vector[0] =  1.0;
                        else       retrieve_vector[0] = -1.0;
                    }



                    Normalize(retrieve_vector);

                    cubemap_theta[i][j][k] = acos(-retrieve_vector[1]); // [0, PI]
                    cubemap_phi[i][j][k]   = atan2(retrieve_vector[2], retrieve_vector[0]) + PI; // [0, 2PI]
                }
            }

        } // for(int i = 0; i < 6; i++){

        init = false;
    }


    static unsigned int filtered_pixel_data[6][CUBEMAP_RESOLUTION][CUBEMAP_RESOLUTION][3];


    if(isNewLighting){    

        for(int i = 0; i < 6; i++){

            for(int j = 0; j < CUBEMAP_RESOLUTION; j++){
                for(int k = 0; k < CUBEMAP_RESOLUTION; k++){

                    double theta = cubemap_theta[i][j][k] + view_theta;
                    double phi   = cubemap_phi[i][j][k]   + view_phi;

                    if(theta > PI){
                        theta = 2.0*PI-theta;
                        phi   += PI;
                    }else if(theta < 0.0){
                        theta = -theta;
                        phi   += PI;
                    }

                    if(phi < 0.0)         phi += 2.0*PI;
                    else if(phi > 2.0*PI) phi -= 2.0*PI;

                    int pixel_coord_x = (phi)*(InvPI*0.5) * png_width;
                    int pixel_coord_y = (theta*InvPI) * png_height;

                    if(pixel_coord_x < 0)                     pixel_coord_x = 0;
                    else if(pixel_coord_x >= (int)png_width)  pixel_coord_x = png_width-1;
                    if(pixel_coord_y < 0)                     pixel_coord_y = 0;
                    else if(pixel_coord_y >= (int)png_height) pixel_coord_y = png_height-1;

                    int pixel_index = (pixel_coord_y*png_width + pixel_coord_x) * 3;

                    for(int m = 0; m < 3; m++) filtered_pixel_data[i][j][k][m] = png_data[pixel_index+m];
                 }
            }

        } // for(int i = 0; i < 6; i++){


        isNewLighting = false;

    } // if(isNewLighting){ 


    // change this from double to char to save memory

    pixel_color[0] = filtered_pixel_data[cubemap_id][ cubemap_pixel_coord[0] ][ cubemap_pixel_coord[1] ][0];
    pixel_color[1] = filtered_pixel_data[cubemap_id][ cubemap_pixel_coord[0] ][ cubemap_pixel_coord[1] ][1];
    pixel_color[2] = filtered_pixel_data[cubemap_id][ cubemap_pixel_coord[0] ][ cubemap_pixel_coord[1] ][2];
}




void GetRandomDirection(double *normal, double *random_direction)
{
    do{
        for(int i = 0; i < 3; i++) random_direction[i] = (random0to1()-0.5)*2.0;
    }while( DotProduct(random_direction, random_direction) > 1.0);

    Normalize(random_direction);

    if( DotProduct(random_direction, normal) < 0.0){ 
        // flip direction
        for(int i = 0; i < 3; i++) random_direction[i] = -random_direction[i];
    }

}

void GetDirectionUsingCosWeightedLambertian(double *point, double *normal, double *random_direction)
{
    // "outgoing" in a physical sense
    // "outgoing_direction" is actually an incoming direction in path tracing

    double theta = acos( sqrt( random0to1() ) );
    double phi   = random0to1()*TwoPI;

    // temporary created rondom direction
    double temp[3] = { sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) };

    // rotate "temp" such that z-axis used for sampling is aligned with "normal"
    double angle_between = acos(normal[2]);
    double rotation_axis[3] = { -normal[1], normal[0], 0.0 };

    Normalize(rotation_axis);
    
    RotateAroundAxis(rotation_axis, angle_between, temp, random_direction);
}


void GetDirectionUsingPhongBRDF(double *point, double *normal, double *outgoing_direction, double shininess, double *random_direction)
{
    // "outgoing" in a physical sense
    // "outgoing_direction" is actually an incoming direction in path tracing

    double reflection_vector[3];

    ComputeReflectionVector(normal, outgoing_direction, reflection_vector);

    do{
        double theta = acos( pow( random0to1(), 1.0/(shininess+1.0) ) );
        double phi   = random0to1()*TwoPI;

        // temporary created rondom direction
        double temp[3] = { sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) };

        // rotate "temp" such that z-axis used for sampling is aligned with "reflection_vector"
        double angle_between = acos(reflection_vector[2]);
        double rotation_axis[3] = { -reflection_vector[1], reflection_vector[0], 0.0 };

        Normalize(rotation_axis);

        RotateAroundAxis(rotation_axis, angle_between, temp, random_direction);

    }while( DotProduct(random_direction, normal) < 0.0); // random_direction must be directed to the direction of normal

}



void TracePathWithImportanceSampling(double *eye, double *ray, double *color, int &cubemap_id, int *cubemap_pixel_coord, double weight_of_path, int depth)
{
    double probability_to_kill     = 0.0;
    double inv_probability_to_kill = 1.0;

    if(weight_of_path < 0.001 || depth > MAX_DEPTH){ // if weight_of_path < 0.001, we stop recursion with probability of 50%
        probability_to_kill     = 0.5;
        inv_probability_to_kill = 1.0 / (1.0-probability_to_kill);

        if(random0to1() < probability_to_kill){
            color[0] = color[1] = color[2] = 0.0;
            return;
        }
    }


    int    nearest_object_index;
    double hit_point[3], object_normal[3], object_diffuse[3], object_specular[3], object_shininess;
    bool   isLight = false;

    double probability_for_diffuse, inv_probability_for_diffuse, inv_probability_for_specular; 

     if( (nearest_object_index = FindNearestIntersection(eye, ray, -1, hit_point, depth)) != -1 ){

        GetMaterialProperty(hit_point, nearest_object_index, isLight, object_normal, object_diffuse, object_specular, object_shininess, 
                            probability_for_diffuse, inv_probability_for_diffuse, inv_probability_for_specular);

 
        if(nearest_object_index < 6){

            cubemap_id = nearest_object_index;

            static double inv_cube_size_x = 1.0 / (cube_bbmax[0]-cube_bbmin[0]);
            static double inv_cube_size_y = 1.0 / (cube_bbmax[1]-cube_bbmin[1]);
            static double inv_cube_size_z = 1.0 / (cube_bbmax[2]-cube_bbmin[2]);

            if(cubemap_id == 0 || cubemap_id == 1){
                cubemap_pixel_coord[0] = ((hit_point[0]-cube_bbmin[0]) * inv_cube_size_x * CUBEMAP_RESOLUTION);
                cubemap_pixel_coord[1] = ((hit_point[2]-cube_bbmin[2]) * inv_cube_size_z * CUBEMAP_RESOLUTION);
            }else if(cubemap_id == 2 || cubemap_id == 4){
                cubemap_pixel_coord[0] = ((hit_point[0]-cube_bbmin[0]) * inv_cube_size_x * CUBEMAP_RESOLUTION);
                cubemap_pixel_coord[1] = ((hit_point[1]-cube_bbmin[1]) * inv_cube_size_y * CUBEMAP_RESOLUTION);
            }else{ // cubemap_id == 3 || cubemap_id == 5
                cubemap_pixel_coord[0] = ((hit_point[2]-cube_bbmin[2]) * inv_cube_size_z * CUBEMAP_RESOLUTION);
                cubemap_pixel_coord[1] = ((hit_point[1]-cube_bbmin[1]) * inv_cube_size_y * CUBEMAP_RESOLUTION);
            }

            if(cubemap_pixel_coord[0] < 0)                        cubemap_pixel_coord[0] = 0;
            else if(cubemap_pixel_coord[0] >= CUBEMAP_RESOLUTION) cubemap_pixel_coord[0] = CUBEMAP_RESOLUTION-1;
            if(cubemap_pixel_coord[1] < 0)                        cubemap_pixel_coord[1] = 0;
            else if(cubemap_pixel_coord[1] >= CUBEMAP_RESOLUTION) cubemap_pixel_coord[1] = CUBEMAP_RESOLUTION-1;


            // test this function
            // GetPixelColor(cubemap_id, cubemap_pixel_coord, object_diffuse, 0, 0);

            // light
            double light_intensity_scale;

            if(depth == 0) light_intensity_scale = 1.0;
            else           light_intensity_scale = 4.0;

            for(int i = 0; i < 3; i++) color[i] = inv_probability_to_kill * light_intensity_scale;// * object_diffuse[i];


            return;

        }else{
            // not light

            double incoming_direction[3], incoming_color[3];

            if( random0to1() <= probability_for_diffuse ){
                // evaluate diffuse

                GetDirectionUsingCosWeightedLambertian(hit_point, object_normal, incoming_direction);

                double weight = inv_probability_to_kill * inv_probability_for_diffuse;

                TracePathWithImportanceSampling(hit_point, incoming_direction, incoming_color, cubemap_id, cubemap_pixel_coord, weight_of_path*weight, depth+1);

                for(int i = 0; i < 3; i++) color[i] = weight * incoming_color[i]*object_diffuse[i];

            }else{
                // evaluate specular

                GetDirectionUsingPhongBRDF(hit_point, object_normal, ray, object_shininess, incoming_direction);

                double weight = inv_probability_to_kill * inv_probability_for_specular * DotProduct(incoming_direction, object_normal);

                TracePathWithImportanceSampling(hit_point, incoming_direction, incoming_color, cubemap_id, cubemap_pixel_coord, weight_of_path*weight, depth+1);

                for(int i = 0; i < 3; i++) color[i] = weight * incoming_color[i]*object_specular[i];
            }


            return;
        }
    
    }else{ // ray does not hit any objects
        color[0] = color[1] = color[2] = 0.0;
        return;
    }

}



void OutputWeightIntoFile(double *eye, GLfloat *image, ofstream **path_tracing_data)
{
    static double deltaW = 512.0/window_width;
    static double deltaH = 512.0/window_height;

    for (int i = 0; i < window_height; i++){
        for (int j = 0; j < window_width; j++) {

            double ray[3] = { deltaW*(j+random0to1())-eye[0], deltaH*(i+random0to1())-eye[1], 512.0-eye[2] };
            double color[3] = { 0.0, 0.0, 0.0, };

            int cubemap_id = -1, cubemap_pixel_coord[2];

            Normalize(ray);

            TracePathWithImportanceSampling(eye, ray, color, cubemap_id, cubemap_pixel_coord, 1.0, 0);

            if(cubemap_id != -1){ // ray hits one of the cubemaps

                unsigned int pixel_color[3];
                GetPixelColor(cubemap_id, cubemap_pixel_coord, pixel_color, 0, 0);

                static float scale = 1.0/255.0;
                image[(i*window_width+j)*3]    += (GLfloat)(color[0]*pixel_color[0]*scale);
                image[(i*window_width+j)*3 +1] += (GLfloat)(color[1]*pixel_color[1]*scale);
                image[(i*window_width+j)*3 +2] += (GLfloat)(color[2]*pixel_color[2]*scale);


                unsigned int   row = i*window_width+j;
                unsigned short col = cubemap_pixel_coord[1]*CUBEMAP_RESOLUTION+cubemap_pixel_coord[0];
                float          weight[3];

                static double inv_N_SAMPLES = 1.0 / N_SAMPLES;

                for(int k = 0; k < 3; k++) weight[k] = color[k]*inv_N_SAMPLES;

                unsigned int   row_id       = row/N_ROWS_IN_BLOCK;
                unsigned short row_in_block = row%N_ROWS_IN_BLOCK;

                path_tracing_data[cubemap_id][row_id].write((char*)&row_in_block, sizeof(short));
                path_tracing_data[cubemap_id][row_id].write((char*)&col,          sizeof(short));
                path_tracing_data[cubemap_id][row_id].write((char*)weight,        sizeof(weight));
            }

        }
    }

}
