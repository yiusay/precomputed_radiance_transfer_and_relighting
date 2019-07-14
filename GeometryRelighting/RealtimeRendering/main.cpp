
#include "monteCarloPathTracing.h"
#include "mesh.h"

#include <windows.h>
#include <gdiplus.h>
using namespace Gdiplus;
#pragma comment (lib,"Gdiplus.lib")

Mesh mesh;

GLint window_width = 512, window_height = 512;

float *png_data;
unsigned int png_width, png_height;

vector<Vertex> table_vertices;
int table_step = 2;



struct LightPixel {
    unsigned char  cubemap_id;
    unsigned short col;
};

struct LightPixelInfo {
    LightPixel light_pixel;
    float      weight;

    LightPixelInfo(int cubemap_id_in, int col_in, float w){
        light_pixel.cubemap_id = (unsigned char)cubemap_id_in;
        light_pixel.col        = (unsigned short)col_in;
        weight                 = w;
    }
};

bool operator>(const LightPixelInfo &a, const LightPixelInfo &b)
{
    return a.weight > b.weight;
}


vector<LightPixel> light_pixels;
unsigned int n_used_light;

vector< vector<SparseMatrixEntry> > TWt[6][3]; // for each cubemap's each color component



double  theta = 0.0, phi = 0.0;

GLfloat startx, starty;
GLfloat model_angle1 = 0.0, model_angle2 = 0.0, scale = 1.0, eye[3] = { 256.0, 256.0, 1228.0 };
bool    left_click = 0, right_click = 0;


bool isNewLighting = true;


void mouse(int button, int state, int x, int y)
{
    if ( button == GLUT_LEFT_BUTTON ) {
        if ( state == GLUT_DOWN ) {
            left_click = true;
            startx   = x;
            starty   = y;
        }else if (state == GLUT_UP) {
            left_click = false;
        }
    }else{ // button == GLUT_RIGHT_BUTTON
        if ( state == GLUT_DOWN ) {  
            right_click = true;
            startx   = x;
            starty   = y;
        }else if (state == GLUT_UP) {
            right_click = false;
        }
    }

}

void motion( int x, int y )
{ 
    if ( left_click && !right_click ) {       // rotating image
        phi   += (x - startx)*0.005;
        theta -= (y - starty)*0.005;

        if(phi < 0.0)         phi += 2.0*PI;
        else if(phi > 2.0*PI) phi -= 2.0*PI;

        if(theta < -PI)     theta = -PI;
        else if(theta > PI) theta = PI;

        isNewLighting = true;
    }else if( !left_click && right_click ) {       // rotating model
        model_angle1 += (x - startx);
   //     model_angle2 += (y - starty);
    }else{ // if( left_click && right_click ) // scaling
        scale -= (y - starty) * 0.01;
    }

    startx = x;
    starty = y;

    glutPostRedisplay();
}


void keyboard( unsigned char c, int x, int y )
{
    switch(c){
    case 'q':
        exit(0);
    case 'z':
        n_used_light = (unsigned int)Max(n_used_light*0.8, 1); 
        break;
    case 'x': 
        n_used_light = (unsigned int)Min(n_used_light/0.8, CUBEMAP_RESOLUTION*CUBEMAP_RESOLUTION*6);
        break;
    case 'a':
        if(n_used_light > 1) n_used_light--; 
        break;
    case 's': 
        if(n_used_light < CUBEMAP_RESOLUTION*CUBEMAP_RESOLUTION*6-1) n_used_light++; 
        break;
    }

    glutPostRedisplay();
}



void GetPixelColor(int cubemap_id, int *cubemap_pixel_coord, float *pixel_color, double view_theta, double view_phi)
{

    static double cubemap_theta[6][CUBEMAP_RESOLUTION][CUBEMAP_RESOLUTION];
    static double   cubemap_phi[6][CUBEMAP_RESOLUTION][CUBEMAP_RESOLUTION];

    static bool init = true;

    if(init){
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


    static float filtered_pixel_data[6][CUBEMAP_RESOLUTION][CUBEMAP_RESOLUTION][3];
    
    if(isNewLighting){    

      //  cerr << "isNewLighting " << endl;

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

                    for(int m = 0; m < 3; m++)  filtered_pixel_data[i][j][k][m] = png_data[pixel_index+m];
                 }
            }

        } // for(int i = 0; i < 6; i++){
        
        isNewLighting = false;

    } // if(isNewLighting){ 


    pixel_color[0] = filtered_pixel_data[cubemap_id][ cubemap_pixel_coord[0] ][ cubemap_pixel_coord[1] ][0];
    pixel_color[1] = filtered_pixel_data[cubemap_id][ cubemap_pixel_coord[0] ][ cubemap_pixel_coord[1] ][1];
    pixel_color[2] = filtered_pixel_data[cubemap_id][ cubemap_pixel_coord[0] ][ cubemap_pixel_coord[1] ][2];
}


void RenderCubemap()
{
    static double retrieve_theta[8], retrieve_phi[8];
    
    static vector<Vertex> sphere_vertices;
    static vector<Face>   sphere_faces;
    static vector< pair<float, float> > sphere_tex_coords;


    static GLuint texName;

    static bool init = true;

    if(init){

        double r =  768;
        int    div = 10;

        for(int i = -80; i <= 80; i+=div){   
            for(int j = 0; j <= 360; j+=div){
                double coord_in[3] = { r*cos(i/180.0*PI)*-cos(j/180.0*PI), r*sin(i/180.0*PI), r*cos(i/180.0*PI)*-sin(j/180.0*PI)};
                sphere_vertices.push_back( Vertex(coord_in, -1) );

                sphere_tex_coords.push_back( make_pair(j/360.0, (i+90.0)/180.0) );
            }
        }

        // south pole
        sphere_vertices.push_back( Vertex(0.0, -r, 0.0, -1) );
        sphere_tex_coords.push_back( make_pair(0.5, 0.0) );

        // north pole
        sphere_vertices.push_back( Vertex(0.0, r, 0.0, -1) );
        sphere_tex_coords.push_back( make_pair(0.5, 1.0) );

        for(int i = 0; i < 180/div; i++){
            for(int j = 0; j < 360/div; j++){
                if(i == 0)       sphere_faces.push_back( Face(sphere_vertices.size()-2, j+1, j, -1) );
                else if(i == 17) sphere_faces.push_back( Face((i-1)*37+j, (i-1)*37+j+1, sphere_vertices.size()-1, -1) );
                else{
                    sphere_faces.push_back( Face((i-1)*37+j+1,   i*37+j+1, i*37+j, -1) );
                    sphere_faces.push_back( Face((i-1)*37+j, (i-1)*37+j+1, i*37+j, -1) );
                }
            }
        }


        glGenTextures(1, &texName);
        glBindTexture(GL_TEXTURE_2D, texName);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, png_width, png_height, 0, GL_RGB, GL_FLOAT, png_data);

        init = false;
    }

    static int cube_face_id[6][4] = { {3,2,6,7}, {0,4,5,1}, {0,1,2,3}, {1,5,6,2}, {7,6,5,4}, {0,3,7,4}, };

    glEnable(GL_TEXTURE_2D);    
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

    glBindTexture(GL_TEXTURE_2D, texName);

    glPushMatrix();
    glTranslatef(256, 0, 256);
    glRotatef(theta*InvPI*180.0, 0, 0, 1);
    glRotatef(phi*InvPI*180.0,   0, 1, 0);


    for(unsigned int i = 0; i < sphere_faces.size(); i++){
        glBegin(GL_TRIANGLES);
        glTexCoord2f(sphere_tex_coords[ sphere_faces[i].v_id[0] ].first, sphere_tex_coords[ sphere_faces[i].v_id[0] ].second);
        glVertex3dv(sphere_vertices[ sphere_faces[i].v_id[0] ].point);
        glTexCoord2f(sphere_tex_coords[ sphere_faces[i].v_id[1] ].first, sphere_tex_coords[ sphere_faces[i].v_id[1] ].second);
        glVertex3dv(sphere_vertices[ sphere_faces[i].v_id[1] ].point);
        glTexCoord2f(sphere_tex_coords[ sphere_faces[i].v_id[2] ].first, sphere_tex_coords[ sphere_faces[i].v_id[2] ].second);
        glVertex3dv(sphere_vertices[ sphere_faces[i].v_id[2] ].point);
        glEnd();
     }

    glPopMatrix();
    glDisable(GL_TEXTURE_2D);  
}

void renderBitmapString(float x, float y, float z, void *font, char *string) 
{  
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glColor3d(1, 1, 1);
    glRasterPos3f(x, y, z);

    //float a[4];
    //glGetFloatv(GL_CURRENT_RASTER_POSITION, a);
    //cerr << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << endl;

    int len = (int) strlen(string);
    for (int i = 0; i < len; i++) {
        glutBitmapCharacter(font, string[i]);
    }

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

}






void display()
{
    static GLfloat *buffer;
    static int k = 0;
    static vector< vector<SparseMatrixEntry> > wavelet_matrix;

    static int total_n_vertices = mesh.vertices.size() + table_vertices.size();
    
    if(k == 0){ cerr << "total_n_vertices " << total_n_vertices << endl;
        buffer = new GLfloat [total_n_vertices*3];
     
        Generate2DHaarMatrix(wavelet_matrix, CUBEMAP_RESOLUTION);
    }

    k++;

    vector<double> wavelet_light[6][3];

    for(int i = 0; i < 6; i++){

        wavelet_light[i][0].resize(CUBEMAP_RESOLUTION*CUBEMAP_RESOLUTION, 0.0);
        wavelet_light[i][1].resize(CUBEMAP_RESOLUTION*CUBEMAP_RESOLUTION, 0.0);
        wavelet_light[i][2].resize(CUBEMAP_RESOLUTION*CUBEMAP_RESOLUTION, 0.0);

        for(int j = 0; j < CUBEMAP_RESOLUTION; j++){
            for(int k = 0; k < CUBEMAP_RESOLUTION; k++){

                int    pixel_coord[2] = { k, j };
                float  pixel_color[3];

                GetPixelColor(i, pixel_coord, pixel_color, theta, phi);

                int col = j*CUBEMAP_RESOLUTION+k;

                for(unsigned int l = 0; l < wavelet_matrix[col].size(); l++){
                    wavelet_light[i][0][wavelet_matrix[col][l].row] += wavelet_matrix[col][l].val * pixel_color[0];
                    wavelet_light[i][1][wavelet_matrix[col][l].row] += wavelet_matrix[col][l].val * pixel_color[1];
                    wavelet_light[i][2][wavelet_matrix[col][l].row] += wavelet_matrix[col][l].val * pixel_color[2];
                }

            }
        }

    } // for(int i = 0; i < 6; i++){

    // sort each light pixel based on area-weighted selection method (i.e. # of covering pixels*light_intensity_coefficient ) 
    vector<LightPixelInfo> light_pixel_info;

    for(int i = 0; i < 6; i++){

        for(int j = 0; j < CUBEMAP_RESOLUTION; j++){
            for(int k = 0; k < CUBEMAP_RESOLUTION; k++){

                int    pixel_coord[2] = { k, j };
                float  pixel_color[3];

                GetPixelColor(i, pixel_coord, pixel_color, theta, phi);

                float weight = 0;

                int col = j*CUBEMAP_RESOLUTION+k;

                weight += TWt[i][0][col].size() * fabs(wavelet_light[i][0][col]);
                weight += TWt[i][1][col].size() * fabs(wavelet_light[i][1][col]); 
                weight += TWt[i][2][col].size() * fabs(wavelet_light[i][2][col]);

                light_pixel_info.push_back( LightPixelInfo(i, col, weight) );
            }
        }

    }

    sort(light_pixel_info.begin(), light_pixel_info.end(), greater<LightPixelInfo>());

    for(int i = 0; i < total_n_vertices*3; i++)  buffer[i] = 0.0;

    int n_mult = 0;

    for(unsigned int i = 0; i < n_used_light; i++){

        LightPixel *light_pixel = &(light_pixel_info[i].light_pixel);

        for(int j = 0; j < 3; j++){
            for(unsigned int k = 0; k < TWt[light_pixel->cubemap_id][j][light_pixel->col].size(); k++){
                float val = TWt[light_pixel->cubemap_id][j][light_pixel->col][k].val;
                int   row = TWt[light_pixel->cubemap_id][j][light_pixel->col][k].row;

                buffer[row*3+j] += val * wavelet_light[light_pixel->cubemap_id][j][light_pixel->col];

                n_mult += 1;
            }
        }

    }


    glClearColor( 0.0, 0.0, 0.0, 0.0 );
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    gluLookAt(eye[0], eye[1], eye[2],
              eye[0], eye[1], 0.0,
                 0.0,    1.0,   0.0 );
    glTranslatef(256, 0, 256);
    glScalef(scale, scale, scale);
    glRotatef(model_angle1, 0, 1, 0);
    glRotatef(model_angle2, 1, 0, 0);
    glTranslatef(-256, 0, -256);


    for(unsigned int i = 0; i < mesh.faces.size(); i++){
        int *v_id = mesh.faces[i].v_id;

        glBegin(GL_TRIANGLES);
        glColor3f(buffer[v_id[0]*3], buffer[v_id[0]*3+1], buffer[v_id[0]*3+2]);
        glVertex3dv(mesh.vertices[ v_id[0] ].point);
        glColor3f(buffer[v_id[1]*3], buffer[v_id[1]*3+1], buffer[v_id[1]*3+2]);
        glVertex3dv(mesh.vertices[ v_id[1] ].point);
        glColor3f(buffer[v_id[2]*3], buffer[v_id[2]*3+1], buffer[v_id[2]*3+2]);
        glVertex3dv(mesh.vertices[ v_id[2] ].point);
        glEnd();
    }

    static int n_rows = (512/table_step)+1, n_cols = (512/table_step)+1;

    // render table
    for(int i = 0; i < n_cols-1; i++){
        for(int j = 0; j < n_rows-1; j++){
            
            unsigned int v_id[4] = { n_rows*i+j, n_rows*i+j+1, n_rows*(i+1)+j+1, n_rows*(i+1)+j };

            glBegin(GL_QUADS);
            glColor3f(buffer[(mesh.vertices.size()+v_id[0])*3]*0.8, 
                      buffer[(mesh.vertices.size()+v_id[0])*3+1], 
                      buffer[(mesh.vertices.size()+v_id[0])*3+2]*0.8);
            glVertex3dv(table_vertices[ v_id[0] ].point);
            glColor3f(buffer[(mesh.vertices.size()+v_id[1])*3]*0.8, 
                      buffer[(mesh.vertices.size()+v_id[1])*3+1], 
                      buffer[(mesh.vertices.size()+v_id[1])*3+2]*0.8);
            glVertex3dv(table_vertices[ v_id[1] ].point);
            glColor3f(buffer[(mesh.vertices.size()+v_id[2])*3]*0.8, 
                        buffer[(mesh.vertices.size()+v_id[2])*3+1], 
                        buffer[(mesh.vertices.size()+v_id[2])*3+2]*0.8);
            glVertex3dv(table_vertices[ v_id[2] ].point);
            glColor3f(buffer[(mesh.vertices.size()+v_id[3])*3]*0.8, 
                       buffer[(mesh.vertices.size()+v_id[3])*3+1], 
                       buffer[(mesh.vertices.size()+v_id[3])*3+2]*0.8);
            glVertex3dv(table_vertices[ v_id[3] ].point);
            glEnd();
        }
    }

    RenderCubemap();


    static double time[10];
    static int ct = 0;

    char info1[256], info2[256], info3[256];

    sprintf(info1, "FPS: %lf", 1.0/(clock()-time[ct%10])*10.0 * CLOCKS_PER_SEC);
    sprintf(info2, "#lights: %d / %d", n_used_light, CUBEMAP_RESOLUTION*CUBEMAP_RESOLUTION*6);
    sprintf(info3, "#multiplications: %d", n_mult);

    renderBitmapString(-0.95, 0.9, 0, GLUT_BITMAP_HELVETICA_18 ,  info1);
    renderBitmapString(-0.95, 0.8, 0, GLUT_BITMAP_HELVETICA_18 , info2);
    renderBitmapString(-0.95, 0.7, 0, GLUT_BITMAP_HELVETICA_18,  info3);

    time[ct%10] = clock();
    ct++;

    glutSwapBuffers();
    glutPostRedisplay();
 }


int main(int argc, char *argv[])
{
    GdiplusStartupInput gdiplusStartupInput;
    ULONG_PTR           gdiplusToken;

    // Initialize GDI+.
    GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);

    png_data      = LoadPngData(L"stars.png", png_width, png_height);

    GdiplusShutdown(gdiplusToken);
    
    if( mesh.ReadOFFFile("../../hw1/models/rocker-arm.off") == false ) return false;


    //////////////////////
    // generate table
    /////////////////////

    for(int i = 0; i <= 512; i += table_step){
        for(int j = 0; j <= 512; j += table_step){
            double coord_in[3] = {i, 0, j};
            table_vertices.push_back( Vertex(coord_in, -1) );
        }
    }


    _setmaxstdio(2048);
    cerr << "max number of files that can be used " << _getmaxstdio() << endl;



    for(int i = 0; i < 6; i++){ // for each cubemap

        for(int j = 0; j < 3; j++){ // for each color component

            cerr << "reading " << i << " " << j << endl;

            TWt[i][j].resize(CUBEMAP_RESOLUTION*CUBEMAP_RESOLUTION);

            char filename[100];

            switch(j){
            case 0: sprintf_s(filename, 100, "data/ltm_%dr", i);  break;
            case 1: sprintf_s(filename, 100, "data/ltm_%dg", i);  break;
            case 2: sprintf_s(filename, 100, "data/ltm_%db", i);  break;
            }

            ifstream fin(filename, ios::binary);

            int col = 0;

            while(true){
                int n_rows_in_column;
                fin.read((char*)&n_rows_in_column, sizeof(int));

                if(fin.eof()) break;

                for(int k = 0; k < n_rows_in_column; k++){
                    SparseMatrixEntry data;
                    fin.read((char*)&data, sizeof(SparseMatrixEntry));

                    //  if(fabs(data.val) > 0.0019) 
                    TWt[i][j][col].push_back(data);
                }

                col++;
            } // while(true){

        } // for(int j = 0; j < 3; j++){

    }  // for(int i = 0; i < 6; i++){


    // sort each light pixel based on their covering area i.e. number of pixles 
    vector<LightPixelInfo> light_pixel_info;

    for(int i = 0; i < 6; i++){

        for(int j = 0; j < CUBEMAP_RESOLUTION*CUBEMAP_RESOLUTION; j++){

            unsigned int n_covering_pixels = 0;

            n_covering_pixels += TWt[i][0][j].size(); 
            n_covering_pixels += TWt[i][1][j].size(); 
            n_covering_pixels += TWt[i][2][j].size();       

            light_pixel_info.push_back( LightPixelInfo(i, j, n_covering_pixels) );
        }

    }

    sort(light_pixel_info.begin(), light_pixel_info.end(), greater<LightPixelInfo>());

    for(unsigned int i = 0; i < light_pixel_info.size(); i++) light_pixels.push_back(light_pixel_info[i].light_pixel);

    n_used_light = light_pixels.size();


 
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowPosition(300, 20);
    glutInitWindowSize(window_width, window_height);
    glutCreateWindow(argv[0]);
    glutDisplayFunc(display);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    GLInit();
    glutMainLoop();

    return 1;
}

