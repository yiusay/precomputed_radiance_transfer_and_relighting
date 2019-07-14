
#include "monteCarloPathTracing.h"
#include "mesh.h"
#include "kdtree.h"

#include <windows.h>
#include <gdiplus.h>
using namespace Gdiplus;
#pragma comment (lib,"Gdiplus.lib")

Mesh   mesh;
Kdtree kdtree;

GLint window_width = 512, window_height = 512;

unsigned char *png_data;
unsigned int   png_width, png_height;

double *table_texture;
int     table_texture_width, table_texture_height;

vector<Vertex> table_vertices;
int table_step = 2;

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
        model_angle2 += (y - starty);
    }else{ // if( left_click && right_click ) // scaling
        scale -= (y - starty) * 0.01;
    }

    startx = x;
    starty = y;

    glutPostRedisplay();
}






void display()
{
    glClearColor( 0.0, 0.0, 0.0, 0.0 );
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    gluLookAt(eye[0], eye[1], eye[2],
              eye[0], eye[1],   0.0,
                 0.0,    1.0,   0.0 );
    glTranslatef(256, 0, 256);
    glScalef(scale, scale, scale);
    glRotatef(model_angle1, 0, 1, 0);
    glRotatef(model_angle2, 1, 0, 0);
    glTranslatef(-256, 0, -256);


    static GLfloat *image, *buffer;
    static int k = 0;

    static ofstream *path_tracing_data[6];

    static int total_n_vertices = mesh.vertices.size() + table_vertices.size();

    if(k == 0){
        // initialize rendering target
        image = new GLfloat [total_n_vertices*3];
        for(unsigned int i = 0; i < total_n_vertices*3; i++) image[i] = 0.0;

        buffer = new GLfloat [total_n_vertices*3];


        int n_blocks = (int)ceil(total_n_vertices / (double)N_ROWS_IN_BLOCK);

        for(int i = 0; i < 6; i++){
            path_tracing_data[i] = new ofstream [n_blocks];

            for(int j = 0; j < n_blocks; j++){
                char filename[100];
                sprintf_s(filename, 100, "data/data_for_cubemap%d_%d", i, j);

                path_tracing_data[i][j].open(filename, ios::binary);

                if(path_tracing_data[i][j].is_open() == false){
                    cerr << "fail to open: " << filename << endl;
                    exit(1);
                }
            }
        }

    } // if(k == 0){

    srand(k++);

    cerr << "# samples: " << k << endl;

    OutputWeightIntoFile(image, path_tracing_data);


    double inv_k = 1.0 / (double)(k); // # of samples should be doubled if we render w/ direct lighting
    for(unsigned int i = 0; i < total_n_vertices*3; i++) buffer[i] = image[i] * inv_k;

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

    for(unsigned int i = 0; i < n_cols-1; i++){
        for(unsigned int j = 0; j < n_rows-1; j++){
            
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


    glutSwapBuffers();

    if(k < N_SAMPLES){
        glutPostRedisplay();
    }else{
        int n_blocks = (int)ceil(total_n_vertices / (double)N_ROWS_IN_BLOCK);

        for(int i = 0; i < 6; i++){
            for(int j = 0; j < n_blocks; j++) path_tracing_data[i][j].close();
            delete [] path_tracing_data[i];
        }

        OutputLightTransportMatrices();

        exit(0);
    }
}


int main(int argc, char *argv[])
{
    GdiplusStartupInput gdiplusStartupInput;
    ULONG_PTR           gdiplusToken;

    // Initialize GDI+.
    GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);

    png_data      = LoadPngData(L"stars.png", png_width, png_height);

    GdiplusShutdown(gdiplusToken);


    if( mesh.ReadOFFFile("cow.off") == false ) return false;

    kdtree.InitKdtree(mesh.vertices, mesh.faces);


    //////////////////////
    // generate table
    /////////////////////

    for(int i = 0; i <= 512; i += table_step){
        for(int j = 0; j <= 512; j += table_step){
            double coord_in[3] = {i, 0, j};
            table_vertices.push_back( Vertex(coord_in, -1) );
        }
    }

    for(unsigned int i = 0; i < table_vertices.size(); i++){
        table_vertices[i].normal[0] = 0.0;
        table_vertices[i].normal[1] = 1.0;
        table_vertices[i].normal[2] = 0.0;
    }


    //if((window_height*window_height)%N_ROWS_IN_BLOCK != 0){
    //    cerr << "[(window_height*window_height)%N_ROWS_IN_BLOCK] must be 0.\n";
    //    exit(1);
    //}

    _setmaxstdio(2048);
    cerr << "max number of files that can be used " << _getmaxstdio() << endl;

    if(_getmaxstdio() < ceil((mesh.vertices.size()+table_vertices.size())/(double)N_ROWS_IN_BLOCK)*18){
        cerr << ceil((mesh.vertices.size()+table_vertices.size())/(double)N_ROWS_IN_BLOCK)*18 <<  " files must be able to be opened at the same time.\n";
        exit(1);
    }


 
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL);
    glutInitWindowPosition(300, 20);
    glutInitWindowSize(window_width, window_height);
    glutCreateWindow(argv[0]);
    glutDisplayFunc(display);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    GLInit();
    glutMainLoop();

    return 1;
}

