
#include "monteCarloPathTracing.h"

#include <windows.h>
#include <gdiplus.h>
using namespace Gdiplus;
#pragma comment (lib,"Gdiplus.lib")



GLint window_width = 512, window_height = 512;

unsigned char *png_data;
unsigned int   png_width, png_height;

double  theta = 0.0, phi = 0.0;
GLfloat startx, starty;
bool    left_click = 0, right_click = 0;

double eye[3] = { 256.0, 256.0, 1228.0 };

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
    if ( left_click && !right_click ) {       // rotating
        phi   += (x - startx)*0.005;
        theta -= (y - starty)*0.005;

        if(phi < 0.0)         phi += 2.0*PI;
        else if(phi > 2.0*PI) phi -= 2.0*PI;

        if(theta < -PI)     theta = -PI;
        else if(theta > PI) theta = PI;

        isNewLighting = true;
    }

    startx = x;
    starty = y;

    glutPostRedisplay();
}





void display()
{
    glClearColor( 0.0, 1.0, 0.0, 0.0 );
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    static GLfloat *image, *buffer;
    static int k = 0;

    static ofstream *path_tracing_data[6];

    if(k == 0){
        // initialize rendering target
        image = new GLfloat [window_width*window_height*3];
        for(int i = 0; i < window_width*window_height*3; i++) image[i] = 0.0;

        buffer = new GLfloat [window_width*window_height*3];


        int n_blocks = (window_height*window_height) / N_ROWS_IN_BLOCK;

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

    OutputWeightIntoFile(eye, image, path_tracing_data);


    double inv_k = 1.0 / (double)(k); // # of samples should be doubled if we render w/ direct lighting
    for(int i = 0; i < window_width*window_height*3; i++) buffer[i] = image[i] * inv_k;

    glDrawPixels(window_width, window_height, GL_RGB, GL_FLOAT, buffer);


    glutSwapBuffers();
 
    if(k < N_SAMPLES){
        glutPostRedisplay();
    }else{
        int n_blocks = (window_height*window_height) / N_ROWS_IN_BLOCK;

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

    png_data      = LoadPngData(L"area_light1.png", png_width, png_height);

    GdiplusShutdown(gdiplusToken);

    if((window_height*window_height)%N_ROWS_IN_BLOCK != 0){
        cerr << "[(window_height*window_height)%N_ROWS_IN_BLOCK] must be 0.\n";
        exit(1);
    }

    _setmaxstdio(2048);
    cerr << "max number of files that can be used " << _getmaxstdio() << endl;

    if(_getmaxstdio() < (window_height*window_height)/N_ROWS_IN_BLOCK*3){
        cerr << (window_height*window_height)/N_ROWS_IN_BLOCK*3 <<  " files must be able to be opened at the same time.\n";
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
   glutMainLoop();

   return 1;
}

