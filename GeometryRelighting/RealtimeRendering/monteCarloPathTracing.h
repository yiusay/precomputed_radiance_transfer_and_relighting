#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <vector>
#include <algorithm>
#include <list>
#include <stack>
#include <utility>
#include <GL/glut.h>
using namespace std;

#include <windows.h>
#include <gdiplus.h>
using namespace Gdiplus;
#pragma comment (lib,"Gdiplus.lib")

class Mesh;
class Kdtree;

struct SparseMatrixEntry {
    int   row;
    float val;

    SparseMatrixEntry(){};
    SparseMatrixEntry(int row_in, float val_in){
        row = row_in;
        val = val_in;
    }

    bool operator<(const int &a){
        return a < row;
    }
    bool operator==(const int &a){
        return a == row;
    }
};

struct SparseMatrixEntry1Byte {
    int  row;
    char val;

    SparseMatrixEntry1Byte(){};
    SparseMatrixEntry1Byte(int row_in, char val_in){
        row = row_in;
        val = val_in;
    }
    bool operator<(const int &a){
        return a < row;
    }
    bool operator==(const int &a){
        return a == row;
    }
};



extern void GLInit();
extern void CrossProduct( double *a, double *b, double *c );
extern void Normalize(double *a);
extern double DotProduct(double *a, double *b);
extern double DotProduct4D(double *a, double *b);
extern void GetEdgeVector(double *v1, double *v2, double *edge_vector);
extern void GetArea(double *normal, double &area);
extern double GetDistance(double *a, double *b);
extern double GetLength(double *a);
extern void Swap(double &a, double &b);
extern double Max(double a, double b);
extern double Min(double a, double b);
extern bool SolveLinearSystem(double (*matrix)[4], double *rhs, double *solution);
extern void RotateAroundAxis(double *rotation_axis, double theta, double *vec_in, double *vec_out);
extern void ComputeReflectionVector(double *normal, double *vector_in, double *reflection_vector);
extern float* LoadPngData(WCHAR *filename, unsigned int &img_width, unsigned int &img_height);


extern void SparseMatrixMultSparseMatrix(vector< vector<SparseMatrixEntry> > &a, 
                                         vector< vector<SparseMatrixEntry> > &b, 
                                         vector< vector<SparseMatrixEntry> > &res);
extern void TransposeSparseMatrix(vector< vector<SparseMatrixEntry> > &a);

extern void OutputWeightIntoFile(GLfloat *image, ofstream **path_tracing_data);
extern void OutputLightTransportMatrices();
extern void Generate2DHaarMatrix(vector< vector<SparseMatrixEntry> > &wavelet_matrix, int image_dim);


extern GLint window_width, window_height;
extern float *png_data;
extern unsigned int   png_width, png_height;


extern bool isNewLighting;

extern Mesh   mesh;


#define CUBEMAP_RESOLUTION 64
#define PI    3.14159265
#define TwoPI 6.28318531
#define InvPI 0.318309886


#define EPSILON 1.0e-6