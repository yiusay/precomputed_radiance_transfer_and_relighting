
#include "monteCarloPathTracing.h"

inline void CrossProduct(double *a, double *b, double *c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

inline void Normalize(double *a)
{
    double length = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

    if(!length) return;

    double inv_length = 1.0/length;

    a[0] *= inv_length;
    a[1] *= inv_length;
    a[2] *= inv_length;
}

inline double DotProduct(double *a, double *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double DotProduct4D(double *a, double *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
}

inline void GetArea(double *normal, double &area)
{
    // assuming that "normal" is not normalized at this point 
    area = 0.5 * sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
}


inline double GetDistance(double *a, double *b)
{
    return sqrt( pow(a[0]-b[0], 2.0) + pow(a[1]-b[1], 2.0) + pow(a[2]-b[2], 2.0) ); 
}

inline double GetLength(double *a)
{
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]); 
}


inline void Swap(double &a, double &b)
{
    double t = a;  a = b;  b = t;
}

inline double Max(double a, double b)
{
    if(a > b) return a;
    else      return b;
}

inline double Min(double a, double b)
{
    if(a < b) return a;
    else      return b;
}


inline bool SolveLinearSystem(double (*matrix)[4], double *rhs, double *solution)
{  
    // perform gaussian elimination

    for(int i = 0; i <= 3; i++){
        double maxAbsoluteValue = -1.0;
        int    pivot_index;

        // choose maximum "matrix[j][i]" as a pivot
        for(int j = i; j <= 3; j++){
            if(fabs(matrix[j][i]) > maxAbsoluteValue){
                maxAbsoluteValue = fabs(matrix[j][i]);
                pivot_index      = j;
            }
        }

        // matrix is singular
        if(maxAbsoluteValue < 1.0e-6) return false;

        for(int j = i; j <= 3; j++) Swap(matrix[i][j], matrix[pivot_index][j]);
        Swap(rhs[i], rhs[pivot_index]);

        double scale = 1.0 / matrix[i][i];

        for(int j = i+1; j <= 3; j++){
            double pivot = -matrix[j][i]*scale;

            for(int k = 0;   k <= i;   k++) matrix[j][k] = 0.0;
            for(int k = i+1; k <= 3; k++){
                if(fabs(matrix[i][k]) > 1.0e-6)
                    matrix[j][k] += matrix[i][k] * pivot;
                else
                    break;
            }

            rhs[j] += rhs[i] * pivot;
        }
    }

   

//    cerr << "diag " << matrix[0][0] << " " << matrix[1][1] << " "  << matrix[2][2] << " " <<  matrix[3][3] << endl;
    

//    if(fabs(matrix[3][3]) < 1.0e-6) return false;

    for(int i = 3; i >= 0; i--){  
        solution[i] = 0.0;
        for(int j = i+1; j <= 3; j++) solution[i] += solution[j]*matrix[i][j];

        solution[i] = (rhs[i]-solution[i])/matrix[i][i];
    }


    return true;
}

inline void ComputeReflectionVector(double *normal, double *vector_in, double *reflection_vector)
{
    double scale = 2.0 * fabs( DotProduct(vector_in, normal) );

    for(int i = 0; i < 3; i++) reflection_vector[i] = vector_in[i] + scale*normal[i];
    Normalize(reflection_vector);
}


void RotateAroundAxis(double *rotation_axis, double theta, double *vec_in, double *vec_out)
{
    // rotate "vec_in" around "rotation_axis" using quarternion

    Normalize(rotation_axis);

    double common_factor = sin(theta*0.5);
    
    double a = cos(theta*0.5);
    double b = rotation_axis[0] * common_factor;
    double c = rotation_axis[1] * common_factor;
    double d = rotation_axis[2] * common_factor;
 
    double mat[9] = { a*a+b*b-c*c-d*d,     2*(b*c-a*d),      2*(b*d+a*c),
                          2*(b*c+a*d), a*a-b*b+c*c-d*d,      2*(c*d-a*b),
                          2*(b*d-a*c),     2*(c*d+a*b),  a*a-b*b-c*c+d*d };

    for(int i = 0; i < 3; i++) vec_out[i] = 0.0;

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++) vec_out[i] += mat[i*3+j] * vec_in[j];
    }

}


bool operator<(SparseMatrixEntry &a, SparseMatrixEntry &b)
{
    return a.row < b.row;
}


inline void AddElementToColumnBasedSparseMatrix(vector<SparseMatrixEntry> &col, int row, float val)
{
   //vector<SparseMatrixEntry>::iterator p = lower_bound(col.begin(), col.end(), row); 

   //if(p != col.end() && p->row == row) p->val += val;
   //else                                col.insert(p, SparseMatrixEntry(row, val));

    col.push_back( SparseMatrixEntry(row, val) );
   if(col.size() > 2000000) cerr << "koeta!\n";
}

    

unsigned char* LoadPngData(WCHAR *filename, unsigned int &img_width, unsigned int &img_height)
{
    Bitmap bitmap(filename);

    img_width  = bitmap.GetWidth();
    img_height = bitmap.GetHeight();

    cerr << img_width << " " << img_height << endl;

    unsigned char *pixel_data = new unsigned char [img_width*img_height*3];

    for(unsigned int i = 0; i < img_height; i++){
        for(unsigned int j = 0; j < img_width; j++){
            Color pixelColor;
            bitmap.GetPixel(j, i, &pixelColor);

            pixel_data[((img_height-1-i)*img_width+j)*3+0] = pixelColor.GetR();
            pixel_data[((img_height-1-i)*img_width+j)*3+1] = pixelColor.GetG();
            pixel_data[((img_height-1-i)*img_width+j)*3+2] = pixelColor.GetB();
        }
    }

    return pixel_data;
}
