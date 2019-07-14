
#include "monteCarloPathTracing.h"


bool operator<(const SparseMatrixEntry &a, const SparseMatrixEntry &b)
{
    return a.row < b.row;
}
  


inline void RefreshSparseColumn(vector<SparseMatrixEntry> &column)
{
    sort(column.begin(), column.end());

    vector<SparseMatrixEntry> refreshed_column;
    unsigned int i = 0;

    while(i < column.size()){
        double sum_entry = 0.0;
        int    current_row = column[i].row;
       
        while(current_row == column[i].row){
            sum_entry += column[i].val;
            i++;
        }

        if(fabs(sum_entry) > 1.0e-12) 
            refreshed_column.push_back( SparseMatrixEntry(current_row, sum_entry) );
    }

    column.swap(refreshed_column); 
}



inline void SparseMatrixMultSparseMatrix(vector< vector<SparseMatrixEntry> > &a, 
                                         vector< vector<SparseMatrixEntry> > &b, 
                                         vector< vector<SparseMatrixEntry> > &res)
{
    res.resize(b.size());

    for(unsigned int i = 0; i < b.size(); i++){
        for(unsigned int j = 0; j < b[i].size(); j++){

            for(unsigned int k = 0; k < a[ b[i][j].row ].size(); k++){
                res[i].push_back( SparseMatrixEntry(a[ b[i][j].row ][k].row, a[ b[i][j].row ][k].val*b[i][j].val) );
            }

        }

        // merge entries with the same row into one 
        RefreshSparseColumn(res[i]);
    }

}


inline void TransposeSparseMatrix(vector< vector<SparseMatrixEntry> > &a)
{
    vector< vector<SparseMatrixEntry> > res(a.size());

    for(unsigned int i = 0; i < a.size(); i++){
        for(unsigned int j = 0; j < a[i].size(); j++)
            res[ a[i][j].row ].push_back( SparseMatrixEntry(i, a[i][j].val) );
    }

     a.swap(res); 
}


void Generate2DHaarMatrix(vector< vector<SparseMatrixEntry> > &wavelet_matrix, int image_dim)
{
    // assuming that 2d image is represented by 1d vector 
    // i.e. i-th entry of the vector is at {i%IMAGE_RESOLUTION,  (int)(i/IMAGE_RESOLUTION)} in the corresponding 2d image
   

    // initialize wavelet_matrix as identity matrix
    wavelet_matrix.resize(image_dim*image_dim);

    for(int i = 0; i < image_dim*image_dim; i++){
        wavelet_matrix[i].push_back( SparseMatrixEntry(i, 1.0) );
    }


    int conversion_dim = image_dim;

    while(conversion_dim > 1){

        if(conversion_dim%2 == 1){
            cerr << "dim must be a power of 2\n"; 
            exit(-1);
        }

        vector< vector<SparseMatrixEntry> > Harr_2d(image_dim*image_dim);

        // fill unchanged part as identity matrix (only [conversion_dim x conversion_dim] part will be changed)
        for(int i = 0; i < conversion_dim; i++){ 
            for(int j = conversion_dim; j < image_dim; j++){ 
                Harr_2d[i*image_dim+j].push_back( SparseMatrixEntry(i*image_dim+j, 1.0) );
            }
        }
        for(int i = conversion_dim; i < image_dim; i++){ 
            for(int j = 0; j < image_dim; j++){ 
                Harr_2d[i*image_dim+j].push_back( SparseMatrixEntry(i*image_dim+j, 1.0) );
            }
        }


        for(int i = 0; i < conversion_dim; i++){
            for(int j = 0; j < conversion_dim; j++){

                // compute contribution of  [conversion_dim x conversion_dim] part@ described as 1d vector

                int contribution[4];

                contribution[0] = 1;
                if(i%2 == 0){
                    contribution[2] = 1;
                    if(j%2 == 0) contribution[1] = contribution[3] =  1;
                    else         contribution[1] = contribution[3] = -1;
                }else{
                    contribution[2] = -1;
                    if(j%2 == 0){
                        contribution[1] =  1; 
                        contribution[3] = -1;
                    }else{
                        contribution[1] = -1; 
                        contribution[3] =  1;
                    }
                }

                int base_row = (int)(i/2),  base_col = (int)(j/2);
                int jump     = (int)(conversion_dim/2);

                int index[4] = {        base_row*image_dim + base_col, 
                                        base_row*image_dim + (base_col+jump),
                                 (base_row+jump)*image_dim + base_col, 
                                 (base_row+jump)*image_dim + (base_col+jump) };  
               
                for(int k = 0; k < 4; k++){
                    Harr_2d[i*image_dim+j].push_back( SparseMatrixEntry(index[k], contribution[k]*0.5) ); 
                }

            }
        }


        vector< vector<SparseMatrixEntry> > mult_result;

        SparseMatrixMultSparseMatrix(Harr_2d, wavelet_matrix, mult_result);

        wavelet_matrix.swap(mult_result);


        conversion_dim /= 2;
        
        static int ct = 0;
   //     if(++ct == 3) break;
    }

}
