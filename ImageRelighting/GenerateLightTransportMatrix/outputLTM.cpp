
#include "monteCarloPathTracing.h"


void OutputLightTransportMatrices()
{
    // construct light transport matrix and output to file (18 light transport matrices (for each rgb component of each cube map))

    vector< vector<SparseMatrixEntry> > wavelet_matrix_transpose;

    Generate2DHaarMatrix(wavelet_matrix_transpose, CUBEMAP_RESOLUTION);
    TransposeSparseMatrix(wavelet_matrix_transpose);


    int n_blocks = (window_height*window_height) / N_ROWS_IN_BLOCK;

    for(int i = 0; i < 6; i++){ // for each cubemap

        for(int j = 0; j < 3; j++){ // for rgb component

            vector< vector<SparseMatrixEntry> > TWt(CUBEMAP_RESOLUTION*CUBEMAP_RESOLUTION);

            // process each row-block one by one
            for(int k = 0; k < n_blocks; k++){

                char filename[100];
                sprintf_s(filename, 100, "data/data_for_cubemap%d_%d", i, k);

                ifstream fin(filename, ios::binary);

                vector< vector<double> > row_block(N_ROWS_IN_BLOCK);

                for(int m = 0; m < N_ROWS_IN_BLOCK; m++)
                    row_block[m].resize(CUBEMAP_RESOLUTION*CUBEMAP_RESOLUTION, 0.0);


                int n_data_read = 0;

                while(true){

                    unsigned short row_in_block;
                    unsigned short col;

                    float weight[3];

                    fin.read((char*)&row_in_block, sizeof(short));
                    fin.read((char*)&col,          sizeof(short));
                    fin.read((char*)weight,        sizeof(weight));

                    if(fin.eof()) break;

                    row_block[row_in_block][col] += weight[j];

                    n_data_read++;
                }

                fin.close();

                cerr << "n_data_read " << n_data_read << endl;

                vector< vector<SparseMatrixEntry> > T_in_block(CUBEMAP_RESOLUTION*CUBEMAP_RESOLUTION);

                int n_entries = 0;

                for(int m = 0; m < N_ROWS_IN_BLOCK; m++){
                    for(int n = 0; n < CUBEMAP_RESOLUTION*CUBEMAP_RESOLUTION; n++){
                        if(fabs(row_block[m][n]) > 1.0e-12){
                            T_in_block[n].push_back( SparseMatrixEntry(m, row_block[m][n]) );
                            n_entries++;
                        }
                    }
                }

                cerr << "T n_entries " << n_entries << endl;
                n_entries = 0;

                vector< vector<SparseMatrixEntry> > TWt_in_block;

                SparseMatrixMultSparseMatrix(T_in_block, wavelet_matrix_transpose, TWt_in_block);
                                

                for(unsigned int m = 0; m < TWt_in_block.size(); m++){
                    for(unsigned int n = 0; n < TWt_in_block[m].size(); n++){
                        if(fabs(TWt_in_block[m][n].val) > 0.001){
                            TWt[m].push_back( SparseMatrixEntry(k*N_ROWS_IN_BLOCK+TWt_in_block[m][n].row, TWt_in_block[m][n].val) );
                            n_entries++;
                        }
                    }
                }

                cerr << "TWt n_entries " << n_entries << endl;

                if((k+1)%10 == 0)
                    cerr << "block " << k+1 << " / " << n_blocks << " done.\n";


            } // for(int k = 0; k < n_blocks; k++){


            int n_entries = 0;

            char filename_out[100];

            switch(j){
            case 0: sprintf_s(filename_out, 100, "data/ltm_%dr", i); break;
            case 1: sprintf_s(filename_out, 100, "data/ltm_%dg", i); break;
            case 2: sprintf_s(filename_out, 100, "data/ltm_%db", i); break;
            }   

            //// should be quantized to 8-bit before outputing to file???
            ofstream ltp_data(filename_out, ios::binary);

            for(unsigned int k = 0; k < TWt.size(); k++){
                unsigned int size = TWt[k].size();
                ltp_data.write((char*)&size, sizeof(int));

                for(unsigned int l = 0; l < TWt[k].size(); l++)
                    ltp_data.write((char*)&TWt[k][l], sizeof(SparseMatrixEntry));

                n_entries += size;
            }

            ltp_data.close();

            cerr << "output lpm " << i << " " << j << " n_entries " << n_entries << " \n";//

        } // for(int j = 0; j < 3; j++){

    } // for(int i = 0; i < 6; i++){



}
