// #include <fstream>
// #include <iostream>
#include "matrix_reader.h"

//using namespace std;

int main(int argc,char** argv){
	printf("./exe Filename\n");
	if (argc != 2) {printf("wrong input!\n");exit(-1);}

	const char* filename = argv[1];
	FILE *fp = fopen(filename,"r");
	idx_t m,n,nonz;
	idx_t* begin_pos;//rowind;
	idx_t* csr;//colptr;
	double *nzval;    
  
	dreadMM_dist(fp, &m, &n, &nonz, &nzval, &csr, &begin_pos);
  int nb;
  double upper[n][n], lower[n][n],mat[n][n];
        
//   for(int outer_row = nb+(nb*num); outer_row<n; outer_row+=nb){
//     for(int outer_col = nb+(nb*num); outer_col<n; outer_col+=nb){
//       for(int row = outer_row; row<outer_row+nb; row++){
//         for(int col = outer_col; col<outer_col+nb; col++){
//           printf("%f ", mat[row][col]);
//         }
//         printf("\n");
//       }
//       printf("\n\n");
//     }
//   }
// }
 
	 for (idx_t i=0; i< n; i++)
   {
     printf("\nNeighbor of vertex %d  begin_pos[%d]:%d  begin_pos[%d]:%d\t",i,i,begin_pos[i],i+1,begin_pos[i+1]);
     for (idx_t j=begin_pos[i]; j< begin_pos[i+1]; j++)
     {
       printf("%d %lf\t",csr[j],nzval[j]);
     }
   }

int col = 0;
_Bool zero = 1; // to check if element is 0
      for(int k = 0; k<n; k++){ //getting elements from row k
          for (idx_t i=0; i< n; i++){ 
            for (idx_t j=begin_pos[i]; j< begin_pos[i+1]; j++) {
              if(csr[j] == k){
                    zero = 0;
                    mat[k][col] = nzval[j];
                    printf("%f ", mat[k][col]);
                    col++; if (col == n) col = 0;
                    j = begin_pos[i+1];
                }  
              }
             if(zero == 1) {
                printf("%d ", 0); 
                upper[k][col] = 0;
                col++; if (col == n) col = 0;
             }
              zero = 1;
          }
               
      }

printf("\n\n MATRIX \n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            
            printf("%f ", mat[i][j]);
        }
        printf("\n");
    }

  printf("How many blocks are there?");
  scanf("%d", &nb);
  printf("%d", nb);

  struct total_blocks {
      int sb_row;
      int sb_col;
      double* nz;
  };
  int row_index = 0;
  int col_index = 0;
  printf("total_blocks is %I64d bytes long ", sizeof(struct total_blocks));
  struct total_blocks* block;
  block =  (struct total_blocks*)malloc(nb * sizeof(struct total_blocks));
  printf("block size: %d", sizeof(block));
    for(int i = 0; i<nb; i++){
        // printf("%d ", i);
        block[i].sb_row = sqrt((n*n)/nb);
        block[i].sb_col = sqrt((n*n)/nb);
        printf("size: %d", block[i].sb_col);
        //block[i] =  (double*)malloc(block[i].sb_row*block[i].sb_col * sizeof(double));
        block[i].nz =  (double*)malloc(block[i].sb_row*block[i].sb_col * sizeof(double));
            if(col_index == n){
                col_index = 0;
                row_index += block[i].sb_row;
            }
        for(int inner_row = 0; inner_row<block[i].sb_row; inner_row++){
            for(int inner_col = 0; inner_col<block[i].sb_col; inner_col++){
                block[i].nz[inner_row*block[i].sb_col+inner_col] = mat[row_index+inner_row][col_index+inner_col];
                //printf("%f ", block[i].nz[inner_row*block[i].sb_col+inner_col]);
            }
        }
        col_index+=block[i].sb_row;
        printf("\nNZ for this block %d:", i);
        for(int j = 0; j<block[i].sb_row*block[i].sb_col; j++){
            //printf("%f ", block[i].nz[j]);
        }
        //printf("\n\n");
    }
    printf("block side %d \n",block[0].sb_col);
int total_row = sqrt(nb);
int s_index = 0; 
int iteration = 0;
double upper_block[4][4];
double lower_block[4][4];

for(int row = 0; row<sqrt((n*n)/nb); row++ ){
    for(int col = 0; col<sqrt((n*n)/nb); col++){
        lower_block[row][col]=0;
        upper_block[row][col]=0;
        if(row == col) lower_block[row][col]=1;
    }
}


printf("total row: %d ",total_row);

for(int s_index = 0; s_index<nb-1; s_index+=total_row+1){
    for(int row = 0; row<block[0].sb_col; row++ ){
        for(int col = 0; col<block[0].sb_col; col++){
            lower[row][col]=0;
            upper[row][col]=0;
            if(row == col) lower[row][col]=1;
        }
    }
    //Category 1 -- STRUCT
    printf("Category 1 index: %d \n", s_index);
    for(int k = 0; k< block[s_index].sb_row; k++){
        for(int i = k; i<block[s_index].sb_row; i++){
            upper[k][i] = block[s_index].nz[k*block[s_index].sb_row+i];
        }
        for(int i = k+1; i<block[s_index].sb_row; i++){
            printf("LOWER IK %f = BLOCK IK %f / UPPER KK %f", lower_block[i][k], (block[s_index].nz[i*block[s_index].sb_row+k]), upper[k][k]);
            lower[i][k] = (block[s_index].nz[i*block[s_index].sb_row+k])/upper[k][k];
             
        }
        for(int j = k+1; j<block[s_index].sb_row; j++){
            for(int i = k+1; i<block[s_index].sb_row; i++){
            block[s_index].nz[i*block[s_index].sb_row+j] -= (lower[i][k]*upper[k][j]);
            printf("new block %d value row %d col %d: %f \n", s_index, i, j, block[s_index].nz[(i*(block[i].sb_col))+j]);

            }
        }

    }
    //adding updated category 1 elements back to block
    for(int k = 0; k< block[s_index].sb_row; k++){
        for(int i = k; i<block[s_index].sb_row; i++){
            block[s_index].nz[k*block[s_index].sb_row+i] = upper[k][i];
        }
        for(int i = k+1; i<block[s_index].sb_row; i++){
            (block[s_index].nz[i*block[s_index].sb_row+k]) = lower[i][k];
             
        }
    }

    //Category 2 -- STRUCT
    printf("\n\n");
    for(int block_index = s_index+total_row; block_index<nb; block_index+=(total_row)){
        printf("Category 2 index: %d \n", block_index);
        for(int k = 0; k<(block[block_index].sb_row); k++){
            for(int i = 0; i<(block[block_index].sb_row); i++){
                block[block_index].nz[i*block[block_index].sb_col+k]/=block[s_index].nz[k*(block[block_index].sb_col)+k];// mat[i][k]/=mat[k][k];
            }
            for(int j = k+1; j<(block[block_index].sb_row);j++){
                for(int i = 0; i <(block[block_index].sb_row);i++){
                    block[block_index].nz[i*(block[block_index].sb_col)+j]-= ((block[block_index].nz[i*(block[block_index].sb_col)+k])*block[s_index].nz[k*block[block_index].sb_col+j]); //mat[i][j]-= (mat[i][k]*mat[k][j]);
                
                }
            } 
            
        }     
    }
    //Category 3 -- STRUCT
    for(int block_index = s_index+1; block_index<(total_row*(iteration+1)); block_index++){
        printf("Category 3 index: %d \n", block_index);
        for(int k = 0; k<block[block_index].sb_row; k++){
            for(int j = 0; j<block[block_index].sb_row; j++){
                for(int i = k+1;i<block[block_index].sb_row; i++){
                    block[block_index].nz[i*block[block_index].sb_col+j]-= ((block[block_index].nz[k*(block[block_index].sb_col)+j])*block[s_index].nz[i*(block[block_index].sb_col)+k]);
                }
            }
        }
    }

    //Category 4 -- STRUCT
    for(int block_index = s_index+total_row+1; block_index<nb; block_index++){
        while((block_index%total_row) <= 0+iteration) block_index++;
        printf("Category 4 index: %d \n", block_index);
        printf("current block: %d \n", block_index);
        printf("left block: %d \n", block_index-(block_index%total_row)+iteration);
        printf("upper block: %d \n\n\n", ((block_index%total_row)+(iteration*total_row)));
            for(int k = 0; k< block[block_index].sb_col; k++){
                for(int j = 0; j< block[block_index].sb_col; j++){
                    for(int i = 0; i< block[block_index].sb_col; i++){
                        if(iteration == 2){
                            printf("current block %d value row %d col %d: %f \n", block_index, i, j, block[block_index].nz[(i*(block[i].sb_col))+j]);
                            printf("- left index %d %d value %f * upper %d %d value %f \n", i, k, block[block_index-(block_index%total_row)+iteration].nz[i*block[i].sb_row+k], k, j, block[(block_index%total_row)+(iteration*total_row)].nz[k*block[i].sb_row+j]);
                        }
                        block[block_index].nz[(i*(block[i].sb_col))+j]-= (block[block_index-(block_index%total_row)+iteration].nz[i*block[i].sb_row+k]*block[(block_index%total_row)+(iteration*total_row)].nz[k*block[i].sb_row+j]);
                        printf("new block %d value row %d col %d: %f \n", block_index, i, j, block[block_index].nz[(i*(block[i].sb_col))+j]);

                    }

                }
            }
    }
    iteration+=1;
}
s_index = nb-1;
    for(int k = 0; k< block[s_index].sb_row; k++){
        for(int i = k; i<block[s_index].sb_row; i++){
            upper[k][i] = block[s_index].nz[k*block[s_index].sb_row+i];
        }
        for(int i = k+1; i<block[s_index].sb_row; i++){
            printf("LOWER IK %f = BLOCK IK %f / UPPER KK %f", lower_block[i][k], (block[s_index].nz[i*block[s_index].sb_row+k]), upper[k][k]);
            lower[i][k] = (block[s_index].nz[i*block[s_index].sb_row+k])/upper[k][k];
             
        }
        for(int j = k+1; j<block[s_index].sb_row; j++){
            for(int i = k+1; i<block[s_index].sb_row; i++){
            block[s_index].nz[i*block[s_index].sb_row+j] -= (lower[i][k]*upper[k][j]);
            printf("new block %d value row %d col %d: %f \n", s_index, i, j, block[s_index].nz[(i*(block[i].sb_col))+j]);

            }
        }

    }
    for(int k = 0; k< block[s_index].sb_row; k++){
        for(int i = k; i<block[s_index].sb_row; i++){
            block[s_index].nz[k*block[s_index].sb_row+i] = upper[k][i];
        }
        for(int i = k+1; i<block[s_index].sb_row; i++){
            (block[s_index].nz[i*block[s_index].sb_row+k]) = lower[i][k];
             
        }
    }
printf("\nBlock values: \n");
for(int i = 0; i<nb; i++){
    printf("\nBlock %d \n", i);
    for(int k = 0; k<(n*n)/nb; k++){
        printf("%f ",block[i].nz[k]);
    }
    printf("\n");
}

// //block format
// // for(int k = 0; k<n; k++){ 
// //         for(int i = k; i < n; i++){
// //             upper[k][i] = mat[k][i];
// //         }
// //         for(int i = k+1; i<n; i++){
// //             lower[i][k] = mat[i][k]/upper[k][k];
// //         }
// //         for(int j = k+1; j<n; j++){
// //             for(int i = k+1; i<n; i++){
// //                 mat[i][j] -= (lower[i][k]*upper[k][j]);
// //             }
// //         }   
// //     }
// nb = 2;
// //Category 1
// for(int num = 0; num<(n/nb); num++){
// //   printf("\nCategory 1\n");
// //   for(int k = 0+(nb*num); k<nb+(nb*num); k++){
// //     for(int i = 0+(nb*num); i<nb+(nb*num); i++){
// //       upper[k][i] = mat[k][i];
// //     }
// //     for(int i = k+1; i<n; i++){
// //         lower[i][k] = mat[i][k]/upper[k][k];
// //     }
// //     for(int j = k+1; j<n; j++){
// //       for(int i = k+1; i<n; i++){
// //           mat[i][j] -= (lower[i][k]*upper[k][j]);
// //         }
// //     } 
// //   }
//   printf("\n\n MATRIX \n");
//     for (int i = 0+(nb*num); i < nb+(nb*num); i++) {
//         for (int j = 0+(nb*num); j < nb+(nb*num); j++) {
            
//             printf("%f ", mat[i][j]);
//         }
//         printf("\n");
//     }


//   // for(int block = 1; block < n/nb; block++){
//   //   for(int k = 1+(nb*num); k<nb+(num*nb); k++){
//   //     for(int i = nb+1+(nb*num); i<block*nb*2+(nb*num); i++){
//   //       mat[i][k]/=mat[k][k];
//   //     }
//   //     for(int j = k+1; j<nb; j++){
//   //       for(int i = block*nb+(nb*num); i<block*nb*2+(nb*num); i++){
//   //         mat[i][j]-= (mat[i][k]*mat[k][j]);
//   //       }
//   //     }
//   //   }
//   // }

//     printf("\nCategory 2 \n");
//   int block = 1;
//   for(int row = nb+(nb*num); row< n; row++){
//     if(row%nb == 0 && row!=nb) printf("\n");
//     for(int col = 0+(nb*num); col < nb+(nb*num); col++){
//       printf("%f ", mat[row][col]);
//     }
//     printf("\n");
//     // block++;
//     // if(block == nb+1){
//     //   block = 1;
//     //   printf("\n");
//     // }
    
//   }

// // for(int block = 1; block<n/nb; block++){
// //   for(int k = 1+(nb*num); k<nb+(nb*num); k++){
// //     for(int j =nb+1+(nb*num); j<block*nb*2+(nb*num); j++){
// //       for(int i = k+1; i<nb; i++){
// //         mat[i][j] -= (mat[i][k]*mat[k][j]);
// //       }
// //     }
// //   }
// // }
//   printf("\nCategory 3\n");
//   for(int col = nb+(nb*num); col < n; col+=nb){
//     for(int row = 0+(nb*num); row<nb+(nb*num); row++){
//       for(int col1 = col; col1<col+nb; col1++){
//         printf("%f ", mat[row][col1]);
//       }
//       printf("\n");
//     }
//     printf("\n");
//   }
// // for(int blockcol = 1; blockcol<((n-(nb*num))/nb); blockcol++){
// //   for(int blockrow = 1; blockrow<((n-(nb*num))/nb); blockrow++){
// //     for(int k = 1+(nb*num); k<nb+(nb*num); k++){
// //       for(int j = 1+(nb*blockcol)+(nb*num); j<nb*2*block+(nb*num); j++){
// //         for(int i = 1+(nb*blockrow)+(nb*num); i<nb*2*block+(nb*num); i++){
// //           mat[i][j] -= mat[i][k]*mat[k][j];
// //         }
// //       }
// //     }
// //   }
// // }
//   printf("\nCategory 4 \n");
//   for(int outer_row = nb+(nb*num); outer_row<n; outer_row+=nb){
//     for(int outer_col = nb+(nb*num); outer_col<n; outer_col+=nb){
//       for(int row = outer_row; row<outer_row+nb; row++){
//         for(int col = outer_col; col<outer_col+nb; col++){
//           printf("%f ", mat[row][col]);
//         }
//         printf("\n");
//       }
//       printf("\n\n");
//     }
//   }
// }








    // for(int block_num = 0; block_num < (n*n/nb*nb); block_num++){ //# block
    //   printf("Block %d \n", block_num);
    //   for(int blockr = block_num*nb; blockr<nb; blockr++){
    //     for(int blockc = block_num*nb; blockc<nb; blockc++){
    //       printf("%f ", mat[blockr][blockc]);
    //     }
    //     printf("\n");
    //   }
    //   printf("\n");
    // }
        // printf("\n");
        //   for(int l = 0; l<=0; l++){//which column
        //       int icsr = k+1; //current row index
        //       for (idx_t a=begin_pos[k]+1; a < begin_pos[k+1]; a++){
        //         printf("\n\n\n %d", csr[a]);
        //          printf("%lf\t", nzval[a]);
        //           if(csr[a]==icsr) {
        //             lower[icsr][k] = nzval[a]/upper[k][k]; 
        //             printf("%d", k);
        //           }
        //           else{
        //             lower[l][k] = 0/upper[k][k]; 
        //           }
        //           icsr++;
        //           // printf("%f ", lower[l][k]);
        //       }
        //     printf("\n ");
        //   }  
          // for(int j = k+1; j<n; j++){
          //   for (idx_t a=begin_pos[k+1]; a < begin_pos[k+2]; a++){
          //       for(int i = k+1; i<n; i++){
          //         if(csr[a] == j){
          //         zero = 0;
          //         nzval[a] -= (lower[i][k]*upper[k][j]);

          //       }
          //       if(zero == 1){
          //           nzval[a] = lower[j][k]*upper[k][]
          //       }
          //       zero = 1;
          //       }
              
          //   }
          // }  
      

    //   for(int k = 0; k<n; k++){ 
    //     for(int i = k; i < n; i++){
    //         upper[k][i] = mat[k][i];
    //     }
    //     for(int i = k+1; i<n; i++){
    //         lower[i][k] = mat[i][k]/upper[k][k];
    //     }
    //     for(int j = k+1; j<n; j++){
    //         for(int i = k+1; i<n; i++){
    //             mat[i][j] -= (lower[i][k]*upper[k][j]);
    //         }
    //     }   
    // }

      

    // printf("\n\n LOWER \n");
    //   for (int i = 0; i < n; i++) {
    //     for (int j = 0; j < n; j++) {
            
    //         printf("%f ", lower[i][j]);
    //     }
    //     printf("\n");
    // }
    //     printf("\n\n UPPER \n");
    //   for (int i = 0; i < n; i++) {
    //     for (int j = 0; j < n; j++) {
            
    //         printf("%f ", upper[i][j]);
    //     }
    //     printf("\n");
    // }

      
    //     for(int i = k+1; i<n; i++){
    //         lower[i][k] = mat[i][k]/upper[k][k];
    //     }
    //     for(int j = k+1; j<n; j++){
    //         for(int i = k+1; i<n; i++){
    //             mat[i][j] -= (lower[i][k]*upper[k][j]);
    //         }
    //     }   
    // }
	//reading each row

	free(begin_pos);
	free(csr);
	free(nzval);
  
}

	// dcreate_matrix_postfix(&A, nrhs, &b, &ldb, &xtrue, &ldx, fp, postfix, &grid);
