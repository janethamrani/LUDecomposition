# LU Decomposition

Factors a matrix as a product of a lower and upper triangular matrix

<img width="500" alt="image" src="https://user-images.githubusercontent.com/26263012/201554361-91bf8b94-2cc9-42d6-aefb-0722e1668298.png">

### Right Looking Algorithm

This algorithm works by computing the partial factorization of a matrix by evaluating each row and column in different iterations in a diagonal fashion and producing an updated submatrix of A consisting of the porition of the original matrix that has not yet been factored. In the following iteration, the proceeding row and column of the submatrix is factored.

<img width="300" alt="image" src="https://user-images.githubusercontent.com/26263012/201565761-be83762a-b447-455e-8508-76f7aa4c6919.png">

The code first creates an empty lower and upper matrix, lower being an identity matrix. Next, a loop that starts with index <i>k = 0</i> to the matrix length is declared;<i> k</i> represents the blocksize of the active submatrix. 
In this loop are three nested loops. The first loop starts with index<i> i = k</i> through the matrix length and sets the upper matrix row equal to its corresponding element in <i>A</i>. It creates the upper triangular matrix as <i>i</i> is initally set to <i>k</i>. In the first outer loop iteration, when<i> k = 0</i>, it will set the entire first row of the upper matrix. 
In the next loop, the index <i>i</i> is set to <i>k+1</i> through the matrix length (<i>k+1</i> to avoid overlapping). This loop sets the block column of the lower matrix to its corresponding value in the <i>A</i> matrix divided by the pivot element (so that the product of the upper and lower matrices can equal <i>A</i>), which is the top left element in the active submatrix.
The last loop is a nested for-loop that updates submatrix for the upcoming iteration using the current block row and column

```sh
for(int k = 0; k<n; k++){ 
        for(int i = k; i < n; i++){
            upper[k][i] = mat[k][i];
        }
        for(int i = k+1; i<n; i++){
            lower[i][k] = mat[i][k]/upper[k][k];
        }
        for(int j = k+1; j<n; j++){
            for(int i = k+1; i<n; i++){
                mat[i][j] -= (lower[i][k]*upper[k][j]);
            }
        }   
    }
```

## CSR - Compressed Sparse Row Format

![image](https://user-images.githubusercontent.com/26263012/203053484-ab0e3165-5480-4360-a494-500bfb734157.png)

<img width="350" alt="image" src="https://user-images.githubusercontent.com/26263012/203054148-e599fc6e-654f-49a8-98b9-0ee0942f05b1.png">

```sh
for (idx_t i=0; i< n; i++)
   {
     printf("\nNeighbor of vertex %d  begin_pos[%d]:%d  begin_pos[%d]:%d\t",i,i,begin_pos[i],i+1,begin_pos[i+1]);
     for (idx_t j=begin_pos[i]; j< begin_pos[i+1]; j++)
     {
       printf("%d %lf\t",csr[j],nzval[j]);
     }
   }
 ```
 
<img width="600" alt="image" src="https://user-images.githubusercontent.com/26263012/203054070-716e7f96-b3d2-44d4-945a-08dcd5984573.png">


Print matrix from CSC Format
```sh
int col = 0;
_Bool zero = 1;
      for (int k = 0; k<3; k++){ //getting elements from row i
          for (idx_t i=0; i< n; i++){ 
            for (idx_t j=begin_pos[i]; j< begin_pos[i+1]; j++)
            {
              if (csr[j] == k){
                    zero = 0;
                    upper[k][col] = nzval[j];
                    printf("%d ", upper[k][col]);
                    col++; if (col == 2) col = 0;
                }
              }
              if (zero == 1) printf("%d ", 0); 
              zero = 1;
            }
          }
```

## Block Format
<img width="100%" alt="image" src="https://user-images.githubusercontent.com/26263012/203840336-27a3a69f-4df4-451d-b2cd-b6d551c996ff.png">

Print blocks from each category
```sh
  printf("\nCategory 1\n");
  for(int i = 0; i<nb; i++){
    for(int k = 0; k<nb; k++){
      printf("%f ", mat[i][k]);
    }
    printf("\n");
  }
  printf("\nCategory 2 \n");
  int block = 1;
  for(int row = nb; row< n; row++){
    for(int col = 0; col < nb; col++){
      printf("%f ", mat[row][col]);
    }
    printf("\n");
    block++;
    if(block == nb+1){
      block = 1;
      printf("\n");
    }
  }

  printf("\nCategory 3\n");
  for(int col = nb; col < n; col+=nb){
    for(int row = 0; row<nb; row++){
      for(int col1 = col; col1<col+nb; col1++){
        printf("%f ", mat[row][col1]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\nCategory 4 \n");
  for(int outer_row = nb; outer_row<n; outer_row+=nb){
    for(int outer_col = nb; outer_col<n; outer_col+=nb){
      for(int row = outer_row; row<outer_row+nb; row++){
        for(int col = outer_col; col<outer_col+nb; col++){
          printf("%f ", mat[row][col]);
        }
        printf("\n");
      }
      printf("\n\n");
    }
  }
```
## Block LU Decomposition Implementation

### Pseudocode
<img width="100%" alt="image" src="https://user-images.githubusercontent.com/26263012/204322770-e2cfe6ea-d7bd-47c6-b077-44b2b89e5640.png">

### My Implementation
```sh 
for(int num = 0; num<(n/nb); num++){
  //Category 1
  for(int k = 0+(nb*num); k<nb+(nb*num); k++){
    for(int i = 0+(nb*num); i<nb+(nb*num); i++){
      upper[k][i] = mat[k][i];
    }
    for(int i = k+1; i<n; i++){
        lower[i][k] = mat[i][k]/upper[k][k];
    }
    for(int j = k+1; j<n; j++){
      for(int i = k+1; i<n; i++){
          mat[i][j] -= (lower[i][k]*upper[k][j]);
        }
    } 
  //Category 2
  for(int block = 1; block < n/nb; block++){
    for(int k = 1+(nb*num); k<nb+(num*nb); k++){
      for(int i = nb+1+(nb*num); i<block*nb*2+(nb*num); i++){
        mat[i][k]/=mat[k][k];
      }
      for(int j = k+1; j<nb; j++){
        for(int i = block*nb+(nb*num); i<block*nb*2+(nb*num); i++){
          mat[i][j]-= (mat[i][k]*mat[k][j]);
        }
      }
    }
  }
//Category 3
 for(int block = 1; block<n/nb; block++){
    for(int k = 1+(nb*num); k<nb+(nb*num); k++){
       for(int j =nb+1+(nb*num); j<block*nb*2+(nb*num); j++){
          for(int i = k+1; i<nb; i++){
             mat[i][j] -= (mat[i][k]*mat[k][j]);
          }
       }
    }
  }
//Category 4  
  for(int blockcol = 1; blockcol<((n-(nb*num))/nb); blockcol++){
     for(int blockrow = 1; blockrow<((n-(nb*num))/nb); blockrow++){
       for(int k = 1+(nb*num); k<nb+(nb*num); k++){
         for(int j = 1+(nb*blockcol)+(nb*num); j<nb*2*block+(nb*num); j++){
           for(int i = 1+(nb*blockrow)+(nb*num); i<nb*2*block+(nb*num); i++){
             mat[i][j] -= mat[i][k]*mat[k][j];
           }
         }
       }
     }
   }
}
```

### Struct Implementation
<img width="50%" alt="image" src="https://user-images.githubusercontent.com/26263012/207111020-b6986dab-5fba-471f-8fb0-c7785fb068c1.png">
```sh 
int total_row = (n*n/nb);
int s_index = 0; 
for(int s_index = 0; s_index<nb-1; s_index+=total_row+1){
    int iteration = 0;
    //Category 1 -- STRUCT
    for(int k = 0; k< block[s_index].sb_row; k++){
        for(int i = k; i<block[s_index].sb_row; i++){
            upper[k][i] = block[s_index].nz[k*block[s_index].sb_row+i];
        }
        for(int i = k+1; i<block[s_index].sb_row; i++){
            lower[i][k] = (block[s_index].nz[i*block[s_index].sb_row+k])/upper[k][k];
        }
        for(int j = k+1; j<block[s_index].sb_row; j++){
            for(int i = k+1; i<block[s_index].sb_row; i++){
            block[s_index].nz[i*block[s_index].sb_row+j] -= (lower[i][k]*upper[k][j]);
            }
        }
    }
    //Category 2 -- STRUCT
    printf("\n\n");
    for(int block_index = s_index+total_row; block_index<nb; block_index+=(total_row)){
        for(int k = 0; k<(block[block_index].sb_row*block[block_index].sb_col); k++){
            for(int i = 0; i<(block[block_index].sb_row*block[block_index].sb_col); i++){
                block[block_index].nz[i*block[block_inde    x].sb_col+k]/=block[s_index].nz[k*block[block_index].sb_col+k]; mat[i][k]/=mat[k][k];
            }
            for(int j = k+1; j<(block[block_index].sb_row*block[block_index].sb_col);j++){
                for(int i = k+1; i <(block[block_index].sb_row*block[block_index].sb_col);i++){
                    block[block_index].nz[i*block[block_index].sb_col+j]-= ((block[block_index].nz[i*block[block_index].sb_col+k])*block[s_index].nz[k*block[block_index].sb_col+j]); //mat[i][j]-= (mat[i][k]*mat[k][j]);
                }
            } 
            
        }
        printf("\n\n");
    }
    //Category 3 -- STRUCT
    for(int block_index = s_index+1; block_index<(total_row*(iteration+1); block_index++){
        for(int k = 1; k<block[block_index].sb_row; k++){
            for(int j = 1; j<block[block_index].sb_row; j++){
                for(int i = k+1;j<block[block_index].sb_row; j++){
                    block[block_index].nz[i*block[block_index].sb_col+j]= ((block[block_index].nz[k*block[block_index].sb_col+j])*block[s_index].nz[i*block[block_index].sb_col+k]);
                }
            }
        }
    }

    //Category 4 -- STRUCT
    for(int block_index = s_index+total_row+1; block_index<nb; block_index++){
        if((block_index%total_row) == 0+iteration) block_index++;
            for(int k = 1; k< block[block_index].sb_col; k++){
                for(int j = 1; j< block[block_index].sb_col; j++){
                    for(int i = 1; i< block[block_index].sb_col; i++){
                        block[block_index].nz[i*block[block_index].sb_col+j]-=block[block_index-(block_index%total_row)+iteration].nz[i*block[i].sb_row+j]*block[(block_index%total_row)+(iteration*total_row)].nz[k*block[i].sb_row+j];
                    }

                }
            }
    }
    iteration+=1;
}
```
