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
      for(int k = 0; k<3; k++){ //getting elements from row i
          for (idx_t i=0; i< n; i++){ 
            for (idx_t j=begin_pos[i]; j< begin_pos[i+1]; j++)
            {
              if(csr[j] == k){
                    zero = 0;
                    upper[k][col] = nzval[j];
                    printf("%d ", upper[k][col]);
                    col++; if (col == 2) col = 0;
                }
              }
              if(zero == 1) printf("%d ", 0); 
              zero = 1;
            }
          }
```
