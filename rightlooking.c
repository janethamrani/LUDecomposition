#include<stdio.h>
void main()
{
    int n = 3;
int lower[n][n], upper[n][n];
int mat[4][4]
        = { { 2, 1, 2 }, { 4, 6, 3 }, { 4, 2, 8 } };
        
for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            // Checking if row is equal to column
            if (row == col)
                lower[row][col] = 1;
            else
                lower[row][col] = 0;
        }
    }
  for (int row = 0; row < n; row++)
        {
            for (int col = 0; col < n; col++)
            {
                upper[row][col] = 0;
            }
        }
 
    
   for(int k = 0; k<n; k++){ // U(k,k:n) = A(k,k:n)
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
     for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            
            printf("%d ", lower[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            
            printf("%d ", upper[i][j]);
        }
        printf("\n");
    }

    
}
