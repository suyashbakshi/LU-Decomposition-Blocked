# LU-Decomposition-Blocked
Implementation of LU decomposition's blocked version of a matrix. 
Feel free to use/modify/distribute the code.
A star wouldn't go amiss if this implementation proved helpful/saved some of your valuable time. :)

The code structure of "LU_Decomposition()" is taken from NERSC website example here: https://www.nersc.gov/users/software/programming-models/openmp/openmp-tasking/openmp-tasking-example-lu/

Implementation of row_func(), col_func(), inner_func() and diag_func() was implemented by me (Suyash Bakshi).

References from NERSC website:

"""
```
Starting with the serial version of LU, shown below, we can see that there is a division with each diagonal element on every element in the rest of the row, and then the changes in this column are then propagated to each row after.
```
```
for(int i=0; i<size; i++) {
    for(int j=i+1; j<size; j++) {
        A[j*size+i] /= A[i*size+i];
        for(int k=i+1; k<size; k++) {
            A[j*size+k] -= A[j*size+i] * A[i*size+k];
        }
    }
}
```
"""

The code LU_decomp_blocked.c implements the following idea:
```
"""
```
This can be broken into 4 functions and applied to a blocked version of the matrix. The first is the same as the original algorithm above, and it's applied to the diagonal blocks. Then the changes from the diagonal block are propegated to the blocks in the same column and the blocks in the same row. A fourth function propegates these changes to the rest of the remaining matrix, and then this is repeated on that remaining inner matrix. Breaking the matrix into blocks not only enables parallelism, but also improves cache usage.
```
for(int i=0; i<num_blocks; i++) {
    diag_func( block_list[i][i] );
    for(int j=i+1; j<num_blocks; j++) {
        row_func( block_list[i][j], block_list[i][i] );
    }
    for(int j=i+1; j<num_blocks; j++) {
        col_func( block_list[j][i], block_list[i][i] );
        for(int k=i+1; k<num_blocks; k++) {
            inner_func( block_list[j][k], block_list[i][k], 
                                          block_list[j][i] );
        }
    }
}
```
"""
