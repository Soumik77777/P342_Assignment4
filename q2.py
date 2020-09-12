'''
Q2. Check whether inverse of the following matrix exists. If yes, find inverse.
[ 0 2 8 6]
[ 0 0 1 2]
[ 0 1 0 1]
[ 3 7 1 0]
'''
#importing library
from library import *

#procedure
matrix_a = read_matrix('mat_2.txt')
n= len(matrix_a)

#Checking if inverse exists
(lower_mat, upper_mat) = lu_decomposition (matrix_a, n)

'''
[A]= [L]*[U]
∴det[A]= det[L]* det[U]
determinant of any triangular matrix= product of diagonals
∴det[L] =1 & det[A] = det[U]
'''
det_a = 1
for i in range(4):
    det_a = det_a*upper_mat[i][i]

if det_a == 0:
    print("The determinant of the matrix is zero. So, the inverse doesn't exist.")
else:
    print("The inverse of the matrix exists. It is:")

    #calculating and printing inverse
    inverse= inverse_by_lu_decomposition(matrix_a, n)

    for i in inverse:
        print(i)


#appending the exact result
'''
The inverse of the matrix exists. It is:
[0.33, -0.25, 1.67, -1.83]
[0.0, 0.08, -0.67, 0.83]
[0.0, 0.17, -0.33, -0.33]
[0.0, -0.08, 0.67, 0.17]
'''
