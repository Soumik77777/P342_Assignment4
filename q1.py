'''
Q1. Solve for xi

x1 + x3 + 2x4 = 6
x2 - 2x3 = -3
x1 + 2x2 - x3 = -2
2x1 + x2 + 3x3 - 2x4 = 0
'''
#importing library
from library import *

#procedure
matrix_a = read_matrix('mat_1.txt')
vector_b = [6, -3, -2, 0]
n= len(matrix_a)

(lower_mat, upper_mat) = lu_decomposition (matrix_a, n)

x = forward_backward_substitution (lower_mat, upper_mat, vector_b, n)


#printing solution for each xi
print("The solution of the set of linear equation is:")
print("x1=", x[0])
print("x2=", x[1])
print("x3=", x[2])
print("x4=", x[3])


#appending the exact result
'''
The solution of the set of linear equation is:
x1= 1.0
x2= -1.0
x3= 1.0
x4= 2.0
'''
