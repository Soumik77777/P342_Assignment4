''' 
Personal library for useful functions; Soumik Bhattacharyya, Roll-1811155

This includes all necessary funcions used in P342,
but the ones that are required in the current week, are placed at the top.
'''

def read_matrix (txt):
    with open(txt, 'r') as a:
        matrix = [[int(num) for num in row.split(' ')] for row in a]

    return matrix


def partial_pivot (m, v, n):
    for i in range (n-1):
        if m[i][i] ==0:
            for j in range (i+1,n):
                if abs(m[j][i]) > abs(m[i][i]):
                    m[i], m[j] = m[j], m[i]
                    v[i], v[j] = v[j], v[i]


def lu_decomposition (matrix, n):

# partial pivoting, because later we've to do a division by diagonal entry
    v=[0, 0, 0, 0]
    partial_pivot(matrix, v, n)

    upper_mat = [[0 for i in range(n)] for j in range(n)]
    lower_mat = [[0 for i in range(n)] for j in range(n)]

    for i in range(n):

        for j in range(i, n): #calculating upper matrix
            sum = 0
            for k in range(i):
                sum += (lower_mat[i][k] * upper_mat[k][j])
            upper_mat[i][j] = matrix[i][j] - sum

        for j in range(i, n): #calculating lower matrix
            if (i == j):
                lower_mat[i][i] = 1
            else:
                sum = 0
                for k in range(i):
                    sum += (lower_mat[j][k] * upper_mat[k][i])

                lower_mat[j][i] = ((matrix[j][i] - sum) / upper_mat[i][i])

    return (lower_mat, upper_mat)


def forward_backward_substitution (lower_mat, upper_mat, vector, n):
    '''
    If we have LUx=B,
    first we solve Ly=B, then Ux=y
    '''
    # forward-substitution
    y = [0] * n
    for i in range(n):
        sum = 0
        for j in range(i):
            sum += lower_mat[i][j] * y[j]

        y[i] = vector[i] - sum

    #backward-substitution
    x = [0] * n
    for i in reversed(range(n)):
        sum = 0
        for j in range(i + 1, n):
            sum+= upper_mat[i][j] * x[j]
        x[i] = (y[i] - sum)/ upper_mat[i][i]

    return (x)


def inverse_by_lu_decomposition (matrix, n):

# partial pivoting, because later we've to do a division by diagonal entry
    v = [0, 0, 0, 0]
    partial_pivot(matrix, v, n)

    inverse = [[0 for i in range(n)] for j in range(n)]

    (lower_mat, upper_mat) = lu_decomposition(matrix, n)

    for i in range (n):
        vector= [0 for j in range(n)]
        vector[i] = 1
        # we're going to calculate each column (j) of the inverse with corresponding vector

        x = forward_backward_substitution(lower_mat, upper_mat, vector, n)
        # x is the vector for each column j in inverse matrix

        for j in range(n):
            inverse [j][i] = round(x[j], 2) #rounding off the result

    return inverse


def matrix_multiplication (matrix_m, matrix_n, n):
    matrix_L = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                matrix_L[i][j] += matrix_n[k][j] * matrix_m[i][k]
    return matrix_L


def gauss_jordan (m, v, n):
    partial_pivot(m, v, n)
    for i in range (n):
        factor1 = m[i][i]
        for j in range (n):
            m[i][j] = m[i][j]/factor1
        for q in range (n):
            v[i][q] = v[i][q]/factor1

        for k in range (n):
            if k!=i and m[k][i]!=0:
                factor2 = m[k][i]
                for l in range (i,n):
                    m[k][l] = m[k][l] - factor2*m[i][l]
                for r in range (n):
                    v[k][r] = v[k][r] - factor2* v[i][r]


