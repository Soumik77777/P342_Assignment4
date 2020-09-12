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
    return (m,v)


def lu_decomposition (matrix, n):

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

    identity = [[0 for i in range(4)] for j in range(4)]
    for i in range(4):
        identity[i][i] = 1
    x = []

    '''The following process could be done easily using a for loop.
    But for those matrices, which are changed by partial pivoting, the columns of final inverse are interchanged.
    But pivoting is important, because at on estep in decomposition, it is divided by diagonal element
    So it is done manually for each row, and result of substitution is appended at the end of each step.'''

    matrix_0 = matrix.copy()
    partial_pivot(matrix_0, identity[0], n)
    (lower_0, upper_0) = lu_decomposition(matrix_0, n)
    x0 = forward_backward_substitution(lower_0, upper_0, identity[0], n)
    x.append(x0)

    matrix_1 = matrix.copy()
    partial_pivot(matrix_1, identity[1], n)
    (lower_1, upper_1) = lu_decomposition(matrix_1, n)
    x1 = forward_backward_substitution(lower_1, upper_1, identity[1], n)
    x.append(x1)

    matrix_2 = matrix.copy()
    partial_pivot(matrix_2, identity[2], n)
    (lower_2, upper_2) = lu_decomposition(matrix_2, n)
    x2 = forward_backward_substitution(lower_2, upper_2, identity[2], n)
    x.append(x2)

    matrix_3 = matrix.copy()
    partial_pivot(matrix_3, identity[3], n)
    (lower_3, upper_3) = lu_decomposition(matrix_3, n)
    x3 = forward_backward_substitution(lower_3, upper_3, identity[3], n)
    x.append(x3)

    inverse = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            inverse[i][j] = round(x[j][i], 3)

    return (inverse)


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


