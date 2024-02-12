import math


def is_symmetric(matrix):
    """
    Checks if a matrix is symmetric.

    Parameters:
        matrix (list of lists): The matrix to check.

    Returns:
        bool: True if the matrix is symmetric, False otherwise.
    """
    n = len(matrix)
    for i in range(n):
        for j in range(i + 1, n):  # Loop only over upper triangular part
            if matrix[i][j] != matrix[j][i]:
                return False
    return True


def cholesky_decomposition(matrix):
    """
    Performs Cholesky decomposition on a symmetric positive-definite matrix.

    Parameters:
        matrix (list of lists): The matrix to decompose.

    Returns:
        list of lists: The lower triangular Cholesky factor.
    """
    n = len(matrix)
    lower = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1):
            if i == j:
                sum_val = sum(lower[i][k] ** 2 for k in range(j))
                lower[i][j] = (matrix[i][i] - sum_val) ** 0.5
            else:
                sum_val = sum(lower[i][k] * lower[j][k] for k in range(j))
                lower[i][j] = (matrix[i][j] - sum_val) / lower[j][j]

    return lower


def doolittle_method(matrix, vector):
    """
    Solves a system of linear equations using the Doolittle method.

    Parameters:
        matrix (list of lists): The coefficient matrix of the system.
        vector (list): The constants vector of the system.

    Returns:
        list: The solution vector.
    """
    n = len(matrix)
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]

    for i in range(n):
        L[i][i] = 1.0
        for j in range(i, n):
            sum_val = sum(L[i][k] * U[k][j] for k in range(i))
            U[i][j] = matrix[i][j] - sum_val

        for j in range(i + 1, n):
            sum_val = sum(L[j][k] * U[k][i] for k in range(i))
            L[j][i] = (matrix[j][i] - sum_val) / U[i][i]

    # Forward substitution to solve Ly = b
    y = [0.0] * n
    for i in range(n):
        y[i] = vector[i]
        for j in range(i):
            y[i] -= L[i][j] * y[j]

    # Back substitution to solve Ux = y
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = y[i]
        for j in range(i + 1, n):
            x[i] -= U[i][j] * x[j]
        x[i] /= U[i][i]

    return x  # Return the solution


# Define the matrix equations for the given problems
problem1_matrix = [
    [1, -2, 3, 2],
    [-1, 5, -5, -2],
    [3, -5, 19, 3],
    [2, -2, 3, 21]
]
problem1_vector = [15, -35, 94, 1]

problem2_matrix = [
    [1, 4, 0, 0],
    [2, 2, 3, 0],
    [4, 3, 6, 3],
    [0, 2, 3, 9]
]
problem2_vector = [20, 36, 122, 60]

# Check if matrices are symmetric
problem1_is_symmetric = is_symmetric(problem1_matrix)
problem2_is_symmetric = is_symmetric(problem2_matrix)

# Solve the matrix equations
if problem1_is_symmetric:
    cholesky_solution = cholesky_decomposition(problem1_matrix)
    print("Solution using Cholesky decomposition:")
    print("X1 =", cholesky_solution[0])
    print("X2 =", cholesky_solution[1])
    print("X3 =", cholesky_solution[2])
    print("X4 =", cholesky_solution[3])
else:
    doolittle_solution = doolittle_method(problem1_matrix, problem1_vector)
    print("Solution using Doolittle method:")
    print("X1 =", doolittle_solution[0])
    print("X2 =", doolittle_solution[1])
    print("X3 =", doolittle_solution[2])
    print("X4 =", doolittle_solution[3])

if problem2_is_symmetric:
    cholesky_solution = cholesky_decomposition(problem2_matrix)
    print("Solution using Cholesky decomposition:")
    print("X1 =", cholesky_solution[0])
    print("X2 =", cholesky_solution[1])
    print("X3 =", cholesky_solution[2])
    print("X4 =", cholesky_solution[3])
else:
    doolittle_solution = doolittle_method(problem2_matrix, problem2_vector)
    print("Solution using Doolittle method:")
    print("X1 =", doolittle_solution[0])
    print("X2 =", doolittle_solution[1])
    print("X3 =", doolittle_solution[2])
    print("X4 =", doolittle_solution[3])
