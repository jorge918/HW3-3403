# Gaussian elimination function
def gaussian_elimination(coefficients, constants):
    n = len(coefficients)
    # Forward Elimination
    for i in range(n):
        max_index = i
        for j in range(i + 1, n):
            if abs(coefficients[j][i]) > abs(coefficients[max_index][i]):
                max_index = j
        coefficients[i], coefficients[max_index] = coefficients[max_index], coefficients[i]
        constants[i], constants[max_index] = constants[max_index], constants[i]
        pivot = coefficients[i][i]
        if pivot == 0:
            raise ValueError("Zero pivot encountered")
        for j in range(i + 1, n):
            factor = coefficients[j][i] / pivot
            for k in range(i, n):
                coefficients[j][k] -= factor * coefficients[i][k]
            constants[j] -= factor * constants[i]
    # Back Substitution
    solution = [0] * n
    for i in range(n - 1, -1, -1):
        solution[i] = constants[i] / coefficients[i][i]
        for j in range(i):
            constants[j] -= coefficients[j][i] * solution[i]
    return solution

# Coefficients of the system of equations
coefficients = [
    [1, -1, 3, 2],
    [-1, 5, -5, -2],
    [3, -5, 19, 3],
    [2, 2, 3, 21]
]

# Constants of the system of equations
constants = [15, -35, 94, 1]

# Solve the system of equations
solution = gaussian_elimination(coefficients, constants)

# Print the solution
print("Solution:")
for i, val in enumerate(solution):
    print(f"x{i+1} = {val}")