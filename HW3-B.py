import math


def gaussian_elimination(coefficients, constants):
    """
    Solves a system of linear equations using Gaussian Elimination method.

    Parameters:
        coefficients (list of lists): The coefficients of the system of equations.
        constants (list): The constants of the system of equations.

    Returns:
        list: The solution to the system of equations.

    Raises:
        ValueError: If a zero pivot is encountered during the elimination process.
    """
    n = len(coefficients)

    # Forward Elimination
    for i in range(n):
        # Find the row with the maximum absolute coefficient for pivot selection
        max_index = i
        for j in range(i + 1, n):
            if abs(coefficients[j][i]) > abs(coefficients[max_index][i]):
                max_index = j

        # Swap rows to make the pivot element the largest in the column
        coefficients[i], coefficients[max_index] = coefficients[max_index], coefficients[i]
        constants[i], constants[max_index] = constants[max_index], constants[i]

        pivot = coefficients[i][i]
        if pivot == 0:
            raise ValueError("Zero pivot encountered")

        # Perform elimination below the pivot element
        for j in range(i + 1, n):
            factor = coefficients[j][i] / pivot
            for k in range(i, n):
                coefficients[j][k] -= factor * coefficients[i][k]
            constants[j] -= factor * constants[i]

    # Back Substitution
    solution = [0] * n
    for i in range(n - 1, -1, -1):
        # Calculate the solution for each variable starting from the bottom
        solution[i] = constants[i] / coefficients[i][i]
        for j in range(i):
            # Update constants to be used in previous equations
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
    print(f"x{i + 1} = {val}")
