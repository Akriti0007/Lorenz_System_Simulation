import numpy as np

# Lorenz system parameters
sigma = 10
rho = 28
beta = 8 / 3

# Equilibrium point 1
E0 = np.array([0, 0, 0])

# Equilibrium points 2 and 3
value = np.sqrt(beta * (rho - 1))

E_plus = np.array([value, value, rho - 1])
E_minus = np.array([-value, -value, rho - 1])

# Display equilibrium points
print("Equilibrium Points")

print("E0 =", E0)
print("E+ =", E_plus)
print("E- =", E_minus)

# Define Jacobian matrix
def jacobian(x, y, z):

    return np.array([
        [-sigma, sigma, 0],
        [rho - z, -1, -x],
        [y, x, -beta]
    ])

# Calculate eigenvalues
points = {
    "E0": E0,
    "E+": E_plus,
    "E-": E_minus
}

print("\nStability Analysis")

for name, point in points.items():

    J = jacobian(*point)

    eigenvalues = np.linalg.eigvals(J)

    print("\n", name)
    print("Jacobian Matrix:")
    print(J)

    print("Eigenvalues:")
    print(eigenvalues)