import numpy as np

sigma = 10
rho = 28
beta = 8 / 3

value = np.sqrt(beta * (rho - 1))

equilibrium_points = {
    "Origin": (0, 0, 0),
    "E+": (value, value, rho - 1),
    "E-": (-value, -value, rho - 1)
}

def jacobian(x, y, z):

    return np.array([
        [-sigma, sigma, 0],
        [rho - z, -1, -x],
        [y, x, -beta]
    ])

for name, point in equilibrium_points.items():

    J = jacobian(*point)

    eigenvalues = np.linalg.eigvals(J)

    print("\nEquilibrium:", name)
    print("Eigenvalues:", eigenvalues)

    if np.all(np.real(eigenvalues) < 0):
        print("Stable equilibrium")

    elif np.any(np.real(eigenvalues) > 0):
        print("Unstable equilibrium")

    else:
        print("Further analysis required")