import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Lorenz parameters
sigma = 10
rho = 28
beta = 8 / 3

# Lorenz equations
def lorenz(t, state):

    x, y, z = state

    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z

    return [dx, dy, dz]


# Calculate equilibrium points
value = np.sqrt(beta * (rho - 1))

E0 = [0, 0, 0]
Eplus = [value, value, rho - 1]
Eminus = [-value, -value, rho - 1]

# Solve the system
t = np.linspace(0, 50, 5001)

solution = solve_ivp(
    lorenz,
    (0, 50),
    [1, 1, 1],
    t_eval=t,
    rtol=1e-9,
    atol=1e-12
)

x, y, z = solution.y

# Plot Lorenz attractor
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

ax.plot(x, y, z, linewidth=0.5,
        label="Lorenz trajectory")

# Mark equilibrium points
ax.scatter(*E0, s=80, marker="o",
           label="E0")

ax.scatter(*Eplus, s=100, marker="^",
           label="E+")

ax.scatter(*Eminus, s=100, marker="^",
           label="E-")

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

ax.set_title("Lorenz Attractor and Equilibrium Points")
ax.legend()

plt.tight_layout()
plt.savefig("equilibrium_visualization.png", dpi=300)
plt.show()