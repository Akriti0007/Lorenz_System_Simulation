import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Lorenz system parameters
sigma = 10
rho = 28
beta = 8 / 3

# Define the system of differential equations
def lorenz(t, state):
    x, y, z = state

    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z

    return [dx, dy, dz]


# Initial conditions
initial_conditions = [1, 1, 1]

# Simulation time
t = np.linspace(0, 50, 5001)

# Solve the differential equations
solution = solve_ivp(
    lorenz,
    (0, 50),
    initial_conditions,
    t_eval=t
)

# Extract results
x = solution.y[0]
y = solution.y[1]
z = solution.y[2]

# Plot the Lorenz attractor
fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111, projection="3d")

ax.plot(x, y, z, linewidth=0.5)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title("Lorenz Attractor")

plt.tight_layout()
plt.savefig("lorenz_attractor.png", dpi=300)
plt.show()

# Experiment 2: Sensitivity to Initial Conditions

initial_1 = [1, 1, 1]
initial_2 = [1.0001, 1, 1]

time = np.linspace(0, 40, 4001)

# First simulation
sol1 = solve_ivp(
    lorenz, (0, 40), initial_1,
    t_eval=time, rtol=1e-9, atol=1e-12
)

# Second simulation
sol2 = solve_ivp(
    lorenz, (0, 40), initial_2,
    t_eval=time, rtol=1e-9, atol=1e-12
)

# Calculate trajectory separation
distance = np.linalg.norm(
    sol1.y - sol2.y, axis=0
)

# Plot trajectory separation
plt.figure(figsize=(9, 5))

plt.semilogy(time, distance)

plt.xlabel("Time")
plt.ylabel("Distance between trajectories")
plt.title("Sensitivity to Initial Conditions")
plt.grid(True)

plt.savefig("sensitivity_analysis.png", dpi=300)
plt.show()

# Experiment 3: Effect of System Parameters

rho_values = [10, 28, 40]

fig = plt.figure(figsize=(16, 5))

for i, rho_value in enumerate(rho_values):

    # Define Lorenz system for each rho value
    def lorenz_parameter(t, state):

        x, y, z = state

        dx = sigma * (y - x)
        dy = x * (rho_value - z) - y
        dz = x * y - beta * z

        return [dx, dy, dz]

    # Solve differential equations
    sol = solve_ivp(
        lorenz_parameter,
        (0, 50),
        [1, 1, 1],
        t_eval=np.linspace(0, 50, 5001),
        rtol=1e-9,
        atol=1e-12
    )

    # Plot 3D attractor
    ax = fig.add_subplot(1, 3, i + 1, projection="3d")

    ax.plot(
        sol.y[0],
        sol.y[1],
        sol.y[2],
        linewidth=0.5
    )

    ax.set_title(f"Rho = {rho_value}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

fig.suptitle("Effect of Rho on Lorenz System")

plt.tight_layout()
plt.savefig("parameter_variation.png", dpi=300)
plt.show()

