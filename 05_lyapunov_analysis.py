
import numpy as np
from scipy.integrate import solve_ivp

# Lorenz system parameters
sigma = 10
rho = 28
beta = 8 / 3

# Lorenz differential equations
def lorenz(t, state):

    x, y, z = state

    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z

    return [dx, dy, dz]


# Initial condition
initial = np.array([1.0, 1.0, 1.0])

# Initial perturbation
epsilon = 1e-8
perturbation = np.array([epsilon, 0, 0])

nearby = initial + perturbation

# Numerical settings
dt = 0.1
total_time = 100
steps = int(total_time / dt)

# Accumulate logarithmic separation
sum_log = 0.0
valid_time = 0.0

for i in range(steps):

    # Integrate both trajectories over one interval
    sol1 = solve_ivp(
        lorenz,
        (0, dt),
        initial,
        rtol=1e-10,
        atol=1e-12
    )

    sol2 = solve_ivp(
        lorenz,
        (0, dt),
        nearby,
        rtol=1e-10,
        atol=1e-12
    )

    initial = sol1.y[:, -1]
    nearby = sol2.y[:, -1]

    # Calculate separation
    difference = nearby - initial
    distance = np.linalg.norm(difference)

    if distance > 0 and np.isfinite(distance):

        sum_log += np.log(distance / epsilon)
        valid_time += dt

        # Renormalise perturbation
        nearby = initial + epsilon * difference / distance

    else:
        raise RuntimeError("Invalid trajectory separation")

# Estimate largest Lyapunov exponent
lyapunov_exponent = sum_log / valid_time

print("Lorenz System Lyapunov Analysis")
print("--------------------------------")
print("Sigma:", sigma)
print("Rho:", rho)
print("Beta:", beta)
print("Total simulation time:", valid_time)
print("Largest Lyapunov exponent:", lyapunov_exponent)