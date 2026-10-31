
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Lorenz system parameters
sigma = 10
rho = 28
beta = 8 / 3

# Initial condition
initial = np.array([1.0, 1.0, 1.0])

# Numerical settings
epsilon = 1e-6
dt = 0.01
burn_time = 50
total_time = 500

# Lorenz differential equations
def lorenz(t, state):

    x, y, z = state

    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z

    return [dx, dy, dz]


# STEP 1: Remove initial transient
burn = solve_ivp(
    lorenz,
    (0, burn_time),
    initial,
    rtol=1e-10,
    atol=1e-12
)

initial = burn.y[:, -1]

# Introduce small perturbation
nearby = initial + np.array([epsilon, 0, 0])


# STEP 2: Calculate Lyapunov exponent
steps = int(total_time / dt)

sum_log = 0.0
elapsed = 0.0

times = []
estimates = []

for i in range(steps):

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

    # Calculate distance between trajectories
    difference = nearby - initial
    distance = np.linalg.norm(difference)

    if distance <= 0 or not np.isfinite(distance):
        raise RuntimeError("Invalid trajectory separation")

    # Accumulate logarithmic growth
    sum_log += np.log(distance / epsilon)

    elapsed += dt

    # Renormalise perturbation
    nearby = initial + epsilon * difference / distance

    # Store running estimate
    times.append(elapsed)
    estimates.append(sum_log / elapsed)


# STEP 3: Display final result
lyapunov_exponent = sum_log / elapsed

print("Lorenz System Lyapunov Analysis")
print("--------------------------------")
print("Sigma:", sigma)
print("Rho:", rho)
print("Beta:", beta)
print("Transient time:", burn_time)
print("Measurement time:", elapsed)
print("Largest Lyapunov exponent:", lyapunov_exponent)


# STEP 4: Plot convergence graph
plt.figure(figsize=(10, 6))

plt.plot(
    times,
    estimates,
    color="blue",
    label="Estimated Lyapunov exponent"
)

plt.axhline(
    0.9,
    color="red",
    linestyle="--",
    label="Approximate reference: 0.9"
)

plt.xlabel("Measurement Time")
plt.ylabel("Estimated Largest Lyapunov Exponent")

plt.title("Lorenz System: Lyapunov Exponent Convergence")

plt.legend()
plt.grid(True)

# Display the convergence after initial fluctuations
plt.xlim(0, 500)
plt.ylim(0.5, 1.3)
plt.tight_layout()

plt.savefig("lyapunov_convergence.png", dpi=300)

plt.show()