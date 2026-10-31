
import numpy as np
import matplotlib.pyplot as plt
import csv

# ---------------------------------
# Lorenz System Parameters
# ---------------------------------

sigma = 10
beta = 8 / 3

rho_values = [10, 20, 23, 28, 35, 40]

# Numerical settings
epsilon = 1e-7
dt = 0.01

burn_time = 100
measurement_time = 200


# ---------------------------------
# Lorenz Differential Equations
# ---------------------------------

def lorenz(state, rho):

    x, y, z = state

    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z

    return np.array([dx, dy, dz])


# ---------------------------------
# Fourth Order Runge-Kutta Method
# ---------------------------------

def rk4_step(state, dt, rho):

    k1 = lorenz(state, rho)

    k2 = lorenz(
        state + 0.5 * dt * k1, rho
    )

    k3 = lorenz(
        state + 0.5 * dt * k2, rho
    )

    k4 = lorenz(
        state + dt * k3, rho
    )

    return state + (dt / 6) * (
        k1 + 2*k2 + 2*k3 + k4
    )


# ---------------------------------
# Largest Lyapunov Exponent
# ---------------------------------

def calculate_lyapunov(rho):

    initial = np.array([1.0, 1.0, 1.0])

    # Remove initial transient
    burn_steps = int(burn_time / dt)

    for _ in range(burn_steps):

        initial = rk4_step(initial, dt, rho)

    # Introduce small perturbation
    nearby = initial + np.array(
        [epsilon, 0, 0]
    )

    steps = int(measurement_time / dt)

    sum_log = 0.0

    estimates = []
    times = []

    for i in range(steps):

        initial = rk4_step(initial, dt, rho)

        nearby = rk4_step(nearby, dt, rho)

        difference = nearby - initial

        distance = np.linalg.norm(difference)

        if distance == 0 or not np.isfinite(distance):
            raise RuntimeError("Invalid separation")

        # Accumulate logarithmic growth
        sum_log += np.log(distance / epsilon)

        # Renormalise perturbation
        nearby = initial + epsilon * difference / distance

        elapsed = (i + 1) * dt

        # Running Lyapunov estimate
        estimates.append(sum_log / elapsed)

        times.append(elapsed)

    exponent = sum_log / measurement_time

    return exponent, times, estimates


# ---------------------------------
# Run Experiments
# ---------------------------------

results = []

print("Lorenz System: Lyapunov Analysis")
print("----------------------------------")

for rho in rho_values:

    exponent, times, estimates = calculate_lyapunov(rho)

    results.append((rho, exponent))

    print(
        f"Rho = {rho:5.1f} | "
        f"Largest Lyapunov Exponent = {exponent:.6f}"
    )


# ---------------------------------
# Save Results as CSV
# ---------------------------------

with open("lyapunov_rho_results.csv", "w", newline="") as file:

    writer = csv.writer(file)

    writer.writerow([
        "Rho",
        "Largest Lyapunov Exponent"
    ])

    writer.writerows(results)


# ---------------------------------
# Plot Lyapunov Exponent vs Rho
# ---------------------------------

rho_plot = [item[0] for item in results]

exponent_plot = [item[1] for item in results]

plt.figure(figsize=(10, 6))

plt.plot(
    rho_plot,
    exponent_plot,
    marker="o",
    linewidth=2,
    label="Estimated Lyapunov exponent"
)

plt.axhline(
    0,
    linestyle="--",
    label="Zero exponent reference"
)

plt.xlabel("Rayleigh Parameter (rho)")

plt.ylabel("Largest Lyapunov Exponent")

plt.title(
    "Effect of Rayleigh Parameter on Chaotic Behaviour"
)

plt.xticks(rho_values)

plt.legend()
plt.grid(True)

plt.tight_layout()

plt.savefig(
    "lyapunov_vs_rho.png",
    dpi=300
)

plt.show()
