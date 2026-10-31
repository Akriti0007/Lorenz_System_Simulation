
import numpy as np
import matplotlib.pyplot as plt
import csv

from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

# ---------------------------------
# Lorenz System Parameters
# ---------------------------------

sigma = 10
beta = 8 / 3

rho_values = np.linspace(0.5, 100, 300)

# Numerical settings
initial = [1.0, 1.0, 1.0]

transient_time = 100
measurement_time = 100
dt = 0.02

# Hopf bifurcation value for classical Lorenz parameters
rho_hopf = sigma * (sigma + beta + 3) / (sigma - beta - 1)


# ---------------------------------
# Lorenz Differential Equations
# ---------------------------------

def lorenz(t, state, rho):

    x, y, z = state

    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z

    return [dx, dy, dz]


# ---------------------------------
# Storage
# ---------------------------------

bifurcation_rho = []
bifurcation_z = []

equilibrium_rho = []
equilibrium_z = []


# ---------------------------------
# Bifurcation Analysis
# ---------------------------------

print("Lorenz System Bifurcation Analysis")
print("-----------------------------------")

for rho in rho_values:

    # Integrate through transient
    transient = solve_ivp(
        lorenz,
        (0, transient_time),
        initial,
        args=(rho,),
        method="DOP853",
        rtol=1e-9,
        atol=1e-11
    )

    final_state = transient.y[:, -1]

    # Integrate long-term behaviour
    t_eval = np.arange(
        0, measurement_time + dt / 2, dt
    )

    solution = solve_ivp(
        lorenz,
        (0, measurement_time),
        final_state,
        args=(rho,),
        t_eval=t_eval,
        method="DOP853",
        rtol=1e-9,
        atol=1e-11
    )

    z = solution.y[2]

    # Find local maxima of z
    peaks, _ = find_peaks(z, prominence=1e-5)

    for peak in peaks:

        bifurcation_rho.append(rho)
        bifurcation_z.append(z[peak])

    # Equilibrium branch for rho below Hopf point
    if rho <= 1:

        equilibrium_rho.append(rho)
        equilibrium_z.append(0.0)

    elif rho < rho_hopf:

        equilibrium_rho.append(rho)
        equilibrium_z.append(rho - 1)


    if int(rho * 10) % 100 == 0:

        print(f"Processing rho = {rho:.2f}")


# ---------------------------------
# Save Numerical Results
# ---------------------------------

with open(
    "bifurcation_results.csv",
    "w",
    newline=""
) as file:

    writer = csv.writer(file)

    writer.writerow(["Rho", "Local Maximum of Z"])

    writer.writerows(
        zip(bifurcation_rho, bifurcation_z)
    )

print("\nBifurcation analysis completed.")

print("Hopf bifurcation parameter:", rho_hopf)


# ---------------------------------
# Plot Bifurcation Diagram
# ---------------------------------

plt.figure(figsize=(12, 7))

plt.scatter(
    bifurcation_rho,
    bifurcation_z,
    s=2,
    label="Long-term local maxima of Z"
)

plt.plot(
    equilibrium_rho,
    equilibrium_z,
    linestyle="--",
    linewidth=1.5,
    label="Stable equilibrium branch"
)

plt.axvline(
    rho_hopf,
    linestyle=":",
    label=f"Hopf bifurcation (rho={rho_hopf:.2f})"
)

plt.xlabel("Rayleigh Parameter (rho)")

plt.ylabel("Local Maxima of Z")

plt.title(
    "Lorenz System: Bifurcation Diagram"
)

plt.legend()
plt.grid(True)

plt.tight_layout()

plt.savefig(
    "bifurcation_diagram.png",
    dpi=300
)

plt.show()

