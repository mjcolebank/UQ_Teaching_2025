import numpy as np
import matplotlib.pyplot as plt

# Parameters
T = 5  # End time for numerical solution
L = 1  # Length of the 1D element
dx = 0.05  # Spatial discretization
dt = 1e-3  # Temporal discretization
x_domain = np.arange(0, L+dx, dx)  # x points
t_domain = np.arange(0, T+dt, dt)  # t points
Nx = len(x_domain)  # Number of x points
Nt = len(t_domain)  # Number of t points

# Initial and boundary conditions
IC = np.ones(Nx)
BC_0 = 4  # BC at x=0
BC_L = 5  # BC at x=L

# Thermal conductivity parameter
alpha = 0.01  # Thermal diffusivity constant (m^2/s)

# Stability check
stab_number = max(alpha) * dt / (dx ** 2)
if stab_number >= 0.5:
    raise Warning('Numerical Instability may occur.')

def heateq_1D(alpha, t_domain, x_domain, IC, BC_0, BC_L):
    Nx = len(x_domain)       # Number of spatial grid points
    Nt = len(t_domain)       # Number of time steps
    dx = x_domain[1] - x_domain[0]  # Spatial step size
    dt = t_domain[1] - t_domain[0]  # Time step size

    r = dt / (dx ** 2)  # Stability parameter

    if np.isscalar(alpha):
        alpha_x = alpha * np.ones(Nx)
    else:
        alpha_x = alpha

    # Initial condition
    u = IC.copy()

    # Boundary conditions
    u[0] = BC_0  # Left boundary (Dirichlet condition)
    u[-1] = BC_L  # Right boundary (Dirichlet condition)

    u_save = np.zeros((Nt, Nx))
    # Time-stepping loop
    for t_i in range(Nt):
        u_new = u.copy()
        for i in range(1, Nx-1):
            u_new[i] = u[i] + alpha_x[i] * r * (u[i+1] - 2*u[i] + u[i-1])
        u = u_new
        u_save[t_i, :] = u

    return u_save

# Run the 1D heat equation solver
u_save = heateq_1D(alpha, t_domain, x_domain, IC, BC_0, BC_L)

# Plot the final temperature distribution
plt.figure()
plt.plot(x_domain, u_save[0, :], '-o', label='t=0')
plt.plot(x_domain, u_save[len(t_domain)//2, :], '-o', label='t=T/2')
plt.plot(x_domain, u_save[-1, :], '-o', label='t=T')
plt.legend()
plt.xlabel('Position (m)')
plt.ylabel('Temperature (C)')
plt.title('Temperature distribution along the rod')
plt.grid(True)
plt.show()
