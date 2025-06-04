import casadi as ca
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(src_dir)
from simulation_parameters import * 


### Initial conditions analytical feasibility check ###

## Initial linear velocity ##
# Linear momentum feasibility check
I_max = (T_max - m*g) * H_max                                      # Maximum impulse available ----------, Ns
v_f   = ca.sqrt(vz_f**2 + vx_f**2)                                 # Final linear velocity --------------, m/s
p_f   = m * v_f                                                    # Allowed final momentum -------------, Ns
p_i   = I_max - p_f                                                # Maximum allowed initial momentum ---, Ns
v_0   = p_i / m                                                    # Maximum initial linear velocity ----, m/s
# Generate angles for circle
theta = np.linspace(0, 2 * np.pi, 500)
print("Maximum initial linear velocity:", v_0, "m/s")

# Parametric equations for a circle of radius v_0 -> v_0** = vz_max**2 + vx_max**2
vx = v_0 * np.cos(theta)
vz = v_0 * np.sin(theta)

# Ignore positive vz values as they would require the vehicle to turn around
for i in vz:
    if i > 0:
        vz[vz == i] = 0

# Create figure
plt.figure(figsize=(8, 6))

# Define plot limits
x_limit = v_0 * 1.2
y_limit = v_0 * 1.2

# Create mesh for shading
x_mesh = np.linspace(-x_limit, x_limit, 100)
y_mesh = np.linspace(-y_limit, y_limit, 100)
X, Y = np.meshgrid(x_mesh, y_mesh)

# Define red shaded area for forbidden velocities
mask_forbidden = (X**2 + Y**2 > v_0**2) | (Y > 0)
plt.contourf(X, Y, mask_forbidden.astype(int), levels=[0.5, 1.5], colors=['red'], alpha=0.3)

# Plot the semicircle
plt.plot(vx, vz, 'b-', linewidth=2, label=f'v₀ = {v_0:.2f} m/s')

# Fill the allowed area in green
plt.fill(vx, vz, alpha=0.2, color='green', label='Allowed Initial Velocities')

# Create custom legend entries
from matplotlib.patches import Patch
legend_elements = [
    plt.Line2D([0], [0], color='blue', linewidth=2, label=f'v₀ = {v_0:.2f} m/s'),
    Patch(facecolor='red', alpha=0.3, label='Forbidden Initial Velocities'),
    Patch(facecolor='green', alpha=0.2, label='Allowed Initial Velocities')

]

# Formatting
plt.xlabel('vx [m/s]')
plt.ylabel('vz [m/s]')
plt.title('Maximum Initial Velocity Envelope')
plt.grid(True)
plt.xlim(-x_limit, x_limit)
plt.ylim(-y_limit, y_limit * 0.1)
plt.gca().set_aspect('equal', adjustable='box')
plt.legend(handles=legend_elements)

## Initial angular velocity ##
# Angular momentum feasibility check
tau_max = T_max * ca.sin(gim_max) * E_z                            # Maximum torque available -----------, Nm
dL_max =  tau_max * H_max                                          # Maximum impulse available ----------, Ns
L_f    = J_y * wy_f                                                # Maximum final momentum -------------, Ns
L_i    = dL_max - L_f                                              # Max allowed initial angular momentum, Ns
w_0    = L_i / J_y * rad                                           # Max initial pitch rate -------------, deg/s
print(f'Maximum initial angular velocity: {w_0:.2f} deg/s')


## Initial position ##
# Position feasibility check
# Create a grid of velocity components within the allowed semicircle
n_grid = 500
vx_range = np.linspace(-v_0, v_0, n_grid//2)
vz_range = np.linspace(-v_0, 0, n_grid//2)  # Only negative vz (downward)

# Only calculate for positive half of vx_range
vx_positive = vx_range[vx_range >= 0]
VX_pos, VZ = np.meshgrid(vx_positive, vz_range)
Z_0_MIN_pos = np.full_like(VX_pos, np.nan)

# Calculate minimum required initial altitude for positive vx values only
for i in range(VX_pos.shape[0]):
    for j in range(VX_pos.shape[1]):
        vx_i = VX_pos[i, j]
        vz_i = VZ[i, j]
            
        v_mag = np.sqrt(vx_i**2 + vz_i**2)
        
        # Check if within velocity envelope
        if v_mag > v_0:
            continue  # Leave as nan
        
        a_T = T_max / m
        
        # Define symbolic variable
        theta = ca.SX.sym('theta')
        
        # Acceleration components
        ax = - a_T * ca.sin(theta)  # Horizontal acceleration (braking)
        az = a_T * ca.cos(theta) - g  # Vertical acceleration (braking + gravity)
        
        eq = vx_i * az - vz_i * ax
        
        # Create function and solver
        f = ca.Function('f', [theta], [eq])
        
        solver = ca.rootfinder('solver', 'newton', f)
        
        guess = [np.pi/4]
        try:
            theta_sol = solver(guess)
        except Exception as e:
            print(f"Solver failed for vx={vx_i}, vz={vz_i}: {e}")
            continue
        theta_val = float(theta_sol)
        
        # Verify solution is valid
        ax_val = - a_T * np.sin(theta_val)
        az_val = a_T * np.cos(theta_val) - g
        
        # Check if times are indeed equal (within tolerance)
        t_x = abs(-vx_i / ax_val)
        t_z = abs(-vz_i / az_val)

        if t_x > H_max:  # Can't have negative time
            print(f"Skipping invalid time: t_x={t_x}, t_z={t_z} for vx={vx_i}, vz={vz_i}")
            continue

        if az_val <= 0:
            continue

        z_0_needed = -(vz_i * t_z + 0.5 * az_val * t_z**2)
            
        Z_0_MIN_pos[i, j] = z_0_needed

# Mirror the results to create the full grid
VX, VZ = np.meshgrid(vx_range, vz_range)
Z_0_MIN = np.full_like(VX, np.nan)

# Fill in the positive half
pos_indices = vx_range >= 0
Z_0_MIN[:, pos_indices] = Z_0_MIN_pos

# Mirror to negative half (excluding the zero column to avoid duplication)
neg_indices = vx_range < 0
Z_0_MIN[:, neg_indices] = np.fliplr(Z_0_MIN_pos)

# Create 3D mesh plot
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection='3d')

# Only plot feasible points
valid_mask = ~np.isnan(Z_0_MIN)
if np.any(valid_mask):
    # Create color-coded surface plot with wireframe
    surf = ax.plot_surface(VX, VZ, Z_0_MIN, 
                          cmap='viridis', 
                          alpha=0.8,
                          antialiased=True,
                          linewidth=0.5,
                          edgecolors='black')
    
    # Add color bar
    fig.colorbar(surf, ax=ax, shrink=0.6, aspect=15, label='Minimum Initial Altitude (m)')
    
else:
    print("No feasible solutions found!")

# Labels and title
ax.set_xlabel('Horizontal Velocity vx (m/s)')
ax.set_ylabel('Vertical Velocity vz (m/s)')
ax.set_zlabel('Minimum Initial Altitude (m)')
ax.set_title('Minimum Required Initial Altitude vs Initial Velocity')

# Set viewing angle
ax.view_init(elev=20, azim=45)

plt.tight_layout()

plt.show()
