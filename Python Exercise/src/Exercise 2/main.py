import casadi as ca
import numpy as np
import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(src_dir)
from lib.dynamics.dynamics import Dynamics
from simulation_parameters import *  # Import simulation parameters from the file
from lib.plotter.plotter import ResultPlotter

# Time parameters
H                 = 50                         # Maximum time horizon-------------------, s
N_steps           = 100                         # Number of time steps-------------------, -
dt                = H/N_steps                  # Time step------------------------------, s
t                 = np.linspace(0, H, N_steps) # Time vector----------------------------, s

# Initial conditions [x, z, theta, vx, vz, wy]
x_i               = 0                          # Initial horizontal position------------, m
h_i               = 0                          # Initial altitude-----------------------, m
theta_i           = 0                          # Initial pitch angle--------------------, rad
vx_i              = 0                          # Initial horizontal velocity------------, m/s
vz_i              = 0                          # Initial vertical velocity--------------, m/s
wy_i              = 0                          # Initial angular velocity---------------, rad/s

# Initial state vector
X_i     = ca.DM([x_i, h_i, theta_i, vx_i, vz_i, wy_i]) 

## Actuator input ##
# Thrust
T = ca.DM.zeros(1, N_steps)       # Engine thrust------------------, N   

# Generate thrust profile
for k in range(N_steps):
    time = t[k]
    if time < 5.0:
        T[0, k] = m * g  # weight = thrust,  so rocket doesnt go underground
    elif time < 10:
        T[0, k] = 1.5 * m * g  # Just enough thrust to liftoft 5 seconds
    elif time < 15:
        T[0, k] = m * g  # Constannt velocity
    elif time < 17.5 :
        T[0, k] = 0  # Decelerate to zero velocity
    else:
        T[0, k] = m * g # Hover

# Gimbal angle
gimbal = ca.DM.zeros(1, N_steps)  # Gimbal angle-------------------, rad

# Generate gimbal profile
for k in range(N_steps):
    time = t[k]
    if time < 25:
        gimbal[0, k] = 0  # Zero gimbal for first 25 seconds
    elif time < 30 + (2*np.pi*J_y)/(m*g * E_z * np.sin(15 * deg)):
        gimbal[0, k] = 15 * deg  # do a flip
    else:
        gimbal[0, k] = 0     

# Initialize control matrix
U_prop = ca.vertcat(T, gimbal)  # 2 controls x N_steps time points

### Dynamics setup ###
dyn   = Dynamics()  # Initialize the dynamics class
f     = dyn.f       # Continuous dynamics function
F     = dyn.F       # Discrete dynamics function
N_X   = dyn.N_X     # Number of states

X_prop = ca.DM.zeros(N_X, N_steps + 1)  # State trajectory
X_prop[:, 0] = X_i                      # Set initial state
X_dot_prop = ca.DM.zeros(N_X, N_steps)  # State derivatives trajectory

dt_prop = np.ones(N_steps) * dt  # Time step trajectory (fixed dt case)

for k in range(N_steps):
    X_prop[:, k+1] = F(X_prop[:, k], U_prop[:, k], dt)  # Next state using explicit RK4 integration
    X_dot_prop[:, k] = f(X_prop[:, k], U_prop[:, k])  # State derivatives
    if X_prop[1, k+1] < 0:
        X_prop[:, k+1] = X_prop[:, k+1]*0
        U_prop[:, k+1:] = U_prop[:, k+1:]*0  # Stop control propagation
        print("Warning: The spacecraft has crashed :/.")
        break

result = {
    'X': X_prop,
    'X_dot': X_dot_prop,
    'U': U_prop,
    'dt': dt_prop
}

### PLOTS ###
ResultPlotter(result, title="Trajectory Propagation")
