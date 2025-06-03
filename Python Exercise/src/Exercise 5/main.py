import casadi as ca
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(src_dir)
from simulation_parameters import * 
from lib.dynamics.dynamics import Dynamics
from lib.controller.setup import OptimizationSetup
from lib.controller.optimizer import Optimizer


### Initial conditions analytical feasibility check ###

# Linear momentum feasibility check
v_0   = ca.sqrt(vz_max**2 + vx_max**2)                             # Maximum initial linear velocity ----, m/s
p_i   = m * v_0                                                    # Initial linear momentum ------------, Ns
v_f   = ca.sqrt(vz_f**2 + vx_f**2)                                 # Final linear velocity --------------, m/s
p_f   = m * v_f                                                    # Maximum final momentum -------------, Ns
I_req = p_i - p_f                                                  # Minimum impulse required -----------, Ns
min_T_max = I_req / H_max  + m * g                                 # Minimum required thrust to land-----, N



# Angular momentum feasibility check
w_0    = 360 * deg                                                 # Maximum initial angular velocity ---, rad/s
L_i    = J_y * w_0                                                 # Initial angular momentum -----------, Ns
w_f    = wy_f                                                      # Final angular velocity -------------, m/s
L_f    = J_y * w_f                                                 # Maximum final momentum -------------, Ns
dL_req = p_i - p_f                                                 # Minimum impulse required -----------, Ns
min_gim_max = ca.asin(dL_req / (H_max * E_z * min_T_max)) * rad    # Minimum required gimbal to land-----, deg

# Actuator dynamics
dyn = Dynamics(implicit=False, act_dyn=True)

### Controller setup ###
setup = OptimizationSetup()  # Initialize the optimization setup

# Optimization problem
opt   = Optimizer(dyn, setup, guess="guess.csv")
results = opt.solve(console_print = True, plot = True, save_traj = True, save_plots=True) # Solve the optimization problem


print(f"Minimum required maximum thrust: {min_T_max:.0f} N")
print(f"Required gimbal: {min_gim_max:.1f} deg")