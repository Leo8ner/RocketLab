import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(src_dir)
from simulation_parameters import * 
from lib.dynamics.dynamics import Dynamics
from lib.controller.setup import OptimizationSetup
from lib.controller.optimizer import Optimizer

### Dynamics setup ###
dyn = Dynamics(implicit=False, act_dyn=False) # Initialize the dynamics class

### Controller setup ###
setup = OptimizationSetup(seed=10)  # Initialize the optimization setup

# Optimization problem
opt   = Optimizer(dyn, setup, guess="guess.csv")  # Initialize the optimizer with the dynamics and setup
results = opt.solve(console_print = True, plot = True, save_traj=True, save_plots=True) # Solve the optimization problem