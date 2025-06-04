import casadi as ca
import sys
import os
import random
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.abspath(os.path.join(current_dir, '../..'))
sys.path.append(src_dir)
from simulation_parameters import *

class OptimizationSetup:
    """
    Handles cost function weights, normalization, constraints.
    """
    
    def __init__(self, seed=None):
        """
        Initialize the optimization setup with simulation parameters.
        
        Args:
            seed (int, optional): Random seed for reproducibility. Defaults to None.
        """  
        self.seed = seed      
        self._setup_cost_matrices()
        self._setup_normalization()
        self._setup_constraints()
    
    def _setup_cost_matrices(self):
        """Setup cost function weight matrices."""
        # Cost function weights
        self.Q_weights = [W_x, W_z, W_theta, W_vx, W_vz, W_wy]
        self.R_weights = [W_T, W_gim]
        self.D_weights = [W_T_rate, W_gim_rate]
        self.M_weights = [W_H]
        
        # Create cost matrices
        self.Q = ca.diag(ca.MX(self.Q_weights))  # State cost matrix
        self.R = ca.diag(ca.MX(self.R_weights))  # Control cost matrix
        self.D = ca.diag(ca.MX(self.D_weights))  # Control rate cost matrix
        self.M = ca.diag(ca.MX(self.M_weights))  # Time horizon cost matrix
    
    def _setup_normalization(self):
        """Setup normalization matrices for states, controls, and time horizon."""
        # State normalization
        self.X_max = ca.DM([x_max, z_max, theta_max, 
                           vx_max, vz_max, wy_max])
        X_norm = ca.diag(1 / self.X_max**2)
        self.Q_norm = ca.mtimes([self.Q, X_norm])
        
        # Control normalization
        self.U_max = ca.DM([T_max, gim_max])
        U_norm = ca.diag(1 / self.U_max**2 / N_steps)
        self.R_norm = ca.mtimes([self.R, U_norm])
        
        # Control rate normalization
        self.U_rate_max = ca.DM([throttle_rate, gim_rate])
        U_rate_norm = ca.diag(1 / self.U_rate_max**2 / N_steps)
        self.D_norm = ca.mtimes([self.D, U_rate_norm])
        
        # Time horizon normalization
        self.H_max = ca.DM([H_max])
        H_norm = ca.diag(1 / self.H_max**2)
        self.M_norm = ca.mtimes([self.M, H_norm])
    
    def _setup_constraints(self):
        """Setup optimization parameters including bounds and conditions."""
        # Initial conditions [x, z, theta, vx, vz, wy]
        self._generate_initial_conditions()

        # Final conditions [x, z, theta, vx, vz, wy]
        self.X_f = ca.DM([x_f, h_f, theta_f, vx_f, vz_f, wy_f])
        
        # Final state error variables [x, z, theta, vx, vz, wy]
        self.slack = ca.DM([x_e, h_e, theta_e, vx_e, vz_e, wy_e])
        
        # State bounds [x, z, theta, vx, vz, wy]
        self.lb_X = ca.DM([x_lb, z_lb, theta_lb, vx_lb, vz_lb, wy_lb])
        self.ub_X = ca.DM([x_ub, z_ub, theta_ub, vx_ub, vz_ub, wy_ub])
        
        # Control bounds [T, gimbal]
        self.lb_U = ca.DM([T_min, gim_min])
        self.ub_U = ca.DM([T_max, gim_max])
    
    def _generate_initial_conditions(self):
        """
        Generate randomized initial conditions with dispersions.
        
        Args:
            seed: Random seed for reproducibility (optional).
            
        Returns:
            ca.DM: Initial state vector with random dispersions
        """
        if self.seed is not None:
            random.seed(self.seed)
        
        # Random dispersions
        x_r     = x_d     * random.uniform(-1, 1)
        h_r     = h_d     * random.uniform(-1, 1)
        theta_r = theta_d * random.uniform(-1, 1)
        vx_r    = vx_d    * random.uniform(-1, 1)
        vz_r    = vz_d    * random.uniform(-1, 1)
        wy_r    = wy_d    * random.uniform(-1, 1)
        
        # Initial conditions with dispersions
        self.X_i = ca.DM([x_i + x_r, h_i + h_r, theta_i + theta_r,
                         vx_i + vx_r, vz_i + vz_r, wy_i + wy_r])
            
    def get_cost_matrices(self):
        """
        Get all cost matrices.
        
        Returns:
            dict: Dictionary containing all cost matrices
        """
        return {
            'Q': self.Q_norm,
            'R': self.R_norm,
            'D': self.D_norm,
            'M': self.M_norm
        }
    
    def get_bounds(self):
        """
        Get state and control bounds.
        
        Returns:
            dict: Dictionary containing bounds
        """
        return {
            'lb_X': self.lb_X,
            'ub_X': self.ub_X,
            'lb_U': self.lb_U,
            'ub_U': self.ub_U,
            'X_max': self.X_max,
            'U_max': self.U_max,
            'U_rate_max': self.U_rate_max
        }
    
    def get_conditions(self):
        """
        Get initial and final conditions.
        
        Returns:
            dict: Dictionary containing conditions
        """
        return {
            'X_f': self.X_f,
            'slack': self.slack,
            'X_i': self.X_i 
        }

if __name__ == "__main__":
    # Initialize the optimization setup
    opt_setup = OptimizationSetup()
    
    # Get cost matrices
    cost_matrices = opt_setup.get_cost_matrices()
    
    # Get bounds
    bounds = opt_setup.get_bounds()

    # Get conditions
    conditions = opt_setup.get_conditions()
    
    print(f"\nInitial conditions: {conditions}")
    print(f"Cost matrices: {cost_matrices}")
    print(f"Bounds: {bounds}")