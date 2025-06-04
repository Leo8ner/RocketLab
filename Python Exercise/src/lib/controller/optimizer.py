import casadi as ca
import numpy as np
import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.abspath(os.path.join(current_dir, '../..'))
sys.path.append(src_dir)
from lib.dynamics.dynamics import Dynamics
from lib.controller.setup import OptimizationSetup
from lib.plotter.plotter import ResultPlotter
from simulation_parameters import *  # Import simulation parameters from the file

class Optimizer:
    def __init__(self, dyn: Dynamics, setup: OptimizationSetup, guess=None):
        
        self.N_X        = dyn.get_state_dim()                              # Number of states----------------- , -
        self.N_U        = dyn.get_control_dim()                            # Number of controls----------------, -
        self.f          = dyn.get_continuous_dynamics()                    # Continuous dynamics function------, X_dot = f(X, U)
        self.F          = dyn.get_discrete_integrator()                    # Discrete dynamics function--------, X_kp1 = F(X_k, U_k, dt)
        
        cost_matrices   = setup.get_cost_matrices()                        # Cost matrices---------------------, Q, R, D, M
        self.Q          = cost_matrices['Q']                               # State cost matrix-----------------, -
        self.R          = cost_matrices['R']                               # Control cost matrix-------------- , -
        self.D          = cost_matrices['D']                               # Control rate cost matrix----------, -
        self.M          = cost_matrices['M']                               # Time horizon cost matrix----------, -
        
        conditions      = setup.get_conditions()                           # Conditions------------------------, X_i, X_f, slack
        self.X_i        = conditions['X_i']                                # Initial state---------------------, [m, m, rad, m/s, m/s, rad/s]
        self.X_f        = conditions['X_f']                                # Final state-----------------------, [m, m, rad, m/s, m/s, rad/s]
        self.slack      = conditions['slack']                              # Slack for final state constraints-, [m, m, rad, m/s, m/s, rad/s]
        
        bounds          = setup.get_bounds()                               # Bounds----------------------------, lb_X, ub_X, lb_U, ub_U, lb_dt, ub_dt
        self.U_rate_max = bounds['U_rate_max']                             # Maximum control rate--------------, [N/s, rad/s]
        self.lb_X       = bounds['lb_X']                                   # State lower bounds----------------, [m, m, rad, m/s, m/s, rad/s]
        self.ub_X       = bounds['ub_X']                                   # State upper bounds----------------, [m, m, rad, m/s, m/s, rad/s]
        self.lb_U       = bounds['lb_U']                                   # Control lower bounds--------------, [N, rad]
        self.ub_U       = bounds['ub_U']                                   # Control upper bounds--------------, [N, rad]
        self.lb_dt      = lb_dt                                            # Time horizon lower bound----------, s
        self.ub_dt      = ub_dt                                            # Time horizon upper bound----------, s
        self.U_max      = bounds['U_max']                                  # Normalized control maximum--------, [N, rad]
        if self.N_X == 8:                                                  # If actuator dynamics are included
            self.ub_X = ca.vertcat(self.ub_X, gim_max, gim_rate)           # Include gimbal angle and rate-----, [m, m, rad, m/s, m/s, rad/s, rad, rad/s]
            self.lb_X = ca.vertcat(self.lb_X, -gim_max, -gim_rate)         # Include gimbal angle and rate-----, [m, m, rad, m/s, m/s, rad/s, rad, rad/s]

        self.opti       = ca.Opti()                                        # Optimization problem
        U_mat           = ca.repmat(self.U_max, 1, N_steps)                # Control normalization matrix------, [N, rad]
        self.U          = U_mat * self.opti.variable(self.N_U, N_steps)    # Controls--------------------------, [N, rad]
        self.dt         = self.opti.variable(N_steps)                      # Initializing time step------------, s
        self.H          = ca.sum1(self.dt)                                 # Total time horizon----------------, s
        self.implicit   = dyn.is_implicit()                                # Implicit dynamics flag-------------, -
        
        self.X_max      = bounds['X_max']                                  # Normalized state maximum----------, [m, m, rad, m/s, m/s, rad/s]
        if self.N_X == 8:                                                  # If actuator dynamics are included
            self.X_max = ca.vertcat(self.X_max, ca.DM([gim_max, gim_rate]))# Include gimbal angle and rate
        X_mat           = ca.repmat(self.X_max, 1, N_steps + 1)            # State normalization matrix--------, [m, m, rad, m/s, m/s, rad/s]
        self.X_all               = X_mat * self.opti.variable(self.N_X, N_steps+1)  # States----------------------------, [m, m, rad, m/s, m/s, rad/s]
        if self.N_X == 8:                                                  # If actuator dynamics are included
            self.X = self.X_all[:-2, :]                                    # Exclude gimbal angle and rate from states
        else:
            self.X = self.X_all                                            # States without actuator dynamics, [m, m, rad, m/s, m/s, rad/s]
        
        self.guess = guess                                                 # Initial guess for the optimization problem
        self.setup_problem()

    def setup_problem(self):

        # State and control bounds
        self.opti.subject_to(self.opti.bounded(self.lb_U, self.U, self.ub_U))     # Control constraints
        self.opti.subject_to(self.opti.bounded(self.lb_X, self.X_all, self.ub_X)) # State constraints
        self.opti.subject_to(self.opti.bounded(self.lb_dt, self.dt, self.ub_dt))  # Time horizon constraints

        # Initial state constraints
        self.opti.subject_to(self.X[:, 0] == self.X_i)

        # Final state constraints
        self.opti.subject_to(self.opti.bounded(-self.slack, self.X[:, -1], self.slack))    # State constraints
        T_max_vec = self.dt * throttle_rate
        objective = 0

        for k in range(N_steps-1):
            # Control rate constraints
            if self.N_X == 8:  # If no actuator dynamics are included
                self.opti.subject_to(
                self.opti.bounded(self.U[0, k] - T_max_vec[k], self.U[0, k + 1], T_max_vec[k] + self.U[0, k])
                )
                U_rate = (self.U[:, k + 1] - self.U[:, k]) / self.dt[k]                   # Control rate
                objective += ca.bilin(self.D, U_rate)                    # Control rate cost
            X_kp1 = self.F(self.X_all[:, k], self.U[:, k], self.dt[k])                        # Next state
            self.opti.subject_to(self.X_all[:, k + 1] == X_kp1)          # Dynamics constraints RK4
            objective += ca.bilin(self.R, self.U[:, k])              # Control cost
        if self.N_X == 8:  # If no actuator dynamics are included
            U_rate = (self.U[:, -1] - self.U[:, -1]) / self.dt[-1]                   # Control rate
            objective += ca.bilin(self.D, U_rate)                    # Control rate cost
        X_kp1 = self.F(self.X_all[:, -2], self.U[:, -1], self.dt[-1])                        # Next state
        self.opti.subject_to(self.X_all[:, -1] == X_kp1)          # Dynamics constraints RK4
        objective += ca.bilin(self.R, self.U[:, -1])              # Control cost
        objective += ca.bilin(self.Q, self.X[:, -1])                 # Final state cost    
        objective += ca.bilin(self.M, self.H)                        # Time horizon cost

        if self.guess:
            self._set_initial_guess()  # Set initial guess if provided

        self.opti.minimize(objective)

    def _set_initial_guess(self):
        """
        Set an initial guess for the optimization problem.
        This is useful for speeding up convergence.
        """
        try:
            dt_guess = np.loadtxt(self.guess, delimiter=",", usecols=range(1), skiprows=1, max_rows=N_steps).T
            X_guess  = np.loadtxt(self.guess, delimiter=",", usecols=range(1, self.N_X+1), skiprows=1).T
            U_guess  = np.loadtxt(self.guess, delimiter=",", usecols=range(self.N_X+1, self.N_X+3), skiprows=1, max_rows=N_steps).T
        except FileNotFoundError:
            print("No initial guess found. Using default values.")
            X_guess  = np.zeros((self.N_X, N_steps + 1))
            U_guess  = np.zeros((self.N_U, N_steps))
            dt_guess = np.ones(N_steps) * (H_0 / N_steps)

        self.opti.set_initial(self.X_all, X_guess)  # Set initial guess for states
        self.opti.set_initial(self.U, U_guess)      # Set initial guess for controls
        self.opti.set_initial(self.dt, dt_guess)    # Set initial guess for time steps

    def solve(self, console_print=False, plot=False, save_traj=False, save_plots=False):
        """
        Solve the optimization problem and return the optimal state, control, and time horizon.

        Parameters:
        print (bool): If True, print the optimization results.
        plot (bool): If True, plot the optimization results.

        Returns:
        opt_X (np.ndarray): Optimal state trajectory.
        opt_U (np.ndarray): Optimal control trajectory.
        """

        # Set solver options
        p_opts = {
            "expand": True,  # Expand the problem to allow for larger problems
        }
        s_opts = {
            "print_level": 5,
            #"mu_strategy": "adaptive",
            "warm_start_init_point": "yes",
            "max_iter": 10000,
            #"tol": 10e-6,
            #"dual_inf_tol": 10e-6,
            #"constr_viol_tol": 10e-6,
            #"compl_inf_tol": 10e-6,
            # "acceptable_tol": self.margin,
            # "acceptable_constr_viol_tol": self.margin,
            # "acceptable_dual_inf_tol": self.margin,
            # "acceptable_compl_inf_tol": self.margin
        }
        self.opti.solver("ipopt", p_opts, s_opts)  # Choose IPOPT as the solver

        try:
            sol    = self.opti.solve()                  # Solve the optimization problem
            opt_X  = self.opti.value(self.X_all)            # Extract the optimal state trajectory
            opt_U  = self.opti.value(self.U)            # Extract the optimal control trajectory
            opt_dt = self.opti.value(self.dt)           # Extract the optimal time horizon
            opt_H  = self.opti.value(self.H)            # Total time horizon

        except RuntimeError as e:
            print(f"Solver failed: {e}")
            opt_X = self.opti.debug.value(self.X_all)  # Extract the debug values of the state trajectory
            opt_U = self.opti.debug.value(self.U)
            opt_dt = self.opti.debug.value(self.dt)
            opt_H  = self.opti.debug.value(self.H)      # Total time horizon

        # Compute state derivatives for each time step
        opt_X_dot = np.zeros((self.N_X, N_steps))
        for k in range(N_steps):
            opt_X_dot[:, k] = np.array(self.f(opt_X[:, k], opt_U[:, k])).flatten()

        result = {
            'X': opt_X,
            'U': opt_U,
            'X_dot': opt_X_dot,
            'dt': opt_dt,
            'H': opt_H
        }

        if console_print:
            print("\nManeuver duration:", round(opt_H, 2), "s")
            print("\nFinal State:")
            print(f"x: {round(opt_X[0, -1], 2)}m, z: {round(opt_X[1, -1], 2)}m, theta: {round(opt_X[2, -1]*180/ca.pi, 2)}deg, vx: {round(opt_X[3, -1], 2)}m/s, vz: {round(opt_X[4, -1], 2)}m/s, wy: {round(opt_X[5, -1]*180/ca.pi, 2)}deg/s\n")
        if plot:
            gimbal_dynanmics = self.N_X == 8
            ResultPlotter(result, "Optimal Trajectory", gim_dyn=gimbal_dynanmics, save=save_plots)  # Plot the results
        if save_traj:
            import csv
            from itertools import zip_longest
            # Save the nominal results
            file = open("results.csv", mode="w", newline="")
            write = csv.writer(file)
            if self.N_X == 6:
                # Write header
                write.writerow(["dt", "x", "z", "theta", "vx", "vz", "wy", "T", "gimbal"])

                # Fill missing values with an empty string
                rows = zip_longest(opt_dt, opt_X[0], opt_X[1], opt_X[2], opt_X[3], opt_X[4], opt_X[5], opt_U[0], opt_U[1])
            else:
                # Write header
                write.writerow(["dt", "x", "z", "theta", "vx", "vz", "wy", "gim_real", "gim_rate", "T", "gimbal"])

                # Fill missing values with an empty string
                rows = zip_longest(opt_dt, opt_X[0], opt_X[1], opt_X[2], opt_X[3], opt_X[4], opt_X[5], 
                                   opt_X[6], opt_X[7], opt_U[0], opt_U[1])
            write.writerows(rows)  # Write transposed data

        

            file.close()            
        return result
