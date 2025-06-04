import casadi as ca
import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.abspath(os.path.join(current_dir, '../..'))
sys.path.append(src_dir)
from simulation_parameters import *

class Dynamics:
    """
    Rocket dynamics class.
    
    Provides both continuous dynamics function f(X, U) -> X_dot
    and discrete integrator F(X, U, dt) -> X_next
    """
    
    def __init__(self, implicit=False, act_dyn=False):
        """Initialize rocket dynamics with physical parameters."""
        
        # Initialize symbolic variables and dynamics
        self.implicit = implicit  # Flag for implicit dynamics
        self.act_dyn = act_dyn    # Flag for actuator dynamics
        self._setup_symbolic_variables()
        if self.act_dyn:
            self._setup_dynamics_with_TVC()
        else:
            self._setup_dynamics()
        self._setup_integrator()
    
    def _setup_symbolic_variables(self):
        """Setup symbolic state and control variables."""
        # States
        self.X = ca.vertcat(
            ca.SX.sym('x'),          # Horizontal position--, m
            ca.SX.sym('z'),          # Vertical position----, m
            ca.SX.sym('theta'),      # Pitch angle----------, rad
            ca.SX.sym('vx'),         # Horizontal velocity--, m/s
            ca.SX.sym('vz'),         # Vertical velocity----, m/s
            ca.SX.sym('wy'),         # Angular velocity-----, rad/s
        )
        if self.act_dyn:
            self.X = ca.vertcat(
                self.X, 
                ca.SX.sym('gim_real'), # Real gimbal angle--, rad
                ca.SX.sym('gim_dot')   # Gimbal angle rate--, rad/s
            )
        self.N_X = self.X.size1()    # Number of states-----, -

        # Controls
        self.U = ca.vertcat(
            ca.SX.sym('T'),          # Engine thrust--------, N
            ca.SX.sym('gimbal'),     # Gimbal angle---------, rad
        )
        self.N_U = self.U.size1()    # Number of controls---, -

        # Time step
        self.dt = ca.SX.sym('dt')    # Time step------------, s
    
    def _setup_dynamics(self):
        """Setup continuous dynamics equations."""
        # Extract state variables
        theta, vx, vz, wy = self.X[2], self.X[3], self.X[4], self.X[5]
        
        # Extract control variables
        T = self.U[0]                # Engine thrust----------------------, N              
        gimbal = self.U[1]           # Gimbal angle-----------------------, rad                 
        
        # Thrust in body frame
        Tz_body = T * ca.cos(gimbal) # Vertical thrust in body frame------, N  
        Tx_body = T * ca.sin(gimbal) # Horizontal thrust in body frame----, N
        T_body = ca.vertcat(Tx_body, Tz_body)  # Total thrust in body frame---------, N

        # Rotation matrix for body-to-inertial frame transformation
        Rot = ca.vertcat(
            ca.horzcat(ca.cos(theta), -ca.sin(theta)), 
            ca.horzcat(ca.sin(theta),  ca.cos(theta))
        )

        # Thrust in inertial frame
        T_inert = ca.mtimes(Rot, T_body)       # Total thrust in inertial frame-----, N
        Tx_inert = T_inert[0]                  # Horizontal thrust in inertial frame, N
        Tz_inert = T_inert[1]                  # Vertical thrust in inertial frame--, N

        # Torque
        tau = T * ca.sin(gimbal) * E_z         # Torque due to gimbal---------------, Nm

        # State derivatives
        x_dot = vx                             # Horizontal velocity----------------, m/s 
        z_dot = vz                             # Vertical velocity------------------, m/s   
        theta_dot = wy                         # Angular velocity-------------------, rad/s  
        vx_dot = Tx_inert / m                  # Horizontal acceleration------------, m/s^2
        vz_dot = Tz_inert / m - g              # Vertical acceleration--------------, m/s^2
        wy_dot = tau / J_y                     # Angular acceleration---------------, rad/s^2

        # Complete state derivative vector
        self.X_dot = ca.vertcat(x_dot, z_dot, theta_dot, vx_dot, vz_dot, wy_dot)

        # Create continuous dynamics function
        self.f = ca.Function('f', [self.X, self.U], [self.X_dot], 
                            ["X", "U"], ["X_dot"])
        
    def _setup_dynamics_with_TVC(self):
        """Setup continuous dynamics equations."""
        # Extract state variables
        theta, vx, vz, wy = self.X[2], self.X[3], self.X[4], self.X[5]
        gim_real, gim_dot = self.X[6], self.X[7]
        # Extract control variables
        T = self.U[0]                # Engine thrust----------------------, N              
        gimbal = self.U[1]           # Gimbal angle-----------------------, rad                 
        
        # Thrust in body frame
        Tz_body = T * ca.cos(gim_real) # Vertical thrust in body frame------, N  
        Tx_body = T * ca.sin(gim_real) # Horizontal thrust in body frame----, N
        T_body = ca.vertcat(Tx_body, Tz_body)  # Total thrust in body frame---------, N

        # Rotation matrix for body-to-inertial frame transformation
        Rot = ca.vertcat(
            ca.horzcat(ca.cos(theta), -ca.sin(theta)), 
            ca.horzcat(ca.sin(theta),  ca.cos(theta))
        )

        # Thrust in inertial frame
        T_inert = ca.mtimes(Rot, T_body)       # Total thrust in inertial frame-----, N
        Tx_inert = T_inert[0]                  # Horizontal thrust in inertial frame, N
        Tz_inert = T_inert[1]                  # Vertical thrust in inertial frame--, N

        # Torque
        tau = T * ca.sin(gim_real) * E_z         # Torque due to gimbal---------------, Nm

        # State derivatives
        x_dot = vx                             # Horizontal velocity----------------, m/s 
        z_dot = vz                             # Vertical velocity------------------, m/s   
        theta_dot = wy                         # Angular velocity-------------------, rad/s  
        vx_dot = Tx_inert / m                  # Horizontal acceleration------------, m/s^2
        vz_dot = Tz_inert / m - g              # Vertical acceleration--------------, m/s^2
        wy_dot = tau / J_y                     # Angular acceleration---------------, rad/s^2
        gim_ddot = -2 * zeta * w_n * gim_dot - w_n**2 * gim_real + w_n**2 * gimbal

        # Complete state derivative vector
        self.X_dot = ca.vertcat(x_dot, z_dot, theta_dot, vx_dot, vz_dot, wy_dot, gim_dot, gim_ddot)

        # Create continuous dynamics function
        self.f = ca.Function('f', [self.X, self.U], [self.X_dot], 
                            ["X", "U"], ["X_dot"])
    
    def _setup_integrator(self):
        """Setup discrete integrator using collocation method."""
        # Define DAE for integrator
        if self.implicit:
            dae = {
                'x': self.X,
                'u': self.U,
                'p': self.dt,
                'ode': self.X_dot * self.dt
            }
            
            # Create integrator
            self.F = ca.integrator('F', 'collocation', dae)
        else:
            # Explicit integrator for non-implicit case
            X_next = self._rk4_explicit(self.f, self.X, self.U, self.dt)
            self.F = ca.Function('F', [self.X, self.U, self.dt], [X_next],
                                ["X", "U", "dt"], ["X_next"])
    
    def _rk4_explicit(self, f, x, u, h):
        """ One step of ode integration with Runge-Kutta of order 4

        Args:
            f:  implementation of an ode, such that dx/dt = f(x, u)
            x:  initial value
            u:  control input
            h: time step
        """
        k1       = f(x, u)
        k2       = f(x + h/2 * k1, u)
        k3       = f(x + h/2 * k2, u)
        k4       = f(x + h * k3, u)
        return x + h/6 * (k1 + 2*k2 + 2*k3 + k4)


    def get_continuous_dynamics(self):
        """
        Get continuous dynamics function.
        
        Returns:
            ca.Function: f(X, U) -> X_dot
        """
        return self.f
    
    def get_discrete_integrator(self):
        """
        Get discrete integrator function.
        
        Returns:
            ca.Function: F(X, U, dt) -> X_next
        """
        return self.F
    
    def get_state_dim(self):
        """Get state dimension."""
        return self.N_X
    
    def get_control_dim(self):
        """Get control dimension."""
        return self.N_U
    
    def is_implicit(self):
        """Check if the dynamics are implicit.
        Returns:
            bool: True if implicit, False otherwise
        """
        return self.implicit


if __name__ == "__main__":
    # Create rocket dynamics instance
    rocket = Dynamics()
    
    # Get dynamics functions
    f = rocket.get_continuous_dynamics()
    F = rocket.get_discrete_integrator()
    
    print(f"Continuous dynamics function: {f}")
    print(f"Discrete integrator: {F}")