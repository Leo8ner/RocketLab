import casadi as ca

### CONSTANTS AND PARAMETERS ###
# Physical constants
g                 = 9.81               # Earth gravity--------------------------, m/s^2
deg               = ca.pi / 180        # Degrees to radians---------------------, rad/deg
rad               = 180 / ca.pi        # Radians to degrees---------------------, deg/rad
percent           = 0.01               # Percentage-----------------------------, -

# Spacecraft parameters
m                 = 5.0e3              # Rocket mass----------------------------, kg
E_z               = 3                  # Engine distance from CoM---------------, m
J_y               = 7.5e5              # Rotational inertia---------------------, kg*m^2

# Initial conditions
x_i               = 0                  # Initial horizontal position------------, m (guess)
h_i               = 150                # Initial altitude-----------------------, m
theta_i           = 0                  # Initial pitch angle--------------------, rad
vx_i              = 0                  # Initial horizontal velocity------------, m/s
vz_i              = -10                # Initial vertical velocity--------------, m/s
wy_i              = 0                  # Initial angular velocity---------------, rad/s

# Dispersions
x_d               = 5                  # Initial horizontal position------------, m (guess)
h_d               = 10                 # Initial altitude-----------------------, m
theta_d           = 5 * deg            # Initial pitch angle--------------------, rad
vx_d              = 1                  # Initial horizontal velocity------------, m/s
vz_d              = 1                  # Initial vertical velocity--------------, m/s
wy_d              = 0                  # Initial angular velocity---------------, rad/s

# Final state conditions
x_f               = 0                  # Final horizontal position--------------, m
h_f               = 0                  # Final altitude-------------------------, m
theta_f           = 0                  # Final pitch angle----------------------, rad
vx_f              = 0                  # Final horizontal velocity--------------, m/s
vz_f              = 0                  # Final vertical velocity----------------, m/s
wy_f              = 0                  # Final angular velocity-----------------, rad/s

# Acceptable landing error
x_e               = 3                  # Final max horizontal position error ---, m 
h_e               = 0                  # Final max altitude error --------------, m
theta_e           = 2 * deg            # Final max pitch angle error -----------, rad 
vx_e              = 0.5                # Final max horizontal velocity error ---, m/s
vz_e              = 0                  # Final max vertical velocity error -----, m/s
wy_e              = 1 * deg            # Final max angular velocity error ------, rad/s

# Maximum expected state variables
x_max             = x_d                # Maximum expected horizontal position---, m
z_max             = h_i + h_d          # Maximum expected vertical position-----, m
theta_max         = theta_d            # Maximum expected pitch angle-----------, rad
vx_max            = vx_d               # Maximum expected horizontal velocity---, m/s
vz_max            = abs(vz_i) + vz_d   # Maximum expected vertical velocity-----, m/s
wy_max            = wy_e               # Maximum expected angular velocity------, rad/s

# State upper bounds - these are just conservative estimates
x_ub              = 2 * x_max          # Upper bound for horizontal position----, m
z_ub              = 2 * z_max          # Upper bound for vertical position------, m
theta_ub          = 90 * deg           # Upper bound for pitch angle------------, rad
vx_ub             = 10 * vx_max        # Upper bound for horizontal velocity----, m/s
vz_ub             = 0                  # Upper bound for vertical velocity------, m/s
wy_ub             = 90 * deg           # Upper bound for angular velocity-------, rad/s

# State lower bounds - these are just conservative estimates
x_lb              = -x_ub              # Lower bound for horizontal position----, m
z_lb              = 0                  # Lower bound for vertical position------, m
theta_lb          = -theta_ub          # Lower bound for pitch angle------------, rad
vx_lb             = -vx_ub             # Lower bound for horizontal velocity----, m/s
vz_lb             = -10 * vz_max       # Lower bound for vertical velocity------, m/s
wy_lb             = -wy_ub             # Lower bound for angular velocity-------, rad/s

# Time & solver parameters
H_0               = 15                 # Time horizon initial guess-------------, s 
H_max             = 35                 # Maximum time horizon-------------------, s
N_steps           = 50                 # Number of time steps-------------------, -
lb_dt             = 0.0001             # Minimum time step----------------------, s
ub_dt             = H_max/N_steps      # Maximum time step----------------------, s
margin            = 1e-6               # Margin for the next state--------------, - 

# Actuator limits
T_max             = 7e4                # Maximum engine thrust------------------, N
throttle          = 40 * percent       # Engine throttle------------------------, -
T_min             = throttle * T_max   # Engine minimum thrust------------------, N
throttle_rate     = 1e4                # Engine maximum throttle rate-----------, N/s
gim_max           = 15 * deg           # Max gimbal angle-----------------------, rad 
gim_min           = -gim_max           # Min gimbal angle-----------------------, rad
gim_rate          = 45 * deg           # Max gimbal rate------------------------, rad/s
w_n               = 10          # Gimbal natural frequency---------------, rad/s
zeta              = 0.9                # Gimbal damping ratio-------------------, -


# Cost function weights
W_x               = 1                  # Weight for horizontal position---------, -
W_z               = 1                  # Weight for vertical position-----------, -
W_theta           = 1                  # Weight for pitch angle-----------------, -
W_vx              = 1                  # Weight for horizontal velocity---------, -
W_vz              = 1                  # Weight for vertical velocity-----------, -
W_wy              = 1                  # Weight for angular velocity------------, -
W_T               = 1                  # Weight for engine thrust---------------, -
W_gim             = 1                  # Weight for gimbal angle----------------, -
W_T_rate          = 1                  # Weight for engine throttle rate--------, -
W_gim_rate        = 0                  # Weight for gimbal rate-----------------, -
W_H               = 10                 # Weight for time horizon----------------, -