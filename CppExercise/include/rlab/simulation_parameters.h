/* This header file contains the definition of the parameters for the optimization problem */
#ifndef SIMULATIONPARAMETERS_H
#define SIMULATIONPARAMETERS_H

/// CONSTANTS AND PARAMETERS ///

// Constants
const double PI         = 3.141592653589793; // Pi------------------------------, -
const double DEG        = PI / 180.0;        // Degrees to radians conversion---, rad/deg
const double RAD        = 180.0 / PI;        // Radians to degrees conversion---, deg/rad
const double g          = 9.81;              // Earth gravity-------------------, m/s^2
const double PERCENT    = 0.01;              // Percentage----------------------, -

// Spacecraft parameters
const double m          = 5.0e3;             // Rocket mass---------------------, kg
const double E_z        = 3.0;               // Engine distance from CoM--------, m
const double J_y        = 7.5e5;             // Rotational inertia--------------, kg*m^2

// Initial conditions
const double x_i        = 0.0;               // Initial horizontal position-----, m
const double h_i        = 150.0;             // Initial altitude----------------, m
const double theta_i    = 0.0;               // Initial pitch angle-------------, rad
const double vx_i       = 0.0;               // Initial horizontal velocity-----, m/s
const double vz_i       = -10.0;             // Initial vertical velocity-------, m/s
const double wy_i       = 0.0;               // Initial angular velocity--------, rad/s

// Dispersions
const double x_d        = 5.0;               // Initial horizontal position-----, m
const double h_d        = 10.0;              // Initial altitude----------------, m
const double theta_d    = 5.0 * DEG;         // Initial pitch angle-------------, rad
const double vx_d       = 1.0;               // Initial horizontal velocity-----, m/s
const double vz_d       = 1.0;               // Initial vertical velocity-------, m/s
const double wy_d       = 0.0;               // Initial angular velocity--------, rad/s

// Final state conditions
const double x_f        = 0.0;               // Final horizontal position-------, m
const double h_f        = 0.0;               // Final altitude------------------, m
const double theta_f    = 0.0;               // Final pitch angle---------------, rad
const double vx_f       = 0.0;               // Final horizontal velocity-------, m/s
const double vz_f       = 0.0;               // Final vertical velocity---------, m/s
const double wy_f       = 0.0;               // Final angular velocity----------, rad/s

// Acceptable landing error
const double x_e        = 3.0;               // Max horizontal error------------, m
const double h_e        = 0.0;               // Max altitude error--------------, m
const double theta_e    = 2.0 * DEG;         // Max pitch angle error-----------, rad
const double vx_e       = 0.5;               // Max horizontal velocity error---, m/s
const double vz_e       = 0.0;               // Max vertical velocity error-----, m/s
const double wy_e       = 1.0 * DEG;         // Max angular velocity error------, rad/s

// Maximum expected state variables
const double x_max      = x_d;               // Max horizontal position---------, m
const double z_max      = h_i + h_d;         // Max vertical position-----------, m
const double theta_max  = theta_d;           // Max pitch angle-----------------, rad
const double vx_max     = vx_d;              // Max horizontal velocity---------, m/s
const double vz_max     = -vz_i + vz_d;      // Max vertical velocity-----------, m/s
const double wy_max     = wy_e;              // Max angular velocity------------, rad/s 

// State upper bounds
const double x_ub       = 2.0 * x_max;       // Upper bound horizontal position-, m
const double z_ub       = 2.0 * z_max;       // Upper bound vertical position---, m
const double theta_ub   = 90.0 * DEG;        // Upper bound pitch angle---------, rad
const double vx_ub      = 10.0 * vx_max;     // Upper bound horizontal velocity-, m/s
const double vz_ub      = 0.0;               // Upper bound vertical velocity---, m/s
const double wy_ub      = 90.0 * DEG;        // Upper bound angular velocity----, rad/s

// State lower bounds
const double x_lb       = -x_ub;             // Lower bound horizontal position-, m
const double z_lb       = 0.0;               // Lower bound vertical position---, m
const double theta_lb   = -theta_ub;         // Lower bound pitch angle---------, rad
const double vx_lb      = -vx_ub;            // Lower bound horizontal velocity-, m/s
const double vz_lb      = -10.0 * vz_max;    // Lower bound vertical velocity---, m/s
const double wy_lb      = -wy_ub;            // Lower bound angular velocity----, rad/s

// Time & solver parameters
const double H_0        = 15.0;              // Initial time horizon------------, s
const double H_max      = 35.0;              // Max time horizon----------------, s
const int    N_steps    = 50;               // Time steps----------------------, -
const double lb_dt      = 0.001;            // Min time step-------------------, s
const double ub_dt      = H_max / N_steps;   // Max time step-------------------, s
const double margin     = 1e-6;              // Margin for next state-----------, -

// Actuator limits
const double T_max      = 7e4;               // Max engine thrust---------------, N
const double throttle   = 40.0 * PERCENT;    // Engine throttle-----------------, -
const double T_min      = throttle * T_max;  // Min engine thrust---------------, N
const double throt_rate = 1e4;               // Max throttle rate---------------, N/s
const double gim_max    = 15.0 * DEG;        // Max gimbal angle----------------, rad
const double gim_min    = -gim_max;          // Min gimbal angle----------------, rad
const double gim_rate   = 45.0 * DEG;        // Max gimbal rate-----------------, rad/s
const double w_n        = 10.0;              // Gimbal natural frequency--------, rad/s
const double zeta       = 0.9;               // Gimbal damping ratio------------, -

// Cost function weights
const double W_x        = 1.0;               // Weight for horizontal position--, -
const double W_z        = 1.0;               // Weight for vertical position----, -
const double W_theta    = 1.0;               // Weight for pitch angle----------, -
const double W_vx       = 1.0;               // Weight for horizontal velocity--, -
const double W_vz       = 1.0;               // Weight for vertical velocity----, -
const double W_wy       = 1.0;               // Weight for angular velocity-----, -
const double W_T        = 1.0;               // Weight for engine thrust--------, -
const double W_gim      = 1.0;               // Weight for gimbal angle---------, -
const double W_T_rate   = 1.0;               // Weight for engine throttle rate-, -
const double W_gim_rate = 0.0;               // Weight for gimbal rate----------, -
const double W_H        = 10.0;              // Weight for time horizon---------, -

#endif // SIMULATIONPARAMETERS_H
