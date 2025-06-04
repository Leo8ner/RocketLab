#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <casadi/casadi.hpp>
#include <rlab/simulation_parameters.h>

using namespace casadi;

class Dynamics {

    bool is_implicit;
    bool actuator_dynamics;
    SX X, U, dt, X_dot;
    Function F, f; 
    int n_X, n_U;

    void setupDynamicsWithoutTVC(); // Setup dynamics without thrust vector control dynamics
    void setupDynamicsWithTVC();    // Setup dynamics with thrust vector control dynamics
    void setupIntegrator();           // Setup the integrator function

public:
    // Constructor
    Dynamics(bool is_implicit = false, bool actuator_dynamics = false); 

    // Getters
    Function getDiscreteDynamics() const;   // Returns the discrete dynamics function F
    Function getContinuousDynamics() const; // Returns the continuous dynamics function f
    int getStateDim() const;                      // Returns the number of states
    int getInputDim() const;                      // Returns the number of controls

};

// RK4 integrator
SX rk4(const SX& x_dot, const SX& x, const SX& dt);

#endif // DYNAMICS_H