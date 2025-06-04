#include <rlab/dynamics.h>

using namespace casadi;

// Constructor implementation
Dynamics::Dynamics(bool is_implicit, bool actuator_dynamics) :
    is_implicit(is_implicit), actuator_dynamics(actuator_dynamics) {

    X = SX::vertcat({
        SX::sym("x"), 
        SX::sym("z"),
        SX::sym("theta"),
        SX::sym("vx"),
        SX::sym("vz"),
        SX::sym("wy")
    });

    n_X = X.size1();

    U = SX::vertcat({
        SX::sym("T"), 
        SX::sym("gimbal")
    });

    n_U = U.size1();
    
    dt = SX::sym("dt");

    if (actuator_dynamics) {
        X = SX::vertcat({
            X, 
            SX::sym("gim_real"), 
            SX::sym("gim_dot")
        });
        n_X = X.size1();

        setupDynamicsWithTVC();
    } else {
        setupDynamicsWithoutTVC();
    }
    setupIntegrator();
}

void Dynamics::setupDynamicsWithTVC() {
    // Define the dynamics with TVC dynamics
    std::cout << "Dynamics with TVC dynamics" << std::endl;
    // State variables
    SX theta    = X(2);
    SX vx       = X(3);
    SX vz       = X(4);
    SX wy       = X(5);
    SX gim_real = X(6);
    SX gim_dot  = X(7);

    // Control variables
    SX T    = U(0); // Thrust
    SX gim  = U(1); // Gimbal angle (input reference)

    // Thrust in body frame
    SX Tz_body = T * cos(gim_real);  // Vertical thrust
    SX Tx_body = T * sin(gim_real);  // Horizontal thrust
    SX T_body  = SX::vertcat({Tx_body, Tz_body});

    // Rotation matrix
    SX Rot = SX::vertcat({
        SX::horzcat({cos(theta), -sin(theta)}),
        SX::horzcat({sin(theta),  cos(theta)})
    });

    // Thrust in inertial frame
    SX T_inert = mtimes(Rot, T_body);
    SX Tx_inert = T_inert(0);
    SX Tz_inert = T_inert(1);

    // Torque from gimbal
    SX tau = T * sin(gim_real) * E_z;

    // Dynamics
    SX x_dot      = vx;
    SX z_dot      = vz;
    SX theta_dot  = wy;
    SX vx_dot     = Tx_inert / m;
    SX vz_dot     = Tz_inert / m - g;
    SX wy_dot     = tau / J_y;
    SX gim_ddot   = -2 * zeta * w_n * gim_dot - w_n * w_n * gim_real + w_n * w_n * gim;

    // State derivative vector
    X_dot = SX::vertcat({
        x_dot,
        z_dot,
        theta_dot,
        vx_dot,
        vz_dot,
        wy_dot,
        gim_dot,
        gim_ddot
    });

    // CasADi function
    f = Function("f", {X, U}, {X_dot}, {"X", "U"}, {"X_dot"});
}

void Dynamics::setupDynamicsWithoutTVC() {
    // Define the dynamics without TVC dynamics

    // State variables
    SX theta    = X(2);
    SX vx       = X(3);
    SX vz       = X(4);
    SX wy       = X(5);

    // Control variables
    SX T    = U(0); // Thrust
    SX gim  = U(1); // Gimbal angle (input reference)

    // Thrust in body frame
    SX Tz_body = T * cos(gim);  // Vertical thrust
    SX Tx_body = T * sin(gim);  // Horizontal thrust
    SX T_body  = SX::vertcat({Tx_body, Tz_body});

    // Rotation matrix
    SX Rot = SX::vertcat({
        SX::horzcat({cos(theta), -sin(theta)}),
        SX::horzcat({sin(theta),  cos(theta)})
    });

    // Thrust in inertial frame
    SX T_inert = mtimes(Rot, T_body);
    SX Tx_inert = T_inert(0);
    SX Tz_inert = T_inert(1);

    // Torque from gimbal
    SX tau = T * sin(gim) * E_z;

    // Dynamics
    SX x_dot      = vx;
    SX z_dot      = vz;
    SX theta_dot  = wy;
    SX vx_dot     = Tx_inert / m;
    SX vz_dot     = Tz_inert / m - g;
    SX wy_dot     = tau / J_y;

    // State derivative vector
    X_dot = SX::vertcat({
        x_dot,
        z_dot,
        theta_dot,
        vx_dot,
        vz_dot,
        wy_dot,
    });

    // CasADi function
    f = Function("f", {X, U}, {X_dot}, {"X", "U"}, {"X_dot"});
}

void Dynamics::setupIntegrator() {

    Function integ;

    SXDict dae = {{"x", X}, {"u", U}, {"p", dt}, {"ode", X_dot*dt}};

    Dict opts;
    if (is_implicit) {
        // Create integrator options
        opts["collocation_scheme"] = "legendre";
        opts["interpolation_order"] = 3;
        opts["simplify"] = true;
        opts["rootfinder"] = "fast_newton";
        
        // Create the integrator function
        integ = integrator("F", "collocation", dae, opts);

    } else {
        //SX X_next = rk4(X_dot, X, dt);
        opts["simplify"] = true;
        integ = integrator("F", "rk", dae, opts);

    }
    // Set the function F to the integrator
    F = integ.map(N_steps, "serial");
}

Function Dynamics::getDiscreteDynamics() const{
    // Returns the integrator function F
    return F;
}

Function Dynamics::getContinuousDynamics() const {
    // Returns the continuous dynamics function f
    return f;
}

int Dynamics::getStateDim() const {
    // Returns the number of states
    return n_X;
}

int Dynamics::getInputDim() const {
    // Returns the number of controls
    return n_U;
}


// RK4 integrator
SX rk4(const SX& x_dot, const SX& x, const SX& dt) {
    SX k1{x_dot};
    SX k2{substitute(x_dot, x, x + dt / 2 * k1)};
    SX k3{substitute(x_dot, x, x + dt / 2 * k2)};
    SX k4{substitute(x_dot, x, x + dt * k3)};
    return x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
}