#include <casadi/casadi.hpp>
#include <rlab/controller.h>
#include <rlab/dynamics.h>
#include <iostream>
#include <chrono>
#include <cstdlib>

using namespace casadi;

int main() {
    // Start the timer
    // This is used to measure the time taken by the optimization process
    auto start = std::chrono::high_resolution_clock::now();

    // Dynamics
    Dynamics dyn(false, true); // Create an instance of the dynamics class
    
    // Constraints
    OptimizationSetup setup; // Create an instance of the Constraints class

    Optimizer opti(dyn, setup, "../output/guess.csv"); // Create an instance of the optimizer class

    // Call the solver
    DMDict res = opti.solve(true, true, true, true); // Solve the optimization problem

    // Stop the timer
    auto end = std::chrono::high_resolution_clock::now();
    // Calculate the elapsed time
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start) / 1000.0;
    // Print the elapsed time
    std::cout << "Computation Time: " << elapsed.count() << " s" << std::endl;
    std::cout << "Maneuver Duration: " << res["H"] << " s" << std::endl;
    
    return 0;
}