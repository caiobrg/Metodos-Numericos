#include "solver_rk4.hpp"
#include <filesystem>
#include <iostream>
#include <string>

int main() {
  // Create the results directory (if it doesn't already exists)
  std::filesystem::create_directory("results/St");
  std::filesystem::create_directory("results/Dt");
  std::filesystem::create_directory("results/Re");

  double stopTime = 10.0;       // amount of time* to simulate for
  double inicialVelocity = 0.0; // Starting velocity* of the sphere

  std::cout << "Starting simulations..." << std::endl;

  // ---------------------------------------------------------
  // 1: Re -> 0 for diferent values of Stokes' number
  // ---------------------------------------------------------
  double standardDeltaT = 0.01;
  double stokes_values[] = {0.1, 0.5, 1.0, 1.5, 2.0};

  for (double St : stokes_values) {
    std::string filename =
        "results/St/resultados_St_" + std::to_string(St) + ".csv";
    solver_rk4(St, 0.0, inicialVelocity, stopTime, standardDeltaT, filename);
  }

  // ---------------------------------------------------------
  // 2: Influence of the size of deltaT (discretization)
  // ---------------------------------------------------------
  double St = 1.0;
  double Re = 0.5;
  double deltaT_values[] = {2.0, 0.5, 0.1, 0.05, 0.01};

  for (double dt : deltaT_values) {
    std::string filename =
        "results/Dt/resultados_dt_" + std::to_string(dt) + ".csv";
    solver_rk4(St, Re, inicialVelocity, stopTime, dt, filename);
  }

  // ---------------------------------------------------------
  // 3, 4 e 5:Inertia effects  (Re != 0)
  // ---------------------------------------------------------
  double reynolds_values[] = {0, 0.1, 0.5, 1.0};

  for (double Re : reynolds_values) {
    std::string filename =
        "results/Re/resultados_Re_" + std::to_string(Re) + ".csv";
    solver_rk4(St, Re, inicialVelocity, stopTime, standardDeltaT, filename);
  }

  std::cout << "Simulations completed" << std::endl;
  return 0;
}
