// Libraries
#pragma once
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

// Obtained through algebraic manipulation
inline double dimensionless_velocity_derivative(double stokesNumber,
                                                double reynoldsNumber,
                                                double velocity) {
  return ((1.0 - velocity - 3.0 / 8.0 * reynoldsNumber * velocity * velocity) /
          stokesNumber);
}

inline int solver_rk4(double stokesNumber, double reynoldsNumber,
                      double initialVelocity, double stopTime, double deltaT,
                      const std::string &filename) {

  std::ofstream outputFile(filename);

  outputFile << std::setprecision(15);

  if (!outputFile.is_open()) {
    std::cerr << "Error: it was not possible to create the file" << std::endl;
    return 0; // Ends the program with an error if it was not possible to create
              // or read the file
  }

  outputFile << "Time(t*),Velocity(v*)\n";

  double velocity = initialVelocity;

  for (double t = 0; t < stopTime; t = t + deltaT) {
    outputFile << t << "," << velocity << "\n";

    double k1 = dimensionless_velocity_derivative(stokesNumber, reynoldsNumber,
                                                  velocity);
    double k2 = dimensionless_velocity_derivative(stokesNumber, reynoldsNumber,
                                                  velocity + (k1 * deltaT / 2));
    double k3 = dimensionless_velocity_derivative(stokesNumber, reynoldsNumber,
                                                  velocity + (k2 * deltaT / 2));
    double k4 = dimensionless_velocity_derivative(stokesNumber, reynoldsNumber,
                                                  velocity + (k3 * deltaT));

    velocity = velocity + (deltaT / 6) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
  }

  outputFile.close();

  return 0;
}
