#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

double f(double x) { return std::pow(x, 3.0) - 13.0 * x - 12.0; }

int main() {
  double x0 = 6.5, passo = 0.5;
  double x1 = x0 - passo;
  double x2 = x1 - passo;
  double x3 = 0.0, tol = 1.0e-6, erro = 1.0;
  int iter = 0, max_iter = 100;

  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << "                              METODO DE MULLER                 "
               "                 "
            << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::fixed << std::setprecision(4)
            << " Chutes iniciais (x0, x1, x2): " << x0 << ", " << x1 << ", "
            << x2 << std::endl;
  std::cout << " Tolerancia:                   " << std::scientific << tol
            << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << std::setw(6) << " Iter " << " | " << std::setw(15)
            << "      x3       "
            << " | " << std::setw(17) << "      f(x3)      " << " | "
            << std::setw(17) << "      Erro       " << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;

  std::ofstream outfile("raizes.dat");

  while (erro > tol && iter < max_iter) {
    iter++;
    double h0 = x1 - x0;
    double h1 = x2 - x1;
    double delta0 = (f(x1) - f(x0)) / h0;
    double delta1 = (f(x2) - f(x1)) / h1;

    double a = (delta1 - delta0) / (h1 + h0);
    double b = a * h1 + delta1;
    double c = f(x2);

    double discriminante = std::pow(b, 2.0) - 4.0 * a * c;

    if (discriminante < 0.0) {
      std::cout << "Aviso: Discriminante negativo. Raízes complexas detectadas."
                << std::endl;
      break;
    }

    double den;
    if (b >= 0.0) {
      den = b + std::sqrt(discriminante);
    } else {
      den = b - std::sqrt(discriminante);
    }

    x3 = x2 - (2.0 * c) / den;
    erro = std::abs(x3 - x2);

    std::cout << std::setw(5) << iter << " | " << std::fixed
              << std::setprecision(8) << std::setw(15) << x3 << " | "
              << std::scientific << std::setprecision(8) << std::setw(17)
              << f(x3) << " | " << std::setw(17) << erro << std::endl;

    x0 = x1;
    x1 = x2;
    x2 = x3;
  }

  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << " Raiz final encontrada: " << std::fixed << std::setprecision(8)
            << x3 << std::endl;
  std::cout << " Iteracoes realizadas:  " << iter << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;

  if (outfile.is_open()) {
    outfile << "Raiz encontrada: " << std::fixed << std::setprecision(8) << x3
            << std::endl;
    outfile.close();
  }

  return 0;
}
