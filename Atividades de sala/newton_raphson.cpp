#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

double f(double x) {
  return std::pow(x, 3.0) - 6.0 * std::pow(x, 2.0) + 9.0 * x - 4.0;
}

double f_derivada1(double x) { return 3.0 * std::pow(x, 2.0) - 12.0 * x + 9.0; }

double f_derivada2(double x) { return 6.0 * x - 12.0; }

double f_NR(double x) { return x - (f(x) / f_derivada1(x)); }

double f_NRmod(double x) {
  return x - ((f(x) * f_derivada1(x)) /
              (std::pow(f_derivada1(x), 2.0) - (f(x) * f_derivada2(x))));
}

int main() {
  double x_atual = 100;
  double tol = 1.0e-6;
  int max_iter = 100, iter = 0;
  double erro = 1.0;

  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << "                       METODO DE NEWTON-RAPHSON                "
               "                 "
            << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << " Chute inicial: " << std::fixed << std::setprecision(4)
            << x_atual << std::endl;
  std::cout << " Tolerancia:    " << std::scientific << tol << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << std::setw(6) << " Iter " << " | " << std::setw(15)
            << "    x_atual     " << " | " << std::setw(17)
            << "      f(x)       "
            << " | " << std::setw(17) << "      Erro       " << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;

  while (erro > tol && iter < max_iter) {
    iter++;
    double fx = f(x_atual);
    double df = f_derivada1(x_atual);

    double x_novo = f_NR(x_atual); // or f_NRmod(x_atual)
    erro = std::abs(x_novo - x_atual);

    std::cout << std::setw(5) << iter << " | " << std::fixed
              << std::setprecision(8) << std::setw(15) << x_atual << " | "
              << std::scientific << std::setprecision(8) << std::setw(17) << fx
              << " | " << std::setw(17) << erro << std::endl;

    if (std::abs(df) < 1e-12) {
      std::cout << "Aviso: Derivada próxima de zero. Parando." << std::endl;
      break;
    }

    x_atual = x_novo;
  }

  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << " Raiz encontrada:      " << std::fixed << std::setprecision(8)
            << x_atual << std::endl;
  std::cout << " Iteracoes realizadas: " << iter << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;

  return 0;
}
