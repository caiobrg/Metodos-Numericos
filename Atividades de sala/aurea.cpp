#include <cmath>
#include <iomanip>
#include <iostream>

double p(double ra, double v, double r1, double r2, double r3) {
  double num = v * r3 * ra;
  double den = r1 * (ra + r2 + r3) + r3 * ra + r3 * r2;
  return std::pow(num / den, 2.0) / ra;
}

int main() {
  double ra_opt, p_max, ra1, ra2;
  double r1 = 8.0, r2 = 12.0, r3 = 10.0;
  double v = 80.0;
  double razao = (std::sqrt(5.0) - 1.0) / 2.0;
  double ra_upper = 100.0;
  double ra_lower = 0.0;
  double length_1;
  double tol = 1.0e-4;
  double erro = 1.0;
  int iter = 0;

  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << "                      OTIMIZACAO: METODO DA RAZAO AUREA        "
               "                 "
            << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << std::setw(6) << " Iter " << " | " << std::setw(14)
            << "   ra_lower    " << " | " << std::setw(14) << "   ra_upper    "
            << " | " << std::setw(14) << "  ra_otimizado " << " | "
            << std::setw(14) << "      Erro      " << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;

  while (std::abs(erro) > tol) {
    iter++;
    length_1 = (ra_upper - ra_lower) * razao;
    ra1 = ra_lower + length_1;
    ra2 = ra_upper - length_1;

    if (p(ra2, v, r1, r2, r3) > p(ra1, v, r1, r2, r3)) {
      ra_upper = ra1;
    } else {
      ra_lower = ra2;
    }

    erro = ra1 - ra2;
    ra_opt = (ra_upper + ra_lower) / 2.0;

    std::cout << std::setw(6) << iter << " | " << std::fixed
              << std::setprecision(8) << std::setw(14) << ra_lower << " | "
              << std::setw(14) << ra_upper << " | " << std::setw(14) << ra_opt
              << " | " << std::scientific << std::setprecision(6)
              << std::setw(14) << erro << std::endl;
  }

  ra_opt = (ra_upper + ra_lower) / 2.0;
  p_max = p(ra_opt, v, r1, r2, r3);

  std::cout << "---------------------------------------------------------------"
               "--------------------"
            << std::endl;
  std::cout << "Ra otimizado: " << std::fixed << std::setprecision(8) << ra_opt
            << std::endl;
  std::cout << "Potência máxima: " << p_max << std::endl;
  std::cout << "---------------------------------------------------------------"
               "--------------------"
            << std::endl;

  return 0;
}
