#include <cmath>
#include <iomanip>
#include <iostream>

double p(double ra, double v, double r1, double r2, double r3) {
  double num = v * r3 * ra;
  double den = r1 * (ra + r2 + r3) + r3 * ra + r3 * r2;
  return std::pow(num / den, 2.0) / ra;
}

int main() {
  double v = 80.0, r1 = 8.0, r2 = 12.0, r3 = 10.0;
  double tol = 1.0e-4, erro = 1.0;
  double ra_0 = 10.0, ra_1 = 39.0, ra_2 = 91.0, ra_3, ra_anterior;
  double p_0, p_1, p_2, p_3;

  if (ra_0 == 0.0)
    ra_0 += 1.0e-4;

  auto calculate_ra3 = [&](double r0, double r1, double r2) {
    double p0 = p(r0, v, r1, r2, r3);
    double p1 = p(r1, v, r1, r2, r3);
    double p2 = p(r2, v, r1, r2, r3);
    double num = p0 * (std::pow(r1, 2.0) - std::pow(r2, 2.0)) +
                 p1 * (std::pow(r2, 2.0) - std::pow(r0, 2.0)) +
                 p2 * (std::pow(r0, 2.0) - std::pow(r1, 2.0));
    double den =
        p0 * (r1 - r2) * 2.0 + p1 * (r2 - r0) * 2.0 + p2 * (r0 - r1) * 2.0;
    return num / den;
  };

  ra_3 = calculate_ra3(ra_0, ra_1, ra_2);
  p_0 = p(ra_0, v, r1, r2, r3);
  p_1 = p(ra_1, v, r1, r2, r3);
  p_2 = p(ra_2, v, r1, r2, r3);
  p_3 = p(ra_3, v, r1, r2, r3);

  while (std::abs(erro) > tol) {
    ra_anterior = ra_3;
    double ra_lower, ra_mid, ra_upper;

    if (ra_3 > ra_1) {
      if (p_3 > p_1) {
        ra_lower = ra_1;
        ra_mid = ra_3;
        ra_upper = ra_2;
      } else {
        ra_lower = ra_0;
        ra_mid = ra_1;
        ra_upper = ra_3;
      }
    } else {
      if (p_3 > p_1) {
        ra_lower = ra_0;
        ra_mid = ra_3;
        ra_upper = ra_1;
      } else {
        ra_lower = ra_3;
        ra_mid = ra_1;
        ra_upper = ra_2;
      }
    }

    ra_0 = ra_lower;
    ra_1 = ra_mid;
    ra_2 = ra_upper;
    p_0 = p(ra_0, v, r1, r2, r3);
    p_1 = p(ra_1, v, r1, r2, r3);
    p_2 = p(ra_2, v, r1, r2, r3);

    ra_3 = calculate_ra3(ra_0, ra_1, ra_2);
    p_3 = p(ra_3, v, r1, r2, r3);

    erro = (ra_3 - ra_anterior) / ra_3;
  }

  std::cout << std::fixed << std::setprecision(8) << "Ra_3: " << ra_3
            << " | P_3: " << p_3 << std::endl;

  return 0;
}
