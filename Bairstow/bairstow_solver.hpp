#pragma once
#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

/**
 * @brief Estrutura para armazenar o resultado da resolução do Método de
 * Bairstow.
 */
struct bairstowResult {
  std::vector<std::complex<double>>
      roots; ///< Vetor contendo as raízes encontradas (reais ou complexas).
  int total_iterations; ///< Número total de iterações realizadas em todo o
                        ///< processo.
  bool success;         ///< Indica se o processo convergiu dentro do limite de
                        ///< iterações.
  double r_converged;   ///< Valor final de r para o primeiro fator quadrático
                        ///< encontrado.
  double s_converged;   ///< Valor final de s para o primeiro fator quadrático
                        ///< encontrado.
  bool
      exact_root; ///< Se encontrou o par de raízes (r,s) em 1 iteração ou menos
};

/**
 * @brief Calcula as raízes de uma equação quadrática ax^2 + bx + c = 0.
 *
 * @param a Coeficiente de x^2.
 * @param b Coeficiente de x.
 * @param c Termo constante.
 * @return std::vector<std::complex<double>> Vetor com as duas raízes.
 */
inline std::vector<std::complex<double>> get_quadratic_roots(double a, double b,
                                                             double c) {
  std::vector<std::complex<double>> roots;

  // Cálculo do discriminante: Delta = b^2 - 4ac
  double discriminant = b * b - 4.0 * a * c;

  std::complex<double> root1, root2;

  if (discriminant >= 0) {
    // Raízes reais
    double sqrt_D = std::sqrt(discriminant);
    root1 = std::complex<double>((-b + sqrt_D) / (2.0 * a), 0.0);
    root2 = std::complex<double>((-b - sqrt_D) / (2.0 * a), 0.0);
  } else {
    // Raízes complexas conjugadas
    double sqrt_abs_D = std::sqrt(std::abs(discriminant));
    double real_part = -b / (2.0 * a);
    double imag_part = sqrt_abs_D / (2.0 * a);

    root1 = std::complex<double>(real_part, imag_part);
    root2 = std::complex<double>(real_part, -imag_part);
  }
  roots.push_back(root1);
  roots.push_back(root2);
  return roots;
}

/**
 * @brief Implementa o Método de Bairstow para encontrar todas as raízes de um
 * polinômio.
 *
 * O método busca fatores quadráticos da forma (x^2 - rx - s).
 *
 * @param aCoeffs Vetores de coeficientes do polinômio (do termo constante até o
 * grau n).
 * @param r_guess Chute inicial para r.
 * @param s_guess Chute inicial para s.
 * @param tolerance Tolerância para o critério de parada.
 * @param max_iterations Número máximo de iterações por fator quadrático.
 * @return bairstowResult Estrutura com raízes e metadados da execução.
 */
inline bairstowResult bairstowSolve(std::vector<double> aCoeffs, double r_guess,
                                    double s_guess, double tolerance,
                                    int max_iterations) {

  // Vetor de raízes final
  std::vector<std::complex<double>> roots;

  // Vetor temporário para raízes de cada fator quadrático
  std::vector<std::complex<double>> localRoots;

  // Inicialização das variáveis r e s com os chutes iniciais
  double r = r_guess;
  double s = s_guess;

  // Grau atual do polinômio (número de coeficientes - 1)
  int n = aCoeffs.size() - 1;

  // Variáveis de incremento para o método de Newton-Raphson
  double delS;
  double delR;

  int totalIterations = 0;
  bool process_success = true;

  // Variáveis para rastrear a primeira bacia de convergência (usado no fractal)
  bool first_factor_found = false;
  double first_r_conv = 0.0;
  double first_s_conv = 0.0;
  bool exact_root = false;

  // Enquanto o grau do polinômio for maior que 2, busca-se um fator quadrático
  while (n > 2) {
    double err = tolerance + 1;
    int iteration = 0;

    // Vetores para a divisão sintética (b) e para as derivadas (c)
    std::vector<double> bCoeffs(aCoeffs.size(), 0.0);
    std::vector<double> cCoeffs(aCoeffs.size(), 0.0);

    // Loop de refinamento de r e s (Método de Newton-Raphson para duas
    // variáveis)
    while (err > tolerance && iteration < max_iterations) {
      totalIterations++;
      iteration++;

      // b_n = a_n
      bCoeffs[n] = aCoeffs[n];

      // b_{n-1} = a_{n-1} + r * b_n
      bCoeffs[n - 1] = aCoeffs[n - 1] + r * bCoeffs[n];

      // c_n = b_n
      cCoeffs[n] = bCoeffs[n];

      // c_{n-1} = b_{n-1} + r * c_n
      cCoeffs[n - 1] = bCoeffs[n - 1] + r * cCoeffs[n];

      for (int i = n - 2; i >= 0; i--) {
        // b_i = a_i + r*b_{i+1} + s*b_{i+2}
        bCoeffs[i] = aCoeffs[i] + r * bCoeffs[i + 1] + s * bCoeffs[i + 2];

        if (i > 0) {
          // c_i = b_i + r*c_{i+1} + s*c_{i+2}
          cCoeffs[i] = bCoeffs[i] + r * cCoeffs[i + 1] + s * cCoeffs[i + 2];
        }
      }

      // O sistema linear para delR e delS é:
      // c2 * delR + c3 * delS = -b1
      // c1 * delR + c2 * delS = -b0
      // Onde ci e bi são os coeficientes calculados acima.
      // Podemos aplicar a regra de cramer para resolver esse sistema

      double det = cCoeffs[2] * cCoeffs[2] - cCoeffs[1] * cCoeffs[3];
      if (std::abs(det) < 1e-12) {
        // Caso de determinante nulo (singularidade), aplicamos um pequeno salto
        delR = 0.01;
        delS = 0.01;
      } else {
        // Solução via Regra de Cramer
        delR = (-bCoeffs[1] * cCoeffs[2] + bCoeffs[0] * cCoeffs[3]) / det;
        delS = (-bCoeffs[0] * cCoeffs[2] + bCoeffs[1] * cCoeffs[1]) / det;
      }

      // Atualização de r e s
      r += delR;
      s += delS;

      // Cálculo do erro relativo para critério de parada
      double r_den = (r != 0.0) ? r : 1e-10;
      double s_den = (s != 0.0) ? s : 1e-10;
      err = std::max(std::abs(delR / r_den), std::abs(delS / s_den));
    }

    // Se atingir o limite de iterações sem convergir, interrompe o processo
    if (iteration >= max_iterations) {
      process_success = false;
      break;
    }

    // Registra o primeiro fator encontrado (para plotagem de bacias no fractal)
    if (!first_factor_found) {
      first_r_conv = r;
      first_s_conv = s;
      if (iteration <= 1) { // Caso a convergência seja imediata
        exact_root = true;
      }
      first_factor_found = true;
    }

    // Extrai as duas raízes do fator quadrático x^2 - rx - s
    localRoots = get_quadratic_roots(1, -r, -s);
    roots.push_back(localRoots[0]);
    roots.push_back(localRoots[1]);

    // Deflação do polinômio: o novo polinômio é dado pelos coeficientes b
    // descartando-se os dois primeiros (que seriam o resto, agora zero)
    aCoeffs = std::vector<double>(bCoeffs.begin() + 2, bCoeffs.end());
    n = aCoeffs.size() - 1;
  }

  // Tratamento do polinômio restante (grau 2 ou 1)
  if (n == 2) {
    // Equação do segundo grau direta
    localRoots = get_quadratic_roots(aCoeffs[2], aCoeffs[1], aCoeffs[0]);
    roots.push_back(localRoots[0]);
    roots.push_back(localRoots[1]);
  } else if (n == 1) {
    // Equação do primeiro grau direta: a1*x + a0 = 0 -> x = -a0/a1
    roots.push_back((-aCoeffs[0]) / aCoeffs[1]);
  }

  // Preenchimento do resultado final
  bairstowResult final_result;
  final_result.roots = roots;
  final_result.total_iterations = totalIterations;
  final_result.success = process_success;
  final_result.r_converged = first_r_conv;
  final_result.s_converged = first_s_conv;
  final_result.exact_root = exact_root;

  return final_result;
}
